function [sys,x0,str,ts] = HIF3(t,x,u,flag)
% HIF3_ECM3  ─  H∞ FILTER  (Level-1 S-Function, ECM-3, discrete)
%   Mengimplementasikan langkah (6)–(13) pada Table-I H-Infinity Filter
%   untuk estimasi SoC & tegangan RC battery Li-ion.
%
% INPUT  u(1:9) :  I  Vt  R0-R3  C1-C3   (arus positif = discharge)
% OUTPUT sys    : [SoC_hat  Vt_hat  V1_hat  V2_hat  V3_hat]
%
% © 2025  Harits_Estimator SoC
% -------------------------------------------------------------------------

% ------------------- parameter tetap -------------------------------------
Ts      = 1;                          % [s]  sample-time
Qbatt   = 1.1*3600;                   % [C]  kapasitas sel (≈1.1 Ah)
eta     = 1;                          % coulombic efficiency
theta0  = 1e-2;                       % H∞ performance bound (bisa adaptif)

Qproc   = diag([1e-10, 1e-5, 1e-5, 1e-5]); % Q_k  (process noise)
Rmeas0  = 1e-4;                       % R_k  (meas noise)  → skalar
Sweight = eye(4);                     % S_k  (state-error weight)
Lk      = eye(4);                     % ε_k = L_k x_k  (monitor semua state)

% ------------------- struktur S-Function ---------------------------------
switch flag
% -------------------------------------------------------------------------
case 0                      % ============ mdlInitializeSizes ============
% -------------------------------------------------------------------------
    sizes = simsizes;
    sizes.NumContStates  = 0;
    sizes.NumDiscStates  = 20;         % 4 state + 16 unsur kovarian
    sizes.NumOutputs     = 5;
    sizes.NumInputs      = 9;
    sizes.DirFeedthrough = 1;
    sizes.NumSampleTimes = 1;
    sys  = simsizes(sizes);

    xhat0 = [0.9; 0; 0; 0];            % SoC≈0.9 , tegangan RC=0
    P0    = diag([1e-4, 1e-2, 1e-2, 1e-2]);
    x0    = [xhat0; P0(:)];

    str = [];
    ts  = [Ts 0];

% -------------------------------------------------------------------------
case 2                      % ============ mdlUpdate (core) ==============
% -------------------------------------------------------------------------
    % ---------- (a) ambil input ------------------------------------------
    I  = u(1);  Vt_meas = u(2);
    R0 = u(3);  R1 = u(4);  R2 = u(5);  R3 = u(6);
    C1 = u(7);  C2 = u(8);  C3 = u(9);

    % ---------- (b) de-vector state --------------------------------------
    xhat = x(1:4);                       % [SoC; V1; V2; V3]
    P    = reshape(x(5:end),4,4);        % P_k⁺ (4×4)

    % =====================================================================
    % I I .   T I M E   U P D A T E   (6)–(7)
    % =====================================================================
    A1 = exp(-Ts/(R1*C1));
    A2 = exp(-Ts/(R2*C2));
    A3 = exp(-Ts/(R3*C3));
    Fk = diag([1, A1, A2, A3]);          % Jacobian f_x

    % (6) state prior
    xhat_m = [ xhat(1) - eta*Ts*I/Qbatt ;          % SoC
               A1*xhat(2) + R1*(1-A1)*I ;
               A2*xhat(3) + R2*(1-A2)*I ;
               A3*xhat(4) + R3*(1-A3)*I ];

    % (7) covariance prior
    Pm = Fk*P*Fk' + Qproc;

    % =====================================================================
    % I I I .   M E A S U R E M E N T   M O D E L
    % =====================================================================
    dVoc  = Voc_grad(xhat_m(1));
    Hk    = [dVoc, -1, -1, -1];         % Jacobian h_x
    Rk    = Rmeas0;                     % (bisa adaptif)
    Rinv  = 1/Rk;

    Voc_m  = open_circuit_voltage(xhat_m(1));
    yhat_m = Voc_m - I*R0 - sum(xhat_m(2:4));   % h(x̂⁻)
    ek     = Vt_meas - yhat_m;                 % inovasi (10)

    % =====================================================================
    % S̄_k  &   K E T A K S A M A A N  (8)–(9)
    % =====================================================================
    Sk_bar = Lk' * Sweight * Lk;        % (8)
    theta  = theta0;                    % reset / bisa buat adaptif

    condMat = inv(Pm) - theta * Sk_bar + Hk' * Rinv * Hk;   % (9)
    while any(eig(condMat) <= 0)
        theta = theta * 2;              % naikkan θ hingga syarat terpenuhi
        condMat = inv(Pm) - theta * Sk_bar + Hk' * Rinv * Hk;
    end

    % =====================================================================
    % G A I N   H - ∞   (11)
    % =====================================================================
    Gamma = eye(4) - theta * Sk_bar * Pm + (Hk' * Rinv * Hk) * Pm;
    Ginv  = Gamma \ eye(4);             % = inv(Gamma) (tanpa instabilitas)
    Kk    = Fk * Pm * Ginv * Hk' * Rinv;

    % =====================================================================
    % S T A T E   &   C O V U P D A T E (12)–(13)
    % =====================================================================
    xhat_p = xhat_m + Kk * ek;          % (12)
    Pp     = Pm * Ginv;                 % (13)

    % ---------- housekeeping ---------------------------------------------
    xhat_p(1) = min(max(xhat_p(1),0),1);  % limit SoC ke [0,1]

    % ---------- simpan kembali -------------------------------------------
    sys = [xhat_p; Pp(:)];

% -------------------------------------------------------------------------
case 3                      % ============ mdlOutputs ====================
% -------------------------------------------------------------------------
    xhat = x(1:4);                       % gunakan x̂_k⁺
    I    = u(1);  R0 = u(3);

    Voc  = open_circuit_voltage(xhat(1));
    Vt_hat = Voc - I*R0 - sum(xhat(2:4));

    sys = [xhat(1); Vt_hat; xhat(2:4)];

% -------------------------------------------------------------------------
case {1,4,9}                % -------- unused flags ----------------------
    sys = [];

otherwise
    error(['HIF3_ECM3 : unhandled flag = ' num2str(flag)]);
end
end
% ========================================================================
%                O C V ( S o C )  &  G R A D I E N T
% ========================================================================
function V = open_circuit_voltage(SoC)
p = [-1.0979e3  4.9095e3  -9.0725e3  8.9523e3 ...
      -5.0861e3  1.6719e3  -3.0330e2  2.7200e1  2.3000];
V = polyval(p, SoC);
end

function dV = Voc_grad(SoC)
p  = [-1.0979e3  4.9095e3  -9.0725e3  8.9523e3 ...
      -5.0861e3  1.6719e3  -3.0330e2  2.7200e1  2.3000];
dV = polyval(polyder(p), SoC);
end
