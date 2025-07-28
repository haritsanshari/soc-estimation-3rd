function [sys,x0,str,ts] = HIF2(t,x,u,flag)
% H∞ FILTER  –  Level-1 S-Function (ECM orde-3, discrete-time)
%   Persamaan mengikuti Table-I: (6)–(13)
%
% INPUT  (width = 9):  I, Vt_meas, R0, R1, R2, R3, C1, C2, C3
% OUTPUT (width = 5):  SoC_hat, Vt_hat, V1_hat, V2_hat, V3_hat
%
% © 2025  Harits_Estimator SoC
% -------------------------------------------------------------------------

% ---------------- parameter tetap ---------------------------------------
Ts      = 1;              % [s]  sample-time
Qbatt   = 1.1*3600;       % [C]  kapasitas sel (≈1.1 Ah)
eta     = 1;              % coulombic efficiency
theta   = 1e-2;           % H∞ performance bound (can adapt later)

Qproc   = diag([1e-8 1e-5 1e-5 1e-5]);   % Qk : process-noise
Rmeas0  = 1e-4;                           % Rk : meas-noise (skalar)
Sweight = eye(4);                        % S_k (state-error weight)
Lk      = eye(4);                        % L_k (ε_k = L_k x_k) – monitor seluruh state

% ---------------- struktur S-Function -----------------------------------
switch flag
% -------------------------------------------------------------------------
case 0          % =============== mdlInitializeSizes ======================
% -------------------------------------------------------------------------
    sizes = simsizes;
    sizes.NumContStates  = 0;
    sizes.NumDiscStates  = 20;         % x̂(4) + vec(P)(16)
    sizes.NumOutputs     = 5;
    sizes.NumInputs      = 9;
    sizes.DirFeedthrough = 1;
    sizes.NumSampleTimes = 1;
    sys  = simsizes(sizes);

    xhat0 = [0.9; 0; 0; 0];            % SoC ≈ 0.9; tegangan RC = 0
    P0    = diag([1e-4 1e-2 1e-2 1e-2]);
    x0    = [xhat0; P0(:)];

    str = [];
    ts  = [Ts 0];

% -------------------------------------------------------------------------
case 2          % =============== mdlUpdate (H∞ core) =====================
% -------------------------------------------------------------------------
    % ---------- ambil input ------------------------------------------------
    I  = u(1);  Vt_meas = u(2);
    R0 = u(3);  R1 = u(4);  R2 = u(5);  R3 = u(6);
    C1 = u(7);  C2 = u(8);  C3 = u(9);

    % ---------- de-vector state -------------------------------------------
    xhat = x(1:4);                          % [SoC; V1; V2; V3]
    P    = reshape(x(5:end),4,4);           % kovarian 4×4

    % =====================================================================
    %                I I .  T I M E  –  U P D A T E  (6)–(7)
    % =====================================================================
    % Dinamika diskrit eksak
    A1 = exp(-Ts/(R1*C1));  A2 = exp(-Ts/(R2*C2));  A3 = exp(-Ts/(R3*C3));
    Fk = diag([1, A1, A2, A3]);              % Jacobian/ state matrix

    % (6) prior state
    xhat_m = [ xhat(1) - eta*Ts*I/Qbatt;         % SoC
               A1*xhat(2) + R1*(1-A1)*I;         % V1
               A2*xhat(3) + R2*(1-A2)*I;         % V2
               A3*xhat(4) + R3*(1-A3)*I ];       % V3

    % (7) prior covariance
    Pm = Fk*P*Fk' + Qproc;

    % =====================================================================
    %             I I I .  M E A S U R E M E N T  M O D E L
    % =====================================================================
    dVoc  = Voc_grad(xhat_m(1));
    Hk    = [ dVoc, -1, -1, -1 ];           % Jacobian output
    Rk    = Rmeas0;                         % (dapat diadaptasi)
    Rinv  = 1/Rk;

    Voc_m  = open_circuit_voltage(xhat_m(1));
    yhat_m = Voc_m - I*R0 - xhat_m(2) - xhat_m(3) - xhat_m(4);   % h(x)
    ek     = Vt_meas - yhat_m;                                  % inovasi  (10)

    % =====================================================================
    %             S̄_k  &  K E T A K S A M A A N   R O B U S T   (8)–(9)
    % =====================================================================
    Sk_bar = Lk' * Sweight * Lk;            % (8)

    % (9) check: (P⁻)⁻¹ − θ S̄ + HᵀR⁻¹H   > 0
    test_mat = inv(Pm) - theta*Sk_bar + Hk'*Rinv*Hk;
    if any(eig(test_mat) <= 0)
        % -- tingkatkan θ atau regularisasi --------------------------------
        theta  = theta * 2;                 % strategi sederhana; bisa diganti adaptif
        warning('HIF:thetaUP','Inequality (9) violated – θ doubled to %g',theta);
        test_mat = inv(Pm) - theta*Sk_bar + Hk'*Rinv*Hk;  % ulang cek
    end

    % =====================================================================
    %                   G A I N   H ∞   (11)
    % =====================================================================
    Gamma = eye(4) - theta*Sk_bar*Pm + Hk'*Rinv*Hk*Pm;  % term dalam (I-…)
    % Hindari inv():  solve  X = Gamma \ (HᵀR⁻¹)
    Kk = Fk*Pm / Gamma * Hk' * Rinv;

    % =====================================================================
    %                   S T A T E   &   P   U P D A T E (12)–(13)
    % =====================================================================
    xhat_p = xhat_m + Kk*ek;                           % (12)
    Pp     = Pm / Gamma;                               % (13)

    % --- housekeeping -----------------------------------------------------
    xhat_p(1) = max(0,min(1,xhat_p(1)));               % saturasi SoC

    % ------------ simpan kembali -----------------------------------------
    sys = [xhat_p; Pp(:)];

% -------------------------------------------------------------------------
case 3          % =============== mdlOutputs ==============================
% -------------------------------------------------------------------------
    xhat = x(1:4);
    I    = u(1);  R0 = u(3);
    Voc  = open_circuit_voltage(xhat(1));
    Vt_hat = Voc - I*R0 - xhat(2) - xhat(3) - xhat(4);

    sys = [xhat(1); Vt_hat; xhat(2:4)];

% -------------------------------------------------------------------------
case {1,4,9}    % Unused flags
    sys = [];
otherwise
    error(['Unhandled flag = ',num2str(flag)]);
end
end
% ========================================================================
%                 V o c ( S o C )   &   g r a d i e n t
% ========================================================================
function V = open_circuit_voltage(SoC)
p = [-1.0979e3  4.9095e3  -9.0725e3  8.9523e3  ...
      -5.0861e3  1.6719e3  -3.0330e2  2.7200e1  2.3000];
V = polyval(p,SoC);
end

function dV = Voc_grad(SoC)
p  = [-1.0979e3  4.9095e3  -9.0725e3  8.9523e3  ...
      -5.0861e3  1.6719e3  -3.0330e2  2.7200e1  2.3000];
dp = polyder(p);
dV = polyval(dp,SoC);
end
