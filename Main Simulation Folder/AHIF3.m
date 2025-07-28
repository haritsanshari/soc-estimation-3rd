function [sys,x0,str,ts] = AHIF3(t,x,u,flag)
% AHIF3_ECM3  ─  Adaptive H-Infinity Filter (Level-1 S-Function, ECM-3)
%   Estimasi SoC dan tegangan RC (V1-V3) dengan adaptasi R, Q, dan theta
%   berdasarkan inovasi error terkini (moving error).
%   Mengimplementasikan modifikasi dari HIF3.m dengan mekanisme adaptif.
% -------------------------------------------------------------------------

% ------------------- parameter tetap -------------------------------------
Ts      = 1;                          % [s]  sample-time
Qbatt   = 1.1*3600;                   % [C]  kapasitas sel (≈1.1 Ah)
eta     = 1;                          % coulombic efficiency
theta0  = 1e-2;                       % initial H∞ performance bound
N       = 120;                        % panjang jendela error

Qproc   = diag([1e-10, 1e-5, 1e-5, 1e-5]); % inisiasi awal Q_k
Rmeas0  = 1e-4;                       % nilai awal R_k
Sweight = eye(4);                    % S_k  (state-error weight)
Lk      = eye(4);                    % ε_k = L_k x_k

persistent ek_win;
if isempty(ek_win)
    ek_win = zeros(N,1);
end

% ------------------- struktur S-Function ---------------------------------
switch flag
case 0
    sizes = simsizes;
    sizes.NumContStates  = 0;
    sizes.NumDiscStates  = 20;         % 4 state + 16 kovarian
    sizes.NumOutputs     = 5;
    sizes.NumInputs      = 9;
    sizes.DirFeedthrough = 1;
    sizes.NumSampleTimes = 1;
    sys = simsizes(sizes);

    xhat0 = [0.9; 0; 0; 0];            % SoC≈0.9 , tegangan RC=0
    P0    = diag([1e-2, 1e-2, 1e-2, 1e-2]);
    x0    = [xhat0; P0(:)];

    str = [];
    ts  = [Ts 0];

case 2
    I  = u(1); Vt_meas = u(2);
    R0 = u(3); R1 = u(4); R2 = u(5); R3 = u(6);
    C1 = u(7); C2 = u(8); C3 = u(9);

    xhat = x(1:4);
    P    = reshape(x(5:end),4,4);

    % Time Update
    A1 = exp(-Ts/(R1*C1));
    A2 = exp(-Ts/(R2*C2));
    A3 = exp(-Ts/(R3*C3));
    Fk = diag([1, A1, A2, A3]);

    xhat_m = [ xhat(1) - eta*Ts*I/Qbatt;
               A1*xhat(2) + R1*(1-A1)*I;
               A2*xhat(3) + R2*(1-A2)*I;
               A3*xhat(4) + R3*(1-A3)*I ];

    Pm = Fk*P*Fk' + Qproc;

    dVoc = Voc_grad(xhat_m(1));
    Hk   = [dVoc, -1, -1, -1];

    Voc_m  = open_circuit_voltage(xhat_m(1));
    yhat_m = Voc_m - I*R0 - sum(xhat_m(2:4));
    ek     = Vt_meas - yhat_m;

    ek_win = [ek_win(2:end); ek];
    Mk = mean(ek_win.^2);

    CPCT = Hk * Pm * Hk';
    Rk = Mk - CPCT;
    Rk = max(min(Rk, 1e-3), 1e-6);
    Rinv = 1/Rk;

    Sk_bar = Lk' * Sweight * Lk;
    theta = theta0;
    condMat = inv(Pm) - theta * Sk_bar + Hk' * Rinv * Hk;

    while any(eig(condMat) <= 0)
        theta = theta * 2;
        condMat = inv(Pm) - theta * Sk_bar + Hk' * Rinv * Hk;
    end

    Gamma = eye(4) - theta * Sk_bar * Pm + (Hk' * Rinv * Hk) * Pm;
    Ginv = inv(Gamma + 1e-6 * eye(4));
    Kk = Fk * Pm * Ginv * Hk' * Rinv;

    dx = Kk * ek;
    limdx = [3*abs(I*Ts/Qbatt);
             3*abs(R1*(1-A1)*I);
             3*abs(R2*(1-A2)*I);
             3*abs(R3*(1-A3)*I)];
    dx = min(max(dx, -limdx), limdx);

    xhat_p = xhat_m + dx;
    xhat_p(1) = min(max(xhat_p(1),0),1);
    xhat_p(2:4) = min(max(xhat_p(2:4), -10), 10);

    Qadapt = Kk * Mk * Kk';
    Qadapt = min(max(Qadapt, -1e-5), 1e-4);
    Qproc = Qadapt;

    Pp = Pm * Ginv;
    sys = [xhat_p; Pp(:)];

case 3
    xhat = x(1:4);
    I = u(1); R0 = u(3);

    Voc = open_circuit_voltage(xhat(1));
    Vt_hat = Voc - I*R0 - sum(xhat(2:4));

    sys = [xhat(1); Vt_hat; xhat(2:4)];

case {1,4,9}
    sys = [];

otherwise
    error(['AHIF3_ECM3 : unhandled flag = ' num2str(flag)]);
end
end

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