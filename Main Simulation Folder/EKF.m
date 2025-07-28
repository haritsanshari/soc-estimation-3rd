function [sys,x0,str,ts] = EKF(t,x,u,flag)
% EKF –  Level-1 S-Function
%  EKF untuk estimasi State-of-Charge (SoC) baterai Li-ion
%  menggunakan Equivalent-Circuit-Model (ECM) orde-3.
%
%  -------- Input port (width = 9) --------
%   u(1) : I      – arus (A), arah ± mengosongkan (+ discharge)
%   u(2) : Vt_meas – tegangan terminal terukur (V)
%   u(3) : R0  (Ohm)
%   u(4) : R1  (Ohm)
%   u(5) : R2  (Ohm)
%   u(6) : R3  (Ohm)
%   u(7) : C1  (F)
%   u(8) : C2  (F)
%   u(9) : C3  (F)
%
%  -------- Output port (width = 5) --------
%   y(1) : SoC_hat          (0–1)
%   y(2) : Vt_hat (prediksi EKF, V)
%   y(3) : V1_hat - tegangan RC1 (V)
%   y(4) : V2_hat - tegangan RC2 (V)
%   y(5) : V3_hat - tegangan RC3 (V)
%
%  Internal discrete states (nx = 4 + 16 = 20):
%     x_hat(4)  : SoC, V1, V2, V3
%     P(4x4)    : kovarians diratakan (column-major)
%
%  © 2025 Harits_Estimator SoC
% -------------------------------------------------------------------------

% ------------------ konfigurasi umum -------------------------------------
Ts      = 1;        % (s)  – selang waktu sampling (ubah sesuai model)
Qbatt   = 1.1*3600; % (C)  – kapasitas nominal sel (mis. 5 Ah)
eta     = 1;      % coulombic efficiency
Qproc   = diag([1e-8 1e-5 1e-5 1e-5]);   % proses-noise (tuned)
Rmeas   = 1e-4;                          % measurement-noise (V^2)

switch flag
% -------------------------------------------------------------------------
case 0      % mdlInitializeSizes
% -------------------------------------------------------------------------
    sizes                = simsizes;
    sizes.NumContStates  = 0;
    sizes.NumDiscStates  = 20;
    sizes.NumOutputs     = 5;
    sizes.NumInputs      = 9;
    sizes.DirFeedthrough = 1;
    sizes.NumSampleTimes = 1;
    sys = simsizes(sizes);

    % initial state & P (SoC=1, tegangan RC = 0, P = diag besar)
    x_hat0 = [0.9; 0; 0; 0];
    P0     = diag([1e-4 1e-2 1e-2 1e-2]);
    x0     = [x_hat0; P0(:)];
    str    = [];
    ts     = [Ts 0];

% -------------------------------------------------------------------------
case 2      % mdlUpdate : EKF core
% -------------------------------------------------------------------------
    % ----------- ambil sinyal masukan ------------------------------------
    I        = u(1);          % arus (A)
    Vt_meas  = u(2);          % tegangan terminal (V)
    R0       = u(3);  R1 = u(4);  R2 = u(5);  R3 = u(6);
    C1       = u(7);  C2 = u(8);  C3 = u(9);

    % ----------- pecah vektor state --------------------------------------
    x_hat    = x(1:4);                     % [SoC; V1; V2; V3]
    P        = reshape(x(5:end),4,4);      % kovarians 4x4

    % ===== 1) PREDIKSI ===================================================
    A1 = exp(-Ts/(R1*C1));
    A2 = exp(-Ts/(R2*C2));
    A3 = exp(-Ts/(R3*C3));

    % fungsi transisi nonlinear
    f = [ x_hat(1) - eta*Ts*I/Qbatt;
          A1*x_hat(2) + R1*(1-A1)*I;
          A2*x_hat(3) + R2*(1-A2)*I;
          A3*x_hat(4) + R3*(1-A3)*I ];

    % Jacobian Fk = df/dx
    Fk = [1,     0,     0,     0;
          0,    A1,     0,     0;
          0,     0,    A2,     0;
          0,     0,     0,    A3];

    % kovarian prediksi
    Pm = Fk*P*Fk' + Qproc;

    % ===== 2) PREDIKSI OUTPUT ============================================
    Voc  = open_circuit_voltage(f(1));
    Vt_hat = Voc - I*R0 - f(2) - f(3) - f(4);

    % ===== 3) KOREKSI ====================================================
    Hk = [ Voc_grad(f(1)), -1, -1, -1 ];   % Jacobian dh/dx
    y_tilde = Vt_meas - Vt_hat;            % inovasi

    S   = Hk*Pm*Hk' + Rmeas;               % kovarian inovasi   (1x1)
    K   = (Pm*Hk')/S;                      % gain Kalman        (4x1)

    x_hat_new = f + K*y_tilde;             % update state
    P_new     = (eye(4) - K*Hk)*Pm;        % update kovarian

    % jaga SoC 0-1
    x_hat_new(1) = min(max(x_hat_new(1),0),1);

    % simpan ke vektor diskret
    sys = [x_hat_new; P_new(:)];

% -------------------------------------------------------------------------
case 3      % mdlOutputs
% -------------------------------------------------------------------------
    % keluarkan estimasi terbaru
    x_hat = x(1:4);
    Voc   = open_circuit_voltage(x_hat(1));
    I     = u(1);
    R0    = u(3);
    Vt_hat = Voc - I*R0 - x_hat(2) - x_hat(3) - x_hat(4);

    sys = [x_hat(1); Vt_hat; x_hat(2:4)];   % [SoC; Vt; V1 V2 V3]

% -------------------------------------------------------------------------
case {1,4,9}   % mdlDerivatives / mdlTerminate / unhandled
% -------------------------------------------------------------------------
    sys = [];

otherwise
    error(['Unhandled flag = ',num2str(flag)]);
end
end

% =====================================================================
% ----  VOC(SoC) & turunannya  ----------------------------------------
function V = open_circuit_voltage(SoC)
% Koefisien polinomial orde-8 (dari data Anda)
% p = [ 0.5180e3  -1.9091e3   2.8322e3  -2.1626e3  0.9016e3  -0.2006e3  0.0217e3  0.0024e3];
p = [-1.0979e3  4.9095e3  -9.0725e3  8.9523e3 ...
         -5.0861e3  1.6719e3  -3.0330e2  2.7200e1  2.3000];
V = polyval(p,SoC);
end

function dV = Voc_grad(SoC)
% p = [ 0.5180e3  -1.9091e3   2.8322e3  -2.1626e3  0.9016e3  -0.2006e3  0.0217e3  0.0024e3];
p = [-1.0979e3  4.9095e3  -9.0725e3  8.9523e3 ...
         -5.0861e3  1.6719e3  -3.0330e2  2.7200e1  2.3000];
dp = polyder(p);
dV = polyval(dp,SoC);
end
