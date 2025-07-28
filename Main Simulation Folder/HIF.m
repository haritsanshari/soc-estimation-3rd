function [sys,x0,str,ts] = HIF(t,x,u,flag)
% H∞ Filter Level-1 S-Function for SoC Estimation (ECM Orde-3)
% INPUT:  I, Vt_meas, R0, R1, R2, R3, C1, C2, C3
% OUTPUT: SoC_hat, Vt_hat, V1_hat, V2_hat, V3_hat

Ts      = 1;             % Sampling time (s)
Qbatt   = 1.1 * 3600;    % Kapasitas baterai (Coulomb)
eta     = 1;           % Efisiensi coulombic
theta   = 0.01;          % Parameter H∞
Qproc   = diag([1e-8 1e-5 1e-5 1e-5]);
Rmeas   = 1e-4;
Sweight = eye(4);        % Pembobotan error state
L       = eye(4);        % Estimasi langsung x

switch flag
    case 0
        sizes = simsizes;
        sizes.NumContStates  = 0;
        sizes.NumDiscStates  = 20;  % xhat(4) + P(4x4)
        sizes.NumOutputs     = 5;
        sizes.NumInputs      = 9;
        sizes.DirFeedthrough = 1;
        sizes.NumSampleTimes = 1;
        sys = simsizes(sizes);

        xhat0 = [0.9; 0; 0; 0];       % SoC = 1; V1–V3 = 0
        P0 = diag([1e-4 1e-2 1e-2 1e-2]);
        x0 = [xhat0; P0(:)];
        str = [];
        ts = [Ts 0];

    case 2  % Update block
        I  = u(1); Vt = u(2); R0 = u(3);
        R1 = u(4); R2 = u(5); R3 = u(6);
        C1 = u(7); C2 = u(8); C3 = u(9);

        xhat = x(1:4);
        P = reshape(x(5:end),4,4);

        A1 = exp(-Ts/(R1*C1));  A2 = exp(-Ts/(R2*C2));  A3 = exp(-Ts/(R3*C3));
        A = diag([1 A1 A2 A3]);
        f = [xhat(1) - eta*Ts*I/Qbatt;
             A1*xhat(2) + R1*(1-A1)*I;
             A2*xhat(3) + R2*(1-A2)*I;
             A3*xhat(4) + R3*(1-A3)*I];

        % Predict
        Pm = A * P * A' + Qproc;

        % Jacobian of output
        dVoc = Voc_grad(f(1));
        Ck = [dVoc, -1, -1, -1];

        % Prediction of output
        Voc_est = open_circuit_voltage(f(1));
        yhat = Voc_est - I*R0 - f(2) - f(3) - f(4);
        ek = Vt - yhat;

        % Sk_bar
        Sk_bar = L' * Sweight * L;

        % Gain H∞
        inv_term = eye(4) - theta*Sk_bar*Pm + Ck'*Rmeas^(-1)*Ck*Pm;
        K = A * Pm * (inv(inv_term)) * Ck' * Rmeas^(-1);

        % Update state
        xhat_new = f + K * ek;

        % Update covariance
        P_new = Pm * inv(inv_term);

        xhat_new(1) = min(max(xhat_new(1), 0), 1);  % clamp SoC
        sys = [xhat_new; P_new(:)];

    case 3  % Output block
        xhat = x(1:4);
        I  = u(1); R0 = u(3);
        Voc_est = open_circuit_voltage(xhat(1));
        Vt_hat = Voc_est - I*R0 - xhat(2) - xhat(3) - xhat(4);
        sys = [xhat(1); Vt_hat; xhat(2:4)];

    case {1,4,9}
        sys = [];
    otherwise
        error(['Unhandled flag = ', num2str(flag)]);
end
end

% --------- Voc Polynomial & Gradient ----------
function V = open_circuit_voltage(SoC)
    % p = [ 0.5180e3 -1.9091e3 2.8322e3 -2.1626e3 0.9016e3 -0.2006e3 0.0217e3 0.0024e3 ];
    p = [-1.0979e3  4.9095e3  -9.0725e3  8.9523e3 ...
         -5.0861e3  1.6719e3  -3.0330e2  2.7200e1  2.3000];
    V = polyval(p, SoC);
end

function dV = Voc_grad(SoC)
    % p = [ 0.5180e3 -1.9091e3 2.8322e3 -2.1626e3 0.9016e3 -0.2006e3 0.0217e3 0.0024e3 ];
    p = [-1.0979e3  4.9095e3  -9.0725e3  8.9523e3 ...
         -5.0861e3  1.6719e3  -3.0330e2  2.7200e1  2.3000];
    dp = polyder(p);
    dV = polyval(dp, SoC);
end
