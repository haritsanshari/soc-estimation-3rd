function [sys,x0,str,ts,simStateCompliance] = TAU_ONLY(t,x,u,flag)
% TAU_ONLY  S-Function Level-1 (MATLAB)
%   Mengonversi theta1..theta8 → tau1..tau3 dan
%   MENCETAK nilai real + imag setiap langkah.
% Kode ini digunakan untuk menguji validitas dari hasil akar pers. kubik

    switch flag
        case 0
            [sys,x0,str,ts,simStateCompliance] = mdlInitializeSizes;
        case 3
            sys = mdlOutputs(t,x,u);     % <-- kirim t ke mdlOutputs
        otherwise
            sys = [];
    end
end

%% ---------------------------------------------------------------------
function [sys,x0,str,ts,simStateCompliance] = mdlInitializeSizes
    sizes               = simsizes;
    sizes.NumContStates = 0;
    sizes.NumDiscStates = 0;
    sizes.NumOutputs    = 3;    % tau1..3
    sizes.NumInputs     = 8;    % theta1..8
    sizes.DirFeedthrough= 1;
    sizes.NumSampleTimes= 1;

    sys = simsizes(sizes);
    x0  = [];
    str = [];
    ts  = [0 0];                % inherit Ts
    simStateCompliance = 'DefaultSimState';
end

%% ---------------------------------------------------------------------
function sys = mdlOutputs(t,~,u)          % <-- t dipakai untuk log
    % --------- theta ----------
    theta = u(:);
    th2 = theta(2);  th3 = theta(3);  th4 = theta(4);

    % --------- koefisien polinom ----------
    T   = 1;
    D   = 1 - th2 - th3 - th4;
    p3  = 1;
    p2  = -(T/2)*(3 - th2 + th3 + 3*th4)/D;
    p1  =  (T^2/4)*(3 + th2 + th3 - 3*th4)/D;
    p0  = -(T^3/8)*(1 + th2 - th3 + th4)/D;

    tau_all = roots([p3 p2 p1 p0]);      % bisa kompleks
    tau_r   = real(tau_all);

    % --------- CETAK real + imag ----------
    fprintf('t = %-6.1f  tau = [%.4g%+.4gi, %.4g%+.4gi, %.4g%+.4gi]\n',...
            t,...
            real(tau_all(1)), imag(tau_all(1)),...
            real(tau_all(2)), imag(tau_all(2)),...
            real(tau_all(3)), imag(tau_all(3)));

    % kembalikan hanya bagian real (3×1)
    sys = tau_r(1:3);
end

