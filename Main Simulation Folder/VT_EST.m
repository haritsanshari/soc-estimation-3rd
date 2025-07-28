function [sys,x0,str,ts,simStateCompliance] = VT_EST(t,x,u,flag)
% VT_EST S-Function Level-1: Estimasi Vt (terminal voltage) dari parameter fisik
%
% INPUT (u):
%   u(1) = Voc
%   u(2) = R0
%   u(3) = R1
%   u(4) = R2
%   u(5) = R3
%   u(6) = C1
%   u(7) = C2
%   u(8) = C3
%   u(9) = I(k)   (arus pada saat k)
%
% STATES diskret (x):
%   x(1) = V_RC1  (tegangan di cabang RC1 pada langkah k)
%   x(2) = V_RC2  (tegangan di cabang RC2 pada langkah k)
%   x(3) = V_RC3  (tegangan di cabang RC3 pada langkah k)
%
% OUTPUT:
%   sys = Vt_est = Voc - R0*I(k) - [V_RC1 + V_RC2 + V_RC3]
%
% Asumsi: Sample time T = 1. Jika perlu diubah, sesuaikan nilai T pada mdlInitializeSizes
% dan di mdlUpdate.
%

switch flag
    case 0
        % Inisialisasi ukuran‐ukuran (sizes)
        [sys,x0,str,ts,simStateCompliance] = mdlInitializeSizes;
    case 2
        % Update state diskret (x_k -> x_{k+1})
        sys = mdlUpdate(t,x,u);
    case 3
        % Hitung output (Vt_est) berdasarkan state x dan input u
        sys = mdlOutputs(t,x,u);
    otherwise
        sys = [];
end
end

% ====================================================
function [sys,x0,str,ts,simStateCompliance] = mdlInitializeSizes
    sizes = simsizes;
    sizes.NumContStates   = 0;    % Tidak ada continuous states
    sizes.NumDiscStates   = 3;    % Tiga state diskret: V_RC1, V_RC2, V_RC3
    sizes.NumOutputs      = 1;    % Hanya satu output: Vt_est
    sizes.NumInputs       = 9;    % Sembilan input: Voc, R0–R3, C1–C3, I(k)
    sizes.DirFeedthrough  = 1;    % OUTPUT tidak langsung feedthrough dari u (karena butuh x)
    sizes.NumSampleTimes  = 1;    % Satu sample time diskret

    sys = simsizes(sizes);
    x0  = [0; 0; 0];   % Inisialisasi V_RC1=0, V_RC2=0, V_RC3=0 pada k=0
    str = [];
    ts  = [1 0];      % Sample time diskret T = 1 (ubah jika perlu)
    simStateCompliance = 'DefaultSimState';
end

% ====================================================
function xkp1 = mdlUpdate(t,x,u)
    % Update state diskret: x_k = [V_RC1(k); V_RC2(k); V_RC3(k)]
    % Kita hitung x_{k+1} = alpha .* x_k + beta .* I(k)
    %
    % Input u:
    %   u(3)=R1, u(4)=R2, u(5)=R3, u(6)=C1, u(7)=C2, u(8)=C3, u(9)=I_k
    %
    T = 1;                            % Langkah diskret (satuan waktu)
    R1 = u(3);  R2 = u(4);  R3 = u(5);
    C1 = u(6);  C2 = u(7);  C3 = u(8);
    I_k = u(9);

    % Hitung tau_i = R_i * C_i
    tao1 = R1 * C1;
    tao2 = R2 * C2;
    tao3 = R3 * C3;

    % Hitung alpha_i = exp(-T/tao_i), beta_i = R_i * (1 - alpha_i)
    alpha1 = exp(- T / (tao1));
    alpha2 = exp(- T / (tao2));
    alpha3 = exp(- T / (tao3));
    beta1  = R1 * (1 - alpha1);
    beta2  = R2 * (1 - alpha2);
    beta3  = R3 * (1 - alpha3);

    V_RC = x(:);  % V_RC pada langkah k
    % Update ke langkah k+1:
    V_RC1_next = alpha1 * V_RC(1) + beta1  * I_k;
    V_RC2_next = alpha2 * V_RC(2) + beta2  * I_k;
    V_RC3_next = alpha3 * V_RC(3) + beta3  * I_k;

    xkp1 = [V_RC1_next; V_RC2_next; V_RC3_next];
end

% ====================================================
function Vt_est = mdlOutputs(t,x,u)
    % Hitung output Vt_est = Voc - R0*I - (V_RC1 + V_RC2 + V_RC3)
    %
    % Input u:
    %   u(1)=Voc, u(2)=R0, u(9)=I_k
    Voc = u(1);
    R0  = u(2);
    I_k = u(9);

    V_RC = x(:);  % V_RC pada langkah k

    Vt_est = Voc - R0 * I_k - sum(V_RC);
end
