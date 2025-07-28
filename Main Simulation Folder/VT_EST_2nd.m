function [sys,x0,str,ts,simStateCompliance] = VT_EST_2nd(t,x,u,flag)
% VT_EST2  –  Estimasi Vt ECM orde-2  (Level-1 S-Function)
%
% INPUT  u(1:7) :
%   1  Voc(k)   – open-circuit voltage            [V]
%   2  R0       – ohmic resistance                [Ω]
%   3  R1       – resistansi cabang RC-1          [Ω]
%   4  R2       – resistansi cabang RC-2          [Ω]
%   5  C1       – kapasitansi cabang RC-1         [F]
%   6  C2       – kapasitansi cabang RC-2         [F]
%   7  I(k)     – arus sel positif = discharge    [A]
%
% STATE diskret x :
%   x(1) = V_RC1(k)   – tegangan elemen RC-1
%   x(2) = V_RC2(k)   – tegangan elemen RC-2
%
% OUTPUT :
%   sys = Vt_hat(k) = Voc − R0·I − (V_RC1 + V_RC2)
% -------------------------------------------------------------------------

switch flag
    case 0,  [sys,x0,str,ts,simStateCompliance] = mdlInitializeSizes;
    case 2,  sys = mdlUpdate(x,u);        % x(k) -> x(k+1)
    case 3,  sys = mdlOutputs(x,u);       % hitung Vt̂(k)
    otherwise, sys = [];                  % flags 1,4,9 tidak dipakai
end
end

%% =============== INIT ====================================================
function [sys,x0,str,ts,simStateCompliance] = mdlInitializeSizes
Tsample = 1;                      %#ok<NASGU> % ubah di sini jika perlu

S = simsizes;
S.NumContStates  = 0;
S.NumDiscStates  = 2;             % V_RC1 V_RC2
S.NumOutputs     = 1;             % Vt_hat
S.NumInputs      = 7;             % Voc R0 R1 R2 C1 C2 I
S.DirFeedthrough = 1;             % output memakai u(k)
S.NumSampleTimes = 1;

sys = simsizes(S);
x0  = [0; 0];                     % inisialisasi V_RC = 0
str = [];
ts  = [1 0];                      % Ts = 1 s (sesuaikan jika perlu)
simStateCompliance = 'DefaultSimState';
end

%% =============== UPDATE (flag==2) ========================================
function x_next = mdlUpdate(x,u)
T = 1;                             % sample-time [s]  – samakan dgn ts(1)

% Ambil parameter fisik dan arus
R1 = u(3);  R2 = u(4);
C1 = u(5);  C2 = u(6);
I  = u(7);

% Lindungi τ = R·C dari nol/negatif
tau1 = max(R1*C1, eps);
tau2 = max(R2*C2, eps);

% Koefisien diskret
alpha1 = exp(-T / tau1);           beta1 = R1 * (1 - alpha1);
alpha2 = exp(-T / tau2);           beta2 = R2 * (1 - alpha2);

% State lama
V1 = x(1);  V2 = x(2);

% Update state k→k+1
V1_next = alpha1 * V1 + beta1 * I;
V2_next = alpha2 * V2 + beta2 * I;

x_next = [V1_next; V2_next];
end

%% =============== OUTPUT (flag==3) ========================================
function Vt_hat = mdlOutputs(x,u)
Voc = u(1);   R0 = u(2);   I = u(7);
Vt_hat = Voc - R0 * I - sum(x);    % x = [V_RC1; V_RC2]
end
