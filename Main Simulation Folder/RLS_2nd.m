function [sys,x0,str,ts,simStateCompliance] = RLS_2nd(t,x,u,flag)
% RLS14  –  Recursive Least-Squares klasik untuk ECM orde-2
%
% INPUT  u(1:6) : [ I(k)  I(k-1)  I(k-2)  Vt(k-1)  Vt(k-2)  Vt(k) ]   (A,V)
% OUTPUT y(1:14):
%   ├─ θ(1:6)       – parameter regresi
%   ├─ y(7)         – error  e(k) = Vt(k) - Vt_hat(k)
%   ├─ y(8)         – Vt_hat(k)
%   └─ y(9:14)      – diag{P(k)}  ≡  P11 … P66
%
%   Model linier diskrét (orde-2):
%     Vt(k) = θ1 + θ2 Vt(k-1) + θ3 Vt(k-2) + θ4 I(k) + θ5 I(k-1) + θ6 I(k-2)

switch flag
    case 0
        [sys,x0,str,ts,simStateCompliance] = mdlInitializeSizes;
    case 2
        sys = mdlUpdate(t,x,u);
    case 3
        sys = mdlOutputs(t,x,u);
    case {1,4,9}
        sys = [];          % Unused flags
    otherwise
        error(['Unhandled flag = ', num2str(flag)]);
end
end

%% =============== INITIALISATION =========================================
function [sys,x0,str,ts,simStateCompliance] = mdlInitializeSizes
nTheta = 6;                       % jumlah parameter
nState = nTheta^2 + nTheta + 1;   % P(:) + θ + simpan error

S = simsizes;
S.NumContStates  = 0;
S.NumDiscStates  = nState;
S.NumOutputs     = nTheta + 2 + nTheta;   % 6 θ + err + ŷ + 6 diag(P) = 14
S.NumInputs      = 6;                      % Ik,Ik-1,Ik-2,Vk-1,Vk-2,Vk
S.DirFeedthrough = 1;
S.NumSampleTimes = 1;

sys = simsizes(S);

P0     = 1e3 * eye(nTheta);       % kovarian awal besar → konvergensi cepat
theta0 = zeros(nTheta,1);
err0   = 0;

x0 = [reshape(P0,[],1); theta0; err0];

str = [];
ts  = [1 0];                      % Ts = 1 s (ubah jika perlu)
simStateCompliance = 'DefaultSimState';
end

%% ================= UPDATE (flag==2) =====================================
function sys = mdlUpdate(~,x,u)
nTheta  = 6;
lambda  = 0.99;                   % faktor pelupa (adjustable)

% -------- Regresor ψ(k) --------------------------------------------------
% ψ = [1; Vt(k-1); Vt(k-2); I(k); I(k-1); I(k-2)];
psi = [1; u(4); u(5); u(1); u(2); u(3)];

% -------- State sebelumnya ----------------------------------------------
P     = reshape(x(1:nTheta^2), nTheta, nTheta);
theta = x(nTheta^2+1 : nTheta^2+nTheta);

% -------- RLS ------------------------------------------------------------
err     = u(6) - psi.' * theta;                 % error prediksi
gainDen = lambda + psi.' * P * psi;             % skalar penyebut
K       = (P * psi) / gainDen;                  % vektor gain RLS

theta_n = theta + K * err;                      % pembaruan parameter
P_n     = (P - K * psi.' * P) / lambda;         % pembaruan kovarian
P_n     = (P_n + P_n.')/2;                      % paksa simetri numerik

% -------- Simpan state ---------------------------------------------------
sys = [reshape(P_n,[],1); theta_n; err];
end

%% ================= OUTPUT (flag==3) =====================================
function sys = mdlOutputs(~,x,u)
nTheta = 6;

% Rekonstruksi-------------------------------------------------------------
psi   = [1; u(4); u(5); u(1); u(2); u(3)];
P     = reshape(x(1:nTheta^2), nTheta, nTheta);
theta = x(nTheta^2+1 : nTheta^2+nTheta);

Vt_hat = psi.' * theta;
err    = x(end);

% Output vector------------------------------------------------------------
sys = [theta; err; Vt_hat; diag(P)];

% Validasi---------------------------------------------------------------
if ~isreal(sys) || any(~isfinite(sys)) || numel(sys) ~= 14
    error('Output invalid: NaN/Inf/non-real atau size ≠ 14');
end
end
