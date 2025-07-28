function [sys,x0,str,ts,simStateCompliance] = RLS18(t,x,u,flag)
% RLS klasik untuk ECM Orde-3 (8 parameter regresi)
% Input  u = [I(k); I(k-1); I(k-2); I(k-3); Vt(k-1); Vt(k-2); Vt(k-3); Vt(k)]
% Output y = [theta(1:8); error; yhat; P(11)-P(88)]

switch flag
    case 0
        [sys,x0,str,ts,simStateCompliance] = mdlInitializeSizes;
    case 2
        sys = mdlUpdate(t,x,u);
    case 3
        sys = mdlOutputs(t,x,u);
    case {1,4,9}
        sys = [];
    otherwise
        error(['Unhandled flag = ', num2str(flag)]);
end
end

%% ========== INIT ==========
function [sys,x0,str,ts,simStateCompliance] = mdlInitializeSizes
nTheta = 8;
nState = nTheta^2 + nTheta + 1;   % P(:) + theta + error

sizes = simsizes;
sizes.NumContStates  = 0;
sizes.NumDiscStates  = nState;
sizes.NumOutputs     = nTheta + 2 + nTheta;   % 8 θ + err + yhat + 8 diag-P
sizes.NumInputs      = 8;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;

sys = simsizes(sizes);

P0     = 1000 * eye(nTheta);
theta0 = zeros(nTheta,1);
err0   = 0;

x0 = [reshape(P0,[],1); theta0; err0];

str = [];
ts  = [1 0];
simStateCompliance = 'DefaultSimState';
end

%% ========== UPDATE ==========
function sys = mdlUpdate(~,x,u)
nTheta = 8;

% — Regresor —
psi = [1; u(5:7); u(1:4)];  % [1; Vt(k-1:-3); I(k:-3)]

% — Ambil state sebelumnya —
P     = reshape(x(1:nTheta^2), nTheta, nTheta);
theta = x(nTheta^2+1 : nTheta^2+nTheta);

% — RLS —
lambda = 0.99;                         % forgetting factor
error  = u(8) - psi' * theta;          % Vt(k) − ŷ
gainDen = lambda + psi' * P * psi;
P_new   = (1/lambda)*(P - (P*(psi*psi')*P)/gainDen);
theta_new = theta + P_new*psi*error;

% — Simpan state —
sys = [reshape(P_new,[],1); theta_new; error];
end

%% ========== OUTPUT ==========
function sys = mdlOutputs(~,x,u)
nTheta = 8;

% — Re-hitung regressor & prediksi —
psi   = [1; u(5:7); u(1:4)];
P     = reshape(x(1:nTheta^2), nTheta, nTheta);
theta = x(nTheta^2+1 : nTheta^2+nTheta);

yhat  = psi' * theta;
error = x(end);

% — Bangun vektor keluaran —
sys = [theta; error; yhat; diag(P)];

% Validasi ukuran & realitas
if ~isreal(sys) || any(~isfinite(sys)) || numel(sys) ~= 18
    error('Output invalid: NaN/Inf/non-real atau ukuran ≠ 18');
end
end
