function [sys,x0,str,ts,simStateCompliance] = RLS(t,x,u,flag)
% RLS klasik untuk ECM Orde-3 dengan 8 parameter regresi
% Urutan input: [I(k); I(k-1); I(k-2); I(k-3); Vt(k-1); Vt(k-2); Vt(k-3); Vt(k)]
% Output: [theta(1:8); error; yhat]

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

% =========== INIT ============
function [sys,x0,str,ts,simStateCompliance] = mdlInitializeSizes
nTheta = 8;
nState = nTheta^2 + nTheta + 1;  % P (64), theta (8), error (1)

sizes = simsizes;
sizes.NumContStates  = 0;
sizes.NumDiscStates  = nState;
sizes.NumOutputs     = nTheta + 2;  % theta(8), error, yhat
sizes.NumInputs      = 8;           % [I(k:k-3); Vt(k-1:k-3); Vt(k)]
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;

sys = simsizes(sizes);

P0 = 10 * eye(nTheta);
theta0 = zeros(nTheta,1);
err0 = 0;

x0 = [reshape(P0,[],1); theta0; err0];

str = [];
ts  = [1 0];
simStateCompliance = 'DefaultSimState';
end

% =========== UPDATE ============
function sys = mdlUpdate(t,x,u)
nTheta = 8;

% === Input unpacking ===
i_vals    = u(1:4);      % I(k), I(k-1), I(k-2), I(k-3)
vt_delays = u(5:7);      % Vt(k-1), Vt(k-2), Vt(k-3)
vt_k      = u(8);        % Vt(k)

% === Construct regressor vector ===
psi = [1; vt_delays(:); i_vals(:)];

% === Retrieve previous state ===
P = reshape(x(1:nTheta^2), nTheta, nTheta);
theta = x(nTheta^2+1 : nTheta^2+nTheta);

% === RLS classical update ===
lambda = 0.99;  % Fixed forgetting factor
error = vt_k - psi' * theta;

gain_denom = lambda + psi' * P * psi;
P_new = (1 / lambda) * (P - (P * (psi * psi') * P) / gain_denom);
theta_new = theta + P_new * psi * error;    

% === Update state ===
sys = [reshape(P_new,[],1); theta_new; error];
end

% =========== OUTPUT ============
function sys = mdlOutputs(t,x,u)
nTheta = 8;

% === Input unpacking ===
i_vals    = u(1:4);      % I(k), I(k-1), ...
vt_delays = u(5:7);      % Vt(k-1), ...
vt_k      = u(8);        % current target Vt

% === Regressor ===
psi = [1; vt_delays(:); i_vals(:)];
theta = x(nTheta^2+1 : nTheta^2+nTheta);
yhat = psi' * theta;
error = x(end);

% === Output ===
sys = [theta; error; yhat];

% Validate output
if ~isreal(sys) || any(isnan(sys)) || any(isinf(sys)) || numel(sys) ~= 10
    error('Output invalid: NaN, Inf, atau bukan vektor real ukuran 10');
end
end
