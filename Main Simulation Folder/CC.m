function [sys,x0,str,ts,simStateCompliance] = CC(t,x,u,flag)
% SOCCounter – Level-1 S-Function (MATLAB)
% Coulomb counting estimator for battery SoC
% Input : u(1) = current (A), positive = discharge
% Output: y(1) = estimated SoC (0 to 1)

switch flag
    case 0
        [sys,x0,str,ts,simStateCompliance] = mdlInitializeSizes;
    case 2
        sys = mdlUpdate(t,x,u);
    case 3
        sys = mdlOutputs(t,x,u);
    otherwise
        sys = [];
end
end

function [sys,x0,str,ts,simStateCompliance] = mdlInitializeSizes
sizes = simsizes;
sizes.NumContStates  = 0;
sizes.NumDiscStates  = 1;   % SoC only
sizes.NumOutputs     = 1;
sizes.NumInputs      = 1;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;

sys = simsizes(sizes);
x0  = 1.0;  % initial SoC = 100%
str = [];
ts  = [1 0]; % 1 second sample time
simStateCompliance = 'DefaultSimState';
end

function sys = mdlUpdate(t,x,u)
% Constants
Ts    = 1;            % sample time [s]
Qbatt = 1.1 * 3600;   % battery capacity [C]
eta   = 1;          % coulombic efficiency

I  = u(1);            % input current (A)
SoC = x(1);

% Update rule: SoC(k+1) = SoC(k) - eta*I*Ts/Q
SoC_new = SoC - (eta * Ts * I) / Qbatt;

sys = SoC_new;
end

function sys = mdlOutputs(t,x,u)
sys = x(1); % output SoC
end
