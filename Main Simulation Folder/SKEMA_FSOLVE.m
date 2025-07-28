function [sys,x0,str,ts,simStateCompliance] = SKEMA_FSOLVE(t,x,u,flag)
% SKEMA_FSOLVE – S‑Function Level‑1 (MATLAB)
% -------------------------------------------------------------------------
% Konversi theta(1..8) ⇒ parameter fisik ECM orde‑3 menggunakan *fsolve*
% untuk 6 variabel:
%       R1, R2, R3   (Ω)
%       C1, C2, C3   (F)
%   Dengan 6 persamaan non‑linier:
%       • 3 persamaan tau (∑, Σ, Π)  —— tau_i = R_i*C_i
%       • 3 persamaan kombinasi R & tau —— ► lihat rumus di bawah.
% -------------------------------------------------------------------------
% INPUT  : u(1..8)  = theta1 .. theta8   (dari blok RLS)
% OUTPUT : [Voc R0 R1 R2 R3 C1 C2 C3 tau1 tau2 tau3]   (11×1)
% -------------------------------------------------------------------------
% Catatan implementasi
%   ▸ Discrete state  = 6 (menyimpan x = [R; C] last‑valid) agar reset otomatis
%   ▸ fsolve dipanggil di *mdlUpdate*; jika gagal → hold‑last‑valid
%   ▸ mdlOutputs hanya membaca state & menghitung Voc, R0, tau
% -------------------------------------------------------------------------

switch flag
    case 0
        [sys,x0,str,ts,simStateCompliance] = mdlInitializeSizes;
    case 2      % UPDATE : jalankan fsolve & simpan jika valid
        sys = mdlUpdate(t,x,u);
    case 3      % OUTPUT : keluarkan param + tau
        sys = mdlOutputs(t,x,u);
    otherwise
        sys = [];
end
end

%% ======================================================================
function [sys,x0,str,ts,simStateCompliance] = mdlInitializeSizes
    sizes = simsizes;
    sizes.NumContStates  = 0;
    sizes.NumDiscStates  = 6;          % R1 R2 R3 C1 C2 C3
    sizes.NumOutputs     = 11;
    sizes.NumInputs      = 8;
    sizes.DirFeedthrough = 1;
    sizes.NumSampleTimes = 1;

    sys = simsizes(sizes);

    % —— seed reasonable (editable) ——
    Rseed = [0.02; 0.01; 0.005];   % Ω
    Cseed = [1000; 400; 40];       % F  (tau ≈ 20‑40 s)
    x0    = [Rseed; Cseed];

    str = [];
    ts  = [1 0];                   % sample time = 1 s
    simStateCompliance = 'DefaultSimState';
end

%% ======================================================================
function sys = mdlOutputs(~,x,u)
    % State → R & C
    R = x(1:3);
    C = x(4:6);
    tau = R.*C;

    % theta & analitik Voc,R0
    th = u(:);
    th1=th(1); th2=th(2); th3=th(3); th4=th(4);
    th5=th(5); th6=th(6); th7=th(7); th8=th(8);

    D  = 1 - th2 - th3 - th4;    D  = safeguardDen(D);
    D0 = 1 + th2 - th3 + th4;    D0 = safeguardDen(D0);

    Voc = th1 / D;
    R0  = (-th5 + th6 - th7 + th8) / D0;

    sys = [ Voc; R0; R(:); C(:); tau(:) ];
end

%% ======================================================================
function newX = mdlUpdate(t,x,u)
    Rprev = x(1:3);
    Cprev = x(4:6);

    th = u(:);
    th1=th(1); th2=th(2); th3=th(3); th4=th(4);
    th5=th(5); th6=th(6); th7=th(7); th8=th(8);

    T = 1;
    D  = 1 - th2 - th3 - th4;    D  = safeguardDen(D);
    D0 = 1 + th2 - th3 + th4;    D0 = safeguardDen(D0);
    R0 = (-th5 + th6 - th7 + th8) / D0;

    % —— Right‑hand constants ————————————
    S1 = (T/2)  * (3 - th2 + th3 + 3*th4) / D;      % tau sum
    S2 = (T^2/4)* (3 + th2 + th3 - 3*th4) / D;      % tau pair‑sum
    S3 = (T^3/8)* (1 + th2 - th3 + th4) / D;        % tau product

    B1 = -(T^2/4)*(3*th5 - th6 - th7 + 3*th8 + R0*(3+th2+th3-3*th4)) / D;
    B2 = -(T/2)  *(3*th5 + th6 - th7 - 3*th8 + R0*(3-th2+th3+3*th4)) / D;
    B3 = (-(th5+th6+th7+th8)/D) - R0;

    % —— fsolve : unknowns [R1 R2 R3 C1 C2 C3] ———————————
    x0 = [Rprev; Cprev];
    opts = optimoptions('fsolve','Display','off','MaxIter',200,'FunctionTolerance',1e-12);
    [xSol,~,exitflag] = fsolve(@eqs,x0,opts);

    if exitflag <= 0  || any(xSol(1:3)<=0) || any(xSol(4:6)<=0)
        % gagal / solusi tak fisik → hold last valid
        newX = x;
        if abs(round(t)-t)<1e-9
            fprintf('t=%g  fsolve fail, hold states\n',t);
        end
        return;
    end
    newX = xSol(:);

    % — nested ————————————————————————————
    function F = eqs(z)
        R = z(1:3);
        C = z(4:6);
        tau = R.*C;
        tau1=tau(1); tau2=tau(2); tau3=tau(3);
        R1=R(1); R2=R(2); R3=R(3);

        F = [ tau1*tau2*tau3                           - S3;
              (tau1*tau2 + tau1*tau3 + tau2*tau3)      - S2;
              (tau1 + tau2 + tau3)                     - S1;
              (R1*tau2*tau3 + R2*tau1*tau3 + R3*tau1*tau2) - B1;
              (R1*(tau2+tau3)+ R2*(tau1+tau3)+ R3*(tau1+tau2)) - B2;
              (R1+R2+R3)                               - B3 ];
    end
end

%% ======================================================================
function D = safeguardDen(D)
    if abs(D) < 1e-10
        D = sign(D)*1e-10;
    end
end
