function [sys,x0,str,ts,simStateCompliance] = SKEMA_A(t,x,u,flag)
% SKEMA_A – S‑Function Level‑1 (MATLAB)
% ---------------------------------------------------------------
%  ▸ Mengonversi parameter regresi theta1..theta8 (orde‑3 ECM)
%    menjadi : Voc, R0, R1‑R3, C1‑C3, tau1‑tau3.
%  ▸ Menyimpan tau & R terakhir yang valid di D‑Work (discrete states),
%    sehingga *reset* otomatis setiap kali simulasi di‑Run ulang.
%  ▸ Bila:
%        – akar polinom tidak 3 buah real‑positif,   **atau**
%        – matriks A sangat ill‑conditioned / singular,
%    maka blok akan *hold‑last‑valid* (tidak update state).
% ---------------------------------------------------------------
%  INPUT  (8×1): theta1..theta8
%  OUTPUT (11×1): [Voc R0 R1 R2 R3 C1 C2 C3 tau1 tau2 tau3]
% ---------------------------------------------------------------

switch flag
    case 0   % =========== INITIALISATION ======================
        [sys,x0,str,ts,simStateCompliance] = mdlInitializeSizes;

    case 2   % =========== UPDATE DISCRETE STATES ==============
        sys = mdlUpdate(t,x,u);

    case 3   % =========== OUTPUT ==============================
        sys = mdlOutputs(t,x,u);

    otherwise
        sys = [];
end
end

%% =============================================================
function [sys,x0,str,ts,simStateCompliance] = mdlInitializeSizes
    sizes = simsizes;
    sizes.NumContStates  = 0;
    sizes.NumDiscStates  = 6;   % [tau1 tau2 tau3 R1 R2 R3]
    sizes.NumOutputs     = 11;  % Voc R0 R1 R2 R3 C1 C2 C3 tau1 tau2 tau3
    sizes.NumInputs      = 8;   % theta1..theta8
    sizes.DirFeedthrough = 1;
    sizes.NumSampleTimes = 1;

    sys = simsizes(sizes);

    % ------- seed awal (dipakai pertama kali) -----------------
    tauSeed = [30; 5; 1];        % detik (boleh diubah)
    RSeed   = [0.02; 0.01; 0.005];

    x0  = [tauSeed; RSeed];      % 6×1
    str = [];
    ts  = [1 0];                 % sample‑time diskret 1 s
    simStateCompliance = 'DefaultSimState';
end

%% =============================================================
function sys = mdlOutputs(~,x,u)
    % ----- ambil state terakhir valid -------------------------
    tauPrev = x(1:3);
    Rprev   = x(4:6);

    % ----- de‑serialize theta ---------------------------------
    th = u(:);
    th1=th(1); th2=th(2); th3=th(3); th4=th(4);
    th5=th(5); th6=th(6); th7=th(7); th8=th(8);
    T = 1;

    %% 1) Voc dan R0 ------------------------------------------
    D  = 1 - th2 - th3 - th4;   D  = safeguardDen(D);
    D0 = 1 + th2 - th3 + th4;   D0 = safeguardDen(D0);
    Voc = th1 / D;
    R0  = (-th5 + th6 - th7 + th8) / D0;

    %% 2) Hitung tau (akar polinom) ----------------------------
    sigma3 = (T^3/8)*(1 + th2 - th3 + th4) / D;
    sigma2 = (T^2/4)*(3 + th2 + th3 - 3*th4) / D;
    sigma1 = (T/2)  *(3 - th2 + th3 + 3*th4) / D;
    r      = roots([1 -sigma1 sigma2 -sigma3]);
    tau    = chooseTau(r);

    %% 3) Hitung R1‑R3 ----------------------------------------
    [R, okR] = solveR(tau, R0, th, T, D);

    % ---- validasi keseluruhan -------------------------------
    if isempty(tau) || ~okR
        tau = tauPrev;   % HOLD
        R   = Rprev;
    end

    % ---- Kapasitansi ----------------------------------------
    C = zeros(3,1);
    tolR = 1e-12;
    for k = 1:3
        if abs(R(k)) > tolR
            C(k) = tau(k) / R(k);
        end
    end

    %% 4) keluaran --------------------------------------------
    sys = [ Voc; R0; R(:); C(:); tau(:) ];
end

%% =============================================================
function newX = mdlUpdate(~,x,u)
    % Re‑hitung tau & R memakai fungsi sama dengan mdlOutputs.
    %  → jika valid, simpan; jika tidak, pertahankan state lama.

    tauPrev = x(1:3);
    Rprev   = x(4:6);

    % theta
    th = u(:);
    th1=th(1); th2=th(2); th3=th(3); th4=th(4);
    th5=th(5); th6=th(6); th7=th(7); th8=th(8);
    T=1;

    D  = 1 - th2 - th3 - th4;    D  = safeguardDen(D);
    D0 = 1 + th2 - th3 + th4;    D0 = safeguardDen(D0);
    R0 = (-th5 + th6 - th7 + th8) / D0;

    sigma3 = (T^3/8)*(1 + th2 - th3 + th4) / D;
    sigma2 = (T^2/4)*(3 + th2 + th3 - 3*th4) / D;
    sigma1 = (T/2)  *(3 - th2 + th3 + 3*th4) / D;

    r   = roots([1 -sigma1 sigma2 -sigma3]);
    tau = chooseTau(r);
    [R, okR] = solveR(tau, R0, th, T, D);

    if isempty(tau) || ~okR
        tau = tauPrev;
        R   = Rprev;
    end

    newX = [tau(:); R(:)];
end

%% ==================== helper ==========================================
function D = safeguardDen(D)
    if abs(D) < 1e-8
        D = sign(D)*1e-8;
    end
end

function tau = chooseTau(r)
    % Ambil akar real‑positif; butuh tiga buah.
    tauReal = r( abs(imag(r))<1e-6 & real(r)>0 );
    if numel(tauReal) < 3
        tau = [];
        return;
    end
    taureal = real(tauReal(:));   % ascending, ambil tiga terkecil
    tau = taureal(1:3);
end

function [R, ok] = solveR(tau, R0, th, T, D)
    ok = false;
    if isempty(tau)
        R = zeros(3,1); return; end

    tau1=tau(1); tau2=tau(2); tau3=tau(3);
    th2=th(2); th3=th(3); th4=th(4);
    th5=th(5); th6=th(6); th7=th(7); th8=th(8);

    B1 = -(T^2/4)*(3*th5-th6-th7+3*th8 + R0*(3+th2+th3-3*th4))/D;
    B2 = -(T/2)*( 3*th5+th6-th7-3*th8 + R0*(3-th2+th3+3*th4))/D;
    B3 = (- (th5+th6+th7+th8)/D) - R0;
    B  = [B1; B2; B3];

    A  = [ tau2*tau3,   tau1*tau3,   tau1*tau2;
           tau2+tau3,   tau1+tau3,   tau1+tau2;
           1,           1,           1         ];

    if cond(A) > 1e12
        R = zeros(3,1); return; end

    R = real(A\B);
    if any(R<=0)
        R = zeros(3,1); return; end

    ok = true;
end
