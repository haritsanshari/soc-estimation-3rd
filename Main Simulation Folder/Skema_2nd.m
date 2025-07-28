    function [sys,x0,str,ts,simStateCompliance] = Skema_2nd(t,x,u,flag)
    % Theta2ECM2  –  Konversi θ(1:6) → [R0 R1 C1 R2 C2 Voc]
    %
    % INPUT  u : [θ1 θ2 θ3 θ4 θ5 θ6]   (hasil RLS14)
    % OUTPUT y : [R0 R1 C1 R2 C2 Voc]
    %
    % Model ARX orde-2:
    %   Vt(k) = θ1 + θ2·Vt(k-1) + θ3·Vt(k-2) + θ4·I(k) + θ5·I(k-1) + θ6·I(k-2)
    %
    % Konversi mengikuti Zhang et al. (2018) dan turunan senior:
    %   • Voc  = θ1 / (1 − θ2 − θ3)
    %   • R0   = − (θ4 − θ5 + θ6)/(θ2 − θ3 + 1)
    %   • R1,R2,C1,C2 dihitung via sistem aljabar (τ = RC)
    % -------------------------------------------------------------------------
    
    switch flag
        case 0, [sys,x0,str,ts,simStateCompliance] = mdlInitSizes;
        case 3, sys = mdlOutputs(u);
        case {1,2,4,9}, sys = [];             % tak dipakai
        otherwise, error('Unhandled flag = %d',flag);
    end
    end
    
    %% =============== INITIALISASI ===========================================
    function [sys,x0,str,ts,simStateCompliance] = mdlInitSizes
    sizes = simsizes;
    sizes.NumContStates  = 0;
    sizes.NumDiscStates  = 0;
    sizes.NumOutputs     = 6;         % R0 R1 C1 R2 C2 Voc
    sizes.NumInputs      = 6;         % θ1 … θ6
    sizes.DirFeedthrough = 1;         % langsung gunakan u pada flag==3
    sizes.NumSampleTimes = 1;         % sample-time tetap
    
    sys = simsizes(sizes);
    
    x0  = [];                         % tak ada state
    str = [];
    ts  = [1 0];                      % Ts = 1 s (ubah jika perlu)
    
    simStateCompliance = 'HasNoSimState';
    end
    
    %% =================== OUTPUTS (flag==3) ==================================
    function y = mdlOutputs(theta)
    % ---- Ambil koefisien ----------------------------------------------------
    theta = double(theta(:));         % pastikan kolom
    t1 = theta(1); t2 = theta(2); t3 = theta(3);
    t4 = theta(4); t5 = theta(5); t6 = theta(6);
    
    % ---- Voc ---------------------------------------------------------------
    denVoc = 1 - t2 - t3;
    Voc = safeDiv(t1 , denVoc);
    
    % ---- Rumus bantu --------------------------------------------------------
    a1 = t2; a2 = t3; a3 = t4; a4 = t5; a5 = t6;
    
    R0 = safeDiv( -(a3 - a4 + a5) , (a1 - a2 + 1) );
    
    eq26 = safeDiv( (a3 + a4 + a5) , (a1 + a2 - 1) ) - R0;
    
    disc = (a2 + 1)^2 - (a1 + a2 - 1)*(a1 - a2 + 1);
    disc = max(disc, 0);                         % hindari sqrt<0
    
    tau1 = (-(a2 + 1) - sqrt(disc)) / (2*(a1 + a2 - 1));
    tau2 = (-(a2 + 1) + sqrt(disc)) / (2*(a1 + a2 - 1));
    
    % pastikan tau1 < tau2
    if tau1 > tau2, tmp = tau1; tau1 = tau2; tau2 = tmp; end
    
    eq23 = safeDiv( (a3 - a5) , (a1 + a2 - 1) ) - ...
           safeDiv( (a3 - a4 + a5) , (a1 - a2 + 1) ) * ...
           safeDiv( (a2 + 1) , (a1 + a2 - 1) );
    
    B = [1 1; tau1 tau2] \ [eq26; eq23];    % [R1; R2]
    R1 = B(1);  R2 = B(2);
    
    C1 = safeDiv(tau1 , R1);
    C2 = safeDiv(tau2 , R2);
    
    % ---- Sanitasi (jaga supaya positif & finite) ---------------------------
    R0 = clean(R0);  R1 = clean(R1);  R2 = clean(R2);
    C1 = clean(C1);  C2 = clean(C2);  Voc = clean(Voc);
    
    y  = [Voc R0 R1 C1 R2 C2]';
    end
    
    %% =================== UTILITAS ===========================================
    function z = safeDiv(num,den)
    % Bagi dengan perlindungan divide-by-zero
    epsd = 1e-12;
    z = num ./ (den + (abs(den)<epsd).*sign(den).*epsd);
    end
    
    function v = clean(v)
    % Ganti NaN, Inf, atau nilai negatif dengan nol
    v(~isfinite(v) | v<0) = 0;
    end
