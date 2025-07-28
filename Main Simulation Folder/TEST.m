function [sys,x0,str,ts,simStateCompliance] = TEST(t,x,u,flag)
% THETA2PHYS  Level-1 S-Function
%   Mengubah parameter regresi (theta1-theta6) hasil RLS ECM-2 menjadi
%   parameter fisik baterai: Voc, R0-R2, C1-C2.
%
% INPUT  u(1:6) : theta1 … theta6          (hasil estimator/RLS)
% OUTPUT sys    : [Voc  R0  R1  R2  C1  C2]
%
% 2025 © Harits_Estimator SoC
% -------------------------------------------------------------------------

switch flag
    %=====================================================================
    case 0      % Initialization
        sizes               = simsizes;
        sizes.NumContStates = 0;
        sizes.NumDiscStates = 0;
        sizes.NumOutputs    = 6;     % Voc R0 R1 R2 C1 C2
        sizes.NumInputs     = 6;     % theta1 … theta6
        sizes.DirFeedthrough= 1;     % langsung pakai input di mdlOutputs
        sizes.NumSampleTimes= 1;
        sys                 = simsizes(sizes);

        x0  = [];            % no states
        str = [];
        ts  = [1 0];         % sample-time 1 s (ubah sesuai kebutuhan)

        simStateCompliance = 'HasNoSimState';

    %=====================================================================
    case 1      % Derivatives  (tak ada state kontinyu)
        sys = [];

    %=====================================================================
    case 2      % Update (tak ada state diskrit)
        sys = [];

    %=====================================================================
    case 3      % Outputs
        theta = u(:).';              % pastikan baris-vektor
        if numel(theta) ~= 6
            error('Input width harus 6 (theta1…theta6).');
        end

        %----------- Ekstraksi koefisien regresi -------------------------
        a0 = theta(1);
        a1 = theta(2);
        a2 = theta(3);
        a3 = theta(4);
        a4 = theta(5);
        a5 = theta(6);

        %----------- Hitung Voc (discrete form, lihat persamaan 20) ------
        Voc = a0 / max(real(1 - a1 - a2), eps);

        %----------- Transformasi θ → parameter fisik --------------------
        % Persamaan (21)-(26) dalam literatur Zhang/Simão - diadaptasi
        % NB: asumsi orde-2 ECM: R0 + dua cabang RC

        % R0
        R0 = -(a3 - a4 + a5) / (a1 - a2 + 1);

        % Variabel bantu
        eq26 = (a3 + a4 + a5)/(a1 + a2 - 1) - R0;
        disc = (a2 + 1)^2 - (a1 + a2 - 1)*(a1 - a2 + 1);
        disc = max(disc, 0);         % hindari imag kecil akibat numerik

        arg1 = (-(a2 + 1) - sqrt(disc)) / (2*(a1 + a2 - 1));
        arg2 = (-(a2 + 1) + sqrt(disc)) / (2*(a1 + a2 - 1));

        eq21 = min(arg1, arg2);
        eq22 = max(arg1, arg2);

        eq23 = (a3 - a5)/(a1 + a2 - 1) - ...
               ((a3 - a4 + a5)/(a1 - a2 + 1))*((a2 + 1)/(a1 + a2 - 1));

        % Selesaikan dua persamaan linear (eq26 & eq23) untuk R1, R2
        B = [1 1; eq21 eq22] \ [eq26; eq23];   % R1, R2
        R1 = B(1);     R2 = B(2);

        % Kapasitansi
        C1 = eq21 / max(R1, eps);
        C2 = eq22 / max(R2, eps);

        %----------- Sanitasi (hindari NaN/Inf/negatif) ------------------
        param = [Voc R0 R1 R2 C1 C2];
        param(~isfinite(param) | param < 0) = 0;

        sys = param;

    %=====================================================================
    case 4      % Next time-hit (tak dipakai karena sample-time tetap)
        sys = t + 1;

    %=====================================================================
    case 9      % Terminate
        sys = [];

    %=====================================================================
    otherwise
        error('Unhandled flag = %d', flag);
end
