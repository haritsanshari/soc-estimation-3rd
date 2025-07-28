function [sys,x0,str,ts] = AHIF(t,x,u,flag)
% AHINF  Level-1 S-Function: 
%   Adaptive H∞ Filter untuk ECM Orde-3 dengan update kovarian secara real time
%
%   INPUT (u):
%     u(1) = I        (arus baterai, A)
%     u(2) = Vt_meas  (tegangan terminal terukur, V)
%     u(3) = R0       (ohmic resistance, Ω)
%     u(4) = R1
%     u(5) = R2
%     u(6) = R3
%     u(7) = C1
%     u(8) = C2
%     u(9) = C3
%
%   OUTPUT (sys):
%     sys(1) = SoC_hat    (perkiraan state-of-charge)
%     sys(2) = Vt_hat     (tegangan terminal terprediksi)
%     sys(3) = V1_hat
%     sys(4) = V2_hat
%     sys(5) = V3_hat
%
%   Adaptasi H∞ Filter sesuai Pers. (24)–(34) di Tabel II [24]:
%     - M̂k = (1/N) ∑ ei ei^T ;  R̂k = M̂k − Hk Pk− Hk^T         (30) 
%     - Kk = Fk Pk− [ I − θ S̄k Pk− + Hk^T Rk^{-1} Hk Pk− ]^{-1} Hk^T Rk^{-1}  (31) 
%     - Q̂k = Kk M̂k Kk^T                                           (32) 
%     - x̂k^+ = x̂k^- + Kk ek                                      (33)
%     - Pk^+ = Pk^- [ I − θ S̄k Pk^- + Hk^T Rk^{-1} Hk Pk^- ]^{-1}   (34)
%
%   BUFFER SIZE untuk inovasi:
   N = 20;                   % panjang jendela geser untuk M̂k (sesuai Tabel II, eks: 120) 

   Ts      = 1;               % Sampling time (s)
   Qbatt   = 1.1 * 3600;      % Kapasitas baterai (Coulomb)
   eta     = 1;               % Efisiensi coulombic
   theta   = 0.01;            % Parameter H∞
   % --- Tidak ada Qproc, Rmeas statis: akan diupdate adaptif di discrete states ---

switch flag
    case 0  % Initialization
        sizes = simsizes;
        sizes.NumContStates  = 0;
        % ----------------------------------------------
        %    Jumlah Discrete States:
        %      4  = xhat (SoC, V1, V2, V3)
        %     16  = P(:)  (kovarian 4×4 disusun kolom-per-kolom)
        %     N   = e2_buffer(1..N) untuk menyimpan e^2 lama
        %      1  = idx_buffer (penunjuk posisi dalam e2_buffer)
        %     16  = Q_flat(:) (kovarian proses 4×4 dari iterasi sebelumnya)
        %      1  = R_meas  (kovarian pengukuran scalar dari iterasi sebelumnya)
        % ----------------------------------------------
        sizes.NumDiscStates  = 4 + 16 + N + 1 + 16 + 1;  
        sizes.NumOutputs     = 5;
        sizes.NumInputs      = 9;
        sizes.DirFeedthrough = 1;
        sizes.NumSampleTimes = 1;
        sys = simsizes(sizes);

        % ----- Inisialisasi x̂₀⁺ dan P₀⁺ -----
        xhat0 = [0.9; 0; 0; 0];                    % Perkiraan awal: SoC=0.9, V1-3=0
        P0    = diag([1e-4 1e-2 1e-2 1e-2]);       % Kovarian awal (P₀⁺)
        % ----- Inisialisasi buffer inovasi (semua nol) -----
        e2_buff0 = zeros(N,1);                     % Belum ada inovasi sebelumnya
        idx0     = 0;                               % Indeks awal
        % ----- Inisialisasi Qproc dan Rmeas -----
        Qproc0 = diag([1e-8 1e-5 1e-5 1e-5]);       % Nilai awal Q₀ (4×4)
        Rmeas0 = 1e-4;                              % Nilai awal R₀ (scalar)

        % Susun vektor x0 (size = 4 +16 + N +1 +16 +1)
        x0 = [ ...
            xhat0;              ... %  1:4
            P0(:);               ... %  5:20
            e2_buff0;            ... % 21:20+N
            idx0;                ... % 21+N
            Qproc0(:);           ... % 22+N : 37+N
            Rmeas0               ... % 38+N
        ];
        str = [];
        ts  = [Ts 0];

    case 2  % Update block (Time + Measurement + Adaptive Update)
        % -----------------------------------
        % Ambil input:
        I    = u(1);
        Vt   = u(2);
        R0   = u(3);
        R1   = u(4);
        R2   = u(5);
        R3   = u(6);
        C1   = u(7);
        C2   = u(8);
        C3   = u(9);

        % -----------------------------------
        % Ambil state lama:
        xhat       = x(1:4);                             % [SoC; V1; V2; V3]
        P          = reshape(x(5:20), 4, 4);              % kovarian P_k-
        e2_buffer  = x(21:20+N);                          % buffer inovasi e^2
        idx_buffer = x(21+N);                             % pointer lama (1..N)
        Qproc      = reshape(x(22+N:37+N), 4, 4);         % Q dari iterasi sebelumnya
        Rmeas      = x(38+N);                             % R dari iterasi sebelumnya (scalar)

        % -----------------------------------
        % 1) **Time Update (Prediksi)** sesuai Pers. (25)–(27)
        %  - Hitung matriks diskretisasi A untuk tiap RC:
        A1 = exp(-Ts/(R1*C1));  
        A2 = exp(-Ts/(R2*C2));  
        A3 = exp(-Ts/(R3*C3));
        Fk = diag([1, A1, A2, A3]);    % F_{k-1} 

        %  - Fungsi non-linier f(xhat,u):
        f = [ ...
            xhat(1) - eta*Ts*I/Qbatt;           ... % SoC⁻
            A1*xhat(2) + R1*(1-A1)*I;           ... % V1⁻
            A2*xhat(3) + R2*(1-A2)*I;           ... % V2⁻
            A3*xhat(4) + R3*(1-A3)*I            ... % V3⁻
        ];

        %  - Prediksi kovarian: Pm = Fk P_k-1 Fk' + Qproc_{k-1}  (Pers. (26)) 
        Pm = Fk * P * Fk' + Qproc;

        %  - Update S̄k = L^T S L ; di sini L = I, S = I ⇒ S̄k = I  (Pers. (27))
        Sk_bar = eye(4);

        % -----------------------------------
        % 2) **Measurement Update (Linearization & Inovasi)**

        % 2a) Hitung Jacobian output Hk = [dVoc/dSoC, -1, -1, -1]
        dVoc = Voc_grad(f(1));
        Hk   = [dVoc, -1, -1, -1];   % (1×4)

        % 2b) Prediksi keluaran: ŷ = Voc(SoC⁻) - I*R0 - Σ Vi⁻
        Voc_est = open_circuit_voltage(f(1));
        yhat    = Voc_est - I*R0 - f(2) - f(3) - f(4);

        % 2c) Hitung inovasi ek = y_meas - ŷ
        ek = Vt - yhat;               % (skalar)

        % -----------------------------------
        % 3) **Adaptive Update untuk Rk** (Pers. (30))
        %   M̂k = (1/N) ∑_{i=1..N} e_i^2
        %   Disini kita simpan e^2 dalam buffer ring.
        idx_new = mod(idx_buffer, N) + 1;         % indeks baru (1..N)
        e2_buffer(idx_new) = ek^2;                % gantikan elemen tertua
        M_hat = sum(e2_buffer) / N;               % rata-rata e^2  (skalar) 

        %   R̂k = M̂k - Hk Pm Hk^T    (Pers. (30))  
        R_hat = M_hat - Hk * Pm * Hk';
        % Pastikan R_hat > 0 (bila perlu clamp):
        if R_hat < 1e-8
            R_hat = 1e-8;
        end

        % -----------------------------------
        % 4) **Kalman Gain Adaptif** (Pers. (31))
        %   inv_term = [ I - θ S̄k Pm + Hk^T Rk^{-1} Hk Pm ]
        inv_term = eye(4) ...
                   - theta * (Sk_bar * Pm) ...
                   + (Hk' * (1/R_hat) * Hk) * Pm;   % Hk'*(1/R_hat)*Hk = Hk^T R^{-1} Hk

        %   Kk = Fk Pm [inv(inv_term)] Hk^T Rk^{-1}
        K = Fk * Pm * (inv(inv_term)) * Hk' * (1/R_hat);  % (4×1) 

        % -----------------------------------
        % 5) **Adaptive Update untuk Qk** (Pers. (32))
        %   Q̂k = Kk M̂k Kk^T    (4×4) 
        Q_hat = K * M_hat * K';  % (4×4)

        % -----------------------------------
        % 6) **Measurement Update** (Pers. (33) & (34))
        %   x̂k⁺ = f + K ek       (4×1) :contentReference[oaicite:1]{index=1}
        xhat_new = f + K * ek;

        %   Pk⁺ = Pm [ inv(inv_term ) ]   (4×4) :contentReference[oaicite:2]{index=2}
        P_new = Pm * inv(inv_term);

        % Clamp SoC di rentang [0,1]
        xhat_new(1) = min(max(xhat_new(1), 0), 1);

        % -----------------------------------
        % Paketkan state diskret baru:
        sys = [ ...
            xhat_new;                   ... % 1:4
            P_new(:);                   ... % 5:20
            e2_buffer;                  ... % 21:20+N
            idx_new;                    ... % 21+N
            Q_hat(:);                   ... % 22+N : 37+N
            R_hat                       ... % 38+N
        ];

    case 3  % Output block
        % Output = [SoC; Vt_hat; V1; V2; V3]
        xhat    = x(1:4);
        I       = u(1);
        R0      = u(3);
        Voc_est = open_circuit_voltage(xhat(1));
        Vt_hat  = Voc_est - I*R0 - xhat(2) - xhat(3) - xhat(4);
        sys     = [ xhat(1); Vt_hat; xhat(2:4) ];

    case {1,4,9}
        sys = [];      % tidak dipakai untuk flag ini

    otherwise
        error(['Unhandled flag = ', num2str(flag)]);
end
end

% --------- Voc Polynomial & Gradient ----------
function V = open_circuit_voltage(SoC)
    % Koefisien polinomial Voc(SoC) (sesuai literatur Zhang et al.)
    p = [-1.0979e3,  4.9095e3,  -9.0725e3,  8.9523e3,  ...
         -5.0861e3,   1.6719e3,  -3.0330e2,  2.7200e1,  2.3000];
    V = polyval(p, SoC);
end

function dV = Voc_grad(SoC)
    % Turunan polinomial Voc'(SoC)
    p  = [-1.0979e3,  4.9095e3,  -9.0725e3,  8.9523e3,  ...
          -5.0861e3,  1.6719e3,  -3.0330e2, 2.7200e1,  2.3000];
    dp = polyder(p);
    dV = polyval(dp, SoC);
end
