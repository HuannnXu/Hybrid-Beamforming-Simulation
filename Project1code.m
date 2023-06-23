%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project 1
%Huan Xu
%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem 1
% parameter setting for result 1
clear all;
N = 64; %BS has N antennas
M = 16;%users equip with M antennas
NRF = 6; % and N^{RF}_t transmit RF chains
K = 1;% serve K users
% and N^{RF}_r transmit RF chains
d = 6;%each user requires d data streams
Ns = K * d;
V_D = zeros(NRF, Ns);
V_RF = ones(N, NRF);
W_Dk = ones(NRF, d);
W_RFk = ones(M, NRF);
% y_tilda = zeros(1, Ns);
L = 15;
iter = 100;
SNR = -10 :2: 6; % dB
P = 0.1;
SNRlin = 10.^(SNR / 10);
gamma2 = P / (N * NRF);
R_avg = zeros(size(SNR));
MonteCarloIterNum = 100;
%% iterate
for snr = 1 : length(SNRlin)
    SNR(snr)
    sigma2 = P / SNRlin(snr);    
    for iter = 1 : MonteCarloIterNum
    %% channel model
        alpha_list = randn(1, L) + 1j * randn(1, L);
        phi_r_list = unifrnd(0, pi * 2, [1, L]);
        phi_t_list = unifrnd(0, pi * 2, [1, L]);
        H = zeros(M, N);
        for l = 1 : L
            a_phi_r = 1 / sqrt(N) * exp(1j * pi * sin(phi_r_list(l)) .* [0 : N - 1])';
            a_phi_t = 1 / sqrt(M) * exp(1j * pi * sin(phi_t_list(l)) .* [0 : M - 1])';
            H = H + alpha_list(l) * a_phi_t * a_phi_r';
        end
        H = sqrt(N * M / L) * H;
        %% RF Precoder Design
        F1 = H' * H;
        C = zeros(NRF - 1, NRF - 1,NRF);
        G = zeros(N, N, NRF);
        flag = 0;
        while flag == 0
            flag = 1;
            for j = 1 : NRF
                VRFj_ba = V_RF(:, [1 : j - 1, j + 1 : end]);
                C(:, :, j) = eye(NRF - 1, NRF - 1) + (gamma2 / sigma2) .* (VRFj_ba' * F1 * VRFj_ba);
                G(:, :, j) = gamma2 / sigma2 * F1 - (gamma2 / sigma2) ^ 2 .* (F1 * VRFj_ba * pinv(C(:, :, j)) * VRFj_ba' * F1);
                for i = 1 : N
                    yeta = 0;
                    for l = 1 : N
                        if l ~= i && i <= N
                            yeta = yeta + G(i, l, j) * V_RF(l, j);
                        end
                        if yeta == 0
                            V_RF(i, j) = 1;
                        else
                            V_RF(i, j) = yeta / abs(yeta);                    
                        end
                        if isnan(V_RF(i, j))
                            flag = 0;
                        end
                    end
                end
            end
        end
        %% Digital Design of precoders - get V_D
        N0 = P / SNRlin(snr);
        H_eff = H * V_RF;
        Q = V_RF' * V_RF;
        Q_invsqrt = inv(sqrtm(Q));
        sizeQ = size(Q);
        [~, ~, U_e] = svd(H_eff * Q_invsqrt);
        Gamma_e = eye(sizeQ(1)) * sqrt(P / NRF);
        V_D = Q_invsqrt * U_e * Gamma_e;
        %% Find W_RFk
        V_t = V_RF * V_D;
        F2 = H * V_t * V_t' * H';
        C = zeros(NRF - 1, NRF - 1, NRF);
        G = zeros(M, M, NRF);
        flag = 0;
        while flag == 0
            flag = 1;
            for j = 1 : NRF
                WRFj_ba = W_RFk(:, [1 : j - 1, j + 1 : end]);
                C(:, :, j) = eye(NRF - 1, NRF - 1) + (1 / (sigma2 * M)) .* (WRFj_ba' * F2 * WRFj_ba);
                G(:, :, j) = 1 / (sigma2 * M) .* F2 - (1 / (sigma2 * M)) ^ 2 .* (F2 * WRFj_ba * inv(C(:, :, j)) * WRFj_ba' * F2);      
                for i = 1 : M
                    yeta = 0;
                    for l = 1 : M
                        if l ~= i && i <= M
                            yeta = yeta + G(i, l, j) * W_RFk(l, j);
                        end
                        if yeta == 0
                            W_RFk(i, j) = 1;
                        else
                            W_RFk(i, j) = yeta / abs(yeta);
                        end
                    end
                end
            end
        end
        %% Design of hybrid combiners get W_Dk
        J = W_RFk' * H * V_t * V_t' * H' * W_RFk + sigma2 * W_RFk' * W_RFk;
        W_Dk = inv(J) * W_RFk' * H * V_t;
        %% Results
        W_t = W_RFk * W_Dk;
        R = log2(det(eye(M) + 1 / sigma2 * ...
                              W_t * inv(W_t' * W_t) * W_t' * H * V_t * V_t' * H'));
        R_avg(snr) = R_avg(snr) + R;
    end
    R_avg(snr) = R_avg(snr) / MonteCarloIterNum;
end
%% Plot
figure;
plot(SNR, R_avg,'-o', "LineWidth", 2);
title("Spectral efficiencies VS SNR, 64x16 MIMO, N^{RF}=N_S=6");
xlabel("SNR(dB)")
ylabel("Spectral Efficiency(bits/s/Hz)")
grid on;
