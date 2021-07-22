function [rate, ee, power, Lt_opt, kappa, bit_alloc, converged] = systemModel(Nt, Nr, Lt, Ns, total_num_of_clusters, total_num_of_rays, rho, noise_variance, maxDBiters, check_rate_approximation_flag, min_bit, max_bit)

    %% Setup
    K = 12;
    rate = zeros(1, K);
    ee = zeros(1, K);
    power = zeros(1, K);
    Lt_opt = zeros(1, K);
    bit_alloc = zeros(Lt, 2);
    p_dac = 0.01;
    Pcp = 1000;
    p_PA = 0.01;
    Pps = 10*10^(-3);
    Pt = 100*10^(-3);
    codebook = fft(eye(Nt));
    F_RF_DFT = 1/sqrt(Nt)*codebook(:, 1:Lt);
    
    %% mmWave channel model
    [H, ~, At] = parametric_mmwave_channel(Nr, Nt, total_num_of_clusters, total_num_of_rays);
    
    %% Optimal beamforming
    [~,F_BB,F_RF,F_DBF,W] = beamformer(H, At, Lt, Ns, 'fft_codebook');
    
    %% Digital beamforming with high resolution DACs (8bits) and no power constraint
    rate(:, 1) = real(log2(det(eye(Ns) + rho*1/Ns*inv(noise_variance^2*(W'*W))*W'*H*(F_DBF*F_DBF')*H'*W)));
    power(:, 1) = ((2*Nt*2^8*p_dac) + Nt*Pt + Nt^2*Pps)/Pcp + rho*p_PA*trace(F_DBF*F_DBF');
    ee(:, 1) = rate(:, 1) / power(:, 1);
    Lt_opt(:, 1) = Nt;

    % Hybrid beamforming
    
    % Random bit allocation
    % 1-8 bits cases
    indx = 2;
    for bit=1:8
        Lt_opt(:, indx) = Ns;
        bit_vec = bit*ones(Lt_opt(:, 2),1);
        Delta = diag(sqrt(1-pi*sqrt(3)/2 * 2.^(-2*bit_vec)));
        [ee(:,indx), rate(:,indx), power(:,indx)] = ee_computation(eye(Lt_opt(:, indx)), noise_variance, bit_vec, Lt_opt(:, indx), Ns, W, H, F_RF(:, 1:Lt_opt(:,2)), F_BB(1:Lt_opt(:,2), 1:Lt_opt(:,2)), Delta, p_dac, Pt, Pps, Pcp, rho);
        indx = indx + 1;
    end
    
    % Random bit allocation
    random_bits_vec = randi([min_bit max_bit], Ns, 1);
    Lt_opt(:, 10) = Ns;
    Delta = diag(sqrt(1-pi*sqrt(3)/2 * 2.^(-2*random_bits_vec)));
    [ee(:,10), rate(:,10), power(:, 10)] = ee_computation(eye(Lt_opt(:, 10)), noise_variance, random_bits_vec, Lt_opt(:, 10), Ns, W, H, F_RF(:, 1:Lt_opt(:,2)), F_BB(1:Lt_opt(:,2), 1:Lt_opt(:,2)), Delta, p_dac, Pt, Pps, Pcp, rho);
         
    % Energy-efficient optimally selected streams (proposed)
    random_bits_vec = randi([min_bit max_bit], Lt, 1);
    Delta = diag(sqrt(1-pi*sqrt(3)/2 * 2.^(-2*random_bits_vec)));
    [S_proposed, Lt_opt(:, 11), kappa, converged, ~] = proposed(H,W,F_RF_DFT,F_BB,Ns,Lt,maxDBiters, random_bits_vec, noise_variance, p_dac, Pt, Pcp, Pps, check_rate_approximation_flag, rho);
    [ee(:,11), rate(:,11), power(:,11)] = ee_computation(S_proposed, noise_variance, random_bits_vec, Lt_opt(:, 11), Ns, W, H, F_RF_DFT, F_BB, Delta, p_dac, Pt, Pps, Pcp, rho);
    bit_alloc(:, 1) = S_proposed*random_bits_vec;
    
    % Exhaustive search with common bit resolution for all RF chains
    if(Ns<=12)
        binary_combinations = dec2bin(0:2^Ns - 1) - '0';
        for lt=1:size(binary_combinations, 1)
            S_bin_comb = zeros(Lt);
            S_bin_comb(1:Ns,1:Ns) = diag(binary_combinations(lt, :));
            for bit = min_bit:max_bit
               bit_vec = bit*ones(Lt,1);
               Delta = diag(sqrt(1-pi*sqrt(3)/2 * 2.^(-2*bit_vec)));
               [ee_bit_bf, rate_bit_bf, power_bit_bf] = ee_computation(S_bin_comb, noise_variance, bit_vec, nnz(diag(S_bin_comb)), Ns, W, H, F_RF, F_BB, Delta, p_dac, Pt, Pps, Pcp, rho);
               if ee_bit_bf>ee(:,12)
                 ee(:,12) = ee_bit_bf;
                 power(:, 12) = power_bit_bf;
                 rate(:, 12) = rate_bit_bf;
                 Lt_opt(:, 12) = nnz(diag(S_bin_comb));
                 bit_bf = bit;
                 bit_alloc(:, 2) = S_bin_comb*bit_vec;
                 S_exhaustive = S_bin_comb;
               end
           end
        end
    end

%     % Active RF minimization, with common bit resolution for all RF chains
%     S_iterative = zeros(Lt, 1);
%     for lt=1:Lt
%         F_BB_lt = F_BB(1:lt, :);
%         F_RF_lt = F_RF_DFT(:, 1:lt); 
%         for bit = min_bit:max_bit
%            bit_vec = bit*ones(lt,1);
%            Delta = diag(sqrt(1-pi*sqrt(3)/2 * 2.^(-2*bit_vec)));
%            [ee_bit_it, rate_bit_it, power_bit_it] = ee_computation(eye(lt), noise_variance, bit_vec, lt, Ns, W, H, F_RF_lt, F_BB_lt, Delta, p_dac, Pt, Pps, Pcp, rho);
%            if ee_bit_it>ee(:,13)
%              ee(:,13) = ee_bit_it;
%              rate(:, 13) = rate_bit_it;
%              power(:, 13) = power_bit_it;
%              Lt_opt(:, 13) = lt;
%              bit_bf = bit;
%              S_iterative(1:lt, :) = ones(lt, 1);
%            end
%        end
%     end
    

end
