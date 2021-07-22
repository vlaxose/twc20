function [H,Ar,At] = parametric_mmwave_channel(Mr, Mt, total_num_of_clusters, total_num_of_rays)
% Implementation of "Simplified Spatial Correlation Models for Clustered
% MIMO Channels With Different Array Configurations", R. Heath, 2007

% Assumptions:
% - all the clusters have equal power (sigma_phi)

% Initialization
H = zeros(Mr, Mt);
Hl = zeros(Mr, Mt);
Ar = zeros(Mr, total_num_of_clusters*total_num_of_rays);
At = zeros(Mt, total_num_of_clusters*total_num_of_rays);
% Main
index = 1;
for tap = 1:total_num_of_clusters
  
    for ray=1:total_num_of_rays
        rayleigh_coeff = 1/sqrt(2)*(randn(1)+1j*randn(1));
        Ar(:, index) = angle(genLaplacianSamples(1), Mr);
        At(:, index) = angle(genLaplacianSamples(1), Mt);
        
        Hl = Hl + rayleigh_coeff*Ar(:, index)*At(:, index)';

        index = index + 1;
    end
    H = H + Hl;
end

    H = sqrt(Mr*Mt)/sqrt(total_num_of_clusters*total_num_of_rays)*H;
end

% Generate the transmit and receive array responces
function vectors_of_angles=angle(phi, M)

    % For Uniform Linear Arrays (ULA) compute the phase shift
    Ghz = 30;
    wavelength = 30/Ghz; % w=c/lambda
    array_element_spacing = 0.5*wavelength;
    wavenumber = 2*pi/wavelength; % k = 2pi/lambda
    phi0 = 0; % mean AOA
    phase_shift = wavenumber*array_element_spacing*sin(phi0-phi)*(0:M-1).';
    vectors_of_angles = 1/sqrt(M)*exp(1j*phase_shift);
end


% Random variable generator based on inverse transform sampling
function x=genLaplacianSamples(N)
    u = rand(N,1);
    % Inverse transform of trancuted Laplacian
    sigma_phi = 50; % standard deviation of the power azimuth spectrum (PAS)
    beta = 1/(1-exp(-sqrt(2)*pi/sigma_phi));
    x = beta*(exp(-sqrt(2)/sigma_phi*pi) - cosh(u));
end