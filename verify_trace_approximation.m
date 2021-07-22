clear
clc
Mt_adc = 128;

R = U + Z;
Cz = 1;
Z = Cz*randn(Mt_adc);


% MC = 100;

% sigma_u_range = [1/6 1/5 1/4 1/3 1/2 1 2];
% H = eye(Mt_adc);
% 
% constellationVector  = [-1 1];
% constellationSize = 2;            
% Es = 1;
% Eb = Es / log2(constellationSize);
% sigma = sqrt(Eb);
% 
% bit_range = [1 16];
% 
% meanMI = zeros(length(sigma_u_range),length(bit_range));
% meanTFC = zeros(length(sigma_u_range),length(bit_range));
% inputSNR = zeros(length(sigma_u_range),length(bit_range));
% for r=1:MC
% 
%   for r=1:length(bit_range)
%     bits = bit_range(r);
%     mi = zeros(length(sigma_u_range),1);
%     trFC = zeros(length(sigma_u_range),1);
%     trCovX = zeros(length(sigma_u_range),1);
%     u =  randsrc(Mt_adc, 1, constellationVector);
%     x = zeros(Mt_adc, length(sigma_u_range));
%     
%     for b=1:length(sigma_u_range)
%       sigma_u = sigma_u_range(b)*sigma;
%       sigma_e = 0;%sqrt(1-pi*sqrt(3)/2 * 2.^(-2*bits));
%       Delta = eye(Mt_adc); %sqrt(1-pi*sqrt(3)/2 * 2.^(-2*bits))*
%       x(:, b) = Delta*u + sigma_e*randn(Mt_adc,1);
%       Cov_x = sigma_u^2*Delta*Delta' + sigma_e^2*eye(Mt_adc);
%       trCovX(b) = trace(Cov_x)/Mt_adc;
% 
% 
%       inputSNR(b, r) = 10*log10(sigma_u^2/(trCovX(b)-sigma_u^2));
% 
%       y = H*x(:,b);
%       mi(b) = mutualinfo(y,x(:, b));
%       FIM = H*Cov_x*H';
%       trFC(b) = trace(FIM*Cov_x)/trace(x(:,b)*x(:,b)')^2;
% 
%     end
% 
%     meanMI(:, r) = meanMI(:, r) + mi;
%     meanTFC(:, r) = meanTFC(:, r) + trFC;
%     
%   end
% end
% figure;
% subplot(1,2,1)
% plot(inputSNR(:, 1), meanMI(:, 1)/MC, 'o-', inputSNR(:, 1), meanTFC(:, 1)/MC, 'x-');
% xlabel('input SNR (dB)')
% ylabel('Mutual Information');
% grid on;
% title('1 bit DACs')
% subplot(1,2,2)
% plot(inputSNR(:, 2), meanMI(:, 2)/MC, 'o-', inputSNR(:, 2), meanTFC(:, 2)/MC, 'x-');
% xlabel('input SNR (dB)')
% ylabel('Mutual Information');
% grid on;
% title('3 bit DACs')
% legend('Mutual Information (theoretical)', 'trace(J C_x) (approximation)')
% 
