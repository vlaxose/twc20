clear
clc

rho_range = ([0.01 0.1 1 10 100 1000]);
noise_variance = sqrt(10^(-5/10));
Nt=128;Lt=128;Ns=8;Nr=Ns;
maxDBiters = 5;
min_bit=1;
max_bit=8;
K = 12;
mean_bitalloc = zeros(length(rho_range), 2);
min_bitalloc = zeros(length(rho_range), 2);
max_bitalloc = zeros(length(rho_range), 2);
for snr_index = 1:length(rho_range)
   disp(['pu: ', num2str(rho_range(snr_index ))])
   [~,~,~,~,~,~,~, mean_bitalloc(snr_index,:), min_bitalloc(snr_index,:), max_bitalloc(snr_index,:)] = monte_carlo_sims(rho_range(snr_index), noise_variance, Nt, Nr, Lt, Ns, K, maxDBiters, min_bit, max_bit);
end

figure;subplot(1,2,1)
p=plot(10*log10(rho_range), mean_bitalloc(:, 1)); hold on;
set(p,'LineWidth',2, 'LineStyle', '-', 'MarkerEdgeColor', 'Green', 'MarkerFaceColor', 'Green', 'Marker', 's', 'MarkerSize', 8, 'Color', 'Green');
p=plot(10*log10(rho_range), min_bitalloc(:,1)); hold on;
set(p,'LineWidth',2, 'LineStyle', '--', 'Color', 'Blue');
p=plot(10*log10(rho_range), max_bitalloc(:,1)); hold on;
set(p,'LineWidth',2, 'LineStyle', '--', 'Color', 'Red');
grid on;
title('Proposed')
ylabel({'Sum of all DACs resolutions'}, 'FontSize', 12)
xlabel('Transmit Power (dBm)', 'FontSize', 12)
% lg = legend('Mean values', 'Minimum values', 'Maximum values', 'Location', 'Best');
lg.FontSize = 8;
subplot(1,2,2)
p=plot(10*log10(rho_range), mean_bitalloc(:, 2)); hold on;
set(p,'LineWidth',2, 'LineStyle', '-', 'MarkerEdgeColor', 'Green', 'MarkerFaceColor', 'Green', 'Marker', 'o', 'MarkerSize', 8, 'Color', 'Green');
p=plot(10*log10(rho_range), min_bitalloc(:,2)); hold on;
set(p,'LineWidth',2, 'LineStyle', '--', 'Color', 'Blue');
p=plot(10*log10(rho_range), max_bitalloc(:,2)); hold on;
set(p,'LineWidth',2, 'LineStyle', '--', 'Color', 'Red');
grid on;
title('common-exhaustive-HBF')
ylabel({'Sum of all DACs resolutions'}, 'FontSize', 12)
xlabel('Transmit Power (dBm)', 'FontSize', 12)
lg.FontSize = 8;
lg = legend('Mean values', 'Minimum values', 'Maximum values', 'Location', 'Best', 'Orientation', 'Horizontal');
savefig(['./results/rate_bitalloc.fig'])
print(['./results/rate_bitalloc.eps'],'-depsc')
