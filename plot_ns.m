clear
clc

ns_range = [1 5 10 15 20 25];
rho = 1;
noise_variance = sqrt(10^(-5/10));
Nt=128;Lt=128;
maxDBiters = 5;
min_bit=1;
max_bit=8;
K = 12;
mean_rate = zeros(length(ns_range), K);
mean_ee = zeros(length(ns_range), K);
mean_Lt_opt = zeros(length(ns_range), K);
for ns_index = 1:length(ns_range)
    Ns = ns_range(ns_index);
    Nr = Ns;
   disp(['ns: ', num2str(ns_range(ns_index ))])
   [mean_rate(ns_index, :), mean_ee(ns_index, :)] = monte_carlo_sims(rho, noise_variance, Nt, Nr, Lt, Ns, K, maxDBiters, min_bit, max_bit);
end

W = 1;
figure;
p=plot(ns_range, W*mean_ee(:, 1)); hold on;
set(p,'LineWidth',2, 'LineStyle', '-', 'Color', 'Black');
p=plot(ns_range, W*mean_ee(:, 2));hold on;
set(p,'LineWidth',0.4, 'LineStyle', '--', 'Color', 'Black');
p=plot(ns_range, W*mean_ee(:, 3));hold on;
set(p,'LineWidth',0.6, 'LineStyle', '--', 'Color', 'Black');
p=plot(ns_range, W*mean_ee(:, 4));hold on;
set(p,'LineWidth',0.8, 'LineStyle', '--', 'Color', 'Black');
p=plot(ns_range, W*mean_ee(:, 5));hold on;
set(p,'LineWidth',1, 'LineStyle', '-', 'Color', 'Black');
p=plot(ns_range, W*mean_ee(:, 6));hold on;
set(p,'LineWidth',1.2, 'LineStyle', ':', 'Color', 'Black');
p=plot(ns_range, W*mean_ee(:, 7));hold on;
set(p,'LineWidth',1.4, 'LineStyle', '-', 'Color', 'Black');
p=plot(ns_range, W*mean_ee(:, 8));hold on;
set(p,'LineWidth',1.6, 'LineStyle', '--', 'Color', 'Black');
p=plot(ns_range, W*mean_ee(:, 9));hold on;
set(p,'LineWidth',1.8, 'LineStyle', '-.', 'Color', 'Black');
p=plot(ns_range, W*mean_ee(:, 10));hold on;
set(p,'LineWidth',2, 'LineStyle', ':', 'MarkerEdgeColor', 'Black', 'MarkerFaceColor', 'Black', 'Marker', 'o', 'MarkerSize', 6, 'Color', 'Black');
p=plot(ns_range, W*mean_ee(:, 11));hold on;
set(p,'LineWidth',2, 'LineStyle', '-', 'MarkerEdgeColor', 'Green', 'MarkerFaceColor', 'Green', 'Marker', 's', 'MarkerSize', 6, 'Color', 'Green');
p=plot(ns_range, W*mean_ee(:, 12));hold on;
set(p,'LineWidth',2, 'LineStyle', '--', 'MarkerEdgeColor', 'Blue', 'MarkerFaceColor', 'Blue', 'Marker', '>', 'MarkerSize', 6, 'Color', 'Blue');
% p=plot(ns_range, W*mean_ee(:, 13));hold on;
% set(p,'LineWidth',2, 'LineStyle', '--', 'MarkerEdgeColor', 'Blue', 'MarkerFaceColor', 'Blue', 'Marker', '<', 'MarkerSize', 6, 'Color', 'Blue');
grid on;
xlabel('Number of transmission streams', 'FontSize', 11)
ylabel({'Energy Efficiency', '(Mbits/Joule)'}, 'FontSize', 11)
% lg = legend('DBF', '1-bit-HBF', '2-bit-HBF', '3-bit-HBF', '4-bit-HBF', '5-bit-HBF', '6-bit-HBF', '7-bit-HBF', '8-bit-HBF', 'random-bits-HBF', 'Proposed', 'common-exhaustive-HBF', 'common-iterative-HBF', 'Location', 'Best');
lg = legend('DBF', '1-bit-HBF', '2-bit-HBF', '3-bit-HBF', '4-bit-HBF', '5-bit-HBF', '6-bit-HBF', '7-bit-HBF', '8-bit-HBF', 'random-bits-HBF', 'Proposed', 'common-exhaustive-HBF', 'Location', 'Best');
lg.FontSize = 8;

savefig(['./results/ee_ns_',num2str(min_bit),'_',num2str(max_bit),'.fig'])
print(['./results/ee_ns_',num2str(min_bit),'_',num2str(max_bit),'.eps'],'-depsc')


W = 1;
figure;
p=plot(ns_range, W*mean_rate(:, 1)); hold on;
set(p,'LineWidth',2, 'LineStyle', '-', 'Color', 'Black');
p=plot(ns_range, W*mean_rate(:, 2));hold on;
set(p,'LineWidth',0.4, 'LineStyle', '--', 'Color', 'Black');
p=plot(ns_range, W*mean_rate(:, 3));hold on;
set(p,'LineWidth',0.6, 'LineStyle', '--', 'Color', 'Black');
p=plot(ns_range, W*mean_rate(:, 4));hold on;
set(p,'LineWidth',0.8, 'LineStyle', '--', 'Color', 'Black');
p=plot(ns_range, W*mean_rate(:, 5));hold on;
set(p,'LineWidth',1, 'LineStyle', '-', 'Color', 'Black');
p=plot(ns_range, W*mean_rate(:, 6));hold on;
set(p,'LineWidth',1.2, 'LineStyle', ':', 'Color', 'Black');
p=plot(ns_range, W*mean_rate(:, 7));hold on;
set(p,'LineWidth',1.4, 'LineStyle', '-', 'Color', 'Black');
p=plot(ns_range, W*mean_rate(:, 8));hold on;
set(p,'LineWidth',1.6, 'LineStyle', '--', 'Color', 'Black');
p=plot(ns_range, W*mean_rate(:, 9));hold on;
set(p,'LineWidth',1.8, 'LineStyle', '-.', 'Color', 'Black');
p=plot(ns_range, W*mean_rate(:, 10));hold on;
set(p,'LineWidth',2, 'LineStyle', ':', 'MarkerEdgeColor', 'Black', 'MarkerFaceColor', 'Black', 'Marker', 'o', 'MarkerSize', 6, 'Color', 'Black');
p=plot(ns_range, W*mean_rate(:, 11));hold on;
set(p,'LineWidth',2, 'LineStyle', '-', 'MarkerEdgeColor', 'Green', 'MarkerFaceColor', 'Green', 'Marker', 's', 'MarkerSize', 6, 'Color', 'Green');
p=plot(ns_range, W*mean_rate(:, 12));hold on;
set(p,'LineWidth',2, 'LineStyle', '--', 'MarkerEdgeColor', 'Blue', 'MarkerFaceColor', 'Blue', 'Marker', '>', 'MarkerSize', 6, 'Color', 'Blue');
% p=plot(ns_range, W*mean_ee(:, 13));hold on;
% set(p,'LineWidth',2, 'LineStyle', '--', 'MarkerEdgeColor', 'Blue', 'MarkerFaceColor', 'Blue', 'Marker', '<', 'MarkerSize', 6, 'Color', 'Blue');
grid on;
xlabel('Number of transmission streams', 'FontSize', 11)
ylabel({'Spectral Efficiency', '(bits/sec)'}, 'FontSize', 11)
% lg = legend('DBF', '1-bit-HBF', '2-bit-HBF', '3-bit-HBF', '4-bit-HBF', '5-bit-HBF', '6-bit-HBF', '7-bit-HBF', '8-bit-HBF', 'random-bits-HBF', 'Proposed', 'common-exhaustive-HBF', 'common-iterative-HBF', 'Location', 'Best');
lg = legend('DBF', '1-bit-HBF', '2-bit-HBF', '3-bit-HBF', '4-bit-HBF', '5-bit-HBF', '6-bit-HBF', '7-bit-HBF', '8-bit-HBF', 'random-bits-HBF', 'Proposed', 'common-exhaustive-HBF', 'Location', 'Best');
lg.FontSize = 8;

savefig(['./results/rate_ns_',num2str(min_bit),'_',num2str(max_bit),'.fig'])
print(['./results/rate_ns_',num2str(min_bit),'_',num2str(max_bit),'.eps'],'-depsc')

