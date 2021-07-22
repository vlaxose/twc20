clear
clc

Nt_range = [48 64 128];
noise_variance = sqrt(10^(-5/10));
rho = 2;
Ns=1;
Nr=8;
Lr=Nr;
min_bit=8;
max_bit=8;
maxDBiters = 5;
kappa = zeros(length(Nt_range), maxDBiters+1);
K=13;

for nt_index = 1:length(Nt_range)
   Nt = Nt_range(nt_index);
   Lt = 0.5*Nt;
   disp(['Lt: ', num2str(Nt_range(nt_index ))])
   [~,~,~,kappa(nt_index,:)] = monte_carlo_sims(rho, noise_variance, Nt, Nr, Lt, Ns, K, maxDBiters, min_bit, max_bit);
end


figure;
p=plot(1:maxDBiters+1, kappa(1, :)); hold on;
set(p,'LineWidth',1, 'LineStyle', '-', 'MarkerEdgeColor', 'Black', 'MarkerFaceColor', 'White', 'Marker', 'o', 'MarkerSize', 8, 'Color', 'Black');
p=plot(1:maxDBiters+1, kappa(2, :));hold on;
set(p,'LineWidth',1, 'LineStyle', '-', 'MarkerEdgeColor', 'Black', 'MarkerFaceColor', 'White', 'Marker', 's', 'MarkerSize', 8, 'Color', 'Black');
p=plot(1:maxDBiters+1, kappa(3, :));hold on;
set(p,'LineWidth',1, 'LineStyle', '-', 'MarkerEdgeColor', 'Black', 'MarkerFaceColor', 'White', 'Marker', 'h', 'MarkerSize', 8, 'Color', 'Black');
xlabel('Dinkelbach iteration', 'FontSize', 12)
ylabel('log(1+ \mu tr(A))', 'FontSize', 12)
grid on;

g = legend(['N_T=',num2str(Nt_range(1))], ['N_T=',num2str(Nt_range(2))], ['N_T=',num2str(Nt_range(3))], 'Location', 'Best');
lg.FontSize = 11;
savefig('./results/DB_convergence.fig')
print('./results/DB_convergence.eps','-depsc')
