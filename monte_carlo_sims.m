function [mean_rate, mean_ee, mean_power, mean_Lt_opt, min_Lt_opt, max_Lt_opt, mean_kappa, mean_bitalloc, min_bitalloc, max_bitalloc] = monte_carlo_sims(rho, noise_variance, Nt, Nr, Lt, Ns, K, maxDBiters, min_bit, max_bit)

  total_num_of_clusters = 3;
  total_num_of_rays = 4;
  
  total_monte_carlo_realizations = 1000;
  total_monte_carlo_realizations_converged = 0;
  mean_rate = zeros(1, K);
  mean_ee = zeros(1, K);
  mean_power = zeros(1, K);
  mean_bitalloc = zeros(1, 2);  
  mean_kappa = zeros(1, maxDBiters+1);
  mean_Lt_opt = zeros(1, K);

  for r=1:total_monte_carlo_realizations

    [rate, ee, power, Lt_opt, kappa, bit_alloc, converged] = systemModel(Nt, Nr, Lt, Ns, total_num_of_clusters, total_num_of_rays, rho, noise_variance, maxDBiters, 0, min_bit, max_bit);

    if(converged)
      mean_rate = mean_rate+rate;
      mean_ee = mean_ee + ee;
      mean_power = mean_power + power;
      for k=1:K
          if r==1
            min_Lt_opt(1,k) = Lt_opt(k);
            max_Lt_opt(1,k) = Lt_opt(k);
          else
            min_Lt_opt(1,k) = min(min_Lt_opt(1,k), Lt_opt(k));
            max_Lt_opt(1,k) = max(max_Lt_opt(1,k), Lt_opt(k));
          end
      end
      mean_Lt_opt = mean_Lt_opt + Lt_opt;
      for k=1:2
          if r==1
            min_bitalloc(1,k) = sum(bit_alloc(:, k));
            max_bitalloc(1,k) = sum(bit_alloc(:, k));
          else
            min_bitalloc(1,k) = min(min_bitalloc(1,k), sum(bit_alloc(:, k)));
            max_bitalloc(1,k) = max(max_bitalloc(1,k), sum(bit_alloc(:, k)));
          end
      end
      mean_bitalloc = mean_bitalloc + sum(bit_alloc, 1);      
      mean_kappa = mean_kappa + kappa;
      total_monte_carlo_realizations_converged = total_monte_carlo_realizations_converged+1;
    end
  end
  
   mean_rate = mean_rate/total_monte_carlo_realizations_converged;
   mean_ee = mean_ee/total_monte_carlo_realizations_converged;
   mean_power = mean_power/total_monte_carlo_realizations_converged;
   mean_Lt_opt = mean_Lt_opt/total_monte_carlo_realizations_converged;
   mean_bitalloc = mean_bitalloc/total_monte_carlo_realizations_converged;
   mean_kappa = mean_kappa/total_monte_carlo_realizations_converged;
  
end