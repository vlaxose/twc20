function [F,F_BB,F_RF,F_DBF,W] = beamformer(H, At, Lt, Ns, type)

    %% Initialization
    [Nr, Nt] = size(H);
    F_RF  = eye(Nt, Lt);
    F_BB = zeros(Lt, Ns);
    
    %% Obtaion the optimal case, i.e., digital beamformer
    [U,D,V] = svd(H);
    F_DBF = V(:, 1:Ns);

    %%
    switch(type)
        case 'angular_codebook'
            D = At;
        case 'fft_codebook'
            D = 1/sqrt(Nt)*fft(eye(Nt));
    end

    
    %%
    % TX
    if(Nt==Lt)
        F_RF = eye(Lt);
        F_BB = F_DBF;
        F = F_RF*F_BB;
    else
        for indx=1:Lt
          % Get the index with the maximum energy
          diff = F_DBF - F_RF*F_BB;
          Psi = D'*diff/norm(diff, 'fro');
          C = diag(Psi*Psi');
          [~, I] = sort(C, 'descend');

          % Update the precoders
          F_RF(:, indx) = D(:,I(1));
          F_BB = F_RF\F_DBF;
        end
        
        F = 1/norm(F_BB)*F_RF*F_BB;
        
        if(norm(F_DBF - F)^2/norm(F_DBF)^2>1e-1)
            warning(['A/D beamformer did not converge, with error:', num2str(norm(F_DBF - F, 'fro'))]);
        end
    end
    
    % RX
    W =  U(:, 1:Ns);
    W = W/ norm(W, 'fro');

        

end

