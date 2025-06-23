%% Batch driver for exercise 7  – parameter sweep using MONODOMAIN_SOLVER
clear; clc;

dt_vec    = [0.10 0.05 0.025];
nCells    = [64 128 256];        % elements per side
Sigma_fac = [10 1 0.1];          % Σ_d / Σ_h
T_final   = 35;                  % [s]
save_steps = [];

results = [];                                             % accumulator
for dt = dt_vec
    for nc = nCells
        for fac = Sigma_fac
            [t_act,isM,inRange] = monodomain_solver(dt,nc,fac,T_final,0, save_steps);
            results = [results; dt nc fac t_act isM inRange];
            fprintf('dt=%6.3f  nc=%3d  Σ_d/Σ_h=%4.1f  t_act=%7.3f  M=%d  in[0,1]=%d\n',...
                    dt,nc,fac,t_act,isM,inRange);
        end
    end
end

T = array2table(results,...
    'VariableNames',{'dt','n_e','Sigma_fac','t_act','M_matrix','u_in_[0,1]'});
disp(' '); disp(T);