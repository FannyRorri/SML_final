%% Visual highlight of three conductivity scenarios (exercise 8)
clear; clc;

dt      = 0.1;          % moderate time step
nc      = 64;           % 128×128 elements
facList = [0.1];    % Σ_d / Σ_h to illustrate
T_final = 35;
plot    = 1;
save_steps = [0,50,100,150,200,250,300,350];

for fac = facList
    figure('Name',sprintf('Final u – Σ_d/Σ_h = %.1f',fac)); clf
    [~,~,~,U,nodes,X,Y] = monodomain_solver(dt,nc,fac,T_final,plot,save_steps);
    Ugrid = reshape(U,nc+1,nc+1);
    contourf(X,Y,Ugrid,linspace(-.1,1.1,20),'LineColor','none');
    axis equal tight; caxis([-.1 1.1]); colorbar
    title(sprintf('Final potential, Σ_d/Σ_h = %.1f',fac));
end