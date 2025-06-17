%% Monodomain FEM solver – final corrected version
%   Scientific Machine Learning – Final Project
%   Author:  Leon Ackermann
%   Updated: 2025-06-08
% -------------------------------------------------------------------------
clear; clc; close all;

%% 1. Physical parameters (from project hand-out)
Sigma_h = 9.5298e-4;              % healthy conductivity
Sigma_d = 10 * Sigma_h;          % diseased conductivity

a  = 18.515;  ft = 0.2383;  fr = 0;  fd = 1;
f_u = @(u) a.*(u-fr).*(u-ft).*(u-fd);   % cubic reaction term

%% 2. Numerical parameters
T   = 200;             % total simulation time  [s]
dt  = 0.1;          % time step
nvx = 64;  nvy = 64;          % vertices in x and y
hx  = 1/(nvx-1);  hy = 1/(nvy-1);

%% 3. Structured mesh -----------------------------------------------------
[X,Y] = meshgrid(linspace(0,1,nvx), linspace(0,1,nvy));
nodes = [X(:)  Y(:)];                    % nv × 2
nv    = size(nodes,1);

id = reshape(1:nv, nvx, nvy);            % node-index grid

% ---- element connectivity (counter-clockwise square) -------------------
n1 = id(1:end-1, 1:end-1);   n1 = n1(:);     % lower-left
n2 = id(2:end  , 1:end-1);   n2 = n2(:);     % lower-right
n3 = id(2:end  , 2:end  );   n3 = n3(:);     % upper-right
n4 = id(1:end-1, 2:end  );   n4 = n4(:);     % upper-left

conn = [n1'; n2'; n3'; n4'];               % 4 × ne
ne   = numel(n1);

% ---- element centres (vectorised) --------------------------------------
xc = (nodes(n1,1) + nodes(n2,1) + nodes(n3,1) + nodes(n4,1))/4;
yc = (nodes(n1,2) + nodes(n2,2) + nodes(n3,2) + nodes(n4,2))/4;

%% 4. Heterogeneous conductivity -----------------------------------------
r2      = [0.12^2  0.152^2  0.12^2];           % r² values from PDF
centres = [0.30 0.70 ;                   % circle 1
           0.70 0.30 ;                   % circle 2
           0.50 0.50 ];                  % circle 3

is_diseased = false(ne,1);
for k = 1:3
    dx = xc - centres(k,1);
    dy = yc - centres(k,2);
    is_diseased = is_diseased | (dx.^2 + dy.^2 <= r2(k));
end

fprintf('Diseased elements detected: %d / %d\n', nnz(is_diseased), ne);
assert(nnz(is_diseased) > 0, 'No diseased elements detected – check mesh.');

sigma_e              = Sigma_h * ones(ne,1);
sigma_e(is_diseased) = Sigma_d;

%% 5. Conductivity sanity plot -------------------------------------------
figure('Name','Conductivity map');
sigma_grid = reshape(sigma_e, nvx-1, nvy-1)';        % row-major view


%% 6. Assemble FEM matrices ----------------------------------------------
fprintf('Assembling M and K_sigma …\n');
M       = assembleMass(nvx, nvy, hx, hy);               % sparse nv × nv
K_sigma = assembleDiffusion_modified(nvx, nvy, hx, hy, sigma_e);
A       = M + dt*K_sigma;                               % left-hand side

%% 7. Initial condition ---------------------------------------------------
U                   = zeros(nv,1);
U(nodes(:,1)>=0.9 & nodes(:,2)>=0.9) = 1;   % corner stimulus

%% 8. Time stepping (IMEX Euler) -----------------------------------------
nSteps    = round(T/dt);
plotEvery = 1;
figure('Name','Monodomain simulation');

for n = 1:nSteps
    U = A \ (M*(U - dt*f_u(U)));          % explicit reaction, implicit diff.
    
    if mod(n,plotEvery)==0 || n==nSteps
        Ugrid = reshape(U, nvy, nvx);
        contourf(X, Y, Ugrid, linspace(-.1,1.1,12), 'LineColor','none');
        hold on
        %contour(X, Y, Ugrid, [ft ft], 'w', 'LineWidth', 2);
        %viscircles(centres, sqrt(r2') , 'LineStyle','--', 'Color','r');
        hold off
        axis equal tight, caxis([-.1 1.1]), colorbar
        title(sprintf('t = %.2f  s  (%d / %d)', n*dt, n, nSteps));
        drawnow
    end
end
disp('Simulation complete.');