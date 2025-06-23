function [t_act,isM,inRange,U,nodes,X,Y] = monodomain_solver( ...
           dt,nc,Sigma_fac,T_final,plotEvery,saveSteps)
%MONODOMAIN_SOLVER  2-D monodomain IMEX FEM on a structured square mesh
%
%   [...] = MONODOMAIN_SOLVER(dt,nc,Sigma_fac,T_final,plotEvery,...
%                             saveSteps,savePrefix)
%
%   New optional inputs
%     saveSteps  – vector of time-step indices whose solution should be
%                  exported as PNG (default [], i.e. none)
%     savePrefix – file-name stem for exported frames (default 'frame')
%
%   Rest of interface unchanged – see original header for details.
% -------------------------------------------------------------------------
if nargin<4, T_final   = 35;       end
if nargin<5, plotEvery = 0;        end
if nargin<6, saveSteps = [];       end
if nargin<7, savePrefix = 'frame'; end

%% 1. Physical parameters -------------------------------------------------
Sigma_h = 9.5298e-4;
Sigma_d = Sigma_fac*Sigma_h;

a  = 18.515;  f_t = 0.2383;  f_r = 0;  f_d = 1;
f_u = @(u) a.*(u-f_r).*(u-f_t).*(u-f_d);       % cubic reaction

%% 2. Mesh ---------------------------------------------------------------
nvx = nc+1;  nvy = nvx;
hx  = 1/(nvx-1);  hy = hx;

[X,Y] = meshgrid(linspace(0,1,nvx), linspace(0,1,nvy));
nodes = [X(:) Y(:)];                 nv  = size(nodes,1);
id    = reshape(1:nv,nvx,nvy);

n1 = id(1:end-1,1:end-1); n1 = n1(:);
n2 = id(2:end  ,1:end-1); n2 = n2(:);
n3 = id(2:end  ,2:end  ); n3 = n3(:);
n4 = id(1:end-1,2:end  ); n4 = n4(:);
conn = [n1';n2';n3';n4'];            ne  = numel(n1);

xc = mean(reshape(nodes(conn,1),4,ne))';
yc = mean(reshape(nodes(conn,2),4,ne))';

%% 3. Conductivity map ---------------------------------------------------
r2  = [0.10^2 0.15^2 0.10^2];
ctr = [0.30 0.70; 0.70 0.30; 0.50 0.50];
is_d = false(ne,1);
for k = 1:3
    dx = xc-ctr(k,1);  dy = yc-ctr(k,2);
    is_d = is_d | (dx.^2+dy.^2 <= r2(k));
end
sigma_e           = Sigma_h*ones(ne,1);
sigma_e(is_d)     = Sigma_d;

%% 4. FEM matrices -------------------------------------------------------
M = assembleMass(nvx,nvy,hx,hy);
K = assembleDiffusion_modified(nvx,nvy,hx,hy,sigma_e);
A = M + dt*K;

offDiag = A - spdiags(diag(A),0,nv,nv);
[~,~,v] = find(offDiag);
isM = all(diag(A)>0) && all(v<=1e-12);

%% 5. Initial condition --------------------------------------------------
U = zeros(nv,1);
U(nodes(:,1)>=0.9 & nodes(:,2)>=0.9) = 1;
stim = (nodes(:,1)>=0.9 & nodes(:,2)>=0.9);  % mask, if needed later

%% 6. Time loop ----------------------------------------------------------
nSteps  = round(T_final/dt);
t_act   = NaN;   active = false;

if plotEvery>0 || ~isempty(saveSteps)
    hFig = figure('Name','Monodomain evolution');
end

for n = 1:nSteps
    RHS = M*(U - dt*f_u(U));
    U   = A \ RHS;

    if ~active && all(U(~stim) > f_t)
        t_act = n*dt;   active = true;
    end

    wantRegularPlot =  plotEvery>0 && (mod(n,plotEvery)==0 || n==nSteps);
    wantSavedFrame  =  ismember(n,saveSteps);

    if wantRegularPlot || wantSavedFrame
        figure(hFig); cla
        Ugrid = reshape(U,nvy,nvx);
        contourf(X,Y,Ugrid,linspace(-.1,1.1,12),'LineColor','none');
        axis equal tight; caxis([-.1 1.1]); colorbar
        title(sprintf('t = %.2f  s  (%d / %d)',n*dt,n,nSteps));
        drawnow

        if wantSavedFrame
            fname = sprintf('step_%05d.png', n);
            exportgraphics(gca, fname, 'Resolution', 300, 'BackgroundColor','none');        end
    end
end

%% 7. Range check --------------------------------------------------------
epsU    = 1e-6;
inRange = (min(U)>=-epsU) && (max(U)<=1+epsU);
end