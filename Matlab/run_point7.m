%% Batch driver for item 7  – Monodomain FEM
%  Runs 3×3×3 parameter sweep  (Δt, mesh, Σ_d) and prints diagnostics
%  Author:  Leon Ackermann   – 2025-06-08
% -----------------------------------------------------------------------
clear; clc;

%% Physical parameters
Sigma_h = 9.5298e-4;
Sigma_fac = [10, 1, 0.1];                % Σ_d / Σ_h
a = 18.515; ft = 0.2383; fr = 0; fd = 1;
f_u = @(u) a.*(u-fr).*(u-ft).*(u-fd);

T_final = 35;

%% Sweeps
dt_vec  = [0.10 0.05 0.025];
nCells  = [64 128 256];                  % elements per side

res = [];                                % results table

for dt = dt_vec
  for nc = nCells
    nvx = nc+1; nvy = nvx;               % vertices
    hx  = 1/(nvx-1);  hy = hx;

    % ---- structured mesh ----------------------------------------------
    [X,Y]  = meshgrid(linspace(0,1,nvx), linspace(0,1,nvy));
    nodes  = [X(:) Y(:)];
    nv     = size(nodes,1);
    id     = reshape(1:nv,nvx,nvy);

    n1 = id(1:end-1,1:end-1); n1 = n1(:);
    n2 = id(2:end  ,1:end-1); n2 = n2(:);
    n3 = id(2:end  ,2:end  ); n3 = n3(:);
    n4 = id(1:end-1,2:end  ); n4 = n4(:);
    conn = [n1'; n2'; n3'; n4'];                % 4×ne
    ne   = numel(n1);

    xc = mean(reshape(nodes(conn,1),4,ne))';
    yc = mean(reshape(nodes(conn,2),4,ne))';

    % disease mask
    r2 = [0.12^2 0.152^2 0.12^2];
    ctr= [0.30 0.70; 0.70 0.30; 0.50 0.50];
    is_d = false(ne,1);
    for k=1:3, dx = xc-ctr(k,1); dy = yc-ctr(k,2);
        is_d = is_d | (dx.^2+dy.^2 <= r2(k));
    end

    % ---- loop over Σ_d values -----------------------------------------
    for fac = Sigma_fac
      Sigma_d = fac*Sigma_h;
      sigma_e           = Sigma_h*ones(ne,1);
      sigma_e(is_d)     = Sigma_d;

      M = assembleMass(nvx,nvy,hx,hy);
      K = assembleDiffusion_modified(nvx,nvy,hx,hy,sigma_e);
      A = M + dt*K;                       % constant system matrix

      % safe M-matrix test
      offDiag = A - spdiags(diag(A),0,nv,nv);
      [~,~,v] = find(offDiag);
      isM = all(diag(A) > 0) && all(v <= 1e-14);

      % initial state
      U = zeros(nv,1);  U(nodes(:,1)>=.9 & nodes(:,2)>=.9) = 1;

      nSteps = round(T_final/dt);  active=false; t_act=NaN;
      for n=1:nSteps
          U = A \ (M*(U - dt*f_u(U)));
          if ~active && any(U>ft), t_act = n*dt; active=true; end
      end
      inRange = min(U)>=0 && max(U)<=1;

      res = [res; dt nc Sigma_d t_act isM inRange];
      fprintf('dt=%6.3f  nc=%3d  Σd=%8.1e  t_act=%7.3f  M=%d  in[0,1]=%d\n',...
              dt,nc,Sigma_d,t_act,isM,inRange);
    end
  end
end

%% pretty table
T = array2table(res, ...
     'VariableNames',{'dt','nCells','Sigma_d','t_act','M_matrix','u_in_[0,1]'});
disp(' '), disp(T)