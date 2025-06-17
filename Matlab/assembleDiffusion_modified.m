function A = assembleDiffusion_modified(nvx, nvy, hx, hy, sigma_per_element)
% ASSEMBLEDIFFUSION_MODIFIED Assembles the global stiffness matrix K for
% a heterogeneous conductivity Sigma.
%
%   This version CORRECTLY uses the sigma_per_element input to scale
%   the local stiffness matrix for each element.

% 1D reference diffusion matrix
Aref = [1 -1; -1 1];
% 1D reference mass matrix
Mref = [1/3 1/6; 1/6 1/3];

% 1D diffusion matrix along x
Ax = 1/hx * Aref;
% 1D diffusion matrix along y
Ay = 1/hy * Aref;

% 1D mass matrix along x
Mx = hx * Mref;
% 1D mass matrix along y
My = hy * Mref;

% Calculate the LOCAL STIFFNESS MATRIX for a REFERENCE element with sigma = 1.
Aloc_base = kron(My, Ax) + kron(Ay, Mx);

% --- Connectivity (using the provided functions' convention) ---
nv = nvx * nvy;
ne = (nvx-1) * (nvy-1);

if length(sigma_per_element) ~= ne
    error('The length of sigma_per_element must match the number of elements.');
end

id = 1:nv;
id = reshape(id, nvx, nvy);
a = id(1:end-1, 1:end-1); a = a(:)';
b = id(2:end, 1:end-1);   b = b(:)';
c = id(1:end-1, 2:end);   c = c(:)';
d = id(2:end, 2:end);     d = d(:)';
conn = [a; b; c; d];

% --- Assembly using sparse ---
ii = (1:4)';
ii = repmat(ii, [1 4]);
jj = ii';
I = conn(ii(:), :);
J = conn(jj(:), :);

% =======================================================================
% --- THIS IS THE CRITICAL FIX ---
% Instead of replicating a single Aloc, we scale the base local matrix
% by the sigma value for each element. The outer product between the
% 16x1 Aloc_base vector and the 1xne sigma_per_element vector creates
% the 16-by-ne values matrix V in one efficient operation.
V = Aloc_base(:) * sigma_per_element';
% =======================================================================

A = sparse(I(:), J(:), V(:), nv, nv);
end