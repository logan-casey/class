% main_solve.m


load('hw_solution_tightgrid_rtildeonly_feb15_s10.mat')
opts = struct();
% % initialize
opts.rtilde_init = eq.rtilde;
opts.wbar_init = eq.wbar;
% Warm starts for V/policy
opts.V_init = eq.V;
opts.pol_init = struct('pol_ik', eq.pol_ik, 'pol_ib', eq.pol_ib);
clearvars -except opts

par = hw.params("full");

% Shock discretization for solve
[zgrid, Pz] = hw.tauchen(par.rho, par.sigma_eps, par.Nz_solve, par.tauchen_m);

% Capital/debt grids
par.k_exponents = [0,0.5,1:15];
grid = hw.grids(par, zgrid);

% bmin = grid.bgrid(1);
% bmax = grid.bgrid(end);
bmin = -600;
bmax = -bmin;

% coarse negative grid
bneg = linspace(bmin, 0, 11)';

% dense positive grid up to bmax (concentrate near 0 with power >1)
Nbpos = 81;
x = linspace(0,1,Nbpos)'; 
bpos = (x.^1.5) * bmax;   % 1.5 concentrates more points near 0

grid.bgrid = unique([bneg; bpos]);  % unique removes duplicate 0
grid.Nb = numel(grid.bgrid);


% Revised net worth grid (debugging starter)
Nw = 120;
wgrid = linspace(-2*grid.kbar, 2*grid.kbar, Nw);

% Solve equilibrium
opts.outer_maxit   = 900;
opts.tol_rtilde    = 1e-3;
opts.outer_verbose = true;

% Damping and caps for rtilde
opts.omega_rtilde  = 0.02;
opts.cap_rtilde    = 0.5;
opts.minPrND       = 1e-4;

% Damping and caps for wbar
opts.damp_wbar = true;
opts.eta_wbar = 0.5;
opts.wbar_s = 10;

% Inner (Howard) controls
opts.inner_maxit    = 5;
opts.inner_tol      = 1e-6;
opts.howard_iters   = 20;
opts.improve_every  = 10;
opts.inner_verbose  = true;


eq = hw.solve_rtilde(wgrid, zgrid, Pz, par, grid, opts);

save('hw_solution_tightgrid_rtildeonly_feb15_s10.mat.mat', 'eq', 'par', 'grid', 'wgrid', 'zgrid', 'Pz');
