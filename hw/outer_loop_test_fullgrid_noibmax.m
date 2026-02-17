% main_solve.m


load('hw_solution_papergrid_feb15_nos.mat')
opts = struct();
% % initialize
opts.rtilde_init = eq.rtilde;
% Warm starts for V/policy
opts.V_init = eq.V;
opts.pol_init = struct('pol_ik', eq.pol_ik, 'pol_ib', eq.pol_ib);
clearvars -except opts

par = hw.params("full");

% Shock discretization for solve
[zgrid, Pz] = hw.tauchen(par.rho, par.sigma_eps, par.Nz_solve, par.tauchen_m);

% Capital/debt grids
grid = hw.grids(par, zgrid);

Nw = 30;
% wmax = (1-par.tau_c_pos)*grid.kbar^par.alpha/par.r;
% wgrid = linspace(-wmax, wmax, Nw);
% Revised net worth grid (debugging starter)
wgrid = linspace(-2*grid.kbar, 2*grid.kbar, Nw)';

% Solve equilibrium
opts.outer_maxit   = 200;
opts.tol_rtilde    = 1e-3;
opts.outer_verbose = true;

% Damping and caps for rtilde
opts.omega_rtilde  = 0.02;%0.02;
opts.cap_rtilde    = 0.5;%0.5;
opts.minPrND       = 1e-4;

% Damping and caps for wbar
opts.damp_wbar = true;
opts.eta_wbar = 0.5;
opts.wbar_s = 0.1;

% Inner (Howard) controls
opts.inner_maxit    = 5;
opts.inner_tol      = 1e-6;
opts.howard_iters   = 20;
opts.improve_every  = 10;
opts.inner_verbose  = true;


eq = hw.solve_equilibrium(wgrid, zgrid, Pz, par, grid, opts);

save('hw_solution_papergrid_feb16_noibmax.mat', 'eq', 'par', 'grid', 'wgrid', 'zgrid', 'Pz');
