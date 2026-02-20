% loads parameters, intitializes, calls solve_equilibrium.m

% initialize
load('')
opts = struct();
opts.rtilde_init = eq.rtilde;
opts.V_init = eq.V;
opts.pol_init = struct('pol_ik', eq.pol_ik, 'pol_ib', eq.pol_ib);
clearvars -except opts

par = hw.params("full");

% Shock discretization
[zgrid, Pz] = hw.tauchen(par.rho, par.sigma_eps, par.Nz_solve, par.tauchen_m);

% Capital/debt grids
grid = hw.grids(par, zgrid);

Nw = 80;
wmax = (1-par.tau_c_pos)*grid.kbar^par.alpha/par.r;
wgrid = linspace(-wmax, wmax, Nw);

% Solve options
opts.outer_maxit   = 1000;
opts.tol_rtilde    = 1e-3;
opts.outer_verbose = true;
% Inner (Howard) controls
opts.inner_maxit    = 5;
opts.inner_tol      = 1e-6;
opts.howard_iters   = 20;
opts.improve_every  = 10;
opts.inner_verbose  = true;

% Damping and caps for rtilde
opts.omega_rtilde  = 0.02;
opts.cap_rtilde    = 1e4;
opts.minPrND       = 1e-4;

% Damping and caps for wbar
opts.damp_wbar = false;
opts.eta_wbar = 0.5;
opts.wbar_s = 0;

eq = hw.solve_equilibrium(wgrid, zgrid, Pz, par, grid, opts);

save('.mat', 'eq', 'par', 'grid', 'wgrid', 'zgrid', 'Pz');
