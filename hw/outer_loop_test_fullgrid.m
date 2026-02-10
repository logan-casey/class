% main_solve.m

par = hw.params("full");

% Shock discretization for solve
[zgrid, Pz] = hw.tauchen(par.rho, par.sigma_eps, par.Nz_solve, par.tauchen_m);

% Capital/debt grids
grid = hw.grids(par, zgrid);

Nw = 120;
% wmax = (1-par.tau_c_pos)*grid.kbar^par.alpha/par.r;
% wgrid = linspace(-wmax, wmax, Nw);
% Revised net worth grid (debugging starter)
wgrid = linspace(-2*grid.kbar, 2*grid.kbar, Nw);

% Solve equilibrium
opts = struct();
opts.outer_maxit   = 60;
opts.tol_rtilde    = 1e-3;
opts.outer_verbose = true;

% Damping and caps for rtilde
opts.omega_rtilde  = 0.02;
opts.cap_rtilde    = 2;
opts.minPrND       = 1e-10;

% Damping and caps for wbar
opts.damp_wbar = true;
opts.eta_wbar = 0.5;

% Inner (Howard) controls
opts.inner_maxit    = 5;
opts.inner_tol      = 1e-6;
opts.howard_iters   = 20;
opts.improve_every  = 10;
opts.inner_verbose  = true;

eq = hw.solve_equilibrium(wgrid, zgrid, Pz, par, grid, opts);

save('hw_solution_papergrid_adaptivedamping.mat', 'eq', 'par', 'grid', 'wgrid', 'zgrid', 'Pz');
