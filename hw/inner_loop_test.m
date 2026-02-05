par = hw.params("full");

% z grid
[zgrid, Pz] = hw.tauchen(par.rho, par.sigma_eps, par.Nz_solve, par.tauchen_m);

% k,b grids
grid = hw.grids(par, zgrid);

% w grid
wgrid = hw.build_wgrid(-2*grid.kbar, 2*grid.kbar, 80);

% rtilde schedule initial guess: flat at r
rtilde = par.r * ones(grid.Nk, grid.Nb, length(zgrid));

opts = struct('maxit', 50, 'tol', 1e-5, 'verbose', true, ...
    'omega', 0.2, 'cap', 10, 'minPrND', 1e-10, ...
    'howard_iters', 50, 'improve_every', 5);

sol = hw.solve_vfi_howard(wgrid, zgrid, Pz, rtilde, par, grid, opts);
