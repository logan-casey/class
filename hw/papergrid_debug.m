
load('hw_solution_papergrid_feb16_nodamping_b')
opts = struct();
opts.outer_maxit = 10;
opts.tol_rtilde = 1e-3;
opts.outer_verbose = true;
opts.omega_rtilde = 0.01;
opts.cap_rtilde = 0.5;
opts.minPrND = 1e-12;
opts.wbar_s = 10;
opts.V_init = eq.V;
opts.pol_init = struct('pol_ik', eq.pol_ik, 'pol_ib', eq.pol_ib);
opts.rtilde_init = eq.rtilde;
opts.inner_verbose = false;
opts.feas_opts = struct('require_PrND', true, 'minPrND', 1e-12, 'eps_mono', 1e-4);  % relax eps_mono

eq_relax_tol = hw.solve_equilibrium(wgrid, zgrid, Pz, par, grid, opts);

% Compare
fprintf('Old mean b-choice: %.3f\n', mean(eq.pol_ib(:)));
fprintf('New mean b-choice: %.3f\n', mean(eq_relax_tol.pol_ib(:)));
fprintf('New mean V: %.4g (vs old %.4g)\n', mean(eq_relax_tol.V(:)), mean(eq.V(:)));
fprintf('New rtilde: mean=%.6g, max=%.6g\n', mean(eq_relax_tol.rtilde(:)), max(eq_relax_tol.rtilde(:)));

save('papergrid_debugs.mat','-mat')