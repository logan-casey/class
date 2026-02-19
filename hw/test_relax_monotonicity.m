% test_relax_monotonicity
% Compare behavior with and without relax_monotonicity in feasible_ibmax
% Usage: run this after loading eq, par, grid, wgrid, zgrid, Pz
load('hw_solution_tightgrid_feb18_2iter_nos')
if ~exist('eq','var') || ~exist('par','var')
    error('Load eq, par, grid, wgrid, zgrid, Pz first.');
end

fprintf('\n=== STEP 1: Compare ibmax with and without relax_monotonicity ===\n');
opts_strict = struct('require_PrND', true, 'relax_monotonicity', false);
opts_relax  = struct('require_PrND', true, 'relax_monotonicity', true);

ib_strict = hw.feasible_ibmax(eq.rtilde, zgrid, Pz, eq.wbar, par, grid, opts_strict);
ib_relax  = hw.feasible_ibmax(eq.rtilde, zgrid, Pz, eq.wbar, par, grid, opts_relax);

fprintf('ibmax mean (strict mono)           : %.3f\n', mean(ib_strict(:)));
fprintf('ibmax mean (relax mono)            : %.3f\n', mean(ib_relax(:)));
fprintf('cells where ib_relax > ib_strict  : %d (%.1f%%)\n', ...
    nnz(ib_relax > ib_strict), 100*nnz(ib_relax > ib_strict)/numel(ib_strict));

% Show max gains
max_gain = max(ib_relax(:) - ib_strict(:));
fprintf('max gain in ibmax per (k,z)       : %d\n', max_gain);

fprintf('\n=== STEP 2: Run short equilibrium with relax_monotonicity ===\n');
% Warm-start with current eq
opts = struct();
opts.outer_maxit = 5;  % just 5 iterations for quick test
opts.tol_rtilde = 1e-3;
opts.outer_verbose = true;
opts.omega_rtilde = 0.05;
opts.cap_rtilde = 1e4;
opts.minPrND = 1e-12;
opts.wbar_s = 0.5;
opts.V_init = eq.V;
opts.pol_init = struct('pol_ik', eq.pol_ik, 'pol_ib', eq.pol_ib);
opts.rtilde_init = eq.rtilde;

% Key: pass relax_monotonicity to inner VFI
opts.inner_verbose = false;
opts.feas_opts = struct('relax_monotonicity', true, 'require_PrND', true, 'minPrND', 1e-12);

fprintf('Starting equilibrium solve with relax_monotonicity=true...\n');
eq_relax = hw.solve_equilibrium(wgrid, zgrid, Pz, par, grid, opts);

fprintf('\n=== STEP 3: Compare borrowing and default stats ===\n');

% Policy frequency: average ibmax used
iw_all = repmat(1:numel(wgrid), numel(zgrid), 1);
iz_all = repmat((1:numel(zgrid))', 1, numel(wgrid));

ib_choices_old = eq.pol_ib(:);
ib_choices_new = eq_relax.pol_ib(:);

fprintf('Old policy: mean b-choice (index)  : %.3f\n', mean(ib_choices_old));
fprintf('New policy: mean b-choice (index)  : %.3f\n', mean(ib_choices_new));
fprintf('Change in mean b-choice            : %.3f\n', mean(ib_choices_new) - mean(ib_choices_old));

% Estimate default frequencies (rough): count states where b > 0 and default is likely
% rough proxy: count where PrND is low
fprintf('\nEstimated state-space borrowing rates:\n');
fprintf('Old: firms with b=0                : %.1f%%\n', 100*mean(eq.pol_ib(:)==1));
fprintf('New: firms with b=0                : %.1f%%\n', 100*mean(eq_relax.pol_ib(:)==1));

% Value function: check if relaxing monotonicity improved welfare
fprintf('\nValue function comparison (rough):\n');
fprintf('Old V: mean=%.4g, std=%.4g\n', mean(eq.V(:)), std(eq.V(:)));
fprintf('New V: mean=%.4g, std=%.4g\n', mean(eq_relax.V(:)), std(eq_relax.V(:)));
fprintf('Change in mean V                   : %.4g\n', mean(eq_relax.V(:)) - mean(eq.V(:)));

% rtilde: compare spreads
fprintf('\nInterest rate comparison:\n');
fprintf('Old rtilde: mean=%.4g, std=%.4g, max=%.4g\n', mean(eq.rtilde(:)), std(eq.rtilde(:)), max(eq.rtilde(:)));
fprintf('New rtilde: mean=%.4g, std=%.4g, max=%.4g\n', mean(eq_relax.rtilde(:)), std(eq_relax.rtilde(:)), max(eq_relax.rtilde(:)));

fprintf('\nDone.\n');

