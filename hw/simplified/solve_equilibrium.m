function eq = solve_equilibrium(wgrid, zgrid, Pz, par, grid, opts)
%HW.SOLVE_EQUILIBRIUM  Full solution with outer iteration on bond yields.
%
% eq = hw.solve_equilibrium(wgrid, zgrid, Pz, par, grid, opts)
%
% Outer loop:
%   1) Given rtilde_old, solve Bellman -> V, policy
%   2) Compute wbar from V
%   3) Update rtilde_new using discrete eq. (20)
%   4) Iterate until convergence

Nz = length(zgrid);

% Initialize rtilde
if isempty(opts.rtilde_init)
    rtilde = par.r * ones(grid.Nk, grid.Nb, Nz);
else
    rtilde = opts.rtilde_init;
end

inner_opts = opts;

for outer_it = 1:opts.outer_maxit

    inner_opts.V_init        = V_init;
    inner_opts.pol_init      = pol_init;

    % SOLVE BELLMAN FOR V, pol
    sol = hw.solve_vfi_howard(wgrid, zgrid, Pz, rtilde, par, grid, inner_opts);

    V      = sol.V;
    pol_ik = sol.pol_ik;
    pol_ib = sol.pol_ib;
    wbar   = sol.wbar;

    % GET Rtilde
    [rtilde_new, ~] = hw.update_rtilde(rtilde, zgrid, Pz, wbar, par, grid, opts);

    diff_r = max(abs(rtilde_new(:) - rtilde(:)));

    % Update warm starts
    rtilde = rtilde_new;
    V_init = V;
    pol_init = struct('pol_ik', pol_ik, 'pol_ib', pol_ib);

    % Stopping rule
    if diff_r < opts.tol_rtilde
        break;
    end

eq = struct();
eq.rtilde   = rtilde;
eq.V        = V;
eq.pol_ik   = pol_ik;
eq.pol_ib   = pol_ib;
eq.wbar     = wbar;
eq.outer_it = outer_it;

end