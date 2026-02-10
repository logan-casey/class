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
%
% Inputs
%   wgrid, zgrid, Pz : grids and transition
%   par, grid        : parameter and grid structs
%   opts             : struct with fields (defaults below)
%       OUTER:
%         .outer_maxit    = 50
%         .tol_rtilde     = 1e-4     % sup norm on rtilde
%         .tol_V          = 1e-5     % optional, sup norm on V
%         .outer_verbose  = true
%         .omega_rtilde   = 0.2      % damping used in update_rtilde
%         .cap_rtilde     = 10       % cap for yields
%         .minPrND        = 1e-10
%       INNER (Howard):
%         .inner_maxit     = 50
%         .inner_tol       = 1e-6
%         .howard_iters    = 50
%         .improve_every   = 5
%         .inner_verbose   = false
%       INIT:
%         .rtilde_init     = []       % if empty, flat r
%         .V_init          = []       % warm start for first outer iter
%         .pol_init        = []       % warm start policy for first outer iter
%
% Output eq struct:
%   .rtilde, .V, .pol_ik, .pol_ib, .wbar
%   .outer_it, .outer_diff_r, .outer_diff_V
%   .history (optional): diffs each iteration

if nargin < 6, opts = struct(); end

% ---- defaults: outer ----
if ~isfield(opts,'outer_maxit'),   opts.outer_maxit = 50; end
if ~isfield(opts,'tol_rtilde'),    opts.tol_rtilde = 1e-4; end
if ~isfield(opts,'tol_V'),         opts.tol_V = 1e-5; end
if ~isfield(opts,'outer_verbose'), opts.outer_verbose = true; end
if ~isfield(opts,'omega_rtilde'),  opts.omega_rtilde = 0.2; end
if ~isfield(opts,'cap_rtilde'),    opts.cap_rtilde = 10; end
if ~isfield(opts,'damp_wbar'),     opts.damp_wbar = false; end
if ~isfield(opts,'eta_wbar'),      opts.eta_wbar = 0; end
if ~isfield(opts,'minPrND'),       opts.minPrND = 1e-10; end

% ---- defaults: inner (Howard) ----
if ~isfield(opts,'inner_maxit'),    opts.inner_maxit = 50; end
if ~isfield(opts,'inner_tol'),      opts.inner_tol = 1e-6; end
if ~isfield(opts,'howard_iters'),   opts.howard_iters = 50; end
if ~isfield(opts,'improve_every'),  opts.improve_every = 5; end
if ~isfield(opts,'inner_verbose'),  opts.inner_verbose = false; end

% ---- init ----
if ~isfield(opts,'rtilde_init'), opts.rtilde_init = []; end
if ~isfield(opts,'V_init'),      opts.V_init = []; end
if ~isfield(opts,'pol_init'),    opts.pol_init = []; end

Nz = length(zgrid);

% Initialize rtilde
if isempty(opts.rtilde_init)
    rtilde = par.r * ones(grid.Nk, grid.Nb, Nz);
else
    rtilde = opts.rtilde_init;
end

% Warm starts for V/policy
V_init   = opts.V_init;
pol_init = opts.pol_init;

history.diff_r = zeros(opts.outer_maxit,1);
history.diff_V = zeros(opts.outer_maxit,1);

V_prev = [];

for outer_it = 1:opts.outer_maxit

    % ---- Inner solve: VFI with Howard ----
    inner_opts = struct();
    inner_opts.maxit         = opts.inner_maxit;
    inner_opts.tol           = opts.inner_tol;
    inner_opts.verbose       = opts.inner_verbose;
    inner_opts.howard_iters  = opts.howard_iters;
    inner_opts.improve_every = opts.improve_every;
    inner_opts.V_init        = V_init;
    inner_opts.pol_init      = pol_init;

    sol = hw.solve_vfi_howard(wgrid, zgrid, Pz, rtilde, par, grid, inner_opts);

    V      = sol.V;
    pol_ik = sol.pol_ik;
    pol_ib = sol.pol_ib;
    wbar_new   = sol.wbar;

    % ---- Outer update of rtilde via equation (20) ----
    upd_opts = struct();
    upd_opts.omega   = opts.omega_rtilde;
    upd_opts.cap     = opts.cap_rtilde;
    upd_opts.minPrND = opts.minPrND;

    % damping for better solution
    if opts.damp_wbar
        eta = opts.eta_wbar;
    else
        eta = 1;
    end
    if ~exist('wbar_old','var')
        wbar_old = wbar_new;
    end
    wbar = (1-eta)*wbar_old + eta*wbar_new;
    wbar_old = wbar;
    rtilde_new = hw.update_rtilde(rtilde, zgrid, Pz, wbar, par, grid, upd_opts);

    diff_r = max(abs(rtilde_new(:) - rtilde(:)));
    history.diff_r(outer_it) = diff_r;

    % Optional V convergence diagnostic across outer iterations
    if isempty(V_prev)
        diff_V = NaN;
    else
        diff_V = max(abs(V(:) - V_prev(:)));
    end
    history.diff_V(outer_it) = diff_V;

    if opts.outer_verbose
        if isnan(diff_V)
            fprintf("Outer it=%d: supdiff rtilde=%.3e\n", outer_it, diff_r);
        else
            fprintf("Outer it=%d: supdiff rtilde=%.3e, supdiff V=%.3e\n", outer_it, diff_r, diff_V);
        end
    end

    % Update iterates + warm starts
    rtilde = rtilde_new;
    V_prev = V;
    V_init = V;
    pol_init = struct('pol_ik', pol_ik, 'pol_ib', pol_ib);

    % Stopping rule: primarily rtilde convergence (matches paper description)
    if diff_r < opts.tol_rtilde
        if opts.outer_verbose
            fprintf("Converged on rtilde: %.3e < %.3e at outer it=%d\n", diff_r, opts.tol_rtilde, outer_it);
        end
        break;
    end

    % Optional secondary stop if V stable too
    if ~isnan(diff_V) && diff_r < 10*opts.tol_rtilde && diff_V < opts.tol_V
        if opts.outer_verbose
            fprintf("Converged on rtilde and V at outer it=%d\n", outer_it);
        end
        break;
    end
end

% Trim history
history.diff_r = history.diff_r(1:outer_it);
history.diff_V = history.diff_V(1:outer_it);

eq = struct();
eq.rtilde   = rtilde;
eq.V        = V;
eq.pol_ik   = pol_ik;
eq.pol_ib   = pol_ib;
eq.wbar     = wbar;
eq.outer_it = outer_it;
eq.outer_diff_r = history.diff_r(end);
eq.outer_diff_V = history.diff_V(end);
eq.history  = history;

end