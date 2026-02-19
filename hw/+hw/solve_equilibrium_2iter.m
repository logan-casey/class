function eq = solve_equilibrium_2iter(wgrid, zgrid, Pz, par, grid, opts)
%HW.SOLVE_EQUILIBRIUM  Nested fixed point:
%   (A) Solve VFI to contraction holding rtilde fixed
%   (B) Update rtilde until pricing errors are small holding (V,wbar) fixed
%   Repeat A-B until joint convergence.
%
% Compared to the previous version, we DO NOT interleave 1 rtilde update per VFI solve.
% Instead we "block solve" each mapping more tightly to reduce oscillations from the
% discontinuous default set.

if nargin < 6, opts = struct(); end

% ---------------- defaults: outer ----------------
if ~isfield(opts,'outer_maxit'),      opts.outer_maxit = 50; end
if ~isfield(opts,'outer_verbose'),    opts.outer_verbose = true; end
if ~isfield(opts,'tol_rtilde'),       opts.tol_rtilde = 1e-4; end
if ~isfield(opts,'tol_V'),            opts.tol_V = 1e-5; end

% ---------------- defaults: VFI block ----------------
if ~isfield(opts,'inner_maxit'),      opts.inner_maxit = 50; end
if ~isfield(opts,'inner_tol'),        opts.inner_tol = 1e-6; end
if ~isfield(opts,'howard_iters'),     opts.howard_iters = 50; end
if ~isfield(opts,'improve_every'),    opts.improve_every = 5; end
if ~isfield(opts,'inner_verbose'),    opts.inner_verbose = false; end

% ---------------- defaults: pricing block ----------------
if ~isfield(opts,'price_maxit'),      opts.price_maxit = 30; end
if ~isfield(opts,'tol_price'),        opts.tol_price = 1e-4; end  % target on pricing error (preferred)
if ~isfield(opts,'omega_rtilde'),     opts.omega_rtilde = 0.2; end
if ~isfield(opts,'cap_rtilde'),       opts.cap_rtilde = 10; end
if ~isfield(opts,'minPrND'),          opts.minPrND = 1e-10; end

% ---------------- defaults: wbar damping ----------------
if ~isfield(opts,'damp_wbar'),        opts.damp_wbar = false; end
if ~isfield(opts,'eta_wbar'),         opts.eta_wbar = 1.0; end

% ---------------- defaults: init ----------------
if ~isfield(opts,'rtilde_init'),      opts.rtilde_init = []; end
if ~isfield(opts,'V_init'),           opts.V_init = []; end
if ~isfield(opts,'pol_init'),         opts.pol_init = []; end

% ---------------- defaults: update_rtilde smoothing ----------------
if ~isfield(opts,'smooth_default'),   opts.smooth_default = false; end
if ~isfield(opts,'smooth_s'),         opts.smooth_s = 10; end  % only used if smooth_default=true

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

% History containers
history.outer.diff_r    = zeros(opts.outer_maxit,1);
history.outer.diff_V    = zeros(opts.outer_maxit,1);
history.outer.price_err = zeros(opts.outer_maxit,1);
history.outer.price_it  = zeros(opts.outer_maxit,1);

V_prev = [];
pol_change = NaN;

for outer_it = 1:opts.outer_maxit

    % ============================================================
    % (A) VFI BLOCK: solve Bellman to contraction holding rtilde fixed
    % ============================================================
    inner_opts = struct();
    inner_opts.maxit           = opts.inner_maxit;
    inner_opts.tol             = opts.inner_tol;
    inner_opts.verbose         = opts.inner_verbose;
    inner_opts.howard_iters    = opts.howard_iters;
    inner_opts.improve_every   = opts.improve_every;
    inner_opts.V_init          = V_init;
    inner_opts.pol_init        = pol_init;
    inner_opts.freeze_feasible = false;
    inner_opts.feas_opts = opts.feas_opts;

    sol = hw.solve_vfi_howard(wgrid, zgrid, Pz, rtilde, par, grid, inner_opts);

    V      = sol.V;
    pol_ik = sol.pol_ik;
    pol_ib = sol.pol_ib;
    wbar_new = sol.wbar(:).';

    % Policy change diagnostic
    if outer_it > 1
        pol_change = mean( (pol_ik(:) ~= pol_ik_prev(:)) | (pol_ib(:) ~= pol_ib_prev(:)) );
    else
        pol_change = NaN;
    end
    pol_ik_prev = pol_ik; pol_ib_prev = pol_ib;

    % wbar damping (optional)
    if ~exist('wbar_old','var')
        wbar_old = wbar_new;
    end
    if opts.damp_wbar
        eta = opts.eta_wbar;
    else
        eta = 1.0;
    end
    wbar = (1-eta)*wbar_old + eta*wbar_new;
    wbar_old = wbar;

    % Monotonicity diagnostics (optional)
    if opts.outer_verbose
        mono_viol = max(max(max(0, V(:,1:end-1) - V(:,2:end))));
        fprintf("max violation V(w,z) nondecreasing in z: %.3e\n", mono_viol);
        wbar_viol = max(max(0, diff(wbar))); % positive diff means increasing -> violation
        fprintf("max violation wbar nonincreasing in z: %.3e\n", wbar_viol);
    end

    % V convergence diagnostic across outer blocks
    if isempty(V_prev)
        diff_V = NaN;
    else
        diff_V = max(abs(V(:) - V_prev(:)));
    end
    history.outer.diff_V(outer_it) = diff_V;

    % ============================================================
    % (B) PRICING BLOCK: iterate update_rtilde holding (V,wbar) fixed
    % ============================================================
    upd_opts = struct();
    upd_opts.cap            = opts.cap_rtilde;
    upd_opts.minPrND        = opts.minPrND;
    upd_opts.smooth_default = opts.smooth_default;
    upd_opts.smooth_s       = opts.smooth_s;

    % You can still adapt omega based on pol_change, but keep it FIXED inside the pricing loop.
    omega = opts.omega_rtilde;
    if ~isnan(pol_change) && pol_change > 0.10
        omega_eff = min(omega, 0.01);
    elseif ~isnan(pol_change) && pol_change > 0.03
        omega_eff = min(omega, 0.02);
    else
        omega_eff = omega;
    end
    upd_opts.omega = omega_eff;

    price_err = NaN;
    diff_r    = NaN;
    def_diag  = struct();

    for pit = 1:opts.price_maxit

        rtilde_old = rtilde;
        [rtilde_new, def_diag] = hw.update_rtilde(rtilde_old, zgrid, Pz, wbar, par, grid, upd_opts);

        diff_r = max(abs(rtilde_new(:) - rtilde_old(:)));

        % Prefer a true "pricing error" if update_rtilde returns it; otherwise fall back to diff_r.
        % Recommended convention: def_diag.rerr is an array of pricing residuals (in bp or yield units).
        if isfield(def_diag,'rerr') && ~isempty(def_diag.rerr)
            price_err = max(abs(def_diag.rerr(:)));
        elseif isfield(def_diag,'price_err') && ~isempty(def_diag.price_err)
            price_err = def_diag.price_err;
        else
            price_err = diff_r;
        end

        rtilde = rtilde_new;

        if opts.outer_verbose
            fprintf("  Pricing it=%d: omega=%.3g, supdiff r=%.3e, price_err=%.3e\n", ...
                pit, omega_eff, diff_r, price_err);
        end

        % Stop pricing block when pricing residuals are small (preferred),
        % otherwise when rtilde is stable.
        if price_err < opts.tol_price || diff_r < opts.tol_rtilde
            break;
        end
    end

    history.outer.diff_r(outer_it)    = diff_r;
    history.outer.price_err(outer_it) = price_err;
    history.outer.price_it(outer_it)  = pit;

    if opts.outer_verbose
        if isnan(diff_V)
            fprintf("Outer block=%d: pol_change=%.3e, price_it=%d, price_err=%.3e, diff_r=%.3e\n", ...
                outer_it, pol_change, pit, price_err, diff_r);
        else
            fprintf("Outer block=%d: pol_change=%.3e, supdiff V=%.3e, price_it=%d, price_err=%.3e, diff_r=%.3e\n", ...
                outer_it, pol_change, diff_V, pit, price_err, diff_r);
        end
    end

    % Warm starts for next outer block
    V_prev  = V;
    V_init  = V;
    pol_init = struct('pol_ik', pol_ik, 'pol_ib', pol_ib);

    % Joint stopping rule
    conv_r = (diff_r < opts.tol_rtilde) || (price_err < opts.tol_price);
    conv_V = (~isnan(diff_V) && diff_V < opts.tol_V);

    if conv_r && conv_V
        if opts.outer_verbose
            fprintf("Converged: (pricing/rtilde) and V at outer block=%d\n", outer_it);
        end
        break;
    end
end

% Trim history
history.outer.diff_r    = history.outer.diff_r(1:outer_it);
history.outer.diff_V    = history.outer.diff_V(1:outer_it);
history.outer.price_err = history.outer.price_err(1:outer_it);
history.outer.price_it  = history.outer.price_it(1:outer_it);

% Output eq
eq = struct();
eq.rtilde   = rtilde;
eq.V        = V;
eq.pol_ik   = pol_ik;
eq.pol_ib   = pol_ib;
eq.wbar     = wbar;
eq.outer_it = outer_it;

eq.outer_diff_r    = history.outer.diff_r(end);
eq.outer_diff_V    = history.outer.diff_V(end);
eq.outer_price_err = history.outer.price_err(end);

eq.history  = history;
eq.def_diag = def_diag;   % last pricing-block diagnostics

end
