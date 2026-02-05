function sol = solve_vfi_howard(wgrid, zgrid, Pz, rtilde, par, grid, opts)
%HW.SOLVE_VFI_HOWARD  VFI with Howard policy improvement.
%
% Enhancements:
%   - Freeze feasibility: compute ibmax ONCE and hold fixed inside the solver
%   - Optionally freeze wbar inside the solver (recommended for stability)
%
% opts additional fields (defaults):
%   .freeze_feasible = true
%   .freeze_wbar     = true
%   .wbar_fixed      = []      % if provided, use this wbar (overrides freeze_wbar)
%   .feas_opts       = struct('require_PrND',true,'minPrND',1e-10,'eps_mono',1e-10)

if nargin < 7, opts = struct(); end
if ~isfield(opts,'maxit'),         opts.maxit = 200; end
if ~isfield(opts,'tol'),           opts.tol = 1e-6; end
if ~isfield(opts,'verbose'),       opts.verbose = true; end
if ~isfield(opts,'howard_iters'),  opts.howard_iters = 20; end
if ~isfield(opts,'improve_every'), opts.improve_every = 1; end
if ~isfield(opts,'V_init'),        opts.V_init = []; end
if ~isfield(opts,'pol_init'),      opts.pol_init = []; end

% new defaults
if ~isfield(opts,'freeze_feasible'), opts.freeze_feasible = true; end
if ~isfield(opts,'freeze_wbar'),     opts.freeze_wbar = true; end
if ~isfield(opts,'wbar_fixed'),      opts.wbar_fixed = []; end
if ~isfield(opts,'feas_opts') || isempty(opts.feas_opts)
    opts.feas_opts = struct('require_PrND', true, 'minPrND', 1e-10, 'eps_mono', 1e-10);
end

Nw = length(wgrid);
Nz = length(zgrid);

% Initialize V
if ~isempty(opts.V_init)
    V = opts.V_init;
else
    V = zeros(Nw, Nz);
end

% Initialize policy
if ~isempty(opts.pol_init) && isfield(opts.pol_init,'pol_ik') && isfield(opts.pol_init,'pol_ib')
    pol_ik = opts.pol_init.pol_ik;
    pol_ib = opts.pol_init.pol_ib;
else
    pol_ik = ones(Nw, Nz);
    pol_ib = ones(Nw, Nz);
end

Vimp = V; % ensure defined

% ------------------------------------------------------------
% Compute (frozen) wbar to be used inside the solver
% ------------------------------------------------------------
if ~isempty(opts.wbar_fixed)
    wbar_used = opts.wbar_fixed(:).';
else
    wbar_used = hw.default_wbar(V, wgrid);
    if ~opts.freeze_wbar
        % if not freezing, we will refresh wbar_used each iteration below
    end
end

% ------------------------------------------------------------
% Compute (frozen) feasibility ibmax (once) if requested
% ------------------------------------------------------------
if opts.freeze_feasible
    ibmax_used = hw.feasible_ibmax(rtilde, zgrid, Pz, wbar_used, par, grid, opts.feas_opts);
else
    ibmax_used = ones(grid.Nk, Nz);
end

for it = 1:opts.maxit

    % Build interpolants for current V
    F = cell(Nz,1);
    for izp = 1:Nz
        F{izp} = griddedInterpolant(wgrid, V(:,izp), 'linear', 'nearest');
    end

    % Update wbar if not frozen and not externally fixed
    if isempty(opts.wbar_fixed) && ~opts.freeze_wbar
        wbar_used = hw.default_wbar(V, wgrid);
    end

    % Update feasibility if not frozen
    if ~opts.freeze_feasible
        ibmax_used = hw.feasible_ibmax(rtilde, zgrid, Pz, wbar_used, par, grid, opts.feas_opts);
    end

    % ------------------------------------------------------------
    % (A) Policy Improvement (expensive)
    % ------------------------------------------------------------
    if mod(it-1, opts.improve_every) == 0
        Vimp = -Inf(Nw, Nz);

        for iz = 1:Nz
            for iw = 1:Nw
                [Vbest, ik_best, ib_best] = hw.bellman_rhs_fast( ...
                    iw, iz, V, wgrid, zgrid, Pz, rtilde, par, grid, wbar_used, F, ibmax_used);

                Vimp(iw, iz)   = Vbest;
                pol_ik(iw, iz) = ik_best;
                pol_ib(iw, iz) = min(ib_best, ibmax_used(ik_best, iz));
            end
        end

        diff_imp = max(abs(Vimp(:) - V(:)));
        V = Vimp;

        if opts.verbose
            fprintf("Howard outer it=%d: after improve supdiff=%.3e\n", it, diff_imp);
        end

        if diff_imp < opts.tol
            break;
        end
    end

    % Rebuild interpolants after improvement (V changed)
    for izp = 1:Nz
        F{izp} = griddedInterpolant(wgrid, V(:,izp), 'linear', 'nearest');
    end

    % If not frozen, refresh wbar_used after the improvement step
    if isempty(opts.wbar_fixed) && ~opts.freeze_wbar
        wbar_used = hw.default_wbar(V, wgrid);
    end

    % If not frozen, refresh feasibility after improvement too
    if ~opts.freeze_feasible
        ibmax_used = hw.feasible_ibmax(rtilde, zgrid, Pz, wbar_used, par, grid, opts.feas_opts);
    end

    % ------------------------------------------------------------
    % (B) Howard Policy Evaluation: fixed policy
    % ------------------------------------------------------------
    for h = 1:opts.howard_iters
        Vnew = zeros(Nw, Nz);

        for iz = 1:Nz
            prob_row = Pz(iz, :).';
            zvec = zgrid(:);
            wbarvec = wbar_used(:);

            for iw = 1:Nw
                wtilde = wgrid(iw);

                ik = pol_ik(iw, iz);
                ib = pol_ib(iw, iz);

                % clamp to frozen/current feasibility set
                ib = min(ib, ibmax_used(ik, iz));

                kp = grid.kgrid(ik);
                bp = grid.bgrid(ib);
                rt = rtilde(ik, ib, iz);

                % one-period flow
                % D = wtilde + bp - kp;
                % ********************************* questionable
                proceeds = bp / (1 + par.r*(1 - par.tau_i));
                D = wtilde + proceeds - kp;

                if D >= 0
                    flow = D - hw.tax_dist(D, par);
                else
                    E = -D;
                    flow = -E - hw.equity_cost(E, par);
                end

                % continuation across z'
                profit_vec = zvec .* (kp^par.alpha);
                taxbase_vec = profit_vec - par.delta*kp - rt*bp;
                Tc_vec = par.tau_c_pos .* max(taxbase_vec, 0) + par.tau_c_neg .* min(taxbase_vec, 0);

                w_real_vec = (1 - par.delta)*kp + profit_vec - Tc_vec - (1 + rt)*bp;
                w_next_vec = max(wbarvec, w_real_vec);

                Vnext_vec = zeros(Nz,1);
                for izp = 1:Nz
                    Vnext_vec(izp) = F{izp}(w_next_vec(izp));
                end

                EV = prob_row' * Vnext_vec;

                Vnew(iw, iz) = flow + (1/(1+par.r)) * EV;
            end
        end

        V = Vnew;
    end

    % check convergence across outer loops (after Howard eval)
    diff = max(abs(V(:) - Vimp(:)));
    if opts.verbose
        fprintf("Howard outer it=%d: end-of-outer diff~%.3e\n", it, diff);
    end

    if diff < opts.tol
        break;
    end
end

sol = struct();
sol.V      = V;
sol.pol_ik = pol_ik;
sol.pol_ib = pol_ib;

% For output: report wbar implied by the final V (even if we froze inside)
sol.wbar = hw.default_wbar(V, wgrid);

sol.it   = it;
sol.diff = 0;

end
