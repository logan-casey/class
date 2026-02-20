function sol = solve_vfi_howard(wgrid, zgrid, Pz, rtilde, par, grid, opts)
% Solve Bellman

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

for it = 1:opts.maxit

    % Build interpolants for current V
    F = cell(Nz,1);
    for izp = 1:Nz
        F{izp} = griddedInterpolant(wgrid, V(:,izp), 'linear', 'nearest');
    end

    % Update wbar
    wbar_used = hw.default_wbar(V, wgrid);

    % Update feasibility
    ibmax_used = hw.feasible_ibmax(rtilde, zgrid, Pz, wbar_used, par, grid, opts.feas_opts);


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

    if diff_imp < opts.tol
        break;
    end
end

sol = struct();
sol.V      = V;
sol.pol_ik = pol_ik;
sol.pol_ib = pol_ib;

% For output: report wbar implied by the final V
sol.wbar = hw.default_wbar(V, wgrid);

sol.it   = it;
sol.diff = 0;

end
