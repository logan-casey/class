function [Vbest, ik_best, ib_best] = bellman_rhs(iw, iz, V, wgrid, zgrid, Pz, rtilde, par, grid, wbar)
%HW.BELLMAN_RHS  One-step Bellman maximization at state (wtilde_i, z_i).
%
% State:
%   wtilde = wgrid(iw), z = zgrid(iz)
%
% Controls (discrete search):
%   kp in grid.kgrid, bp in grid.bgrid (no additional feasibility filter yet)
%
% Inputs
%   iw, iz   : indices for current state
%   V        : current value function guess [Nw x Nz]
%   wgrid    : revised net worth grid [Nw x 1]
%   zgrid    : z grid [Nz x 1]
%   Pz       : transition matrix [Nz x Nz]
%   rtilde   : bond yield schedule [Nk x Nb x Nz] (indexed by kp, bp, iz)
%   par      : parameters
%   grid     : struct with kgrid, bgrid
%   wbar     : default boundary [1 x Nz] computed from current V
%
% Outputs
%   Vbest, ik_best, ib_best : maximizing value and argmax indices

wtilde = wgrid(iw);
Nz = length(zgrid);

Vbest   = -Inf;
ik_best = NaN;
ib_best = NaN;

for ik = 1:grid.Nk
    kp = grid.kgrid(ik);

    for ib = 1:grid.Nb
        bp = grid.bgrid(ib);

        % Net distribution to shareholders (can be negative = issuance)
        D = wtilde + bp - kp;

        if D >= 0
            payout = D;
            flow = payout - hw.tax_dist(payout, par);
        else
            E = -D; % issuance amount
            flow = -E - hw.equity_cost(E, par);
        end

        % continuation value
        EV = 0;
        rt = rtilde(ik, ib, iz);

        for izp = 1:Nz
            prob = Pz(iz, izp);
            zprime = zgrid(izp);

            % realized net worth
            w_real = hw.realized_networth(kp, bp, zprime, rt, par);

            % revised net worth (default reset)
            w_next = hw.revised_networth(w_real, wbar(izp));

            % interpolate V(w_next, izp)
            Vnext = interp1_clamped(wgrid, V(:, izp), w_next);

            EV = EV + prob * Vnext;
        end

        Vcand = flow + (1 / (1 + par.r)) * EV;

        if Vcand > Vbest
            Vbest   = Vcand;
            ik_best = ik;
            ib_best = ib;
        end
    end
end

end

function v = interp1_clamped(x, y, xq)
% Linear interpolation on a grid, clamping queries to [min(x), max(x)].
xq = min(max(xq, x(1)), x(end));
v = interp1(x, y, xq, 'linear');
end
