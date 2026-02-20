function ibmax = feasible_ibmax(rtilde, zgrid, Pz, wbar, par, grid, opts)
%HW.FEASIBLE_IBMAX  Compute max feasible b-index for each (k,z).
%
% Implements a discrete version of the paper's restriction that we only
% consider (rtilde,b) pairs where debt value is increasing in promised yield.
% Under the maintained mapping b -> rtilde(k,b,z), this corresponds to
% restricting to the (weakly) nondecreasing part of rtilde as b increases.
%
% Additionally, we can require PrND > minPrND based on the default set.
%
% Inputs
%   rtilde : [Nk x Nb x Nz] current yield schedule
%   zgrid, Pz, wbar, par, grid : as usual
%   opts : struct (optional)
%       .eps_mono       tolerance for monotonicity, default 1e-10
%       .minPrND        default 1e-10
%       .require_PrND   (true/false), default true
%       .relax_monotonicity (true/false), default false
%           if true: skip monotonicity check, only enforce PrND if require_PrND=true
%
% Output
%   ibmax : [Nk x Nz] integer, maximum feasible ib index (>=1)

if nargin < 7, opts = struct(); end
if ~isfield(opts,'eps_mono'),           opts.eps_mono = 1e-10; end
if ~isfield(opts,'minPrND'),            opts.minPrND = 1e-10; end
if ~isfield(opts,'require_PrND'),       opts.require_PrND = true; end
if ~isfield(opts,'relax_monotonicity'), opts.relax_monotonicity = false; end

eps_mono = opts.eps_mono;
minPrND  = opts.minPrND;
relax_mono = opts.relax_monotonicity;

Nz = length(zgrid);
ibmax = ones(grid.Nk, Nz);

wbar = wbar(:);
zvec = zgrid(:);

for iz = 1:Nz
    prob_row = Pz(iz,:).';

    for ik = 1:grid.Nk
        kp = grid.kgrid(ik);

        % We'll scan b from low to high and stop when constraint breaks
        last_r = -Inf;
        last_ok = 1;

        for ib = 1:grid.Nb
            bp = grid.bgrid(ib);
            rt = rtilde(ik, ib, iz);

            % Check monotonicity (unless relaxed)
            if ~relax_mono
                if rt + eps_mono < last_r
                    break;
                end
            end
            last_r = rt;

            if opts.require_PrND && bp > 0
                % Compute PrND under current rt for default classification
                profit = zvec .* (kp^par.alpha);
                taxbase = profit - par.delta*kp - rt*bp;
                Tc = par.tau_c_pos.*max(taxbase,0) + par.tau_c_neg.*min(taxbase,0);
                wreal = (1-par.delta)*kp + profit - Tc - (1+rt)*bp;

                def = (wreal < wbar);
                PrND = sum(prob_row .* (~def));

                if PrND <= minPrND
                    break;
                end
            end

            last_ok = ib;
        end

        ibmax(ik, iz) = last_ok;
    end
end

end
