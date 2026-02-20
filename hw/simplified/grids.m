function grid = grids(par, zgrid)
%HW.GRIDS  Build capital and debt grids given parameters and z grid.
%
% grid = hw.grids(par, zgrid)
%
% Outputs:
%   grid.kbar, grid.kgrid
%   grid.bmin, grid.bmax, grid.bgrid
%   grid.Nk, grid.Nb

if ~isfield(par, "alpha") || isnan(par.alpha)
    error("par.alpha must be set (not NaN).");
end
if ~isfield(par, "delta")
    error("par.delta missing.");
end
if isempty(zgrid)
    error("zgrid is empty.");
end

grid = struct();
zbar = max(zgrid);

% kbar solves: zbar * alpha * kbar^(alpha-1) = delta
grid.kbar = ((zbar * par.alpha) / par.delta)^(1/(1 - par.alpha));

% Capital grid as in paper: kbar*(1-delta)^[0, 0.5, 1, ..., 15]
if ~isfield(par, "k_exponents") || isempty(par.k_exponents)
    error("Need par.k_exponents")
else
    exponents = par.k_exponents;
end
grid.kgrid = grid.kbar .* (1 - par.delta) .^ exponents(:);
grid.Nk = numel(grid.kgrid);

% Debt grid: half as many points; symmetric endpoints
% bmax = (1 - tau_c_pos) * kbar^alpha / r
if ~isfield(par, "tau_c_pos") || ~isfield(par, "r")
    error("Need par.tau_c_pos and par.r.");
end

bmax = (1 - par.tau_c_pos) * (grid.kbar^par.alpha) / par.r;
grid.bmin = -bmax;
grid.bmax =  bmax;

% Ensure Nb integer
if ~isfield(par, "Nb") || isempty(par.Nb)
    Nb = grid.Nk/2;
else
    Nb = par.Nb;
end
if abs(Nb - round(Nb)) > 1e-12
    Nb = round(Nb);
end

grid.Nb = Nb;
grid.bgrid = linspace(grid.bmin, grid.bmax, grid.Nb)';

end