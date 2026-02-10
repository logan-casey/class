function sol = solve_vfi(wgrid, zgrid, Pz, rtilde, par, grid, opts)
%HW.SOLVE_VFI  Solve equity value function given a fixed bond yield schedule rtilde.
%
% sol = hw.solve_vfi(wgrid, zgrid, Pz, rtilde, par, grid, opts)
%
% Inputs
%   wgrid, zgrid, Pz, rtilde, par, grid as before
%   opts: struct with fields
%       .maxit   (default 500)
%       .tol     (default 1e-6)
%       .verbose (default true)
%
% Output sol struct:
%   .V         [Nw x Nz]
%   .pol_ik    [Nw x Nz] (index of kp)
%   .pol_ib    [Nw x Nz] (index of bp)
%   .wbar      [1 x Nz]
%   .it        iterations used
%   .diff      last sup norm

if nargin < 7, opts = struct(); end
if ~isfield(opts,'maxit'),   opts.maxit = 500; end
if ~isfield(opts,'tol'),     opts.tol = 1e-6; end
if ~isfield(opts,'verbose'), opts.verbose = true; end

Nw = length(wgrid);
Nz = length(zgrid);

V = zeros(Nw, Nz);          % initial guess (you can improve later)
pol_ik = ones(Nw, Nz);
pol_ib = ones(Nw, Nz);

for it = 1:opts.maxit

    % default boundary implied by current V
    wbar = hw.default_wbar(V, wgrid);

    Vnew = -Inf(Nw, Nz);

    for iz = 1:Nz
        for iw = 1:Nw
            [Vbest, ik_best, ib_best] = hw.bellman_rhs(iw, iz, V, wgrid, zgrid, Pz, rtilde, par, grid, wbar);

            Vnew(iw, iz) = Vbest;
            pol_ik(iw, iz) = ik_best;
            pol_ib(iw, iz) = ib_best;
        end
    end

    diff = max(abs(Vnew(:) - V(:)));
    V = Vnew;

    if opts.verbose && (mod(it,10)==0 || it==1)
        fprintf("VFI it=%d, supdiff=%.3e\n", it, diff);
    end

    if diff < opts.tol
        break;
    end
end

sol = struct();
sol.V      = V;
sol.pol_ik = pol_ik;
sol.pol_ib = pol_ib;
sol.wbar   = hw.default_wbar(V, wgrid);
sol.it     = it;
sol.diff   = diff;

end
