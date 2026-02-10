function diag = hw_find_interior_default_point(eq, grid, zgrid, Pz, par)
% Returns a point (ik,ib,iz) with interior default: 0<PrND<1, if any.

Nk = grid.Nk; Nb = grid.Nb; Nz = length(zgrid);
wbar = eq.wbar(:);

diag = struct('found',false,'ik',NaN,'ib',NaN,'iz',NaN,'PrND',NaN,'shareDef',NaN);

for iz = 1:Nz
    prob = Pz(iz,:)';   % over z'
    for ik = 1:Nk
        kp = grid.kgrid(ik);
        profit_vec = zgrid(:) .* (kp^par.alpha);

        for ib = 1:Nb
            bp = grid.bgrid(ib);
            rt = eq.rtilde(ik,ib,iz);

            taxbase_vec = profit_vec - par.delta*kp - rt*bp;
            Tc_vec = par.tau_c_pos .* max(taxbase_vec,0) + par.tau_c_neg .* min(taxbase_vec,0);

            wreal_vec = (1-par.delta)*kp + profit_vec - Tc_vec - (1+rt)*bp;

            def = (wreal_vec < wbar);
            PrND = sum(prob .* (~def));
            shareDef = sum(def)/Nz;     % purely across z' grid points

            if PrND > 1e-6 && PrND < 1-1e-6
                diag.found = true;
                diag.ik = ik; diag.ib = ib; diag.iz = iz;
                diag.PrND = PrND;
                diag.shareDef = shareDef;
                return
            end
        end
    end
end
end
