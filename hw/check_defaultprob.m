iz = ceil(length(zgrid)/2);
ik = grid.Nk;              % max k
kp = grid.kgrid(ik);

wbar = eq.wbar(:);
prob = Pz(iz,:)';

for ib = 1:grid.Nb
    bp = grid.bgrid(ib);
    r_guess = eq.rtilde(ik,ib,iz);

    % compute w_realized across z'
    profit = zgrid(:) .* (kp^par.alpha);
    taxbase = profit - par.delta*kp - r_guess*bp;
    Tc = par.tau_c_pos.*max(taxbase,0) + par.tau_c_neg.*min(taxbase,0);
    wreal = (1-par.delta)*kp + profit - Tc - (1+r_guess)*bp;

    def = (wreal < wbar);
    PrND = sum(prob .* (~def));

    % recovery ratios
    R = zeros(length(zgrid),1);
    for izp=1:length(zgrid)
        R(izp) = hw.recovery_R(kp, zgrid(izp), wbar(izp), par);
    end
    ED = sum(prob .* def .* (R ./ max(bp,1e-12)));

    fprintf("ib=%d, b=%.2f, r_guess=%.3f, PrND=%.3g, ED=%.3g\n", ...
            ib, bp, r_guess, PrND, ED);
end
