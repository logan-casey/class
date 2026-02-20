% iz = ceil(length(zgrid)/2);
ib = 50;
ik = grid.Nk;              % max k
% ik = ;
kp = grid.kgrid(ik);


for iz = 1:par.Nz_solve
    bp = grid.bgrid(ib);
    % r_guess = eq.rtilde(ik,ib,iz);
    r_guess = 0.025;
    wbar = eq.wbar(iz);

    % compute w_realized
    profit = zgrid(iz) .* (kp^par.alpha);
    taxbase = profit - par.delta*kp - r_guess*bp;
    Tc = par.tau_c_pos.*max(taxbase,0) + par.tau_c_neg.*min(taxbase,0);
    wreal = (1-par.delta)*kp + profit - Tc - (1+r_guess)*bp;

    def = (wreal < wbar);
    fprintf("iz=%.3f, wreal=%.3f, wbar=%.3f, w-wbar=%.3f, def =%i\n", iz, wreal, wbar, wreal-wbar, def)
end
