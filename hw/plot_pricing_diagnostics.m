% function plot_pricing_diagnostics(eq, grid, zgrid, Pz, par)

fprintf('\nGenerating pricing diagnostics...\n');

Nz = length(zgrid);

iz = ceil(Nz/2);            % middle productivity
ik = ceil(grid.Nk/2);       % mid capital

bvec = grid.bgrid(:);

PrND = zeros(length(bvec),1);

kp = grid.kgrid(ik);
wbar = eq.wbar(:);

prob_row = Pz(iz,:)';

profit_vec = zgrid(:) .* (kp^par.alpha);

for ib = 1:length(bvec)

    bp = bvec(ib);
    rt = eq.rtilde(ik,ib,iz);

    taxbase = profit_vec - par.delta*kp - rt*bp;
    Tc = par.tau_c_pos.*max(taxbase,0) + par.tau_c_neg.*min(taxbase,0);

    wreal = (1-par.delta)*kp + profit_vec - Tc - (1+rt)*bp;

    def = (wreal < wbar);

    PrND(ib) = sum(prob_row .* (~def));

end

figure;
subplot(3,1,1)
plot(bvec, eq.rtilde(ik,:,iz),'LineWidth',1.5)
title('Bond yield rtilde vs b')

subplot(3,1,2)
plot(bvec, PrND,'LineWidth',1.5)
title('Non-default probability vs b')

subplot(3,1,3)
plot(bvec, eq.pol_ib(:,iz),'LineWidth',1.5)
title('Chosen borrowing index vs b')

sgtitle('Hennessy-Whited pricing diagnostics')

