% function hw_diag_pricing_curve(eq, grid, zgrid, Pz, par)

Nz = length(zgrid);
Nb = grid.Nb;

% pick representative indices
diag = hw_find_interior_default_point(eq, grid, zgrid, Pz, par);
if ~diag.found
    warning('No interior-default points found on current grid / rtilde. Diagnostics will be degenerate.');
else
    fprintf('Interior default at (ik,ib,iz)=(%d,%d,%d): PrND=%.4f, shareDef(zprime grid)=%.2f\n', ...
        diag.ik, diag.ib, diag.iz, diag.PrND, diag.shareDef);
end

ik = diag.ik; ib = diag.ib; iz = diag.iz;
kp = grid.kgrid(ik);

prob = Pz(iz,:)';         % distribution over z' given current z
wbar = eq.wbar(:);

bvec = grid.bgrid(:);

PrND = zeros(Nb,1);
ED   = zeros(Nb,1);
rerr = zeros(Nb,1);

profit_vec = zgrid(:) .* (kp^par.alpha);

% lender discount factor if pricing is risk-free after taxes (match your model)
betaL = 1/(1 + par.r*(1-par.tau_i));

for ib = 1:Nb
    bp = bvec(ib);
    rt = eq.rtilde(ik,ib,iz);

    % realized w across z'
    taxbase_vec = profit_vec - par.delta*kp - rt*bp;
    Tc_vec = par.tau_c_pos .* max(taxbase_vec,0) + par.tau_c_neg .* min(taxbase_vec,0);

    wreal_vec = (1-par.delta)*kp + profit_vec - Tc_vec - (1+rt)*bp;

    def = (wreal_vec < wbar);
    PrND(ib) = sum(prob .* (~def));

    % recovery in default states
    Rvec = zeros(Nz,1);
    for izp = 1:Nz
        Rvec(izp) = hw.recovery_R(kp, zgrid(izp), wbar(izp), par);
    end

    % expected discounted payoff to lender (in units of face value)
    % payoff = bp if repay; = R(z') if default
    payoff = (~def).*bp + def.*Rvec;
    PV = betaL * (prob' * payoff);

    % pricing error: PV should equal price paid today.
    % If you price at q = bp/(1 + rt*(1-par.tau_i)) (your “proceeds”),
    % then error = PV - q.
    q = bp;
    if bp > 0
        q = bp / (1 + rt*(1-par.tau_i));
    end

    rerr(ib) = PV - q;

    % optional: expected recovery ratio conditional on default
    denom = max(bp,1e-12);
    ED(ib) = sum(prob .* def .* (Rvec ./ denom));
end

figure;
subplot(3,1,1)
plot(bvec, eq.rtilde(ik,:,iz), 'LineWidth',1.5)
title(sprintf('rtilde(b) at ik=%d, iz=%d', ik, iz))

subplot(3,1,2)
plot(bvec, PrND, 'LineWidth',1.5)
title('Pr(non-default) vs b')

subplot(3,1,3)
plot(bvec, rerr, 'LineWidth',1.5); yline(0,'--');
title('Bond pricing error PV - q (should cross 0 once)')

sgtitle('HW diagnostic: pricing curve / root behavior')

% print the worst pricing error
[~, ibmax] = max(abs(rerr));
fprintf('Worst pricing error at b=%.4f: rerr=%.4e, rtilde=%.4f, PrND=%.4f\n', ...
    bvec(ibmax), rerr(ibmax), eq.rtilde(ik,ibmax,iz), PrND(ibmax));

% end


fprintf('--- tail b diagnostics (last 3 b points) ---\n');
for j = max(1,Nb-2):Nb
    bp = grid.bgrid(j);
    rt = eq.rtilde(ik, j, iz);
    prob = Pz(iz,:)';
    kp = grid.kgrid(ik);

    profit_vec = zgrid(:) .* (kp^par.alpha);
    taxbase_vec = profit_vec - par.delta*kp - rt*bp;
    Tc_vec = par.tau_c_pos .* max(taxbase_vec,0) + par.tau_c_neg .* min(taxbase_vec,0);
    wreal_vec = (1-par.delta)*kp + profit_vec - Tc_vec - (1+rt)*bp;

    def = (wreal_vec < eq.wbar(:));
    PrND = sum(prob .* (~def));

    % expected recovery / default payoff term (match update_rtilde!)
    Rvec = zeros(length(zgrid),1);
    for izp=1:length(zgrid)
        Rvec(izp) = hw.recovery_R(kp, zgrid(izp), eq.wbar(izp), par);
    end
    ER = sum(prob .* def .* Rvec);

    fprintf('j=%d bp=%.4g rt=%.4g PrND=%.6f ER=%.6g min(wreal-wbar)=%.6g\n', ...
        j, bp, rt, PrND, ER, min(wreal_vec - eq.wbar(:)));
end
