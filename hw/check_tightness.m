ik = round(rand(1)*par.Nk)+1;
ib = round(rand(1)*par.Nz_solve)+1;

kp = grid.kgrid(eq.pol_ik(ik,ib)); bp = grid.bgrid(eq.pol_ib(ik,ib)); rt = 0.025;  % use whatever you are testing
zvec = zgrid(:);

profit_vec = zvec .* (kp^par.alpha);

taxbase_vec = profit_vec - par.delta*kp - rt*bp;
Tc_vec = par.tau_c_pos.*max(taxbase_vec,0) + par.tau_c_neg.*min(taxbase_vec,0);

wreal_vec = (1-par.delta)*kp + profit_vec - Tc_vec - (1+rt)*bp;

fprintf("kp=%.4f, kp^alpha=%.4f\n", kp, kp^par.alpha);
fprintf("z range = %.4f to %.4f (range %.4f)\n", min(zvec), max(zvec), max(zvec)-min(zvec));
fprintf("profit range = %.4f\n", max(profit_vec)-min(profit_vec));
fprintf("wreal range  = %.4f\n", max(wreal_vec)-min(wreal_vec));
fprintf("Tc range     = %.4f\n", max(Tc_vec)-min(Tc_vec));
