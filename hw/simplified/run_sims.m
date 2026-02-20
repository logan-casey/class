% load or compute eq, par, grid, wgrid, zgrid, Pz (from your solve_equilibrium output)
load('hw_solution_papergrid_feb16_nodamping_b.mat')
stats = hw.simulate_from_eq(eq, par, grid, wgrid, zgrid, Pz, 20000, 200, 14, 1);
save('sims1.mat', 'stats','-mat')
disp('Returned stats:');
disp(stats);