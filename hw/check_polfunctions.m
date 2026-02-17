[x, y] = meshgrid(zgrid, wgrid);
zb = zeros(size(eq.pol_ib));
for i = 1:size(eq.pol_ib,1)
    for j = 1:size(eq.pol_ib,2)
        zb(i,j) = grid.bgrid(eq.pol_ib(i,j));
    end
end
figure
surf(x,y,zb);
zk = zeros(size(eq.pol_ik));
for i = 1:size(eq.pol_ik,1)
    for j = 1:size(eq.pol_ik,2)
        zk(i,j) = grid.kgrid(eq.pol_ik(i,j));
    end
end
figure
surf(x,y,zk);

figure
surf(x,y,eq.V);

