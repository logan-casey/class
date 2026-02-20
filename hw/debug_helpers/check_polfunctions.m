[x, y] = meshgrid(zgrid, wgrid);
zb = zeros(size(eq.pol_ib));
for i = 1:size(eq.pol_ib,1)
    for j = 1:size(eq.pol_ib,2)
        zb(i,j) = grid.bgrid(eq.pol_ib(i,j));
    end
end
figure
surf(x,y,zb);
saveas(gcf, 'pol_ib.png','png')
zk = zeros(size(eq.pol_ik));
for i = 1:size(eq.pol_ik,1)
    for j = 1:size(eq.pol_ik,2)
        zk(i,j) = grid.kgrid(eq.pol_ik(i,j));
    end
end
figure
surf(x,y,zk);
saveas(gcf, 'pol_ik.png','png')

figure
surf(x,y,eq.V);
saveas(gcf, 'val.png','png')

