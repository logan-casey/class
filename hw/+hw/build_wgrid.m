function wgrid = build_wgrid(wmin, wmax, Nw)
%HW.BUILD_WGRID  Simple linear grid for revised net worth.
if Nw < 2
    error("Nw must be >= 2");
end
wgrid = linspace(wmin, wmax, Nw)';
end
