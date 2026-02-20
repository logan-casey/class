function wbar = default_wbar(V, wgrid)
%HW.DEFAULT_WBAR
%
% Extract default threshold wbar(z) from value function.
%
% Inputs
%   V      : value function array [Nw x Nz]
%   wgrid  : net worth grid [Nw x 1]
%
% Output
%   wbar   : [1 x Nz] vector of default thresholds

[Nw, Nz] = size(V);

wbar = zeros(1, Nz);

for iz = 1:Nz

    Vslice = V(:, iz);

    % find first index where value >= 0
    idx = find(Vslice >= 0, 1, 'first');

    if isempty(idx)
        % equity always negative -> choose upper bound
        wbar(iz) = wgrid(end);
        continue;

    elseif idx == 1
        % already nonnegative at lowest point
        wbar(iz) = wgrid(1);
        continue;
   
    else
        % linear interpolation between idx-1 and idx
        w_low  = wgrid(idx-1);
        w_high = wgrid(idx);
    
        V_low  = Vslice(idx-1);
        V_high = Vslice(idx);
    
        % solve V=0 linearly
        weight = -V_low / (V_high - V_low);
    
        wbar(iz) = w_low + weight * (w_high - w_low);
    end

end
end