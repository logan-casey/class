function wtilde = revised_networth(w, wbar)
%HW.REVISED_NETWORTH
%
% Revised net worth:
%   wtilde = max( wbar(z'), w )
%
% Inputs
%   w     : realized net worth (scalar or array)
%   wbar  : default threshold corresponding to same z'
%
% Output
%   wtilde : revised net worth

wtilde = max(wbar, w);

end