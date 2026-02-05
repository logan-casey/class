function is_default = default_indicator(w_realized, wbar_zprime)
%HW.DEFAULT_INDICATOR
%
% Returns logical indicator of default.

is_default = (w_realized < wbar_zprime);

end
