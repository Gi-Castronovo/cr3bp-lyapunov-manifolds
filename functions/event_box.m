function [value, isterminal, direction] = event_box(~, Y)
%EVENT_BOX Stops the integration when the trajectory leaves the plot domain
%   (|x| > 2.5 or |y| > 1.5, nondimensional units).

x = Y(1);  y = Y(2);
value      = [x - 2.5;  -2.5 - x;  abs(y) - 1.5];
isterminal = [1; 1; 1];
direction  = [0; 0; 0];
end
