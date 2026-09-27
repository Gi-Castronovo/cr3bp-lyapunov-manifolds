function [value, isterminal, direction] = event_vertical_section(~, Y, x_section, mu)
%EVENT_VERTICAL_SECTION Stops the integration on the Poincare section x = x_section.
%   Additional terminal events act as safety limits: the trajectory leaves
%   the Moon region (|y| > 0.5 or x outside [-2.5, 2.4]) or comes closer
%   than 0.005 (about 1900 km) to the Moon, where the point-mass model is
%   no longer meaningful.

x = Y(1);  y = Y(2);
dist_moon = sqrt((x - (1 - mu))^2 + y^2);

value = [x - x_section;          % 1: Poincare section (target event)
         abs(y) - 0.5;           % 2: vertical limit
         x - 2.4;                % 3: right limit
         -2.5 - x;               % 4: left limit
         dist_moon - 0.005];     % 5: Moon proximity
isterminal = [1; 1; 1; 1; 1];
direction  = [0; 0; 0; 0; 0];
end
