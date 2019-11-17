% calculates distance from the water surface where pressure is p0 to the point, where pressure is p
% pressure, mBar
% pressure at water surface, mBar
% water density, kg/m^3
% gravity acceleration at sea level, m/s^2
function [depth_m] = PHX_DepthByPressure_Calc(p, p0, rho, g)
   depth_m = (100.0 * (p - p0) / (rho * g));
end
