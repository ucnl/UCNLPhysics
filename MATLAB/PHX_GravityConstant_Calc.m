% calculates gravity at sea level vs latitude
% WGS84 ellipsoid gravity formula
% latitude, signed from -90 to 90
function [gravity_acc] = PHX_GravityConstant_Calc(latitude)
    PHX_GE = 9.7803253359;
    PHX_K  = 0.00193185265241;
    PHX_E  = 0.00669437999013;

    phi_sq = (sin(latitude * pi / 180))^2;    
    gravity_acc = PHX_GE * ((1 + PHX_K * phi_sq) / sqrt(1 - PHX_E * phi_sq));
end