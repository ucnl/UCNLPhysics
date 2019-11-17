% calculates in situ density of water
% millero et al 1980, deep-sea res.,27a,255-264
% jpots ninth report 1978,tenth report 1980
% t - temperature in °C
% p - pressure in mBar
% s - water salinity in PSU
function [water_density] = PHX_WaterDensity_Calc(t, p, s)

    p = p / 1000.0;
    sr = sqrt(abs(s));

    sig = 4.8314E-4 * s;
    sig = sig + ((-1.6546E-6 * t + 1.0227E-4) * t - 5.72466E-3) * sr;
    sig = sig + (((5.3875E-9 * t - 8.2467E-7) * t + 7.6438E-5) * t - 4.0899E-3) * t + 0.824493;
    sig = sig * s;
    sig = sig + ((((6.536332E-9 * t - 1.120083E-6) * t + 1.001685E-4) * t - 9.095290E-3) * t + 6.793952E-2) * t - 0.157406;

    b = ((9.1697E-10 * t + 2.0816E-8) * t - 9.9348E-7) * s + (5.2787E-8 * t - 6.12293E-6) * t + 8.50935E-5;

    k0 = ((-5.3009E-4 * t + 1.6483E-2) * t + 7.944E-2) * sr;
    k0 = k0 + ((-6.1670E-5 * t + 1.09987E-2) * t - 0.603459) * t + 54.6746;
    k0 = k0 * s;
    k0 = k0 + (((-5.155288E-5 * t + 1.360477E-2) * t - 2.327105) * t + 148.4206) * t + 19652.21;

    a = 1.91075E-4 * sr + (-1.6078E-6 * t - 1.0981E-5) * t + 2.2838E-3;
    a = a* s;
    a = a + ((-5.77905E-7 * t + 1.16092E-4) * t + 1.43713E-3) * t + 3.239908;

    k = (b * p + a) * p + k0;

    water_density = 1000.0 + (k * sig + 1000.0 * p) / (k - p);
end