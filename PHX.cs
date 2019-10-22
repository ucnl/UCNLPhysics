using System;

namespace UCNLPhysics
{
    public static class PHX
    {
        public static readonly double PHX_FWTR_DENSITY_KGM3        = 998.02;  // Fresh water density at 20°C
        public static readonly double PHX_FWTR_SOUND_SPEED_MPS     = 1500.0;  // Default speed of sound in water
        public static readonly double PHX_FWTR_SOUND_SPEED_MPS_MIN = 1300.0;  // Min value for speed of sound
        public static readonly double PHX_FWTR_SOUND_SPEED_MPS_MAX = 1600.0;  // Max value for speed of sound
        public static readonly double PHX_FWTR_SALINITY_PSU        = 0.0;     // Default water salinity, PSU
        public static readonly double PHX_GRAVITY_ACC_MPS2         = 9.80665; // ISO 80000-3:2006
        public static readonly double PHX_ATM_PRESSURE_MBAR        = 1013.25; // Average at sea level

        static readonly double PHX_GE = 9.7803253359;
        static readonly double PHX_K  = 0.00193185265241;
        static readonly double PHX_E  = 0.00669437999013;

        static readonly double C00 = 1402.388;
        static readonly double C01 = 5.03830;
        static readonly double C02 = -5.81090E-2;
        static readonly double C03 = 3.3432E-4;
        static readonly double C04 = -1.47797E-6;
        static readonly double C05 = 3.1419E-9;
        static readonly double C10 = 0.153563;
        static readonly double C11 = 6.8999E-4;
        static readonly double C12 = -8.1829E-6;
        static readonly double C13 = 1.3632E-7;
        static readonly double C14 = -6.1260E-10;
        static readonly double C20 = 3.1260E-5;
        static readonly double C21 = -1.7111E-6;
        static readonly double C22 = 2.5986E-8;
        static readonly double C23 = -2.5353E-10;
        static readonly double C24 = 1.0415E-12;
        static readonly double C30 = -9.7729E-9;
        static readonly double C31 = 3.8513E-10;
        static readonly double C32 = -2.3654E-12;
        static readonly double A00 = 1.389;
        static readonly double A01 = -1.262E-2;
        static readonly double A02 = 7.166E-5;
        static readonly double A03 = 2.008E-6;
        static readonly double A04 = -3.21E-8;
        static readonly double A10 = 9.4742E-5;
        static readonly double A11 = -1.2583E-5;
        static readonly double A12 = -6.4928E-8;
        static readonly double A13 = 1.0515E-8;
        static readonly double A14 = -2.0142E-10;
        static readonly double A20 = -3.9064E-7;
        static readonly double A21 = 9.1061E-9;
        static readonly double A22 = -1.6009E-10;
        static readonly double A23 = 7.994E-12;
        static readonly double A30 = 1.100E-10;
        static readonly double A31 = 6.651E-12;
        static readonly double A32 = -3.391E-13;
        static readonly double B00 = -1.922E-2;
        static readonly double B01 = -4.42E-5;
        static readonly double B10 = 7.3637E-5;
        static readonly double B11 = 1.7950E-7;
        static readonly double D00 = 1.727E-3;
        static readonly double D10 = -7.9836E-6;


        /// <summary>
        /// calculates in situ density of water
        /// millero et al 1980, deep-sea res.,27a,255-264
        /// jpots ninth report 1978,tenth report 1980
        /// </summary>
        /// <param name="t">temperature, Celsius degree</param>
        /// <param name="p">pressure, mBar</param>
        /// <param name="s">Salinity, PSU</param>
        /// <returns>water density in kg/m^3</returns>
        public static double PHX_WaterDensity_Calc(double t, double p, double s)
        {
            double sig, k0, b, a, k, sr;
            p = p / 1000.0;
            sr = Math.Sqrt(Math.Abs(s));
            sig = 4.8314E-4 * s;
            sig += ((-1.6546E-6 * t + 1.0227E-4) * t - 5.72466E-3) * sr;
            sig += (((5.3875E-9 * t - 8.2467E-7) * t + 7.6438E-5) * t - 4.0899E-3) * t + 0.824493;
            sig *= s;
            sig += ((((6.536332E-9 * t - 1.120083E-6) * t + 1.001685E-4) * t - 9.095290E-3) * t + 6.793952E-2) * t - 0.157406;
            b = ((9.1697E-10 * t + 2.0816E-8) * t - 9.9348E-7) * s + (5.2787E-8 * t - 6.12293E-6) * t + 8.50935E-5;
            k0 = ((-5.3009E-4 * t + 1.6483E-2) * t + 7.944E-2) * sr;
            k0 += ((-6.1670E-5 * t + 1.09987E-2) * t - 0.603459) * t + 54.6746;
            k0 *= s;
            k0 += (((-5.155288E-5 * t + 1.360477E-2) * t - 2.327105) * t + 148.4206) * t + 19652.21;
            a = 1.91075E-4 * sr + (-1.6078E-6 * t - 1.0981E-5) * t + 2.2838E-3;
            a *= s;
            a += ((-5.77905E-7 * t + 1.16092E-4) * t + 1.43713E-3) * t + 3.239908;
            k = (b * p + a) * p + k0;
            sig = 1000.0 + (k * sig + 1000.0 * p) / (k - p);

            return sig;
        }

        /// <summary>
        /// The UNESCO equation: Chen and Millero (1977)
        /// </summary>
        /// <param name="t">temperature, Celsius degree</param>
        /// <param name="p_mBar">pressure, mBar</param>
        /// <param name="s">salinity, PSU</param>
        /// <returns>Speed of sound in m/s</returns>
        public static double PHX_SpeedOfSound_Calc(double t, double p_mBar, double s)
        {
            double t2 = t * t;
            double t3 = t2 * t;
            double t4 = t3 * t;
            double p = p_mBar / 1000.0f;

            double Cw = (C00 + C01 * t + C02 * t2 + C03 * t3 + C04 * t4 + C05 * t4 * t) +
                       (C10 + C11 * t + C12 * t2 + C13 * t3 + C14 * t4) * p +
                       (C20 + C21 * t + C22 * t2 + C23 * t3 + C24 * t4) * p * p +
                       (C30 + C31 * t + C32 * t2) * p * p * p;

            double A = (A00 + A01 * t + A02 * t2 + A03 * t3 + A04 * t4) +
                      (A10 + A11 * t + A12 * t2 + A13 * t3 + A14 * t4) * p +
                      (A20 + A21 * t + A22 * t2 + A23 * t3) * p * p +
                      (A30 + A31 * t + A32 * t2) * p * p * p;

            double B = B00 + B01 * t + (B10 + B11 * t) * p;

            double D = D00 + D10 * p;

            double c = Cw + A * s + B * Math.Sqrt(s * s * s) + D * s * s;

            return c;
        }

        /// <summary>
        /// calculates gravity at sea level vs latitude
        /// WGS84 ellipsoid gravity formula
        /// </summary>
        /// <param name="latitude">latitude, signed from -90 to 90</param>
        /// <returns>Gravity acceleration at sea level, m/s^2</returns>
        public static double PHX_GravityConstant_Calc(double latitude)
        {
            double phi_sq = Math.Sin(latitude);
            phi_sq *= phi_sq;
            return PHX_GE * ((1 + PHX_K * phi_sq) / Math.Sqrt(1 - PHX_E * phi_sq));
        }

        /// <summary>
        /// calculates distance from the water surface where pressure is p0 to the point, where pressure is p
        /// </summary>
        /// <param name="p">pressure, mBar</param>
        /// <param name="p0">pressure at water surface, mBar</param>
        /// <param name="rho">water density, kg/m^3</param>
        /// <param name="g">gravity acceleration at sea level, m/s^2</param>
        /// <returns>depth (distance from water surface)</returns>
        public static double PHX_DepthByPressure_Calc(double p, double p0, double rho, double g)
        {
            return (100.0 * (p - p0) / (rho * g));
        }
    }
}
