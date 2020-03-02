using System;
using System.Collections.Generic;

namespace UCNLPhysics
{
    public struct TSProfilePoint
    {
        public double Z;
        public double Pa;
        public double T;
        public double S;
    }

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
   
        private static double Linterp(double x1, double y1, double x2, double y2, double x)
        {
            return y1 + (x - x1)*(y2 - y1)/(x2 - x1);
        }

        public static double PHX_Approx_Pressure_By_Depth(double h, double p0, double rho, double g)
        {
            return h * rho * g / 100.0 + p0;
        }

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
            p = p / 1000.0;
            double sr = Math.Sqrt(Math.Abs(s));

            double sig = (4.8314E-4 * s +
                       ((-1.6546E-6 * t + 1.0227E-4) * t - 5.72466E-3) * sr +
                       (((5.3875E-9 * t - 8.2467E-7) * t + 7.6438E-5) * t - 4.0899E-3) * t + 0.824493) * s +
                      ((((6.536332E-9 * t - 1.120083E-6) * t + 1.001685E-4) * t - 9.095290E-3) * t + 6.793952E-2) * t - 0.157406;

            double b = ((9.1697E-10 * t + 2.0816E-8) * t - 9.9348E-7) * s + (5.2787E-8 * t - 6.12293E-6) * t + 8.50935E-5;

            double k0 = (((-5.3009E-4 * t + 1.6483E-2) * t + 7.944E-2) * sr +
                      ((-6.1670E-5 * t + 1.09987E-2) * t - 0.603459) * t + 54.6746) * s +
                     (((-5.155288E-5 * t + 1.360477E-2) * t - 2.327105) * t + 148.4206) * t + 19652.21;

            double a = (1.91075E-4 * sr + (-1.6078E-6 * t - 1.0981E-5) * t + 2.2838E-3) * s +
                    ((-5.77905E-7 * t + 1.16092E-4) * t + 1.43713E-3) * t + 3.239908;

            double k = (b * p + a) * p + k0;
            return 1000.0 + (k * sig + 1000.0 * p) / (k - p);
        }

        /// <summary>
        /// The UNESCO equation: Chen and Millero (1977)
        /// </summary>
        /// <param name="t">temperature, Celsius degree</param>
        /// <param name="p">pressure, mBar</param>
        /// <param name="s">salinity, PSU</param>
        /// <returns>Speed of sound in m/s</returns>
        [Obsolete]
        public static double PHX_SpeedOfSound_Calc(double t, double p, double s)
        {
            return PHX_SpeedOfSoundInWater_UNESCO(t, p, s);
        }

        /// <summary>
        /// The UNESCO equation: Chen and Millero (1977)
        /// </summary>
        /// <param name="t">temperature, Celsius degree</param>
        /// <param name="p">pressure, mBar</param>
        /// <param name="s">salinity, PSU</param>
        /// <returns>Speed of sound in m/s</returns>
        public static double PHX_SpeedOfSoundInWater_UNESCO(double t, double p, double s)
        {
            var t2 = t * t;
            var t3 = t2 * t;
            var t4 = t3 * t;
            p = p / 1000.0;
            var p2 = p * p;
            var p3 = p2 * p;

            var Cw = (1402.388 + 5.03830 * t + -5.81090E-2 * t2 + 3.3432E-4 * t3 + -1.47797E-6 * t4 + 3.1419E-9 * t4 * t) +
                       (0.153563 + 6.8999E-4 * t + -8.1829E-6 * t2 + 1.3632E-7 * t3 + -6.1260E-10 * t4) * p +
                       (3.1260E-5 + -1.7111E-6 * t + 2.5986E-8 * t2 + -2.5353E-10 * t3 + 1.0415E-12 * t4) * p2 +
                       (-9.7729E-9 + 3.8513E-10 * t + -2.3654E-12 * t2) * p3;

            var A = (1.389 + -1.262E-2 * t + 7.166E-5 * t2 + 2.008E-6 * t3 + -3.21E-8 * t4) +
                      (9.4742E-5 + -1.2583E-5 * t + -6.4928E-8 * t2 + 1.0515E-8 * t3 + -2.0142E-10 * t4) * p +
                      (-3.9064E-7 + 9.1061E-9 * t + -1.6009E-10 * t2 + 7.994E-12 * t3) * p2 +
                      (1.100E-10 + 6.651E-12 * t + -3.391E-13 * t2) * p3;

            var B = -1.922E-2 + -4.42E-5 * t + (7.3637E-5 + 1.7950E-7 * t) * p;

            var D = 1.727E-3 + -7.9836E-6 * p;

            return Cw + A * s + B * Math.Sqrt(s * s * s) + D * s * s;
        }

        /// <summary>
        /// Calculate speed of sound (according to Francois & Garrison, JASA 72 (6) p1886) 
        /// </summary>
        /// <param name="t">temperature in °C</param>
        /// <param name="d">depth in meters</param>
        /// <param name="s">salinity in PSU</param>
        /// <returns></returns>
        public static double PHX_SpeedOfSoundInWater_FrancoisGarrison(double t, double d, double s)
        {
           return 1412 + 3.21 * t + 1.19 * s + 0.0167 * d;
        }

        /// <summary>
        /// calculates gravity at sea level vs latitude
        /// WGS84 ellipsoid gravity formula
        /// </summary>
        /// <param name="latitude">latitude, signed from -90 to 90</param>
        /// <returns>Gravity acceleration at sea level, m/s^2</returns>
        public static double PHX_GravityConstant_Calc(double latitude)
        {
            double phi_sq = Math.Sin(latitude * Math.PI / 180.0);
            phi_sq *= phi_sq;
            return PHX_GE * ((1 + PHX_K * phi_sq) / Math.Sqrt(1 - PHX_E * phi_sq));
        }

        /// <summary>
        /// Calculates distance from the water surface where pressure is p0 to the point, where pressure is p
        /// To take into account compression of the water column the better way to use water density
        /// estimated for the point with Pm = (P-P0)/2
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

        /// <summary>
        /// Calculated vertical distance from the water surface with atmispheric pressure p0
        /// by measured pressure p considering compression of water column with constant temperature t
        /// and salinity s.
        /// </summary>
        /// <param name="p">Measured pressure, mBar</param>
        /// <param name="p0">Measured atmospheric pressure, mBar</param>
        /// <param name="t">Water temperature, °C</param>
        /// <param name="s">Salinity, PSU</param>
        /// <param name="g">Gravity acceleration, m/s^2</param>
        /// <param name="Np">Number of pressure intervals for integration</param>
        /// <returns>depth (distance from water surface)</returns>
        public static double PHX_DepthByPressure_CompressionTSConst_Calc(double p, double p0, double t, double s, double g, int Np)
        {
            double dp = (p - p0) / Np;
            double he = 0;

            for (int n = 0; n < Np; n++)
            {
                he += 1.0/PHX_WaterDensity_Calc(t, p0 + n * dp, s);
            }
            
            return 100.0 * he * dp / g;
        }

        public static double PHX_DepthBy_Pressure_TSProfile_Calc(double pm, double p0, double g, int Np, 
            TSProfilePoint[] tsProfile)
        {
            if ((tsProfile == null) && (tsProfile.Length < 1))
                throw new ArgumentException("tsProfile should have more than 1 elements");

            if ((pm >= tsProfile[0].Pa) && (pm <= tsProfile[tsProfile.Length - 1].Pa))
            {
                double p1 = tsProfile[0].Pa;
                double t1 = tsProfile[0].T;
                double s1 = tsProfile[0].S;
                double p2 = tsProfile[1].Pa;
                double t2 = tsProfile[1].T;
                double s2 = tsProfile[1].S;

                double dp = (pm - p0) / Np;
                int pointIdx = 1;
                double h = 0;
                double t, p, s;

                for (int n = 0; n < Np; n++)
                {
                    p = p0 + n * dp;
                    
                    if (p > p2) 
                    {
                         p1 = p2;
                         t1 = t2;
                         s1 = s2;
                         pointIdx++;
                         p2 = tsProfile[pointIdx].Pa;
                         t2 = tsProfile[pointIdx].T;
                         s2 = tsProfile[pointIdx].S;
                    }
                    t = Linterp(p1, t1, p2, t2, p);
                    s = Linterp(p1, s1, p2, s2, p);                    
                    h += 1.0 / PHX_WaterDensity_Calc(t, p, s);
                }

                h *= 100.0 * dp / g;
                return h;
            }
            else
            {
                throw new ArgumentOutOfRangeException("pm is beyond the specified tsProfile");
            }
        }

        /// <summary>
        /// Estimates vertical distance traveled by sound during time tof considering specified temperature and salinity profile.
        /// By dividing whole time in Nt intervals and integrating the distance with layers of water with different density.
        /// </summary>
        /// <param name="tof">A half of traveling time (one way traveling time), sec</param>
        /// <param name="g">gravity acceleration, m/s^2</param>
        /// <param name="Nt">Number of intervals to time</param>
        /// <param name="tsProfile">temperature and salinity profile</param>
        /// <returns>Depth as a traveled distance through medium with layers with different temperature, pressure and salinity, m</returns>
        public static double PHX_VerticalSoundPath_By_Time_TSProfile_Calc(double tof, double g, int Nt, TSProfilePoint[] tsProfile)
        {
            if ((tsProfile == null) && (tsProfile.Length < 1))
                throw new ArgumentException("tsProfile should have more than 1 elements");

            double approxDist = PHX_FWTR_SOUND_SPEED_MPS * tof;
            
            if ((approxDist >= tsProfile[0].Z) && (approxDist <= tsProfile[tsProfile.Length - 1].Z))
            {
                double p1 = tsProfile[0].Pa;
                double t1 = tsProfile[0].T;
                double s1 = tsProfile[0].S;
                double z1 = tsProfile[0].Z;
                double p2 = tsProfile[1].Pa;
                double t2 = tsProfile[1].T;
                double s2 = tsProfile[1].S;
                double z2 = tsProfile[1].Z;

                double dt = tof / Nt;
                int pointIdx = 1;
                double h = 0;                
                double v, t, p, s;

                v = PHX_SpeedOfSoundInWater_UNESCO(t1, p1, s1);

                for (int n = 0; n < Nt; n++)
                {
                    h += dt * v;

                    if (h > z2)
                    {
                        p1 = p2;
                        t1 = t2;
                        s1 = s2;
                        z1 = z2;
                        pointIdx++;
                        p2 = tsProfile[pointIdx].Pa;
                        t2 = tsProfile[pointIdx].T;
                        s2 = tsProfile[pointIdx].S;
                        z2 = tsProfile[pointIdx].Z;
                    }

                    t = Linterp(z1, t1, z2, t2, h);
                    p = Linterp(z1, p1, z2, p2, h);
                    s = Linterp(z1, s1, z2, s2, h);
                    v = PHX_SpeedOfSoundInWater_UNESCO(t, p, s);
                }

                return h;
            }
            else
            {
                throw new ArgumentOutOfRangeException("tof is beyond the specified tsProfile");
            }
        }
    }
}
