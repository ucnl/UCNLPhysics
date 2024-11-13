using System;
using System.Collections.Generic;

namespace UCNLPhysics
{
    public struct TSProfilePoint
    {
        public double Z;
        public double T;
        public double S;

        public TSProfilePoint(double z, double t, double s)
        {
            Z = z;
            T = t;
            S = s;
        }
    }

    public static class PHX
    {
        public static readonly double PHX_FWTR_DENSITY_KGM3        = 998.02;  // Fresh water density at 20°C
        public static readonly double PHX_FWTR_SOUND_SPEED_MPS     = 1500.0;  // Default speed of sound in water
        public static readonly double PHX_FWTR_SOUND_SPEED_MPS_MIN = 1300.0;  // Min value for speed of sound
        public static readonly double PHX_FWTR_SOUND_SPEED_MPS_MAX = 1600.0;  // Max value for speed of sound
        public static readonly double PHX_FWTR_SALINITY_PSU        = 0.0;     // Default water salinity, PSU
        public static readonly double PHX_SALINITY_PSU_MIN         = 0.0;
        public static readonly double PHX_SALINITY_PSU_MAX         = 44.0;
        public static readonly double PHX_GRAVITY_ACC_MPS2         = 9.80665; // ISO 80000-3:2006
        public static readonly double PHX_GRAVITY_ACC_MPS2_MIN     = 9.7639;
        public static readonly double PHX_GRAVITY_ACC_MPS2_MAX     = 9.8337;

        public static readonly double PHX_ATM_PRESSURE_MBAR        = 1013.25; // Average at sea level
        
        #region Obsolete
        [Obsolete]
        public static double PHX_WaterDensity_Calc(double t, double p, double s)
        {
            return Water_density_calc(t, p, s);
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
            return Speed_of_sound_UNESCO_calc(t, p, s);
        }              
       
        /// <summary>
        /// calculates gravity at sea level vs latitude
        /// WGS84 ellipsoid gravity formula
        /// </summary>
        /// <param name="latitude">latitude, signed from -90 to 90</param>
        /// <returns>Gravity acceleration at sea level, m/s^2</returns>
        [Obsolete]
        public static double PHX_GravityConstant_Calc(double latitude)
        {
            return Gravity_constant_wgs84_calc(latitude);
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
        [Obsolete]
        public static double PHX_DepthByPressure_Calc(double p, double p0, double rho, double g)
        {
            return Depth_by_pressure_calc(p, p0, rho, g);
        }
        #endregion

        // Interpolates a value with given x coordinate by two given points (x1,y1) and (x2,y2)
        public static double Linterp(double x1, double y1, double x2, double y2, double x)
        {
            return y1 + (x - x1) * (y2 - y1) / (x2 - x1);
        }

        /// calculates in situ density of water
        /// millero et al 1980, deep-sea res.,27a,255-264
        /// jpots ninth report 1978,tenth report 1980
        public static double Water_density_calc(double t, double p, double s)
        {
            p = p / 1000.0;
            double sr = Math.Sqrt(Math.Abs(s));
            double sig = (((4.8314E-4 * s) +
                         ((-1.6546E-6 * t + 1.0227E-4) * t - 5.72466E-3) * sr +
                         (((5.3875E-9 * t - 8.2467E-7) * t + 7.6438E-5) * t - 4.0899E-3) * t + 0.824493) * s) +
                         ((((6.536332E-9 * t - 1.120083E-6) * t + 1.001685E-4) * t - 9.095290E-3) * t + 6.793952E-2) * t - 0.157406;

            double b = ((9.1697E-10 * t + 2.0816E-8) * t - 9.9348E-7) * s + (5.2787E-8 * t - 6.12293E-6) * t + 8.50935E-5;

            double k0 = (((((-5.3009E-4 * t + 1.6483E-2) * t + 7.944E-2) * sr) +
                        ((-6.1670E-5 * t + 1.09987E-2) * t - 0.603459) * t + 54.6746) * s) +
                        (((-5.155288E-5 * t + 1.360477E-2) * t - 2.327105) * t + 148.4206) * t + 19652.21;

            double a = (1.91075E-4 * sr + (-1.6078E-6 * t - 1.0981E-5) * t + 2.2838E-3) * s +
                      ((-5.77905E-7 * t + 1.16092E-4) * t + 1.43713E-3) * t + 3.239908;

            double k = (b * p + a) * p + k0;

            return 1000.0 + (k * sig + 1000.0 * p) / (k - p);
        }

        /// The UNESCO equation: Chen and Millero (1977)
        public static double Speed_of_sound_UNESCO_calc(double t, double p, double s)
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

        /// Calculates gravity at sea level vs latitude
        /// WGS84 ellipsoid gravity formula
        /// phi in degrees
        public static double Gravity_constant_wgs84_calc(double phi)
        {
            double phi_sq = Math.Sin(phi * Math.PI / 180.0);
            phi_sq *= phi_sq;
            return (9.7803253359 * ((1.0 + 0.00193185265241 * phi_sq) / Math.Sqrt(1.0 - 0.00669437999013 * phi_sq)));
        }

        /// calculates distance from the water surface where pressure is p0 to the point, where pressure is p
        public static double Depth_by_pressure_calc(double p, double p0, double rho, double g)
        {
            return 100.0 * (p - p0) / (rho * g);                        
        }

        // Calculates pressure of a water column with given height (distance between
        // the water surface and the given point) assuming constant water density
        // h - depth, m
        // p0 - atmospheric pressure, mBar
        // rho - water density, kg/m^3
        // g - gravity acceleration, m/s^2
        public static double Pressure_by_depth_calc(double h, double p0, double rho, double g)
        {
            return h * rho * g / 100.0 + p0;
        }

        // Calculates depth (as a distance between the water surface and a point with
        // the given pressure) by the specified TS-profile
        // pm - pressure measured at the point, mBar
        // p0 - atmospheric pressure, mBar
        // g - gravity acceleration, m/s^2
        // tsProfile - vertical Temperature-Salinity profile at the given point
        // as an array of (f64, f64, f64)
        //   z - vertical coordinate, m (positive, 0 - water surface)
        //   t - temperature, °C
        //   s - salinity, PSU
        // Np - number of pressure intervals for integration
        public static double Depth_by_pressure_ts_profile(double pm, double p0, double g, int n_p, TSProfilePoint[] ts_profile)
        {
            if (n_p <= 0) 
            {
                throw new ArgumentOutOfRangeException("Specified number of time intervals Nt should be greater than zero");
            }

            if (ts_profile.Length < 2)
            {
                throw new ArgumentOutOfRangeException("tsProfile has to contain at least two points");
            }

            double t1 = ts_profile[0].T;
            double s1 = ts_profile[0].S;
            double rho0 = Water_density_calc(t1, p0, s1);
            double p1 = Pressure_by_depth_calc(ts_profile[0].Z, p0, rho0, g);
            double pe = Pressure_by_depth_calc(ts_profile[ts_profile.Length - 1].Z, p0, rho0, g);

            if ((pm < p1) || (pm > pe))
            {
                throw new ArgumentOutOfRangeException("Specified pressure is beyond the specified TS-profile");
            }
            
            int p_idx = 1;
            double t2 = ts_profile[p_idx].T;
            double s2 = ts_profile[p_idx].S;
            double p2 = Pressure_by_depth_calc(ts_profile[p_idx].Z, p0, rho0, g);
             
            double dp = (pm - p0) / n_p;
            double h = 0.0, rho, t, p = p0, s;
    
            while (p < pm) 
            {
                p += dp;

                if (p > p2)
                {

                    p1 = p2;
                    t1 = t2;
                    s1 = s2;
                    p_idx += 1;

                    t2 = ts_profile[p_idx].T;
                    s2 = ts_profile[p_idx].S;
                    p2 = Pressure_by_depth_calc(ts_profile[p_idx].Z, p0, rho0, g);
                }

                t = Linterp(p1, t1, p2, t2, p);
                s = Linterp(p1, s1, p2, s2, p);

                rho = Water_density_calc(t, p, s);
                h += 1.0 / rho;
            }

            return h * 100.0 * dp / g;
        }

        // Calculates the path, which sound traveled in vertical direction
        // between the water surface and the deepest point during
        // a given time of flight considering given temperature and salinity profile
        // tof - time of flight, sec
        // Nt - number of time intervals for integration
        // tsProfile - vertical Temperature-Salinity profile at the given point
        // as an array of TSPoint
        //   z - vertical coordinate, m (positive, 0 - water surface)
        //   t - temperature, °C
        //   s - salinity, PSU
        public static double Vertical_sound_path_ts_profile(double tof, double g, int n_t, TSProfilePoint[] ts_profile)
        {

	        if (ts_profile.Length < 2)
            {
                throw new ArgumentOutOfRangeException("tsProfile has to contain at least two points");
            }
  
            if (n_t <= 0)
            {
                throw new ArgumentOutOfRangeException("Specified number of time intervals Nt should be greater than zero");
            }
        
            double z1 = ts_profile[0].Z;
            double t1 = ts_profile[0].T;
            double s1 = ts_profile[0].S;
            double rho0 = Water_density_calc(t1, PHX_ATM_PRESSURE_MBAR, s1);
            double p1 = Pressure_by_depth_calc(z1, PHX_ATM_PRESSURE_MBAR, rho0, g);
  
            double v = Speed_of_sound_UNESCO_calc(t1, p1, s1);
  
            if (v * tof > ts_profile[ts_profile.Length - 1].Z)
            {
                throw new ArgumentOutOfRangeException("Specified time of flight is beyond the specified TS-profile");
            }
  
            int p_idx = 1;
            double z2 = ts_profile[p_idx].Z;
            double t2 = ts_profile[p_idx].T;
            double s2 = ts_profile[p_idx].S;
            double p2 = Pressure_by_depth_calc(z2, PHX_ATM_PRESSURE_MBAR, rho0, g);
  
            double dt = tof / n_t;
            double h = 0.0;
            double t;
            double p;
            double s;
            double tt = 0.0;
  
            while (tt < tof)
            {
  
                tt += dt;
                h = h + dt * v;
  
                if (h > z2)
                {  
                    p1 = p2;
                    t1 = t2;
                    s1 = s2;
                    z1 = z2;
                    p_idx = p_idx + 1;
  
                    z2 = ts_profile[p_idx].Z;
                    t2 = ts_profile[p_idx].T;
                    s2 = ts_profile[p_idx].S;
                    p2 = Pressure_by_depth_calc(z2, PHX_ATM_PRESSURE_MBAR + p1, rho0, g);
                }
  
                t = Linterp(z1, t1, z2, t2, h);
                p = Linterp(z1, p1, z2, p2, h);
                s = Linterp(z1, s1, z2, s2, h);
                v = Speed_of_sound_UNESCO_calc(t, p, s);
            }
  
            return h;
        }

        // Calculated the freezing temperature of seawater (in °C) with specified pressure and salinity.
        // According to:
        // Algorithms for computation of fundamental properties of seawater.
        // Unesco technical papers in marine science vol. 44, 1983, pp. 30
        // https://darchive.mblwhoilibrary.org/bitstream/handle/1912/2470/059832eb.pdf
        // p - pressure, mBar
        // s - PSU
        public static double Water_fpoint_calc(double p, double s)
        {
           return (-0.0575 + 1.710523E-3 * Math.Sqrt(s) - 2.154996E-4 * s) * s - 7.53E-6 * p;
        }

        // calculation of absorption according to:
        // Francois & Garrison, J. Acoust. Soc. Am., Vol. 72, No. 6, December 1982
        // f frequency (kHz)
        // T Temperature (degC)
        // S Salinity (ppt)
        // D Depth (m)
        // pH Acidity
        public static double Alpha_e_FrancoisGarrison_calc(double f, double t, double s, double h, double pH)
        {
            // For f = 1..500 kHz:
            // -2 < T < 22 °C
            // 30 < S < 35 PSU
            // 0 < D < 3.5 km

            // For f > 500 kHz:

            // 0 < T < 30 °C
            // 0 < S < 40 PSU
            // 0 < D < 10 km            

            // Total absorption = Boric Acid Contrib. + Magnesium Sulphate Contrib. + Pure Water Contrib.

            // Measured ambient temp
            double t_kel = 273.15 + t;
            double fsq = f * f;            

            // Calculate speed of sound (according to Francois & Garrison, JASA 72 (6) p1886)
            double c = 1412 + 3.21 * t + 1.19 * s + 0.0167 * h;

            // Boric acid contribution
            double A1 = (8.86 / c) * Math.Pow(10, 0.78 * pH - 5.0);
            double P1 = 1;
            double f1 = 2.8 * Math.Sqrt(s / 35) * Math.Pow(10, 4.0 - 1245 / t_kel);
            double Boric = (A1 * P1 * f1 * fsq) / (fsq + f1 * f1);

            // MgSO4 contribution
            double A2 = 21.44 * (s / c) * (1 + 0.025 * t);
            double P2 = 1 - 1.37E-4 * h + 6.2E-9 * h * h;
            double f2 = (8.17 * Math.Pow(10, 8 - 1990 / t_kel)) / (1 + 0.0018 * (s - 35));
            double MgSO4 = (A2 * P2 * f2 * fsq) / (fsq + f2 * f2);

            // Pure water contribution
            double A3;
            if (t <= 20)
            {
                A3 = 4.937E-4 - 2.59E-5 * t + 9.11E-7 * t * t - 1.5E-8 * t * t * t;
            }
            else
            {
                A3 = 3.964E-4 - 1.146E-5 * t + 1.45E-7 * t * t - 6.5E-10 * t * t * t;
            }

            double P3 = 1 - 3.83E-5 * h + 4.9E-10 * h * h;
            double H2O = A3 * P3 * fsq;

            // Total absorption
            return Boric + MgSO4 + H2O;
        }

    }
}
