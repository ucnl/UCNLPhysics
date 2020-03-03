// (C) Aleksandr Dikarev, 2016-2020
// physx.c

#include <core_physx.h>

// Interpolates a value with given x coordinate by two given points (x1,y1) and (x2,y2)
float PHX_linterp(float x1, float y1, float x2, float y2, float x)
{
    return y1 + (x - x1) * (y2 - y1) / (x2 - x1);
}

/// calculates in situ density of water
/// millero et al 1980, deep-sea res.,27a,255-264
/// jpots ninth report 1978,tenth report 1980
float PHX_water_density_calc(float t, float p, float s)
{
    p = p / 1000.0f;
    float sr = _sqrt_f32(s);
    float sig = (((4.8314E-4f * s) +
                 ((-1.6546E-6f * t + 1.0227E-4f) * t - 5.72466E-3f) * sr +
                 (((5.3875E-9f * t - 8.2467E-7f) * t + 7.6438E-5f) * t - 4.0899E-3f) * t + 0.824493f) * s) +
                 ((((6.536332E-9f * t - 1.120083E-6f) * t + 1.001685E-4f) * t - 9.095290E-3f) * t + 6.793952E-2f) * t - 0.157406f;

    float b = ((9.1697E-10f * t + 2.0816E-8f) * t - 9.9348E-7f) * s + (5.2787E-8f * t - 6.12293E-6f) * t + 8.50935E-5f;

    float k0 = (((((-5.3009E-4f * t + 1.6483E-2f) * t + 7.944E-2f) * sr) +
                ((-6.1670E-5f * t + 1.09987E-2f) * t - 0.603459f) * t + 54.6746f) * s) +
                (((-5.155288E-5f * t + 1.360477E-2f) * t - 2.327105f) * t + 148.4206f) * t + 19652.21f;

    float a = (1.91075E-4f * sr + (-1.6078E-6f * t - 1.0981E-5f) * t + 2.2838E-3f) * s +
              ((-5.77905E-7f * t + 1.16092E-4f) * t + 1.43713E-3f) * t + 3.239908f;

    float k = (b * p + a) * p + k0;

    return 1000.0f + (k * sig + 1000.0 * p) / (k - p);
}

/// The UNESCO equation: Chen and Millero (1977)
float PHX_speed_of_sound_UNESCO_calc(float t, float p, float s)
{
    p = p / 1000.0;
    float sr = _sqrt_f32(s);

    float d = 1.727E-3f - 7.9836E-6f * p;

    float b_1 = 7.3637E-5f + 1.7945E-7f * t;
    float b_0 = -1.922E-2f - 4.42E-5f * t;
    float b = b_0 + b_1 * p;

    float a_3 = (-3.389E-13f * t + 6.649E-12f)  * t + 1.100E-10f;
    float a_2 = ((7.988E-12f * t - 1.6002E-10f) * t + 9.1041E-9f) * t - 3.9064E-7f;
    float a_1 = (((-2.0122E-10f * t + 1.0507E-8f)  * t - 6.4885E-8f) * t - 1.2580E-5f) * t + 9.4742E-5f;
    float a_0 = (((-3.21E-8f * t + 2.006E-6f) * t + 7.164E-5f) * t - 1.262E-2f) * t + 1.389f;
    float a = ((a_3 * p + a_2) * p + a_1) * p + a_0;

    float c_3 = (-2.3643E-12f * t + 3.8504E-10f) * t - 9.7729E-9f;
    float c_2 = (((1.0405E-12f * t - 2.5335E-10f) * t + 2.5974E-8f) * t - 1.7107E-6f)  * t + 3.1260E-5f;
    float c_1 = (((-6.1185E-10f * t + 1.3621E-7f)  * t - 8.1788E-6f) * t + 6.8982E-4f)  * t + 0.153563f;
    float c_0 = ((((3.1464E-9f  * t - 1.47800E-6f) * t + 3.3420E-4f) * t - 5.80852E-2f) * t + 5.03711f) * t + 1402.388f;
    float c  = ((c_3 * p + c_2) * p + c_1) * p + c_0;

    return c + (a + b * sr + d * s) * s;
}

/// Calculates gravity at sea level vs latitude
/// WGS84 ellipsoid gravity formula
float PHX_gravity_constant_wgs84_calc(float phi)
{
    float phi_sq = _sin_f32(phi);
    phi_sq *= phi_sq;
    return (9.7803253359f * ((1.0f + 0.00193185265241f * phi_sq) / _sqrt_f32(1.0f - 0.00669437999013f * phi_sq)));
}

/// calculates distance from the water surface where pressure is p0 to the point, where pressure is p
float PHX_depth_by_pressure_calc(float p, float p0, float rho, float g)
{
    return 100.0f * (p - p0) / (rho * g);
}

// Calculates pressure of a water column with given height (distance between
// the water surface and the given point) assuming constant water density
// h - depth, m
// p0 - atmospheric pressure, mBar
// rho - water density, kg/m^3
// g - gravity acceleration, m/s^2
float PHX_pressure_by_depth_calc(float h, float p0, float rho, float g)
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
bool PHX_depth_by_pressure_ts_profile(float pm, float p0, float g, int n_p, TSProfileItem_Struct* ts_profile, int ts_profile_size,  float* dpt)
{
    if ((n_p <= 0) ||
    	(pm < p0) || (p0 < 0) ||
    	(ts_profile_size < 2))
    {
    	return false;
    }

    float t1 = ts_profile[0].t;
    float s1 = ts_profile[0].s;
    float rho0 = PHX_water_density_calc(t1, p0, s1);
    float p1 = PHX_pressure_by_depth_calc(ts_profile[0].z, p0, rho0, g);

    float pe = PHX_pressure_by_depth_calc(ts_profile[ts_profile_size - 1].z, p0, rho0, g);

    if ((pm < p1) || (pm > pe))
    {
        return false;
    }

    int p_idx = 1;
    float t2 = ts_profile[p_idx].t;
    float s2 = ts_profile[p_idx].s;
    float p2 = PHX_pressure_by_depth_calc(ts_profile[p_idx].z, p0, rho0, g);

    float dp = (pm - p0) / n_p;
    *dpt = 0.0f;

    float rho;
    float t;
    float p = p0;
    float s;

    while (p < pm) {

        p += dp;

        if (p > p2) {

            p1 = p2;
            t1 = t2;
            s1 = s2;
            p_idx++;

            t2 = ts_profile[p_idx].t;
            s2 = ts_profile[p_idx].s;
            p2 = PHX_pressure_by_depth_calc(ts_profile[p_idx].z, p0, rho0, g);
        }

        t = PHX_linterp(p1, t1, p2, t2, p);
        s = PHX_linterp(p1, s1, p2, s2, p);

        rho = PHX_water_density_calc(t, p, s);
        *dpt += 1.0f / rho;
    }

    *dpt = (*dpt * 100.0f * dp / g);
    return true;
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
bool PHX_vertical_sound_path_ts_profile(float tof, float g, int n_t, TSProfileItem_Struct* ts_profile, int ts_profile_size,  float* h)
{

	if ((n_t <= 0) ||
		(ts_profile_size < 2))
	{
		return false;
	}

    float z1 = ts_profile[0].z;
    float t1 = ts_profile[0].t;
    float s1 = ts_profile[0].s;
    float rho0 = PHX_water_density_calc(t1, PHX_ATM_PRESSURE_MBAR, s1);
    float p1 = PHX_pressure_by_depth_calc(z1, PHX_ATM_PRESSURE_MBAR, rho0, g);

    float v = PHX_speed_of_sound_UNESCO_calc(t1, p1, s1);

    if (v * tof > ts_profile[ts_profile_size - 1].z)
    {
    	return false;
    }

    int p_idx = 1;
    float z2 = ts_profile[p_idx].z;
    float t2 = ts_profile[p_idx].t;
    float s2 = ts_profile[p_idx].s;
    float p2 = PHX_pressure_by_depth_calc(z2, PHX_ATM_PRESSURE_MBAR, rho0, g);

    float dt = tof / n_t;
    *h = 0.0;
    float t, p, s, tt = 0.0f;

    while (tt < tof)
    {
        tt += dt;
        *h = *h + dt * v;

        if (*h > z2)
        {
            p1 = p2;
            t1 = t2;
            s1 = s2;
            z1 = z2;
            p_idx++;

            z2 = ts_profile[p_idx].z;
            t2 = ts_profile[p_idx].t;
            s2 = ts_profile[p_idx].s;
            p2 = PHX_pressure_by_depth_calc(z2, PHX_ATM_PRESSURE_MBAR + p1, rho0, g);
        }

        t = PHX_linterp(z1, t1, z2, t2, *h);
        p = PHX_linterp(z1, p1, z2, p2, *h);
        s = PHX_linterp(z1, s1, z2, s2, *h);
        v = PHX_speed_of_sound_UNESCO_calc(t, p, s);
    }

    return true;
}

// Calculated the freezing temperature of seawater (in °C) with specified pressure and salinity.
// According to:
// Algorithms for computation of fundamental properties of seawater.
// Unesco technical papers in marine science vol. 44, 1983, pp. 30
// https://darchive.mblwhoilibrary.org/bitstream/handle/1912/2470/059832eb.pdf
// p - pressure, mBar
// s - PSU
float PHX_water_fpoint_calc(float p, float s)
{
   return (-0.0575f + 1.710523E-3f * _sqrt_f32(s) - 2.154996E-4f * s) * s - 7.53E-6f * p;
}
