// (C) Aleksandr Dikarev, 2016

#include <physx.h>

#define PHYSX_GE (9.7803253359f)
#define PHYSX_K  (0.00193185265241f)
#define PHYSX_E  (0.00669437999013f)

#define C00	1402.388f
#define C01	5.03830f
#define C02	-5.81090E-2f
#define C03	3.3432E-4f
#define C04	-1.47797E-6f
#define C05	3.1419E-9f
#define C10	0.153563f
#define C11	6.8999E-4f
#define C12	-8.1829E-6f
#define C13	1.3632E-7f
#define C14	-6.1260E-10f
#define C20	3.1260E-5f
#define C21	-1.7111E-6f
#define C22	2.5986E-8f
#define C23	-2.5353E-10f
#define C24	1.0415E-12f
#define C30	-9.7729E-9f
#define C31	3.8513E-10f
#define C32	-2.3654E-12f
#define A00	1.389f
#define A01	-1.262E-2f
#define A02	7.166E-5f
#define A03	2.008E-6f
#define A04	-3.21E-8f
#define A10	9.4742E-5f
#define A11	-1.2583E-5f
#define A12	-6.4928E-8f
#define A13	1.0515E-8f
#define A14	-2.0142E-10f
#define A20	-3.9064E-7f
#define A21	9.1061E-9f
#define A22	-1.6009E-10f
#define A23	7.994E-12f
#define A30	1.100E-10f
#define A31	6.651E-12f
#define A32	-3.391E-13f
#define B00	-1.922E-2f
#define B01	-4.42E-5f
#define B10	7.3637E-5f
#define B11	1.7950E-7f
#define D00	1.727E-3f
#define D10	-7.9836E-6f


// calculates in situ density.
// millero et al 1980, deep-sea res.,27a,255-264
// jpots ninth report 1978,tenth report 1980
// t: temperature, Celsius degree
// p: pressure, mBar
// s: salinity, ppm
// result: kg/m3
float PHYSX_WaterDensity_Calc(float t, float p, float s)
{
	float sig, k0, b, a, k, sr;
	p = p / 1000.0f;
	sr = _sqrt_f32(_abs_f32(s));
	sig = 4.8314E-4f * s;
	sig += ((-1.6546E-6f * t + 1.0227E-4f) * t - 5.72466E-3f) * sr;
	sig += (((5.3875E-9f * t - 8.2467E-7f) * t + 7.6438E-5f) * t - 4.0899E-3f) * t + 0.824493f;
	sig *= s;
	sig += ((((6.536332E-9f * t - 1.120083E-6f) * t + 1.001685E-4f) * t - 9.095290E-3f) * t + 6.793952E-2f) * t - 0.157406f;
	b = ((9.1697E-10f * t + 2.0816E-8f) * t - 9.9348E-7f) * s + (5.2787E-8f * t - 6.12293E-6f) * t + 8.50935E-5f;
	k0 = ((-5.3009E-4f * t + 1.6483E-2f) * t + 7.944E-2f) * sr;
	k0 += ((-6.1670E-5f * t + 1.09987E-2f) * t - 0.603459f) * t + 54.6746f;
	k0 *= s;
	k0 += (((-5.155288E-5f * t + 1.360477E-2f) * t - 2.327105f) * t + 148.4206f) * t + 19652.21f;
	a = 1.91075E-4f * sr + (-1.6078E-6f * t - 1.0981E-5f) * t + 2.2838E-3f;
	a *= s;
	a += ((-5.77905E-7f * t + 1.16092E-4f) * t + 1.43713E-3f) * t + 3.239908f;
	k = (b * p + a) * p + k0;
	sig = 1000.0f + (k * sig + 1000.0f * p) / (k - p);

	return sig;
}

// The UNESCO equation: Chen and Millero (1977)
// t: temperature, Celsius degree
// p: pressure, mBar
// s: salinity, ppm
// result: m/s
float PHYSX_SpeedOfSound_Calc(float t, float p_mbar, float s)
{
	float t2 = t * t;
	float t3 = t2 * t;
	float t4 = t3 * t;
	float p = p_mbar / 1000.0f;

	float Cw = (C00 + C01 * t + C02 * t2 + C03 * t3 + C04 * t4 + C05 * t4 * t) +
			   (C10 + C11 * t + C12 * t2 + C13 * t3 + C14 * t4) * p +
			   (C20 + C21 * t + C22 * t2 + C23 * t3 + C24 * t4) * p * p +
			   (C30 + C31 * t + C32 * t2) * p * p * p;

	float A = (A00 + A01 * t + A02 * t2 + A03 * t3 + A04 * t4) +
			  (A10 + A11 * t + A12 * t2 + A13 * t3 + A14 * t4) * p +
			  (A20 + A21 * t + A22 * t2 + A23 * t3) * p * p +
			  (A30 + A31 * t + A32 * t2) * p * p * p;

	float B = B00 + B01 * t + (B10 + B11 * t) * p;

	float D = D00 + D10 * p;

	float c = Cw + A * s + B * _sqrt_f32(s * s * s) + D * s * s;

	return c;
}

// calculates gravity at sea level vs latitude
// WGS84 ellipsoid gravity formula
// latitude: signed float [-90..90]
// result: m/s2
float PHYSX_GravityConstant_Calc(float latitude)
{
	float phi_sq = _sin_f32(latitude);
	phi_sq *= phi_sq;
	return PHYSX_GE * ((1 + PHYSX_K * phi_sq) / _sqrt_f32(1 - PHYSX_E * phi_sq));
}

// Hydrostatic pressure
// p: pressure in mBar
// p0: pressure on the water surface
// rho: water density, kg/m3
// result: depth in m
float PHYSX_DepthByPressure_Calc(float p, float p0, float rho, float g)
{
	return (100.0f * (p - p0) / (rho * g));
}
