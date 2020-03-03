// (C) Aleksandr Dikarev, 2016-2020
// physx.h

#ifndef _CORE_PHX_H_
#define _CORE_PHX_H_

#if !defined(bool)
typedef unsigned char bool;

#define true  ((bool)(1u))
#define false ((bool)(0u))
#endif

#define PHX_FWTR_DENSITY_KGM3        (998.02f)       // Fresh water density at 20°C
#define PHX_FWTR_SOUND_SPEED_MPS     (1500.0f)       //
#define PHX_FWTR_SALINITY_PSU        (0.0f)          //
#define PHX_GRAVITY_ACC_MPS2         (9.80665f)      // ISO 80000-3:2006
#define PHX_ATM_PRESSURE_MBAR        (1013.25f)      // Average at sea level

extern float _abs_f32(float value);
extern float _sqrt_f32(float value);
extern float _sin_f32(float alpha);

typedef struct
{
	float z;
	float t;
	float s;
} TSProfileItem_Struct;

// Interpolates a value with given x coordinate by two given points (x1,y1) and (x2,y2)
float PHX_linterp(float x1, float y1, float x2, float y2, float x);

/// calculates in situ density of water
/// millero et al 1980, deep-sea res.,27a,255-264
/// jpots ninth report 1978,tenth report 1980
float PHX_water_density_calc(float t, float p, float s);

/// The UNESCO equation: Chen and Millero (1977)
float PHX_speed_of_sound_UNESCO_calc(float t, float p, float s);

/// Calculates gravity at sea level vs latitude
/// WGS84 ellipsoid gravity formula
float PHX_gravity_constant_wgs84_calc(float phi);

/// calculates distance from the water surface where pressure is p0 to the point, where pressure is p
float PHX_depth_by_pressure_calc(float p, float p0, float rho, float g);

// Calculates pressure of a water column with given height (distance between
// the water surface and the given point) assuming constant water density
// h - depth, m
// p0 - atmospheric pressure, mBar
// rho - water density, kg/m^3
// g - gravity acceleration, m/s^2
float PHX_pressure_by_depth_calc(float h, float p0, float rho, float g);

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
bool PHX_depth_by_pressure_ts_profile(float pm, float p0, float g, int n_p, TSProfileItem_Struct* ts_profile, int ts_profile_size,  float* dpt);

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
bool PHX_vertical_sound_path_ts_profile(float tof, float g, int n_t, TSProfileItem_Struct* ts_profile, int ts_profile_size,  float* h);

// Calculated the freezing temperature of seawater (in °C) with specified pressure and salinity.
// According to:
// Algorithms for computation of fundamental properties of seawater.
// Unesco technical papers in marine science vol. 44, 1983, pp. 30
// https://darchive.mblwhoilibrary.org/bitstream/handle/1912/2470/059832eb.pdf
// p - pressure, mBar
// s - PSU
float PHX_water_fpoint_calc(float p, float s);


#endif
