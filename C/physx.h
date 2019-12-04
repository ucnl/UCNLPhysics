// (C) Aleksandr Dikarev, 2016

#ifndef _PHYSX_H_
#define _PHYSX_H_

#define PHX_FWTR_DENSITY_KGM3        (998.02f)       // Fresh water density at 20°C
#define PHX_FWTR_SOUND_SPEED_MPS     (1500.0f)       //
#define PHX_FWTR_SALINITY_PSU        (0.0f)          //
#define PHX_GRAVITY_ACC_MPS2         (9.80665f)      // ISO 80000-3:2006
#define PHX_ATM_PRESSURE_MBAR        (1013.25f)      // Average at sea level

extern float _abs_f32(float value);
extern float _sqrt_f32(float value);
extern float _sin_f32(float alpha);

// calculates in situ density.
// millero et al 1980, deep-sea res.,27a,255-264
// jpots ninth report 1978,tenth report 1980
// t: temperature, Celsius degree
// p: pressure, mBar
// s: salinity, ppm
// result: kg/m3
float PHYSX_WaterDensity_Calc(float t, float p, float s);

// The UNESCO equation: Chen and Millero (1977)
// t: temperature, Celsius degree
// p: pressure, mBar
// s: salinity, ppm
// result: m/s
float PHYSX_SpeedOfSound_Calc(float t, float p, float s);

// calculates gravity at sea level vs latitude
// WGS84 ellipsoid gravity formula
// latitude: signed float [-90..90]
// result: m/s2
float PHYSX_GravityConstant_Calc(float latitude);

// Hydrostatic pressure
// p: pressure in mBar
// p0: pressure on the water surface
// rho: water density, kg/m3
// result: depth in m
float PHYSX_DepthByPressure_Calc(float p, float p0, float rho, float g);

#endif
