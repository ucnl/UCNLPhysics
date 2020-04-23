use std::f64;

pub const PHX_FWTR_DENSITY_KGM3: f64        = 998.02;  // Fresh water density at 20°C
pub const PHX_FWTR_SOUND_SPEED_MPS: f64     = 1500.0;  // Default speed of sound in water
pub const PHX_FWTR_SOUND_SPEED_MPS_MIN: f64 = 1300.0;  // Min value for speed of sound
pub const PHX_FWTR_SOUND_SPEED_MPS_MAX: f64 = 1600.0;  // Max value for speed of sound
pub const PHX_FWTR_SALINITY_PSU: f64        = 0.0;     // Default water salinity, PSU
pub const PHX_GRAVITY_ACC_MPS2: f64         = 9.80665; // ISO 80000-3:2006
pub const PHX_ATM_PRESSURE_MBAR: f64        = 1013.25; // Average at sea level

const PHX_GE: f64 = 9.7803253359;
const PHX_K: f64  = 0.00193185265241;
const PHX_E: f64  = 0.00669437999013;

// Interpolates a value with given x coordinate by two given points (x1,y1) and (x2,y2)
pub fn phx_linterp(x1: f64, y1: f64, x2: f64, y2: f64, x: f64) -> f64 {
    y1 + (x - x1)*(y2 - y1)/(x2 - x1)
}

/// calculates in situ density of water
/// millero et al 1980, deep-sea res.,27a,255-264
/// jpots ninth report 1978,tenth report 1980
pub fn phx_water_density_calc(t_c: f64, p_mbar: f64, s_psu: f64) -> f64 {

    let p = p_mbar / 1000.0;
    let sr = s_psu.sqrt();
    let sig = (((4.8314E-4 * s_psu) +
               ((-1.6546E-6 * t_c + 1.0227E-4) * t_c - 5.72466E-3) * sr +
               (((5.3875E-9 * t_c - 8.2467E-7) * t_c + 7.6438E-5) * t_c - 4.0899E-3) * t_c + 0.824493) * s_psu) +
               ((((6.536332E-9 * t_c - 1.120083E-6) * t_c + 1.001685E-4) * t_c - 9.095290E-3) * t_c + 6.793952E-2) * t_c - 0.157406;

    let b = ((9.1697E-10 * t_c + 2.0816E-8) * t_c - 9.9348E-7) * s_psu + (5.2787E-8 * t_c - 6.12293E-6) * t_c + 8.50935E-5;

    let k0 = (((((-5.3009E-4 * t_c + 1.6483E-2) * t_c + 7.944E-2) * sr) + 
              ((-6.1670E-5 * t_c + 1.09987E-2) * t_c - 0.603459) * t_c + 54.6746) * s_psu) +
               (((-5.155288E-5 * t_c + 1.360477E-2) * t_c - 2.327105) * t_c + 148.4206) * t_c + 19652.21;

    let a = (1.91075E-4 * sr + (-1.6078E-6 * t_c - 1.0981E-5) * t_c + 2.2838E-3) * s_psu +
             ((-5.77905E-7 * t_c + 1.16092E-4) * t_c + 1.43713E-3) * t_c + 3.239908;

    let k = (b * p + a) * p + k0;

    (1000.0 + (k * sig + 1000.0 * p) / (k - p))
}

/// The UNESCO equation: Chen and Millero (1977)
pub fn phx_speed_of_sound_unesco_calc(t: f64, p: f64, s: f64) -> f64 {

    /*
    let t2 = t * t;
    let t3 = t2 * t;
    let t4 = t3 * t;
    let p = p / 1000.0;
    let p2 = p * p;
    let p3 = p2 * p;

    let c_w = (1402.388 + 5.03830 * t + -5.81090E-2 * t2 + 3.3432E-4 * t3 + -1.47797E-6 * t4 + 3.1419E-9 * t4 * t) +
             (0.153563 + 6.8999E-4 * t + -8.1829E-6 * t2 + 1.3632E-7 * t3 + -6.1260E-10 * t4) * p +
             (3.1260E-5 +  -1.7111E-6 * t +  2.5986E-8 * t2 + -2.5353E-10 * t3 + 1.0415E-12 * t4) * p2 +
             (-9.7729E-9 + 3.8513E-10 * t + -2.3654E-12 * t2) * p3;

    let a = (1.389 + -1.262E-2 * t + 7.166E-5 * t2 + 2.008E-6 * t3 + -3.21E-8 * t4) +
            (9.4742E-5 + -1.2583E-5 * t + -6.4928E-8 * t2 + 1.0515E-8 * t3 + -2.0142E-10 * t4) * p +
            (-3.9064E-7 + 9.1061E-9 * t + -1.6009E-10 * t2 + 7.994E-12 * t3) * p2 +
            (1.100E-10 + 6.651E-12 * t + -3.391E-13 * t2) * p3;

    let b = -1.922E-2 + -4.42E-5 * t + (7.3637E-5 + 1.7950E-7 * t) * p;

    let d = 1.727E-3 + -7.9836E-6 * p;
  
    (c_w + a * s + b * s.powi(3).sqrt() + d * s * s)
    */

    let p = p / 1000.0;
    let sr = s.abs().sqrt();

    let d = 1.727E-3 - 7.9836E-6 * p;

    let b_1 = 7.3637E-5 + 1.7945E-7 * t;
    let b_0 = -1.922E-2 - 4.42E-5 * t;
    let b = b_0 + b_1 * p;

    let a_3 = (-3.389E-13 * t + 6.649E-12)  * t + 1.100E-10;
    let a_2 = ((7.988E-12 * t - 1.6002E-10) * t + 9.1041E-9) * t - 3.9064E-7;
    let a_1 = (((-2.0122E-10 * t + 1.0507E-8)  * t - 6.4885E-8) * t - 1.2580E-5) * t + 9.4742E-5;
    let a_0 = (((-3.21E-8 * t + 2.006E-6) * t + 7.164E-5) * t - 1.262E-2) * t + 1.389;
    let a = ((a_3 * p + a_2) * p + a_1) * p + a_0;

    let c_3 = (-2.3643E-12 * t + 3.8504E-10) * t - 9.7729E-9;
    let c_2 = (((1.0405E-12 * t - 2.5335E-10) * t + 2.5974E-8) * t - 1.7107E-6)  * t + 3.1260E-5;
    let c_1 = (((-6.1185E-10 * t + 1.3621E-7)  * t - 8.1788E-6) * t + 6.8982E-4)  * t + 0.153563;
    let c_0 = ((((3.1464E-9  * t - 1.47800E-6) * t + 3.3420E-4) * t - 5.80852E-2) * t + 5.03711) * t + 1402.388;
    let c  = ((c_3 * p + c_2) * p + c_1) * p + c_0;

    (c + (a + b * sr + d * s) * s)
}

/// Calculates gravity at sea level vs latitude
/// WGS84 ellipsoid gravity formula
pub fn phx_gravity_constant_wgs84_calc(lat_rad: f64) -> f64 {

    let phi_sq = lat_rad.sin().powi(2);    
    (PHX_GE * ((1.0 + PHX_K * phi_sq) / (1.0 - PHX_E * phi_sq).sqrt()))
}

/// calculates distance from the water surface where pressure is p0 to the point, where pressure is p
pub fn phx_depth_by_pressure_calc(p_mbar: f64, p0_mbar: f64, rho_kg_by_m3: f64, g_m_by_s2: f64) -> f64 {

    (100.0 * (p_mbar - p0_mbar) / (rho_kg_by_m3 * g_m_by_s2))
}


// Calculates pressure of a water column with given height (distance between
// the water surface and the given point) assuming constant water density
// h - depth, m
// p0 - atmospheric pressure, mBar
// rho - water density, kg/m^3
// g - gravity acceleration, m/s^2
pub fn phx_pressure_by_depth_calc(h: f64, p0: f64, rho: f64, g: f64) -> f64 {
    h * rho * g / 100.0 + p0
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
pub fn phx_depth_by_pressure_ts_profile(pm: f64, p0: f64, g: f64, ts_profile: &[(f64, f64, f64)], n_p: i32) -> f64 {

    if n_p <= 0 {
        panic!("Specified number of time intervals Nt should be greater than zero");
    }

    if ts_profile.len() < 2 {
        panic!("tsProfile has to contain at least two points");
    }

    let mut t1 = ts_profile[0].1;
    let mut s1 = ts_profile[0].2;
    let rho0 = phx_water_density_calc(t1, p0, s1);
    let mut p1 = phx_pressure_by_depth_calc(ts_profile[0].0, p0, rho0, g);

    let pe = phx_pressure_by_depth_calc(ts_profile[ts_profile.len() - 1].0, p0, rho0, g);

    if (pm < p1) || (pm > pe) {
        panic!("Specified pressure is beyond the specified TS-profile");
    }

    let mut p_idx = 1;
    let mut t2 = ts_profile[p_idx].1;
    let mut s2 = ts_profile[p_idx].2;
    let mut p2 = phx_pressure_by_depth_calc(ts_profile[p_idx].0, p0, rho0, g);

    let dp = (pm - p0) / (n_p as f64);
    let mut h = 0.0;
    
    let mut rho;
    let mut t;
    let mut p = p0;
    let mut s;    
    
    while p < pm {

        p += dp;

        if p > p2 {

            p1 = p2;
            t1 = t2;
            s1 = s2;
            p_idx += 1;

            t2 = ts_profile[p_idx].1;
            s2 = ts_profile[p_idx].2;
            p2 = phx_pressure_by_depth_calc(ts_profile[p_idx].0, p0, rho0, g);
        }

        t = phx_linterp(p1, t1, p2, t2, p);
        s = phx_linterp(p1, s1, p2, s2, p);

        rho = phx_water_density_calc(t, p, s);
        h += 1.0 / rho;
    }

    (h * 100.0 * dp / g)
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
pub fn phx_vertical_sound_path_ts_profile(tof: f64, n_t: i32, g: f64, ts_profile: &[(f64, f64, f64)]) -> f64 {
      
    if ts_profile.len() < 2 {
      panic!("tsProfile has to contain at least two points");
    }
  
    if n_t <= 0 {
      panic!("Specified number of time intervals Nt should be greater than zero");
    }
        
    let mut z1 = ts_profile[0].0;
    let mut t1 = ts_profile[0].1;
    let mut s1 = ts_profile[0].2;
    let rho0 = phx_water_density_calc(t1, PHX_ATM_PRESSURE_MBAR, s1);
    let mut p1 = phx_pressure_by_depth_calc(z1, PHX_ATM_PRESSURE_MBAR, rho0, g);
  
    let mut v = phx_speed_of_sound_unesco_calc(t1, p1, s1);
  
    if v * tof > ts_profile[ts_profile.len() - 1].0 {
      panic!("Specified time of flight is beyond the specified TS-profile");
    }
  
    let mut p_idx = 1;
    let mut z2 = ts_profile[p_idx].0;
    let mut t2 = ts_profile[p_idx].1;
    let mut s2 = ts_profile[p_idx].2;
    let mut p2 = phx_pressure_by_depth_calc(z2, PHX_ATM_PRESSURE_MBAR, rho0, g);
  
    let dt = tof / (n_t as f64);
    let mut h = 0.0;
    let mut t;
    let mut p;
    let mut s;
    let mut tt = 0.0;
  
    while tt < tof {
  
        tt += dt;
        h = h + dt * v;
  
        if h > z2 {
  
            p1 = p2;
            t1 = t2;
            s1 = s2;
            z1 = z2;
            p_idx = p_idx + 1;
  
            z2 = ts_profile[p_idx].0;
            t2 = ts_profile[p_idx].1;
            s2 = ts_profile[p_idx].2;
            p2 = phx_pressure_by_depth_calc(z2, PHX_ATM_PRESSURE_MBAR + p1, rho0, g);
        }
  
        t = phx_linterp(z1, t1, z2, t2, h);
        p = phx_linterp(z1, p1, z2, p2, h);
        s = phx_linterp(z1, s1, z2, s2, h);
        v = phx_speed_of_sound_unesco_calc(t, p, s);
    }
  
    h
}
  
// Calculated the freezing temperature of seawater (in °C) with specified pressure and salinity.
// According to:
// Algorithms for computation of fundamental properties of seawater. 
// Unesco technical papers in marine science vol. 44, 1983, pp. 30
// https://darchive.mblwhoilibrary.org/bitstream/handle/1912/2470/059832eb.pdf
// p - pressure, mBar
// s - PSU
pub fn phx_water_fpoint_calc(p: f64, s: f64) -> f64 {
   (-0.0575 + 1.710523E-3 * s.abs().sqrt() - 2.154996E-4 * s) * s - 7.53E-6 * p
}

// calculation of absorption according to:
// Francois & Garrison, J. Acoust. Soc. Am., Vol. 72, No. 6, December 1982
// f frequency (kHz)
// T Temperature (degC)
// S Salinity (ppt)
// D Depth (m)
// pH Acidity
pub fn alpha_e_francois_garrison_calc(f: f64, t: f64, s: f64, h: f64, ph: f64) -> f64 {
    // Total absorption = Boric Acid Contrib. + Magnesium Sulphate Contrib. + Pure Water Contrib.

    // Measured ambient temp
    let t_kel: f64 = 273.15 + t;
    let fsq: f64 = f * f;            

    // Calculate speed of sound (according to Francois & Garrison, JASA 72 (6) p1886)
    let c: f64 = 1412.0 + 3.21 * t + 1.19 * s + 0.0167 * h;

    // Boric acid contribution
    let a1 = (8.86 / c) * 10f64.powf(0.78 * ph - 5.0);
    let p1 = 1.0;
    let f1 = 2.8 * (s / 35.0).sqrt() * 10f64.powf(4.0 - 1245.0 / t_kel);
    let boric = (a1 * p1 * f1 * fsq) / (fsq + f1 * f1);

    // MgSO4 contribution
    let a2 = 21.44 * (s / c) * (1.0 + 0.025 * t);
    let p2 = 1.0 - 1.37E-4 * h + 6.2E-9 * h * h;
    let f2 = (8.17 * 10f64.powf(8.0 - 1990.0 / t_kel)) / (1.0 + 0.0018 * (s - 35.0));
    let mgso4 = (a2 * p2 * f2 * fsq) / (fsq + f2 * f2);

    let a3: f64;
    if t <= 20.0 {
        a3 = 4.937E-4 - 2.59E-5 * t + 9.11E-7 * t * t - 1.5E-8 * t * t * t;
    }
    else {
        a3 = 3.964E-4 - 1.146E-5 * t + 1.45E-7 * t * t - 6.5E-10 * t * t * t;
    }

    let p3 = 1.0 - 3.83E-5 * h + 4.9E-10 * h * h;
    let h2o = a3 * p3 * fsq;

    // Total absorption
    (boric + mgso4 + h2o)
}






#[macro_export]
macro_rules! assert_approx_eq {
    ($a:expr, $b:expr, $eps:expr) => {{
        let (a, b) = (&$a, &$b);
        let eps = $eps;
        assert!(
            (*a - *b).abs() < eps,
            "assertion failed \
             (left: `{:?}`, right: `{:?}`, eps: `{:?}`, real diff: `{:?}`)",
            *a,
            *b,
            eps,
            (*a - *b).abs()
        );
    }};
}

#[cfg(test)]
mod tests {

    use super::*;   

    #[test]
    fn phx_water_density_calc_test() {
        // Test according to reference fresh water density from:
        // Tanaka, M., Girard, G., Davis, R., Peuto, A., & Bignell, N. (2001).
        // Recommended table for the density of water between 0°C and 40°C based on
        // recent experimental reports. Metrologia, 38(4), 301–309. doi:10.1088/0026-1394/38/4/3 
        // 
        let ref_salinity: f64 = 0.0;
        let ref_temperatures: [f64; 41] = [  0.0,  1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0, 
                                            10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0,
                                            20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0,
                                            30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0,
                                            40.0 ];

        let ref_densities: [f64; 41] = [ 999.8428, 999.9017, 999.9429, 999.9672, 999.9749, 999.9668, 999.9431, 999.9045,
                                         999.8513, 999.7839, 999.7027, 999.6081, 999.5005, 999.3801, 999.2474, 999.1026,
                                         998.9459, 998.7778, 998.5984, 998.4079, 998.2067, 997.9950, 997.7730, 997.5408,
                                         997.2988, 997.0470, 996.7857, 996.5151, 996.2353, 995.9465, 995.6488, 995.3424,
                                         995.0275, 994.7041, 994.3724, 994.0326, 993.6847, 993.3290, 992.9654, 992.5941,
                                         992.2152 ];

        for p_idx in 0..ref_densities.len() {
            assert_approx_eq!(phx_water_density_calc(ref_temperatures[p_idx], PHX_ATM_PRESSURE_MBAR, ref_salinity), 
                              ref_densities[p_idx], 5.2E-2);
        }        

        // Testing according to fresh & seawater density from:
        // ITTC – Recommended Procedures 7.5-02-01-03. Fresh Water and Seawater Properties 
        // Effective Date 2011 Revision 02, Updated / Edited by Approved Quality Systems Group of the 28 th ITTC 26 th ITTC 2011 
        // Date 09/2016
        // https://ittc.info/media/7503/75-02-01-03.pdf

        // From pp. 4-5, Table 1:
        let ref_salinity: f64 = 0.0;
        let ref_temperatures: [f64; 31] = [ 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0,
                                            20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0,
                                            30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0 ];
                                            
        let ref_densities: [f64; 31] = [ 999.7025, 999.6079, 999.5004, 999.3801, 999.2474, 999.1026, 998.9461,
                                         998.7780, 998.5986, 998.4083, 998.2072, 997.9955, 997.7735, 997.5414,
                                         997.2994, 997.0476, 996.7864, 996.5158, 996.2360, 995.9471, 995.6495,
                                         995.3431, 995.0281, 994.7048, 994.3731, 994.0333, 993.6855, 993.3298,
                                         992.9663, 992.5951, 992.2164 ];
        
        for p_idx in 0..ref_densities.len() {
            assert_approx_eq!(phx_water_density_calc(ref_temperatures[p_idx], PHX_ATM_PRESSURE_MBAR, ref_salinity), 
                            ref_densities[p_idx], 5.2E-2);
        }     


        // From pp. 7-8, Table 3
        let ref_salinity: f64 = 35.16; // +/- 0.007
        let ref_temperatures: [f64; 30] = [  1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0, 10.0, 
                                            11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0,
                                            21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0 ];
                                            
        let ref_densities: [f64; 30] = [ 1028.0941, 1028.0197, 1027.9327, 1027.8336, 1027.7225, 1027.6000,
                                         1027.4662, 1027.3214, 1027.1659, 1027.0000, 1026.8238, 1026.6376,
                                         1026.4416, 1026.2360, 1026.0210, 1025.7967, 1025.5633, 1025.3210,
                                         1025.0700, 1024.8103, 1024.5421, 1024.2656, 1023.9808, 1023.6881,
                                         1023.3873, 1023.0788, 1022.7626, 1022.4389, 1022.1078, 1021.7694 ];
        
        for p_idx in 0..ref_densities.len() {
            assert_approx_eq!(phx_water_density_calc(ref_temperatures[p_idx], PHX_ATM_PRESSURE_MBAR, ref_salinity), 
                            ref_densities[p_idx], 0.5);
        }     


        // From pp. 10-11, Table 4
        let ref_temperature: f64 = 15.0; 
        let ref_salinities: [f64; 31] = [ 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0,
                                          20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0,
                                          30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0 ];
                                            
        let ref_densities: [f64; 31] = [ 1006.7950, 1007.5571, 1008.3191, 1009.0812, 1009.8434, 1010.6056,
                                         1011.3680, 1012.1305, 1012.8932, 1013.6561, 1014.4192, 1015.1824,
                                         1015.9459, 1016.7097, 1017.4736, 1018.2379, 1019.0023, 1019.7670,
                                         1020.5320, 1021.2973, 1022.0628, 1022.8286, 1023.5946, 1024.3609,
                                         1025.1275, 1025.8944, 1026.6615, 1027.4289, 1028.1966, 1028.9646,
                                         1029.7328 ];
        
        
        for p_idx in 0..ref_densities.len() {
            assert_approx_eq!(phx_water_density_calc(ref_temperature, PHX_ATM_PRESSURE_MBAR, ref_salinities[p_idx]), 
                            ref_densities[p_idx], 0.2);
        }     


        // Don't remember where these values from ), but it should be ok
        let ref_temperature: f64 = 0.0;
        let ref_salinity: f64 = 35.0;
        let ref_pressures: [f64; 7] = [    0.0,  100000.0,  200000.0,  400000.0,  600000.0,  800000.0,  1000000.0 ];
        let ref_densities: [f64; 7] = [ 1028.13,   1032.85,   1037.47,   1046.40,   1054.95,   1063.15,    1071.04 ];

        for p_idx in 0..ref_densities.len() {
            assert_approx_eq!(phx_water_density_calc(ref_temperature, PHX_ATM_PRESSURE_MBAR + ref_pressures[p_idx], ref_salinity), 
                            ref_densities[p_idx], 0.1);
        }   
    }

    #[test]
    fn phx_speed_of_sound_calc_test() {
        // Reference values according to 
        // Algorithms for computation of fundamental properties of seawater. 
        // Unesco technical papers in marine science vol. 44, 1983, pp. 28
        // https://darchive.mblwhoilibrary.org/bitstream/handle/1912/2470/059832eb.pdf
        
        let ref_t: [f64; 5] = [ 0.0, 10.0, 20.0, 30.0, 40.0 ];
        let ref_p: [f64; 11] = [ 0.0, 1E5, 2E5, 3E5, 4E5, 5E5, 6E5, 7E5, 8E5, 9E5, 1E6 ];

        let ref_v_s25: [[f64; 5]; 11] = [ [ 1435.8, 1477.7, 1510.3, 1535.2, 1553.4 ],
                                          [ 1452.0, 1494.1, 1527.0, 1552.1, 1570.6 ],
                                          [ 1468.6, 1510.7, 1543.6, 1569.0, 1587.6 ],
                                          [ 1485.6, 1527.5, 1560.3, 1585.7, 1604.5 ],
                                          [ 1502.8, 1544.3, 1576.9, 1602.4, 1621.3 ],
                                          [ 1520.4, 1561.3, 1593.6, 1619.0, 1638.0 ],
                                          [ 1538.1, 1578.4, 1610.3, 1635.5, 1654.6 ],
                                          [ 1556.0, 1595.6, 1626.9, 1651.9, 1671.0 ],
                                          [ 1574.1, 1612.8, 1643.5, 1668.2, 1687.2 ],
                                          [ 1592.2, 1630.1, 1660.2, 1684.5, 1703.3 ],
                                          [ 1610.4, 1647.4, 1676.8, 1700.6, 1719.2 ]
                                        ];

        let ref_v_s30: [[f64; 5]; 11] = [ [ 1442.5, 1483.7, 1515.9, 1540.4, 1558.3 ],
                                          [ 1458.8, 1500.2, 1532.5, 1557.3, 1575.4 ],
                                          [ 1475.4, 1516.8, 1549.2, 1574.1, 1592.3 ],
                                          [ 1492.4, 1533.6, 1565.8, 1590.8, 1609.2 ],
                                          [ 1509.7, 1550.4, 1582.4, 1607.4, 1626.0 ],
                                          [ 1527.2, 1567.4, 1599.1, 1624.0, 1642.7 ],
                                          [ 1544.9, 1584.4, 1615.7, 1640.4, 1659.2 ],
                                          [ 1562.7, 1601.5, 1632.3, 1656.8, 1675.6 ],
                                          [ 1580.7, 1618.7, 1648.9, 1673.1, 1691.8 ],
                                          [ 1598.8, 1636.0, 1665.5, 1689.3, 1707.8 ],
                                          [ 1616.8, 1653.3, 1682.1, 1705.4, 1723.5 ]
                                        ];

        let ref_v_s35: [[f64; 5]; 11] = [ [ 1449.1, 1489.8, 1521.5, 1545.6, 1563.2 ],
                                          [ 1465.5, 1506.3, 1538.1, 1562.4, 1580.2 ],
                                          [ 1482.3, 1523.0, 1554.7, 1579.2, 1597.1 ],
                                          [ 1499.3, 1539.7, 1571.3, 1595.9, 1613.9 ],
                                          [ 1516.5, 1556.5, 1587.9, 1612.5, 1630.7 ],
                                          [ 1534.0, 1573.4, 1604.5, 1629.0, 1647.3 ],
                                          [ 1551.6, 1590.4, 1621.0, 1645.4, 1663.8 ],
                                          [ 1569.4, 1607.5, 1637.6, 1661.7, 1680.1 ],
                                          [ 1587.2, 1624.6, 1654.1, 1677.9, 1696.2 ],
                                          [ 1605.2, 1641.8, 1670.6, 1694.0, 1712.2 ],
                                          [ 1623.2, 1659.0, 1687.2, 1710.1, 1727.8 ]
                                        ];                                        


        for t_idx in 0..ref_t.len() {
            for p_idx in 0..ref_p.len() {                
                assert_approx_eq!(phx_speed_of_sound_unesco_calc(ref_t[t_idx], ref_p[p_idx], 25.0), ref_v_s25[p_idx][t_idx], 0.1);                
                assert_approx_eq!(phx_speed_of_sound_unesco_calc(ref_t[t_idx], ref_p[p_idx], 30.0), ref_v_s30[p_idx][t_idx], 0.1);                
                assert_approx_eq!(phx_speed_of_sound_unesco_calc(ref_t[t_idx], ref_p[p_idx], 35.0), ref_v_s35[p_idx][t_idx], 0.1);
            }
        }

        assert_approx_eq!(phx_speed_of_sound_unesco_calc(40.0, 1000000.0, 40.0), 1731.995, 0.001);
    }    

    #[test]
    fn phx_gravity_constant_wgs84_calc_test() {

        let ref_lat_deg: [f64; 10] = [ 0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0 ];
        let ref_g: [f64; 10] = [ 9.7801, 9.7819, 9.7863, 9.7932, 9.8016, 9.8100, 9.8200, 9.8261, 9.8300, 9.8322 ];

        for p_idx in 0..ref_lat_deg.len() {
            assert_approx_eq!(phx_gravity_constant_wgs84_calc(ref_lat_deg[p_idx].to_radians()), ref_g[p_idx], 1E-3);
        }
        
        for lt in 0..90 {
            let lt_rad = (lt as f64).to_radians();
            assert_approx_eq!(phx_gravity_constant_wgs84_calc(lt_rad), 
                              phx_gravity_constant_wgs84_calc(-lt_rad),
                              1E-6);
        }
    }

    #[test]
    fn phx_pressure_by_depth_calc_test() {

        // Reference values according to 
        // Algorithms for computation of fundamental properties of seawater. 
        // Unesco technical papers in marine science vol. 44, 1983, pp. 28
        // https://darchive.mblwhoilibrary.org/bitstream/handle/1912/2470/059832eb.pdf

        let ref_salinity = 35.0;
        let ref_temperature = 0.0;
        let ref_pressure: [f64; 11] = [ 5E4, 1E5, 2E5, 3E5, 4E5, 5E5, 6E5, 7E5, 8E5, 9E5, 1E6 ];
        let ref_lats_deg: [f64; 5] = [ 0.0, 30.0, 45.0, 60.0, 90.0 ];
        let ref_dpt_lat: [[f64; 11]; 5] = [ [ 496.65, 992.12, 1979.55, 2962.43, 3940.88, 4915.04, 5885.03, 6850.95, 7812.93, 8771.07, 9725.47 ],
                                            [ 496.00, 990.81, 1976.94, 2958.52, 3935.68, 4908.56, 5877.27, 6841.92, 7802.63, 8759.51, 9712.65 ],
                                            [ 495.34, 989.50, 1974.33, 2954.61, 3930.49, 4902.08, 5869.51, 6832.89, 7792.33, 8747.95, 9699.84 ],
                                            [ 494.69, 988.19, 1971.72, 2950.71, 3925.30, 4895.60, 5861.76, 6823.86, 7782.04, 8736.40, 9687.84 ],
                                            [ 494.03, 986.88, 1969.11, 2946.81, 3920.10, 4889.13, 5854.01, 6814.84, 7771.76, 8724.85, 9674.23 ]
                                          ];

        for l_idx in 0..ref_lats_deg.len() {

            for p_idx in 0..ref_pressure.len() {

                let g = phx_gravity_constant_wgs84_calc(ref_lats_deg[l_idx].to_radians());
                // taking the water density in the midpoint to consider the compression of water
                let rho = phx_water_density_calc(ref_temperature, ref_pressure[p_idx] / 2.0, ref_salinity);
                let h_est = phx_depth_by_pressure_calc(ref_pressure[p_idx], 0.0, rho, g);

                // calculated values should deviate less than 0.07% of reference values
                assert_approx_eq!(h_est, ref_dpt_lat[l_idx][p_idx], ref_dpt_lat[l_idx][p_idx] * 0.0007);
            }
        }
    }

    #[test]
    fn phx_water_fpoint_calc_test() {
        
        // Reference values according to 
        // Algorithms for computation of fundamental properties of seawater. 
        // Unesco technical papers in marine science vol. 44, 1983, pp. 30
        // https://darchive.mblwhoilibrary.org/bitstream/handle/1912/2470/059832eb.pdf

        assert_approx_eq!(phx_water_fpoint_calc(50000.0, 40.0), -2.588567, 1E-6);

        let ref_pressure: [f64; 6] = [ 0.0, 1E4, 2E4, 3E4, 4E4, 5E4 ];
        let ref_salinity: [f64; 8] = [ 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0 ];
        let ref_tf: [[f64; 6]; 8] = [ [ -0.274, -0.349, -0.424, -0.500, -0.575, -0.650 ],
                                      [ -0.542, -0.618, -0.693, -0.768, -0.844, -0.919 ],
                                      [ -0.812, -0.887, -0.962, -1.038, -1.113, -1.188 ],
                                      [ -1.083, -1.159, -1.234, -1.309, -1.384, -1.460 ],
                                      [ -1.358, -1.434, -1.509, -1.584, -1.660, -1.735 ],
                                      [ -1.638, -1.713, -1.788, -1.864, -1.939, -2.014 ],
                                      [ -1.922, -1.998, -2.073, -2.148, -2.224, -2.299 ],
                                      [ -2.212, -2.287, -2.363, -2.438, -2.513, -2.589 ]
                                    ];

        for s_idx in 0..ref_salinity.len() {

            for p_idx in 0..ref_pressure.len() {

                let tf_est = phx_water_fpoint_calc(ref_pressure[p_idx], ref_salinity[s_idx]);                
                assert_approx_eq!(tf_est, ref_tf[s_idx][p_idx], 0.001);
            }
        }
    }
}
