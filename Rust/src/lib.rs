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

const C00: f64 = 1402.388;
const C01: f64 = 5.03830;
const C02: f64 = -5.81090E-2;
const C03: f64 = 3.3432E-4;
const C04: f64 = -1.47797E-6;
const C05: f64 = 3.1419E-9;
const C10: f64 = 0.153563;
const C11: f64 = 6.8999E-4;
const C12: f64 = -8.1829E-6;
const C13: f64 = 1.3632E-7;
const C14: f64 = -6.1260E-10;
const C20: f64 = 3.1260E-5;
const C21: f64 = -1.7111E-6;
const C22: f64 = 2.5986E-8;
const C23: f64 = -2.5353E-10;
const C24: f64 = 1.0415E-12;
const C30: f64 = -9.7729E-9;
const C31: f64 = 3.8513E-10;
const C32: f64 = -2.3654E-12;
const A00: f64 = 1.389;
const A01: f64 = -1.262E-2;
const A02: f64 = 7.166E-5;
const A03: f64 = 2.008E-6;
const A04: f64 = -3.21E-8;
const A10: f64 = 9.4742E-5;
const A11: f64 = -1.2583E-5;
const A12: f64 = -6.4928E-8;
const A13: f64 = 1.0515E-8;
const A14: f64 = -2.0142E-10;
const A20: f64 = -3.9064E-7;
const A21: f64 = 9.1061E-9;
const A22: f64 = -1.6009E-10;
const A23: f64 = 7.994E-12;
const A30: f64 = 1.100E-10;
const A31: f64 = 6.651E-12;
const A32: f64 = -3.391E-13;
const B00: f64 = -1.922E-2;
const B01: f64 = -4.42E-5;
const B10: f64 = 7.3637E-5;
const B11: f64 = 1.7950E-7;
const D00: f64 = 1.727E-3;
const D10: f64 = -7.9836E-6;


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
pub fn phx_speed_of_sound_calc(t_c: f64, p_mbar: f64, s_psu: f64) -> f64 {

    let t2 = t_c * t_c;
    let t3 = t2 * t_c;
    let t4 = t3 * t_c;
    let p = p_mbar / 1000.0;

    let c_w = (C00 + C01 * t_c + C02 * t2 + C03 * t3 + C04 * t4 + C05 * t4 * t_c) +
             (C10 + C11 * t_c + C12 * t2 + C13 * t3 + C14 * t4) * p +
             (C20 + C21 * t_c + C22 * t2 + C23 * t3 + C24 * t4) * p * p +
             (C30 + C31 * t_c + C32 * t2) * p * p * p;

    let a = (A00 + A01 * t_c + A02 * t2 + A03 * t3 + A04 * t4) +
            (A10 + A11 * t_c + A12 * t2 + A13 * t3 + A14 * t4) * p +
            (A20 + A21 * t_c + A22 * t2 + A23 * t3) * p * p +
            (A30 + A31 * t_c + A32 * t2) * p * p * p;

    let b = B00 + B01 * t_c + (B10 + B11 * t_c) * p;

    let d = D00 + D10 * p;

    (c_w + a * s_psu + b * s_psu.powi(3).sqrt() + d * s_psu * s_psu)
}

/// Calculates gravity at sea level vs latitude
/// WGS84 ellipsoid gravity formula
pub fn phx_gravity_constant_wgs84_calc(lat_rad: f64) -> f64 {

    let phi_sq = lat_rad.sin().powi(2);    
    (PHX_GE * ((1.0 + PHX_K * phi_sq) / (1.0 - PHX_E * phi_sq).sqrt()))
}

/// calculates distance from the water surface where pressure is p0 to the point, where pressure is p
pub fn phx_depth_dy_pressure_calc(p_mbar: f64, p0_mbar: f64, rho_kg_by_m3: f64, g_m_by_s2: f64) -> f64 {

    (100.0 * (p_mbar - p0_mbar) / (rho_kg_by_m3 * g_m_by_s2))
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

        assert_approx_eq!(phx_water_density_calc(0.0, PHX_ATM_PRESSURE_MBAR, 0.0), 999.9, 1E-1);
        assert_approx_eq!(phx_water_density_calc(4.0, PHX_ATM_PRESSURE_MBAR, 0.0), 1000.0, 1E-1);
        assert_approx_eq!(phx_water_density_calc(20.0, PHX_ATM_PRESSURE_MBAR, 0.0), 998.2, 1E-1);
        assert_approx_eq!(phx_water_density_calc(40.0, PHX_ATM_PRESSURE_MBAR, 0.0), 992.2, 1E-1);
        
        assert_approx_eq!(phx_water_density_calc(0.0, PHX_ATM_PRESSURE_MBAR, 35.0), 1028.13, 1E-1);
        assert_approx_eq!(phx_water_density_calc(0.0, PHX_ATM_PRESSURE_MBAR + 100000.0, 35.0), 1032.85, 1E-1);
        assert_approx_eq!(phx_water_density_calc(0.0, PHX_ATM_PRESSURE_MBAR + 200000.0, 35.0), 1037.47, 1E-1);
        assert_approx_eq!(phx_water_density_calc(0.0, PHX_ATM_PRESSURE_MBAR + 400000.0, 35.0), 1046.40, 1E-1);
        assert_approx_eq!(phx_water_density_calc(0.0, PHX_ATM_PRESSURE_MBAR + 600000.0, 35.0), 1054.95, 1E-1);
        assert_approx_eq!(phx_water_density_calc(0.0, PHX_ATM_PRESSURE_MBAR + 800000.0, 35.0), 1063.15, 1E-1);
        assert_approx_eq!(phx_water_density_calc(0.0, PHX_ATM_PRESSURE_MBAR + 1000000.0, 35.0), 1071.04, 1E-1);
    }

    #[test]
    fn phx_speed_of_sound_calc_test() {
        
        assert_approx_eq!(phx_speed_of_sound_calc(0.0, PHX_ATM_PRESSURE_MBAR, 0.0), 1403.0, 1.0);
        assert_approx_eq!(phx_speed_of_sound_calc(5.0, PHX_ATM_PRESSURE_MBAR, 0.0), 1427.0, 1.0);
        assert_approx_eq!(phx_speed_of_sound_calc(10.0, PHX_ATM_PRESSURE_MBAR, 0.0), 1447.0, 1.0);
        assert_approx_eq!(phx_speed_of_sound_calc(20.0, PHX_ATM_PRESSURE_MBAR, 0.0), 1482.0, 1.0);
        assert_approx_eq!(phx_speed_of_sound_calc(30.0, PHX_ATM_PRESSURE_MBAR, 0.0), 1509.0, 1.0);
        assert_approx_eq!(phx_speed_of_sound_calc(40.0, PHX_ATM_PRESSURE_MBAR, 0.0), 1529.0, 1.0);

        assert_approx_eq!(phx_speed_of_sound_calc(0.0, PHX_ATM_PRESSURE_MBAR, 35.0), 1449.0, 1.0);
        assert_approx_eq!(phx_speed_of_sound_calc(10.0, PHX_ATM_PRESSURE_MBAR, 35.0), 1490.0, 1.0);
        assert_approx_eq!(phx_speed_of_sound_calc(20.0, PHX_ATM_PRESSURE_MBAR, 35.0), 1522.0, 1.0);
        assert_approx_eq!(phx_speed_of_sound_calc(30.0, PHX_ATM_PRESSURE_MBAR, 35.0), 1546.0, 1.0);
        assert_approx_eq!(phx_speed_of_sound_calc(40.0, PHX_ATM_PRESSURE_MBAR, 35.0), 1563.0, 1.0);
    }    

    #[test]
    fn phx_gravity_constant_wgs84_calc_test() {

        assert_approx_eq!(phx_gravity_constant_wgs84_calc((0.0 as f64).to_radians()), 9.7801, 1E-3);
        assert_approx_eq!(phx_gravity_constant_wgs84_calc((10.0 as f64).to_radians()), 9.7819, 1E-3);
        assert_approx_eq!(phx_gravity_constant_wgs84_calc((20.0 as f64).to_radians()), 9.7863, 1E-3);
        assert_approx_eq!(phx_gravity_constant_wgs84_calc((30.0 as f64).to_radians()), 9.7932, 1E-3);
        assert_approx_eq!(phx_gravity_constant_wgs84_calc((40.0 as f64).to_radians()), 9.8016, 1E-3);
        assert_approx_eq!(phx_gravity_constant_wgs84_calc((50.0 as f64).to_radians()), 9.8100, 1E-3);
        assert_approx_eq!(phx_gravity_constant_wgs84_calc((60.0 as f64).to_radians()), 9.8200, 1E-3);
        assert_approx_eq!(phx_gravity_constant_wgs84_calc((70.0 as f64).to_radians()), 9.8261, 1E-3);
        assert_approx_eq!(phx_gravity_constant_wgs84_calc((80.0 as f64).to_radians()), 9.8300, 1E-3);
        assert_approx_eq!(phx_gravity_constant_wgs84_calc((90.0 as f64).to_radians()), 9.8322, 1E-3);

        for lt in 0..90 {
            let lt_rad = (lt as f64).to_radians();
            assert_approx_eq!(phx_gravity_constant_wgs84_calc(lt_rad), 
                              phx_gravity_constant_wgs84_calc(-lt_rad),
                              1E-6);
        }
    }
}
