using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using UCNLPhysics;

namespace UCNLPhysics_tests
{
    [TestClass]
    public class PHX_Test
    {
        [TestMethod]
        public void Water_density_calc_test()
        {
            // Test according to reference fresh water density from:
            // Tanaka, M., Girard, G., Davis, R., Peuto, A., & Bignell, N. (2001).
            // Recommended table for the density of water between 0°C and 40°C based on
            // recent experimental reports. Metrologia, 38(4), 301–309. doi:10.1088/0026-1394/38/4/3 
            // 
            double ref_salinity = 0.0;
            double[] ref_temperatures = new double[] {  0.0,  1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0, 
                                                       10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0,
                                                       20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0,
                                                       30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0,
                                                       40.0 };

            double[] ref_densities = new double[] { 999.8428, 999.9017, 999.9429, 999.9672, 999.9749, 999.9668, 999.9431, 999.9045,
                                                    999.8513, 999.7839, 999.7027, 999.6081, 999.5005, 999.3801, 999.2474, 999.1026,
                                                    998.9459, 998.7778, 998.5984, 998.4079, 998.2067, 997.9950, 997.7730, 997.5408,
                                                    997.2988, 997.0470, 996.7857, 996.5151, 996.2353, 995.9465, 995.6488, 995.3424,
                                                    995.0275, 994.7041, 994.3724, 994.0326, 993.6847, 993.3290, 992.9654, 992.5941,
                                                    992.2152 };

            for (int p_idx = 0; p_idx < ref_densities.Length; p_idx++)
            {
                Assert.AreEqual(PHX.Water_density_calc(ref_temperatures[p_idx], PHX.PHX_ATM_PRESSURE_MBAR, ref_salinity),
                                ref_densities[p_idx],
                                5.2E-2);
            }

            // Testing according to fresh & seawater density from:
            // ITTC – Recommended Procedures 7.5-02-01-03. Fresh Water and Seawater Properties 
            // Effective Date 2011 Revision 02, Updated / Edited by Approved Quality Systems Group of the 28 th ITTC 26 th ITTC 2011 
            // Date 09/2016
            // https://ittc.info/media/7503/75-02-01-03.pdf

            // From pp. 4-5, Table 1:
            ref_salinity = 0.0;
            ref_temperatures = new double[] { 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0,
                                              20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0,
                                              30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0 };

            ref_densities = new double[] { 999.7025, 999.6079, 999.5004, 999.3801, 999.2474, 999.1026, 998.9461,
                                           998.7780, 998.5986, 998.4083, 998.2072, 997.9955, 997.7735, 997.5414,
                                           997.2994, 997.0476, 996.7864, 996.5158, 996.2360, 995.9471, 995.6495,
                                           995.3431, 995.0281, 994.7048, 994.3731, 994.0333, 993.6855, 993.3298,
                                           992.9663, 992.5951, 992.2164 };

            for (int p_idx = 0; p_idx < ref_densities.Length; p_idx++)
            {
                Assert.AreEqual(PHX.Water_density_calc(ref_temperatures[p_idx], PHX.PHX_ATM_PRESSURE_MBAR, ref_salinity),
                                ref_densities[p_idx],
                                5.2E-2);
            }


            // From pp. 7-8, Table 3
            ref_salinity = 35.16; // +/- 0.007
            ref_temperatures = new double[] {  1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0, 10.0, 
                                              11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0,
                                              21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0 };

            ref_densities = new double[] { 1028.0941, 1028.0197, 1027.9327, 1027.8336, 1027.7225, 1027.6000,
                                           1027.4662, 1027.3214, 1027.1659, 1027.0000, 1026.8238, 1026.6376,
                                           1026.4416, 1026.2360, 1026.0210, 1025.7967, 1025.5633, 1025.3210,
                                           1025.0700, 1024.8103, 1024.5421, 1024.2656, 1023.9808, 1023.6881,
                                           1023.3873, 1023.0788, 1022.7626, 1022.4389, 1022.1078, 1021.7694 };

            for (int p_idx = 0; p_idx < ref_densities.Length; p_idx++)
            {
                Assert.AreEqual(PHX.Water_density_calc(ref_temperatures[p_idx], PHX.PHX_ATM_PRESSURE_MBAR, ref_salinity),
                                ref_densities[p_idx],
                                0.5);
            }


            // From pp. 10-11, Table 4
            double ref_temperature = 15.0;
            double[] ref_salinities = new double[] { 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0,
                                                     20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0,
                                                     30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0 };

            ref_densities = new double[] { 1006.7950, 1007.5571, 1008.3191, 1009.0812, 1009.8434, 1010.6056,
                                           1011.3680, 1012.1305, 1012.8932, 1013.6561, 1014.4192, 1015.1824,
                                           1015.9459, 1016.7097, 1017.4736, 1018.2379, 1019.0023, 1019.7670,
                                           1020.5320, 1021.2973, 1022.0628, 1022.8286, 1023.5946, 1024.3609,
                                           1025.1275, 1025.8944, 1026.6615, 1027.4289, 1028.1966, 1028.9646,
                                           1029.7328 };

            for (int p_idx = 0; p_idx < ref_densities.Length; p_idx++)
            {
                Assert.AreEqual(PHX.Water_density_calc(ref_temperature, PHX.PHX_ATM_PRESSURE_MBAR, ref_salinities[p_idx]),
                                ref_densities[p_idx],
                                0.2);
            }

            // Don't remember where are these values from ), but it should be ok
            ref_temperature = 0.0;
            ref_salinity = 35.0;
            double[] ref_pressures = new double[] { 0.0, 100000.0, 200000.0, 400000.0, 600000.0, 800000.0, 1000000.0 };
            ref_densities = new double[] { 1028.13, 1032.85, 1037.47, 1046.40, 1054.95, 1063.15, 1071.04 };

            for (int p_idx = 0; p_idx < ref_densities.Length; p_idx++)
            {
                Assert.AreEqual(PHX.Water_density_calc(ref_temperature, PHX.PHX_ATM_PRESSURE_MBAR + ref_pressures[p_idx], ref_salinity),
                                ref_densities[p_idx],
                                0.1);
            }
        }

        [TestMethod]
        public void Speed_of_sound_calc_test()
        {
            // Reference values according to 
            // Algorithms for computation of fundamental properties of seawater. 
            // Unesco technical papers in marine science vol. 44, 1983, pp. 28
            // https://darchive.mblwhoilibrary.org/bitstream/handle/1912/2470/059832eb.pdf

            double[] ref_t = new double[] { 0.0, 10.0, 20.0, 30.0, 40.0 };
            double[] ref_p = new double[] { 0.0, 1E5, 2E5, 3E5, 4E5, 5E5, 6E5, 7E5, 8E5, 9E5, 1E6 };

            double[,] ref_v_s25 = new double[,] { { 1435.8, 1477.7, 1510.3, 1535.2, 1553.4 },
                                              { 1452.0, 1494.1, 1527.0, 1552.1, 1570.6 },
                                              { 1468.6, 1510.7, 1543.6, 1569.0, 1587.6 },
                                              { 1485.6, 1527.5, 1560.3, 1585.7, 1604.5 },
                                              { 1502.8, 1544.3, 1576.9, 1602.4, 1621.3 },
                                              { 1520.4, 1561.3, 1593.6, 1619.0, 1638.0 },
                                              { 1538.1, 1578.4, 1610.3, 1635.5, 1654.6 },
                                              { 1556.0, 1595.6, 1626.9, 1651.9, 1671.0 },
                                              { 1574.1, 1612.8, 1643.5, 1668.2, 1687.2 },
                                              { 1592.2, 1630.1, 1660.2, 1684.5, 1703.3 },
                                              { 1610.4, 1647.4, 1676.8, 1700.6, 1719.2 },
                                            };

            double[,] ref_v_s30 = new double[,] { { 1442.5, 1483.7, 1515.9, 1540.4, 1558.3 },
                                              { 1458.8, 1500.2, 1532.5, 1557.3, 1575.4 },
                                              { 1475.4, 1516.8, 1549.2, 1574.1, 1592.3 },
                                              { 1492.4, 1533.6, 1565.8, 1590.8, 1609.2 },
                                              { 1509.7, 1550.4, 1582.4, 1607.4, 1626.0 },
                                              { 1527.2, 1567.4, 1599.1, 1624.0, 1642.7 },
                                              { 1544.9, 1584.4, 1615.7, 1640.4, 1659.2 },
                                              { 1562.7, 1601.5, 1632.3, 1656.8, 1675.6 },
                                              { 1580.7, 1618.7, 1648.9, 1673.1, 1691.8 },
                                              { 1598.8, 1636.0, 1665.5, 1689.3, 1707.8 },
                                              { 1616.8, 1653.3, 1682.1, 1705.4, 1723.5 }
                                            };

            double[,] ref_v_s35 = new double[,] { { 1449.1, 1489.8, 1521.5, 1545.6, 1563.2 },
                                              { 1465.5, 1506.3, 1538.1, 1562.4, 1580.2 },
                                              { 1482.3, 1523.0, 1554.7, 1579.2, 1597.1 },
                                              { 1499.3, 1539.7, 1571.3, 1595.9, 1613.9 },
                                              { 1516.5, 1556.5, 1587.9, 1612.5, 1630.7 },
                                              { 1534.0, 1573.4, 1604.5, 1629.0, 1647.3 },
                                              { 1551.6, 1590.4, 1621.0, 1645.4, 1663.8 },
                                              { 1569.4, 1607.5, 1637.6, 1661.7, 1680.1 },
                                              { 1587.2, 1624.6, 1654.1, 1677.9, 1696.2 },
                                              { 1605.2, 1641.8, 1670.6, 1694.0, 1712.2 },
                                              { 1623.2, 1659.0, 1687.2, 1710.1, 1727.8 }
                                            };

            for (int t_idx = 0; t_idx < ref_t.Length; t_idx++)
            {
                for (int p_idx = 0; p_idx < ref_p.Length; p_idx++)
                {
                    Assert.AreEqual(PHX.Speed_of_sound_UNESCO_calc(ref_t[t_idx], ref_p[p_idx], 25.0), ref_v_s25[p_idx, t_idx], 0.1);
                    Assert.AreEqual(PHX.Speed_of_sound_UNESCO_calc(ref_t[t_idx], ref_p[p_idx], 30.0), ref_v_s30[p_idx, t_idx], 0.1);
                    Assert.AreEqual(PHX.Speed_of_sound_UNESCO_calc(ref_t[t_idx], ref_p[p_idx], 35.0), ref_v_s35[p_idx, t_idx], 0.1);
                }
            }

            Assert.AreEqual(PHX.Speed_of_sound_UNESCO_calc(40.0, 1000000.0, 40.0), 1731.995, 0.001);
        }

        [TestMethod]
        public void Gravity_constant_wgs84_calc_test()
        {
            double[] ref_lat_deg = new double[] { 0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0 };
            double[] ref_g = new double[] { 9.7801, 9.7819, 9.7863, 9.7932, 9.8016, 9.8100, 9.8200, 9.8261, 9.8300, 9.8322 };

            for (int p_idx = 0; p_idx < ref_lat_deg.Length; p_idx++)
            {
                Assert.AreEqual(PHX.Gravity_constant_wgs84_calc(ref_lat_deg[p_idx]), ref_g[p_idx], 1E-3);
            }

            for (int lt = 0; lt < 90; lt++)
            {
                Assert.AreEqual(PHX.Gravity_constant_wgs84_calc(lt),
                                PHX.Gravity_constant_wgs84_calc(-lt),
                                1E-6);
            }
        }

        [TestMethod]
        public void Pressure_by_depth_calc_test()
        {
            // Reference values according to 
            // Algorithms for computation of fundamental properties of seawater. 
            // Unesco technical papers in marine science vol. 44, 1983, pp. 28
            // https://darchive.mblwhoilibrary.org/bitstream/handle/1912/2470/059832eb.pdf

            double ref_salinity = 35.0;
            double ref_temperature = 0.0;
            double[] ref_pressures = new double[] { 5E4, 1E5, 2E5, 3E5, 4E5, 5E5, 6E5, 7E5, 8E5, 9E5, 1E6 };
            double[] ref_lats_deg = new double[] { 0.0, 30.0, 45.0, 60.0, 90.0 };
            double[,] ref_dpt_lat = { { 496.65, 992.12, 1979.55, 2962.43, 3940.88, 4915.04, 5885.03, 6850.95, 7812.93, 8771.07, 9725.47 },
                                      { 496.00, 990.81, 1976.94, 2958.52, 3935.68, 4908.56, 5877.27, 6841.92, 7802.63, 8759.51, 9712.65 },
                                      { 495.34, 989.50, 1974.33, 2954.61, 3930.49, 4902.08, 5869.51, 6832.89, 7792.33, 8747.95, 9699.84 },
                                      { 494.69, 988.19, 1971.72, 2950.71, 3925.30, 4895.60, 5861.76, 6823.86, 7782.04, 8736.40, 9687.84 },
                                      { 494.03, 986.88, 1969.11, 2946.81, 3920.10, 4889.13, 5854.01, 6814.84, 7771.76, 8724.85, 9674.23 }
                                    };

            for (int l_idx = 0; l_idx < ref_lats_deg.Length; l_idx++)
            {
                for (int p_idx = 0; p_idx < ref_pressures.Length; p_idx++)
                {
                    double g = PHX.Gravity_constant_wgs84_calc(ref_lats_deg[l_idx]);

                    // taking the water density in the midpoint to consider the compression of water
                    double rho = PHX.Water_density_calc(ref_temperature, ref_pressures[p_idx] / 2.0, ref_salinity);
                    double h_est = PHX.Depth_by_pressure_calc(ref_pressures[p_idx], 0.0, rho, g);

                    // calculated values should deviate less than 0.07% of reference values
                    Assert.AreEqual(h_est, ref_dpt_lat[l_idx, p_idx], ref_dpt_lat[l_idx, p_idx] * 0.0007);
                }
            }
        }

        [TestMethod]
        public void Water_fpoint_calc_test()
        {
            // Reference values according to 
            // Algorithms for computation of fundamental properties of seawater. 
            // Unesco technical papers in marine science vol. 44, 1983, pp. 30
            // https://darchive.mblwhoilibrary.org/bitstream/handle/1912/2470/059832eb.pdf

            Assert.AreEqual(PHX.Water_fpoint_calc(50000.0, 40.0), -2.588567, 1E-6);

            double[] ref_pressures = new double[] { 0.0, 1E4, 2E4, 3E4, 4E4, 5E4 };
            double[] ref_salinities = new double[] { 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0 };
            double[,] ref_tf = new double[,] { { -0.274, -0.349, -0.424, -0.500, -0.575, -0.650 },
                                           { -0.542, -0.618, -0.693, -0.768, -0.844, -0.919 },
                                           { -0.812, -0.887, -0.962, -1.038, -1.113, -1.188 },
                                           { -1.083, -1.159, -1.234, -1.309, -1.384, -1.460 },
                                           { -1.358, -1.434, -1.509, -1.584, -1.660, -1.735 },
                                           { -1.638, -1.713, -1.788, -1.864, -1.939, -2.014 },
                                           { -1.922, -1.998, -2.073, -2.148, -2.224, -2.299 },
                                           { -2.212, -2.287, -2.363, -2.438, -2.513, -2.589 }
                                         };

            for (int s_idx = 0; s_idx < ref_salinities.Length; s_idx++)
            {
                for (int p_idx = 0; p_idx < ref_pressures.Length; p_idx++)
                {
                    double tf_est = PHX.Water_fpoint_calc(ref_pressures[p_idx], ref_salinities[s_idx]);
                    Assert.AreEqual(tf_est, ref_tf[s_idx, p_idx], 0.001);
                }
            }
        }
    }
}
