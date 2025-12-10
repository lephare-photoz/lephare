#include <valarray>
#include <array>

namespace {

#define NUM_EMISSION_LINES 65
#define NUM_NEBULAR_DATA 71
// proportion of Helium compared to hydrogen = n(HeII)/n(HI)
#define n_heII 0.1
// [cm^3 s^-1] : total recombination coeff for
// hydrogen in case B (except to groundstate), for Te = 10kK
//  Different from Schearer, use Osterbrock
#define alpha_B 2.59e-13

/*! Rest frame wavelength of the emission lines considered in LePHARE
'H_Lya', '[CII]', '[OIII]', '[CIII]', '[NII]', '[CII]',
'[MgII]', '[OII]a', '[OII]b', '[NeIII]', 'H_zeta+[HeI]', 'H_eps+[NeIII]',
'[HeI]', '[SII]', '[SII]', 'H_delta', '[OIII]', 'H_gamma', '[HeI]',
'[ArIV]+[HeI]', 'H_beta', '[OIII]', '[OIII]', '[NI]', '[NII]', '[HeI]',
'[OI]', '[SIII]', '[NII]', 'H_alpha', '[NII]', '[HeI]', '[SII]', '[SII]',
'[HeI]', '[ArIII]', '[OII]', '[OII]', '[ArIII]', '[SIII]', '[SIII]',
'[SII]', '[SII]', '[SII]', '[HeI]', 'Pa_Pa', 'Pb_Pa', 'Pg_Pa', 'Pd_Pa',
'Pe_Pa', 'P6_HI', 'P7_HI', 'P8_HI', 'P9_HI', 'P10_HI', 'P11_HI', 'P12_HI',
'Bra', 'Brb', 'Brg', 'Brd', 'Bre', 'Br6', 'Br7', 'Br8'
*/
const std::array<double, NUM_EMISSION_LINES> emission_lines = {
    1215.67,   1335.708,  1663.48,   1908.73,   2140.896,  2325.726,  2799.116,
    3726.04,   3728.8,    3869.063,  3888.856,  3968.771,  4026.191,  4068.601,
    4076.349,  4101.734,  4363.21,   4340.472,  4471.48,   4712.249,  4861.352,
    4958.911,  5006.843,  5199.079,  5754.594,  5875.62,   6300.304,  6312.064,
    6548.041,  6562.787,  6583.461,  6678.151,  6716.44,   6730.812,  7065.19,
    7135.803,  7319.458,  7330.196,  7751.057,  9068.611,  9530.586,  10286.731,
    10320.491, 10336.407, 10827.334, 18750.976, 12818.072, 10938.174, 10049.368,
    9545.969,  9229.667,  9014.909,  8862.782,  8750.473,  8665.018,  8598.392,
    8545.383,  40511.742, 26251.517, 21655.267, 19445.581, 18174.237, 17362.143,
    16806.509, 16407.21};

// In the empirical GalSED::generateEmEmpUV ratio estimation, we only
// consider Lyman_alpha, OIIa and OIIb; Hbeta, OIII doublet et Halpha.
// Hbeta is taken as normalization.
// published in a slightly different form in Ilbert et al. 2009
const std::array<double, NUM_EMISSION_LINES> empirical_ratio = {
    22.2, 0, 0, 0, 0,    0,    0, 0.81, 0.81, 0, 0, 0, 0,    0,  0, 0, 0,
    0,    0, 0, 1, 0.21, 0.59, 0, 0,    0,    0, 0, 0, 2.90, 0., 0, 0, 0,
    0,    0, 0, 0, 0,    0,    0, 0,    0,    0, 0, 0, 0,    0,  0, 0, 0,
    0,    0, 0, 0, 0,    0,    0, 0,    0,    0, 0, 0, 0,    0};

// conversion factor between Halpa and other lines
// Use a ratio 1 for Halpha/OII, according to Iary. Drop the 1.77 since
// attenuation seems to be alreay included I took OIII5007/OIII4959=3 I
// took 4.081 betwen 5007 and Hbeta, as the physical recipes.
// This ratio will be modified later as a function of redshift.
// For NII/Halpha, take 0.3 according to the BPT diagram.
// For Lyman_alpha, use the same factor as for
// physical recipes since fescape fraction is applied.
const std::array<double, NUM_EMISSION_LINES> empirical_ratio2 = {
    22.2, 0, 0, 0, 0,    0,     0, 1.425, 1.425, 0, 0, 0, 0,    0,    0, 0, 0,
    0,    0, 0, 1, 1.36, 4.081, 0, 0,     0,     0, 0, 0, 2.85, 0.86, 0, 0, 0,
    0,    0, 0, 0, 0,    0,     0, 0,     0,     0, 0, 0, 0,    0,    0, 0, 0,
    0,    0, 0, 0, 0,    0,     0, 0,     0,     0, 0, 0, 0,    0};

/* Line ratios for 3 different metallicities:
 * Z1 (Z < 0.02),
 * Z2 (0.02 < Z <= 0.2),
 * Z3-Z5 (0.2 < Z)
 * BC03 convention: 62 solar 0.02, 52 subsolar 0.008, 42 subsolar 0.004
 * For Z < 0.03 in solar unit, Z1=0.0004
 */
const std::array<double, NUM_EMISSION_LINES> Z1_line = {
    22.2,    0.0,     0.0,     0.0,    0.0,    0.0,    0.0,    0.245,  0.245,
    0.295,   0.203,   0.27,    0.015,  0.005,  0.002,  0.256,  0.109,  0.466,
    0.036,   0.01,    1.0,     1.097,  3.159,  0.003,  0.0,    0.096,  0.008,
    0.009,   0.005,   2.87,    0.015,  0.026,  0.037,  0.029,  0.028,  0.027,
    0.012,   0.007,   0.067,   0.0,    0.0,    0.0,    0.0,    0.0,    0.232,
    0.352,   0.165,   0.0906,  0.0368, 0.0366, 0.0254, 0.0184, 0.0138, 0.0106,
    0.0185,  0.00666, 0.00541, 0.0884, 0.0471, 0.0281, 0.0186, 0.0127, 0.00914,
    0.00683, 0.00524};
// 0.03 < Z <= 0.3 in solar unit, Z2=0.004
const std::array<double, NUM_EMISSION_LINES> Z2_line = {
    22.2,    0.0,     0.058,   0.0,    0.0,    0.0,    0.31,   0.896,  0.896,
    0.416,   0.192,   0.283,   0.015,  0.017,  0.007,  0.256,  0.066,  0.466,
    0.036,   0.014,   1.0,     1.617,  4.752,  0.01,   0.0,    0.108,  0.041,
    0.017,   0.059,   2.87,    0.175,  0.03,   0.188,  0.138,  0.023,  0.071,
    0.027,   0.014,   0.176,   0.2,    0.51,   0.0,    0.0,    0.0,    0.225,
    0.352,   0.165,   0.0906,  0.0368, 0.0366, 0.0254, 0.0184, 0.0138, 0.0106,
    0.0185,  0.00666, 0.00541, 0.0884, 0.0471, 0.0281, 0.0186, 0.0127, 0.00914,
    0.00683, 0.00524};
// 0.3 < Z in solar unit, Z3-Z5=0.008-0.05
const std::array<double, NUM_EMISSION_LINES> Z3_line = {
    22.2,    0.11,    0.01,    0.18,   0.01,   0.29,   0.07,   1.505,  1.505,
    0.3,     0.107,   0.159,   0.015,  0.029,  0.011,  0.256,  0.01,   0.466,
    0.05,    0.0,     1.0,     1.399,  4.081,  0.03,   0.01,   0.14,   0.13,
    0.03,    0.136,   2.87,    0.404,  0.03,   0.3,    0.21,   0.04,   0.035,
    0.026,   0.014,   0.086,   0.365,  0.945,  0.048,  0.058,  0.054,  0.229,
    0.352,   0.165,   0.0906,  0.0368, 0.0366, 0.0254, 0.0184, 0.0138, 0.0106,
    0.0185,  0.00666, 0.00541, 0.0884, 0.0471, 0.0281, 0.0186, 0.0127, 0.00914,
    0.00683, 0.00524};

// to be used in SED::applyExtLines. Calzetti value would be 0.44 everywhere
const std::array<double, NUM_EMISSION_LINES> ebvFac = {
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

// Data for the contribution of the continuum from the nebular regions
//
//  Wavelength grid for wavelength & emission coefficient :
//  Points from Aller (1984)   for lambda=1000 to 10000.1 (54 pts)
//  Points from Ferland (1980) for lambda > 10000.1     (17 pts)
const std::valarray<double> ga_lamb = {
    1000.,  1200.,  1300.,  1400.,  1500.,  1600.,   1800.,   2051.,  2053.,
    2200.,  2400.1, 2599.4, 2600.8, 2725.4, 2855.2,  2997.9,  3121.4, 3122.,
    3331.,  3421.4, 3422.,  3527.,  3642.,  3648.,   3679.2,  3679.9, 4000.,
    4282.8, 4499.9, 4996.5, 5096.,  5450.8, 5695.8,  5700.,   5995.9, 6633.8,
    6635.8, 6999.9, 7438.7, 7441.1, 7848.,  7850.4,  8193.3,  8196.2, 8196.7,
    8198.5, 8207.,  8209.,  8265.4, 8268.1, 8499.9,  9000.,   9500.1, 10000.1,
    14591., 14592., 22799., 22800., 32805., 32806.,  44705.,  44706., 58460.,
    58461., 74145., 74146., 91199., 91200., 109878., 109879., 132173.};
/*  Emission coefficients from Aller (1984, Table 4-9, p.102)
    for lambda <=10000.1 Ang,
    and Ferland (1980, PASP 92, 506, Table 1) for lambda >=14591 Ang
    all for Te = 10000 K: gamma_H in [1.d-40 ergs cm^3 s^-1 Hz^-1]
*/
// Hydrogen recombination
const std::valarray<double> ga_H_val = {
    0.001,  0.009,  0.023,  0.051,  0.100,  0.181,  0.486,  1.272,  1.280,
    2.026,  3.453,  5.404,  5.419,  6.928,  8.740,  11.020, 13.231, 13.242,
    17.477, 19.490, 19.505, 21.976, 24.841, 1.387,  1.434,  1.435,  1.950,
    2.461,  2.883,  3.929,  4.150,  4.959,  5.534,  5.544,  6.252,  7.801,
    7.806,  8.693,  9.754,  9.760,  10.729, 10.735, 11.537, 11.544, 11.545,
    11.549, 11.569, 4.317,  4.368,  4.371,  4.578,  5.020,  5.449,  5.867,
    9.04,   5.93,   8.51,   6.90,   8.50,   7.56,   8.66,   8.06,   8.87,
    8.47,   9.11,   8.82,   9.34,   9.14,   9.58,   9.42,   9.80};

/*   gamma_2q in [1.d-40 ergs cm^3 s^-1 Hz^-1]
     assume ga_2q=0 for lambda >=10000 A
*/
// To be check
const std::valarray<double> ga_2q_val = {
    0.000, 0.000, 5.611, 8.236, 9.210, 9.478, 9.196, 8.433, 8.426, 7.947, 7.314,
    6.747, 6.743, 6.423, 6.096, 5.770, 5.508, 5.507, 5.206, 4.947, 4.946, 4.771,
    4.588, 4.579, 4.532, 4.531, 4.088, 3.746, 3.512, 3.055, 2.975, 2.715, 2.554,
    2.552, 2.375, 2.064, 2.063, 1.913, 1.751, 1.750, 1.616, 1.615, 1.515, 1.514,
    1.514, 1.514, 1.511, 1.511, 1.496, 1.496, 1.438, 1.324, 1.221, 1.129, 0.,
    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,
    0.,    0.,    0.,    0.,    0.};
/*  gamma_HeI in [1.d-40 ergs cm^3 s^-1 Hz^-1]
    assume ga_hei=ga_h for lambda >=10000 A
*/
// Helium recombination
const std::valarray<double> ga_HeI_val = {
    0.046,  0.079,  0.126,  0.210,  0.400,  0.869,  1.650,  3.895,  3.926,
    6.087,  9.014,  11.861, 7.156,  9.976,  12.651, 15.326, 17.443, 15.873,
    20.979, 22.989, 4.321,  4.902,  5.501,  5.532,  5.687,  1.450,  2.070,
    2.807,  3.310,  4.943,  5.231,  6.176,  6.759,  6.771,  7.583,  9.086,
    8.661,  9.391,  10.176, 10.000, 10.607, 8.081,  8.497,  5.880,  5.880,
    5.041,  5.048,  5.050,  5.098,  4.360,  4.596,  5.065,  5.483,  5.860,
    9.04,   5.93,   8.51,   6.90,   8.50,   7.56,   8.66,   8.06,   8.87,
    8.47,   9.11,   8.82,   9.34,   9.14,   9.58,   9.42,   9.80};

const std::valarray<double>
ga_total = ga_H_val + ga_2q_val + 0.1 * ga_HeI_val;
//ga_total = ga_H_val + ga_2q_val + n_heII * ga_HeI_val;

}  // namespace
