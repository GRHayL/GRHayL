
//********************************************
// Basic definitions for module finite_difference:

// Set the number of ghost zones
// Note that upwinding in e.g., BSSN requires that NGHOSTS = FD_CENTDERIVS_ORDER/2 + 1 <- Notice the +1.
#define NGHOSTS 3
//********************************************

// Step 0b: Set up numerical grid structure, first in space...
static const int Nxx = 8;
// const REAL domain_size = 10.0;
static const double dxx = 2.5;
static const double invdx = 1.0/dxx;
static const int Nxx0 = Nxx;
static const int Nxx1 = Nxx;
static const int Nxx2 = Nxx;
static const int Nxx_plus_2NGHOSTS0 = Nxx + 2*NGHOSTS;
static const int Nxx_plus_2NGHOSTS1 = Nxx + 2*NGHOSTS;
static const int Nxx_plus_2NGHOSTS2 = Nxx + 2*NGHOSTS;


//********************************************
// Basic definitions for module reference_metric:
//********************************************


//********************************************
// Basic definitions for module grid:
// EVOLVED VARIABLES:
#define NUM_EVOL_GFS 5
#define STILDED0GF	0
#define STILDED1GF	1
#define STILDED2GF	2
#define RHO_STARGF	3
#define TAU_TILDEGF	4


// AUXILIARY VARIABLES:
#define NUM_AUX_GFS 0


// AUXEVOL VARIABLES:
#define NUM_AUXEVOL_GFS 102
#define BU0GF	0
#define BU1GF	1
#define BU2GF	2
#define BLU0GF	3
#define BLU1GF	4
#define BLU2GF	5
#define BRU0GF	6
#define BRU1GF	7
#define BRU2GF	8
#define HLLE_FLUX_STILDED0GF	9
#define HLLE_FLUX_STILDED1GF	10
#define HLLE_FLUX_STILDED2GF	11
#define HLLE_FLUX_RHO_STARGF	12
#define HLLE_FLUX_TAU_TILDEGF	13
#define KDD00GF	14
#define KDD01GF	15
#define KDD02GF	16
#define KDD11GF	17
#define KDD12GF	18
#define KDD22GF	19
#define PGF	20
#define P_LGF	21
#define P_RGF	22
#define VU0GF	23
#define VU1GF	24
#define VU2GF	25
#define VLU0GF	26
#define VLU1GF	27
#define VLU2GF	28
#define VRU0GF	29
#define VRU1GF	30
#define VRU2GF	31
#define ALPHAGF	32
#define ALPHA_DD0GF	33
#define ALPHA_DD1GF	34
#define ALPHA_DD2GF	35
#define ALPHA_FACEGF	36
#define BETAU0GF	37
#define BETAU1GF	38
#define BETAU2GF	39
#define BETAU_DD00GF	40
#define BETAU_DD01GF	41
#define BETAU_DD02GF	42
#define BETAU_DD10GF	43
#define BETAU_DD11GF	44
#define BETAU_DD12GF	45
#define BETAU_DD20GF	46
#define BETAU_DD21GF	47
#define BETAU_DD22GF	48
#define BETA_FACEU0GF	49
#define BETA_FACEU1GF	50
#define BETA_FACEU2GF	51
#define CS2_LGF	52
#define CS2_RGF	53
#define GAMMADD00GF	54
#define GAMMADD01GF	55
#define GAMMADD02GF	56
#define GAMMADD11GF	57
#define GAMMADD12GF	58
#define GAMMADD22GF	59
#define GAMMADD_DD000GF	60
#define GAMMADD_DD001GF	61
#define GAMMADD_DD002GF	62
#define GAMMADD_DD010GF	63
#define GAMMADD_DD011GF	64
#define GAMMADD_DD012GF	65
#define GAMMADD_DD020GF	66
#define GAMMADD_DD021GF	67
#define GAMMADD_DD022GF	68
#define GAMMADD_DD110GF	69
#define GAMMADD_DD111GF	70
#define GAMMADD_DD112GF	71
#define GAMMADD_DD120GF	72
#define GAMMADD_DD121GF	73
#define GAMMADD_DD122GF	74
#define GAMMADD_DD220GF	75
#define GAMMADD_DD221GF	76
#define GAMMADD_DD222GF	77
#define GAMMA_FACEDD00GF	78
#define GAMMA_FACEDD01GF	79
#define GAMMA_FACEDD02GF	80
#define GAMMA_FACEDD11GF	81
#define GAMMA_FACEDD12GF	82
#define GAMMA_FACEDD22GF	83
#define HGF	84
#define H_LGF	85
#define H_RGF	86
#define RHOBGF	87
#define RHOB_LGF	88
#define RHOB_RGF	89
#define U4U0GF	90
#define U4U1GF	91
#define U4U2GF	92
#define U4U3GF	93
#define U4LU0GF	94
#define U4LU1GF	95
#define U4LU2GF	96
#define U4LU3GF	97
#define U4RU0GF	98
#define U4RU1GF	99
#define U4RU2GF	100
#define U4RU3GF	101

// Declare the IDX4S(gf,i,j,k) macro, which enables us to store 4-dimensions of
//   data in a 1D array. In this case, consecutive values of "i"
//   (all other indices held to a fixed value) are consecutive in memory, where
//   consecutive values of "j" (fixing all other indices) are separated by
//   Nxx_plus_2NGHOSTS0 elements in memory. Similarly, consecutive values of
//   "k" are separated by Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1 in memory, etc.
#define IDX4S(g,i,j,k)                                                  \
  ( (i) + Nxx_plus_2NGHOSTS0 * ( (j) + Nxx_plus_2NGHOSTS1 * ( (k) + Nxx_plus_2NGHOSTS2 * (g) ) ) )
#define IDX4ptS(g,idx) ( (idx) + (Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2) * (g) )
#define IDX3S(i,j,k) ( (i) + Nxx_plus_2NGHOSTS0 * ( (j) + Nxx_plus_2NGHOSTS1 * ( (k) ) ) )
#define LOOP_REGION(i0min,i0max, i1min,i1max, i2min,i2max)              \
  for(int i2=i2min;i2<i2max;i2++) for(int i1=i1min;i1<i1max;i1++) for(int i0=i0min;i0<i0max;i0++)


static const int kronecker_delta[4][3] = { { 0,0,0 },
                                    { 1,0,0 },
                                    { 0,1,0 },
                                    { 0,0,1 } };

void initialize_structs(const int idx,
                        const double *restrict auxevol_gfs, 
                        primitive_quantities *restrict prims, 
                        primitive_quantities *restrict prims_r, 
                        primitive_quantities *restrict prims_l, 
                        metric_quantities *restrict metric, 
                        metric_quantities *restrict metric_face,
                        extrinsic_curvature *restrict curv, 
                        metric_derivatives *restrict metric_derivs);

void read_from_binary_file_all(const char *restrict binary_file, 
                               double *restrict auxevol_gfs);
void read_from_binary_file_recons(const char *restrict binary_file, 
                                  double *restrict auxevol_gfs);

