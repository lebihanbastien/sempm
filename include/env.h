#ifndef DEFINE_ENV_H_INCLUDED
#define DEFINE_ENV_H_INCLUDED

/**
 * \file env.h
 * \brief Routines for the definition of the working environment via several C structure
 *        that represent Three or Four-Body models of the Sun-Earth-Moon system.
 *        The main routine is init_fbp_lib, which initializes a FBPL structure, i.e. a
 *        structure that contains all the necessary constants and coefficients
 *        (equations of motion, Hamiltonian) for a given Four-Body problem
 *        around a given libration point. The coefficients are retrieved from txt files
 *        stored in subfolders of ./data/
 *        that must have been computed with the routine qbtbp().
 *
 *        The physical constants hard-coded in this file have been taken from:
 *          - Andreu 1998 for the Sun-Earth-Moon system,
 *          - JPL/Goddard Space Flight Center websites for the rest of the solar system.
 *
 *        Note that most of the routines initialize C structures with C routines,
 *        because they were taken from original C code made in 2014.
 *        Most of these C structures could be implemented again as C++ classes with
 *        proper constructors. This would requires quite a bit of work because
 *        most of these structures are used on the fly in the current code.
 *        This file could also be merged with init.cpp in the long run, in order
 *        to have one single file for the initialization routines.
 *
 * \author BLB
 */

//Std
#include <iostream>
#include <fstream>
#include <sstream>

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <string.h>

//Gsl
#include <gsl_complex_math.h>

//Custom
#include "ofs.h"
#include "config.h"
#include "constants.h"

//----------------------------------------------------------------------------------------
//            Structures
//----------------------------------------------------------------------------------------
/**
 * \struct LibPoint
 * \brief Structure to describe a Libration point.
 **/
typedef struct LibPoint LibPoint;
struct LibPoint
{
    //number
    int number;
    //position [adim] in the corresponding CR3BP
    double position[3];
    //distance to closest primary (only for l1,l2,l3);
    double gamma_i;
    //energy in the proper CR3BP
    double Ei;
    //jacobi constant in the proper CR3BP
    double Ci;
};

/**
 * \struct Body
 * \brief Characteristic constants of a celestial body (Mass, orbital period...)
 **/
typedef struct Body Body;
struct Body
{
    //Mass [kg]
    double M;
    //Gravitationnal parameter [km^3/s^2]
    double GM;
    //Equatorial Radius [km]
    double Req;
    //Mean radius [km]
    double Rm;
    //Sidereal orbit period [s]
    double T;
    //Semi-major axis [km]
    double a;
    //Name
    char name[50];
};

/**
 * \struct CR3BP
 * \brief  Environment of a given Circular Restricted Three-Body Problem.
 **/
typedef struct CR3BP CR3BP;
struct CR3BP
{
    Body m1;   //First primary
    Body m2;   //Second primary

    double mu; // Gravitational constant [-]
    double L;  // Distance parameter [km]
    double T;  // Time parameter [s]
    double R1; // Radius of the first primary  [km]
    double R2; // Radius of the second primary [km]
    double rh; // Hill's radius [-]

    //Libration points (adim)
    LibPoint l1;
    LibPoint l2;
    LibPoint l3;
    LibPoint l4;
    LibPoint l5;

    //Name
    char name[50];
};

/**
 * \struct USYS
 * \brief  Unit system structure containing a given set of parameters in a given unit system.
 **/
typedef struct USYS USYS;
struct USYS
{
    //Note: a name should be added!
    int label;      //type of unit system
    double mu_EM;   //mass ratio (EM)
    double mu_SEM;  //mass ratio (SEM)
    double mu_SE;   //mass ratio (SE)
    double as;      //radius of the circular motion approximating the quasi-circular exterior (Earth+Moon-Sun) motion
    double ns;      //Mean angular motion of this circular motion
    double ai;      //radius of the circular motion approximating the quasi-circular interior (Earth-Moon) motion
    double ni;      //Mean angular motion of this circular motion
    double n;       //difference of mean angular motion
    double T;       //Period
    double ms;      //Sun mass
    double me;      //Earth mass
    double mm;      //Moon mass
    double lecc;    //lunar eccentricity
};

/**
 * \struct CSYS
 * \brief  Coordinates system structure for a given Four-Body problem about a given
 *         libration point. It mainly contains some constants that depend on the libration
 *         point (gamma, c1) and the coefficients of the equations of motion.
 *         It also contains some useful strings, such as the subfolders where the data
 *         is stored. Finally, for simplicity, the CSYS structure also contains
 *         information on the type of parameterization of a certain manifold,
 *         as currently desired by the user. For the latter data, a dedicated structure
 *         could be designed.
 **/
typedef struct CSYS CSYS;
struct CSYS
{
    //Name and label
    string name;    //name of the unit system
    int label;      //type of unit cs

    //Associated CR3BP
    CR3BP cr3bp;

    //default unit system
    USYS us;

    //Simple coefficients
    double mu;      //mass ratio of the associated CR3BP
    double c1;      //c1 coefficient, defined wrt to a given libration point li
    double c2;      //c2 coefficient, defined wrt to a given libration point li
    double gamma;   //gamma coefficient, defined wrt to a given libration point li
    int li;         //associated libration point
    int man_type;    //associated manifold type

    //Arrays
    double* coeffs; //Default set of vector field coefficients
    double* Ps;     //Sun   position in EM coordinates
    double* Pe;     //Earth position in EM coordinates
    double* Pm;     //Moon  position in EM coordinates
    double* ps;     //Sun   position in NC coordinates
    double* pe;     //Earth position in NC coordinates
    double* pm;     //Moon  position in NC coordinates

    //Main folders
    string F_COC;   //Change of coordinates Translation+Floquet+Complexification
    string F_PRINT; //For printing results
    string F_COEF;  //Coefficients alphas for integration in the Earth-Moon system.

    //Parameterization folders
    string F_PMS;   //Global Parameterization folder, chosen within the folders below:
    string F_GS;    //Graph style           |
    string F_NF;    //Normal form style     |  all for Center Manifold
    string F_MS;    //Mixed style           |

    string F_CS;    //Center-Stable
    string F_CU;    //Center-Unstable
    string F_CUS;   //Center-Hyperbolic

    //QBTBP: solution of the Quasi-Bicircular Three-Body Problem, in the proper units
    Ofsc zt;    //inner motion
    Ofsc Zt;    //outer motion
    Ofsc ztdot; //time derivative of the inner motion
    Ofsc Ztdot; //time derivative of the outer motion
};

/**
 * \struct FBP
 * \brief  Structure that describes a given Four-Body Problem.
 **/
typedef struct FBP FBP;
struct FBP
{
    Body m1;   //First  primary
    Body m2;   //Second primary
    Body m3;   //Third  primary

    CR3BP cr3bp1; //First associated CR3BP
    CR3BP cr3bp2; //Second associated CR3BP
};


/**
 * \struct FBPL
 * \brief  Structure for a given Four-Body Problem focused on one libration point.
 *         The libration point must be L1, L2 (should also work for L3 in some cases).
 *
 *         The FBPL structure contains several variables, including:
 *
 *          - li, the current libration point at hand (1 or 2).
 *
 *          - us_em, a USYS (unit system) structure. It contains several
 *          constants (distances, frequencies) in Earth-Moon  normalized units.
 *          - us_sem, a USYS (unit system) structure. It contains several
 *          constants (distances, frequencies) in Sun-Earth normalized units.
 *
 *          - cs_em_l1, a CSYS (coordinate system) structure. It contains several
 *          constants in Earth-Moon normalized units that are associated with the EML1
 *          point (gamma, c1), as well as the coefficients of the vector field of
 *          the Four-Body Problem, both in EM and EMNC coordinates.
 *          - cs_em_l2, cs_sem_l1, cs_sem_l2, equivalents structures
 *          for the other points.
 *
 *          Moreover, the USYS and CSYS structures are duplicated for easy access:
 *          - us, equal one of the previous USYS structures:
 *                  * us_em if coord_sys == Csts::EM,
 *                  * us_Sem if coord_sys == Csts::SEM.
 *          - cs_em, equal to cs_em_l1,2 if li_EM == 1,2.
 *          - cs_sem, equal to cs_sem_l1,2 if li_SEM == 1,2.
 *          - cs, equal to
 *                  * cs_em_l1,2 if li_EM == 1,2 & coord_sys == Csts::EM,
 *                  * cs_sem_l1,2 if li_SEM == 1,2 & coord_sys == Csts::SEM.
 *
 *         See code below for the other variables.
 *
 *         Note that there are obviously some unecessary variables that have
 *         been gathered during the iterative implementation of the present code. A proper
 *         cleaning of this mega structure is probably required.
 *
 **/
typedef struct FBPL FBPL;
struct FBPL
{
    //Libration points at hand
    int li_EM;     //number of the EM  libration point considered (1 or 2)
    int li_SEM;    //number of the SEM libration point considered (1 or 2)

    //Model
    int model;     // Current model Csts::QBCP, Csts::CRTBP, etc.
    int coord_sys;  // Current framework: Either Csts:EM or Csts::SEM
    int li;        // Current libration point: Either li_EM or li_SEM
    int param_style;       // Current style of parameterization

    //Unit systems
    USYS us_em;   //EM unit system: constants in EM normalized units
    USYS us_sem;  //SEM unit system: constants in SEM normalized units
    USYS us;      //default unit system

    //Coordinate systems (constants and coefficients of the equations of motion)
    CSYS cs_em_l1;      //EM around L1
    CSYS cs_em_l2;      //EM around L2
    CSYS cs_em_l3;      //EM around L3
    CSYS cs_sem_l1;     //SEM around L1
    CSYS cs_sem_l2;     //SEM around L2
    CSYS cs_sem_l3;     //SEM around L3

    //Default coordinate systems
    CSYS cs_em;     //Default EM csys  (for COCs)
    CSYS cs_sem;    //Default SEM csys (for COCs)
    CSYS cs;        //Default csys     (for integration)

    // Misc
    int n_order_fourier;            //Order of the Fourier expansions
    int eff_nf;        //Effective order of the Fourier expansions, for some specific computations
    int is_norm;        //Are the equations of motion normalized?
    int numberOfCoefs; //number of coefficients in the equations of motion

    //6*6 matrix B used for Floquet transformation
    double* B;
};

/**
 * \struct QBCP_I
 * \brief  Hybrid environment for contination procedures between two models (model1 and model2).
 *         In practice, the effective environment is
 *                     "model_eff =  (1.0 - epsilon)*model1 + epsilon*model2"
 *
 **/
typedef struct QBCP_I QBCP_I;
struct QBCP_I
{
    FBPL model1;
    FBPL model2;
    double epsilon;
};

//----------------------------------------------------------------------------------------
//            Initialization routines
//----------------------------------------------------------------------------------------

/**
 * \brief Initialize a FBPL structure, i.e. a FBP (Four-Body Problem) focused on two
 *        libration points: one the EM system and one of the SEM system.
 *        As for now, the libration point must be L1 or L2 for both systems.
 *        The FBPL structure fbpl contains several variables, including:
 *
 *          - fbpl.li, the current libration point at hand (1 or 2).
 *
 *          - fbpl.us_em, a USYS (unit system) structure. It contains several
 *          constants (distances, frequencies) in Earth-Moon  normalized units.
 *          - fbpl.us_sem, a USYS (unit system) structure. It contains several
 *          constants (distances, frequencies) in Sun-Earth normalized units.
 *
 *          - fbpl.cs_em_l1, a CSYS (coordinate system) structure. It contains several
 *          constants in Earth-Moon normalized units that are associated with the EML1
 *          point (gamma, c1), as well as the coefficients of the vector field of
 *          the Four-Body Problem, both in EM and EMNC coordinates.
 *          - fbpl.cs_em_l2, cs_sem_l1, cs_sem_l2, equivalents structures
 *          for the other points.
 *
 *          Moreover, the USYS and CSYS structures are duplicated for easy access:
 *          - fbpl.us, equal one of the previous USYS structures:
 *                  * us_em if coord_sys == Csts::EM,
 *                  * us_Sem if coord_sys == Csts::SEM.
 *          - fbpl.cs_em, equal to cs_em_l1,2 if li_EM == 1,2.
 *          - fbpl.cs_sem, equal to cs_sem_l1,2 if li_SEM == 1,2.
 *          - fbpl.cs, equal to
 *                  * cs_em_l1,2 if li_EM == 1,2 & coord_sys == Csts::EM,
 *                  * cs_sem_l1,2 if li_SEM == 1,2 & coord_sys == Csts::SEM.
 *
 *        In the end, we mainly use fbpl.li (libration point), fbpl.us (constants in the
 *        suitable normalized units), and fbpl.cs (constants and coefficients about the
 *        proper libration point).
 *        See the definition of the FBPL structure to see the other variables contained
 *        in it.
 *
 * Note that the FBPL is focused on one libration point but contains the constants and
 * coefficients for all four points (EML1, EML2, SEL1, SEL2).
 * The focus of the FBPL structure can be changed from one point to another via
 * the routines change_coord and change_li_coord. Moreover, the coefficients for each
 * point can be accessed via calls of the form fbpl->cs_em_l1, fbpl->cs_em_l2, etc.
 *
 * The inputs of this routines are:
 * \param fbpl a pointer on the FBPL structure to initialize.
 * \param fbp a pointer on the FBP parent structure.
 * \param li_EM number of the libration point for the EM system.
 * \param li_SEM number of the libration point for the SEM system.
 * \param model: either Csts::QBCP, Csts::BCP, Csts::CRTBP, etc.
 * \param coord_sys: default coordinate system for this structure:
 *           - if coord_sys == Csts::EM,  the fbpl is focused on the li_EM  point of the
 *             EM system.
 *           - if coord_sys == Csts::SEM, the fbpl is focused on the li_SEM point of the
 *             SEM system.
 *        The focus can be change dynamically during computation, via the routines
 *        change_coord and change_li_coord.
 * \param param_style: type of parameterization of the manifolds (Csts::GRAPH, etc). Note that
 *        the param_style influences the number of coefficients taken into account in the
 *        Fourier series. Indeed, for graph method, the reduced vector field is non
 *        autonomous, and full Fourier series are used. For normal form, the reduced
 *        vector field is quasi autonomous and we can safely reduce the order of the
 *        series to 5 (11 coefficients taken into account).
 * \param man_type_EM: type of manifold about li_EM (Csts::MAN_CENTER, etc).
 * \param man_type_SEM: type of manifold about li_SEM.
 * \param is_new an integer: equal to 1 if no solution has been previously
 *        computed with the routine qbtbp(), 0 otherwise.
 * \param is_norm: are the equations of motion normalized?  Probably deprecated,
 *        should be always true. Kept for consistency with older code.
 *
 *
 * Note that the FBP structure fbp is used only for the initialization of the coordinate
 * systems. More precisely, it contains some parameters specific to each libration point
 * (gamma, c1, etc), via its CR3BP structures
 * (see init_fbp, the routine that initializes the QBCP structures).
 *
 * Note also that there is no speficic need to use two libration points at the
 * same time (one of the EM and one of the SEM systems). This choice was made with the
 * study of li_EM-li_SEM connections in mind, which are now computed in another package.
 *
 **/
void init_fbp_lib(FBPL* fbpl, FBP* fbp, int li_EM, int li_SEM,
               int model, int coord_sys, int param_style, int man_type_EM, int man_type_SEM,
               int is_new, int is_norm);

/**
 * \brief Initializes a unit system in the form of a USYS structure,
 *        such as Earth-Moon or Sun-(Earth+Moon) unit systems.
 *        Consistent with appendices of the PhD manuscript (October 2017).
 *        The Earth-Moon constants are taken from Andreu 1998.
 * \param usys pointer on the USYS structure to initialize.
 * \param label the type of unit system (Csts:EM, Csts:SEM, etc).
 * \param model the type of model (Csts::QBCP, etc).
 **/
void init_usys(USYS* usys, int label, int model);

/**
 * \brief Initializes a coordinate systems (CSYS structure), with associated
 *        vector field coefficients, data folder names, and unit system.
 * \param csys pointer on the CSYS structure to initialize.
 * \param fbpl pointer on the FBPL structure that contains csys.
 * \param fbp pointer on the FBP structure that contains parameters specific to
 *        each libration points (namely, the value of gamma)
 * \param li number of the libration point to focus on (L1, L2).
 * \param coord_sys indix of the coordinate system to use (Csts::EM, Csts::SEM).
 * \param man_type: type of manifold about li (Csts::MAN_CENTER, Csts::MAN_CENTER_S...).
 * \param is_new boolean. if true, the qbtbp has not been computed via the qbtbp() routine,
 *        so the vector field coefficients is not initialized.
 *
 *   Note that the FBP structure is used only for the initialization of the coordinate
 *   systems. More precisely, it contains some parameters specific to each libration point
 *  (gamma), via its CR3BP structures.
 **/
void init_csys(CSYS* csys, FBPL* fbpl, FBP* fbp, int li,
               int coord_sys, int man_type, int is_new);

/**
 * \brief Initialize a Four-Body Problem in the form of a FBP structure.
 * \param fbp pointer on the FBP structure to init.
 * \param first name of the first primary (e.g Sun)
 * \param second name of the second primary (e.g. Earth)
 * \param third name of the third primary  (e.g. Moon)
 *
 * At the end of this routine:
 *    - fbp->cr3bp1 is initialized as the CRTBP of (second, third).
 *    - fbp->cr3bp2 is initialized as the CRTBP of (first, second).
 *
 * WARNING: in the case of the Sun-Earth-Moon system (in fact, the only interesting case
 * for us...), we need to set (SUN, EARTH_AND_MOON), instead of (SUN, EARTH), in
 * fbp->cr3bp2. So keep in mind that if the configuration is SUN/EARTH/MOON,
 * the result will NOT be (first, second) in fbp->cr3bp2, but (first, second+third).
 *
 **/
void init_fbp(FBP* fbp, int first, int second, int third);

/**
 * \brief Initializes a certain Circular Restricted 3-Body Problem as a CR3BP structure.
 * \param cr3bp pointer on the CR3BP structure
 * \param name_1 name of the first primary
 * \param name_2 name of the second primary
 **/
void init_cr3bp(CR3BP* cr3bp, int name_1, int name_2);

/**
 * \brief Initializes a libration point.
 * \param libp a pointer towards the LibPoint structure to initialize.
 * \param cr3bp a CR3BP structure that contains useful coefficients.
 * \param number the indix of the libration point to init.
 **/
void init_lib_point(LibPoint* libp, CR3BP cr3bp, int number);

/**
* \brief Initialize one celestial body
* \param body a pointer on the Body structure to init.
* \param name the name of the body in integer format (consistent with SPICE numerotation)
**/
void init_body(Body* body, int name);

/**
 *  \brief Initializes two FBPL (model1 and model2) with two different models for
 *         continuation process from one model (epsilon = 0.0) to the other
 *        (epsilon = 1.0).
 **/
void init_fbp_cont(QBCP_I* model, FBPL* model1, FBPL* model2,
                int name_1, int name_2, int name_3, int is_norm, int li_EM, int li_SEM,
                int is_new, int mod1, int mod2, int coord_sys, int param_style);

//----------------------------------------------------------------------------------------
//            Change of coordinate systems
//----------------------------------------------------------------------------------------
/**
 *  \brief Change the default coordinate system of the FBPL structure to coord_sys.
 **/
void change_coord(FBPL& fbpl, int coord_sys);

/**
 *  \brief Change the default coordinate system to coord_sys and
 *         the libration point to li in the FBPL structure.
 **/
void change_li_coord(FBPL& fbpl, int coord_sys, int li);

//----------------------------------------------------------------------------------------
//            Subroutines - I/O
//----------------------------------------------------------------------------------------
/**
 *  \brief Return the string corresponding to the libration point number provided
 *         (e.g. "L1" if li == 1).
 **/
string init_F_LI(int li);

/**
 *  \brief Return the string corresponding to the model indix provided
 *         (e.g. "QBCP" if model == Csts::QBCP).
 **/
string init_F_MODEL(int model);

/**
 *  \brief Return the string corresponding to the coordinate system
 *         provided (e.g. "EM" if coord_sys == Csts::EM).
 **/
string init_F_COORDSYS(int coord_sys);

/**
 *  \brief Return the folder name corresponding to the prefix/model/framework/libration
 *         point number combination provided (e.g. "prefix/QBCP/EM/L1").
 **/
string init_F_FOLDER(string prefix, int model, int coord_sys, int li);

/**
 * \brief Retrieve a set of coefficients, given as Fourier series from a txt file.
 * \param filename the name of the txt file.
 * \param params a pointer toward the array to update.
 * \param n_order_fourier the order of the Fourier series.
 * \param shift the indix from which to start the storage of the coefficients in params.
 * \param flag: if flag == 1, the coefficients computed via Fast Fourier Transform (FFT)
 *        are used. Otherwise, the expansions obtained through Fourier series algebraic
 *        manipulations are used.
 **/
void read_fourier_coef(string filename, double* params, int n_order_fourier, int shift, int flag, int number);

//----------------------------------------------------------------------------------------
//            Subroutines - Computation
//----------------------------------------------------------------------------------------
/**
 * \brief Compute the CR3BP potential energy for the given state and mu.
 * \param y  the state array.
 * \param mu the mass ratio of the current CR3BP.
 **/
double crtbp_energy(double y[], double mu);

/**
 * \brief Compute the coefficient cn for a given libration point (L1, L2, and L3 for now)
 * \param fbpl a reference to the FBPL initialized around the selected libration point.
 * \param n the indix of the coefficient to compute.
 *  See double cn_coeff(int li, double gamma, double mu, int n).
 **/
double cn_coeff(FBPL& fbpl, int n);

/**
 * \brief Compute the coefficient cn for a given libration point (L1, L2, and L3 for now)
 * \param li the number of the current libration point (1,2,3)
 * \param gamma the gamma parameter associated to the current libration point
 * \param mu the mass ratio of the current TBP system
 * \param n the indix of the coefficient to compute.
 *
 * We recall that the cn coefficient are given by the following equation:
 * \f$  c_n = \frac{1}{\gamma_j^3} \left( (\pm1)^n \mu + (-1)^n \frac{(1-\mu) \gamma_j^{(n+1)}}{(1 -\mp \gamma_j)^{(n+1)}} \right) \f$ for \f$ L_j, j = 1,2 \f$, and where:
 * - The upper (resp. lower) signs stand for \f$ L_1 \f$ (resp. \f$ L_2 \f$).
 * - \f$ \gamma_j \f$ is the distance from \f$ L_j \f$ to the smallest primary.
 *
 * The L3 value is taken from Richardson, 1980:
 *
 * \f$  c_n = \frac{1}{\gamma_j^3} \left( 1 - \mu  + \frac{\mu \gamma_j^{n+1}}{(1 + \gamma_j)^{n+1}} \right) \f$.
 **/
double cn_coeff(int li, double gamma, double mu, int n);

/**
 * \brief Using the Newton-Raphson method, find the root of a function known to lie close
 *        to the first guess x1. The root will be refined until its accuracy is known
 *        within ± xacc. funcd is a user-supplied routine that returns both the
 *        function value and the first derivative of the function at the point x.
 **/
double rtnewt(void (*funcd)(double, int, double, double*, double*),
              double x1, double xacc, double mu, int number);

/**
 * \brief Provides the function value and its first derivative for Newton's method.
 *        f corresponds to the quintic equation satisfied by the Li-m2 distance for
 *        the L1/L2 cases and by 1-(Li-m1 distance) for the L3 case.
 **/
void quintic_eq(double mu, int number, double y, double* f, double* df);

/**
 *  \brief This routine performs a change of coordinates:
 *         From SYSTEM (EM or SEM) to NORMALIZED-CENTERED (NC) coordinates.
 *         It is just used here for the primaries.
 */
void sys_to_nc_prim(double Zc[3], double zc[3], double c1, double gamma);


#endif // DEFINE_ENV_H_INCLUDED
