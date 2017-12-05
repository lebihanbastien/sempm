#ifndef POINCARE_H_INCLUDED
#define POINCARE_H_INCLUDED

extern "C"
{
    #include "nrutil.h"
}


#include "pmode.h"
#include "init.h"

//Integration methods
#define DUAL_INT          1
#define DUAL_INT_NO_RESET 2
#define DUAL_INT_STEPPED  3
#define SINGLE_INT        4

//Types of map
#define PMAP        1
#define TMAP        2
#define EMAP        3
#define IMAP        4
#define HMAP        5
#define IMAPPLANAR  6


/**
 * \file poincare.h
 * \brief Computation of Poincare, Error, and Period map for the QBCP.
 * \author BLB.
 * \date May 2015
 * \version 1.0
 *
 */

/**
 *  \struct Pmap
 *  \brief  Poincare map parameters.
 *
 *  Contains all necessary parameters and constants to define a suitable Poincare map
 *  (energy level, order of expansions, maximum number of events, initial time...).
 *
 **/
typedef struct Pmap Pmap;
struct Pmap
{
    //Expansions
    int ofts_order;      //the order of the expansions
    int ofs_order;       //the order of the expansions

    //Integration & event
    double tmax_nc;      //final time
    double t0_nc;        //init time
    double t_proj_nc;    //test time
    double f_proj_T;     //test projection frequency
    double T;            //period
    int max_events;      //maximum z=0 events allowed
    double eproj;        //the threshold above which a reset of z(t) is needed
    double pabs;         //absolute tolerance during integration
    double prel;         //relative tolerance during integration
    double proot;        //the precision on the roots z = 0

    //Energy
    double h_nc;           //value of the Hamiltonian
    double h_nc_li;        //value of the Hamiltonian @Li
    double dh_nc_t0;          //h_nc-h_nc_li at t = t0

    //Grid
    int  n_sol;          //number of solutions
    double si_min;       //min of the grid
    double si_max;       //max of the grid

    //Checking divergence
    double max_rad_div;

    //Dimension
    int var_dim_h;

    //Model bool
    bool graph_style;

    //Type
    int type;

    //Plot
    int is_plot;

    //Method of integration
    int int_method;

    //Parallel computation of not
    int is_par;
};

/**
 *  \struct Orbit
 *  \brief  Defines a given orbit with proper arrays to store results.
 **/
typedef struct Orbit Orbit;
struct Orbit
{
    //------------------------------------------------------------------------------------
    //Parent
    //------------------------------------------------------------------------------------
    Pmap   *pmap;       //Poincare map (parent)
    FBPL   *fbpl;       //QBCP around a given Li point (parent)

    //------------------------------------------------------------------------------------
    //Parameterization (common to all orbits)
    //------------------------------------------------------------------------------------
    vector<Oftsc>*  W;      //z(t) = W(s(t), t)
    vector<Oftsc>*  Wh;     //zh(t) = Wh(s(t), t)
    matrix<Oftsc>*  DW;     //Jacobian of W
    Ofsc* ofs;              //Auxiliary Ofs object
    double  n;              //Pulsation of the QBCP


    //------------------------------------------------------------------------------------
    //COC (common to all orbits)
    //------------------------------------------------------------------------------------
    matrix<Ofsc>* PC;       //COC matrix
    vector<Ofsc>* V;        //COC vector

    //------------------------------------------------------------------------------------
    //For event detection
    //------------------------------------------------------------------------------------
    value_params *val_par;  //Event parameters
    value_function *fvalue; //Value function for event detection
    double** z0_mat;        //Pointer towards the stored position of events (NC)
    double** zh0_mat;       //Pointer towards the stored position of events (TFR)
    double** s0_mat;        //Pointer towards the stored position of events (RCM)
    double*  te_mat;        //Pointer towards the stored time of events
    double*  hz;            //Energy z(t) at each event
    double*  hw;            //Energy W(s(t),t) at each event
    double*  ePm;           //Projection error at each event
    int last_indix;         //Index of the last computed event in z0_mat
    int reset_number;       //number of reset during computation (for dual integration)
    int      *nevent;       //the label of the events

    //------------------------------------------------------------------------------------
    //Characteristics
    //------------------------------------------------------------------------------------
    double   *z0;           //Initial position in NC coordinates dim = 6
    double   *zh0;          //Initial position in TFR coordinates dim = 6
    double   *si;           //Initial RCM configuration dim = REDUCED_NV
    double   *s0d;          //Initial position in CCM8 coordinates (real+imag part) dim = 2*REDUCED_NV
    cdouble  *s0;           //Initial position in CCM8 coordinates (real+imag part) dim = 4
    double   *xf;           //Final position dim = 6
    double    tf;           //final time after computation
    double    eOm;          //Orbit error
    int       int_method;   //integration method used to compute the orbit; -1 if not computed
    int       label;        //label of the orbit
    int       var_dim_h;    //the dimension along wich we guarantee a certain initial energy value


    //------------------------------------------------------------------------------------
    //ODE integration
    //------------------------------------------------------------------------------------
    OdeStruct *ode_s_6;     //NC ode structure
    OdeStruct *ode_s_6_root; //NC ode structure
};


//----------------------------------------------------------------------------------------
//
//  Poincare map
//
//----------------------------------------------------------------------------------------
/**
 *   \brief Computes a Poincare map
 *   \param pmap a reference to the Poincare map parameters
 *
 *    Requires init_inv_man and init_coc
 **/
void pmap_build_random(Pmap& pmap);


//----------------------------------------------------------------------------------------
//
//  Init ode routines
//
//----------------------------------------------------------------------------------------

/**
 *   \brief Initialize the ode structure for the NC vector field
 **/
void pmap_init_ode(OdeStruct& ode_s_6, Pmap& pmap);


//----------------------------------------------------------------------------------------
//
//  Print & Read
//
//----------------------------------------------------------------------------------------
/**
 *  \brief Returns the filename associated with the current Pmap and method of computation.
 **/
string pmap_filename(Pmap& pmap);

/**
 *   \brief Reset a given binary file
 **/
void header_fprint_bin(string filename);

/**
 *   \brief Print the poincare map of and orbit in a txt file. Binary version. Only the points for which pz > 0 are kept.
 **/
void pmap_orbit_fprint_bin(Orbit* orbit, string filename, int append);

//----------------------------------------------------------------------------------------
//
//  Init of one orbit
//
//----------------------------------------------------------------------------------------
/**
 *   \brief Update the array st0 with values from orbit.s0 and the value sr. Used in orbit_init_pmap and delta_h_nc_orbit.
 *
 *          The following patterns are followed (examples):
 *              - if SEML.model = CRTBP and orbit.dim = 2: s2 = sr, s4 = 0
 *                   st0[0] = orbit->si[0];;
 *                   st0[1] = sr;
 *                   st0[2] = orbit->si[2];
 *                   st0[3] = 0.0;
 *
 *             - if SEML.model = QBCP and orbit.dim = 2: s2 = s4 = sr
 *                   st0[0] = sr;
 *                   st0[1] = orbit->si[1];
 *                   st0[2] = sr;
 *                   st0[3] = orbit->si[3];
 **/
void orbit_update_s0(Orbit* orbit, double st0[], double sr);

/**
    \brief Initialize an orbit wrt a Poincare map so that H(orbit.s0) = H(Pmap)
 **/
int orbit_init_pmap(Orbit& orbit, double st0[]);

/**
    \brief Update the initial conditions (si, s0, z0 and s0d) of the orbit given an array of initial TFC conditions si
 **/
void orbit_update_ic(Orbit& orbit, const double si[], double t0);

/**
 *   \brief Computes the hamiltonian at the position st0, in system coordinates and units.
 *   \param pmap a reference to the pmap that carries a set of useful parameters
 *   \param st0 the input state
 *
 *   WARNING: Direct use of CM inside this
 **/
double h_nc_pmap(Pmap& pmap, double st0[]);

/**
    \brief Computes the difference between an given Ham value and the state configuration defined in the routine orbit_update_s0 (see comments therein).
    \param  sr a double to complete the current tested configuration
    \param  params a pointer to the orbit with a given Ham value (why void? the idea is to a have a generic function, but might be useless at this point)
    \return the difference between the two hamiltonians
 **/
double delta_h_nc_orbit(double sr, void* params);

//----------------------------------------------------------------------------------------
//
// Computation of one orbit
//
//----------------------------------------------------------------------------------------
/**
 *   \brief Computes the poincare map of one given orbit, with the following method:
 *
 *          - integration of the NC vector field and reset of the state with the
 *            semi-analytical approximation of the center manifold each time the
 *            xy-plane is crossed.
 *               -
 *   \param orbit a pointer to the orbit
 **/
int orbit_compute_pmap(Orbit* orbit, int method);

/**
 *   \brief Computes the poincare map of one given orbit, with a single method: only the Nc vector field is computed. The distance with respect to the
 *          central manifold is evaluated each time that z = 0.
 *   \param orbit a pointer to the orbit
 **/
int single_pmap_proj(Orbit* orbit);

/**
 *   \brief Computes the poincare map of one given orbit, with a single method: only the Nc vector field is computed. The distance with respect to the
 *          central manifold is evaluated each time that z = 0.
 *   \param orbit a pointer to the orbit
 **/
int single_pmap_proj(Orbit* orbit);


/**
 *   \brief Refine the root z = 0 for a single integration NC+rvf. Used in single_pmap_proj.
 **/
int refine_root(Orbit* orbit,
                double* yv,
                double* t,
                double* s1,
                double* yv_mat,
                double* t_mat,
                double omega1,
                double omega3,
                int events);

//----------------------------------------------------------------------------------------
//
// Steppers with projection
//
//----------------------------------------------------------------------------------------
/**
 *   \brief Integrates one step the current state yv using projection on the CM if necessary
 **/
int gslc_proj_step(Orbit* orbit,
                   double yv[],
                   double* t,
                   double t0,
                   double t1,
                   double* ePm,
                   int* nreset,
                   double omega1,
                   double omega3,
                   int is_reset_on);

/**
 *   \brief Integrates the current state yv up to t = t1, using projection on the CM if necessary
 **/
int gslc_proj_evolve(Orbit* orbit,
                     double yv[],
                     double* t,
                     double t0,
                     double t1,
                     double* ePm,
                     int* nreset,
                     double omega1,
                     double omega3,
                     int is_reset_on);



//----------------------------------------------------------------------------------------
//
// Plotting
//
//----------------------------------------------------------------------------------------
void orbit_poincare_plot(Orbit* orbit, gnuplot_ctrl* h1, int color);

//----------------------------------------------------------------------------------------
//
//  Orbit C structure handling
//
//----------------------------------------------------------------------------------------
/**
    \brief Initialize one orbit structure
 **/
void init_orbit(Orbit* orbit, vector<Oftsc>*  W, vector<Oftsc>*  Wh,
                matrix<Oftsc>* DW, matrix<Ofsc>*  PC, vector<Ofsc>*  V,
                value_params*   val_par, value_function* fvalue,
                OdeStruct* ode_s_6, OdeStruct* ode_s_6_root,
                Pmap* pmap, FBPL* fbpl, Ofsc* orbit_ofs, int var_dim_h, int label);


/**
    \brief Free one orbit
 **/
void free_orbit(Orbit* orbit);




#endif // POINCARE_H_INCLUDED
