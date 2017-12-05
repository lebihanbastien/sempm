#ifndef PMODE_H_INCLUDED
#define PMODE_H_INCLUDED

/**
 * \file pmode.h
 * \brief Integration and test of the outputs of the parameterization method.
 * \author BLB
 */

#include "pmt.h"
#include "pmeval.h"

/**
 *  \brief Structure for the reduced vector field
 */
typedef struct RVF RVF;
struct RVF
{
    vector<Oftsc>* fh;   //Reduced vector field
    Ofsc* ofs;           //Auxiliary Ofs object
    int order;           //Order of the evaluation     (<= OFTS_ORDER)
    int ofs_order;       //Order of the ofs evaluation (<= OFS_ORDER)
    double n;            //Pulsation of the QBCP
};


//----------------------------------------------------------------------------------------
//
//          Integration of the pm
//
//----------------------------------------------------------------------------------------
/**
 *  \brief Computes the reduce vector field dot(s) = f(s,t)
 *  \param t the current time
 *  \param y the complex state s (CCM), with real and imag parts separated (CCM8 form)
 *  \param f the complex reduced vector field f(s,t), with real and imag parts separated (CCM8 form)
 *  \param params_void a (void) pointer towards a RVF structure that contains:
 *           (1) the Fourier-Taylor expansions of the rvf
 *           (2) the order of the evaluation
 *           (3) the pulsation n of the system
 *           (4) an OFS temp variable.
 *
 *     Note that the use of void *params_void instead of the direct use of an
 *     RVF structure is to comply with GSL ODE integration tools,
 *     that require a vector field of the exact form :
 *          int vf(double t, const double y[], double f[], void *params_void)
 *
 **/
int qbcp_fh(double t, const double y[], double f[], void* params_void);

//----------------------------------------------------------------------------------------
//
//          Tests
//
//----------------------------------------------------------------------------------------
/**
 *  \brief Computations of the orbital and invariance errors along an orbit initialized by
 *         the parameterization of the center manifold W(s,t) = PC(t)*Wh(s,t) + V(t), at
 *         order ofts_order for the Taylor series and ofs_order for the Fourier series.
 *         The errors are computed on the interval tvec and plotted.
 **/
int orb_inv_error(const double st0[],      //RCM initial conditions
                  vector<Oftsc>& Wh,       //TFC manifold
                  matrix<Ofsc>& PC,        //COC matrix: z = PC*zh+V
                  vector<Ofsc>& V,         //COC vector: z = PC*zh+V
                  vector<Oftsc>& FW,       //NC vector field
                  vector<Oftsc>& DWf,      //Jacobian of FW
                  vector<Oftsc>& Wdot,     //Partial derivative of W wrt time
                  double tvec[2],          //Integration interval
                  OdeStruct* ode_s_nc,     //ode structure for NC integration
                  OdeStruct* ode_s_rvf,    //ode structure for RVF integration (reduced vector field)
                  int Npoints,             //Number of points on which the errors are estimated
                  FBPL& fbpl,              //current QBCP
                  int ofts_order,          //Order for the eval of the OFTS objects
                  int ofs_order,           //Order for the eval of the OFS objects
                  gnuplot_ctrl**  ht,      //Gnuplot handlers
                  int color);              //Color of the plots

/**
 *  \brief Test of the parameterization of the central manifold on a
 *         given orbit, through the computation of various errors (eO, eI) and the
 *         relative energy.
 *              - The orbit is initialized by the reduced state vector si[].
 *              - The errors are computed on the interval tvec.
 *              - The errors are computed for the n_ofts_order Taylor orders stored in v_ofts_order.
 *              - The errors are computed for the n_ofs_order Fourier orders stored in v_ofs_order.
 *
 *   Note that: the expansions FW, DWf and Wdot are taken from files,
 *   whereas the parameterization itself CM and CMh, and the reduced vector field Fh
 *   are taken from global objects, defined in init.cpp.
 *
 *   Requires init_inv_man().
 *
 **/
void pm_error_vs_orders_test(int n_ofts_order, int v_ofts_order[],
                             int n_ofs_order, int v_ofs_order[],
                             double si[], double tvec[]);

/**
 *  \brief Computations of the filenames used in the routine orb_inv_error (orbital error, invariance error, etc).
 **/
string orb_inv_error_output_name(string f_plot, string title, const double *st0, int ofts_order, int ofs_order);

//----------------------------------------------------------------------------------------
//
//          DEPRECATED (to be cleaned up)
//
//----------------------------------------------------------------------------------------
/**
 *  \brief Evaluates the contributions of each order in W to the computation of W(s,t), with an arbitrary state (s,t)
 **/
void pmContributions();

/**
 *  \brief Evaluates the l1/linf norm of each order in W.
 **/
void pmNorms();

/**
 *  \brief Small divisors under a certain value
 **/
void pmSmallDivisors(double sdmax);

/**
 *  \brief Evaluates the variations of the initial conditions (IC) wrt to the order in W(s,t), with an arbitrary state (s,t)
 **/
void pmTestIC();

#endif // PMODE_H_INCLUDED
