#ifndef QBCP_H_INCLUDED
#define QBCP_H_INCLUDED

/**
 * \file qbcp.h
 * \brief Implementation of various vector fields in the Sun-Earth-Moon QBCP/CRTBPs,
 *        as well as additional routines for Hamiltonian & tests.
 *
 *        Note that most of the routines are prefixed by qbcp_ despite the fact that they
 *        handle both the QBCP and CRTBP cases (could be changed in future release).
 *
 * \author BLB.
 */

#include <iostream>
#include <time.h>
#include <stdio.h>
#include <complex.h>
#include <gsl_complex_math.h>
#include <fstream>
#include <sstream>

//Custom
#include "ode.h"
#include "qbtbp.h"
#include "diffcorr.h"
#include "env.h"
#include "odezero.h"
#include "init.h"
#include "em_se_in.h"

extern "C" {
#include "gnuplot_i.h"
}

//----------------------------------------------------------------------------------------
//
// Derivatives for ODE integration
//
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
// Vector fields
//----------------------------------------------------------------------------------------
/**
 *  \brief Vector field of the Sun-Earth-Moon QBCP/CRTBP in NC(EM or SE) coordinates
 *         (no variationnal equations).
 **/
int qbcp_vfn_novar(double t, const double z[], double f[], void *params_void);

/**
 *  \brief Vector field of the Sun-Earth-Moon QBCP/CRTBP in NC(EM or SE) coordinates,
 *         with full nonlinear variational equations.
 **/
int qbcp_vfn_varnonlin(double t, const double z[], double f[], void *params_void);

/**
 *  \brief Vector field of the Sun-Earth-Moon QBCP/CRTBP in NC(EM or SE) coordinates,
 *         with linear variational equations, and taking into account the time-dependent
 *         translation so that the dynamical equivalent to the libration point becomes
 *         a fixed point (see nfo2.cpp for practical use).
 **/
int qbcp_vfn_varlin_trans(double t, const double z[], double f[], void *params_void);

/**
 *  \brief Full nonlinear variational equations of the Sun-Earth-Moon QBCP/CRTBP
 *         in NC(EM or SE) coordinates.
 **/
int qbcp_Dfn_varnonlin(double t, const double z[], double **Df, void *params_void);

/**
 *  \brief Vector field of the Sun-Earth-Moon QBCP/CRTBP in EM or SE coordinates.
 **/
int qbcp_vf(double t, const double z[], double f[], void *params_void);

//----------------------------------------------------------------------------------------
// Inner routines for the computation of the vector fields
//----------------------------------------------------------------------------------------
/**
 *  \brief Build the Variational Equation Matrix Q from the arrays b and alpha.
 *         Note that alpha[14] (alpha15) is zero for the QBCP
 **/
void set_vareq_matrix(gsl_matrix *Q, double b[], double alpha[]);

/**
 *  \brief Update the vector field of the state in system (EM or SE) coordinates.
 *         Note that alpha[14] (alpha15) is zero for the QBCP
 **/
int vf_state( const double z[], double f[], double alpha[],
              double ps[], double pe[], double pm[],
              double qps2, double qpe2, double qpm2,
              double ms, double me, double mm );

/**
 *  \brief Update the vector field of the state in NC(EM or SE) coordinates.
 *         Note that alpha[14] (alpha15) is zero for the QBCP
 **/
int vfn_state(const double z[], double f[], double alpha[],
              double ps[], double pe[], double pm[],
              double qps2, double qpe2, double qpm2,
              double ms, double me, double mm,
              double gamma);


/**
 *  \brief Update the variational equation matrix in NC(EM or SE) coordinates.
 **/
int vfn_stm(const double z[], gsl_matrix *Q, double alpha[],
            double ps[], double pe[], double pm[],
            double qps2, double qpe2, double qpm2,
            double ms, double me, double mm,
            double gamma);

/**
 *  \brief Update the variational equation matrix in NC(EM or SE) coordinates.
 *         Case of a double **Df, instead of a gsl_matrix.
 **/
int vfn_stm(const double z[], double **Df, double alpha[],
            double ps[], double pe[], double pm[],
            double qps2, double qpe2, double qpm2,
            double ms, double me, double mm,
            double gamma);

/**
 *  \brief Update the variational equation matrix in system coordinates (non-normalized)
 **/
int vf_stm(const double z[], gsl_matrix *Q, double alpha[],
           double ps[], double pe[], double pm[],
           double qps2, double qpe2, double qpm2,
           double ms, double me, double mm);

/**
 *  \brief Update the linearized variational equation matrix in NC(EM or SE) coordinates.
 **/
int vfn_stm_lin_trans(const double z[], gsl_matrix *Q, double alpha[],
                     double ps[], double pe[], double pm[],
                     double ms, double me, double mm,
                     double gamma);


/**
 *  \brief Computes the second-order derivatives of the linearized potential for one
 *         given primary
 **/
double Uijlin(double pc[], double qpc2, double mc, double factor, int i, int j);

/**
 *  \brief Computes the second-order derivatives of the full nonlinear potential for one
 *         given primary
 **/
double Uij(const double z[], double pc[], double qpc2, double mc, double factor, int i, int j);

//----------------------------------------------------------------------------------------
//
// Hamiltonians
//
//----------------------------------------------------------------------------------------
/**
 *  \brief Hamiltonian of the Sun-Earth-Moon QBCP/CRTBP in EM or SE coordinates
 **/
double qbcp_H(double t, const double z[], void *params_void);

/**
 *  \brief Hamiltonian of the Sun-Earth-Moon QBCP/CRTBP in NC(EM or SE) coordinates
 **/
double qbcp_Hn(double t, const double z[], void *params_void);

//----------------------------------------------------------------------------------------
//
// Floquet analysis
//
//----------------------------------------------------------------------------------------
//
/**
 *  \brief Same as the routine qbcp_vfn_varlin_trans + Integration of the matrix Pbar
 *         that appears in the Floquet analysis of the hamiltonian of order 2 in the
 *         neighborhood of L1,2 (see nfo2.h, and BLB thesis manuscript).
 **/
int qbcp_vfn_varlin_trans_P(double t, const double z[], double f[], void *params_void);

/**
 *  \brief QBCP vector field in the form of the matrix Q1(t), Q2(t) and Q3(t) defined by:
 *              x' =  Q2^T*x + 2Q3*Y
 *              z' = -2Q1 *x -  Q2*z
 *          Used in the Floquet analysis (see nfo2.h)
 *
 *        The vector field is in NC(EM or SE) coordinates.
 **/
int qbcp_vfn_Q(double t, const double z[], gsl_matrix *Q1, gsl_matrix *Q2, gsl_matrix *Q3,
void *params_void);

//----------------------------------------------------------------------------------------
// COC: tests & plots
//----------------------------------------------------------------------------------------
/**
 *  \brief Plot the Primaries' Three-Body motion of the QBCP.
 **/
void QBTBP_IN();


#endif // QBCP_H_INCLUDED
