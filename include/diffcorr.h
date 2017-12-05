#ifndef DIFFCORR_H_INCLUDED
#define DIFFCORR_H_INCLUDED

/**
 * \file  diffcorr.h
 * \brief Contains all the routines that perform differential correction procedures.
 * \author BLB.
 * \date   2016
 * \version 1.0
 */

#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_odeiv2.h>

#include "qbtbp.h"
#include "vf.h"
#include "odezero.h"

extern "C"{
#include "gnuplot_i.h"
#include "nrutil.h"
}

//----------------------------------------------------------------------------------------
//  Inner routines
//----------------------------------------------------------------------------------------
/**
 *  \brief Solve a definite-positive linear system:
 *         - Decompose a ncs x ncs matrix M that is definite-positive, using GSL routines.
 *         - Then inverse the system err_vec = M*corr_vec.
 **/
int ftc_inv_dfls(gsl_matrix* M, gsl_vector* err_vec, gsl_vector* corr_vec, int ncs);

/**
 *  \brief Computes the correction vector associated to the minimum norm solution.
 *         Given:
 *              - an ncs x 1   error vector err_vec
 *              - an nfv x ncs Jacobian err_jac,
 *         This routine computes the correction vector associated to
 *         the minimum norm solution:
 *
 *              corr_vec = err_jac^T x (err_jac x err_jac^T)^{-1} err_vec.
 **/
int ftc_corrvec_mn(gsl_vector* mn_vec, gsl_vector *err_vec, gsl_matrix* err_jac, int nfv, int ncs);

/**
 *  \brief Computes the correction vector for a square system.
 *         Given:
 *              - an ncs x 1   error vector err_vec
 *              - an ncs x ncs Jacobian err_jac,
 *         This routine computes the correction vector given by:
 *
 *              corr_vec = err_jac^{-1} x err_vec.
 **/
int ftc_corrvec_square(gsl_vector* corr_vec, gsl_vector *err_vec, gsl_matrix* err_jac, int ncs);


//----------------------------------------------------------------------------------------
//  Single shooting
//----------------------------------------------------------------------------------------
/**
 *  \brief Performs a differential correction procedure on zv0 in order to get a
 *         periodic orbit of period t_period.
 *         The algorithm assumes that the orbit is in the xy plane, is symmetric wrt to
 *         the x-axis, and has a period T = t_period.
 *
 *         The number of variable is 42: 6 for the state, 36 for the
 *         State Transition Matrix (STM).
 *
 *  \param zv0 the initial conditions that needs to be corrected.
 *  \param t_period the final integration time (period of the orbit).
 *  \param eps_diff: the desired precision (magnitude of the final error).
 *  \param d:  a gsl_odeiv2_driver object used to integrate the equations of motion.
 *  \param n_var: the dimension of zv0.
 *  \param is_plot: if true, the steps of the diffcorr procedure are plotted in a
 *         temporary gnuplot window.
 *
 *   In brief: the algorithm corrects  [ zv0[0] zv0[4] ] in order to get
 *   [y[1] y[3]](t1) = 0, i.e. when the trajectory crosses the line t = t1.
 **/
int single_shoot_sym(double zv0[], double t_period, double eps_diff,
                     gsl_odeiv2_driver *d, int is_plot);


//----------------------------------------------------------------------------------------
//  Multiple shooting
//----------------------------------------------------------------------------------------
/**
 * \brief Multiple shooting scheme with periodicity conditions, on n_patch+1 patch points.
 *
 *   \param (zmd[0:41][0:n_patch], tmd[0:n_patch]) is the first guess (input).
 *   \param (zmdn[0:41][0:n_patch], tmdn[0:n_patch]) is the corrected output.
 *   \param d:  a gsl_odeiv2_driver object used to integrate the equations of motion.
 *   \param n_patch: the number of patch points
 **/
int mult_shoot_period(double** zmd, double* tmd, double** zmdn, double* tmdn,
                      gsl_odeiv2_driver* d, int n_patch, double eps_diff,
                      int is_plot, gnuplot_ctrl* h1, int strong_conv);



//----------------------------------------------------------------------------------------
//  Plot
//----------------------------------------------------------------------------------------
/**
 *  \brief Integrate the state zv[] up to t = t1 on a n_points grid, in either NC or SYS
 *         coordinates, and plot the corresponding result in SYS coordinates
 *         (e.g. EM or SEM coord.).
 *         Then the results are plotted in the xy-plane on a temporary gnuplot window
 *         via the handle *h1.
 *         Print in txt files is included via the integer is_stored.
 **/
int ode_plot_xy(const double zv[], int n_var, double t1, gsl_odeiv2_driver *d,
                gnuplot_ctrl  *h1, int n_points, int color, int is_norm, int is_stored,
                string legend, string filename);

/**
 *  \brief Integrate the states (zmdn**, tmdn*), in either NC or SYS coordinates,
 *         and plot the corresponding result in SYS coordinates (e.g. EM or SEM coord.).
 *         Then the results are plotted in the xy-plane  on a temporary gnuplot window
 *         via the handle *h1.
 **/
int ode_plot_vec(double **zmdn, double *tmdn, int n_var, int n_patch, gsl_odeiv2_driver *d,
               gnuplot_ctrl  *h1, int n_points, int color, int is_norm, string legend);

#endif // DIFFCORR_H_INCLUDED
