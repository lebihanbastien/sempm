#ifndef LPDYNEQ_H_INCLUDED
#define LPDYNEQ_H_INCLUDED

#include "init.h"
#include "em_se_in.h"
extern "C" {
#include "gnuplot_i.h"
}

//========================================================================================
//         Dynamical equivalents to the libration points
//========================================================================================
/**
 *  \brief Computation of the dynamical equivalents of the libration points. Test function.
 *
 *         - This function makes use of the subroutine lpdyneq_single_shooting.
 *           The latter contains the initialization of the first guess
 *           (the geometrical positions of the CRTBP libration points)
 *           and the differential corrector (single shooting).
 *         - After lpdyneq_single_shooting, the resulting initial conditions are
 *           integrated and plotted on a full orbit.
 *         - Then, a test of the periodicity of the orbit is performed, via the
 *           computation of the error |zv(0) - zv(T)|.
 *         - Finally, the dynamical equivalent is computed again via a multiple shooting
 *           approach, for consistency.
 *         - If is_stored is true, the results (x(t), y(t)) in synodical coordinates are
 *           stored in txt files of the form: "./plot/QBCP/DYNEQ_QBCP_EM_L1_NC.txt".
 **/
void compute_dyn_eq_lib_point(FBPL &fbpl, int is_stored);

//========================================================================================
// Subroutines
//========================================================================================
/**
 *  \brief Main routine for the computation of the dynamical equivalent to the libration
 *         points via a single shooting procedure.
 *         The results are given in the form of the initial conditions (42 dimensions)
 *         updated in the array z0v[].
 **/
void lpdyneq_single_shooting(gsl_odeiv2_driver *d, double z0v[]);

/**
 *   \brief Computes M+1 patch points, indexed between 0 and M,
 *          along the T-periodic orbit starting at z0v = [x0, y0, z0, vx0, vy0, vz0, eye(6)]
 *          Each STM is initialized with the identity matrix.
 **/
void lpdyneq_patch_points(const double *z0v, gsl_odeiv2_driver *d, double t_period,
                          double **zm, double *tm, int M);


//========================================================================================
// Periodicity Condition
//========================================================================================
/**
 *  \brief Test the periodicity condition zv(0) = zv(t_period), with initial condition z0v,
 *         for the vector field contained in the driver d.
 *
 *  Note that the integer N is the number of variables associated to the driver d,
 *  and should be also the number of variables in zv.
 *  However, the periodicity condition is tested only on n_var_test variables (a usual
 *  example is n_var = 42 but n_var_test = 6).
 **/
int periodicity_condition(const double zv0[], int n_var, int n_var_test, double t_period,
                          gsl_odeiv2_driver *d, int is_norm);


#endif // LPDYNEQ_H_INCLUDED
