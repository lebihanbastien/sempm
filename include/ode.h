#ifndef CUSTOM_ODE_H_INCLUDED
#define CUSTOM_ODE_H_INCLUDED

/**
 * \file ode.h
 * \brief Additional structures for ODE integrationd, based on GSL implementation.
 * \author BLB
 */

//C++
#include <iostream>

//GSL
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

//Custom
#include "config.h"

using namespace std;


//----------------------------------------------------------------------------------------
// Structure
//----------------------------------------------------------------------------------------
/**
 *   \struct OdeStruct. A super structure for ODE integration based on GSL objects.
 *           Contains a stepper, a controller, a driver, a root solver, and the
 *           associated precisions.
 **/
typedef struct OdeStruct OdeStruct;
struct OdeStruct
{
    //Stepper
    const gsl_odeiv2_step_type *T;
    gsl_odeiv2_step *s;
    //Control
    gsl_odeiv2_control *c;
    //Evolution
    gsl_odeiv2_evolve *e;
    //System
    gsl_odeiv2_system sys;
    //Driver
    gsl_odeiv2_driver * d;
    //Solver for root finding
    gsl_root_fsolver *s_root;

    //Precisions
    double eps_int_rel;  //for integration (relative precision)
    double eps_int_abs;  //for integration (absolute precision)
    double eps_root;     //for root finding

    //Initial step
    double h;

    //Dimension
    int dim;
};

//----------------------------------------------------------------------------------------
// Initialization
//----------------------------------------------------------------------------------------
/**
 *  \brief Initializes an ode structure of the type OdeStruct, including:
 *           - a step type T,
 *           - a root solver type T_root,
 *           - a root solver s_root,
 *           - a root precision eps_root
 *           - an absolute integration precision eps_int_abs
 *           - a  relative integration precision eps_int_rel
 *           - a diffcoor precision esp_diff
 *           - a number of dimension dim,
 *           - a default step h,
 *           - a vector field func,
 *           - its corresponding jacobian,
 *           - some additional parameters in params.
 *
 *        During this initialisation, the following elements are initialized:
 *           - a stepper s,
 *           - a controler c,
 *           - a system sys,
 *           - a driver d containing all the previous objects.
 **/
void init_ode_structure(OdeStruct *ode_s,
                        const gsl_odeiv2_step_type *T,
                        const gsl_root_fsolver_type *T_root,
                        double eps_int_abs,
                        double eps_int_rel,
                        double eps_root,
                        size_t dim,
                        double h,
                        int (* func) (double t, const double y[], double dydt[], void *params),
                        int (* jacobian) (double t, const double y[], double * dfdy, double dfdt[], void * params),
                        void *params);

/**
 *  \brief Initializes an ode structure of the type OdeStruct, including:
 *           - a step type T,
 *           - a root solver type T_root,
 *           - a root solver s_root,
 *           - a number of dimension dim,
 *           - a vector field func,
 *           - its corresponding jacobian,
 *           - some additional parameters in params.
 *
 *        During this initialisation, the following elements are initialized:
 *           - a stepper s,
 *           - a controler c,
 *           - a system sys,
 *           - a driver d containing all the previous objects.
 *
 *        The default precision of the propagator are taken in Config::configManager().
 *        Moreover, the jacobian is always NULL.
 **/
void init_ode_structure(OdeStruct *ode_s,
                        const gsl_odeiv2_step_type *T,
                        const gsl_root_fsolver_type *T_root,
                        size_t dim,
                        int (* func) (double t, const double y[], double dydt[], void *params),
                        void *ode_params);

//----------------------------------------------------------------------------------------
// Free
//----------------------------------------------------------------------------------------
/**
 *  \brief Free the space allocated to an ode structure.
 *         Namely, the stepper s, the controller c, and the driver d.
 **/
void free_ode_structure(OdeStruct *ode_s);

//----------------------------------------------------------------------------------------
// Update, reset
//----------------------------------------------------------------------------------------
/**
 *  \brief Reset an ode structure.
 *         Namely, the stepper s, the controller c, and the driver d.
 **/
void reset_ode_structure(OdeStruct *ode_s);

/**
 *  \brief Update an ode structure of the type OdeStruct, including:
 *           - a step type T,
 *           - a stepper s,
 *           - a controler c,
 *           - a system sys,
 *           - a driver d,
 *           - a root solver s_root,
 *           - a root precision eps_root
 *           - a diffcoor precision esp_diff
 **/
void update_ode_structure(OdeStruct *ode_s,
                         const gsl_odeiv2_step_type *T,
                         gsl_odeiv2_step *s,
                         gsl_odeiv2_control *c,
                         gsl_odeiv2_evolve *e,
                         gsl_odeiv2_system sys,
                         gsl_odeiv2_driver * d,
                         gsl_root_fsolver *s_root,
                         double eps_root);

/**
 * \brief Set the right direction of integration in OdeStruct ode_s
 **/
void set_dir_ode_structure(OdeStruct *ode_s, double t0, double t1);


/**
 * \brief Set the right direction of integration in a driver d
 **/
void set_dir_driver(gsl_odeiv2_driver *d, double t0, double t1);

#endif // CUSTOM_ODE_H_INCLUDED
