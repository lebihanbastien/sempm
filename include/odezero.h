#ifndef CUSTOM_ODEZERO_H_INCLUDED
#define CUSTOM_ODEZERO_H_INCLUDED

/**
 * \file odezero.h
 * \brief Additional structures for ODE integrationd with event detection, based on GSL
 *        implementation.
 * \author BLB
 */

//GSL
#include <gsl/gsl_odeiv2.h>
#include <gsl_interp.h>
#include <gsl/gsl_roots.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <complex.h>
#include <gsl_complex_math.h>
#include <gsl/gsl_blas.h>

//Custom
#include "env.h"
#include "ode.h"
#include <sys/types.h>
#include <sys/wait.h>
#include "ofts.h"

//Gnuplot
extern "C" {
#include "gnuplot_i.h"
}


//----------------------------------------------------------------------------------------
//     Structures
//----------------------------------------------------------------------------------------
/**
 *  \struct value_params odezero.h
 *  \brief  Contains the parameters to trigger an event of the form
 *          f = (dim - value) = 0, with dim = x, y or z.
 **/
typedef struct value_params value_params;
struct value_params
{
    //Maximum of events allowed during integration
    int max_events;
    //Dimension along which the event is triggered: either 0 (x), 1 (y), or 2 (z).
    int dim;
    //The event is of the form f = (dim - value) = 0
    double value;
    //Direction of events: if 1, only crosses from negative to positive values of f are
    //kept. if -1, only crosses from positive to negative values of f are kept.
    //if 0, all values are kept.
    int direction;
};


/**
 *  \struct value_output odezero.h
 *  \brief  Contains the outputs of a value_function for an event procedure f(t, z, params) = 0.
 *
 *  When the function value of a value_function structure is called, a value_output structure is returned, containing:
 *   - The current value f(t, z, params)
 *   - The maximum number of events.
 *   - The direction of events.
 *   Note that both the maximum number of events and the direction of events are already
 *   present in the structure value_params of the \c value_function. However, they are
 *   also returned in the output in order to reproduce the MATLAB procedure for event detection.
 **/
typedef struct value_output value_output;
struct value_output
{
    //The current value f(t, z, params)
    double val;
    //The maximum number of events.
    int max_events;
    //Direction of events: if 1, only crosses from negative to positive values of f are
    //kept. if -1, only crosses from positive to negative values of f are kept.
    //if 0, all values are kept.
    int direction;
};

/**
 *  \struct value_function odezero.h
 *  \brief  Contains both the parameters and the function for an event procedure f(t, z, params) = 0.
 *
 *  When the function value of a value_function structure is called, a value_output structure is returned, containing:
 *   - The current value f(t, z, params)
 *   - The maximum number of events.
 *   - The direction of events.
 *   Note that both the maximum number of events and the direction of events are
 *   already present in the structure \c value_params of the value_function.
 *   However, they are also returned in the output in order to reproduce the MATLAB
 *    procedure for event detection.
 *   Note also that, given the structure value_params, only events of the form
 *   f(t, z, params)  = dim - value = 0, with dim = x, y, or z are possible.
 **/
typedef struct value_function value_function;
struct value_function
{
    //Pointer to the parameters for this event
    value_params *val_par;
    //Pointer to the function f(t, z, params)
    value_output (*value)(double, double [], void *);
};

/**
 *  \struct OdeParams odezero.h
 *  \brief  Structure dedicated to ode integration with an event detection via a
 *          value_function structure.
 **/
struct OdeParams
{
    //Initial time
    double t0;
    //Initial state
    double* y0;
    //Driver for ode integration
    gsl_odeiv2_driver *d;
    //Structure for event detection
    value_function *fvalue;
};


//----------------------------------------------------------------------------------------
//      Routines
//----------------------------------------------------------------------------------------
/**
 * \brief return a structure value_output inspired by MATLAB odezero implementation,
 *        containing:
 *         - val = y[desired dimension] - desired value
 *         - direction: the direction (increasing, decreasing or both)
 *         - max_events: the maximum events allowed by the user before the end of the
 *           integration
 **/
value_output linear_intersection(double t, double yv[], void *params);

/**
 * \brief return the value y[1] = y for a given integration t from t0 = 0.0
 *        and for given initial conditions provided in the structure pointer *params.
 **/
double odezero_event (double t, void *params);

/**
 * \brief return the value y[1] = y for a given integration t from t0 = 0.0
 *        and for given initial conditions provided in the structure pointer *params.
 **/
double ode_y_zero (double t, void *params);




#endif // CUSTOM_ODEZERO_H_INCLUDED

