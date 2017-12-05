#include "odezero.h"

/**
 * \file odezero.cpp
 * \brief Additional structures for ODE integration with event (zero value) detection,
 *        based on GSL implementation.
 * \author BLB
 */

/**
 * \brief return a structure value_output inspired by MATLAB odezero implementation,
 *        containing:
 *         - val = y[desired dimension] - desired value
 *         - direction: the direction (increasing, decreasing or both)
 *         - max_events: the maximum events allowed by the user before the end of the
 *           integration
 **/
value_output linear_intersection(double t, double yv[], void *params)
{
    value_params *p = (value_params *) params;
    value_output val_par;
    val_par.val        = yv[p->dim] - p->value;
    val_par.direction  = p->direction;
    val_par.max_events = p->max_events;

    return val_par;
}

/**
 * \brief return the value y[1] = y for a given integration t from t0 = 0.0
 *        and for given initial conditions provided in the structure pointer *params.
 **/
double odezero_event (double t, void *params)
{
    //------------------------------------------------------------------------------------
    //Init
    //------------------------------------------------------------------------------------
    int i;
    struct OdeParams *p = (struct OdeParams *) params;
    double ystart[42];
    for(i=0; i<42; i++) ystart[i] = (p->y0)[i];
    double t0 = p->t0;

    //------------------------------------------------------------------------------------
    //Integration
    //------------------------------------------------------------------------------------
    //Starting in the right direction
    p->d->h = (t>t0) ? fabs(p->d->h) : -fabs(p->d->h);
    gsl_odeiv2_driver_apply (p->d, &t0, t, ystart);

    //------------------------------------------------------------------------------------
    //Reset
    //------------------------------------------------------------------------------------
    gsl_odeiv2_driver_reset(p->d);

    return p->fvalue->value(t0, ystart, p->fvalue->val_par).val;
}


/**
 * \brief return the value y[1] = y for a given integration t from t0 = 0.0
 *        and for given initial conditions provided in the structure pointer *params.
 **/
double ode_y_zero (double t, void *params)
{
    //------------------------------------------------------------------------------------
    //Init
    //------------------------------------------------------------------------------------
    int i;
    struct OdeParams *p = (struct OdeParams *) params;
    double ystart[42];
    for(i=0; i<42; i++) ystart[i] = (p->y0)[i];
    double t0 = p->t0;

    //------------------------------------------------------------------------------------
    //Integration
    //------------------------------------------------------------------------------------
    gsl_odeiv2_driver_apply (p->d, &t0, t, ystart);

    //------------------------------------------------------------------------------------
    //Reset
    //------------------------------------------------------------------------------------
    gsl_odeiv2_driver_reset(p->d);
    return ystart[1];
}
