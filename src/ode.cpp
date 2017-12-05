#include "ode.h"

/**
 * \file ode.cpp
 * \brief Additional structures for ODE integrationd, based on GSL implementation.
 * \author BLB
 */


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
                        void *ode_params)
{
    //Precisions and initial step
    ode_s->eps_root = eps_root;
    ode_s->h = h;
    ode_s->eps_int_abs = eps_int_abs;
    ode_s->eps_int_rel = eps_int_rel;
    ode_s->dim = dim;

    //Stepper
    ode_s->T = T;
    ode_s->s = gsl_odeiv2_step_alloc(T, dim);
    if (ode_s->s == NULL)
    {
        cout << "failed to allocate stepper object in init_ode_structure" << endl;
    }

    //Control
    ode_s->c = gsl_odeiv2_control_y_new(eps_int_abs, eps_int_rel);
    if (ode_s->c == NULL)
    {
        cout << "failed to allocate control object in init_ode_structure"  << endl;
    }

    //Evolution
    ode_s->e = gsl_odeiv2_evolve_alloc(dim);
    if (ode_s->e == NULL)
    {
        cout << "failed to allocate evolution object in init_ode_structure"  << endl;
    }


    //System
    ode_s->sys.function = func;
    ode_s->sys.jacobian = jacobian;
    ode_s->sys.dimension = dim;
    ode_s->sys.params = ode_params;

    ode_s->d = gsl_odeiv2_driver_alloc_y_new (&ode_s->sys, ode_s->T, ode_s->h, ode_s->eps_int_abs, ode_s->eps_int_rel);
    if (ode_s->d == NULL)
    {
        cout << "failed to allocate driver object in init_ode_structure"  << endl;
    }

    ode_s->s_root = gsl_root_fsolver_alloc (T_root);
    if (ode_s->s_root == NULL)
    {
        cout << "failed to allocate root_fsolver object in init_ode_structure"  << endl;
    }

}

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
                        void *ode_params)
{
    init_ode_structure(ode_s, T, T_root,
                       Config::configManager().G_PREC_ABS(),
                       Config::configManager().G_PREC_REL(),
                       Config::configManager().G_PREC_ROOT(),
                       dim,
                       Config::configManager().G_PREC_HSTART(),
                       func, NULL, ode_params);
}


//----------------------------------------------------------------------------------------
// Free
//----------------------------------------------------------------------------------------
/**
 *  \brief Free the space allocated to an ode structure.
 *         Namely, the stepper s, the controller c, and the driver d.
 **/
void free_ode_structure(OdeStruct *ode_s)
{
    gsl_odeiv2_step_free(ode_s->s);
    gsl_odeiv2_evolve_free(ode_s->e);
    gsl_odeiv2_driver_free(ode_s->d);
}

//----------------------------------------------------------------------------------------
// Update, reset
//----------------------------------------------------------------------------------------
/**
 *  \brief Reset an ode structure.
 *         Namely, the stepper s, the controller c, and the driver d.
 **/
void reset_ode_structure(OdeStruct *ode_s)
{
    gsl_odeiv2_step_reset(ode_s->s);
    gsl_odeiv2_evolve_reset(ode_s->e);
    gsl_odeiv2_driver_reset(ode_s->d);
}


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
                         double eps_root)
{
    ode_s->dim = d->sys->dimension;
    ode_s->T = T;
    ode_s->s = s;
    ode_s->c = c;
    ode_s->e = e;
    ode_s->sys = sys;
    ode_s->d = d;
    ode_s->s_root = s_root;

    ode_s->eps_root = eps_root;
    //ode_s->eps_diff = eps_diff;
    ode_s->h = d->h;

    //Integration tolerances are set to -1 because they are hidden in d and e (be careful during initialization!)
    ode_s->eps_int_abs = -1;
    ode_s->eps_int_rel = -1;
}

/**
 * \brief Set the right direction of integration in OdeStruct ode_s
 **/
void set_dir_ode_structure(OdeStruct *ode_s, double t0, double t1)
{
    int dir = (t0 < t1)? +1:-1;
    ode_s->d->h = dir*fabs(ode_s->d->h);
    ode_s->h = dir*fabs(ode_s->h);
}


/**
 * \brief Set the right direction of integration in a driver d
 **/
void set_dir_driver(gsl_odeiv2_driver *d, double t0, double t1)
{
    int dir = (t0 < t1)? +1:-1;
    d->h = dir*fabs(d->h);
}


