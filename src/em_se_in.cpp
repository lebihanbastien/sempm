/**
 * \file  eminsem.cpp
 * \brief Contains all the routines to perform changes of coordinates between:
 *          - the Earth-Moon (EM) synodical coordinates.
 *          - the Sun-Earth (SEM or SE) synodical coordinates.
 *          - the Normalized-Centered Earth-Moon (NCEM) coordinates.
 *          - the Normalized-Centered  Sun-Earth (NCSE or NCSEM) coordinates.
 *          - the Inertial (IN) coordinates, in EM units, SE units,
 *            or non-normalized units.
 *
 *        This files allows to handle either (position, velocity) or (position, momentum)
 *        format for the coordinates.
 *
 * \author BLB
 */

#include "em_se_in.h"


//----------------------------------------------------------------------------------------
// COC: Velocities <--> Momenta
//----------------------------------------------------------------------------------------
/**
 *  \brief Change the SE velocities into SE momenta
 **/
void se_v_to_se_m(double t, const double zv_se_v[], double zv_se_m[], void *params_void)
{
    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    FBPL* qbp = (FBPL*) params_void;
    double n  =  qbp->us_sem.n;

    //------------------------------------------------------------------------------------
    //evaluate the deltas
    //------------------------------------------------------------------------------------
    double delta[3];
    eval_array_coef(delta, t, n, qbp->n_order_fourier, qbp->cs_sem.coeffs, 3);

    //------------------------------------------------------------------------------------
    //Position to position
    //------------------------------------------------------------------------------------
    for(int i =0; i<3; i++) zv_se_m[i] = zv_se_v[i];

    //------------------------------------------------------------------------------------
    //Velocity to Momenta
    //------------------------------------------------------------------------------------
    //px
    zv_se_m[3] = 1.0/delta[0] * (zv_se_v[3] - delta[1]*zv_se_v[0] - delta[2]*zv_se_v[1]);
    //py
    zv_se_m[4] = 1.0/delta[0] * (zv_se_v[4] - delta[1]*zv_se_v[1] + delta[2]*zv_se_v[0]);
    //pz
    zv_se_m[5] = 1.0/delta[0] * (zv_se_v[5] - delta[1]*zv_se_v[2] );
}

/**
 *  \brief Change the SE momenta into SE velocities
 **/
void se_m_to_se_v(double t, const double zv_se_m[], double zv_se_v[], void *params_void)
{
    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    FBPL* qbp = (FBPL*) params_void;
    double n    =  qbp->us_sem.n;

    //------------------------------------------------------------------------------------
    //evaluate the deltas
    //------------------------------------------------------------------------------------
    double delta[3];
    eval_array_coef(delta, t, n, qbp->n_order_fourier, qbp->cs_sem.coeffs, 3);

    //------------------------------------------------------------------------------------
    //Position to position
    //------------------------------------------------------------------------------------
    for(int i =0; i<3; i++) zv_se_v[i] = zv_se_m[i];

    //------------------------------------------------------------------------------------
    //Momenta to velocities
    //------------------------------------------------------------------------------------
    //vx
    zv_se_v[3] = delta[0]*zv_se_m[3] + delta[1]*zv_se_m[0] + delta[2]*zv_se_m[1];
    //vy
    zv_se_v[4] = delta[0]*zv_se_m[4] + delta[1]*zv_se_m[1] - delta[2]*zv_se_m[0];
    //vz
    zv_se_v[5] = delta[0]*zv_se_m[5] + delta[1]*zv_se_m[2];
}

/**
 *  \brief Change the EM velocities into EM momenta
 **/
void em_v_to_em_m(double t, const double zv_em_v[], double zv_em_m[], void *params_void)
{
    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    FBPL* qbp  = (FBPL*) params_void;
    double n     =  qbp->us_em.n;

    //------------------------------------------------------------------------------------
    //Evaluate the alphas and their derivatives @ t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //------------------------------------------------------------------------------------
    double alpha[3];
    eval_array_coef(alpha, t, n, qbp->n_order_fourier, qbp->cs_em.coeffs, 3);

    //------------------------------------------------------------------------------------
    //Position to position
    //------------------------------------------------------------------------------------
    for(int i =0; i<3; i++) zv_em_m[i] = zv_em_v[i];


    //------------------------------------------------------------------------------------
    //Velocity to Momenta
    //------------------------------------------------------------------------------------
    //px
    zv_em_m[3] = 1.0/alpha[0] * (zv_em_v[3] - alpha[1]*zv_em_v[0] - alpha[2]*zv_em_v[1]);
    //py
    zv_em_m[4] = 1.0/alpha[0] * (zv_em_v[4] - alpha[1]*zv_em_v[1] + alpha[2]*zv_em_v[0]);
    //pz
    zv_em_m[5] = 1.0/alpha[0] * (zv_em_v[5] - alpha[1]*zv_em_v[2] );
}

/**
 *  \brief Change the EM momenta into EM velocities
 **/
void em_m_to_em_v(double t, const double zv_em_m[], double zv_em_v[], void *params_void)
{
    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    FBPL* qbp  = (FBPL*) params_void;
    double n     =  qbp->us_em.n;

    //------------------------------------------------------------------------------------
    //Evaluate the alphas and their derivatives @ t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //------------------------------------------------------------------------------------
    double alpha[3];
    eval_array_coef(alpha, t, n, qbp->n_order_fourier, qbp->cs_em.coeffs, 3);

    //------------------------------------------------------------------------------------
    //Position to position
    //------------------------------------------------------------------------------------
    for(int i =0; i<3; i++) zv_em_v[i] = zv_em_m[i];

    //------------------------------------------------------------------------------------
    //Momenta to velocities
    //------------------------------------------------------------------------------------
    //vx
    zv_em_v[3] = alpha[0]*zv_em_m[3] + alpha[1]*zv_em_m[0] + alpha[2]*zv_em_m[1];
    //vy
    zv_em_v[4] = alpha[0]*zv_em_m[4] + alpha[1]*zv_em_m[1] - alpha[2]*zv_em_m[0];
    //vz
    zv_em_v[5] = alpha[0]*zv_em_m[5] + alpha[1]*zv_em_m[2];
}

//----------------------------------------------------------------------------------------
// Change of unit system
//----------------------------------------------------------------------------------------
/**
 *   \brief From SE units to EM units for a (position, velocity) state vector and time
 *          given in IN coordinates.
 **/
void units_se_to_em(double *tc, double zv_in[], FBPL *fbpl)
{
    //IN[SE] to IN[EM]
    zv_in[0] *= fbpl->us_em.as;
    zv_in[1] *= fbpl->us_em.as;
    zv_in[2] *= fbpl->us_em.as;
    zv_in[3] *= fbpl->us_em.as*fbpl->us_em.ns;
    zv_in[4] *= fbpl->us_em.as*fbpl->us_em.ns;
    zv_in[5] *= fbpl->us_em.as*fbpl->us_em.ns;
    *tc     /= fbpl->us_em.ns;
}

/**
 *   \brief From EM units to SE units for a (position, velocity) state vector and time
 *          given in IN coordinates.
 **/
void units_em_to_se(double *tc, double zv_in[], FBPL *fbpl)
{
    //IN[EM] to IN[SE]
    zv_in[0] /= fbpl->us_em.as;
    zv_in[1] /= fbpl->us_em.as;
    zv_in[2] /= fbpl->us_em.as;
    zv_in[3] /= fbpl->us_em.as*fbpl->us_em.ns;
    zv_in[4] /= fbpl->us_em.as*fbpl->us_em.ns;
    zv_in[5] /= fbpl->us_em.as*fbpl->us_em.ns;
    *tc     *= fbpl->us_em.ns;
}

//----------------------------------------------------------------------------------------
// COC: IN <--> EM
//----------------------------------------------------------------------------------------
/**
 * \brief From EM to IN (in EM units)
 **/
void em_v_to_in(double t, const double zv_em_v[], double zv_in[], FBPL *fbpl)
{
    //------------------------------------------------------------------------------------
    //Init
    //------------------------------------------------------------------------------------
    //Param
    double n  = fbpl->us_em.n;
    double ms = fbpl->us_em.ms;
    double ns = fbpl->us_em.ns;
    double as = fbpl->us_em.as;
    double ni = fbpl->us_em.ni;
    double ai = fbpl->us_em.ai;
    //r
    double r1 = creal(eval_z(fbpl->cs_em.zt, t, n, ni, ai));
    double r2 = cimag(eval_z(fbpl->cs_em.zt, t, n, ni, ai));
    double r  = sqrt(r1*r1 + r2*r2);

    //R
    double R1 = creal(eval_z(fbpl->cs_em.Zt, t, n, ns, as));
    double R2 = cimag(eval_z(fbpl->cs_em.Zt, t, n, ns, as));
    //rdot
    double r1dot = creal(eval_zdot(fbpl->cs_em.zt, fbpl->cs_em.ztdot, t, n, ni, ai));
    double r2dot = cimag(eval_zdot(fbpl->cs_em.zt, fbpl->cs_em.ztdot, t, n, ni, ai));
    double rdot  = 1.0/r*(r1*r1dot + r2*r2dot);
    //Rdot
    double R1dot = creal(eval_zdot(fbpl->cs_em.Zt, fbpl->cs_em.Ztdot, t, n, ns, as));
    double R2dot = cimag(eval_zdot(fbpl->cs_em.Zt, fbpl->cs_em.Ztdot, t, n, ns, as));

    //------------------------------------------------------------------------------------
    // EM to IN: Position
    //------------------------------------------------------------------------------------
    zv_in[0] = r1*zv_em_v[0] - r2*zv_em_v[1] - ms/(1.0+ms)*R1;
    zv_in[1] = r2*zv_em_v[0] + r1*zv_em_v[1] - ms/(1.0+ms)*R2;
    zv_in[2] = r *zv_em_v[2];

    //------------------------------------------------------------------------------------
    // EM to IN: Velocity
    //------------------------------------------------------------------------------------
    zv_in[3] = r1dot*zv_em_v[0] - r2dot*zv_em_v[1] + r1*zv_em_v[3] - r2*zv_em_v[4]- ms/(1.0+ms)*R1dot;
    zv_in[4] = r2dot*zv_em_v[0] + r1dot*zv_em_v[1] + r2*zv_em_v[3] + r1*zv_em_v[4]- ms/(1.0+ms)*R2dot;
    zv_in[5] = rdot *zv_em_v[2] + r *zv_em_v[5];
}

/**
 * \brief From IN to EM (in EM units)
 **/
void in_to_em_v(double t, const double zv_in[], double zv_em_v[], FBPL *fbpl)
{
    //------------------------------------------------------------------------------------
    //Init
    //------------------------------------------------------------------------------------
    //Param
    double n  = fbpl->us_em.n;
    double ms = fbpl->us_em.ms;
    double ns = fbpl->us_em.ns;
    double as = fbpl->us_em.as;
    double ni = fbpl->us_em.ni;
    double ai = fbpl->us_em.ai;

    //r
    double r1 = creal(eval_z(fbpl->cs_em.zt, t, n, ni, ai));
    double r2 = cimag(eval_z(fbpl->cs_em.zt, t, n, ni, ai));
    double r = sqrt(r1*r1 + r2*r2);
    //R
    double R1 = creal(eval_z(fbpl->cs_em.Zt, t, n, ns, as));
    double R2 = cimag(eval_z(fbpl->cs_em.Zt, t, n, ns, as));
    //rdot
    double r1dot = creal(eval_zdot(fbpl->cs_em.zt, fbpl->cs_em.ztdot, t, n, ni, ai));
    double r2dot = cimag(eval_zdot(fbpl->cs_em.zt, fbpl->cs_em.ztdot, t, n, ni, ai));
    //Rdot
    double R1dot = creal(eval_zdot(fbpl->cs_em.Zt, fbpl->cs_em.Ztdot, t, n, ns, as));
    double R2dot = cimag(eval_zdot(fbpl->cs_em.Zt, fbpl->cs_em.Ztdot, t, n, ns, as));

    //Additional parameters
    double a = +pow(r, -4.0)*(r1dot*r*r - 2*r1*(r1*r1dot + r2*r2dot)); //dot(r1/(r*r))
    double b = +pow(r, -4.0)*(r2dot*r*r - 2*r2*(r1*r1dot + r2*r2dot)); //dot(r2/(r*r))
    double c = -pow(r, -3.0)*(r1*r1dot + r2*r2dot);                    //dot(1/r)

    //------------------------------------------------------------------------------------
    // EM to IN: Position
    //------------------------------------------------------------------------------------
    zv_em_v[0] = 1.0/(r*r) * ( +r1*(zv_in[0] + ms/(1.0+ms)*R1) +  r2*(zv_in[1] + ms/(1.0+ms)*R2) );
    zv_em_v[1] = 1.0/(r*r) * ( -r2*(zv_in[0] + ms/(1.0+ms)*R1) +  r1*(zv_in[1] + ms/(1.0+ms)*R2) );
    zv_em_v[2] = 1.0/r * zv_in[2];

    //------------------------------------------------------------------------------------
    // EM to IN: Velocity
    //------------------------------------------------------------------------------------
    zv_em_v[3] = +a*(zv_in[0] + ms/(1.0+ms)*R1) + b*(zv_in[1] + ms/(1.0+ms)*R2)
             + 1.0/(r*r) * ( +r1*(zv_in[3] + ms/(1.0+ms)*R1dot) +  r2*(zv_in[4] + ms/(1.0+ms)*R2dot) );
    zv_em_v[4] = -b*(zv_in[0] + ms/(1.0+ms)*R1) + a*(zv_in[1] + ms/(1.0+ms)*R2)
             + 1.0/(r*r) * ( -r2*(zv_in[3] + ms/(1.0+ms)*R1dot) +  r1*(zv_in[4] + ms/(1.0+ms)*R2dot) );
    zv_em_v[5] = c*zv_in[2] + 1.0/r * zv_in[5];
}

//----------------------------------------------------------------------------------------
// COC: IN <--> SE
//----------------------------------------------------------------------------------------
/**
 * \brief From SE to IN (in SE units)
 **/
void se_v_to_in(double t, const double zv_se_v[], double zv_in[], FBPL *fbpl)
{
    //------------------------------------------------------------------------------------
    //Init
    //------------------------------------------------------------------------------------
    //Param
    double n  = fbpl->us_sem.n;
    double ns = fbpl->us_sem.ns;
    double as = fbpl->us_sem.as;

    //------------------------------------------------------------------------------------
    //r & R
    //------------------------------------------------------------------------------------
    //R
    double R1 = creal(eval_z(fbpl->cs_sem.Zt, t, n, ns, as));
    double R2 = cimag(eval_z(fbpl->cs_sem.Zt, t, n, ns, as));
    double R = sqrt(R1*R1 + R2*R2);

    //------------------------------------------------------------------------------------
    //Derivatives
    //------------------------------------------------------------------------------------
    //Rdot
    double R1dot = creal(eval_zdot(fbpl->cs_sem.Zt, fbpl->cs_sem.Ztdot, t, n, ns, as));
    double R2dot = cimag(eval_zdot(fbpl->cs_sem.Zt, fbpl->cs_sem.Ztdot, t, n, ns, as));
    double Rdot  = 1.0/R*(R1*R1dot + R2*R2dot);

    //------------------------------------------------------------------------------------
    //Position & Velocity of the SE barycenter in inertial coordinates (in SE units)
    //------------------------------------------------------------------------------------
    double zv_in_bary[6];
    zv_in_bary[0] = 0.0;
    zv_in_bary[1] = 0.0;
    zv_in_bary[2] = 0.0;
    zv_in_bary[3] = 0.0;
    zv_in_bary[4] = 0.0;
    zv_in_bary[5] = 0.0;


    //------------------------------------------------------------------------------------
    // SE to IN: Position
    //------------------------------------------------------------------------------------
    zv_in[0] = R1*zv_se_v[0] - R2*zv_se_v[1] + zv_in_bary[0];
    zv_in[1] = R2*zv_se_v[0] + R1*zv_se_v[1] + zv_in_bary[1];
    zv_in[2] = R *zv_se_v[2];

    //------------------------------------------------------------------------------------
    // SE to IN: Velocity
    //------------------------------------------------------------------------------------
    zv_in[3] = R1dot*zv_se_v[0] - R2dot*zv_se_v[1] + R1*zv_se_v[3] - R2*zv_se_v[4] + zv_in_bary[3];
    zv_in[4] = R2dot*zv_se_v[0] + R1dot*zv_se_v[1] + R2*zv_se_v[3] + R1*zv_se_v[4] + zv_in_bary[4];
    zv_in[5] = Rdot *zv_se_v[2] + R *zv_se_v[5];
}

/**
 * \brief From IN to SE (in SE units)
 **/
void in_to_se_v(double t, const double zv_in[], double zv_se_v[], FBPL *fbpl)
{
    //------------------------------------------------------------------------------------
    //Init
    //------------------------------------------------------------------------------------
    //Param
    double n  = fbpl->us_sem.n;
    double ns = fbpl->us_sem.ns;
    double as = fbpl->us_sem.as;

    //------------------------------------------------------------------------------------
    //r & R
    //------------------------------------------------------------------------------------
    //R
    double R1 = creal(eval_z(fbpl->cs_sem.Zt, t, n, ns, as));
    double R2 = cimag(eval_z(fbpl->cs_sem.Zt, t, n, ns, as));
    //h
    double h1 = R1;
    double h2 = R2;
    double h  = sqrt(h1*h1 + h2*h2);

    //------------------------------------------------------------------------------------
    //Derivatives
    //------------------------------------------------------------------------------------
    //Rdot
    double R1dot = creal(eval_zdot(fbpl->cs_sem.Zt, fbpl->cs_sem.Ztdot, t, n, ns, as));
    double R2dot = cimag(eval_zdot(fbpl->cs_sem.Zt, fbpl->cs_sem.Ztdot, t, n, ns, as));
    //hdot
    double h1dot = R1dot;
    double h2dot = R2dot;

    //------------------------------------------------------------------------------------
    //Position & Velocity of the SE barycenter in inertial coordinates and SE units
    //------------------------------------------------------------------------------------
    double zv_in_bary[6];
    zv_in_bary[0] = 0.0;
    zv_in_bary[1] = 0.0;
    zv_in_bary[2] = 0.0;
    zv_in_bary[3] = 0.0;
    zv_in_bary[4] = 0.0;
    zv_in_bary[5] = 0.0;

    //Additional parameters
    double a = +pow(h, -4.0)*(h1dot*h*h - 2*h1*(h1*h1dot + h2*h2dot)); //dot(h1/(h*h))
    double b = +pow(h, -4.0)*(h2dot*h*h - 2*h2*(h1*h1dot + h2*h2dot)); //dot(h2/(h*h))
    double c = -pow(h, -3.0)*(h1*h1dot  + h2*h2dot);                   //dot(1/h)

    //------------------------------------------------------------------------------------
    // SE to IN: Position
    //------------------------------------------------------------------------------------
    zv_se_v[0] = 1.0/(h*h) * ( +h1*(zv_in[0] - zv_in_bary[0]) +  h2*(zv_in[1] - zv_in_bary[1]) );
    zv_se_v[1] = 1.0/(h*h) * ( -h2*(zv_in[0] - zv_in_bary[0]) +  h1*(zv_in[1] - zv_in_bary[1]) );
    zv_se_v[2] = 1.0/h * zv_in[2];

    //------------------------------------------------------------------------------------
    // SE to IN: Velocity
    //------------------------------------------------------------------------------------
    zv_se_v[3] = +a*(zv_in[0] - zv_in_bary[0]) + b*(zv_in[1] - zv_in_bary[1])
             + 1.0/(h*h) * ( +h1*(zv_in[3] - zv_in_bary[3]) +  h2*(zv_in[4] - zv_in_bary[4]) );
    zv_se_v[4] = -b*(zv_in[0] - zv_in_bary[0]) + a*(zv_in[1] - zv_in_bary[1])
             + 1.0/(h*h) * ( -h2*(zv_in[3] - zv_in_bary[3]) +  h1*(zv_in[4] - zv_in_bary[4]) );
    zv_se_v[5] = +c*zv_in[2] + 1.0/h * zv_in[5];
}

//----------------------------------------------------------------------------------------
// COC: SE <--> EM
//----------------------------------------------------------------------------------------
/**
 * \brief From SE to EM (both in position/momenta form)
 **/
void se_m_to_em_m(double t, const double zv_se_m[], double zv_em_m[], FBPL *fbpl)
{
    double tc = t;

    //Momenta to velocities
    double zv_se_v[6];
    se_m_to_se_v(tc, zv_se_m, zv_se_v, fbpl);

    //SE --> IN
    double zv_in[6];
    se_v_to_in(tc, zv_se_v, zv_in, fbpl);

    //IN[SE] to IN[EM]
    units_se_to_em(&tc, zv_in, fbpl);

    //IN --> EM
    double zv_em_v[6];
    in_to_em_v(tc, zv_in, zv_em_v, fbpl);

    //Velocities to momenta
    em_v_to_em_m(tc, zv_em_v, zv_em_m, fbpl);
}

/**
 * \brief From EM to SE (both in position/momenta form)
 **/
void em_m_to_se_m(double t, const double zv_em_m[], double zv_se_m[], FBPL *fbpl)
{
    double tc = t;
    //Momenta to velocities
    double zv_em_v[6];
    em_m_to_em_v(tc, zv_em_m, zv_em_v, fbpl);

    //EM-->IN
    double zv_in[6];
    em_v_to_in(tc, zv_em_v, zv_in, fbpl);

    //IN[EM] to IN[SE]
    units_em_to_se(&tc, zv_in, fbpl);

    //IN-->SE
    double zv_se_v[6];
    in_to_se_v(tc, zv_in, zv_se_v, fbpl);

    //Velocities to momenta
    se_v_to_se_m(tc, zv_se_v, zv_se_m, fbpl);
}

//----------------------------------------------------------------------------------------
// COC: SE <--> NCEM
//----------------------------------------------------------------------------------------
/**
 * \brief From NCEM to SE (both in position/momenta form)
 **/
void ncem_m_to_se_m(double t, const double zv_ncem_m[], double zv_se_m[], FBPL *fbpl)
{
    double zv_em_m[6];
    //NC to EM
    ncem_m_to_em_m(t, zv_ncem_m, zv_em_m, fbpl);
    //EM to SE
    em_m_to_se_m(t, zv_em_m, zv_se_m, fbpl);
}


//----------------------------------------------------------------------------------------
// COC: NCSEM <--> NCEM
//----------------------------------------------------------------------------------------
/**
 * \brief From NCSE to  NCEM (both in position/momenta form)
 **/
void ncse_m_to_ncem_m(double t, const double zv_ncse_m[], double zv_ncem_m[], FBPL *fbpl)
{
    double zv_se_m[6], zv_em_m[6];
    //NC to SE
    ncse_m_to_se_m(t, zv_ncse_m, zv_se_m, fbpl);
    //SE to EM
    se_m_to_em_m(t, zv_se_m, zv_em_m, fbpl);
    //EM to NC (careful, the time should be set into EM units!)
    em_m_to_ncem_m(t/fbpl->us_em.ns, zv_em_m, zv_ncem_m, fbpl);
}

/**
 * \brief From NCEM to  NCSE (both in position/momenta form)
 **/
void ncem_m_to_ncse_m(double tEM, const double zv_ncem_m[], double zv_ncse_m[], FBPL *fbpl)
{
    double zv_em_m[6], zv_se_m[6];
    //NC to EM
    ncem_m_to_em_m(tEM, zv_ncem_m, zv_em_m, fbpl);
    //EM to SE
    em_m_to_se_m(tEM, zv_em_m, zv_se_m, fbpl);
    //SE to NC (careful, the time should be set into SE units!)
    se_m_to_ncse_m(tEM*fbpl->us_em.ns, zv_se_m, zv_ncse_m, fbpl);
}


//----------------------------------------------------------------------------------------
// COC: NCSYS <--> SYS
//----------------------------------------------------------------------------------------
/**
 *  \brief COC: NC(EM or SE) coordinates to SYS (EM or SE) coordinates.
 **/
void ncsys_m_to_sys_m(double t, const double zv_ncsys_m[], double zv_sys_m[], FBPL *qbp)
{
    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    double n     =  qbp->us.n;
    double gamma =  qbp->cs.gamma;
    double c1    =  qbp->cs.c1;
    //------------------------------------------------------------------------------------
    //Evaluate the alphas and their derivatives @ t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //------------------------------------------------------------------------------------
    double alpha[8];
    eval_array_coef(alpha, t, n, qbp->n_order_fourier, qbp->cs.coeffs, 8);

    //------------------------------------------------------------------------------------
    //CoC
    //------------------------------------------------------------------------------------
    //X = -gamma*(x-c1)
    zv_sys_m[0] = -gamma*(zv_ncsys_m[0] - c1);
    //Y = -gamma*y
    zv_sys_m[1] = -gamma*zv_ncsys_m[1];
    //Z = +gamma*z
    zv_sys_m[2] = +gamma*zv_ncsys_m[2];

    //PX = -gamma(px+a2/a1*c1)
    zv_sys_m[3] = -gamma*(zv_ncsys_m[3] + alpha[1]/alpha[0]*c1);
    //PY = -gamma(py-a3/a1*c1)
    zv_sys_m[4] = -gamma*(zv_ncsys_m[4] - alpha[2]/alpha[0]*c1);
    //PZ = +gamma*pz
    zv_sys_m[5] = +gamma*zv_ncsys_m[5];
}

/**
 *  \brief COC: from SYS (EM or SE) coordinates to NC(EM or SE) coordinates.
 **/
void sys_m_to_ncsys_m(double t, const double zv_sys_m[], double zv_ncsys_m[], FBPL *qbp)
{
    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    double n     =  qbp->us.n;
    double gamma =  qbp->cs.gamma;
    double c1    =  qbp->cs.c1;

    //------------------------------------------------------------------------------------
    //Evaluate the alphas and their derivatives @ t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //------------------------------------------------------------------------------------
    double alpha[8];
    eval_array_coef(alpha, t, n, qbp->n_order_fourier, qbp->cs_em.coeffs, 8);

    //------------------------------------------------------------------------------------
    //CoC
    //------------------------------------------------------------------------------------
    //x = -X/gamma + c1
    zv_ncsys_m[0] = -zv_sys_m[0]/gamma  + c1;
    //y = -Y/gamma
    zv_ncsys_m[1] = -zv_sys_m[1]/gamma;
    //z = +Z/gamma
    zv_ncsys_m[2] = +zv_sys_m[2]/gamma;
    //px = -PX/gamma - a2/a1*c1
    zv_ncsys_m[3] = -zv_sys_m[3]/gamma - alpha[1]/alpha[0]*c1;
    //py = -PY/gamma + a3/a1*c1
    zv_ncsys_m[4] = -zv_sys_m[4]/gamma + alpha[2]/alpha[0]*c1;
    //pz = +PZ/gamma
    zv_ncsys_m[5] = +zv_sys_m[5]/gamma;
}


//----------------------------------------------------------------------------------------
// COC: NCEM <--> EM
//----------------------------------------------------------------------------------------
/**
 *  \brief COC: from NCEM coordinates to EM coordinates
 **/
void ncem_m_to_em_m(double t, const double zv_ncem_m[], double zv_em_m[], FBPL *qbp)
{
    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    double n     =  qbp->us_em.n;
    double gamma =  qbp->cs_em.gamma;
    double c1    =  qbp->cs_em.c1;

    //------------------------------------------------------------------------------------
    //Evaluate the alphas and their derivatives @ t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //------------------------------------------------------------------------------------
    double alpha[8];
    eval_array_coef(alpha, t, n, qbp->n_order_fourier, qbp->cs_em.coeffs, 8);

    //------------------------------------------------------------------------------------
    //CoC
    //------------------------------------------------------------------------------------
    //X = -gamma*(x-c1)
    zv_em_m[0] = -gamma*(zv_ncem_m[0] - c1);
    //Y = -gamma*y
    zv_em_m[1] = -gamma*zv_ncem_m[1];
    //Z = +gamma*z
    zv_em_m[2] = +gamma*zv_ncem_m[2];

    //PX = -gamma(px+a2/a1*c1)
    zv_em_m[3] = -gamma*(zv_ncem_m[3] + alpha[1]/alpha[0]*c1);
    //PY = -gamma(py-a3/a1*c1)
    zv_em_m[4] = -gamma*(zv_ncem_m[4] - alpha[2]/alpha[0]*c1);
    //PZ = +gamma*pz
    zv_em_m[5] = +gamma*zv_ncem_m[5];

}

/**
 *  \brief COC: from EM coordinates to NCEM coordinates
 **/
void em_m_to_ncem_m(double t, const double zv_em_m[], double zv_ncem_m[], FBPL *qbp)
{
    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    double n     =  qbp->us_em.n;
    double gamma =  qbp->cs_em.gamma;
    double c1    =  qbp->cs_em.c1;

    //------------------------------------------------------------------------------------
    //Evaluate the alphas and their derivatives @ t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //------------------------------------------------------------------------------------
    double alpha[8];
    eval_array_coef(alpha, t, n, qbp->n_order_fourier, qbp->cs_em.coeffs, 8);

    //------------------------------------------------------------------------------------
    //CoC
    //------------------------------------------------------------------------------------
    //x = -X/gamma + c1
    zv_ncem_m[0] = -zv_em_m[0]/gamma  + c1;
    //y = -Y/gamma
    zv_ncem_m[1] = -zv_em_m[1]/gamma;
    //z = +Z/gamma
    zv_ncem_m[2] = +zv_em_m[2]/gamma;
    //px = -PX/gamma - a2/a1*c1
    zv_ncem_m[3] = -zv_em_m[3]/gamma - alpha[1]/alpha[0]*c1;
    //py = -PY/gamma + a3/a1*c1
    zv_ncem_m[4] = -zv_em_m[4]/gamma + alpha[2]/alpha[0]*c1;
    //pz = +PZ/gamma
    zv_ncem_m[5] = +zv_em_m[5]/gamma;
}


//----------------------------------------------------------------------------------------
// COC: NCSE <--> SE
//----------------------------------------------------------------------------------------
/**
 *  \brief COC: from SE coordinates to NCSE coordinates
 **/
void se_m_to_ncse_m(double t, const double zv_se_m[], double zv_ncse_m[], FBPL *qbp)
{
    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    double n     =  qbp->us_sem.n;
    double gamma =  qbp->cs_sem.gamma;
    double c1    =  qbp->cs_sem.c1;

    //------------------------------------------------------------------------------------
    //Evaluate the alphas and their derivatives @ t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //------------------------------------------------------------------------------------
    double alpha[8];
    eval_array_coef(alpha, t, n, qbp->n_order_fourier, qbp->cs_sem.coeffs, 8);

    //------------------------------------------------------------------------------------
    //CoC
    //------------------------------------------------------------------------------------
    //x = -X/gamma + c1
    zv_ncse_m[0] = -zv_se_m[0]/gamma  + c1;
    //y = -Y/gamma
    zv_ncse_m[1] = -zv_se_m[1]/gamma;
    //z = +Z/gamma
    zv_ncse_m[2] = +zv_se_m[2]/gamma;
    //px = -PX/gamma - a2/a1*c1
    zv_ncse_m[3] = -zv_se_m[3]/gamma - alpha[1]/alpha[0]*c1;
    //py = -PY/gamma + a3/a1*c1
    zv_ncse_m[4] = -zv_se_m[4]/gamma + alpha[2]/alpha[0]*c1;
    //pz = +PZ/gamma
    zv_ncse_m[5] = +zv_se_m[5]/gamma;
}

/**
 *  \brief COC: from NCSE coordinates to SE coordinates
 **/
void ncse_m_to_se_m(double t, const double zv_ncse_m[], double zv_se_m[], FBPL *qbp)
{
    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    double n     =  qbp->us_sem.n;
    double gamma =  qbp->cs_sem.gamma;
    double c1    =  qbp->cs_sem.c1;

    //------------------------------------------------------------------------------------
    //Evaluate the alphas and their derivatives @ t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //------------------------------------------------------------------------------------
    double alpha[8];
    eval_array_coef(alpha, t, n, qbp->n_order_fourier, qbp->cs_sem.coeffs, 8);

    //------------------------------------------------------------------------------------
    //CoC
    //------------------------------------------------------------------------------------
    //X = -gamma*(x-c1)
    zv_se_m[0] = -gamma*(zv_ncse_m[0] - c1);
    //Y = -gamma*y
    zv_se_m[1] = -gamma*zv_ncse_m[1];
    //Z = +gamma*z
    zv_se_m[2] = +gamma*zv_ncse_m[2];

    //PX = -gamma(px+a2/a1*c1)
    zv_se_m[3] = -gamma*(zv_ncse_m[3] + alpha[1]/alpha[0]*c1);
    //PY = -gamma(py-a3/a1*c1)
    zv_se_m[4] = -gamma*(zv_ncse_m[4] - alpha[2]/alpha[0]*c1);
    //PZ = +gamma*pz
    zv_se_m[5] = +gamma*zv_ncse_m[5];
}

