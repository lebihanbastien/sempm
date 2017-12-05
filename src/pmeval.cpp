#include "pmeval.h"

/**
 * \file pmcoc.cpp
 * \brief Evaluation of the high-ofts_order semi-analytical approximations,
 *        along with inner change of coordinates between different coordinates system
 *        for the reduced state s. I This files also deals with change of variables
 *        between Translated-Floquet-Complexified (TFC) and Translated-Floquet (TF)
 *        coordinates for the 6-dim state.
 * \author BLB.
 *
 *
 *  The reduced state s is a d-dim vector (d = 4, 5, or 6), which can be expressed
 *  in :
 *
 *  - RCM (Real Center Manifold) coordinates, as an array s[0:d-1] of double;
 *  - CCM (Complex Center Manifold) coordinates, as an array s[0:d-1] of double complex;
 *  - CCM8 (Complex Center Manifold) coordinates, as an array s[0:2*d-1] of double, which
 * correspond to CCM coordinates for which the real and imaginary parts are separated, for
 * real numerical integration purposes.
 *
 */

//----------------------------------------------------------------------------------------
//
//          Change of coordinates for the reduced state
//
//----------------------------------------------------------------------------------------
/**
 *  \brief from CCM to RCM coordinates
 **/
void ccm_to_rcm(const cdouble sv_ccm[], double sv_rcm[], int nv)
{
    sv_rcm[0] = creal(1.0/sqrt(2)*(sv_ccm[0]   + sv_ccm[2]*I));
    sv_rcm[2] = creal(1.0/sqrt(2)*(sv_ccm[0]*I + sv_ccm[2]));
    sv_rcm[1] = creal(1.0/sqrt(2)*(sv_ccm[1]   + sv_ccm[3]*I));
    sv_rcm[3] = creal(1.0/sqrt(2)*(sv_ccm[1]*I + sv_ccm[3]));
    if(nv == 5)
    {
        sv_rcm[4] = creal(sv_ccm[4]);
    }
    else if(nv == 6)
    {
        sv_rcm[4] = creal(sv_ccm[4]);
        sv_rcm[5] = creal(sv_ccm[5]);
    }
}

/**
 *  \brief from RCM to CCM coordinates
 **/
void rcm_to_ccm(const double sv_rcm[], cdouble sv_ccm[], int nv)
{
    //From real to complex TFC
    sv_ccm[0] = 1.0/sqrt(2)*(sv_rcm[0] - sv_rcm[2]*I);
    sv_ccm[2] = 1.0/sqrt(2)*(sv_rcm[2] - sv_rcm[0]*I);
    sv_ccm[1] = 1.0/sqrt(2)*(sv_rcm[1] - sv_rcm[3]*I);
    sv_ccm[3] = 1.0/sqrt(2)*(sv_rcm[3] - sv_rcm[1]*I);
    if(nv == 5)
    {
        sv_ccm[4] = sv_rcm[4]+I*0.0;
    }
    else if(nv == 6)
    {
        sv_ccm[4] = sv_rcm[4]+I*0.0;
        sv_ccm[5] = sv_rcm[5]+I*0.0;
    }
}

/**
 *  \brief from RCM to CCM coordinates, with real and imag part stored separately
 **/
void rcm_to_cmm8(const double sv_rcm[], double sv_ccm8[])
{
    //From real to complex TFC
    cdouble sv_ccm[REDUCED_NV];
    rcm_to_ccm(sv_rcm, sv_ccm, REDUCED_NV);

    //Store real and imag part separately
    ccm_to_cmm8(sv_ccm, sv_ccm8);
}

/**
 *  \brief from CCM coordinates, with real and imag part stored separately, to RCM coordinates
 **/
void ccm8_to_rcm(const double sv_ccm8[], double sv_rcm[])
{
    //CCM8 to CCM
    cdouble sv_ccm[REDUCED_NV];
    ccm8_to_cmm(sv_ccm8, sv_ccm);
    //CCM to RCM
    ccm_to_rcm(sv_ccm, sv_rcm, REDUCED_NV);
}

/**
 *  \brief from CCM coordinates, with real and imag part stored separately, to CCM coordinates
 **/
void ccm8_to_cmm(const double sv_ccm8[], cdouble sv_ccm[])
{
    int p2 = 0;
    for(int p = 0; p < REDUCED_NV; p++)
    {
        sv_ccm[p]  =   sv_ccm8[p2++];
        sv_ccm[p] += I*sv_ccm8[p2++];
    }
}

/**
 *  \brief from CCM coordinates to CCM coordinates, with real and imag part stored separately.
 **/
void ccm_to_cmm8(const cdouble sv_ccm[], double sv_ccm8[])
{
    //Store real and imag part separately
    int p2 = 0;
    for(int p = 0; p < REDUCED_NV; p++)
    {
        sv_ccm8[p2++] = creal(sv_ccm[p]);
        sv_ccm8[p2++] = cimag(sv_ccm[p]);
    }
}

//----------------------------------------------------------------------------------------
//
//          Change of coordinates for the 6-dimensional state
//
//----------------------------------------------------------------------------------------
/**
 *  \brief from TFC to TF coordinates
 **/
void tfc_to_tf(const cdouble zv_ccm[6], double zv_tf[6])
{
    //First center
    zv_tf[0] = creal(1.0/sqrt(2)*(zv_ccm[0]   + zv_ccm[3]*I));
    zv_tf[3] = creal(1.0/sqrt(2)*(zv_ccm[0]*I + zv_ccm[3]));
    //Second center
    zv_tf[2] = creal(1.0/sqrt(2)*(zv_ccm[2]   + zv_ccm[5]*I));
    zv_tf[5] = creal(1.0/sqrt(2)*(zv_ccm[2]*I + zv_ccm[5]));
    //Hyperbolic dir
    zv_tf[1] = creal(zv_ccm[1]);
    zv_tf[4] = creal(zv_ccm[4]);
}

//----------------------------------------------------------------------------------------
//
//          Evaluation of the high-order pm
//
//----------------------------------------------------------------------------------------
/**
 *   \brief Evaluate the TFC configuration zv_tfc(t) = Wh(sv_ccm, t)
 *   \param sv_ccm an array of complex which gives the configuration to input in CCM coordinates
 *   \param ofts_order the order of the Taylor expansions to use in the evaluation
 *   \param ofs_order the order of the Fourier expansions to use in the evaluation
 *   \param Wh the Fourier-Taylor expansion that contains the parameterization Wh(s,t), in TFC coordinates
 *   \param zv_tfc the 6-dim output vector in TFC coordinates, given as vector<Ofsc> object
 *   \param graph_style if true, the special case of the QBCP is used to compute Wh.
 *          Namely: only Wh[1] and Wh[4] are complete polynomials, the rest of the components are of order one.
 **/
void ccm_to_tfc(cdouble sv_ccm[],
                const int ofts_order,
                const int ofs_order,
                vector<Oftsc>& Wh,
                vector<Ofsc>& zv_tfc,
                bool graph_style)
{
    //------------------------------------------------------------------------------------
    // 2. Update zv_tfc
    //------------------------------------------------------------------------------------
    if(graph_style)
    {
        //--------------------------------------------------------------------------------
        // Using particular geometry of the QBCP
        //--------------------------------------------------------------------------------
        cdouble temp;
        // zv_tfc[0]
        //---------------
        temp = Wh[0].get_coef(1,0)->ofs_get_coef(0);
        zv_tfc[0].set_coef(temp*sv_ccm[0],0);
        // zv_tfc[1]
        //---------------
        Wh[1].evaluate(sv_ccm, zv_tfc[1], ofts_order, ofs_order);
        // zv_tfc[2]
        //---------------
        temp = Wh[2].get_coef(1,1)->ofs_get_coef(0);
        zv_tfc[2].set_coef(temp*sv_ccm[1],0);
        // zv_tfc[3]
        //---------------
        temp = Wh[3].get_coef(1,2)->ofs_get_coef(0);
        zv_tfc[3].set_coef(temp*sv_ccm[2],0);
        // zv_tfc[4]
        //---------------
        Wh[4].evaluate(sv_ccm, zv_tfc[4], ofts_order, ofs_order);
        // zv_tfc[5]
        //---------------
        temp = Wh[5].get_coef(1,3)->ofs_get_coef(0);
        zv_tfc[5].set_coef(temp*sv_ccm[3],0);
    }
    else
    {
        //--------------------------------------------------------------------------------
        // General computation
        //--------------------------------------------------------------------------------
        for(int p = 0; p < 6; p++)
        {
            Wh[p].evaluate(sv_ccm, zv_tfc[p], ofts_order, ofs_order);
        }
    }
}


/**
 *   \brief Evaluate the TFC configuration zv_tfc(t) = Wh(g(sv_rcm), t)
 *   \param sv_rcm an array of double which gives the configuration to input in real CM coordinates
 *   \param ofts_order the order of the Taylor expansions to use in the evaluation
 *   \param ofs_order the order of the Fourier expansions to use in the evaluation
 *   \param Wh the Fourier-Taylor expansion that contains the parameterization Wh(s,t), in TFC coordinates
 *   \param zv_tfc the 6-dim output vector in TFC coordinates, given as vector<Ofsc> object
 *   \param graph_style if true, the special case of the QBCP is used to compute Wh.
 *          Namely: only Wh[1] and Wh[4] are complete polynomials, the rest of the components are of order one.
 **/
void rcm_to_tfc(const double sv_rcm[],
                const int ofts_order,
                const int ofs_order,
                vector<Oftsc>& Wh,
                vector<Ofsc>& zv_tfc,
                bool graph_style)
{
    //------------------------------------------------------------------------------------
    // Inner variables (CCM)
    //------------------------------------------------------------------------------------
    cdouble sv_ccm[REDUCED_NV];

    //------------------------------------------------------------------------------------
    // 1. RCM to CCM
    //------------------------------------------------------------------------------------
    rcm_to_ccm(sv_rcm, sv_ccm, REDUCED_NV);

    //------------------------------------------------------------------------------------
    // 2. Update zv_tfc
    //------------------------------------------------------------------------------------
    ccm_to_tfc(sv_ccm, ofts_order, ofs_order, Wh, zv_tfc, graph_style);
}



/**
 *   \brief Evaluate the TFC configuration zv_tfc(t) = Wh(g(s_ccm8), t)
 *   \param s_ccm8 an array of double which gives the configuration to input in CCM8 coordinates
 *   \param ofts_order the order of the Taylor expansions to use in the evaluation
 *   \param ofs_order the order of the Fourier expansions to use in the evaluation
 *   \param Wh the Fourier-Taylor expansion that contains the parameterization Wh(s,t), in TFC coordinates
 *   \param zv_tfc the 6-dim output vector in TFC coordinates, given as vector<Ofsc> object
 *   \param graph_style if true, the special case of the QBCP is used to compute Wh.
 *          Namely: only Wh[1] and Wh[4] are complete polynomials, the rest of the components are of order one.
 **/
void ccm8_to_tfc(const double s_ccm8[],
                 const int ofts_order,
                 const int ofs_order,
                 vector<Oftsc>& Wh,
                 vector<Ofsc>& zv_tfc,
                 bool graph_style)
{
    //------------------------------------------------------------------------------------
    // Inner variables (CCM)
    //------------------------------------------------------------------------------------
    cdouble sv_ccm[REDUCED_NV];

    //------------------------------------------------------------------------------------
    // 1. CCM8 to CCM
    //------------------------------------------------------------------------------------
    ccm8_to_cmm(s_ccm8, sv_ccm);

    //------------------------------------------------------------------------------------
    // 2. Update zv_tfc
    //------------------------------------------------------------------------------------
    ccm_to_tfc(sv_ccm, ofts_order, ofs_order, Wh, zv_tfc, graph_style);
}

/**
 *   \brief Evaluate the configuration zh1 = Wh(g(sv_rcm), n*t) in TF coordinates (not complex!)
 *   \param sv_rcm an array of double which gives the configuration to input in real CM coordinates
 *   \param t the current time
 *   \param n the pulsation of the system
 *   \param ofts_order the order of the Taylor expansions to use in the evaluation
 *   \param ofs_order the order of the Fourier expansions to use in the evaluation
 *   \param Wh the Fourier-Taylor expansion that contains the parameterization Wh(s,t), in TFC coordinates
 *   \param ofs a side variable to compute OFS objects
 *   \param zv_tf the 6-dim output vector in TF coordinates, given as a vector of double
 *   \param graph_style if true, the special case of the QBCP is used to compute W.
 *          Namely: only Wh[1] and Wh[4] are complete polynomials, the rest of the components are of order one.
 **/
void rcm_to_tf(const double sv_rcm[],
               const double t,
               const double n,
               const int ofts_order,
               const int ofs_order,
               vector<Oftsc>& Wh,
               Ofsc& ofs,
               double z_tf[],
               bool graph_style)
{
    //------------------------------------------------------------------------------------
    // Inner variables (TFC)
    //------------------------------------------------------------------------------------
    cdouble zv_tfc[6];
    //------------------------------------------------------------------------------------
    // Update z0
    //------------------------------------------------------------------------------------
    rcm_to_tfc(sv_rcm, t, n, ofts_order, ofs_order, Wh, ofs, zv_tfc, graph_style);
    //------------------------------------------------------------------------------------
    // TFC to TF
    //------------------------------------------------------------------------------------
    tfc_to_tf(zv_tfc, z_tf);
}



/**
 *   \brief Evaluate the configuration zh1 = Wh(sv_ccm, n*t)
 *   \param sv_ccm an array of cdouble which gives the configuration to input in complex CM coordinates
 *   \param t the current time
 *   \param n the pulsation of the system
 *   \param ofts_order the order of the Taylor expansions to use in the evaluation
 *   \param ofs_order the order of the Fourier expansions to use in the evaluation
 *   \param Wh the Fourier-Taylor expansion that contains the parameterization Wh(s,t), in TFC coordinates
 *   \param ofs a side variable to compute OFS objects
 *   \param zv_tfc the output array to update, in TFC coordinates
 *   \param graph_style if true, the special case of the QBCP is used to compute W.
 *          Namely: only Wh[1] and Wh[4] are complete polynomials, the rest of the components are of order one.
 **/
void ccm_to_tfc(cdouble sv_ccm[],
                const double t,
                const double n,
                const int ofts_order,
                const int ofs_order,
                vector<Oftsc>& Wh,
                Ofsc& ofs,
                cdouble zv_tfc[],
                bool graph_style)
{
    //------------------------------------------------------------------------------------
    // 1. Update zIn
    //------------------------------------------------------------------------------------
    if(graph_style)
    {
        //--------------------------------------------------------------------------------
        // Using particular geometry of the QBCP
        //--------------------------------------------------------------------------------
        for(int p = 0; p < 6; p++)
        {
            if(p == 1 || p == 4)
            {
                //For p = 1,4 normal computation
                Wh[p].evaluate(sv_ccm, ofs, ofts_order, ofs_order);
                zv_tfc[p] = ofs.evaluate(n*t, ofs_order);
            }
            else
            {
                //ofts_order 1 is sufficient for p = 0,2,3,5
                Wh[p].evaluate(sv_ccm, ofs, 1, ofs_order);
                zv_tfc[p] = ofs.evaluate(n*t, ofs_order);
            }
        }
    }
    else
    {
        //--------------------------------------------------------------------------------
        // General computation
        //--------------------------------------------------------------------------------
        for(int p = 0; p < 6; p++)
        {
            Wh[p].evaluate(sv_ccm, ofs, ofts_order, ofs_order);
            zv_tfc[p] = ofs.evaluate(n*t, ofs_order);
        }
    }
}


/**
 *   \brief Evaluate the configuration zh1 = Wh(g(sv_rcm), n*t)
 *   \param sv_rcm an array of 4 double which gives the configuration to input in real CM coordinates
 *   \param t the current time
 *   \param n the pulsation of the system
 *   \param ofts_order the order of the Taylor expansions to use in the evaluation
 *   \param ofs_order the order of the Fourier expansions to use in the evaluation
 *   \param Wh the Fourier-Taylor expansion that contains the parameterization Wh(s,t), in TFC coordinates
 *   \param ofs a side variable to compute OFS objects
 *   \param zv_tfc the output array to update, in TFC coordinates
 *   \param graph_style if true, the special case of the QBCP is used to compute W.
 *          Namely: only Wh[1] and Wh[4] are complete polynomials, the rest of the components are of order one.
 **/
void rcm_to_tfc(const double sv_rcm[],
                const double t,
                const double n,
                const int ofts_order,
                const int ofs_order,
                vector<Oftsc>& Wh,
                Ofsc& ofs,
                cdouble zv_tfc[],
                bool graph_style)
{
    //------------------------------------------------------------------------------------
    // Inner variables (CCM)
    //------------------------------------------------------------------------------------
    cdouble sv_ccm[REDUCED_NV];

    //------------------------------------------------------------------------------------
    // RCM to CCM
    //------------------------------------------------------------------------------------
    rcm_to_ccm(sv_rcm, sv_ccm, REDUCED_NV);

    //------------------------------------------------------------------------------------
    // 2. Update zIn
    //------------------------------------------------------------------------------------
    ccm_to_tfc(sv_ccm, t, n, ofts_order, ofs_order, Wh, ofs, zv_tfc, graph_style);
}


/**
 *   \brief Evaluate the configuration zh1 = Wh(g(sv_ccm8), n*t)
 *   \param sv_ccm8 an array of 8 double which gives the configuration to input in complex CM coordinates
 *   \param t the current time
 *   \param n the pulsation of the system
 *   \param ofts_order the order of the Taylor expansions to use in the evaluation
 *   \param ofs_order the order of the Fourier expansions to use in the evaluation
 *   \param Wh the Fourier-Taylor expansion that contains the parameterization Wh(s,t), in TFC coordinates
 *   \param ofs a side variable to compute OFS objects
 *   \param zv_tfc the output array to update, in TFC coordinates
 *   \param graph_style if true, the special case of the QBCP is used to compute W.
 *          Namely: only Wh[1] and Wh[4] are complete polynomials, the rest of the components are of order one.
 **/
void ccm8_to_tfc(const double sv_ccm8[],
                 const double t,
                 const double n,
                 const int ofts_order,
                 const int ofs_order,
                 vector<Oftsc>& Wh,
                 Ofsc& ofs,
                 cdouble zv_tfc[],
                 bool graph_style)
{
    //------------------------------------------------------------------------------------
    // Inner variables (CCM)
    //------------------------------------------------------------------------------------
    cdouble sv_ccm[REDUCED_NV];

    //------------------------------------------------------------------------------------
    // CCM8 to CCM
    //------------------------------------------------------------------------------------
    ccm8_to_cmm(sv_ccm8, sv_ccm);

    //------------------------------------------------------------------------------------
    // 2. Update zIn
    //------------------------------------------------------------------------------------
    ccm_to_tfc(sv_ccm, t, n, ofts_order, ofs_order, Wh, ofs, zv_tfc, graph_style);
}



/**
 *   \brief Evaluate the configuration zv_nc = W(g(sv_rcm), n*t) with the use of an
 *          intermediate TFC configuration zv_tfc(t) = Wh(g(sv_rcm), t).
 *   \param sv_rcm an array of double which gives the configuration to input in RCM coordinates
 *   \param t the current time
 *   \param n the pulsation of the system
 *   \param ofts_order the order of the Taylor expansions to use in the evaluation
 *   \param ofs_order the order of the Fourier expansions to use in the evaluation
 *   \param Wh the Fourier-Taylor expansion that contains the parameterization Wh(s,t), in TFC coordinates
 *   \param ofs a side variable to compute OFS objects
 *   \param PC COC matrix: z = PC*zh+V
 *   \param V COC vector: z = PC*zh+V
 *   \param zv_nc the output array to update, in NC coordinates
 *   \param graph_style if true, the special case of the QBCP is used to compute Wh.
 *          Namely: only Wh[1] and Wh[4] are complete polynomials, the rest of the
 *          components are of order one.
 **/
void rcm_to_nc_by_tfc(const double sv_rcm[],
                      const double t,
                      const double n,
                      const int ofts_order,
                      const int ofs_order,
                      vector<Oftsc>& Wh,
                      matrix<Ofsc>& PC,
                      vector<Ofsc>& V,
                      double zv_nc[],
                      bool graph_style)
{
    //------------------------------------------------------------------------------------
    // Inner variables (CCM, TFC, NC)
    //------------------------------------------------------------------------------------
    cdouble z_nc_comp[6];
    vector<Ofsc> zv_tfc(6), zv_nc_ofs(6);

    //------------------------------------------------------------------------------------
    // RCM to TFC
    //------------------------------------------------------------------------------------
    rcm_to_tfc(sv_rcm, ofts_order, ofs_order, Wh, zv_tfc, graph_style);

    //------------------------------------------------------------------------------------
    // TFC to NC
    //------------------------------------------------------------------------------------
    apply_coc_ofs(PC, V, zv_tfc, zv_nc_ofs);
    for(int p = 0; p < 6; p++)
    {
        z_nc_comp[p] = zv_nc_ofs[p].evaluate(n*t, ofs_order);
        zv_nc[p] = creal(z_nc_comp[p]);
    }
}

/**
 *   \brief Evaluate the configuration zv_nc = W(g(sv_ccm8), n*t) with the use of an
 *          intermediate TFC configuration zv_tfc(t) = Wh(g(sv_ccm8), t)
 *   \param sv_ccm8 an array of double which gives the configuration to input in CCM8 coordinates
 *   \param t the current time
 *   \param n the pulsation of the system
 *   \param ofts_order the order of the Taylor expansions to use in the evaluation
 *   \param ofs_order the order of the Fourier expansions to use in the evaluation
 *   \param Wh the Fourier-Taylor expansion that contains the parameterization Wh(s,t), in TFC coordinates
 *   \param ofs a side variable to compute OFS objects
 *   \param PC COC matrix: z = PC*zh+V
 *   \param V COC vector: z = PC*zh+V
 *   \param zv_nc the output array to update, in NC coordinates
 *   \param graph_style if true, the special case of the QBCP is used to compute Wh.
 *          Namely: only Wh[1] and Wh[4] are complete polynomials, the rest of the components are of order one.
 **/
void ccm8_to_nc_by_tfc(const double sv_ccm8[],
                       const double t,
                       const double n,
                       const int ofts_order,
                       const int ofs_order,
                       vector<Oftsc>& Wh,
                       matrix<Ofsc>& PC,
                       vector<Ofsc>& V,
                       double zv_nc[],
                       bool graph_style)
{
    //------------------------------------------------------------------------------------
    // Inner variables (CCM, TFC, NC)
    //------------------------------------------------------------------------------------
    cdouble z_nc_comp[6];
    cdouble sv_ccm[REDUCED_NV];
    vector<Ofsc> zv_tfc(6), zv_nc_ofs(6);

    //------------------------------------------------------------------------------------
    // CCM8 to TFC
    //------------------------------------------------------------------------------------
    ccm8_to_cmm(sv_ccm8, sv_ccm);
    ccm_to_tfc(sv_ccm, ofts_order, ofs_order, Wh, zv_tfc, graph_style);

    //------------------------------------------------------------------------------------
    // TFC to NC
    //------------------------------------------------------------------------------------
    apply_coc_ofs(PC, V, zv_tfc, zv_nc_ofs);
    for(int p = 0; p < 6; p++)
    {
        z_nc_comp[p] = zv_nc_ofs[p].evaluate(n*t, ofs_order);
        zv_nc[p] = creal(z_nc_comp[p]);
    }
}


/**
 *   \brief Evaluate the configuration zv_nc = W(sv_ccm, n*t) with the use of an
 *          intermediate TFC configuration zv_tfc(t) = Wh(sv_ccm, t)
 *   \param sv_ccm an array of 4 complex double which gives the configuration to input in CCM coordinates
 *   \param t the current time
 *   \param n the pulsation of the system
 *   \param ofts_order the order of the Taylor expansions to use in the evaluation
 *   \param ofs_order the order of the Fourier expansions to use in the evaluation
 *   \param Wh the Fourier-Taylor expansion that contains the parameterization Wh(s,t), in TFC coordinates
 *   \param ofs a side variable to compute OFS objects
 *   \param PC COC matrix: z = PC*zh+V
 *   \param V COC vector: z = PC*zh+V
 *   \param zv_nc the output array to update
 *   \param graph_style if true, the special case of the QBCP is used to compute Wh.
 *          Namely: only Wh[1] and Wh[4] are complete polynomials, the rest of the components are of order one.
 **/
void ccm_to_nc_by_tfc(cdouble sv_ccm[],
                      const double t,
                      const double n,
                      const int ofts_order,
                      const int ofs_order,
                      vector<Oftsc>& Wh,
                      matrix<Ofsc>& PC,
                      vector<Ofsc>& V,
                      double zv_nc[],
                      bool graph_style)
{
    //------------------------------------------------------------------------------------
    // Inner variables (CCM, TFC, NC)
    //------------------------------------------------------------------------------------
    cdouble z_nc_comp[6];
    vector<Ofsc> zv_tfc(6), zv_nc_ofs(6);

    //------------------------------------------------------------------------------------
    // CCM to TFC
    //------------------------------------------------------------------------------------
    ccm_to_tfc(sv_ccm, ofts_order, ofs_order, Wh, zv_tfc, graph_style);

    //------------------------------------------------------------------------------------
    // TFC to NC
    //------------------------------------------------------------------------------------
    apply_coc_ofs(PC, V, zv_tfc, zv_nc_ofs);
    for(int p = 0; p < 6; p++)
    {
        z_nc_comp[p] = zv_nc_ofs[p].evaluate(n*t, ofs_order);
        zv_nc[p] = creal(z_nc_comp[p]);
    }
}

//----------------------------------------------------------------------------------------
//
//          Reduced vector field
//
//----------------------------------------------------------------------------------------
/**
 *   \brief Evaluate the reduced vector field (RVF) dot(sv_ccm8) = fv_ccm(sv_ccm8, t)
 *   \param sv_ccm8 an array of doubles which gives the current CCM state, with separated real and imag parts.
 *   \param t the current time
 *   \param n the pulsation of the system
 *   \param ofts_order the order of the Taylor expansions to use in the evaluation
 *   \param ofs_order the order of the Fourier expansions to use in the evaluation
 *   \param fh the Fourier-Taylor expansion that contains the RVF.
 *   \param ofs a side variable to compute OFS objects
 *   \param fv_ccm the RVF output array of cdouble to update, in CCM coordinates
 **/
void rvf_eval_cmm8(const double sv_ccm8[],
                   const double t,
                   const double n,
                   const int ofts_order,
                   const int ofs_order,
                   vector<Oftsc>& fh,
                   Ofsc& ofs,
                   cdouble fv_ccm[])
{
    //------------------------------------------------------------------------------------
    //CCM8 to CCM
    //------------------------------------------------------------------------------------
    cdouble s[REDUCED_NV];
    ccm8_to_cmm(sv_ccm8, s);

    //------------------------------------------------------------------------------------
    //Evaluation of fh
    //------------------------------------------------------------------------------------
    for(int p = 0; p < REDUCED_NV; p++)
    {
        fh[p].evaluate(s, ofs, ofts_order, ofs_order);
        fv_ccm[p] = ofs.evaluate(n*t, ofs_order);
    }
}

/**
 *   \brief Evaluate the reduced vector field (RVF) dot(s4) = fv_ccm8(s4, t)
 *   \param sv_ccm8 an array of 8 doubles which gives the current CCM state, with separated real and imag parts.
 *   \param t the current time
 *   \param n the pulsation of the system
 *   \param ofts_order the order of the Taylor expansions to use in the evaluation
 *   \param ofs_order the order of the Fourier expansions to use in the evaluation
 *   \param fh the Fourier-Taylor expansion that contains the RVF.
 *   \param ofs a side variable to compute OFS objects
 *   \param fv_ccm8 the RVF output array of double to update, with separated real and imag parts.
 **/
void rvf8_eval_cmm8(const double sv_ccm8[],
                    const double t,
                    const double n,
                    const int ofts_order,
                    const int ofs_order,
                    vector<Oftsc>& fh,
                    Ofsc& ofs,
                    double fv_ccm8[])
{
    //------------------------------------------------------------------------------------
    //CCM8 to RVF
    //------------------------------------------------------------------------------------
    cdouble sd[REDUCED_NV];
    rvf_eval_cmm8(sv_ccm8, t, n, ofts_order, ofs_order, fh, ofs, sd);

    //------------------------------------------------------------------------------------
    //Separation of real and imag parts (RFV to RVF8)
    //------------------------------------------------------------------------------------
    int p2 = 0;
    for(int p = 0; p < REDUCED_NV; p++)
    {
        fv_ccm8[p2++] = creal(sd[p]);
        fv_ccm8[p2++] = cimag(sd[p]);
    }
}


//----------------------------------------------------------------------------------------
//
//          Projection on the high-order center manifold
//
//----------------------------------------------------------------------------------------
/**
 *  \brief Projection of the current NC state on the central manifold in CCM coordinates
 **/
void nc_proj_cmm(const double z[], const double t, const double n, const int ofs_order,
                 matrix<Ofsc>& CQ, vector<Ofsc>& V, double omega1, double omega3,
                 cdouble sc[], int nv)
{
    //------------------------------------------------------------------------------------
    //TFC & TF coordinates
    //------------------------------------------------------------------------------------
    cdouble zv_tfc[6], zv_tf[6];

    //------------------------------------------------------------------------------------
    //z - V
    //------------------------------------------------------------------------------------
    for(int p = 0; p < 6; p++) zv_tf[p] = z[p] - V[p].evaluate(n*t, ofs_order);

    //------------------------------------------------------------------------------------
    //Update Wh = CQ*(z - V)
    //------------------------------------------------------------------------------------
    for(int k = 0; k <6; k++)
    {
        zv_tfc[k] = 0.0+0.0*I;
        for(int p = 0; p <6; p++)
        {
            zv_tfc[k] += zv_tf[p]* CQ.get_coef(k, p).evaluate(n*t, ofs_order);
        }
    }

    //------------------------------------------------------------------------------------
    //Projection on the center manifold
    //------------------------------------------------------------------------------------
    tfc_to_ccm(zv_tfc, omega1, omega3, sc, nv);
}

/**
 *  \brief Projection of the current TFC state on the central manifold in CCM coordinates.
 *         If the dimension of the reduced vector set is greater than 4, the remaining
 *         coordinates (s5 and possibly s6) are set to zero.
 **/
void tfc_to_ccm(const cdouble zv_tfc[], double omega1, double omega3,
                cdouble sc[], int nv)
{
    //------------------------------------------------------------------------------------
    //Projection on the center manifold
    //------------------------------------------------------------------------------------
    //Wh1 = i*w1*s1 => s1 = -i/w1*Wh1
    sc[0] = -1.0*I/omega1*zv_tfc[0];

    //Wh3 = i*w3*s2 => s2 = -i/w3*Wh3
    sc[1] = -1.0*I/omega3*zv_tfc[2];

    //Wh4 = -i*w1*s3 => s3 = +i/w1*Wh4
    sc[2] = +1.0*I/omega1*zv_tfc[3];

    //Wh6 = -i*w3*s4 => s4 = +i/w3*Wh6
    sc[3] = +1.0*I/omega3*zv_tfc[5];

    if(nv > 4) sc[4] = 0.0;
    if(nv > 5) sc[5] = 0.0;
}
