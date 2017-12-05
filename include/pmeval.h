#ifndef PMCOC_H_INCLUDED
#define PMCOC_H_INCLUDED

#include "coc.h"
#include "matrix.h"

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
void ccm_to_rcm(const cdouble sv_ccm[], double sv_rcm[], int nv);

/**
 *  \brief from RCM to CCM coordinates
 **/
void rcm_to_ccm(const double sv_rcm[], cdouble sv_ccm[], int nv);

/**
 *  \brief from RCM to CCM coordinates, with real and imag part stored separately
 **/
void rcm_to_cmm8(const double sv_rcm[], double sv_ccm8[]);

/**
 *  \brief from CCM coordinates, with real and imag part stored separately, to RCM coordinates
 **/
void ccm8_to_rcm(const double sv_ccm8[], double sv_rcm[]);

/**
 *  \brief from CCM coordinates, with real and imag part stored separately, to CCM coordinates
 **/
void ccm8_to_cmm(const double sv_ccm8[], cdouble sv_ccm[]);

/**
 *  \brief from CCM coordinates to CCM coordinates, with real and imag part stored separately.
 **/
void ccm_to_cmm8(const cdouble sv_ccm[], double sv_ccm8[]);

//----------------------------------------------------------------------------------------
//
//          Change of coordinates for the 6-dimensional state
//
//----------------------------------------------------------------------------------------
/**
 *  \brief from TFC to TF coordinates
 **/
void tfc_to_tf(const cdouble zv_ccm[6], double zv_tf[6]);

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
                bool graph_style);

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
                bool graph_style);


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
                 bool graph_style);

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
               bool graph_style);

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
                bool graph_style);


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
                bool graph_style);


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
                 bool graph_style);


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
                      bool graph_style);

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
                       bool graph_style);


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
                      bool graph_style);

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
                   cdouble fv_ccm[]);

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
                    double fv_ccm8[]);

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
                 cdouble sc[], int nv);

/**
 *  \brief Projection of the current TFC state on the central manifold in CCM coordinates.
 *         If the dimension of the reduced vector set is greater than 4, the remaining
 *         coordinates (s5 and possibly s6) are set to zero.
 **/
void tfc_to_ccm(const cdouble zv_tfc[], double omega1, double omega3,
                cdouble sc[], int nv);
#endif // PMCOC_H_INCLUDED
