#ifndef COC_H_INCLUDED
#define COC_H_INCLUDED


/**
 * \file coc.h
 * \brief Implements the complete change of coordinates (coc) from the
 *        Translated-Floquet-Complexified (TFC) to the Normalized-Centered (NC)
 *        coordinates of the QBCP.
 * \author BLB
 *
 */

//Custom
#include "vf.h"
#include "env.h"
#include "timec.h"
#include "matrix.h"

#include "init.h"
#include "nfo2.h"

#include <stdio.h>
#include <iostream>



//----------------------------------------------------------------------------------------
//
//          OFTS version of the coc
//
//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------
//          Initialization
//----------------------------------------------------------------------------------------
/**
 *  \brief Retrieves the Fouriers series coefficients of the matrices P and Q, the vectors
 *         V, Xe, Xm and Xs and the scalars ILe, ILs and ILm.
 *
 *  - We recall that the change of coordinates between the NC coordinates and the TFC
 *  coordinates is of the form:
 *                              zv_nc = P(t) * C *  zv_tfc(t) + V(t)          (1)
 *  where:
 *
 *
 *  - zv_tfc is the 6-dim state vector in TFC coordinates
 *  - zv_nc is the 6-dim state vector in NC coordinates
 *  - P(t) is a 6 x 6 periodic matrix, whose coefficients are Fourier series
 *  - V(t) is a 6 x 1 periodic vector, whose coefficients are Fourier series
 *  - C is a constant 6 x 6 complex matrix.
 *
 *  In practice, we define and use the matrix PC(t) = P(t)*C.
 *  The inverse of (1) is denoted:
 *                                zv_tfc = CQ(t) * (zv_nc(t) - V(t)           (2)
 *
 *  In this routine,
 *
 *  1. First step
 *      * P, Q and V are retrieved from txt files in the folder \c fbpl.F_COC.
 *      * dot(V) is computed and stored in Vdot.
 *      * The matrices PC = P(t)*C and CQ = inv(PC) are computed and stored in PC and CQ.
 *      * dot(PC) is computed and stored in PCdot.
 *
 *  2. Second step
 *     * The vectors Xe[0:2], Xm[0:2], and Xs[0:2] contain the true
 *       (non-shifted) time-dependent positions of the primaries in the xy-plane. They are
 *       retrieved from txt files in the folder \c fbpl.F_COC.
 *     * The Fourier series ILc = 1/sqrt(Xc[0]^2 + Xc[1]^2) are retrieved from txt files,
 *       for c = e, m, s.
 *
 *
 *  Note that, in the CRTBP case, the Fourier series become constant variables, and are
 *  computed from CRTBP constants rather than retrieved from txt files.
 *
 **/
void init_coc_ofs(matrix<Ofsc>& P, matrix<Ofsc>& Q, matrix<Ofsc>& PC, matrix<Ofsc>& PCdot,
                  matrix<Ofsc>& CQ, vector<Ofsc>& Xe, vector<Ofsc>& Xm, vector<Ofsc>& Xs,
                  vector<Ofsc>& V, vector<Ofsc>& Vdot, Ofsc& ILe, Ofsc& ILm, Ofsc& ILs,
                  FBPL& fbpl);

//----------------------------------------------------------------------------------------
//          Applying
//----------------------------------------------------------------------------------------
/**
 *  \brief Apply the direct change of variables (1) at every order.
 *         The Fourier series are supposed to be in ofs format.
 **/
void apply_coc_ofts(matrix<Ofsc>& PC,
                    vector<Ofsc>& V,
                    vector<Oftsc>& zv_tfc_in,
                    vector<Oftsc>& zv_nc_out);


/**
 *  \brief Apply the change of variables (1) at order m, with or without the zero order
 *         V(t), depending on the boolean is_zero_order_shift.
 *         The Fourier series are supposed to be in ofs format.
 **/
void apply_coc_ofts(matrix<Ofsc>& PC,
                    vector<Ofsc>& V,
                    vector<Oftsc>& zv_tfc_in,
                    vector<Oftsc>& zv_nc_out,
                    int m,
                    int is_zero_order_shift);


/**
 *  \brief Apply the change of variables (1) at order m, with or without the zero order
 *         V(t), depending on the boolean is_zero_order_shift.
 *         The Fourier series are supposed to be in tfs format.
 **/
void apply_coc_tfts(matrix<Ofsc>& PC,
                    vector<Ofsc>& V,
                    vector<Oftsc>& zv_tfc_in,
                    vector<Oftsc>& zv_nc_out,
                    int m,
                    int flag);

/**
 *  \brief Apply the inverse (2) of the change of variables (1) at every order.
 *         The Fourier series are supposed to be in ofs format.
 **/
void apply_inc_coc_ofts(matrix<Ofsc>& CQ,
                        vector<Ofsc>& V,
                        vector<Oftsc>& zv_nc_in,
                        vector<Oftsc>& zv_tfc_out,
                        vector<Oftsc>& zv_tfr_tp);
/**
 *  \brief Apply the derivative of the direct change of variables (1) at every order.
 *         The Fourier series are supposed to be in ofs format.
 *
 *         The change of variables is of the form:
 *
 *                  zv_nc_dot = PCdot * zv_tfc + PC * zv_tfc_dot + Vdot
 **/
void apply_coc_dot_ofts(matrix<Ofsc>& PC,
                        matrix<Ofsc>& PCdot,
                        vector<Ofsc>& Vdot,
                        vector<Oftsc>& zv_tfc,
                        vector<Oftsc>& zv_tfc_dot,
                        vector<Oftsc>& zv_nc_dot);

/**
 *  \brief Apply the derivative of the direct change of variables (1) at order m.
 *         The Fourier series are supposed to be in ofs format.
 *
 *         The change of variables is of the form:
 *
 *                  zv_nc_dot = PCdot * zv_tfc + PC * zv_tfc_dot + Vdot
 **/
void apply_coc_dot_ofts(matrix<Ofsc>& PC,
                        matrix<Ofsc>& PCdot,
                        vector<Ofsc>& Vdot,
                        vector<Oftsc>& zv_tfc,
                        vector<Oftsc>& zv_tfc_dot,
                        vector<Oftsc>& zv_nc_dot,
                        int m);

/**
 *  \brief Apply the derivative of the inverse of the change of variables at every orders.
 *         The Fourier series are supposed to be in ofs format.
 *         The change of variables is of the form:
 *
 *                  zv_tfc_dot = CQ * (zv_nc_dot - PCdot * zv_tfc - Vdot)
 **/
void apply_inv_coc_dot_ofts(matrix<Ofsc>&  CQ,
                            matrix<Ofsc>&  PCdot,
                            vector<Ofsc>&  Vdot,
                            vector<Oftsc>& zv_tfc,
                            vector<Oftsc>& zv_nc_dot,
                            vector<Oftsc>& zv_tfc_dot,
                            vector<Oftsc>& ztp_1,
                            vector<Oftsc>& ztp_2);

/**
 *  \brief Apply the derivative of the inverse of the change of variables at order m.
 *         The Fourier series are supposed to be in ofs format.
 *         The change of variables is of the form:
 *
 *                  zv_tfc_dot = CQ * (zv_nc_dot - PCdot * zv_tfc - Vdot)
 **/
void apply_inv_coc_dot_ofts(matrix<Ofsc>& CQ,
                            matrix<Ofsc>& PCdot,
                            vector<Ofsc>& Vdot,
                            vector<Oftsc>& zv_tfc,
                            vector<Oftsc>& zv_nc_dot,
                            vector<Oftsc>& zv_tfc_dot,
                            vector<Oftsc>& ztp_1,
                            vector<Oftsc>& ztp_2,
                            int m);

/**
 *  \brief Apply the derivative of the inverse of the change of variables at order m.
 *         The Fourier series are supposed to be in tfs format.
 *         The change of variables is of the form:
 *
 *                  zv_tfc_dot = CQ * (zv_nc_dot - PCdot * zv_tfc - Vdot)
 **/
void apply_inv_coc_dot_tfts(matrix<Ofsc>&  CQ,
                            matrix<Ofsc>&  PCdot,
                            vector<Ofsc>&  Vdot,
                            vector<Oftsc>& zv_tfc,
                            vector<Oftsc>& zv_nc_dot,
                            vector<Oftsc>& zv_tfc_dot,
                            vector<Oftsc>& ztp_1,
                            vector<Oftsc>& ztp_2,
                            int m);

//----------------------------------------------------------------------------------------
//
//          Switch tfs/ofs format
//
//----------------------------------------------------------------------------------------
/**
 *  \brief From ofs to tfs format for the whole coc.
 **/
void tfs_from_ofs(matrix<Ofsc>& P,
                  matrix<Ofsc>& Q,
                  matrix<Ofsc>& PC,
                  matrix<Ofsc>& PCdot,
                  matrix<Ofsc>& CQ,
                  vector<Ofsc>& Xe,
                  vector<Ofsc>& Xm,
                  vector<Ofsc>& Xs,
                  vector<Ofsc>& V,
                  vector<Ofsc>& Vdot,
                  Ofsc& ILe,
                  Ofsc& ILm,
                  Ofsc& ILs);

//----------------------------------------------------------------------------------------
//
//          OFS version of the coc
//
//----------------------------------------------------------------------------------------
/**
 *  \brief Apply the direct change of variables (1) for zv_tfc_in/zv_nc_out in
 *         vector<Ofsc> format.
 *         The Fourier series are supposed to be in ofs format.
 **/
void apply_coc_ofs(matrix<Ofsc>& PC, vector<Ofsc>& V, vector<Ofsc>& zv_tfc_in,
                   vector<Ofsc>& zv_nc_out);

#endif // COC_H_INCLUDED
