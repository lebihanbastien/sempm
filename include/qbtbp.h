#ifndef QBTBP_OFS_H_INCLUDED
#define QBTBP_OFS_H_INCLUDED

/**
 * \file qbtbp.h
 * \brief Routines of the Quasi-Bicircular Three-Body Problem (QBTBP). This file
 *        contains the routine to compute and store the QBTBP in Fourier series format.
 *        It also contains the routines to compute and store the coefficients of the
 *        Hamiltonian and the vector field of the Quasi-Bicircular Problem (QBCP).
 *
 *        Note that is also contains routines to do the same for the Bicircular Problem
 *        (BCP), with the same Fourier format.
 *
 * \author BLB
 */

#include <iostream>
#include <time.h>
#include <stdio.h>
#include <complex.h>

//GSL
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_fft_complex.h>

#include <gsl/gsl_sf_bessel.h>

//Custom
#include "ofts.h"
#include "ode.h"
#include "gslc.h"
#include "env.h"
#include "init.h"
#include "em_se_in.h"

extern "C" {
#include "nrutil.h"
}


//----------------------------------------------------------------------------------------
// Main routine: computation of the QBTBP & the QBCP
//----------------------------------------------------------------------------------------
/**
 *  \brief Main routine to compute the Quasi-Bicircular Three-Body Problem (QBTBP).
 *         The iteratives equation of the QBTBP are solved using qbtbp_ofts.
 *         The output is the Earth-Moon inner motion z(t) and Sun(Earth+Moon) outer motion
 *         Z(t) in Fourier series, up to the order OFS_ORDER.
 *         The coefficients of the vector field of the Four-Body Quasi-Bicircular Problem
 *        (QBCP) are also computed in stored if is_stored == true.
 *
 *  \param is_stored: if true, the following data are stored:
 *              - The inner and outer motion of the QBTBP are saved in data/qbtbp, using
 *                the routine qbtbp_ofs_alg_bj_cj.
 *              - The Fourier coefficients of the equations of motion of the QBCP
 *                are computed from the inner and outer motion of the QBTBP and saved in
 *                data/VF/QBCP using qbtbp_ofs_fft_alpha, and qbtbp_ofs_fft_delta. These
 *                equations of motion are computed about EML1,2,3 and SEL1,2,3.
 *
 *  \param is_test_on: if true, the following tests are performed:
 *                   - a test on a full period is performed vs the numerical integration of the QBTBP.
 *                   - a test of the results in Earth-Moon/Inertial/Sun-Earth frameworks.
 *                   - a comparison between the coefficients obtained with FFT and the
 *                     ones obtained via algebraic manipulations on Ofs objects.
 **/
void qbtbp_and_qbcp(int is_test_on, int is_stored);

//----------------------------------------------------------------------------------------
// Computing the QBTBP
//----------------------------------------------------------------------------------------
/**
 *  \brief Computes the Sun-Earth-Moon Quasi-Bicircular Three-Body Problem in Ofts format.
 *         The iteratives equation of the QBTBP are solved. See appendix A of BLB 2017.
 *         The outputs are zr_ofts and Zr_ofts, two Fourier-Taylor series (Ofts) that
 *         describe the inner and outer motions of the QBTBP. The USYS structure us_em
 *         contains the constants of the QBCP in Earth-Moon normalized units.
 */
void qbtbp_ofts(Oftsd& zr_ofts, Oftsd& Zr_ofts, USYS& us_em, CSYS& cs, int is_stored);

//----------------------------------------------------------------------------------------
// Storing the results using FFT
//----------------------------------------------------------------------------------------
/**
 *  \brief FFT and Storage in txt files of the alpha_i Fourier series. Makes use
 *         of FFT routines via qbtbp_ofs_fft_unpack.
 *         The following coefficients are computed and stored:
 *              - The alpha series (1 to 18).
 *              - The position of the primaries in EM coordinates
 *                (redundant with alpha7-12): Xe, Ye, Ze
 *              - The position of the primaries in EMNC coordinates xe, ye, ze
 */
void qbtbp_ofs_fft_alpha( Oftsd& zt,     //zt = normalized Earth-Moon motion
                          Oftsd& Zt,     //Zt = normalized Sun-(Earth+Moon) motion
                          int n_order_fourier,        //Order of the Fourier expansions
                          USYS& us_em,   //Constants in EM units
                          CSYS& cs_em);  //Coefficients in EM units

/**
 *  \brief FFT and Storage in txt files of the delta_i Fourier series. Makes use of
 *         FFT routines via qbtbp_ofs_fft_unpack.
 *         The following coefficients are computed and stored:
 *              - The delta series (1 to 18).
 *              - The position of the primaries in SE coordinates system
 *                (redundant with delta7-12): Xe, Ye, Ze
 *              - The position of the primaries in SENC coordinates system: xe, ye, ze
 */
void qbtbp_ofs_fft_delta(Oftsd& zt,    //zt = normalized Earth-Moon motion
                         Oftsd& Zt,    //Zt = normalized Sun-(Earth+Moon) motion
                         int n_order_fourier,       //Order of the Fourier expansions
                         USYS& us_em,  //Constants in EM units
                         USYS& us_sem, //Constants in SEM units
                         CSYS& cs_sem);//Coefficients in SEM units

/**
 *  \brief Computes one specific FFT and stores it txt file.
 *  \param x_gsl_0: a GSL vector of size N which contains the discrete evaluations of the
 *         function f on which to perform the FFT.
 *  \param filename: a string. The final result is stored in the file filename+"_fft.txt"
 *  \param n_order_fourier: an integer. The order of the Fourier series obtained after FFT.
 *  \param flag: a boolean. If true, the function f is supposed even (sum of cosinus).
 *         If false, f is supposed odd (sum of sinus).
 */
void qbtbp_ofs_fft_unpack(gsl_vector* x_gsl_0, string filename, int n_order_fourier, int N, int flag);

//----------------------------------------------------------------------------------------
// Storing the results algebraic manipulations
//----------------------------------------------------------------------------------------
/**
 *  \brief Computation of the Fourier form of the inner and outer motion of the QBTBP,
 *   stored in data/qbtbp/bj and cj for the double form, bjc and cjc for the complex form
 *  These computations are performed via algebraic manipulation on the Fourier series,
 *  contrary to other similar routines such as qbtbp_ofs_fft_alpha that makes use of FFT.
 *
 *  Note that the computation and storage of the bj, cj, bjc and cjc should probably be
 *  pur in a separate routine: indeed, these Fourier series are still used in other
 *  computations, while the alpha functions computed in this routine are only used in test
 *  routines.
 **/
void qbtbp_ofs_alg_bj_cj(Oftsd& zr_ofts,//zr_ofts = normalized Earth-Moon motion
                     Oftsd& Zr_ofts,//Zr_ofts = normalized Sun-(Earth+Moon) motion
                     int n_order_fourier,        //Order of the Fourier expansions
                     USYS& us_em);  //Constants in EM units


/**
 *  \brief Computation of:
 *         The first 8 alpha functions of the QBCP vector field in EM coordinates.
 *         They are saved in txt files of the form cs.F_COEF+"alpha??.txt"
 *
 *  These computations are performed via algebraic manipulation on the Fourier series,
 *  contrary to other similar routines such as qbtbp_ofs_fft_alpha that makes use of FFT.
 **/
void qbtbp_ofs_alg_alpha(Oftsd& zr_ofts,//zr_ofts = normalized Earth-Moon motion
                         Oftsd& Zr_ofts,//Zr_ofts = normalized Sun-(Earth+Moon) motion
                         Oftsd& z1,     //z1 = \bar{zr_ofts}
                         Oftsd& z2,     //z2 = \bar{zr_ofts}^(-3/2)
                         Oftsd& z3,     //z3 = zr_ofts^(-1/2)
                         Oftsd& z4,     //z4 = z2*z3 = zr_ofts^(-1/2)*\bar{zr_ofts}^(-3/2)
                         Oftsd& z5,     //z5 = \bar{zr_ofts}^(-1/2)
                         Ofsd& sigma1,  //sigma1 = exp(-itheta)
                         Ofsd& sigma2,  //sigma2 = exp(+itheta)
                         int n_order_fourier,        //Order of the Fourier expansions
                         USYS& us_em,   //Constants in EM units
                         CSYS& cs);     //Coefficients in a certain unit system


//----------------------------------------------------------------------------------------
// Reccurence for computing the QBTBP
//----------------------------------------------------------------------------------------
/**
 *  \brief One step of the recurrence scheme of the QBTBP.
 */
void qbtbp_ofts_recurrence(Oftsd& zr_ofts, //zr_ofts = normalized Earth-Moon motion
                           Oftsd& Zr_ofts, //Zr_ofts = normalized Sun-(Earth+Moon) motion
                           Oftsd& z1,      //z1 = \bar{zr_ofts}
                           Oftsd& z2,      //z2 = \bar{zr_ofts}^(-3/2)
                           Oftsd& z3,      //z3 = zr_ofts^(-1/2)
                           Oftsd& z4,      //z4 = z2*z3 = zr_ofts^(-1/2)*\bar{zr_ofts}^(-3/2)
                           Oftsd& z5,      //z5 = \bar{zr_ofts}^(-1/2)
                           Oftsd& z6,      //z6 = z3*z5 = (zr_ofts*\bar{zr_ofts})^(-1/2)
                           Oftsd& a1,      //a1 = mu*exp(itheta)*zr_ofts
                           Oftsd& a2,      //a2 = epsilon*a1
                           Oftsd& a3,      //a3 = Zr_ofts-mu*exp(itheta)*epsilon*zr_ofts
                           Oftsd& a4,      //a4 = a3^(-1/2).
                           Oftsd& a5,      //a5 = \bar{a3}
                           Oftsd& a6,      //a6 = a5^(-3/2).
                           Oftsd& a7,      //a7 = a4*a6,
                           Oftsd& b1,      //b1 = (1-mu)*exp(ithetb)*zr_ofts
                           Oftsd& b2,      //b2 = epsilon*b1
                           Oftsd& b3,      //b3 = Zr_ofts+(1-mu)*exp(ithetb)*epsilon*zr_ofts
                           Oftsd& b4,      //b4 = b3^(-1/2)
                           Oftsd& b5,      //b5 = \bar{b3}
                           Oftsd& b6,      //b6 = b5^(-3/2)
                           Oftsd& b7,      //b7 = b4*b6,
                           Oftsd& Pm,      //Pm = -ms/as^2*exp(-itheta)*b7 + ms/as^2*exp(-itheta)*a7 - z4
                           Oftsd& Qm,      //Qm = -ns^2*mu*b7 - ns^2*(1-mu)*a7
                           Oftsd& epsilon, //epsilon = Ofts with just 1.0 at order 1
                           Ofsd&  Pfm,     //Pfm = Pm(j) at order n
                           Ofsd&  Qfm,     //Qfm = Qm(j) at order n
                           Ofsd&  ufm,     //ufm = zr_ofts(j) at order n
                           Ofsd&  vfm,     //vfm = Zr_ofts(j) at order n
                           Ofsd& sigma1,   //sigma1 = exp(-itheta)
                           Ofsd& sigma2,   //sigma2 = exp(+itheta)
                           int m,          //order
                           int n_order_fourier,         //order of the Fourier expansions
                           USYS& us_em);   //EM units

//----------------------------------------------------------------------------------------
// Integrating the QBTBP
//----------------------------------------------------------------------------------------
/**
 *  \brief Derivatives of the QBTBP. To plug into GSL integrator.
 *
 * Note: the use of GSL library forces us to use double variables
 * As a consequence, in the case of z and Z, we need to use real and imaginary parts
 * as separate variables.
 */
int qbtbp_derivatives(double t, const double y[], double f[], void* params);

//----------------------------------------------------------------------------------------
// Testing the QBTBP
//----------------------------------------------------------------------------------------
/**
 *  \brief Test function to compare the analytical solution of the QBTBP to the
 *         numerical integration of the equations of motion.
 */
void qbtbp_test(double t1, Ofsc& bjc, Ofsc& cjc, OdeStruct ode_s, FBPL& fbpl);

/**
 *  \brief Test function to check the changes of coordinates between EM, Inertial and SEM
 *         coordinates, using the state vector of the Sun, the Earth, and the Moon.
 */
void qbtbp_test_in_em_sem(double t1, Ofsc& bjc, Ofsc& cjc);

/**
 *  \brief Test function used inside qbtbp_test_in_em_sem.
 */
void qbtbp_prim_test(const double xe_0_sem[],
             const double xe_0_em[],
             const double ye_0_in[],
             double t0,
             double t0c);

/**
 *  \brief Comparison of the QBTBP computed via FFT or via OFS (algebraic manipulations).
 */
void qbtbp_test_FFT_vs_OFS(Ofsc& bjc,       //zt(t)
                           Ofsc& cjc,       //Zt(t)
                           int n_order_fourier,          //order of the Fourier expansions
                           int N,           //Number of points
                           int type,        //Type of reference
                           OdeStruct ode_s, //ode structure
                           FBPL& fbpl);     //QBCP
//----------------------------------------------------------------------------------------
// Evaluating the QBTBP
//----------------------------------------------------------------------------------------
/**
 *  \brief Evaluate z(t), with \f$ z(t) = e^{it} z_r(t) \f$ in Earth-Moon units.
 */
cdouble eval_z(Ofsc& zt, double t, double n, double ni, double ai);

/**
 *  \brief Evaluate dz(t)/dt, with \f$ z(t) = e^{it} z_r(t) \f$ in Earth-Moon units.
 */
cdouble eval_zdot(Ofsc& zt, Ofsc& ztdot, double t, double n, double ni, double ai);

/**
 *  \brief Evaluate d2z(t)/dt2, with \f$ z(t) = e^{it} z_r(t) \f$ in Earth-Moon units.
 */
cdouble eval_zddot(Ofsc& zt, Ofsc& ztdot, Ofsc& ztddot, double t, double n, double ni, double ai);

//----------------------------------------------------------------------------------------
//Evaluation of parts or all the alpha/beta routines at a given time
//----------------------------------------------------------------------------------------
/**
 *  \brief Evaluation of Fourier series given as an array of coefficients,
 *         at a given time t.
 */
void eval_array_coef(double* alpha, double t, double omega, int order, double* params, int number);

/**
 *  \brief Evaluation of the time derivatives Fourier series given as an array of
 *         coefficients, at a given time t.
 */
void eval_array_coef_der(double* alpha, double t, double omega, int order, double* params, int number);

//----------------------------------------------------------------------------------------
//Evalution of individual alpha/beta
//or derivative of alpha/beta of a given type (Odd or Even)
//----------------------------------------------------------------------------------------
/**
 *  \brief Evaluate the sum \f$ \sum_{k = 0}^N coef(k) cos(k \omega t)  \f$.
 */
double eval_even_trigo(int order, double* coef, double* cR);

/**
 *  \brief Evaluate the sum \f$ \sum_{k = 0}^N - k \omega coef(k) sin(k \omega t)  \f$.
 */
double eval_even_trigo_der(double omega,  int order, double* coef, double* sR);

/**
 *  \brief Evaluate the sum \f$ \sum_{k = 0}^N coef(k) sin(k \omega t)  \f$.
 */
double eval_odd_trigo(int order, double* coef, double* sR);

/**
 *  \brief Evaluate the sum \f$ \sum_{k = 0}^N  k \omega coef(k) cos(k \omega t)  \f$.
 */
double eval_odd_trigo_der( double omega,  int order, double* coef, double* cR);

//----------------------------------------------------------------------------------------
// Computation of the BCP
//----------------------------------------------------------------------------------------
/**
 *  \brief Main routine to compute the Bicircular Problem in Fourier series format.
 *         Note that this routine computes and stores the coefficients of the BCP with the
 *         same format as the QBCP. Since the time-dependent coefficients of the BCP are
 *         only of order one in theta, this leads to very sparse Fourier series.
 *         However, the QBCP format is kept so that the BCP can be used in the same
 *         routines.
 */
void bcp();

/**
 *  \brief Computes and stores the coefficients of the vector field of the BCP from the
 *         Earth-Moon point of view.
 **/
void bcp_alpha(int n_order_fourier, USYS &us_em, CSYS &cs_em);

/**
 *  \brief Computes and stores the coefficients of the vector field of the BCP from the
 *         Sun-Earth point of view.
 **/
void bcp_delta(int n_order_fourier, USYS &us_sem, CSYS &cs_sem);


#endif // QBTBP_OFS_H_INCLUDED
