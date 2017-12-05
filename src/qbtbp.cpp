/**
 * \file qbtbp.cpp
 * \brief Routines of the Quasi-Bicircular Three-Body Problem (QBTBP) (src). This file
 *        contains the routine to compute and store the QBTBP in Fourier series format.
 *        It also contains the routines to compute and store the coefficients of the
 *        Hamiltonian and the vector field of the Quasi-Bicircular Problem (QBCP).
 *        The main routine is
 *                              void qbtbp_and_qbcp(int is_test_on, int is_stored)
 *
 *        Note that is also contains routines to do the same for the Bicircular Problem
 *        (BCP), with the same Fourier format.
 *
 * \author BLB
 */

#include "qbtbp.h"
using namespace std;


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
 *  \param is_stored: if true, the following data is stored:
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
void qbtbp_and_qbcp(int is_test_on, int is_stored)
{
    //------------------------------------------------------------------------------------
    // Splash
    //------------------------------------------------------------------------------------
    cout << "---------------------------------------------------" << endl;
    cout << "                                                   " << endl;
    cout << "       Resolution of the Sun-Earth-Moon            " << endl;
    cout << "   Quasi-Bicircular Three-Body Problem (QBTBP)     " << endl;
    cout << "    Quasi-Bicircular Four-Body Problem (QBCP)      " << endl;
    cout << "                                                   " << endl;
    cout << "---------------------------------------------------" << endl;
    cout << " qbtbp. The computation will be made up to order " << OFS_ORDER << endl;
    Config::configManager().coutmp();


    //------------------------------------------------------------------------------------
    // Initialization of environement via an FBPL structure
    //------------------------------------------------------------------------------------
    cout << " qbtbp. Initializing the environment..." << endl;
    //Init the Sun-Earth-Moon FBP structure
    FBP fbp;
    init_fbp(&fbp, Csts::SUN, Csts::EARTH, Csts::MOON);

    //Init the Sun-Earth-Moon FBPL structure, with arbitrary configuration
    FBPL fbpl;
    init_fbp_lib(&fbpl, &fbp, 1, 1, Csts::QBCP, Csts::EM, Csts::GRAPH,
              Csts::MAN_CENTER, Csts::MAN_CENTER, true, true);

    //Init of the internal/external motions
    Oftsd ofts_z(1,OFS_ORDER,2,OFS_ORDER);
    Oftsd ofts_Z(1,OFS_ORDER,2,OFS_ORDER);

    //------------------------------------------------------------------------------------
    // Computation of the QBTBP
    //------------------------------------------------------------------------------------
    cout << " qbtbp. Solving the QBTBP..." << endl;
    qbtbp_ofts(ofts_z, ofts_Z, fbpl.us_em, fbpl.cs, is_stored);
    cout << " qbtbp. Done." << endl;

    //------------------------------------------------------------------------------------
    // Storage (if desired)
    //------------------------------------------------------------------------------------
    if(is_stored)
    {
        //Storage of the inner and outer motions in double and complex form
        cout << " qbtbp. Storing the QBTBP inner and outer motions in data/qbtbp/" << endl;
        qbtbp_ofs_alg_bj_cj(ofts_z, ofts_Z, OFS_ORDER, fbpl.us_em);

        //Computation of the coeffs of the vector field for EML1
        cout << " qbtbp. Storing the vector field about EML1 in " << fbpl.cs_em_l1.F_COEF << endl;
        qbtbp_ofs_fft_alpha(ofts_z, ofts_Z, OFS_ORDER, fbpl.us_em, fbpl.cs_em_l1);

        //Computation of the coeffs of the vector field for EML2
        cout << " qbtbp. Storing the vector field about EML2 in " << fbpl.cs_em_l2.F_COEF << endl;
        qbtbp_ofs_fft_alpha(ofts_z, ofts_Z, OFS_ORDER, fbpl.us_em, fbpl.cs_em_l2);

        //Computation of the coeffs of the vector field for EML3
        cout << " qbtbp. Storing the vector field about EML3 in " << fbpl.cs_em_l3.F_COEF << endl;
        qbtbp_ofs_fft_alpha(ofts_z, ofts_Z, OFS_ORDER, fbpl.us_em, fbpl.cs_em_l3);

        //Computation of the coeffs of the vector field for SEL1
        cout << " qbtbp. Storing the vector field about SEL1 in " << fbpl.cs_sem_l1.F_COEF << endl;
        qbtbp_ofs_fft_delta(ofts_z, ofts_Z, OFS_ORDER, fbpl.us_em, fbpl.us_sem, fbpl.cs_sem_l1);

        //Computation of the coeffs of the vector field for SEL2
        cout << " qbtbp. Storing the vector field about SEL2 in " << fbpl.cs_sem_l2.F_COEF << endl;
        qbtbp_ofs_fft_delta(ofts_z, ofts_Z, OFS_ORDER, fbpl.us_em, fbpl.us_sem, fbpl.cs_sem_l2);

        //Computation of the coeffs of the vector field for SEL3
        cout << " qbtbp. Storing the vector field about SEL3 in " << fbpl.cs_sem_l3.F_COEF << endl;
        qbtbp_ofs_fft_delta(ofts_z, ofts_Z, OFS_ORDER, fbpl.us_em, fbpl.us_sem, fbpl.cs_sem_l3);
    }

    //------------------------------------------------------------------------------------
    // Test (if desired)
    //------------------------------------------------------------------------------------
    if(is_test_on)
    {
        //--------------------------------------------------------------------------------
        //Init the integration tools
        //--------------------------------------------------------------------------------
        //Stepper
        const gsl_odeiv2_step_type* T = gsl_odeiv2_step_rk8pd;
        OdeStruct ode_s;
        //Brent-Dekker root finding method
        const gsl_root_fsolver_type* T_root = gsl_root_fsolver_brent;
        //Ode solver parameters
        double param[2];
        param[0] = fbpl.us_em.mu_EM; //note that the qbtbp is computed in EM framework
        param[1] = fbpl.us_em.ms;    //note that the qbtbp is computed in EM framework
        //General structures
        init_ode_structure(&ode_s, T, T_root, 8, qbtbp_derivatives, param);

        //--------------------------------------------------------------------------------
        //Init the inner/outer motions in Ofs format (Fourier series)
        //--------------------------------------------------------------------------------
        Ofsd bj(OFS_ORDER);
        Ofsd cj(OFS_ORDER);
        Ofsc bjc(OFS_ORDER);
        Ofsc cjc(OFS_ORDER);

        //The epsilon paramater
        double eps = 1.0/fbpl.us_em.as;  //note that the qbtbp is computed in EM framework

        //From ots to ofs
        fts2fs(&bj, ofts_z, eps);
        fts2fs(&cj, ofts_Z, eps);

        //From double to complex double
        ofs_double_to_complex(bj, bjc);
        ofs_double_to_complex(cj, cjc);

        //--------------------------------------------------------------------------------
        // Splash
        //--------------------------------------------------------------------------------
        cout << "---------------------------------------------------" << endl;
        cout << "                Test of the QBTBP                  " << endl;
        cout << "---------------------------------------------------" << endl;
        cout << "Three consecutive tests are performed:             " << endl;

        //--------------------------------------------------------------------------------
        //Analytical Vs Numerical results
        //--------------------------------------------------------------------------------
        cout << " 1. The inner and outer motion z(t) and Z(t) are   " << endl;
        cout << " integrated on one period T of the QBTBP, and the  " << endl;
        cout << " consistency between the integration and the       " << endl;
        cout << " parameterization (Fourier series) is checked.     " << endl;
        qbtbp_test(2*M_PI, bjc, cjc, ode_s, fbpl);
        pressEnter(true, "Press Enter to proceed with the tests.");

        cout << " 2. The change of coordinates between Earth-Moon,     " << endl;
        cout << " Inertial, and Sun-Earth coordinates are tested,      " << endl;
        cout << " by looking at how state vectors of the primaries,    " << endl;
        cout << " which are functions of z(t) and Z(t), are transposed " << endl;
        cout << " by these COCs.                                       " << endl;
        qbtbp_test_in_em_sem(0.0*M_PI, bjc, cjc);
        pressEnter(true, "Press Enter to proceed with the tests.");

        cout << " 3. The resulting coefficients of the Earth-Moon equations of motion " << endl;
        cout << " (alpha functions) are tested. These coefficients have been computed " << endl;
        cout << " from z(t) and Z(t) using two different strategies:                  " << endl;
        cout << "  - Using algebraic operations of Fourier series (OFS in the sequel)." << endl;
        cout << "  - Using Fast Fourier Transform (FFT in the sequel).                " << endl;
        cout << " Both are compared.                                                  " << endl;
        cout << " This test is performed only if qbtbp has been called                " << endl;
        cout << " with is_stored == true...                                            " << endl;

        if(is_stored)
        {
            cout << "... which is the case. " << endl;
            qbtbp_test_FFT_vs_OFS(bjc, cjc, OFS_ORDER, 500, 1, ode_s, fbpl);
        }
        else
        {
            cout << "... which is not the case. " << endl;
        }

        pressEnter(true, "Press Enter to end the test session.");
    }
}


//----------------------------------------------------------------------------------------
// Solving the QBTBP
//----------------------------------------------------------------------------------------
/**
 *  \brief Computes the Sun-Earth-Moon Quasi-Bicircular Three-Body Problem in Ofts format.
 *         The iteratives equation of the QBTBP are solved. See appendix A of BLB 2017.
 *         The outputs are zr_ofts and Zr_ofts, two Fourier-Taylor series (Ofts) that
 *         describe the inner and outer motions of the QBTBP. The USYS structure us_em
 *         contains the constants of the QBCP in Earth-Moon normalized units.
 */
void qbtbp_ofts(Oftsd& zr_ofts, Oftsd& Zr_ofts, USYS& us_em, CSYS& cs, int is_stored)
{
    //------------------------------------------------------------------------------------
    //Parameters
    //------------------------------------------------------------------------------------
    int n_order_taylor      = zr_ofts.get_order();         //order of the Taylor expansion
    int n_var_taylor   = zr_ofts.get_nvar();            //number of variables of the Taylor expansion
    int n_order_fourier     = zr_ofts.get_coef_order();        //order of the Fourier coefficient
    int n_var_fourier  = zr_ofts.get_coef_nvar();    //number of variables of the Fourier coefficient

    //Physical param
    double n  = us_em.n;  //mean angular motion in EM units
    double as = us_em.as;           //Sun-(Earth+Moon barycenter) average distance in EM units
    double ms = us_em.ms;           //Sun mass in EM units

    //------------------------------------------------------------------------------------
    //Declaration of the variables
    //------------------------------------------------------------------------------------
    //Order 1
    Ofsd u1(n_order_fourier);
    Ofsd v1(n_order_fourier);
    //Order n of zr and Zr
    Oftsd z1(n_var_taylor, n_order_taylor, n_var_fourier, n_order_fourier);     //z1 = \bar{zr_ofts}
    Oftsd z2(n_var_taylor, n_order_taylor, n_var_fourier, n_order_fourier);     //z2 = \bar{zr_ofts}^(-3/2)
    Oftsd z3(n_var_taylor, n_order_taylor, n_var_fourier, n_order_fourier);     //z3 = zr_ofts^(-1/2)
    Oftsd z4(n_var_taylor, n_order_taylor, n_var_fourier, n_order_fourier);     //z4 = z2*z3 = zr_ofts^(-1/2)*\bar{zr_ofts}^(-3/2)
    Oftsd z5(n_var_taylor, n_order_taylor, n_var_fourier, n_order_fourier);     //z5 = \bar{zr_ofts}^(-1/2)
    Oftsd z6(n_var_taylor, n_order_taylor, n_var_fourier, n_order_fourier);     //z6 = z3*z5 = (zr_ofts*\bar{zr_ofts})^(-1/2)

    Oftsd a1(n_var_taylor, n_order_taylor, n_var_fourier, n_order_fourier);     //a1 = mu*exp(itheta)*zr_ofts
    Oftsd a2(n_var_taylor, n_order_taylor, n_var_fourier, n_order_fourier);     //a2 = epsilon*a1
    Oftsd a3(n_var_taylor, n_order_taylor, n_var_fourier, n_order_fourier);     //a3 = Zr_ofts-mu*exp(itheta)*epsilon*zr_ofts
    Oftsd a4(n_var_taylor, n_order_taylor, n_var_fourier, n_order_fourier);     //a4 = a3^(-1/2).
    Oftsd a5(n_var_taylor, n_order_taylor, n_var_fourier, n_order_fourier);     //a5 = \bar{a3}
    Oftsd a6(n_var_taylor, n_order_taylor, n_var_fourier, n_order_fourier);     //a6 = a5^(-3/2).
    Oftsd a7(n_var_taylor, n_order_taylor, n_var_fourier, n_order_fourier);     //a7 = a4*a6;

    Oftsd b1(n_var_taylor, n_order_taylor, n_var_fourier, n_order_fourier);     //b1 = (1-mu)*exp(ithetb)*zr_ofts
    Oftsd b2(n_var_taylor, n_order_taylor, n_var_fourier, n_order_fourier);     //b2 = epsilon*b1
    Oftsd b3(n_var_taylor, n_order_taylor, n_var_fourier, n_order_fourier);     //b3 = Zr_ofts+(1-mu)*exp(ithetb)*epsilon*zr_ofts
    Oftsd b4(n_var_taylor, n_order_taylor, n_var_fourier, n_order_fourier);     //b4 = b3^(-1/2)
    Oftsd b5(n_var_taylor, n_order_taylor, n_var_fourier, n_order_fourier);     //b5 = \bar{b3}
    Oftsd b6(n_var_taylor, n_order_taylor, n_var_fourier, n_order_fourier);     //b6 = b5^(-3/2)
    Oftsd b7(n_var_taylor, n_order_taylor, n_var_fourier, n_order_fourier);     //b7 = b4*b6;

    Oftsd Pm(n_var_taylor, n_order_taylor, n_var_fourier, n_order_fourier);     //Pm = -ms/as^2*exp(-itheta)*b7 + ms/as^2*exp(-itheta)*a7 - z4
    Oftsd Qm(n_var_taylor, n_order_taylor, n_var_fourier, n_order_fourier);     //Qm = -ns^2*mu*b7 - ns^2*(1-mu)*a7

    //Ofs used to solve the system at each order
    Ofsd  Pfm(n_order_fourier);
    Ofsd  Qfm(n_order_fourier);
    Ofsd  ufm(n_order_fourier);
    Ofsd  vfm(n_order_fourier);

    //------------------------------------------------------------------------------------
    //Specific Ofs and Ofts objects
    //------------------------------------------------------------------------------------
    //sigma1 = exp(-itheta)
    Ofsd sigma1(n_order_fourier);
    sigma1.set_coef(1.0, -1);
    //sigma2 = exp(+itheta)
    Ofsd sigma2(n_order_fourier);
    sigma2.set_coef(1.0, 1);
    //epsilon = Ofts with just 1.0 at order 1
    Oftsd epsilon(n_var_taylor, n_order_taylor, n_var_fourier, n_order_fourier);
    epsilon.set_coef0(1.0, 1, 0);

    //------------------------------------------------------------------------------------
    //Order 0 of zr and Zr
    //------------------------------------------------------------------------------------
    zr_ofts.set_coef0(1.0, 0, 0);
    Zr_ofts.set_coef0(1.0, 0, 0);

    //------------------------------------------------------------------------------------
    //Order 1 (optionnal)
    //------------------------------------------------------------------------------------
    u1.set_coef(-ms/(as*as)/6.0, 0);
    u1.set_coef(-ms/(as*as)*(24.0*n*n+24*n+9.0)/(64*pow(n,4.0)-16*n*n), -2);
    u1.set_coef(+ms/(as*as)*9.0/(64*pow(n,4.0)-16*n*n), 2);

    //------------------------------------------------------------------------------------
    //Order 0 of the recurrence
    //------------------------------------------------------------------------------------
    int m = 0;
    qbtbp_ofts_recurrence(zr_ofts, Zr_ofts, z1, z2, z3, z4, z5, z6,
                          a1, a2, a3, a4, a5, a6, a7,
                          b1, b2, b3, b4, b5, b6, b7,
                          Pm, Qm, epsilon, Pfm, Qfm, ufm, vfm,
                          sigma1, sigma2, m, n_order_fourier, us_em);

    //------------------------------------------------------------------------------------
    //Recurrence: m = 1 to zr_ofts.get_order()
    //------------------------------------------------------------------------------------
    for(m = 1; m <= zr_ofts.get_order(); m++)
    {
        qbtbp_ofts_recurrence(zr_ofts, Zr_ofts, z1, z2, z3, z4, z5, z6,
                              a1, a2, a3, a4, a5, a6, a7,
                              b1, b2, b3, b4, b5, b6, b7,
                              Pm, Qm, epsilon, Pfm, Qfm, ufm, vfm,
                              sigma1, sigma2, m, n_order_fourier, us_em);
    }


    //------------------------------------------------------------------------------------
    // Storage: we compute and store here the first 8 coefficients of the vector field
    // of the QBCP in EM coordinates. These series are obtained through algebraic
    // manipulation of Ofsc (Fourier series) objects. Such results are not used anymore
    // for the evaluation of the vector field. They have been replaced by the coefficients
    // obtained via FFT. However, they are kept here for consistency and for testing
    // purposes. This part could probably be avoided in the future.
    //------------------------------------------------------------------------------------
    if(is_stored)
    {
        qbtbp_ofs_alg_alpha(zr_ofts, Zr_ofts, z1, z2, z3, z4, z5, sigma1, sigma2, n_order_fourier, us_em, cs);
    }
}


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
void qbtbp_ofs_fft_alpha( Oftsd& zt,             //zt = normalized Earth-Moon motion
                          Oftsd& Zt,             //Zt = normalized Sun-(Earth+Moon) motion
                          int n_order_fourier,   //Order of the Fourier expansions
                          USYS& us_em,           //Constants in EM units
                          CSYS& cs_em)           //Coefficients in EM units
{
    //------------------------------------------------------------------------------------
    //Physical params specific to the FBPL
    //------------------------------------------------------------------------------------
    double c1    = cs_em.c1;
    double gamma = cs_em.gamma;

    //------------------------------------------------------------------------------------
    //Mass ratio
    //------------------------------------------------------------------------------------
    double mu_EM = us_em.mu_EM;

    //------------------------------------------------------------------------------------
    //Physical params in EM units
    //------------------------------------------------------------------------------------
    double ns = us_em.ns;  //Sun-(Earth+Moon) mean angular motion
    double ni = us_em.ni;  //Earth-Moon mean angular motion
    double n  = us_em.n;   //n = ni - ns
    double as = us_em.as;  //Sun-(Earth+Moon) mean distance
    double ai = us_em.ai;  //Earth-Moon mean distance
    double ms = us_em.ms;  //Sun mass

    //------------------------------------------------------------------------------------
    //FFT params
    //------------------------------------------------------------------------------------
    int N      = pow(2,12);
    double tf  = 2*M_PI/n;
    double t   = 0;
    double eps = 1.0/as;

    //------------------------------------------------------------------------------------
    //QBTBP init
    //------------------------------------------------------------------------------------
    Ofsd bj(n_order_fourier);
    Ofsd cj(n_order_fourier);
    Ofsc ztc(n_order_fourier);
    Ofsc Ztc(n_order_fourier);

    //------------------------------------------------------------------------------------
    //Building z(t) and Z(t)
    //------------------------------------------------------------------------------------
    //From ots to ofs
    fts2fs(&bj, zt, eps);
    fts2fs(&cj, Zt, eps);

    //From double to cdouble in ztc and Ztc
    ofs_double_to_complex(bj, ztc);
    ofs_double_to_complex(cj, Ztc);

    //Derivatives
    Ofsc ztcdot(ztc);
    Ofsc Ztcdot(Ztc);
    //Derivation
    ztcdot.dot(n);
    Ztcdot.dot(n);

    //Double derivatives
    Ofsc ztcddot(ztcdot);
    Ofsc Ztcddot(Ztcdot);
    //Derivation
    ztcddot.dot(n);
    Ztcddot.dot(n);

    //z(t) and Z(t)
    cdouble zi;
    cdouble Zi;
    cdouble zidot;
    cdouble ziddot;
    cdouble Ziddot;

    //------------------------------------------------------------------------------------
    // GSL objects
    //------------------------------------------------------------------------------------
    gsl_vector* gsl_alpha1  = gsl_vector_calloc(N);
    gsl_vector* gsl_alpha2  = gsl_vector_calloc(N);
    gsl_vector* gsl_alpha3  = gsl_vector_calloc(N);
    gsl_vector* gsl_alpha4  = gsl_vector_calloc(N);
    gsl_vector* gsl_alpha5  = gsl_vector_calloc(N);
    gsl_vector* gsl_alpha6  = gsl_vector_calloc(N);
    gsl_vector* gsl_alpha7  = gsl_vector_calloc(N);
    gsl_vector* gsl_alpha8  = gsl_vector_calloc(N);
    gsl_vector* gsl_alpha9  = gsl_vector_calloc(N);
    gsl_vector* gsl_alpha10 = gsl_vector_calloc(N);
    gsl_vector* gsl_alpha11 = gsl_vector_calloc(N);
    gsl_vector* gsl_alpha12 = gsl_vector_calloc(N);
    gsl_vector* gsl_alpha13 = gsl_vector_calloc(N);
    gsl_vector* gsl_alpha14 = gsl_vector_calloc(N);
    //rq: alpha15 is let equal to zero because of ERTBP vector field
    gsl_vector* gsl_alpha16 = gsl_vector_calloc(N);
    gsl_vector* gsl_alpha17 = gsl_vector_calloc(N);
    gsl_vector* gsl_alpha18 = gsl_vector_calloc(N);

    //Redundancy for the positions of the primaries
    gsl_vector* gsl_Xe = gsl_vector_calloc(N);
    gsl_vector* gsl_Ye = gsl_vector_calloc(N);
    gsl_vector* gsl_Ze = gsl_vector_calloc(N);
    gsl_vector* gsl_Xm = gsl_vector_calloc(N);
    gsl_vector* gsl_Ym = gsl_vector_calloc(N);
    gsl_vector* gsl_Zm = gsl_vector_calloc(N);
    gsl_vector* gsl_Xs = gsl_vector_calloc(N);
    gsl_vector* gsl_Ys = gsl_vector_calloc(N);
    gsl_vector* gsl_Zs = gsl_vector_calloc(N);


    //Positions of the primaries in NC coordinates
    gsl_vector* gsl_xe = gsl_vector_calloc(N);
    gsl_vector* gsl_ye = gsl_vector_calloc(N);
    gsl_vector* gsl_ze = gsl_vector_calloc(N);
    gsl_vector* gsl_xm = gsl_vector_calloc(N);
    gsl_vector* gsl_ym = gsl_vector_calloc(N);
    gsl_vector* gsl_zm = gsl_vector_calloc(N);
    gsl_vector* gsl_xs = gsl_vector_calloc(N);
    gsl_vector* gsl_ys = gsl_vector_calloc(N);
    gsl_vector* gsl_zs = gsl_vector_calloc(N);

    //------------------------------------------------------------------------------------
    //Inner double variables
    //------------------------------------------------------------------------------------
    double rv2;
    double alpha1i;
    double alpha2i;
    double alpha3i;
    double alpha4i;
    double alpha5i;
    double alpha6i;
    double alpha7i;
    double alpha8i;
    double alpha9i;
    double alpha10i;
    double alpha11i;
    double alpha12i;
    double alpha13i;
    double alpha14i;
    //rq: alpha15 is let equal to zero because of ERTBP vector field
    double alpha16i;
    double alpha17i;
    double alpha18i;

    double Ze[3], Zm[3], Zs[3];
    double ze[3], zm[3], zs[3];

    //Temp
    double alphati;

    //Derivatives
    double alpha1doti;
    double alpha2doti;
    double alpha3doti;

    //------------------------------------------------------------------------------------
    //Time loop to prepare the FFT
    //------------------------------------------------------------------------------------
    for(int i = 0 ; i < N; i++)
    {
        //--------------------------------------------------------------------------------
        //update the time
        //--------------------------------------------------------------------------------
        t = tf*i/N;

        //--------------------------------------------------------------------------------
        //z(t), zdot(t), Z(t) and Zdot(t)
        //--------------------------------------------------------------------------------
        zi     =  eval_z(ztc, t, n, ni, ai);
        Zi     =  eval_z(Ztc, t, n, ns, as);
        zidot  =  eval_zdot(ztc, ztcdot, t, n, ni, ai);
        ziddot =  eval_zddot(ztc, ztcdot, ztcddot, t, n, ni, ai);
        Ziddot =  eval_zddot(Ztc, Ztcdot, Ztcddot, t, n, ns, as);
        //rv2= 1/(z*zb)
        rv2    =  creal(1.0/(zi*conj(zi)));

        //--------------------------------------------------------------------------------
        // The alphas
        //--------------------------------------------------------------------------------
        alpha1i  = +rv2;
        alpha2i  = -rv2*creal(zidot*conj(zi));
        alpha3i  = +rv2*cimag(zidot*conj(zi));
        alpha4i  = -ms/(1.0+ms)*creal(Ziddot*conj(zi));
        alpha5i  = -ms/(1.0+ms)*cimag(Ziddot*conj(zi));
        alpha6i  = +creal(cpow(zi*conj(zi), -1.0/2+0.0*I));
        //Sun
        alpha7i  = +rv2*creal(Zi*conj(zi));
        alpha8i  = +rv2*cimag(Zi*conj(zi));
        //Earth
        alpha9i  = +mu_EM;
        alpha10i = +0.0;
        //Moon
        alpha11i = +mu_EM-1.0;
        alpha12i = +0.0;
        //Temp
        alphati = +creal(zi*conj(zi));

        //Derivatives of alpha1,2,3
        alpha1doti = creal(+0.0*I-( zidot*conj(zi)+zi*conj(zidot) ) / cpow(zi*conj(zi), 2.0+0.0*I));
        alpha2doti = creal(+0.0*I-( creal(ziddot*conj(zi) + zidot*conj(zidot))*zi*conj(zi)
                                    - (zidot*conj(zi)+zi*conj(zidot))*creal(zidot*conj(zi)) ) / cpow(zi*conj(zi), 2.0+0.0*I));
        alpha3doti = creal(+0.0*I+( cimag(ziddot*conj(zi) + zidot*conj(zidot))*zi*conj(zi)
                                    - (zidot*conj(zi)+zi*conj(zidot))*cimag(zidot*conj(zi)) ) / cpow(zi*conj(zi), 2.0+0.0*I));

        //Alpha13 and 14, for NC computation
        alpha13i = -(alphati*alphati*(alpha2doti*alpha1i - alpha1doti*alpha2i)+ alphati*(alpha2i*alpha2i+alpha3i*alpha3i))*c1 + alpha4i/gamma;
        alpha14i =   alphati*alphati*(alpha3doti*alpha1i - alpha1doti*alpha3i)*c1 + alpha5i/gamma;

        //Alpha16i
        alpha16i = alpha2doti - alpha2i*alpha1doti/alpha1i + alpha2i*alpha2i + alpha3i*alpha3i;
        //Alpha17i
        alpha17i = alpha3doti - alpha3i*alpha1doti/alpha1i;
        //Alpha18i
        alpha18i = alpha1doti/alpha1i;

        //--------------------------------------------------------------------------------
        // The primaries, again
        //--------------------------------------------------------------------------------
        //Sun, EM coordinates
        Zs[0] = alpha7i;
        Zs[1] = alpha8i;
        Zs[2] = 0.0;
        //Sun, NC coordinates
        sys_to_nc_prim(Zs, zs, c1, gamma);
        //Earth, EM coordinates
        Ze[0] = alpha9i;
        Ze[1] = alpha10i;
        Ze[2] = 0.0;
        //Earth, NC coordinates
        sys_to_nc_prim(Ze, ze, c1, gamma);
        //Moon, EM coordinates
        Zm[0] = alpha11i;
        Zm[1] = alpha12i;
        Zm[2] = 0.0;
        //Moon, NC coordinates
        sys_to_nc_prim(Zm, zm, c1, gamma);

        //--------------------------------------------------------------------------------
        //Storage in GSL objects at each time
        //--------------------------------------------------------------------------------
        gsl_vector_set(gsl_alpha1,  i, alpha1i);
        gsl_vector_set(gsl_alpha2,  i, alpha2i);
        gsl_vector_set(gsl_alpha3,  i, alpha3i);
        gsl_vector_set(gsl_alpha4,  i, alpha4i);
        gsl_vector_set(gsl_alpha5,  i, alpha5i);
        gsl_vector_set(gsl_alpha6,  i, alpha6i);
        gsl_vector_set(gsl_alpha7,  i, alpha7i);
        gsl_vector_set(gsl_alpha8,  i, alpha8i);
        gsl_vector_set(gsl_alpha9,  i, alpha9i);
        gsl_vector_set(gsl_alpha10, i, alpha10i);
        gsl_vector_set(gsl_alpha11, i, alpha11i);
        gsl_vector_set(gsl_alpha12, i, alpha12i);
        gsl_vector_set(gsl_alpha13, i, alpha13i);
        gsl_vector_set(gsl_alpha14, i, alpha14i);

        gsl_vector_set(gsl_alpha16, i, alpha16i);
        gsl_vector_set(gsl_alpha17, i, alpha17i);
        gsl_vector_set(gsl_alpha18, i, alpha18i);

        //--------------------------------------------------------------------------------
        //Primary, EM coordinates
        //--------------------------------------------------------------------------------
        gsl_vector_set(gsl_Xs,  i, Zs[0]);
        gsl_vector_set(gsl_Ys,  i, Zs[1]);
        gsl_vector_set(gsl_Zs,  i, Zs[2]);

        gsl_vector_set(gsl_Xe,  i, Ze[0]);
        gsl_vector_set(gsl_Ye,  i, Ze[1]);
        gsl_vector_set(gsl_Ze,  i, Ze[2]);

        gsl_vector_set(gsl_Xm,  i, Zm[0]);
        gsl_vector_set(gsl_Ym,  i, Zm[1]);
        gsl_vector_set(gsl_Zm,  i, Zm[2]);

        //--------------------------------------------------------------------------------
        //Primary, NC coordinates
        //--------------------------------------------------------------------------------
        gsl_vector_set(gsl_xs,  i, zs[0]);
        gsl_vector_set(gsl_ys,  i, zs[1]);
        gsl_vector_set(gsl_zs,  i, zs[2]);

        gsl_vector_set(gsl_xe,  i, ze[0]);
        gsl_vector_set(gsl_ye,  i, ze[1]);
        gsl_vector_set(gsl_ze,  i, ze[2]);

        gsl_vector_set(gsl_xm,  i, zm[0]);
        gsl_vector_set(gsl_ym,  i, zm[1]);
        gsl_vector_set(gsl_zm,  i, zm[2]);

    }

    //------------------------------------------------------------------------------------
    //Unpack the FFT and put in data file
    //------------------------------------------------------------------------------------
    qbtbp_ofs_fft_unpack(gsl_alpha1, cs_em.F_COEF+"alpha1", n_order_fourier, N, 1);
    qbtbp_ofs_fft_unpack(gsl_alpha2, cs_em.F_COEF+"alpha2", n_order_fourier, N, 0);
    qbtbp_ofs_fft_unpack(gsl_alpha3, cs_em.F_COEF+"alpha3", n_order_fourier, N, 1);
    qbtbp_ofs_fft_unpack(gsl_alpha4, cs_em.F_COEF+"alpha4", n_order_fourier, N, 1);
    qbtbp_ofs_fft_unpack(gsl_alpha5, cs_em.F_COEF+"alpha5", n_order_fourier, N, 0);
    qbtbp_ofs_fft_unpack(gsl_alpha6, cs_em.F_COEF+"alpha6", n_order_fourier, N, 1);
    //Sun
    qbtbp_ofs_fft_unpack(gsl_alpha7, cs_em.F_COEF+"alpha7", n_order_fourier, N, 1);
    qbtbp_ofs_fft_unpack(gsl_alpha8, cs_em.F_COEF+"alpha8", n_order_fourier, N, 0);
    //Earth
    qbtbp_ofs_fft_unpack(gsl_alpha9,  cs_em.F_COEF+"alpha9",  n_order_fourier, N, 1);
    qbtbp_ofs_fft_unpack(gsl_alpha10, cs_em.F_COEF+"alpha10", n_order_fourier, N, 0);
    //Moon
    qbtbp_ofs_fft_unpack(gsl_alpha11, cs_em.F_COEF+"alpha11", n_order_fourier, N, 1);
    qbtbp_ofs_fft_unpack(gsl_alpha12, cs_em.F_COEF+"alpha12", n_order_fourier, N, 0);
    //NC additional coef
    qbtbp_ofs_fft_unpack(gsl_alpha13, cs_em.F_COEF+"alpha13",  n_order_fourier, N, 1);
    qbtbp_ofs_fft_unpack(gsl_alpha14, cs_em.F_COEF+"alpha14",  n_order_fourier, N, 0);
    //VF with velocity additional coef
    qbtbp_ofs_fft_unpack(gsl_alpha16, cs_em.F_COEF+"alpha16",  n_order_fourier, N, 1);
    qbtbp_ofs_fft_unpack(gsl_alpha17, cs_em.F_COEF+"alpha17",  n_order_fourier, N, 0);
    qbtbp_ofs_fft_unpack(gsl_alpha18, cs_em.F_COEF+"alpha18",  n_order_fourier, N, 0);


    //--------------------------------------------------------------------------------
    //Primary, EM coordinates
    //Note that, at this step, the vertical motion of the primaries is undefined,
    //so we can put either Even or Odd in the qbtbp_ofs_fft_unpack option of
    //gsl_Zc without much difference
    //--------------------------------------------------------------------------------
    qbtbp_ofs_fft_unpack(gsl_Xs, cs_em.F_COEF+"Ps1", n_order_fourier, N, 1);
    qbtbp_ofs_fft_unpack(gsl_Ys, cs_em.F_COEF+"Ps2", n_order_fourier, N, 0);
    qbtbp_ofs_fft_unpack(gsl_Zs, cs_em.F_COEF+"Ps3", n_order_fourier, N, 1);

    qbtbp_ofs_fft_unpack(gsl_Xe, cs_em.F_COEF+"Pe1", n_order_fourier, N, 1);
    qbtbp_ofs_fft_unpack(gsl_Ye, cs_em.F_COEF+"Pe2", n_order_fourier, N, 0);
    qbtbp_ofs_fft_unpack(gsl_Ze, cs_em.F_COEF+"Pe3", n_order_fourier, N, 1);

    qbtbp_ofs_fft_unpack(gsl_Xm, cs_em.F_COEF+"Pm1", n_order_fourier, N, 1);
    qbtbp_ofs_fft_unpack(gsl_Ym, cs_em.F_COEF+"Pm2", n_order_fourier, N, 0);
    qbtbp_ofs_fft_unpack(gsl_Zm, cs_em.F_COEF+"Pm3", n_order_fourier, N, 1);


    //--------------------------------------------------------------------------------
    //Primary, NC coordinates
    //Note that, at this step, the vertical motion of the primaries is undefined,
    //so we can put either Even or Odd in the qbtbp_ofs_fft_unpack option of
    //gsl_zc without much difference
    //--------------------------------------------------------------------------------
    qbtbp_ofs_fft_unpack(gsl_xs, cs_em.F_COEF+"ps1", n_order_fourier, N, 1);
    qbtbp_ofs_fft_unpack(gsl_ys, cs_em.F_COEF+"ps2", n_order_fourier, N, 0);
    qbtbp_ofs_fft_unpack(gsl_zs, cs_em.F_COEF+"ps3", n_order_fourier, N, 1);

    qbtbp_ofs_fft_unpack(gsl_xe, cs_em.F_COEF+"pe1", n_order_fourier, N, 1);
    qbtbp_ofs_fft_unpack(gsl_ye, cs_em.F_COEF+"pe2", n_order_fourier, N, 0);
    qbtbp_ofs_fft_unpack(gsl_ze, cs_em.F_COEF+"pe3", n_order_fourier, N, 1);

    qbtbp_ofs_fft_unpack(gsl_xm, cs_em.F_COEF+"pm1", n_order_fourier, N, 1);
    qbtbp_ofs_fft_unpack(gsl_ym, cs_em.F_COEF+"pm2", n_order_fourier, N, 0);
    qbtbp_ofs_fft_unpack(gsl_zm, cs_em.F_COEF+"pm3", n_order_fourier, N, 1);


    //------------------------------------------------------------------------------------
    //Umemory release
    //------------------------------------------------------------------------------------
    gsl_vector_free(gsl_alpha1);
    gsl_vector_free(gsl_alpha2);
    gsl_vector_free(gsl_alpha3);
    gsl_vector_free(gsl_alpha4);
    gsl_vector_free(gsl_alpha5);
    gsl_vector_free(gsl_alpha6);
    gsl_vector_free(gsl_alpha7);
    gsl_vector_free(gsl_alpha8);
    gsl_vector_free(gsl_alpha9);
    gsl_vector_free(gsl_alpha10);
    gsl_vector_free(gsl_alpha11);
    gsl_vector_free(gsl_alpha12);
    gsl_vector_free(gsl_alpha13);
    gsl_vector_free(gsl_alpha14);
    gsl_vector_free(gsl_alpha16);
    gsl_vector_free(gsl_alpha17);
    gsl_vector_free(gsl_alpha18);
    gsl_vector_free(gsl_Xe);
    gsl_vector_free(gsl_Ye);
    gsl_vector_free(gsl_Ze);
    gsl_vector_free(gsl_Xm);
    gsl_vector_free(gsl_Ym);
    gsl_vector_free(gsl_Zm);
    gsl_vector_free(gsl_Xs);
    gsl_vector_free(gsl_Ys);
    gsl_vector_free(gsl_Zs);
    gsl_vector_free(gsl_xe);
    gsl_vector_free(gsl_ye);
    gsl_vector_free(gsl_ze);
    gsl_vector_free(gsl_xm);
    gsl_vector_free(gsl_ym);
    gsl_vector_free(gsl_zm);
    gsl_vector_free(gsl_xs);
    gsl_vector_free(gsl_ys);
    gsl_vector_free(gsl_zs);
}


/**
 *  \brief FFT and Storage in txt files of the delta_i Fourier series. Makes use of
 *         FFT routines via qbtbp_ofs_fft_unpack.
 *         The following coefficients are computed and stored:
 *              - The delta series (1 to 18).
 *              - The position of the primaries in SE coordinates system
 *                (redundant with delta7-12): Xe, Ye, Ze
 *              - The position of the primaries in SENC coordinates system: xe, ye, ze
 */
void qbtbp_ofs_fft_delta(Oftsd& zt,            //zt = normalized Earth-Moon motion
                         Oftsd& Zt,            //Zt = normalized Sun-(Earth+Moon) motion
                         int n_order_fourier,  //Order of the Fourier expansions
                         USYS& us_em,          //Constants in EM units
                         USYS& us_sem,         //Constants in SEM units
                         CSYS& cs_sem)         //Coefficients in SEM units
{
    //------------------------------------------------------------------------------------
    //Physical params in EM units
    //------------------------------------------------------------------------------------
    double as_EM  = us_em.as;     //Sun-(Earth+Moon) mean distance
    double ms_EM  = us_em.ms;     //Sun mass
    double mu_EM  = us_em.mu_EM;  //Mass ratio
    double mu_SEM = us_em.mu_SEM; //Mass ratio

    //------------------------------------------------------------------------------------
    //Physical params in SE units
    //------------------------------------------------------------------------------------
    double ns = us_sem.ns;    //Sun-(Earth+Moon) mean angular motion
    double ni = us_sem.ni;    //Earth-Moon mean angular motion
    double n  = us_sem.n;     //n = ni - ns
    double as = us_sem.as;    //Sun-(Earth+Moon) mean distance
    double ai = us_sem.ai;    //Earth-Moon mean distance

    //------------------------------------------------------------------------------------
    //Physical params specific to the FBPL
    //------------------------------------------------------------------------------------
    double c1    = cs_sem.c1;
    double gamma = cs_sem.gamma;

    //------------------------------------------------------------------------------------
    //FFT params
    //------------------------------------------------------------------------------------
    int N      = pow(2,12);
    double tf  = 2*M_PI/n;
    double t   = 0;
    double eps = 1.0/as_EM;

    //------------------------------------------------------------------------------------
    // Building z(t) and Z(t)
    //------------------------------------------------------------------------------------
    Ofsd bj(n_order_fourier);
    Ofsd cj(n_order_fourier);
    Ofsc ztc(n_order_fourier);
    Ofsc Ztc(n_order_fourier);

    //From ots to ofs
    fts2fs(&bj, zt, eps);
    fts2fs(&cj, Zt, eps);

    //From double to cdouble in ztc and Ztc
    ofs_double_to_complex(bj, ztc);
    ofs_double_to_complex(cj, Ztc);

    //Derivatives.
    Ofsc ztcdot(ztc);
    Ofsc Ztcdot(Ztc);
    //Derivation
    ztcdot.dot(n);
    Ztcdot.dot(n);

    //Double derivatives.
    Ofsc ztcddot(ztcdot);
    Ofsc Ztcddot(Ztcdot);
    //Derivation
    ztcddot.dot(n);
    Ztcddot.dot(n);

    //z(t) and Z(t) and derivatives
    cdouble zi;
    cdouble Zi;
    cdouble Zidot;
    cdouble Ziddot;

    //------------------------------------------------------------------------------------
    // GSL objects
    //------------------------------------------------------------------------------------
    gsl_vector* gsl_delta1  = gsl_vector_calloc(N);
    gsl_vector* gsl_delta2  = gsl_vector_calloc(N);
    gsl_vector* gsl_delta3  = gsl_vector_calloc(N);
    gsl_vector* gsl_delta4  = gsl_vector_calloc(N);
    gsl_vector* gsl_delta5  = gsl_vector_calloc(N);
    gsl_vector* gsl_delta6  = gsl_vector_calloc(N);
    gsl_vector* gsl_delta7  = gsl_vector_calloc(N);
    gsl_vector* gsl_delta8  = gsl_vector_calloc(N);
    gsl_vector* gsl_delta9  = gsl_vector_calloc(N);
    gsl_vector* gsl_delta10 = gsl_vector_calloc(N);
    gsl_vector* gsl_delta11 = gsl_vector_calloc(N);
    gsl_vector* gsl_delta12 = gsl_vector_calloc(N);
    gsl_vector* gsl_delta13 = gsl_vector_calloc(N);
    gsl_vector* gsl_delta14 = gsl_vector_calloc(N);
    //rq: delta15 is let equal to zero because of BCP vector field
    gsl_vector* gsl_delta16 = gsl_vector_calloc(N);
    gsl_vector* gsl_delta17 = gsl_vector_calloc(N);
    gsl_vector* gsl_delta18 = gsl_vector_calloc(N);


    //Redundancy for the positions of the primaries
    gsl_vector* gsl_Xe = gsl_vector_calloc(N);
    gsl_vector* gsl_Ye = gsl_vector_calloc(N);
    gsl_vector* gsl_Ze = gsl_vector_calloc(N);
    gsl_vector* gsl_Xm = gsl_vector_calloc(N);
    gsl_vector* gsl_Ym = gsl_vector_calloc(N);
    gsl_vector* gsl_Zm = gsl_vector_calloc(N);
    gsl_vector* gsl_Xs = gsl_vector_calloc(N);
    gsl_vector* gsl_Ys = gsl_vector_calloc(N);
    gsl_vector* gsl_Zs = gsl_vector_calloc(N);


    //Positions of the primaries in NC coordinates
    gsl_vector* gsl_xe = gsl_vector_calloc(N);
    gsl_vector* gsl_ye = gsl_vector_calloc(N);
    gsl_vector* gsl_ze = gsl_vector_calloc(N);
    gsl_vector* gsl_xm = gsl_vector_calloc(N);
    gsl_vector* gsl_ym = gsl_vector_calloc(N);
    gsl_vector* gsl_zm = gsl_vector_calloc(N);
    gsl_vector* gsl_xs = gsl_vector_calloc(N);
    gsl_vector* gsl_ys = gsl_vector_calloc(N);
    gsl_vector* gsl_zs = gsl_vector_calloc(N);


    //------------------------------------------------------------------------------------
    //Inner double variables
    //------------------------------------------------------------------------------------
    double delta1i;
    double delta2i;
    double delta3i;
    double delta4i;
    double delta5i;
    double delta6i;
    double delta7i;
    double delta8i;
    double delta9i;
    double delta10i;
    double delta11i;
    double delta12i;
    double delta13i;
    double delta14i;
    //rq: delta15 is let equal to zero because of BCP vector field
    double delta16i;
    double delta17i;
    double delta18i;

    double Ze[3], Zm[3], Zs[3];
    double ze[3], zm[3], zs[3];

    //Temp
    double deltati;

    //Derivatives
    double delta1doti;
    double delta2doti;
    double delta3doti;


    //Temp variables
    //double g1;
    //double g2;
    double Rv2;

    //Time loop to prepare the FFT
    for(int i = 0 ; i < N; i++)
    {
        //--------------------------------------------------------------------------------
        //update time
        //--------------------------------------------------------------------------------
        t = tf*i/N;

        //--------------------------------------------------------------------------------
        //z(t),  Z(t) and derivatives
        //--------------------------------------------------------------------------------
        zi     =  eval_z(ztc, t, n, ni, ai);
        Zi     =  eval_z(Ztc, t, n, ns, as);
        Zidot  =  eval_zdot(Ztc, Ztcdot, t, n, ns, as);
        Ziddot =  eval_zddot(Ztc, Ztcdot, Ztcddot, t, n, ns, as);
        //rv2= 1/(Z*Zb)
        Rv2 = creal(1.0/(Zi*conj(Zi)));

        //--------------------------------------------------------------------------------
        // The deltas
        //--------------------------------------------------------------------------------
        delta1i  = +Rv2;
        delta2i  = -Rv2*creal(Zidot * conj(Zi));
        delta3i  = +Rv2*cimag(Zidot * conj(Zi));
        delta4i  = +0.0;
        delta5i  = +0.0;
        delta6i  = sqrt(Rv2);
        //Sun
        delta7i  = +mu_SEM;
        delta8i  = +0.0;
        //Earth
        delta9i  = (mu_EM)*Rv2*creal(zi*conj(Zi)) - ms_EM/(1.0+ms_EM);
        delta10i = (mu_EM)*Rv2*cimag(zi*conj(Zi));
        //Moon
        delta11i = (mu_EM-1)*Rv2*creal(zi*conj(Zi)) - ms_EM/(1.0+ms_EM);
        delta12i = (mu_EM-1)*Rv2*cimag(zi*conj(Zi));


        //Temp
        deltati = +creal(Zi*conj(Zi));

        //Derivatives of delta1,2,3
        delta1doti = creal(+0.0*I-( Zidot*conj(Zi)+Zi*conj(Zidot) ) / cpow(Zi*conj(Zi), 2.0+0.0*I));
        delta2doti = creal(+0.0*I-( creal(Ziddot*conj(Zi) + Zidot*conj(Zidot))*Zi*conj(Zi)
                                    - (Zidot*conj(Zi)+Zi*conj(Zidot))*creal(Zidot*conj(Zi)) ) / cpow(Zi*conj(Zi), 2.0+0.0*I));
        delta3doti = creal(+0.0*I+( cimag(Ziddot*conj(Zi) + Zidot*conj(Zidot))*Zi*conj(Zi)
                                    - (Zidot*conj(Zi)+Zi*conj(Zidot))*cimag(Zidot*conj(Zi)) ) / cpow(Zi*conj(Zi), 2.0+0.0*I));

        //Delta13 and 14, for NC computation
        delta13i = -(deltati*deltati*(delta2doti*delta1i - delta1doti*delta2i)+ deltati*(delta2i*delta2i+delta3i*delta3i))*c1 + delta4i/gamma;
        delta14i =   deltati*deltati*(delta3doti*delta1i - delta1doti*delta3i)*c1 + delta5i/gamma;

        //Delta16i
        delta16i = delta2doti - delta2i*delta1doti/delta1i + delta2i*delta2i + delta3i*delta3i;
        //Delta17i
        delta17i = delta3doti - delta3i*delta1doti/delta1i;
        //Delta18i
        delta18i = delta1doti/delta1i;

        //--------------------------------------------------------------------------------
        // The primaries, again
        //--------------------------------------------------------------------------------
        //Sun, EM coordinates
        Zs[0] = delta7i;
        Zs[1] = delta8i;
        Zs[2] = 0.0;
        //Sun, NC coordinates
        sys_to_nc_prim(Zs, zs, c1, gamma);
        //Earth, EM coordinates
        Ze[0] = delta9i;
        Ze[1] = delta10i;
        Ze[2] = 0.0;
        //Earth, NC coordinates
        sys_to_nc_prim(Ze, ze, c1, gamma);
        //Moon, EM coordinates
        Zm[0] = delta11i;
        Zm[1] = delta12i;
        Zm[2] = 0.0;
        //Moon, NC coordinates
        sys_to_nc_prim(Zm, zm, c1, gamma);

        //--------------------------------------------------------------------------------
        //Storage in GSL objects at each time
        //--------------------------------------------------------------------------------
        gsl_vector_set(gsl_delta1,  i, delta1i);
        gsl_vector_set(gsl_delta2,  i, delta2i);
        gsl_vector_set(gsl_delta3,  i, delta3i);
        gsl_vector_set(gsl_delta4,  i, delta4i);
        gsl_vector_set(gsl_delta5,  i, delta5i);
        gsl_vector_set(gsl_delta6,  i, delta6i);
        gsl_vector_set(gsl_delta7,  i, delta7i);
        gsl_vector_set(gsl_delta8,  i, delta8i);
        gsl_vector_set(gsl_delta9,  i, delta9i);
        gsl_vector_set(gsl_delta10, i, delta10i);
        gsl_vector_set(gsl_delta11, i, delta11i);
        gsl_vector_set(gsl_delta12, i, delta12i);
        gsl_vector_set(gsl_delta13, i, delta13i);
        gsl_vector_set(gsl_delta14, i, delta14i);

        gsl_vector_set(gsl_delta16, i, delta16i);
        gsl_vector_set(gsl_delta17, i, delta17i);
        gsl_vector_set(gsl_delta18, i, delta18i);

        //--------------------------------------------------------------------------------
        //Primary, EM coordinates
        //--------------------------------------------------------------------------------
        gsl_vector_set(gsl_Xs,  i, Zs[0]);
        gsl_vector_set(gsl_Ys,  i, Zs[1]);
        gsl_vector_set(gsl_Zs,  i, Zs[2]);

        gsl_vector_set(gsl_Xe,  i, Ze[0]);
        gsl_vector_set(gsl_Ye,  i, Ze[1]);
        gsl_vector_set(gsl_Ze,  i, Ze[2]);

        gsl_vector_set(gsl_Xm,  i, Zm[0]);
        gsl_vector_set(gsl_Ym,  i, Zm[1]);
        gsl_vector_set(gsl_Zm,  i, Zm[2]);

        //--------------------------------------------------------------------------------
        //Primary, NC coordinates
        //--------------------------------------------------------------------------------
        gsl_vector_set(gsl_xs,  i, zs[0]);
        gsl_vector_set(gsl_ys,  i, zs[1]);
        gsl_vector_set(gsl_zs,  i, zs[2]);

        gsl_vector_set(gsl_xe,  i, ze[0]);
        gsl_vector_set(gsl_ye,  i, ze[1]);
        gsl_vector_set(gsl_ze,  i, ze[2]);

        gsl_vector_set(gsl_xm,  i, zm[0]);
        gsl_vector_set(gsl_ym,  i, zm[1]);
        gsl_vector_set(gsl_zm,  i, zm[2]);

    }

    //------------------------------------------------------------------------------------
    //Unpack the FFT and put in data file
    //------------------------------------------------------------------------------------
    qbtbp_ofs_fft_unpack(gsl_delta1, cs_sem.F_COEF+"alpha1", n_order_fourier, N, 1);
    qbtbp_ofs_fft_unpack(gsl_delta2, cs_sem.F_COEF+"alpha2", n_order_fourier, N, 0);
    qbtbp_ofs_fft_unpack(gsl_delta3, cs_sem.F_COEF+"alpha3", n_order_fourier, N, 1);
    qbtbp_ofs_fft_unpack(gsl_delta4, cs_sem.F_COEF+"alpha4", n_order_fourier, N, 1);
    qbtbp_ofs_fft_unpack(gsl_delta5, cs_sem.F_COEF+"alpha5", n_order_fourier, N, 0);
    qbtbp_ofs_fft_unpack(gsl_delta6, cs_sem.F_COEF+"alpha6", n_order_fourier, N, 1);
    //Sun
    qbtbp_ofs_fft_unpack(gsl_delta7, cs_sem.F_COEF+"alpha7", n_order_fourier, N, 1);
    qbtbp_ofs_fft_unpack(gsl_delta8, cs_sem.F_COEF+"alpha8", n_order_fourier, N, 0);
    //Earth
    qbtbp_ofs_fft_unpack(gsl_delta9,  cs_sem.F_COEF+"alpha9",  n_order_fourier, N, 1);
    qbtbp_ofs_fft_unpack(gsl_delta10, cs_sem.F_COEF+"alpha10", n_order_fourier, N, 0);
    //Moon
    qbtbp_ofs_fft_unpack(gsl_delta11, cs_sem.F_COEF+"alpha11", n_order_fourier, N, 1);
    qbtbp_ofs_fft_unpack(gsl_delta12, cs_sem.F_COEF+"alpha12", n_order_fourier, N, 0);
    //NC additional coef
    qbtbp_ofs_fft_unpack(gsl_delta13, cs_sem.F_COEF+"alpha13",  n_order_fourier, N, 1);
    qbtbp_ofs_fft_unpack(gsl_delta14, cs_sem.F_COEF+"alpha14",  n_order_fourier, N, 0);
    //VF with velocity additionnal coef
    qbtbp_ofs_fft_unpack(gsl_delta16, cs_sem.F_COEF+"alpha16",  n_order_fourier, N, 1);
    qbtbp_ofs_fft_unpack(gsl_delta17, cs_sem.F_COEF+"alpha17",  n_order_fourier, N, 0);
    qbtbp_ofs_fft_unpack(gsl_delta18, cs_sem.F_COEF+"alpha18",  n_order_fourier, N, 0);

    //--------------------------------------------------------------------------------
    //Primary, EM coordinates
    //Note that, at this step, the vertical motion of the primaries is undefined,
    //so we can put either Even or Odd in the qbtbp_ofs_fft_unpack option of gsl_Zc without much difference
    //--------------------------------------------------------------------------------
    qbtbp_ofs_fft_unpack(gsl_Xs, cs_sem.F_COEF+"Ps1", n_order_fourier, N, 1);
    qbtbp_ofs_fft_unpack(gsl_Ys, cs_sem.F_COEF+"Ps2", n_order_fourier, N, 0);
    qbtbp_ofs_fft_unpack(gsl_Zs, cs_sem.F_COEF+"Ps3", n_order_fourier, N, 1);

    qbtbp_ofs_fft_unpack(gsl_Xe, cs_sem.F_COEF+"Pe1", n_order_fourier, N, 1);
    qbtbp_ofs_fft_unpack(gsl_Ye, cs_sem.F_COEF+"Pe2", n_order_fourier, N, 0);
    qbtbp_ofs_fft_unpack(gsl_Ze, cs_sem.F_COEF+"Pe3", n_order_fourier, N, 1);

    qbtbp_ofs_fft_unpack(gsl_Xm, cs_sem.F_COEF+"Pm1", n_order_fourier, N, 1);
    qbtbp_ofs_fft_unpack(gsl_Ym, cs_sem.F_COEF+"Pm2", n_order_fourier, N, 0);
    qbtbp_ofs_fft_unpack(gsl_Zm, cs_sem.F_COEF+"Pm3", n_order_fourier, N, 1);


    //--------------------------------------------------------------------------------
    //Primary, NC coordinates
    //Note that, at this step, the vertical motion of the primaries is undefined,
    //so we can put either Even or Odd in the qbtbp_ofs_fft_unpack option of gsl_zc without much difference
    //--------------------------------------------------------------------------------
    qbtbp_ofs_fft_unpack(gsl_xs, cs_sem.F_COEF+"ps1", n_order_fourier, N, 1);
    qbtbp_ofs_fft_unpack(gsl_ys, cs_sem.F_COEF+"ps2", n_order_fourier, N, 0);
    qbtbp_ofs_fft_unpack(gsl_zs, cs_sem.F_COEF+"ps3", n_order_fourier, N, 1);

    qbtbp_ofs_fft_unpack(gsl_xe, cs_sem.F_COEF+"pe1", n_order_fourier, N, 1);
    qbtbp_ofs_fft_unpack(gsl_ye, cs_sem.F_COEF+"pe2", n_order_fourier, N, 0);
    qbtbp_ofs_fft_unpack(gsl_ze, cs_sem.F_COEF+"pe3", n_order_fourier, N, 1);

    qbtbp_ofs_fft_unpack(gsl_xm, cs_sem.F_COEF+"pm1", n_order_fourier, N, 1);
    qbtbp_ofs_fft_unpack(gsl_ym, cs_sem.F_COEF+"pm2", n_order_fourier, N, 0);
    qbtbp_ofs_fft_unpack(gsl_zm, cs_sem.F_COEF+"pm3", n_order_fourier, N, 1);


    //------------------------------------------------------------------------------------
    //Memory release
    //------------------------------------------------------------------------------------
    gsl_vector_free(gsl_delta1);
    gsl_vector_free(gsl_delta2);
    gsl_vector_free(gsl_delta3);
    gsl_vector_free(gsl_delta4);
    gsl_vector_free(gsl_delta5);
    gsl_vector_free(gsl_delta6);
    gsl_vector_free(gsl_delta7);
    gsl_vector_free(gsl_delta8);
    gsl_vector_free(gsl_delta9);
    gsl_vector_free(gsl_delta10);
    gsl_vector_free(gsl_delta11);
    gsl_vector_free(gsl_delta12);
    gsl_vector_free(gsl_delta13);
    gsl_vector_free(gsl_delta14);
    gsl_vector_free(gsl_delta16);
    gsl_vector_free(gsl_delta17);
    gsl_vector_free(gsl_delta18);
    gsl_vector_free(gsl_Xe);
    gsl_vector_free(gsl_Ye);
    gsl_vector_free(gsl_Ze);
    gsl_vector_free(gsl_Xm);
    gsl_vector_free(gsl_Ym);
    gsl_vector_free(gsl_Zm);
    gsl_vector_free(gsl_Xs);
    gsl_vector_free(gsl_Ys);
    gsl_vector_free(gsl_Zs);
    gsl_vector_free(gsl_xe);
    gsl_vector_free(gsl_ye);
    gsl_vector_free(gsl_ze);
    gsl_vector_free(gsl_xm);
    gsl_vector_free(gsl_ym);
    gsl_vector_free(gsl_zm);
    gsl_vector_free(gsl_xs);
    gsl_vector_free(gsl_ys);
    gsl_vector_free(gsl_zs);
}

/**
 *  \brief Computes one specific FFT and stores it txt file.
 *  \param x_gsl_0: a GSL vector of size N which contains the discrete evaluations of the
 *         function f on which to perform the FFT.
 *  \param filename: a string. The final result is stored in the file filename+"_fft.txt"
 *  \param n_order_fourier: an integer. The order of the Fourier series obtained after FFT.
 *  \param flag: a boolean. If true, the function f is supposed even (sum of cosinus).
 *         If false, f is supposed odd (sum of sinus).
 */
void qbtbp_ofs_fft_unpack(gsl_vector* x_gsl_0, string filename, int n_order_fourier, int array_size, int flag)
{
    //Copy of x_gsl_0 into temp vector
    gsl_vector* xGSL = gsl_vector_calloc(array_size);
    gsl_vector_memcpy(xGSL, x_gsl_0);

    gsl_vector_complex* data_complex = gsl_vector_complex_calloc(array_size);
    gsl_fft_real_wavetable* wavetable = gsl_fft_real_wavetable_alloc (array_size);
    gsl_fft_real_workspace* workspace = gsl_fft_real_workspace_alloc (array_size);
    Ofsc xFFT(n_order_fourier);
    ofstream curentStream;

    gsl_fft_real_transform (xGSL->data, 1, array_size, wavetable, workspace);
    gsl_fft_halfcomplex_unpack(xGSL->data , data_complex->data ,  xGSL->stride ,xGSL->size);
    //Order 0
    if(flag) //even case
        xFFT.set_coef(+GSL_REAL(gsl_vector_complex_get(data_complex, 0))/array_size+I*GSL_IMAG(gsl_vector_complex_get(data_complex, 0))/array_size,  0);
    else //odd case
        xFFT.set_coef(0.0,  0);

    //Order n
    for(int i = 1; i<= n_order_fourier; i++)
    {
        if(flag) //even case
        {
            //Negative frequencies
            xFFT.set_coef(+GSL_REAL(gsl_vector_complex_get(data_complex, array_size-i))/array_size, -i);
            //Positive frequencies
            xFFT.set_coef(+GSL_REAL(gsl_vector_complex_get(data_complex, i))/array_size,  i);
        }
        else //odd case
        {
            //Negative frequencies
            xFFT.set_coef(I*GSL_IMAG(gsl_vector_complex_get(data_complex, array_size-i))/array_size, -i);
            //Positive frequencies
            xFFT.set_coef(I*GSL_IMAG(gsl_vector_complex_get(data_complex, i))/array_size,  i);
        }

    }

    //Storage in txt file
    curentStream.open((filename+"_fft.txt").c_str());
    curentStream <<  xFFT << endl;
    curentStream.close();


    if(flag) //even case
    {
        //Cosinus expansion version
        curentStream.open((filename+"c_fft.txt").c_str());
        curentStream << 0 << " " << creal(xFFT.ofs_get_coef(0)) << endl;
        for(int l = 1; l<=n_order_fourier; l++) curentStream << l << " " << creal(xFFT.ofs_get_coef(-l) + xFFT.ofs_get_coef(l))  <<  endl;
        curentStream.close();
    }
    else //odd case
    {
        //Sinus expansion version
        curentStream.open((filename+"c_fft.txt").c_str());
        for(int l = 0; l<=n_order_fourier; l++)
            curentStream << l << " " <<  cimag(xFFT.ofs_get_coef(-l) - xFFT.ofs_get_coef(l)) << endl;
        curentStream.close();
    }


    gsl_vector_complex_free(data_complex);
    gsl_vector_free(xGSL);
    gsl_fft_real_wavetable_free(wavetable);
    gsl_fft_real_workspace_free(workspace);
}


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
                         USYS& us_em)   //Constants in EM units
{
    //------------------------------------------------------------------------------------
    //Physical parameters
    //------------------------------------------------------------------------------------
    double eps = 1.0/us_em.as; //length epsilon

    //------------------------------------------------------------------------------------
    // Initialization
    //------------------------------------------------------------------------------------
    Ofsd bj(n_order_fourier), cj(n_order_fourier);
    Ofsc bjc(n_order_fourier), cjc(n_order_fourier);

    //------------------------------------------------------------------------------------
    // bj and cj
    //------------------------------------------------------------------------------------
    //From ots to ofs
    fts2fs(&bj, zr_ofts, eps);
    fts2fs(&cj, Zr_ofts, eps);
    //Store in txt files
    ofstream curentStream;
    curentStream.open("data/qbtbp/bj.txt");
    curentStream << "bj: \n" << bj << endl;
    curentStream.close();
    curentStream.open("data/qbtbp/cj.txt");
    curentStream << "cj: \n" << cj << endl;
    curentStream.close();

    //------------------------------------------------------------------------------------
    // bj and cj in complex form
    //------------------------------------------------------------------------------------
    //From double to cdouble
    ofs_double_to_complex(bj, bjc);
    ofs_double_to_complex(cj, cjc);
    //Store in txt files
    curentStream.open("data/qbtbp/bjc.txt");
    curentStream <<  bjc << endl;
    curentStream.close();
    curentStream.open("data/qbtbp/cjc.txt");
    curentStream <<  cjc << endl;
    curentStream.close();
}



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
                         int n_order_fourier, //Order of the Fourier expansions
                         USYS& us_em,   //Constants in EM units
                         CSYS& cs)      //Coefficients in a certain unit system
{
    //------------------------------------------------------------------------------------
    //Physical parameters
    //------------------------------------------------------------------------------------
    double n  = us_em.n;
    double ns = us_em.ns;
    double as = us_em.as;
    double ms = us_em.ms;

    double eps = 1.0/as;                  //length epsilon
    int n_order_taylor        = zr_ofts.get_order();      //order of the Taylor expansion
    int n_var_taylor     = zr_ofts.get_nvar();         //number of variables of the Taylor expansion
    int n_var_fourier    = zr_ofts.get_coef_nvar(); //number of variables of the Fourier coefficient

    //------------------------------------------------------------------------------------
    // Initialization
    //------------------------------------------------------------------------------------
    Ofsd bj(n_order_fourier), cj(n_order_fourier);
    Ofsc bjc(n_order_fourier), cjc(n_order_fourier);

    //r^2  = 1/alpha1
    Ofsd r2(n_order_fourier);                      //Ofs double version
    Oftsd  zr2(n_var_taylor, n_order_taylor, n_var_fourier, n_order_fourier);   //Ofts version

    //The pre computations
    Ofsc d1(n_order_fourier), d2(n_order_fourier), d3(n_order_fourier), d31(n_order_fourier),
    d4(n_order_fourier), d41(n_order_fourier), d5(n_order_fourier), d51(n_order_fourier),
    d6(n_order_fourier), d61(n_order_fourier);
    Ofsc d7(n_order_fourier), d71(n_order_fourier), d8(n_order_fourier), d81(n_order_fourier),
    d9(n_order_fourier), d91(n_order_fourier), d10(n_order_fourier), d101(n_order_fourier), d11(n_order_fourier);
    Ofsc d111(n_order_fourier), d12(n_order_fourier), d121(n_order_fourier),
    d1r(n_order_fourier), d1i(n_order_fourier), d3r(n_order_fourier), d3i(n_order_fourier),
    d4r(n_order_fourier), d4i(n_order_fourier);
    Ofsc d5r(n_order_fourier), d5i(n_order_fourier), d6r(n_order_fourier),
    d6i(n_order_fourier), d7r(n_order_fourier), d7i(n_order_fourier), d8r(n_order_fourier),
    d8i(n_order_fourier), d9r(n_order_fourier);
    Ofsc d9i(n_order_fourier), d10r(n_order_fourier), d10i(n_order_fourier),
    d11r(n_order_fourier), d11i(n_order_fourier), d12r(n_order_fourier), d12i(n_order_fourier);

    //The alphas
    Ofsd alpha1(n_order_fourier), alpha6(n_order_fourier);
    Ofs<cdouble > alpha1c(n_order_fourier), alpha2c(n_order_fourier), alpha3c(n_order_fourier),
    alpha4c(n_order_fourier), alpha5c(n_order_fourier);
    Ofs<cdouble > alpha6c(n_order_fourier), alpha7c(n_order_fourier), alpha8c(n_order_fourier);

    //sigma1c & sigma2c
    Ofs<cdouble > sigma1c(n_order_fourier), sigma2c(n_order_fourier);
    ofs_double_to_complex(sigma1, sigma1c);
    ofs_double_to_complex(sigma2, sigma2c);

    //------------------------------------------------------------------------------------
    // bj and cj
    //------------------------------------------------------------------------------------
    //From ots to ofs
    fts2fs(&bj, zr_ofts, eps);
    fts2fs(&cj, Zr_ofts, eps);
    //From double to cdouble
    ofs_double_to_complex(bj, bjc);
    ofs_double_to_complex(cj, cjc);

    //------------------------------------------------------------------------------------
    // zr and Zr in complex form and derivatives. We recall that:
    //sigma1 = exp(-itheta)
    //sigma2 = exp(+itheta
    //------------------------------------------------------------------------------------
    Ofsc zr(bjc);
    Ofsc Zr(cjc);
    //dot
    Ofsc zrdot(zr);
    zrdot.dot(n);
    Ofsc Zrdot(Zr);
    Zrdot.dot(n);
    //ddot
    Ofsc zrddot(zrdot);
    zrddot.dot(n);
    Ofsc Zrddot(Zrdot);
    Zrddot.dot(n);
    //conj
    Ofsc zrc(zr);
    zrc.conjugate();
    Ofsc Zrc(Zr);
    Zrc.conjugate();

    //Order 0:
    //------------------------------------------------------------------------------------
    //d1 = zr*conj(zr)
    d1.ofs_prod(zr, zrc);
    //d2 = Z*conj(Z) = as^2*Zr*conj(Zr)
    d2.ofs_mprod(Zr, Zrc, as*as+0.0*I);
    //d31 = zr*conj(Zr)
    d31.ofs_prod(zr, Zrc);
    //d3 = z*conj(Z) = as*sigma2*zr*conj(Zr)
    d3.ofs_mprod(sigma2c, d31, as+0.0*I);
    //d41 = Zr*conj(zr)
    d41.ofs_prod(Zr, zrc);
    //d4 = Z*conj(z) = as*sigma1*Zr*conj(zr)
    d4.ofs_mprod(sigma1c, d41, as+0.0*I);
    //d1r = real(d1), d1i = imag(d1)
    ofs_real_part(d1, d1r);
    ofs_imag_part(d1, d1i);
    //d3r = real(d3), d3i = imag(d3)
    ofs_real_part(d3, d3r);
    ofs_imag_part(d3, d3i);
    //d4r = real(d4), d4i = imag(d4)
    ofs_real_part(d4, d4r);
    ofs_imag_part(d4, d4i);


    //First order:
    //------------------------------------------------------------------------------------
    //d51 = dot(zr) + i*zr
    d51.ofs_fsum(zrdot, 1.0+0.0*I, zr, I);
    //d5 = dot(z)*conj(z)
    d5.ofs_prod(d51, zrc);

    //d61 = sigma2*d51 = exp(int)*(dot(zr) + i*zr)
    d61.ofs_prod(sigma2c, d51);
    //d6 = dot(z)*conj(Z) = as*exp(int)*(dot(zr) + i*zr)*conj(Zr)
    d6.ofs_mprod(d61, Zrc, as+0.0*I);

    //d71 = dot(Zr) + i*ns*Zr
    d71.ofs_fsum(Zrdot, 1.0+0.0*I, Zr, I*ns);
    //d7 = dot(Z)*conj(Z) = as*as*(dot(Zr) + i*ns*Zr)*conj(Zr)
    d7.ofs_mprod(d71, Zrc, as*as+0.0*I);

    //d81 = sigma1*d71 = exp(-int)*(dot(Zr) + i*ns*Zr)
    d81.ofs_sprod(sigma1c, d71);
    //d8 = dot(Z)*conj(z) = as*exp(-int)*(dot(Zr) + i*ns*Zr)*conj(z)
    d8.ofs_mprod(d81, zrc, as+0.0*I);

    //d5r = real(d5), d5i = imag(d5)
    ofs_real_part(d5, d5r);
    ofs_imag_part(d5, d5i);
    //d6r = real(d6), d6i = imag(d6)
    ofs_real_part(d6, d6r);
    ofs_imag_part(d6, d6i);
    //d7r = real(d7), d7i = imag(d7)
    ofs_real_part(d7, d7r);
    ofs_imag_part(d7, d7i);
    //d8r = real(d8), d8i = imag(d8)
    ofs_real_part(d8, d8r);
    ofs_imag_part(d8, d8i);

    //Second order:
    //------------------------------------------------------------------------------------
    //d91 = ddot(zr) +2*i*dot(zr) - zr
    d91.ofs_fsum(zrddot, 1.0+0.0*I, zrdot, 2.0*I);
    d91.ofs_fsum(d91, 1.0+0.0*I, zr, -1.0+0.0*I);
    //d9 = ddot(z)*conj(z)
    d9.ofs_prod(d91, zrc);

    //d101 = sigma2*d91 = exp(int)*(ddot(zr) +2*i*dot(zr) - zr)
    d101.ofs_prod(sigma2c, d91);
    //d10 = ddot(z)*conj(Z) = as*exp(int)*(ddot(zr) +2*i*dot(zr) - zr)*conj(Zr)
    d10.ofs_mprod(d101, Zrc, as+0.0*I);

    //d111 = ddot(Zr) +2*i*ns*dot(Zr) - ns*ns*Zr
    d111.ofs_fsum(Zrddot, 1.0+0.0*I, Zrdot, 2*ns*I);
    d111.ofs_fsum(d111, 1.0+0.0*I, Zr, -ns*ns+0.0*I);
    //d11 = ddot(Z)*conj(Z) = as*as*(ddot(Zr) +2*i*ns*dot(Zr) - ns*ns*Zr)*conj(Zr)
    d11.ofs_mprod(d111, Zrc, as*as+0.0*I);

    //d121 = sigma1*d111 = exp(-int)*(ddot(Zr) +2*i*ns*dot(Zr) - ns*ns*Zr)
    d121.ofs_prod(sigma1c, d111);
    //d12 = ddot(Z)*conj(z) = as*exp(-int)*(ddot(Zr) +2*i*ns*dot(Zr) - ns*ns*Zr)*conj(zr)
    d12.ofs_mprod(d121, zrc, as+0.0*I);


    //d9r = real(d9), d9i = imag(d9)
    ofs_real_part(d9, d9r);
    ofs_imag_part(d9, d9i);

    //d10r = real(d10), d10i = imag(d10)
    ofs_real_part(d10, d10r);
    ofs_imag_part(d10, d10i);

    //d11r = real(d11), d11i = imag(d11)
    ofs_real_part(d11, d11r);
    ofs_imag_part(d11, d11i);

    //d12r = real(d12), d12i = imag(d12)
    ofs_real_part(d12, d12r);
    ofs_imag_part(d12, d12i);


    //------------------------------------------------------------------------------------
    // The alphas
    //------------------------------------------------------------------------------------
    //alpha6
    //-----------------------------
    //z4 = 1/r2 = z^(-1/2)*zbar^(-1/2)
    z4.zero();
    z4.ofts_sprod(z3, z5);
    //alpha6
    fts2fs(&alpha6, z4, eps);
    //alpha6c
    ofs_double_to_complex(alpha6, alpha6c);
    //Storage in txt file
    ofs_sst(alpha6c, cs.F_COEF+"alpha6", 1, "");
    //-----------------------------


    //alpha1
    //-----------------------------
    //alpha1 = alpha6*alpha6
    alpha1.ofs_sprod(alpha6,alpha6);
    //alpha1c
    ofs_double_to_complex(alpha1, alpha1c);
    //Storage in txt file
    ofs_sst(alpha1c, cs.F_COEF+"alpha1", 1, "");
    //-----------------------------

    //alpha2
    //-----------------------------
    //alpha2c = -real(dot(z)*conj(z))*1/r2
    alpha2c.ofs_mprod(alpha1c, d5r, -1.0+0.0*I);
    //force the zero order to zero (odd function)
    alpha2c.set_coef(0.0, 0);
    //force the odd nature
    for(int l = -n_order_fourier; l< 0; l++) alpha2c.set_coef(+0.0*I-alpha2c.ofs_get_coef(-l), l);
    //Storage in txt file
    ofs_sst(alpha2c, cs.F_COEF+"alpha2", 0, "");
    //-----------------------------

    //alpha3
    //-----------------------------
    //alpha3c = +imag(dot(z)*conj(z))*1/r2
    alpha3c.ofs_sprod(alpha1c, d5i);
    //Storage in txt file
    ofs_sst(alpha3c, cs.F_COEF+"alpha3", 1, "");
    //-----------------------------

    //alpha4
    //-----------------------------
    alpha4c+= (+0.0*I-ms/(1+ms))*d12r;
    //Storage in txt file
    ofs_sst(alpha4c, cs.F_COEF+"alpha4", 1, "");
    //-----------------------------

    //alpha5
    //-----------------------------
    alpha5c+= (+0.0*I-ms/(1+ms))*d12i;
    //Storage in txt file
    ofs_sst(alpha5c, cs.F_COEF+"alpha5", 0, "");
    //-----------------------------

    //alpha7
    //-----------------------------
    //alpha7
    alpha7c.ofs_sprod(alpha1c, d4r);
    //Storage in txt file
    ofs_sst(alpha7c, cs.F_COEF+"alpha7", 1, "");
    //-----------------------------

    //alpha8
    //-----------------------------
    alpha8c.ofs_sprod(alpha1c, d4i);
    alpha8c.set_coef(0.0, 0);
    //Storage in txt file
    ofs_sst(alpha8c, cs.F_COEF+"alpha8", 0, "");
    //-----------------------------
}


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
                           USYS& us_em)    //EM units

{
    //Physical param
    double n  = us_em.n;
    double mu = us_em.mu_EM;
    double ns = us_em.ns;
    double as = us_em.as;
    double ms = us_em.ms;

    if(m == 0)
    {
        //--------------------------------------------------------------------------------
        //Order 0 of the recurrence
        //--------------------------------------------------------------------------------
        z1.ccopy(zr_ofts, m);
        //z1 = \bar{zr_ofts}
        z1.conjugate(m);
        //z2 = \bar{zr_ofts}^(-3/2). At order 0, z1 = z2
        z2.ccopy(z1, m);
        //z3 = zr_ofts^(-1/2). At order 0, z3 = z
        z3.ccopy(zr_ofts, m);
        //z4 = z2*z3;
        z4.ofts_sprod(z2, z3, m);
        //z5 = \bar{zr_ofts}^(-1/2). At order 0, z5 = z1
        z5.ccopy(z1, m);
        //z6 = z3*z5
        z6.ofts_sprod(z3, z5, m);

        //a1 = mu*exp(itheta)*z
        a1.ofts_smult_tu(zr_ofts, sigma2, mu, m);
        //a2 = epsilon*a1
        a2.ofts_sprod(a1, epsilon, m);
        //a3+= Z
        a3.ofts_smult_u(Zr_ofts, 1.0, m);
        //a3+= -mu*exp(itheta)*epsilon*z
        a3.ofts_smult_u(a2, -1.0, m);
        //a4 = a3^(-1/2). At order 0, a4 = a3
        a4.ccopy(a3, m);

        //a5 = \bar{a3}
        a5.ccopy(a3, m);
        a5.conjugate(m);
        //a6 = a5^(-3/2). At order 0, a6 = a5
        a6.ccopy(a5, m);
        //a7 = a4*a6;
        a7.ofts_sprod(a4,a6, m);


        //b1 = (1-mu)*exp(itheta)*z
        b1.ofts_smult_tu(zr_ofts, sigma2, (1.0-mu), m);
        //b2 = epsilon*b1
        b2.ofts_sprod(b1, epsilon, m);
        //b3+= Z
        b3.ofts_smult_u(Zr_ofts, 1.0, m);
        //b3+= (1-mu)*exp(ithetb)*epsilon*z
        b3.ofts_smult_u(b2, 1.0, m);
        //b4 = b3^(-1/2). At order 0, b4 = b3
        b4.ccopy(b3, m);
        //b5 = \bar{b3}
        b5.ccopy(b3, m);
        b5.conjugate(m);
        //b6 = b5^(-3/2). At order 0, b6 = b5
        b6.ccopy(b5, m);
        //b7 = b4*b6;
        b7.ofts_sprod(b4,b6, m);


        //Pm += -ms/as^2*exp(-itheta)*b7
        Pm.ofts_smult_tu(b7, sigma1, -ms/(as*as),m);
        //Pm +=  ms/as^2*exp(-itheta)*a7
        Pm.ofts_smult_tu(a7, sigma1,  ms/(as*as),m);
        //Pm +=  -z4
        Pm.ofts_smult_u(z4, -1.0, m);


        //Qm += -ns^2*mu*b7
        Qm.ofts_smult_u(b7, -ns*ns*mu, m);
        //Qm += -ns^2*(1-mu)*a7
        Qm.ofts_smult_u(a7, -ns*ns*(1-mu), m);
    }
    else
    {
        //--------------------------------------------------------------------------------
        //Order m of the recurrence
        //--------------------------------------------------------------------------------
        //z1 = zr_ofts at order m-1
        z1.ccopy(zr_ofts, m-1);
        //z1 = \bar{z1} at order m-1
        z1.conjugate(m-1);

        //z2 = \bar{z}^^(-3/2)
        //CAN BE USED BECAUSE z1[0] = 1.0 !!
        if(m>1) z2.ofts_pows(z1, -3.0/2, m-1);
        z2.ofts_pows(z1, -3.0/2, m);

        //z5 = \bar{z}^^(-1/2)
        //CAN BE USED BECAUSE z1[0] = 1.0 !!
        if(m>1) z5.ofts_pows(z1, -1.0/2, m-1);
        z5.ofts_pows(z1, -1.0/2, m);

        //z3 = z^(-1/2)
        //CAN BE USED BECAUSE z[0] = 1.0 !!
        if(m>1) z3.ofts_pows(zr_ofts, -1.0/2, m-1);
        z3.ofts_pows(zr_ofts, -1.0/2, m);

        //z4 = z2*z3;
        z4.ofts_sprod(z2, z3, m);

        //z6 = z3*z5
        z6.ofts_sprod(z3, z5, m);

        //a1 = mu*exp(itheta)*z at order m-1
        if(m>1) a1.ofts_smult_tu(zr_ofts, sigma2, mu, m-1);
        //a2 = epsilon*a1 at order m
        a2.ofts_sprod(a1, epsilon, m);
        //a3+= Z at order m-1
        if(m>1) a3.ofts_smult_u(Zr_ofts, 1.0, m-1);
        //a3+= -mu*exp(itheta)*epsilon*z
        a3.ofts_smult_u(a2, -1.0, m);

        //a4 = a3^(-1/2)
        //CAN BE USED BECAUSE a3[0] = 1.0 !!
        if(m>1) a4.ofts_pows(a3, -1.0/2, m-1);
        a4.ofts_pows(a3, -1.0/2, m);

        //a5 = \bar{a3}
        if(m>1)
        {
            a5.ccopy(a3, m-1);
            a5.conjugate(m-1);
        }
        a5.ccopy(a3, m);
        a5.conjugate(m);

        //a6 = a5^(-3/2)
        //CAN BE USED BECAUSE a5[0] = 1.0 !!
        if(m>1) a6.ofts_pows(a5, -3.0/2, m-1);
        a6.ofts_pows(a5, -3.0/2, m);

        //a7 = a4*a6;
        a7.ofts_sprod(a4,a6, m);

        //b1 = (1-mu)*exp(itheta)*z
        if(m>1) b1.ofts_smult_tu(zr_ofts, sigma2, (1.0-mu), m-1);
        //b2 = epsilon*b1
        b2.ofts_sprod(b1, epsilon, m);
        //b3+= Z
        if(m>1) b3.ofts_smult_u(Zr_ofts, 1.0, m-1);
        //b3+= (1-mu)*exp(itheta)*epsilon*z
        b3.ofts_smult_u(b2, 1.0, m);

        //b4 = b3^(-1/2)
        //CAN BE USED BECAUSE b3[0] = 1.0 !!
        if(m>1) b4.ofts_pows(b3, -1.0/2, m-1);
        b4.ofts_pows(b3, -1.0/2, m);

        //b5 = \bar{b3}
        if(m>1)
        {
            b5.ccopy(b3, m-1);
            b5.conjugate(m-1);
        }
        b5.ccopy(b3, m);
        b5.conjugate(m);

        //b6 = b5^(-3/2)
        //CAN BE USED BECAUSE b5[0] = 1.0 !!
        if(m>1) b6.ofts_pows(b5, -3.0/2, m-1);
        b6.ofts_pows(b5, -3.0/2, m);

        //b7 = b4*b6;
        b7.ofts_sprod(b4,b6, m);

        //Pm += -ms/as^2*exp(-itheta)*b7
        Pm.ofts_smult_tu(b7, sigma1, -ms/(as*as),m);
        //Pm +=  ms/as^2*exp(-itheta)*a7
        Pm.ofts_smult_tu(a7, sigma1,  ms/(as*as),m);
        //Pm +=  -z4
        Pm.ofts_smult_u(z4, -1.0, m);


        //Qm += -ns^2*mu*b7
        Qm.ofts_smult_u(b7, -ns*ns*mu, m);
        //Qm += -ns^2*(1-mu)*a7
        Qm.ofts_smult_u(a7, -ns*ns*(1-mu), m);

        //-------------------------
        // Solving equations
        //-------------------------
        Pfm.ccopy(Pm.get_term(m)->get_coef(0));
        Qfm.ccopy(Qm.get_term(m)->get_coef(0));


        //Order 0
        ufm.set_coef(-1.0/3*Pfm.ofs_get_coef(0), 0);          //u0 = -1/3*p0
        vfm.set_coef(-1.0/(3*ns*ns)*Qfm.ofs_get_coef(0), 0);  //v0 = -1/(3*ns^2)*p0

        //Order 1 to n_order_fourier
        //solving a 2*2 system
        double k1, k2, k3, l1, l2, l3, uj, vj;
        for(int j = 1; j<= n_order_fourier; j++)
        {
            k1 = -pow(j*n,2.0) + 2*j*n - 3.0/2;
            k2 = -pow(j*n,2.0) - 2*j*n - 3.0/2;
            k3 = -3.0/2;

            l1 = -pow(j*n,2.0) + 2*j*n*ns - 3.0/2*ns*ns;
            l2 = -pow(j*n,2.0) - 2*j*n*ns - 3.0/2*ns*ns;
            l3 = -3.0/2*ns*ns;

            //u-j
            uj = ( k2*Pfm.ofs_get_coef(-j) - k3*Pfm.ofs_get_coef(j))/(k1*k2-k3*k3);
            ufm.set_coef(uj, -j);
            //uj
            uj = (-k3*Pfm.ofs_get_coef(-j) + k1*Pfm.ofs_get_coef(j))/(k1*k2-k3*k3);
            ufm.set_coef(uj, j);
            //v-j
            vj = ( l2*Qfm.ofs_get_coef(-j) - l3*Qfm.ofs_get_coef(j))/(l1*l2-l3*l3);
            vfm.set_coef(vj, -j);
            //vj
            vj = (-l3*Qfm.ofs_get_coef(-j) + l1*Qfm.ofs_get_coef(j))/(l1*l2-l3*l3);
            vfm.set_coef(vj, j);

        }

        //Update z and Z
        zr_ofts.get_coef(m,0)->ccopy(ufm);
        Zr_ofts.get_coef(m,0)->ccopy(vfm);

    }

    if(m == zr_ofts.get_order())
    {
        //Last update to account for final order in all series
        //---------------------------------
        //z1 = zr_ofts at order m-1
        z1.ccopy(zr_ofts, m);
        //z1 = \bar{z1} at order m-1
        z1.conjugate(m);
        //z2 = \bar{z}^^(-3/2)
        if(m>1) z2.ofts_pows(z1, -3.0/2, m);
        //z5 = \bar{z}^^(-1/2)
        if(m>1) z5.ofts_pows(z1, -1.0/2, m);
        //z3 = z^(-1/2)
        if(m>1) z3.ofts_pows(zr_ofts, -1.0/2, m);
    }

}

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
int qbtbp_derivatives(double t, const double y[], double f[], void* params)
{
    //------------------------------------------------------------------------------------
    //Parameters for the qbtbp
    //------------------------------------------------------------------------------------
    double mu = * (double* ) params;
    double ms = * ((double* ) (params)+1);

    //------------------------------------------------------------------------------------
    //Reconstruction of z and Z
    //------------------------------------------------------------------------------------
    cdouble z = y[0] + I*y[1];
    cdouble Z = y[2] + I*y[3];

    cdouble temp1 = Z-mu*z;
    temp1 = temp1/pow(cabs(temp1), 3.0);

    cdouble temp2 = Z+(1-mu)*z;
    temp2 = temp2/pow(cabs(temp2), 3.0);

    cdouble zdd = +0.0*I-z/pow(cabs(z), 3.0) + ms*(temp1-temp2);
    cdouble Zdd = -(1+ms)*(mu*temp2 + (1-mu)*temp1);

    //------------------------------------------------------------------------------------
    //Phase space derivatives
    //------------------------------------------------------------------------------------
    f[0] = y[4];
    f[1] = y[5];
    f[2] = y[6];
    f[3] = y[7];
    f[4] = creal(zdd);
    f[5] = cimag(zdd);
    f[6] = creal(Zdd);
    f[7] = cimag(Zdd);

    return GSL_SUCCESS;
}

//----------------------------------------------------------------------------------------
// Testing the QBTBP
//----------------------------------------------------------------------------------------
/**
 *  \brief Test function to compare the analytical solution of the QBTBP to the
 *         numerical integration of the equations of motion.
 */
void qbtbp_test(double t1, Ofsc& bjc, Ofsc& cjc, OdeStruct ode_s, FBPL& fbpl)
{
    //Initialization
    double as = fbpl.us_em.as;
    double n  = fbpl.us_em.n;
    double ns = fbpl.us_em.ns;
    int n_order_fourier    = bjc.get_order();

    //z(0) and Z(0)
    cdouble z0 = bjc.evaluate(0.0);
    cdouble Z0 = as*cjc.evaluate(0.0);

    //zdot(0) and Zdot(0)
    Ofsc zdot(n_order_fourier);
    for(int l = -n_order_fourier; l<=n_order_fourier; l++) zdot.set_coef(I*(1.0+l*n)*bjc.ofs_get_coef(l), l);
    cdouble zdot0 = zdot.evaluate(0.0);

    Ofsc Zdot(n_order_fourier);
    for(int l = -n_order_fourier; l<=n_order_fourier; l++) Zdot.set_coef(I*(ns+l*n)*cjc.ofs_get_coef(l), l);
    cdouble Zdot0 = as*Zdot.evaluate(0.0);

    //Initial conditions
    double yv[8];
    yv[0] = creal(z0);
    yv[1] = cimag(z0);
    yv[2] = creal(Z0);
    yv[3] = cimag(Z0);
    yv[4] = creal(zdot0);
    yv[5] = cimag(zdot0);
    yv[6] = creal(Zdot0);
    yv[7] = cimag(Zdot0);

    cout << "-------------------------------------------" << endl;
    cout << "Initial z and Z, at t0 = " << 0.0 << endl;
    cout << "z(t0) = " << creal(z0) << " + " << cimag(z0) << "i" <<  endl;
    cout << "Z(t0) = " << creal(Z0) << " + " << cimag(Z0) << "i" <<  endl;
    cout << "-------------------------------------------" << endl;

    //Loop
    double h = Config::configManager().G_PREC_HSTART(); //first guess for stepper
    int status;
    double t = 0.0;
    while (t < t1)
    {
        //Evolve one step
        status = gsl_odeiv2_evolve_apply (ode_s.e, ode_s.c, ode_s.s, &ode_s.sys, &t, t1, &h, yv);
        //Break if evolution has gone wrong
        if (status != GSL_SUCCESS) break;
    }

    //Final state
    cout << "Integrated z and Z, at t1 = " << t << endl;
    cout << "z(t1) = " << yv[0] << " + " << yv[1] << "i" << endl;
    cout << "Z(t1) = " << yv[2] << " + " << yv[3] << "i" << endl;
    cout << "-------------------------------------------" << endl;

    //Analytical final state
    cdouble zfinal = (cos(t1)+I*sin(t1))*bjc.evaluate(n*t1);
    cdouble Zfinal = as*(cos(ns*t1)+I*sin(ns*t1))*cjc.evaluate(n*t1);

    cdouble zdotfinal = (cos(t1)+I*sin(t1))*zdot.evaluate(n*t1);
    cdouble Zdotfinal = as*(cos(ns*t1)+I*sin(ns*t1))*Zdot.evaluate(n*t1);


    cout << "Semi-analytical z and Z, at t1 = " << t << endl;
    cout << "z(t1) = " << creal(zfinal) << " + " << cimag(zfinal) << "i" <<  endl;
    cout << "Z(t1) = " << creal(Zfinal) << " + " << cimag(Zfinal) << "i" <<  endl;
    cout << "-------------------------------------------" << endl;


    cout << "Absolute delta between analytical and numerical results: " << endl;
    cout << "dz(t1)    = " << cabs(zfinal - yv[0]-I*yv[1])/*/cabs(zfinal)*100 << " %" */<< endl;
    cout << "dzdot(t1) = " << cabs(zdotfinal - yv[4]-I*yv[5])/*/cabs(zdotfinal)*100 << " %" */<< endl;
    cout << "dZ(t1)    = " << cabs(Zfinal - yv[2]-I*yv[3])/*/cabs(zfinal)*100 << " %" */<< endl;
    cout << "dZdot(t1) = " << cabs(Zdotfinal - yv[6]-I*yv[7])/*/cabs(Zdotfinal)*100 << " %" */<< endl;
    cout << "-------------------------------------------" << endl;
}

/**
 *  \brief Test function to check the changes of coordinates between EM, Inertial and SEM
 *         coordinates, using the state vector of the Sun, the Earth, and the Moon.
 */
void qbtbp_test_in_em_sem(double t1, Ofsc& bjc, Ofsc& cjc)
{
    //Initialization
    USYS us_em  = SEML.us_em;
    USYS us_sem = SEML.us_sem;
    int n_order_fourier      = SEML.n_order_fourier;
    //Derivatives of z and Z
    Ofsc zdot(bjc);
    Ofsc Zdot(cjc);
    zdot.dot(us_em.n);
    Zdot.dot(us_em.n);

    cout << "-------------------------------------------" << endl;
    cout << "           Test in SEM units               " << endl;
    cout << "-------------------------------------------" << endl;
    double t0 = t1;             //new t0 in EM units
    double t0c =t0*us_em.ns;    //new t0 in SEM units

    //z(t0) and Z(t0)
    cdouble z0 = eval_z(bjc, t0, us_em.n, us_em.ni, us_em.ai);
    cdouble Z0 = eval_z(cjc, t0, us_em.n, us_em.ns, us_em.as);
    cdouble zdot0 = eval_zdot(bjc, zdot, t0, us_em.n, us_em.ni, us_em.ai);
    cdouble Zdot0 = eval_zdot(cjc, Zdot, t0, us_em.n, us_em.ns, us_em.as);

    //------------------------------------------------------------------------------------
    //Evaluate the deltas @t = t0c
    //------------------------------------------------------------------------------------
    double delta[15];
    eval_array_coef(delta, t0c, us_sem.n, n_order_fourier, SEML.cs_sem.coeffs, 15);
    double deltad[15];
    eval_array_coef_der(deltad, t0c, us_sem.n, n_order_fourier, SEML.cs_sem.coeffs, 15);

    //------------------------------------------------------------------------------------
    //Evaluate the alphas @t = t0
    //------------------------------------------------------------------------------------
    double alpha[15];
    eval_array_coef(alpha, t0, us_em.n, n_order_fourier, SEML.cs_em.coeffs, 15);
    double alphad[15];
    eval_array_coef_der(alphad, t0, us_em.n, n_order_fourier, SEML.cs_em.coeffs, 15);

    cout << "State of the Earth  " << endl;
    //Earth position & momenta in SEM ref and SEM units
    double xe_0_sem[6];
    xe_0_sem[0] = delta[8];
    xe_0_sem[1] = delta[9];
    xe_0_sem[2] = 0.0;
    xe_0_sem[3] = deltad[8];
    xe_0_sem[4] = deltad[9];
    xe_0_sem[5] = 0.0;
    //Earth position & momenta in EM ref and EM units
    double xe_0_em[6];
    xe_0_em[0] = us_em.mu_EM;
    xe_0_em[1] = 0.0;
    xe_0_em[2] = 0.0;
    xe_0_em[3] = 0.0;
    xe_0_em[4] = 0.0;
    xe_0_em[5] = 0.0;
    //Earth position & momenta in IN ref and EM units
    double ye_0_in[6];
    ye_0_in[0] = -us_em.ms/(1.0+us_em.ms)*creal(Z0) + us_em.mu_EM*creal(z0);
    ye_0_in[1] = -us_em.ms/(1.0+us_em.ms)*cimag(Z0) + us_em.mu_EM*cimag(z0);
    ye_0_in[2] = 0.0;
    ye_0_in[3] = -us_em.ms/(1.0+us_em.ms)*creal(Zdot0) + us_em.mu_EM*creal(zdot0);
    ye_0_in[4] = -us_em.ms/(1.0+us_em.ms)*cimag(Zdot0) + us_em.mu_EM*cimag(zdot0);
    ye_0_in[5] = 0.0;
    //Comparison
    qbtbp_prim_test(xe_0_sem, xe_0_em, ye_0_in, t0, t0c);

    cout << "-------------------------------------------" << endl;
    cout << "State of the Moon  " << endl;
    //Moon position & velocity in SEM ref and SEM units
    double xm0_SEM[6];
    xm0_SEM[0] = delta[10];
    xm0_SEM[1] = delta[11];
    xm0_SEM[2] = 0.0;
    xm0_SEM[3] = deltad[10];
    xm0_SEM[4] = deltad[11];
    xm0_SEM[5] = 0.0;
    //Moon position & velocity in EM ref and EM units
    double xm0_EM[6];
    xm0_EM[0] = us_em.mu_EM-1;
    xm0_EM[1] = 0.0;
    xm0_EM[2] = 0.0;
    xm0_EM[3] = 0.0;
    xm0_EM[4] = 0.0;
    xm0_EM[5] = 0.0;
    //Moon position & momenta in IN ref and EM units
    double ym0_IN[6];
    ym0_IN[0] = -us_em.ms/(1.0+us_em.ms)*creal(Z0) - (1-us_em.mu_EM)*creal(z0);
    ym0_IN[1] = -us_em.ms/(1.0+us_em.ms)*cimag(Z0) - (1-us_em.mu_EM)*cimag(z0);
    ym0_IN[2] = 0.0;
    ym0_IN[3] = -us_em.ms/(1.0+us_em.ms)*creal(Zdot0) - (1-us_em.mu_EM)*creal(zdot0);
    ym0_IN[4] = -us_em.ms/(1.0+us_em.ms)*cimag(Zdot0) - (1-us_em.mu_EM)*cimag(zdot0);
    ym0_IN[5] = 0.0;
    //Comparison
    qbtbp_prim_test(xm0_SEM, xm0_EM, ym0_IN, t0, t0c);


    cout << "-------------------------------------------" << endl;
    cout << "State of the Sun  " << endl;
    //Sun position & velocity in SEM ref and SEM units
    double xs0_SEM[6];
    xs0_SEM[0] = us_sem.mu_SEM;
    xs0_SEM[1] = 0.0;
    xs0_SEM[2] = 0.0;
    xs0_SEM[3] = 0.0;
    xs0_SEM[4] = 0.0;
    xs0_SEM[5] = 0.0;
    //Sun position & velocity in EM ref and EM units
    double xs0_EM[6];
    xs0_EM[0] = alpha[6];
    xs0_EM[1] = alpha[7];
    xs0_EM[2] = 0.0;
    xs0_EM[3] = alphad[6];
    xs0_EM[4] = alphad[7];
    xs0_EM[5] = 0.0;
    //Sun position & momenta in IN ref and EM units
    double ys0_IN[6];
    ys0_IN[0] = 1.0/(1.0+us_em.ms)*creal(Z0);
    ys0_IN[1] = 1.0/(1.0+us_em.ms)*cimag(Z0);
    ys0_IN[2] = 0.0;
    ys0_IN[3] = 1.0/(1.0+us_em.ms)*creal(Zdot0);
    ys0_IN[4] = 1.0/(1.0+us_em.ms)*cimag(Zdot0);
    ys0_IN[5] = 0.0;
    //Comparison
    qbtbp_prim_test(xs0_SEM, xs0_EM, ys0_IN, t0, t0c);
    cout << "-------------------------------------------" << endl;

}

/**
 *  \brief Test function used inside qbtbp_test_in_em_sem.
 */
void qbtbp_prim_test(const double xe_0_sem[],
                     const double xe_0_em[],
                     const double ye_0_in[],
                     double t0,
                     double t0c)
{
    //State vectors
    double xem_sem[6];
    double xem_em[6];
    double xe_0_in_from_sem[6];
    double xe_0_in_from_em[6];
    double xe_0_em_2[6];
    double xe_0_sem_from_em[6];
    double xe0_SEM_from_IN[6];

    //------------------------------------------------------------------------------------
    // SEM -> EM -> IN
    //------------------------------------------------------------------------------------
    //Planet position & momenta in SEM ref and SEM units
    se_v_to_se_m(t0c, xe_0_sem, xem_sem, &SEML);
    //Planet position & momenta in EM ref and EM units
    se_m_to_em_m(t0c, xem_sem, xem_em, &SEML);
    //Planet position & velocity in EM ref and EM units
    em_m_to_em_v(t0, xem_em, xe_0_em_2, &SEML);
    //Planet position & velocity in IN ref and EM units
    em_v_to_in(t0, xe_0_em_2, xe_0_in_from_sem,  &SEML);

    //------------------------------------------------------------------------------------
    // EM -> IN
    //------------------------------------------------------------------------------------
    //Planet position & velocities in IN ref and EM units
    em_v_to_in(t0, xe_0_em, xe_0_in_from_em,  &SEML);

    cout << " The state vector of the primary is computed in SEM (Sun-Earth), " << endl;
    cout << " EM (Earth-Moon) and IN (Inertial) coordinates, at time t0 = " << t0 << endl;
    cout << " Then, the SEM and EM states are transposed in IN coordinates " << endl;
    cout << " in order to compare them with the original IN state vector." << endl;
    cout << "-------------------------------------------" << endl;
    cout << "(SEM -> IN) vs IN     (EM->IN) vs IN " << endl;
    for(int i = 0; i < 6; i++)
    {
        cout <<  fabs(xe_0_in_from_sem[i]-ye_0_in[i]) << "          " << fabs(xe_0_in_from_em[i]-ye_0_in[i]) << endl;
    }
    cout << "-------------------------------------------" << endl;

    //------------------------------------------------------------------------------------
    // EM -> IN -> SE
    //------------------------------------------------------------------------------------
    //Planet position & momenta in EM ref and EM units
    em_v_to_em_m(t0, xe_0_em, xem_em, &SEML);
    //Planet position & momenta in SEM ref and SEM units
    em_m_to_se_m(t0, xem_em, xem_sem, &SEML);
    //Planet position & velocity in SEM ref and SEM units
    se_m_to_se_v(t0c, xem_sem, xe_0_sem_from_em, &SEML);
    //Planet position & velocity in SEM ref and SEM units
    in_to_se_v(t0c, ye_0_in, xe0_SEM_from_IN, &SEML);


    cout << "The same is done with EM vs SEM" << endl;
    cout << "-------------------------------------------" << endl;
    cout << "(EM -> SEM) vs SEM " << endl;

    for(int i = 0; i < 6; i++)
    {
        cout <<  fabs(xe_0_sem[i]-xe_0_sem_from_em[i]) << endl;
    }
    cout << "-------------------------------------------" << endl;
}


/**
 *  \brief Comparison of the QBTBP computed via FFT or via OFS (algebraic manipulations).
 */
void qbtbp_test_FFT_vs_OFS(Ofsc& bjc,             //zt(t)
                           Ofsc& cjc,             //Zt(t)
                           int n_order_fourier,   //order of the Fourier expansions
                           int n_points,          //Number of points
                           int type,              //Type of reference
                           OdeStruct ode_s,       //ode structure
                           FBPL& fbpl)            //QBCP
{
    reset_ode_structure(&ode_s);

    //Creation of the matrices of results
    double** alpha_INT_VS_OFS = dmatrix(0, 7, 0, n_points-1);
    double** alpha_INT_VS_FFT = dmatrix(0, 7, 0, n_points-1);
    double** alpha_INT        = dmatrix(0, 7, 0, n_points-1);

    //Retrieving from txt files
    double* alpha_OFS = dvector(0, 8*(n_order_fourier+1)-1);
    double* alpha_FFT = dvector(0, 8*(n_order_fourier+1)-1);


    cout << "-------------------------------------------" << endl;
    cout << "Retrieving the coefficients computed via alebraic operations (OFS) " << endl;
    read_fourier_coef(fbpl.cs.F_COEF+"alpha", alpha_OFS, n_order_fourier, 0, 0, 8);
    cout << "Retrieving the coefficients computed via FFT " << endl;
    read_fourier_coef(fbpl.cs.F_COEF+"alpha", alpha_FFT, n_order_fourier, 0, 1, 8);

    //------------------------------------------------------------------------------------
    //Physical params in EM units
    //------------------------------------------------------------------------------------
    double ns = fbpl.us_em.ns;  //Sun-(Earth+Moon) mean angular motion
    double ni = fbpl.us_em.ni;  //Earth-Moon mean angular motion
    double n  = fbpl.us_em.n;   //n = ni - ns
    double as = fbpl.us_em.as;  //Sun-(Earth+Moon) mean distance
    double ai = fbpl.us_em.ai;  //Earth-Moon mean distance
    double ms = fbpl.us_em.ms;  //Sun mass

    //------------------------------------------------------------------------------------
    //QBTBP init
    //------------------------------------------------------------------------------------
    Ofsc ztc(bjc);
    Ofsc Ztc(cjc);
    //Derivatives
    Ofsc ztcdot(ztc);
    Ofsc Ztcdot(Ztc);
    //Derivation
    ztcdot.dot(n);
    Ztcdot.dot(n);

    //Double derivatives
    Ofsc ztcddot(ztcdot);
    Ofsc Ztcddot(Ztcdot);
    //Derivation
    ztcddot.dot(n);
    Ztcddot.dot(n);

    //z(t) and Z(t)
    cdouble zi;
    cdouble Zi;
    cdouble zidot;
    cdouble Ziddot;

    double rv2;
    double alpha1i;
    double alpha2i;
    double alpha3i;
    double alpha4i;
    double alpha5i;
    double alpha6i;
    double eval_OFS[8];
    double eval_FFT[8];

    //z(0) and Z(0)
    cdouble z0 = bjc.evaluate(0.0);
    cdouble Z0 = as*cjc.evaluate(0.0);

    //zdot(0) and Zdot(0)
    Ofsc zdot(n_order_fourier);
    for(int l = -n_order_fourier; l<=n_order_fourier; l++) zdot.set_coef(I*(1.0+l*n)*bjc.ofs_get_coef(l), l);
    cdouble zdot0 = zdot.evaluate(0.0);

    Ofsc Zdot(n_order_fourier);
    for(int l = -n_order_fourier; l<=n_order_fourier; l++) Zdot.set_coef(I*(ns+l*n)*cjc.ofs_get_coef(l), l);
    cdouble Zdot0 = as*Zdot.evaluate(0.0);

    //Initial conditions
    double yv[8];
    yv[0] = creal(z0);
    yv[1] = cimag(z0);
    yv[2] = creal(Z0);
    yv[3] = cimag(Z0);
    yv[4] = creal(zdot0);
    yv[5] = cimag(zdot0);
    yv[6] = creal(Zdot0);
    yv[7] = cimag(Zdot0);

    //Loop
    double ti;
    double tf  = 2*M_PI/n;
    double t   = 0.0;
    int status;
    double f[8];
    for(int i = 0 ; i < n_points; i++)
    {
        //update
        ti = tf*i/n_points;

        //Computing the reference values
        if(type == 0)
        {
            //z(t), zdot(t), Z(t) and Zdot(t) if the type of reference is the numerical integration  (type = 0)
            //Evolve one step
            status = gsl_odeiv2_driver_apply (ode_s.d, &t, ti, yv);
            //Compute z(t) and Z(t)
            zi     =  yv[0] + yv[1]*I;
            Zi     =  yv[2] + yv[3]*I;
            zidot  =  yv[4] + yv[5]*I;
            qbtbp_derivatives(t, yv, f, ode_s.sys.params);
            Ziddot =  f[6] + f[7]*I;
            rv2    =  creal(1.0/(zi*conj(zi)));
        }
        else
        {
            //z(t), zdot(t), Z(t) and Zdot(t) if the type of reference is the Fourier expansion (type = 1)
            zi     =  eval_z(ztc, t, n, ni, ai);
            Zi     =  eval_z(Ztc, t, n, ns, as);
            zidot  =  eval_zdot(ztc, ztcdot, t, n, ni, ai);
            Ziddot =  eval_zddot(Ztc, Ztcdot, Ztcddot, t, n, ns, as);
            rv2    =  creal(1.0/(zi*conj(zi)));
        }


        //Evalutation of FFT and OFS results
        eval_array_coef(eval_OFS, t, n, n_order_fourier, alpha_OFS, 8);
        eval_array_coef(eval_FFT, t, n, n_order_fourier, alpha_FFT, 8);

        alpha1i  = rv2;
        alpha2i  = -rv2*creal(zidot*conj(zi));
        alpha3i  = +rv2*cimag(zidot*conj(zi));
        alpha4i  = -ms/(1.0+ms)*creal(Ziddot*conj(zi));
        alpha5i  = -ms/(1.0+ms)*cimag(Ziddot*conj(zi));
        alpha6i  = creal(cpow(zi*conj(zi), -1.0/2+0.0*I));

        //Results in matrices
        alpha_INT_VS_OFS[0][i] = fabs(alpha1i - eval_OFS[0]);
        alpha_INT_VS_OFS[1][i] = fabs(alpha2i - eval_OFS[1]);
        alpha_INT_VS_OFS[2][i] = fabs(alpha3i - eval_OFS[2]);
        alpha_INT_VS_OFS[3][i] = fabs(alpha4i - eval_OFS[3]);
        alpha_INT_VS_OFS[4][i] = fabs(alpha5i - eval_OFS[4]);
        alpha_INT_VS_OFS[5][i] = fabs(alpha6i - eval_OFS[5]);
        alpha_INT_VS_OFS[6][i] = fabs(rv2*creal(Zi*conj(zi)) - eval_OFS[6]);
        alpha_INT_VS_OFS[7][i] = fabs(rv2*cimag(Zi*conj(zi)) - eval_OFS[7]);

        alpha_INT_VS_FFT[0][i] = fabs(alpha1i - eval_FFT[0]);
        alpha_INT_VS_FFT[1][i] = fabs(alpha2i - eval_FFT[1]);
        alpha_INT_VS_FFT[2][i] = fabs(alpha3i - eval_FFT[2]);
        alpha_INT_VS_FFT[3][i] = fabs(alpha4i - eval_FFT[3]);
        alpha_INT_VS_FFT[4][i] = fabs(alpha5i - eval_FFT[4]);
        alpha_INT_VS_FFT[5][i] = fabs(alpha6i - eval_FFT[5]);
        alpha_INT_VS_FFT[6][i] = fabs(rv2*creal(Zi*conj(zi)) - eval_FFT[6]);
        alpha_INT_VS_FFT[7][i] = fabs(rv2*cimag(Zi*conj(zi)) - eval_FFT[7]);

        alpha_INT[0][i] = alpha1i;
        alpha_INT[1][i] = alpha2i;
        alpha_INT[2][i] = alpha3i;
        alpha_INT[3][i] = alpha4i;
        alpha_INT[4][i] = alpha5i;
        alpha_INT[5][i] = alpha6i;
        alpha_INT[6][i] = rv2*creal(Zi*conj(zi));
        alpha_INT[7][i] = rv2*cimag(Zi*conj(zi));

        //Break if evolution has gone wrong
        if (status != GSL_SUCCESS) break;
    }

    //Maximum
    double max_INT_VS_OFS[8];
    double max_INT_VS_FFT[8];
    double max_INT[8];
    for(int j = 0; j < 8; j ++)
    {
        max_INT_VS_OFS[j] = alpha_INT_VS_OFS[j][0];
        max_INT_VS_FFT[j] = alpha_INT_VS_FFT[j][0];
        max_INT[j]        = alpha_INT[j][0];
    }

    for(int i = 1 ; i < n_points; i++)
    {
        for(int j = 0; j < 8; j++)
        {
            if(max_INT_VS_OFS[j] < alpha_INT_VS_OFS[j][i]) max_INT_VS_OFS[j] = alpha_INT_VS_OFS[j][i];
            if(max_INT_VS_FFT[j] < alpha_INT_VS_FFT[j][i]) max_INT_VS_FFT[j] = alpha_INT_VS_FFT[j][i];
            if(max_INT[j] < alpha_INT[j][i]) max_INT[j] = alpha_INT[j][i];
        }
    }

    //Print
    cout << "-------------------------------------------" << endl;
    cout << "      Comparison: FFT vs OFS               " << endl;
    cout << "-------------------------------------------" << endl;
    cout << " - The first 8 coefficients (alpha functions) of the         " << endl;
    cout << " QBCP vector field are compared against a common reference." << endl;
    cout << " - Only the maxima are compared." << endl;

    if(type == 0)
    {
        cout << "The numerical integration is used to compute the reference values (REF) " << endl;
    }
    else
    {
        cout << "Numerical evaluations of the Fourier expansions of z(t) and Z(t) " << endl;
        cout << "are used to compute the reference values (REF) " << endl;
    }
    cout << "i         max(REF)        |max(REF) - max(OFS)|    |max(REF) - max(FFT)|       " << endl;
    for(int j = 0; j < 8; j++)
    {
        cout << j << "      " << max_INT[j] << "           " << max_INT_VS_OFS[j]  << "            " <<  max_INT_VS_FFT[j]  << endl;
    }
    cout << "-------------------------------------------" << endl;

    //Conclusion
    cout << "Conclusion: " << endl;
    cout << "- OFS and FFT coefficients both represent quite well " << endl;
    cout << " the real integrated system  (see comparison with type == 0). " << endl;
    cout << "- However, the FFT coefficients are better suited " << endl;
    cout << " to represent the approximated system given by " << endl;
    cout << " the Fourier expansions computed in qbtbp (see comparison with type == 1). " << endl;
    cout << "- This is ought to the approximations made to compute z(OFS)^{-a}, a > 0, " << endl;
    cout << " when computing the OFS coefficients." << endl;
    cout << "-------------------------------------------" << endl;

    //Memory release
    free_dvector(alpha_OFS, 0, 8*(n_order_fourier+1)-1);
    free_dvector(alpha_FFT, 0, 8*(n_order_fourier+1)-1);
    free_dmatrix(alpha_INT_VS_OFS, 0, 7, 0, n_points-1);
    free_dmatrix(alpha_INT       , 0, 7, 0, n_points-1);
    free_dmatrix(alpha_INT_VS_FFT, 0, 7, 0, n_points-1);
}

//----------------------------------------------------------------------------------------
// Evaluating the QBTBP
//----------------------------------------------------------------------------------------
/**
 *  \brief Evaluate z(t), with \f$ z(t) = e^{it} z_r(t) \f$ in Earth-Moon units.
 */
cdouble eval_z(Ofsc& zt, double t, double n, double ni, double ai)
{
    return ai*(cos(ni*t)+I*sin(ni*t))*zt.evaluate(n*t);
}

/**
 *  \brief Evaluate dz(t)/dt, with \f$ z(t) = e^{it} z_r(t) \f$ in Earth-Moon units.
 */
cdouble eval_zdot(Ofsc& zt, Ofsc& ztdot, double t, double n, double ni, double ai)
{
    return ai*(cos(ni*t)+I*sin(ni*t))*(ztdot.evaluate(n*t) + I*ni*zt.evaluate(n*t));
}

/**
 *  \brief Evaluate d2z(t)/dt2, with \f$ z(t) = e^{it} z_r(t) \f$ in Earth-Moon units.
 */
cdouble eval_zddot(Ofsc& zt, Ofsc& ztdot, Ofsc& ztddot, double t, double n, double ni, double ai)
{
    return ai*(cos(ni*t)+I*sin(ni*t))*( 2*I*ni*ztdot.evaluate(n*t) - ni*ni*zt.evaluate(n*t) + ztddot.evaluate(n*t));
}


//----------------------------------------------------------------------------------------
//Evaluation of parts or all the alpha/beta routines at a given time
//----------------------------------------------------------------------------------------
/**
 *  \brief Evaluation of Fourier series given as an array of coefficients,
 *         at a given time t.
 */
void eval_array_coef(double* alpha, double t, double omega, int n_order, double* params, int number)
{
    double* header;
    int l;
    double cR[n_order];
    double sR[n_order];

    cR[0] = cos(omega*t);
    sR[0] = sin(omega*t);

    for(int i = 1; i< n_order; i++)
    {
        //From trigo formulae (fastest)
        cR[i] =  cR[i-1]*cR[0] - sR[i-1]*sR[0];
        sR[i] =  sR[i-1]*cR[0] + cR[i-1]*sR[0];
        //From GSL
        //cR[i] =  gsl_sf_cos((i+1)*omega*t);
        //sR[i] =  gsl_sf_sin((i+1)*omega*t);
        //Native C
        //cR[i] =  cos((i+1)*omega*t);
        //sR[i] =  sin((i+1)*omega*t);
    }

    for(l = 0, header = params; l < number ; l++, header+=(n_order+1))
    {
        alpha[l] = 0.0;
        if(l==1 || l== 4 || l==7 || l==9 || l==11 || l==13) alpha[l] = eval_odd_trigo(n_order, header, sR);    //Odd funtions (alpha_2,5,8,10,12,14)
        else  alpha[l] = eval_even_trigo(n_order, header, cR);                                                 //Even functions
    }
}

/**
 *  \brief Evaluation of the time derivatives Fourier series given as an array of
 *         coefficients, at a given time t.
 */
void eval_array_coef_der(double* alpha, double t, double omega, int n_order, double* params, int number)
{
    double* header;
    int l;
    double cR[n_order];
    double sR[n_order];

    cR[0] = cos(omega*t);
    sR[0] = sin(omega*t);

    for(int i = 1; i< n_order; i++)
    {
        cR[i] =  cR[i-1]*cR[0] - sR[i-1]*sR[0];
        sR[i] =  sR[i-1]*cR[0] + cR[i-1]*sR[0];
    }

    for(l = 0, header = params; l < number ; l++, header+=(n_order+1))
    {
        alpha[l] = 0.0;
        if(l==1 || l== 4 || l==7 || l==9 || l==11 || l==13)
        {
            alpha[l] = eval_odd_trigo_der(omega , n_order, header, cR); //Odd funtions (alpha_2,5,8,10,12,14)
        }
        else  alpha[l] = eval_even_trigo_der(omega, n_order, header, sR);   //Even functions
    }
}


//----------------------------------------------------------------------------------------
//Evalution of individual alpha/beta
//or derivative of alpha/beta of a given type (Odd or Even)
//----------------------------------------------------------------------------------------
/**
 *  \brief Evaluate the sum \f$ \sum_{k = 0}^N coef(k) cos(k \omega t)  \f$.
 */
double eval_even_trigo(int n_order, double* coef, double* cR)
{
    double result = 0.0;
    for(int i= n_order; i>=1; i--) result += coef[i]*cR[i-1];//even type
    result += coef[0];
    return result;
}

/**
 *  \brief Evaluate the sum \f$ \sum_{k = 0}^N - k \omega coef(k) sin(k \omega t)  \f$.
 */
double eval_even_trigo_der(double omega,  int n_order, double* coef, double* sR)
{
    double result = 0.0;
    for(int i= n_order; i>=1; i--) result += -omega*i*coef[i]*sR[i-1];//even type
    return result;
}

/**
 *  \brief Evaluate the sum \f$ \sum_{k = 0}^N coef(k) sin(k \omega t)  \f$.
 */
double eval_odd_trigo(int n_order, double* coef, double* sR)
{
    double result = 0.0;
    for(int i= n_order; i>=1; i--) result += coef[i]*sR[i-1]; //odd type
    return result;
}

/**
 *  \brief Evaluate the sum \f$ \sum_{k = 0}^N  k \omega coef(k) cos(k \omega t)  \f$.
 */
double eval_odd_trigo_der( double omega,  int n_order, double* coef, double* cR)
{
    double result = 0.0;
    for(int i= n_order; i>=1; i--) result += omega*i*coef[i]*cR[i-1];//odd type
    return result;
}


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
void bcp()
{
    //------------------------------------------------------------------------------------
    // Splash
    //------------------------------------------------------------------------------------
    cout << "---------------------------------------------------" << endl;
    cout << "                                                   " << endl;
    cout << "       Resolution of the Sun-Earth-Moon            " << endl;
    cout << "      Bicircular Four-Body Problem (BCP)           " << endl;
    cout << "                                                   " << endl;
    cout << "---------------------------------------------------" << endl;
    Config::configManager().coutmp();


    //------------------------------------------------------------------------------------
    // Initialization of environement via an FBPL structure
    //------------------------------------------------------------------------------------
    cout << " bcp. Initializing the environment..." << endl;
    //Init the Sun-Earth-Moon FBP structure
    FBP fbp;
    init_fbp(&fbp, Csts::SUN, Csts::EARTH, Csts::MOON);

    //Init the Sun-Earth-Moon FBPL structure, with arbitrary configuration
    FBPL fbpl;
    init_fbp_lib(&fbpl, &fbp, 1, 1, Csts::BCP, Csts::EM, Csts::GRAPH,
              Csts::MAN_CENTER, Csts::MAN_CENTER, true, true);


    //------------------------------------------------------------------------------------
    // 3. Alphas
    //------------------------------------------------------------------------------------
    //Computation of the coeffs of the vector field for EML1
    cout << " bcp. Storing the vector field about EML1 in " << fbpl.cs_em_l1.F_COEF << endl;
    bcp_alpha(OFS_ORDER, fbpl.us_em, fbpl.cs_em_l1);

    //Computation of the coeffs of the vector field for EML2
    cout << " bcp. Storing the vector field about EML2 in " << fbpl.cs_em_l2.F_COEF << endl;
    bcp_alpha(OFS_ORDER, fbpl.us_em, fbpl.cs_em_l2);

    //------------------------------------------------------------------------------------
    // 4. Deltas
    //------------------------------------------------------------------------------------
    //Computation of the coeffs of the vector field for SEL1
    cout << " bcp. Storing the vector field about SEL1 in " << fbpl.cs_sem_l1.F_COEF << endl;
    bcp_delta(OFS_ORDER, fbpl.us_sem, fbpl.cs_sem_l1); //libration points on the Sun-Bem line

    //Computation of the coeffs of the vector field for SEL2
    cout << " bcp. Storing the vector field about SEL2 in " << fbpl.cs_sem_l2.F_COEF << endl;
    bcp_delta(OFS_ORDER, fbpl.us_sem, fbpl.cs_sem_l2); //libration points on the Sun-Bem line
}

/**
 *  \brief Computes and stores the coefficients of the vector field of the BCP from the
 *         Earth-Moon point of view.
 **/
void bcp_alpha(int n_order_fourier, USYS &us_em, CSYS &cs_em)
{
    //------------------------------------------------------------------------------------
    //Parameters for one specific libration point
    //------------------------------------------------------------------------------------
    double mu_EM = us_em.mu_EM;
    double gamma = cs_em.gamma;
    double c1    = cs_em.c1;
    double ms    = us_em.ms;
    double as    = us_em.as;

    //------------------------------------------------------------------------------------
    // 3. Set all alpha functions
    //------------------------------------------------------------------------------------
    Ofs<cdouble > alpha1c(n_order_fourier);
    Ofs<cdouble > alpha2c(n_order_fourier);
    Ofs<cdouble > alpha3c(n_order_fourier);
    Ofs<cdouble > alpha4c(n_order_fourier);
    Ofs<cdouble > alpha5c(n_order_fourier);
    Ofs<cdouble > alpha6c(n_order_fourier);
    Ofs<cdouble > alpha7c(n_order_fourier);
    Ofs<cdouble > alpha8c(n_order_fourier);
    Ofs<cdouble > alpha9c(n_order_fourier);
    Ofs<cdouble > alpha10c(n_order_fourier);
    Ofs<cdouble > alpha11c(n_order_fourier);
    Ofs<cdouble > alpha12c(n_order_fourier);
    Ofs<cdouble > alpha13c(n_order_fourier);
    Ofs<cdouble > alpha14c(n_order_fourier);

    //Redundancy for the positions of the primaries
    Ofs<cdouble > Xe(n_order_fourier);
    Ofs<cdouble > Ye(n_order_fourier);
    Ofs<cdouble > Ze(n_order_fourier);
    Ofs<cdouble > Xm(n_order_fourier);
    Ofs<cdouble > Ym(n_order_fourier);
    Ofs<cdouble > Zm(n_order_fourier);
    Ofs<cdouble > Xs(n_order_fourier);
    Ofs<cdouble > Ys(n_order_fourier);
    Ofs<cdouble > Zs(n_order_fourier);

    //Positions of the primaries in NC coordinates
    Ofs<cdouble > xe(n_order_fourier);
    Ofs<cdouble > ye(n_order_fourier);
    Ofs<cdouble > ze(n_order_fourier);
    Ofs<cdouble > xm(n_order_fourier);
    Ofs<cdouble > ym(n_order_fourier);
    Ofs<cdouble > zm(n_order_fourier);
    Ofs<cdouble > xs(n_order_fourier);
    Ofs<cdouble > ys(n_order_fourier);
    Ofs<cdouble > zs(n_order_fourier);

    //-------------------------
    // The alphas
    //-------------------------
    //alpha1 = 1
    alpha1c.set_coef(1.0+0.0*I, 0);
    //alpha2 = 0
    //alpha3 = 1
    alpha3c.set_coef(1.0+0.0*I, 0);
    //alpha4 = ms/as^2*cos(theta) = ms/(2*as^2)*(exp(itheta) + exp(-itheta))
    alpha4c.set_coef(0.5*ms/(as*as)+0.0*I, -1);
    alpha4c.set_coef(0.5*ms/(as*as)+0.0*I, +1);
    //alpha5 = -ms/as^2*sin(theta) = +I*ms/(2*as^2)*(exp(itheta) - exp(-itheta))
    alpha5c.set_coef(-0.5*ms/(as*as)*I, -1);
    alpha5c.set_coef(+0.5*ms/(as*as)*I, +1);
    //alpha6 = 1
    alpha6c.set_coef(1.0+0.0*I, 0);
    //alpha13 = alpha4/gamma-c1
    alpha13c.ofs_mult(alpha4c, 1.0/gamma+0.0*I);
    alpha13c.add_coef(-c1+0.0*I, 0);
    //alpha14 = alpha5/gamma
    alpha14c.ofs_mult(alpha5c, 1.0/gamma+0.0*I);

    //-------------------------
    // Primaries
    //-------------------------
    //Earth
    alpha9c.set_coef(+mu_EM+0.0*I, 0);
    //Moon
    alpha11c.set_coef(+mu_EM-1.0+0.0*I, 0);
    //Sun
    //alpha7 = as*cos(theta) = as/2*(exp(itheta) + exp(-itheta))
    alpha7c.set_coef(0.5*as+0.0*I, -1);
    alpha7c.set_coef(0.5*as+0.0*I, +1);
    //alpha8 = -as*sin(theta) = +I*as/2*(exp(itheta) - exp(-itheta))
    alpha8c.set_coef(-0.5*as*I, -1);
    alpha8c.set_coef(+0.5*as*I, +1);

    //-------------------------
    // The primaries, again
    //-------------------------
    //Earth, EM coordinates
    Xe.set_coef(+mu_EM+0.0*I, 0);
    //Earth, NC coordinates
    xe.set_coef(c1 - mu_EM/gamma + 0.0*I, 0);
    //Moon, EM coordinates
    Xm.set_coef(+mu_EM-1.0+0.0*I, 0);
    //Moon, NC coordinates
    xm.set_coef(c1 - (mu_EM-1.0)/gamma + 0.0*I, 0);

    //Sun, EM coordinates
    //--------------------------------------------------------------------------------
    //Xs = as*cos(theta) = as/2*(exp(itheta) + exp(-itheta))
    Xs.set_coef(0.5*as+0.0*I, -1);
    Xs.set_coef(0.5*as+0.0*I, +1);
    //Ys = -as*sin(theta) = +I*as/2*(exp(itheta) - exp(-itheta))
    Ys.set_coef(-0.5*as*I, -1);
    Ys.set_coef(+0.5*as*I, +1);
    //Sun, NC coordinates;
    //--------------------------------------------------------------------------------
    xs.ofs_smult(Xs, -1.0/gamma);
    xs.add_coef(c1, 0);
    ys.ofs_smult(Ys, -1.0/gamma);

    //------------------------------------------------------------------------------------
    //Put in data file
    //------------------------------------------------------------------------------------
    ofs_sst(alpha1c, cs_em.F_COEF+"alpha1", 1, "_fft");
    ofs_sst(alpha2c, cs_em.F_COEF+"alpha2", 0, "_fft");
    ofs_sst(alpha3c, cs_em.F_COEF+"alpha3", 1, "_fft");
    ofs_sst(alpha4c, cs_em.F_COEF+"alpha4", 1, "_fft");
    ofs_sst(alpha5c, cs_em.F_COEF+"alpha5", 0, "_fft");
    ofs_sst(alpha6c, cs_em.F_COEF+"alpha6", 1, "_fft");
    //Sun
    ofs_sst(alpha7c, cs_em.F_COEF+"alpha7", 1, "_fft");
    ofs_sst(alpha8c, cs_em.F_COEF+"alpha8", 0, "_fft");
    //Earth
    ofs_sst(alpha9c,  cs_em.F_COEF+"alpha9",  1, "_fft");
    ofs_sst(alpha10c, cs_em.F_COEF+"alpha10", 0, "_fft");
    //Moon
    ofs_sst(alpha11c, cs_em.F_COEF+"alpha11", 1, "_fft");
    ofs_sst(alpha12c, cs_em.F_COEF+"alpha12", 0, "_fft");
    //NC additional coeffs
    ofs_sst(alpha13c, cs_em.F_COEF+"alpha13", 1, "_fft");
    ofs_sst(alpha14c, cs_em.F_COEF+"alpha14", 0, "_fft");

    //--------------------------------------------------------------------------------
    //Primary, EM coordinates
    //Note that, at this step, the vertical motion of the primaries is undefined,
    //so we can put either Even or Odd in the ofs_sst option of gsl_Zc without much difference
    //--------------------------------------------------------------------------------
    ofs_sst(Xs, cs_em.F_COEF+"Ps1", 1, "_fft");
    ofs_sst(Ys, cs_em.F_COEF+"Ps2", 0, "_fft");
    ofs_sst(Zs, cs_em.F_COEF+"Ps3", 1, "_fft");

    ofs_sst(Xm, cs_em.F_COEF+"Pm1", 1, "_fft");
    ofs_sst(Ym, cs_em.F_COEF+"Pm2", 0, "_fft");
    ofs_sst(Zm, cs_em.F_COEF+"Pm3", 1, "_fft");

    ofs_sst(Xe, cs_em.F_COEF+"Pe1", 1, "_fft");
    ofs_sst(Ye, cs_em.F_COEF+"Pe2", 0, "_fft");
    ofs_sst(Ze, cs_em.F_COEF+"Pe3", 1, "_fft");


    //--------------------------------------------------------------------------------
    //Primary, NC coordinates
    //Note that, at this step, the vertical motion of the primaries is undefined,
    //so we can put either Even or Odd in the ofs_sst option of gsl_zc without much difference
    //--------------------------------------------------------------------------------
    ofs_sst(xs, cs_em.F_COEF+"ps1", 1, "_fft");
    ofs_sst(ys, cs_em.F_COEF+"ps2", 0, "_fft");
    ofs_sst(zs, cs_em.F_COEF+"ps3", 1, "_fft");

    ofs_sst(xm, cs_em.F_COEF+"pm1", 1, "_fft");
    ofs_sst(ym, cs_em.F_COEF+"pm2", 0, "_fft");
    ofs_sst(zm, cs_em.F_COEF+"pm3", 1, "_fft");

    ofs_sst(xe, cs_em.F_COEF+"pe1", 1, "_fft");
    ofs_sst(ye, cs_em.F_COEF+"pe2", 0, "_fft");
    ofs_sst(ze, cs_em.F_COEF+"pe3", 1, "_fft");
}

/**
 *  \brief Computes and stores the coefficients of the vector field of the BCP from the
 *         Sun-Earth point of view.
 **/
void bcp_delta(int n_order_fourier, USYS &us_sem, CSYS &cs_sem)
{
    //------------------------------------------------------------------------------------
    //Parameters for one specific libration point
    //------------------------------------------------------------------------------------
    double mu_EM  = us_sem.mu_EM;
    double mu_SE  = us_sem.mu_SEM;
    double gamma  = cs_sem.gamma;
    double c1     = cs_sem.c1;
    //double mm     = us_sem.mm;
    //double me     = us_sem.me;
    double ai     = us_sem.ai;

    //Distance of Earth & Moon from their barycenter
    double am = (1- mu_EM)*ai;
    double ae = mu_EM*ai;
    //Factor that appears in the definition of the deltas.
    //double fem = mm/(am*am) - me/(ae*ae);

    //------------------------------------------------------------------------------------
    // 3. Set all delta functions
    //------------------------------------------------------------------------------------
    Ofs<cdouble > delta1c(n_order_fourier);
    Ofs<cdouble > delta2c(n_order_fourier);
    Ofs<cdouble > delta3c(n_order_fourier);
    Ofs<cdouble > delta4c(n_order_fourier);
    Ofs<cdouble > delta5c(n_order_fourier);
    Ofs<cdouble > delta6c(n_order_fourier);
    Ofs<cdouble > delta7c(n_order_fourier);
    Ofs<cdouble > delta8c(n_order_fourier);
    Ofs<cdouble > delta9c(n_order_fourier);
    Ofs<cdouble > delta10c(n_order_fourier);
    Ofs<cdouble > delta11c(n_order_fourier);
    Ofs<cdouble > delta12c(n_order_fourier);
    Ofs<cdouble > delta13c(n_order_fourier);
    Ofs<cdouble > delta14c(n_order_fourier);

    //Redundancy for the positions of the primaries
    Ofs<cdouble > Xe(n_order_fourier);
    Ofs<cdouble > Ye(n_order_fourier);
    Ofs<cdouble > Ze(n_order_fourier);
    Ofs<cdouble > Xm(n_order_fourier);
    Ofs<cdouble > Ym(n_order_fourier);
    Ofs<cdouble > Zm(n_order_fourier);
    Ofs<cdouble > Xs(n_order_fourier);
    Ofs<cdouble > Ys(n_order_fourier);
    Ofs<cdouble > Zs(n_order_fourier);

    //Positions of the primaries in NC coordinates
    Ofs<cdouble > xe(n_order_fourier);
    Ofs<cdouble > ye(n_order_fourier);
    Ofs<cdouble > ze(n_order_fourier);
    Ofs<cdouble > xm(n_order_fourier);
    Ofs<cdouble > ym(n_order_fourier);
    Ofs<cdouble > zm(n_order_fourier);
    Ofs<cdouble > xs(n_order_fourier);
    Ofs<cdouble > ys(n_order_fourier);
    Ofs<cdouble > zs(n_order_fourier);

    //-------------------------
    // The deltas
    //-------------------------
    //delta1 = 1
    delta1c.set_coef(1.0+0.0*I, 0);
    //delta2 = 0
    //delta3 = 1
    delta3c.set_coef(1.0+0.0*I, 0);
    //delta4 = 0; // OR -fem*cos(theta) = -0.5*fem*(exp(itheta) + exp(-itheta))
    //delta4c.set_coef(-0.5*fem+0.0*I, -1);
    //delta4c.set_coef(-0.5*fem+0.0*I, +1);
    //delta5 = 0; // OR -fem*sin(theta) = +I*0.5*fem*(exp(itheta) - exp(-itheta))
    //delta5c.set_coef(-0.5*fem*I, -1);
    //delta5c.set_coef(+0.5*fem*I, +1);
    //delta6 = 1
    delta6c.set_coef(1.0+0.0*I, 0);
    //delta13 = delta4/gamma-c1
    delta13c.ofs_mult(delta4c, 1.0/gamma+0.0*I);
    delta13c.add_coef(-c1+0.0*I, 0);
    //delta14 = delta5/gamma
    delta14c.ofs_mult(delta5c, 1.0/gamma+0.0*I);

    //-------------------------
    // Primaries
    //-------------------------
    //Earth
    //delta9 = +ae*cos(theta) + mu_SE-1.0  = ae/2*(exp(itheta) + exp(-itheta)) + mu_SE-1.0
    delta9c.set_coef(mu_SE-1.0 +0.0*I, 0);
    delta9c.set_coef(0.5*ae+0.0*I, -1);
    delta9c.set_coef(0.5*ae+0.0*I, +1);
    //delta10 = +ae*sin(theta) = -I*ae/2*(exp(itheta) - exp(-itheta))
    delta10c.set_coef(+0.5*ae*I, -1);
    delta10c.set_coef(-0.5*ae*I, +1);

    //Sun
    delta7c.set_coef(+mu_SE+0.0*I, 0);

    //Moon
    //delta11 = -am*cos(theta) + mu_SE-1.0  = -am/2*(exp(itheta) + exp(-itheta)) + mu_SE-1.0
    delta11c.set_coef(mu_SE-1.0 +0.0*I, 0);
    delta11c.set_coef(-0.5*am+0.0*I, -1);
    delta11c.set_coef(-0.5*am+0.0*I, +1);
    //delta12 = -am*sin(theta) = +I*am/2*(exp(itheta) - exp(-itheta))
    delta12c.set_coef(-0.5*am*I, -1);
    delta12c.set_coef(+0.5*am*I, +1);

    //-------------------------
    // The primaries, again
    //-------------------------
    //Earth, EM coordinates
    //--------------------------------------------------------------------------------
    //Xe = +ae*cos(theta) + mu_SE-1.0  = ae/2*(exp(itheta) + exp(-itheta)) + mu_SE-1.0
    Xe.set_coef(mu_SE-1.0 +0.0*I, 0);
    Xe.set_coef(0.5*ae+0.0*I, -1);
    Xe.set_coef(0.5*ae+0.0*I, +1);
    //Ye = +ae*sin(theta) = -I*ae/2*(exp(itheta) - exp(-itheta))
    Ye.set_coef(+0.5*ae*I, -1);
    Ye.set_coef(-0.5*ae*I, +1);
    //Earth, NC coordinates;
    //--------------------------------------------------------------------------------
    xe.ofs_smult(Xe, -1.0/gamma);
    xe.add_coef(c1, 0);
    ye.ofs_smult(Ye, -1.0/gamma);

    //Sun, EM coordinates
    Xs.set_coef(+mu_SE+0.0*I, 0);
    //Sun, NC coordinates
    xs.set_coef(c1 - mu_SE/gamma + 0.0*I, 0);

    //Moon, EM coordinates
    //--------------------------------------------------------------------------------
    //Xm = -am*cos(theta) + mu_SE-1.0  = am/2*(exp(itheta) + exp(-itheta)) + mu_SE-1.0
    Xm.set_coef(mu_SE-1.0 +0.0*I, 0);
    Xm.set_coef(-0.5*am+0.0*I, -1);
    Xm.set_coef(-0.5*am+0.0*I, +1);
    //Ym = -am*sin(theta) = +I*am/2*(exp(itheta) - exp(-itheta))
    Ym.set_coef(-0.5*am*I, -1);
    Ym.set_coef(+0.5*am*I, +1);
    //Moon, NC coordinates;
    //--------------------------------------------------------------------------------
    xm.ofs_smult(Xm, -1.0/gamma);
    xm.add_coef(c1, 0);
    ym.ofs_smult(Ym, -1.0/gamma);

    //------------------------------------------------------------------------------------
    //Put in data file
    //------------------------------------------------------------------------------------
    ofs_sst(delta1c, cs_sem.F_COEF+"alpha1", 1, "_fft");
    ofs_sst(delta2c, cs_sem.F_COEF+"alpha2", 0, "_fft");
    ofs_sst(delta3c, cs_sem.F_COEF+"alpha3", 1, "_fft");
    ofs_sst(delta4c, cs_sem.F_COEF+"alpha4", 1, "_fft");
    ofs_sst(delta5c, cs_sem.F_COEF+"alpha5", 0, "_fft");
    ofs_sst(delta6c, cs_sem.F_COEF+"alpha6", 1, "_fft");
    //Sun
    ofs_sst(delta7c, cs_sem.F_COEF+"alpha7", 1, "_fft");
    ofs_sst(delta8c, cs_sem.F_COEF+"alpha8", 0, "_fft");
    //Earth
    ofs_sst(delta9c,  cs_sem.F_COEF+"alpha9",  1, "_fft");
    ofs_sst(delta10c, cs_sem.F_COEF+"alpha10", 0, "_fft");
    //Moon
    ofs_sst(delta11c, cs_sem.F_COEF+"alpha11", 1, "_fft");
    ofs_sst(delta12c, cs_sem.F_COEF+"alpha12", 0, "_fft");
    //NC additional coeffs
    ofs_sst(delta13c, cs_sem.F_COEF+"alpha13", 1, "_fft");
    ofs_sst(delta14c, cs_sem.F_COEF+"alpha14", 0, "_fft");

    //--------------------------------------------------------------------------------
    //Primary, EM coordinates
    //Note that, at this step, the vertical motion of the primaries is undefined,
    //so we can put either Even or Odd in the ofs_sst option of gsl_Zc without much difference
    //--------------------------------------------------------------------------------
    ofs_sst(Xs, cs_sem.F_COEF+"Ps1", 1, "_fft");
    ofs_sst(Ys, cs_sem.F_COEF+"Ps2", 0, "_fft");
    ofs_sst(Zs, cs_sem.F_COEF+"Ps3", 1, "_fft");

    ofs_sst(Xm, cs_sem.F_COEF+"Pm1", 1, "_fft");
    ofs_sst(Ym, cs_sem.F_COEF+"Pm2", 0, "_fft");
    ofs_sst(Zm, cs_sem.F_COEF+"Pm3", 1, "_fft");

    ofs_sst(Xe, cs_sem.F_COEF+"Pe1", 1, "_fft");
    ofs_sst(Ye, cs_sem.F_COEF+"Pe2", 0, "_fft");
    ofs_sst(Ze, cs_sem.F_COEF+"Pe3", 1, "_fft");


    //--------------------------------------------------------------------------------
    //Primary, NC coordinates
    //Note that, at this step, the vertical motion of the primaries is undefined,
    //so we can put either Even or Odd in the ofs_sst option of gsl_zc without much difference
    //--------------------------------------------------------------------------------
    ofs_sst(xs, cs_sem.F_COEF+"ps1", 1, "_fft");
    ofs_sst(ys, cs_sem.F_COEF+"ps2", 0, "_fft");
    ofs_sst(zs, cs_sem.F_COEF+"ps3", 1, "_fft");

    ofs_sst(xm, cs_sem.F_COEF+"pm1", 1, "_fft");
    ofs_sst(ym, cs_sem.F_COEF+"pm2", 0, "_fft");
    ofs_sst(zm, cs_sem.F_COEF+"pm3", 1, "_fft");

    ofs_sst(xe, cs_sem.F_COEF+"pe1", 1, "_fft");
    ofs_sst(ye, cs_sem.F_COEF+"pe2", 0, "_fft");
    ofs_sst(ze, cs_sem.F_COEF+"pe3", 1, "_fft");
}

