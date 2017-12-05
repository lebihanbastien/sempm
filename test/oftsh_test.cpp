/**
 * \file oftsh_test.cpp
 * \brief Test file for the Oftsh (Homogeneous Fourier-Taylor polynomials) class.
 * \author BLB
 */

#include "timec.h"
#include "oftsh_test.h"
#define VAR 2.

using namespace std;


/**
 * \fn void oftsh_test()
 * \brief Main routine for Oftsh class test.
 */
void oftsh_test()
{
    //Check the the OFTS_ORDER parameter is even, to ease the computations
    if(OFTS_ORDER %2 != 0)
    {
        cout << "Error in oftsh_test: the parameter OFTS_ORDER must be even." << endl;
        return;
    }

    char ch;
    //------------------------------------------------------------------------------------
    // Splash
    //------------------------------------------------------------------------------------
    cout << "---------------------------------------------------" << endl;
    cout << "                                                   " << endl;
    cout << "         Test routine for the Oftsh class          " << endl;
    cout << "        specialized in the form Oftsh<Ofs>         " << endl;
    cout << "      (Fourier-Taylor homogeneous polynomials)     " << endl;
    cout << "                                                   " << endl;
    cout << "---------------------------------------------------" << endl;
    cout << "The following public functions are NOT tested:     " << endl;
    cout << "---------------------------------------------------" << endl;
    cout << "- The constructors" << endl;
    cout << "- The destructor " << endl;
    cout << "- Oftsh<T>& Oftsh<T>::lcopy (Oftsh<T> const& b)" << endl;
    cout << "- Oftsh<T>& Oftsh<T>::ccopy (Oftsh<T> const& b)" << endl;
    cout << "- Oftsh<T>::link_coefs(T *coef0)" << endl;
    cout << "- Oftsh<T>::set_coef(T value, int pos)" << endl;
    cout << "- Oftsh<T>::add_coef(T value, int pos)" << endl;
    cout << "- Oftsh<T>::set_sub_coef(U value, int pos)" << endl;
    cout << "- Oftsh<T>::set_sub_coef(U value, int pos, int i)" << endl;
    cout << "- Oftsh<T>::set_random_coefs()" << endl;
    cout << "- Oftsh<T>::get_term() const" << endl;
    cout << "- Oftsh<T>::get_term(int i) const" << endl;
    cout << "- Oftsh<T>::get_ptr_first_coef() const" << endl;
    cout << "- Oftsh<T>::get_coef(int i) const" << endl;
    cout << "- Oftsh<T>::get_order() const" << endl;
    cout << "- Oftsh<T>::get_nvar() const" << endl;
    cout << "- Oftsh<T>::zero() and its specializations." << endl;

    cout << "---------------------------------------------------" << endl;
    cout << "The following public functions are tested:         " << endl;
    cout << "---------------------------------------------------" << endl;
    cout << "- Oftsh<T>& Oftsh<T>::conjugate() " << endl;
    cout << "- Oftsh<T>& Oftsh<T>::oftsh_smult_t(Oftsh<T> const& a, T const& c)" << endl;
    cout << "- Oftsh<T>& Oftsh<T>::oftsh_mult_t(Oftsh<T> const& a, T const& c)" << endl;
    cout << "- Oftsh< Ofs<U> >& Oftsh<T>::oftsh_smult_tu(Oftsh< Ofs<U> > const& a, Ofs<U> const& ra, U const& c)" << endl;
    cout << "- Oftsh< Ofs<U> >& Oftsh<T>::oftsh_smult_u(Oftsh< Ofs<U> > const& a, U const& c)" << endl;
    cout << "- Oftsh< Ofs<U> >& Oftsh<T>::oftsh_mult_u(Oftsh< Ofs<U> > const& a, U const& c)" << endl;
    cout << "- Oftsh<T>& Oftsh<T>::oftsh_sprod(Oftsh<T> const& a, Oftsh<T> const& b)" << endl;
    cout << "- Oftsh<T>::oftsh_smprod_u(Oftsh< Ofs<U> > const& a, Oftsh< Ofs<U> > const& b, U const& m)" << endl;
    cout << "- Oftsh<T>& Oftsh<T>::derh(Oftsh< T > const& a, int ni)" << endl;
    cout << "- Oftsh<T>& Oftsh<T>::sderh(Oftsh< T > const& a, int ni)" << endl;
    cout << "- Oftsh<T>& Oftsh<T>::dot(Oftsh<T> const& a, double const&  n)" << endl;

    cout <<  setw(5) << setprecision(1) << setiosflags(ios::scientific);
    cout << "---------------------------------------------------" << endl;
    cout << " 1. Tests are made with polynomials of order " << OFTS_ORDER << "." << endl;
    cout << " 2. Tests are made with number of variables: " << REDUCED_NV << "." << endl;
    cout << " 3. Each coefficient is a series of order " << OFS_ORDER << "." << endl;
    cout << " 4. The series are evaluated at arbitrary coordinates." << endl;
    cout << " 5. Random coefficients are set as input into every Ofs coefficient. " << endl;
    cout << "However, an arbitrary decreasing of the coefficients is set: " << endl;
    cout << " around " << (double)rand()/(10.0*(pow(0,7.0)+1)*RAND_MAX) << " @ order 0." << endl;
    cout << " around " << (double)rand()/(10.0*(pow(OFS_ORDER,7.0)+1)*RAND_MAX) << " @ order "<< OFS_ORDER << "." << endl;
    cout << " 6. Whenever it is possible, expected error is displayed. " << endl;
    cout << "---------------------------------------------------" << endl;
    cout <<  setw(5) << setprecision(15) << std::showpos  <<  setiosflags(ios::scientific);
    cout << "Press Enter to proceed with the tests." << endl;
    scanf("%c",&ch);


    //------------------------------------------------------------------------------------
    // Tests.
    //------------------------------------------------------------------------------------
    oftsh_test_conjugate();
    cout << "Press Enter to proceed with the tests." << endl;
    scanf("%c",&ch);
    cout << endl;

    oftsh_test_smult_t();
    cout << "Press Enter to proceed with the tests." << endl;
    scanf("%c",&ch);
    cout << endl;

    oftsh_test_smult_tt();
    cout << "Press Enter to proceed with the tests." << endl;
    scanf("%c",&ch);
    cout << endl;

    oftsh_test_smult_u();
    cout << "Press Enter to proceed with the tests." << endl;
    scanf("%c",&ch);
    cout << endl;

    oftsh_test_sprod();
    cout << "Press Enter to proceed with the tests." << endl;
    scanf("%c",&ch);
    cout << endl;

    oftsh_test_smprod();
    cout << "Press Enter to proceed with the tests." << endl;
    scanf("%c",&ch);
    cout << endl;

    oftsh_test_sderh();
    cout << "Press Enter to proceed with the tests." << endl;
    scanf("%c",&ch);
    cout << endl;

    oftsh_test_derh();
    cout << "Press Enter to proceed with the tests." << endl;
    scanf("%c",&ch);
    cout << endl;

    oftsh_test_dot();
    cout << "Press Enter to proceed with the tests." << endl;
    scanf("%c",&ch);
    cout << endl;



    cout << "---------------------------------------------------" << endl;
    cout << "End of test session.                               " << endl;
    cout << "---------------------------------------------------" << endl;
}

/**
 * \fn void oftsh_test_conjugate()
 * \brief Test of the routine: Oftsh<T>::conjugate() and its specializations.
 */
void oftsh_test_conjugate()
{
    cout << "---------------------------------------------------" << endl;
    cout << "    Test of the routine: Oftsh<T>::conjugate().    " << endl;
    cout << "---------------------------------------------------" << endl;

    //Initialization of OFTS objects, used to initialize the OFTSH objects
    //--------------
    Oftsd xd_ofts(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc_ofts(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);

    //Initialization
    //--------------
    Oftsh< Ofsd > *xd = xd_ofts.get_term(OFTS_ORDER);
    Oftsh< Ofsc > *xdc = xdc_ofts.get_term(OFTS_ORDER);

    double Xvar[REDUCED_NV];
    cdouble Xvard[REDUCED_NV];

    Ofsd ofs(OFS_ORDER);
    Ofsc ofsd(OFS_ORDER);

    //result storage
    cdouble res1, res2;
    cdouble resd1, resd2;
    //array of coordinates (randomly init)
    for(int i = 0; i < REDUCED_NV; i++) Xvar[i] = 1.0*rand()/RAND_MAX;
    for(int i = 0; i < REDUCED_NV; i++) Xvard[i] = 1.0*rand()/RAND_MAX + I*rand()/RAND_MAX ;

    //Set random coefs in xd and xdc
    xd->set_random_coefs();
    xdc->set_random_coefs();

    //Set results
    //--------------
    xd->evaluate(Xvar, ofs);
    xdc->evaluate_conjugate(Xvard, ofsd);

    res1  = conj(ofs.evaluate(VAR));
    resd1 = conj(ofsd.evaluate(VAR));

    //Take conjugate
    //--------------
    xd->conjugate();
    xdc->conjugate();

    //Set results
    //--------------
    xd->evaluate(Xvar, ofs);
    xdc->evaluate(Xvard, ofsd);

    res2  = ofs.evaluate(VAR);
    resd2 = ofsd.evaluate(VAR);

    //Comparison. double case
    //--------------
    cout << "1. double case.         " << endl;
    cout << "Expected result: " << creal(res1) << cimag(res1) << "i" << endl;
    cout << "Obtained result: " << creal(res2) << cimag(res2) << "i" << endl;
    res2 -= res1;
    cout << "Delta: " << cabs(res2) << endl;

    //Comparison. cdouble case
    //--------------
    cout << "2. cdouble case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;

}

/**
 * \fn void oftsh_test_smult_t()
 * \brief Test of the routine: Oftsh<T>& Oftsh<T>::oftsh_smult_t(Oftsh<T> const& a, T const& c).
 */
void oftsh_test_smult_t()
{
    cout << "-----------------------------------------------------------" << endl;
    cout << " Test of the routines:                                     " << endl;
    cout << " 1. Oftsh<T>::oftsh_smult_t(Oftsh<T> const& a, T const& c) " << endl;
    cout << " 2. Oftsh<T>::oftsh_mult_t(Oftsh<T> const& a, T const& c)  " << endl;
    cout << " 3. Oftsh<T>::oftsh_smult_tu(Oftsh< Ofs<U> > const& a, ...)" << endl;
    cout << "-----------------------------------------------------------" << endl;

    //Initialization of OFTS objects, used to initialize the OFTSH objects
    //--------------
    Oftsd xd_ofts(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc_ofts(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsd xd2_ofts(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc2_ofts(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);

    //Initialization
    //--------------
    Oftsh< Ofsd > *xd = xd_ofts.get_term(OFTS_ORDER);
    Oftsh< Ofsc > *xdc = xdc_ofts.get_term(OFTS_ORDER);
    Oftsh< Ofsd > *xd2 = xd2_ofts.get_term(OFTS_ORDER);
    Oftsh< Ofsc > *xdc2 = xdc2_ofts.get_term(OFTS_ORDER);

    double Xvar[REDUCED_NV];
    cdouble Xvard[REDUCED_NV];


    Ofsd ofs(OFS_ORDER);
    Ofsc ofsd(OFS_ORDER);
    Ofsd ofs2(OFS_ORDER);
    Ofsc ofsd2(OFS_ORDER);
    Ofsd ofs3(OFS_ORDER);
    Ofsc ofsd3(OFS_ORDER);

    //result storage
    cdouble res1, res2;
    cdouble resd1, resd2;
    //array of coordinates (randomly init)
    for(int i = 0; i < REDUCED_NV; i++) Xvar[i] = 1.0*rand()/RAND_MAX;
    for(int i = 0; i < REDUCED_NV; i++) Xvard[i] = 1.0*rand()/RAND_MAX + I*rand()/RAND_MAX ;

    //Set random coefs in xd and xdc
    xd->set_random_coefs();
    xdc->set_random_coefs();
    xd2->set_random_coefs();
    xdc2->set_random_coefs();

    //---------------------------------------------------------
    // 1. Oftsh<T>::oftsh_smult_t(Oftsh<T> const& a, T const& c)
    //---------------------------------------------------------

    cout << "---------------------------------------------------------  " << endl;
    cout << " 1. Oftsh<T>::oftsh_smult_t(Oftsh<T> const& a, T const& c) " << endl;
    cout << "---------------------------------------------------------  " << endl;

    //Set results
    //--------------
    xd->evaluate(Xvar, ofs);
    xdc->evaluate(Xvard, ofsd);

    xd2->evaluate(Xvar, ofs2);
    xdc2->evaluate(Xvard, ofsd2);

    //Only here to set random subcoefficients in ofs3 and ofsd3
    xd->evaluate(Xvar, ofs3);
    xdc->evaluate(Xvard, ofsd3);

    res1  = ofs.evaluate(VAR)  + ofs3.evaluate(VAR)*ofs2.evaluate(VAR);
    resd1 = ofsd.evaluate(VAR) + ofsd3.evaluate(VAR)*ofsd2.evaluate(VAR);

    //Take conjugate
    //--------------
    xd->oftsh_smult_t(*xd2, ofs3);
    xdc->oftsh_smult_t(*xdc2, ofsd3);

    //Set results
    //--------------
    xd->evaluate(Xvar, ofs);
    xdc->evaluate(Xvard, ofsd);

    res2  = ofs.evaluate(VAR);
    resd2 = ofsd.evaluate(VAR);

    //Comparison. double case
    //--------------
    cout << "1. double case.         " << endl;
    cout << "Expected result: " << creal(res1) << cimag(res1) << "i" << endl;
    cout << "Obtained result: " << creal(res2) << cimag(res2) << "i" << endl;
    res2 -= res1;
    cout << "Delta: " << cabs(res2) << endl;
    cout << "Expected error : " << cabs(smult_expct_error(*xd2, ofs3, Xvar, VAR)) << endl;
    cout << "Corrected delta: " << cabs(res2 + smult_expct_error(*xd2, ofs3, Xvar, VAR)) << endl;
    cout << "Higher bound on relative accumulated roundoff errors: ~" << 7.0/2*OFS_ORDER*(OFS_ORDER+1)*Manip::nmon(REDUCED_NV, OFTS_ORDER)*1.1e-16 << endl;

    //Comparison. cdouble case
    //--------------
    cout << "2. cdouble case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;
    cout << "Expected error : " << cabs(smult_expct_error(*xdc2, ofsd3, Xvard, VAR)) << endl;
    cout << "Corrected delta: " << cabs(resd2 + smult_expct_error(*xdc2, ofsd3, Xvard, VAR)) << endl;
    cout << "Higher bound on relative accumulated roundoff errors: ~" << 7.0/2*OFS_ORDER*(OFS_ORDER+1)*Manip::nmon(REDUCED_NV, OFTS_ORDER)*1.1e-16 << endl;

    //---------------------------------------------------------
    // 2. Oftsh<T>::oftsh_mult_t(Oftsh<T> const& a, T const& c)
    //---------------------------------------------------------

    //Set random coefs in xd and xdc
    xd->set_random_coefs();
    xdc->set_random_coefs();
    xd2->set_random_coefs();
    xdc2->set_random_coefs();

    cout << "---------------------------------------------------------" << endl;
    cout << " 2. Oftsh<T>::oftsh_mult_t(Oftsh<T> const& a, T const& c)" << endl;
    cout << "---------------------------------------------------------" << endl;

    //Set results
    //--------------
    xd->evaluate(Xvar, ofs);
    xdc->evaluate(Xvard, ofsd);

    xd2->evaluate(Xvar, ofs2);
    xdc2->evaluate(Xvard, ofsd2);

    //Only here to set random subcoefficients in ofs3 and ofsd3
    xd->evaluate(Xvar, ofs3);
    xdc->evaluate(Xvard, ofsd3);

    res1  = ofs3.evaluate(VAR)*ofs2.evaluate(VAR);
    resd1 = ofsd3.evaluate(VAR)*ofsd2.evaluate(VAR);

    //Take conjugate
    //--------------
    xd->oftsh_mult_t(*xd2, ofs3);
    xdc->oftsh_mult_t(*xdc2, ofsd3);

    //Set results
    //--------------
    xd->evaluate(Xvar, ofs);
    xdc->evaluate(Xvard, ofsd);

    res2  = ofs.evaluate(VAR);
    resd2 = ofsd.evaluate(VAR);

    //Comparison. double case
    //--------------
    cout << "1. double case.         " << endl;
    cout << "Expected result: " << creal(res1) << cimag(res1) << "i" << endl;
    cout << "Obtained result: " << creal(res2) << cimag(res2) << "i" << endl;
    res2 -= res1;
    cout << "Delta: " << cabs(res2) << endl;
    cout << "Expected error : " << cabs(smult_expct_error(*xd2, ofs3, Xvar, VAR)) << endl;
    cout << "Corrected delta: " << cabs(res2 + smult_expct_error(*xd2, ofs3, Xvar, VAR)) << endl;
    cout << "Higher bound on relative accumulated roundoff errors: ~" << 7.0/2*OFS_ORDER*(OFS_ORDER+1)*Manip::nmon(REDUCED_NV, OFTS_ORDER)*1.1e-16 << endl;

    //Comparison. cdouble case
    //--------------
    cout << "2. cdouble case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;
    cout << "Expected error : " << cabs(smult_expct_error(*xdc2, ofsd3, Xvard, VAR)) << endl;
    cout << "Corrected delta: " << cabs(resd2 + smult_expct_error(*xdc2, ofsd3, Xvard, VAR)) << endl;
    cout << "Higher bound on relative accumulated roundoff errors: ~" << 7.0/2*OFS_ORDER*(OFS_ORDER+1)*Manip::nmon(REDUCED_NV, OFTS_ORDER)*1.1e-16 << endl;


    //---------------------------------------------------------
    // 3. Oftsh<T>::oftsh_smult_tu(Oftsh< Ofs<U> > const& a, Ofs<U> const& ra, U const& m)
    //---------------------------------------------------------
    double c = 2.0;
    cdouble cd = 2.0 + I;

    //Set random coefs in xd and xdc
    xd->set_random_coefs();
    xdc->set_random_coefs();
    xd2->set_random_coefs();
    xdc2->set_random_coefs();

    cout << "---------------------------------------------------------" << endl;
    cout << " 3. Oftsh<T>::oftsh_smult_tu(Oftsh< Ofs<U> > const& a, Ofs<U> const& ra, U const& m)" << endl;
    cout << "---------------------------------------------------------" << endl;

    //Set results
    //--------------
    xd->evaluate(Xvar, ofs);
    xdc->evaluate(Xvard, ofsd);

    xd2->evaluate(Xvar, ofs2);
    xdc2->evaluate(Xvard, ofsd2);

    //Only here to set random subcoefficients in ofs3 and ofsd3
    xd->evaluate(Xvar, ofs3);
    xdc->evaluate(Xvard, ofsd3);

    res1  = ofs.evaluate(VAR)  + c*ofs3.evaluate(VAR)*ofs2.evaluate(VAR);
    resd1 = ofsd.evaluate(VAR) + cd*ofsd3.evaluate(VAR)*ofsd2.evaluate(VAR);

    //Take conjugate
    //--------------
    xd->oftsh_smult_tu(*xd2, ofs3, c);
    xdc->oftsh_smult_tu(*xdc2, ofsd3, cd);

    //Set results
    //--------------
    xd->evaluate(Xvar, ofs);
    xdc->evaluate(Xvard, ofsd);

    res2  = ofs.evaluate(VAR);
    resd2 = ofsd.evaluate(VAR);

    //Comparison. double case
    //--------------
    cout << "1. double case.         " << endl;
    cout << "Expected result: " << creal(res1) << cimag(res1) << "i" << endl;
    cout << "Obtained result: " << creal(res2) << cimag(res2) << "i" << endl;
    res2 -= res1;
    cout << "Delta: " << cabs(res2) << endl;
    cout << "Expected error : " << cabs(c*smult_expct_error(*xd2, ofs3, Xvar, VAR)) << endl;
    cout << "Corrected delta: " << cabs(res2 + c*smult_expct_error(*xd2, ofs3, Xvar, VAR)) << endl;
    cout << "Higher bound on relative accumulated roundoff errors: ~" << 7.0/2*OFS_ORDER*(OFS_ORDER+1)*Manip::nmon(REDUCED_NV, OFTS_ORDER)*1.1e-16 << endl;

    //Comparison. cdouble case
    //--------------
    cout << "2. cdouble case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;
    cout << "Expected error : " << cabs(cd*smult_expct_error(*xdc2, ofsd3, Xvard, VAR)) << endl;
    cout << "Corrected delta: " << cabs(resd2 + cd*smult_expct_error(*xdc2, ofsd3, Xvard, VAR)) << endl;
    cout << "Higher bound on relative accumulated roundoff errors: ~" << 7.0/2*OFS_ORDER*(OFS_ORDER+1)*Manip::nmon(REDUCED_NV, OFTS_ORDER)*1.1e-16 << endl;
}

/**
 * \fn void oftsh_test_smult_tt()
 * \brief Test of the routine: Oftsh<T>& oftsh_smult_tt(Oftsh<T> const& a, T const& c1, T const& c2).
 */
void oftsh_test_smult_tt()
{
    cout << "-----------------------------------------------------------" << endl;
    cout << " Test of the routines:                                     " << endl;
    cout << " oftsh_smult_tt(Oftsh<T> const& a, T const& c1, T const& c2) " << endl;
    cout << "-----------------------------------------------------------" << endl;

    //Initialization of OFTS objects, used to initialize the OFTSH objects
    //--------------
    Oftsd xd_ofts(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc_ofts(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsd xd2_ofts(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc2_ofts(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);

    //Initialization
    //--------------
    Oftsh< Ofsd > *xd = xd_ofts.get_term(OFTS_ORDER);
    Oftsh< Ofsc > *xdc = xdc_ofts.get_term(OFTS_ORDER);
    Oftsh< Ofsd > *xd2 = xd2_ofts.get_term(OFTS_ORDER);
    Oftsh< Ofsc > *xdc2 = xdc2_ofts.get_term(OFTS_ORDER);

    double Xvar[REDUCED_NV];
    cdouble Xvard[REDUCED_NV];


    Ofsd ofs(OFS_ORDER);
    Ofsc ofsd(OFS_ORDER);
    Ofsd ofs2(OFS_ORDER);
    Ofsc ofsd2(OFS_ORDER);
    Ofsd ofs3(OFS_ORDER);
    Ofsc ofsd3(OFS_ORDER);
    Ofsd ofs4(OFS_ORDER);
    Ofsc ofsd4(OFS_ORDER);

    Ofsd temp(OFS_ORDER);
    Ofsc tempc(OFS_ORDER);

    //result storage
    cdouble res1, res2;
    cdouble resd1, resd2;
    //array of coordinates (randomly init)
    for(int i = 0; i < REDUCED_NV; i++) Xvar[i] = 1.0*rand()/RAND_MAX;
    for(int i = 0; i < REDUCED_NV; i++) Xvard[i] = 1.0*rand()/RAND_MAX + I*rand()/RAND_MAX ;

    //Set random coefs in xd and xdc
    xd->set_random_coefs();
    xdc->set_random_coefs();
    xd2->set_random_coefs();
    xdc2->set_random_coefs();

    //Set results
    //--------------
    xd->evaluate(Xvar, ofs);
    xdc->evaluate(Xvard, ofsd);

    xd2->evaluate(Xvar, ofs2);
    xdc2->evaluate(Xvard, ofsd2);

    //Only here to set random subcoefficients in ofs3/4 and ofsd3/4
    xd->evaluate(Xvar, ofs3);
    xdc->evaluate(Xvard, ofsd3);
    xd2->evaluate(Xvar, ofs4);
    xdc2->evaluate(Xvard, ofsd4);

    res1  = ofs.evaluate(VAR)  + ofs4.evaluate(VAR)*ofs3.evaluate(VAR)*ofs2.evaluate(VAR);
    resd1 = ofsd.evaluate(VAR) + ofsd4.evaluate(VAR)*ofsd3.evaluate(VAR)*ofsd2.evaluate(VAR);

    //Take conjugate
    //--------------
    xd->oftsh_smult_tt(*xd2, ofs3, ofs4, temp);
    xdc->oftsh_smult_tt(*xdc2, ofsd3, ofsd4, tempc);

    //Set results
    //--------------
    xd->evaluate(Xvar, ofs);
    xdc->evaluate(Xvard, ofsd);

    res2  = ofs.evaluate(VAR);
    resd2 = ofsd.evaluate(VAR);

    //Comparison. double case
    //--------------
    cout << "1. double case.         " << endl;
    cout << "Expected result: " << creal(res1) << cimag(res1) << "i" << endl;
    cout << "Obtained result: " << creal(res2) << cimag(res2) << "i" << endl;
    res2 -= res1;
    cout << "Delta: " << cabs(res2) << endl;
    cout << "Note: the expected error is trickier to obtain since we do two products in a row in ofs_smprod_t " << endl;

    //Comparison. cdouble case
    //--------------
    cout << "2. cdouble case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;
    cout << "Note: the expected error is trickier to obtain since we do two products in a row in ofs_smprod_t " << endl;
}

/**
 * \fn void oftsh_test_smult_u()
 * \brief Test of the routine: Oftsh<T>& Oftsh<T>::oftsh_smult_t(Oftsh<T> const& a, T const& c).
 */
void oftsh_test_smult_u()
{
    cout << "------------------------------------------------------------------" << endl;
    cout << " Test of the routines:                                            " << endl;
    cout << " 1. Oftsh<T>::oftsh_smult_u(Oftsh< Ofs<U> > const& a, U const& c) " << endl;
    cout << " 2. Oftsh<T>::oftsh_mult_u(Oftsh< Ofs<U> > const& a, U const& c)  " << endl;
    cout << "------------------------------------------------------------------" << endl;

    //Initialization of OFTS objects, used to initialize the OFTSH objects
    //--------------
    Oftsd xd_ofts(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc_ofts(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsd xd2_ofts(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc2_ofts(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);

    //Initialization
    //--------------
    Oftsh< Ofsd > *xd = xd_ofts.get_term(OFTS_ORDER);
    Oftsh< Ofsc > *xdc = xdc_ofts.get_term(OFTS_ORDER);
    Oftsh< Ofsd > *xd2 = xd2_ofts.get_term(OFTS_ORDER);
    Oftsh< Ofsc > *xdc2 = xdc2_ofts.get_term(OFTS_ORDER);

    double Xvar[REDUCED_NV];
    cdouble Xvard[REDUCED_NV];


    Ofsd ofs(OFS_ORDER);
    Ofsc ofsd(OFS_ORDER);
    Ofsd ofs2(OFS_ORDER);
    Ofsc ofsd2(OFS_ORDER);
    Ofsd ofs3(OFS_ORDER);
    Ofsc ofsd3(OFS_ORDER);

    double c = 2.0;
    cdouble cd = 2.0 + I;

    //result storage
    cdouble res1, res2;
    cdouble resd1, resd2;
    //array of coordinates (randomly init)
    for(int i = 0; i < REDUCED_NV; i++) Xvar[i] = 1.0*rand()/RAND_MAX;
    for(int i = 0; i < REDUCED_NV; i++) Xvard[i] = 1.0*rand()/RAND_MAX + I*rand()/RAND_MAX ;

    //Set random coefs in xd and xdc
    xd->set_random_coefs();
    xdc->set_random_coefs();
    xd2->set_random_coefs();
    xdc2->set_random_coefs();

    //---------------------------------------------------------
    // 1. Oftsh<T>::oftsh_smult_t(Oftsh<T> const& a, T const& c)
    //---------------------------------------------------------
    cout << "------------------------------------------------------------------" << endl;
    cout << " 1. Oftsh<T>::oftsh_smult_u(Oftsh< Ofs<U> > const& a, U const& c) " << endl;
    cout << "------------------------------------------------------------------" << endl;

    //Set results
    //--------------
    xd->evaluate(Xvar, ofs);
    xdc->evaluate(Xvard, ofsd);

    xd2->evaluate(Xvar, ofs2);
    xdc2->evaluate(Xvard, ofsd2);

    //Only here to set random subcoefficients in ofs3 and ofsd3
    xd->evaluate(Xvar, ofs3);
    xdc->evaluate(Xvard, ofsd3);

    res1  = ofs.evaluate(VAR)  + c*ofs2.evaluate(VAR);
    resd1 = ofsd.evaluate(VAR) + cd*ofsd2.evaluate(VAR);

    //Take conjugate
    //--------------
    xd->oftsh_smult_u(*xd2, c);
    xdc->oftsh_smult_u(*xdc2, cd);

    //Set results
    //--------------
    xd->evaluate(Xvar, ofs);
    xdc->evaluate(Xvard, ofsd);

    res2  = ofs.evaluate(VAR);
    resd2 = ofsd.evaluate(VAR);

    //Comparison. double case
    //--------------
    cout << "1. double case.         " << endl;
    cout << "Expected result: " << creal(res1) << cimag(res1) << "i" << endl;
    cout << "Obtained result: " << creal(res2) << cimag(res2) << "i" << endl;
    res2 -= res1;
    cout << "Delta: " << cabs(res2) << endl;

    //Comparison. cdouble case
    //--------------
    cout << "2. cdouble case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;

    //---------------------------------------------------------
    // 2. Oftsh<T>::oftsh_mult_t(Oftsh<T> const& a, T const& c)
    //---------------------------------------------------------

    //Set random coefs in xd and xdc
    xd->set_random_coefs();
    xdc->set_random_coefs();
    xd2->set_random_coefs();
    xdc2->set_random_coefs();

    cout << "------------------------------------------------------------------" << endl;
    cout << " 2. Oftsh<T>::oftsh_mult_u(Oftsh< Ofs<U> > const& a, U const& c)  " << endl;
    cout << "------------------------------------------------------------------" << endl;

    //Set results
    //--------------
    xd->evaluate(Xvar, ofs);
    xdc->evaluate(Xvard, ofsd);

    xd2->evaluate(Xvar, ofs2);
    xdc2->evaluate(Xvard, ofsd2);

    //Only here to set random subcoefficients in ofs3 and ofsd3
    xd->evaluate(Xvar, ofs3);
    xdc->evaluate(Xvard, ofsd3);

    res1  = c*ofs2.evaluate(VAR);
    resd1 = cd*ofsd2.evaluate(VAR);

    //Take conjugate
    //--------------
    xd->oftsh_mult_u(*xd2, c);
    xdc->oftsh_mult_u(*xdc2, cd);

    //Set results
    //--------------
    xd->evaluate(Xvar, ofs);
    xdc->evaluate(Xvard, ofsd);

    res2  = ofs.evaluate(VAR);
    resd2 = ofsd.evaluate(VAR);

    //Comparison. double case
    //--------------
    cout << "1. double case.         " << endl;
    cout << "Expected result: " << creal(res1) << cimag(res1) << "i" << endl;
    cout << "Obtained result: " << creal(res2) << cimag(res2) << "i" << endl;
    res2 -= res1;
    cout << "Delta: " << cabs(res2) << endl;

    //Comparison. cdouble case
    //--------------
    cout << "2. cdouble case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;

}

/**
 * \fn void oftsh_test_sprod()
 * \brief Test of the routine: Oftsh<T>::oftsh_sprod(Oftsh<T> const& a, Oftsh<T> const& b) and its specializations.
 */
void oftsh_test_sprod()
{
    cout << "-----------------------------------------------------------------" << endl;
    cout << " Test of the routines:                                           " << endl;
    cout << " 1. Oftsh<T>::oftsh_sprod(Oftsh<T> const& a, Oftsh<T> const& b)  " << endl;
    cout << " 2. its specializations.                                         " << endl;
    cout << "-----------------------------------------------------------------" << endl;

    tic();
    //Initialization of OFTS objects, used to initialize the OFTSH objects
    //--------------
    Oftsd xd_ofts(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc_ofts(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsd xd2_ofts(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc2_ofts(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsd xd3_ofts(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc3_ofts(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);

    //Initialization
    //--------------
    Oftsh< Ofsd > *xd = xd_ofts.get_term(OFTS_ORDER);
    Oftsh< Ofsc > *xdc = xdc_ofts.get_term(OFTS_ORDER);
    Oftsh< Ofsd > *xd2 = xd2_ofts.get_term(OFTS_ORDER/2);
    Oftsh< Ofsc > *xdc2 = xdc2_ofts.get_term(OFTS_ORDER/2);
    Oftsh< Ofsd > *xd3 = xd3_ofts.get_term(OFTS_ORDER/2);
    Oftsh< Ofsc > *xdc3 = xdc3_ofts.get_term(OFTS_ORDER/2);

    cout << "Initialization in : " << toc() << "s." << endl;


    double Xvar[REDUCED_NV];
    cdouble Xvard[REDUCED_NV];

    Ofsd ofs(OFS_ORDER);
    Ofsc ofsd(OFS_ORDER);
    Ofsd ofs2(OFS_ORDER);
    Ofsc ofsd2(OFS_ORDER);
    Ofsd ofs3(OFS_ORDER);
    Ofsc ofsd3(OFS_ORDER);

    //result storage
    cdouble res1, res2;
    cdouble resd1, resd2;
    //array of coordinates (randomly init)
    for(int i = 0; i < REDUCED_NV; i++) Xvar[i] = 1.0*rand()/RAND_MAX;
    for(int i = 0; i < REDUCED_NV; i++) Xvard[i] = 1.0*rand()/RAND_MAX + I*rand()/RAND_MAX ;



    //Set random coefs in xd and xdc
    xd->set_random_coefs();
    xdc->set_random_coefs();
    xd2->set_random_coefs();
    xdc2->set_random_coefs();
    xd3->set_random_coefs();
    xdc3->set_random_coefs();

    //Set results
    //--------------
    xd->evaluate(Xvar, ofs);
    xdc->evaluate(Xvard, ofsd);
    xd2->evaluate(Xvar, ofs2);
    xdc2->evaluate(Xvard, ofsd2);
    xd3->evaluate(Xvar, ofs3);
    xdc3->evaluate(Xvard, ofsd3);


    res1  = ofs.evaluate(VAR)  + ofs3.evaluate(VAR)*ofs2.evaluate(VAR);
    resd1 = ofsd.evaluate(VAR) + ofsd3.evaluate(VAR)*ofsd2.evaluate(VAR);

    //Take conjugate
    //--------------
    tic();
    xd->oftsh_sprod(*xd2, *xd3);
    cout << "Product in : " << toc() << "s." << endl;
    tic();
    xdc->oftsh_sprod(*xdc2, *xdc3);
    cout << "Product in : " << toc() << "s." << endl;


    //Set results
    //--------------
    xd->evaluate(Xvar, ofs);
    xdc->evaluate(Xvard, ofsd);

    res2  = ofs.evaluate(VAR);
    resd2 = ofsd.evaluate(VAR);

    //Comparison. double case
    //--------------
    double maxC = fabs(pow(Xvar[0], 1.0*OFTS_ORDER));
    for(int i = 1; i < REDUCED_NV ; i++) if(fabs(pow(Xvar[i], 1.0*OFTS_ORDER)) > maxC) maxC = fabs(pow(Xvar[i], 1.0*OFTS_ORDER));

    cout << "1. double case.         " << endl;
    cout << "Expected result: " << creal(res1) << cimag(res1) << "i" << endl;
    cout << "Obtained result: " << creal(res2) << cimag(res2) << "i" << endl;
    res2 -= res1;
    cout << "Delta: " << cabs(res2) << endl;
    //How do we estimate the error made on the product? Compute the number of single product in Oftsh*Oftsh at given order/nv = nmon(nv, order/2)^2
    //Then an estimation of this error is \pm Ofs*Ofs error times this number times xvar(0)^order.
    cout <<  setw(5) << setprecision(2) << std::showpos  <<  setiosflags(ios::scientific);
    cout << "Expected error : up to ";
    cout << maxC*cabs(sprod_expct_error(xd2->get_coef(0), xd3->get_coef(0), VAR))*pow(Manip::nmon(REDUCED_NV, OFTS_ORDER/2), 2.0) << endl;
    cout <<  setw(5) << setprecision(15) << std::showpos  <<  setiosflags(ios::scientific);

    //Comparison. cdouble case
    //--------------
    maxC = cabs(cpow(Xvard[0], 1.0*OFTS_ORDER+0.0*I));
    for(int i = 1; i < REDUCED_NV ; i++) if(cabs(cpow(Xvard[i], 1.0*OFTS_ORDER+0.0*I)) > maxC) maxC = cabs(cpow(Xvard[i], 1.0*OFTS_ORDER+0.0*I));

    cout << "2. cdouble case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;
    cout <<  setw(5) << setprecision(2) << std::showpos  <<  setiosflags(ios::scientific);
    cout << "Expected error : up to ";
    cout <<  maxC*cabs(sprod_expct_error(xdc2->get_coef(0), xdc3->get_coef(0), VAR))*pow(Manip::nmon(REDUCED_NV, OFTS_ORDER/2), 2.0) << endl;
    cout <<  setw(5) << setprecision(15) << std::showpos  <<  setiosflags(ios::scientific);
}

/**
 * \fn void oftsh_test_smprod()
 * \brief Test of the routine: Oftsh<T>::oftsh_smprod_u(Oftsh< Ofs<U> > const& a, Oftsh< Ofs<U> > const& b, U const& m).
 */
void oftsh_test_smprod()
{
    cout << "-----------------------------------------------------------------" << endl;
    cout << " Test of the routines:                                           " << endl;
    cout << " 1. Oftsh<T>::oftsh_smprod_u(Oftsh< Ofs<U> > const& a, Oftsh< Ofs<U> > const& b, U const& m)  " << endl;
    cout << " 2. its specializations.                                         " << endl;
    cout << "-----------------------------------------------------------------" << endl;

    tic();
    //Initialization of OFTS objects, used to initialize the OFTSH objects
    //--------------
    Oftsd xd_ofts(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc_ofts(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsd xd2_ofts(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc2_ofts(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsd xd3_ofts(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc3_ofts(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);

    //Initialization
    //--------------
    Oftsh< Ofsd > *xd = xd_ofts.get_term(OFTS_ORDER);
    Oftsh< Ofsc > *xdc = xdc_ofts.get_term(OFTS_ORDER);
    Oftsh< Ofsd > *xd2 = xd2_ofts.get_term(OFTS_ORDER/2);
    Oftsh< Ofsc > *xdc2 = xdc2_ofts.get_term(OFTS_ORDER/2);
    Oftsh< Ofsd > *xd3 = xd3_ofts.get_term(OFTS_ORDER/2);
    Oftsh< Ofsc > *xdc3 = xdc3_ofts.get_term(OFTS_ORDER/2);

    cout << "Initialization in : " << toc() << "s." << endl;

    double Xvar[REDUCED_NV];
    cdouble Xvard[REDUCED_NV];

    double c = 2.0;
    cdouble cd = 2.0+I;

    Ofsd ofs(OFS_ORDER);
    Ofsc ofsd(OFS_ORDER);
    Ofsd ofs2(OFS_ORDER);
    Ofsc ofsd2(OFS_ORDER);
    Ofsd ofs3(OFS_ORDER);
    Ofsc ofsd3(OFS_ORDER);

    //result storage
    cdouble res1, res2;
    cdouble resd1, resd2;
    //array of coordinates (randomly init)
    for(int i = 0; i < REDUCED_NV; i++) Xvar[i] = 1.0*rand()/RAND_MAX;
    for(int i = 0; i < REDUCED_NV; i++) Xvard[i] = 1.0*rand()/RAND_MAX + I*rand()/RAND_MAX ;

    //Set random coefs in xd and xdc
    xd->set_random_coefs();
    xdc->set_random_coefs();
    xd2->set_random_coefs();
    xdc2->set_random_coefs();
    xd3->set_random_coefs();
    xdc3->set_random_coefs();

    //Set results
    //--------------
    xd->evaluate(Xvar, ofs);
    xdc->evaluate(Xvard, ofsd);
    xd2->evaluate(Xvar, ofs2);
    xdc2->evaluate(Xvard, ofsd2);
    xd3->evaluate(Xvar, ofs3);
    xdc3->evaluate(Xvard, ofsd3);


    res1  = ofs.evaluate(VAR)  + c*ofs3.evaluate(VAR)*ofs2.evaluate(VAR);
    resd1 = ofsd.evaluate(VAR) + cd*ofsd3.evaluate(VAR)*ofsd2.evaluate(VAR);

    //Take conjugate
    //--------------
    tic();
    xd->oftsh_smprod_u(*xd2, *xd3, c);
    cout << "Product in : " << toc() << "s." << endl;
    tic();
    xdc->oftsh_smprod_u(*xdc2, *xdc3, cd);
    cout << "Product in : " << toc() << "s." << endl;


    //Set results
    //--------------
    xd->evaluate(Xvar, ofs);
    xdc->evaluate(Xvard, ofsd);

    res2  = ofs.evaluate(VAR);
    resd2 = ofsd.evaluate(VAR);

    //Comparison. double case
    //--------------
    double maxC = fabs(pow(Xvar[0], 1.0*OFTS_ORDER));
    for(int i = 1; i < REDUCED_NV ; i++) if(fabs(pow(Xvar[i], 1.0*OFTS_ORDER)) > maxC) maxC = fabs(pow(Xvar[i], 1.0*OFTS_ORDER));

    cout << "1. double case.         " << endl;
    cout << "Expected result: " << creal(res1) << cimag(res1) << "i" << endl;
    cout << "Obtained result: " << creal(res2) << cimag(res2) << "i" << endl;
    res2 -= res1;
    cout << "Delta: " << cabs(res2) << endl;
    //How do we estimate the error made on the product? Compute the number of single product in Oftsh*Oftsh at given order/nv = nmon(nv, order/2)^2
    //Then an estimation of this error is \pm Ofs*Ofs error times this number times xvar(0)^order.
    cout <<  setw(5) << setprecision(2) << std::showpos  <<  setiosflags(ios::scientific);
    cout << "Expected error : up to ";
    cout <<  maxC*cabs(smprod_expct_error(xd2->get_coef(0), xd3->get_coef(0), c, VAR))*pow(Manip::nmon(REDUCED_NV, OFTS_ORDER/2), 2.0) << endl;
    cout <<  setw(5) << setprecision(15) << std::showpos  <<  setiosflags(ios::scientific);


    //Comparison. cdouble case
    //--------------
    maxC = cabs(cpow(Xvard[0], 1.0*OFTS_ORDER+0.0*I));
    for(int i = 1; i < REDUCED_NV ; i++) if(cabs(cpow(Xvard[i], 1.0*OFTS_ORDER+0.0*I)) > maxC) maxC = cabs(cpow(Xvard[i], 1.0*OFTS_ORDER+0.0*I));

    cout << "2. cdouble case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;
    cout <<  setw(5) << setprecision(2) << std::showpos  <<  setiosflags(ios::scientific);
    cout << "Expected error : up to ";
    cout << maxC*cabs(smprod_expct_error(xdc2->get_coef(0), xdc3->get_coef(0), cd, VAR))*pow(Manip::nmon(REDUCED_NV, OFTS_ORDER/2), 2.0) << endl;
    cout <<  setw(5) << setprecision(15) << std::showpos  <<  setiosflags(ios::scientific);
}

/**
 * \fn void oftsh_test_smprod_t()
 * \brief Test of the routine: Oftsh<T>::oftsh_smprod_u(Oftsh< Ofs<U> > const& a, Oftsh< Ofs<U> > const& b, U const& m).
 */
void oftsh_test_smprod_t()
{
    cout << "-----------------------------------------------------------------" << endl;
    cout << " Test of the routines:                                           " << endl;
    cout << " 1. Oftsh<T>::oftsh_smprod_t(Oftsh< T > const& a, Oftsh< T> const& b, T const& c)  " << endl;
    cout << " 2. its specializations.                                         " << endl;
    cout << "-----------------------------------------------------------------" << endl;

    tic();
    //Initialization of OFTS objects, used to initialize the OFTSH objects
    //--------------
    Oftsd xd_ofts(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc_ofts(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsd xd2_ofts(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc2_ofts(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsd xd3_ofts(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc3_ofts(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);

    //Initialization
    //--------------
    Oftsh< Ofsd > *xd = xd_ofts.get_term(OFTS_ORDER);
    Oftsh< Ofsc > *xdc = xdc_ofts.get_term(OFTS_ORDER);
    Oftsh< Ofsd > *xd2 = xd2_ofts.get_term(OFTS_ORDER/2);
    Oftsh< Ofsc > *xdc2 = xdc2_ofts.get_term(OFTS_ORDER/2);
    Oftsh< Ofsd > *xd3 = xd3_ofts.get_term(OFTS_ORDER/2);
    Oftsh< Ofsc > *xdc3 = xdc3_ofts.get_term(OFTS_ORDER/2);

    Ofsd temp(OFS_ORDER);
    Ofsc tempc(OFS_ORDER);

    cout << "Initialization in : " << toc() << "s." << endl;

    double Xvar[REDUCED_NV];
    cdouble Xvard[REDUCED_NV];

    double c = 2.0;
    cdouble cd = 2.0+I;


    Ofsd ofs(OFS_ORDER);
    Ofsc ofsd(OFS_ORDER);
    Ofsd ofs2(OFS_ORDER);
    Ofsc ofsd2(OFS_ORDER);
    Ofsd ofs3(OFS_ORDER);
    Ofsc ofsd3(OFS_ORDER);
    Ofsd ofs4(OFS_ORDER);
    Ofsc ofsd4(OFS_ORDER);

    //result storage
    cdouble res1, res2;
    cdouble resd1, resd2;
    //array of coordinates (randomly init)
    for(int i = 0; i < REDUCED_NV; i++) Xvar[i] = 1.0*rand()/RAND_MAX;
    for(int i = 0; i < REDUCED_NV; i++) Xvard[i] = 1.0*rand()/RAND_MAX + I*rand()/RAND_MAX ;



    //Set random coefs in xd and xdc
    xd->set_random_coefs();
    xdc->set_random_coefs();
    xd2->set_random_coefs();
    xdc2->set_random_coefs();
    xd3->set_random_coefs();
    xdc3->set_random_coefs();

    //Set results
    //--------------
    xd->evaluate(Xvar, ofs);
    xdc->evaluate(Xvard, ofsd);
    xd2->evaluate(Xvar, ofs2);
    xdc2->evaluate(Xvard, ofsd2);
    xd3->evaluate(Xvar, ofs3);
    xdc3->evaluate(Xvard, ofsd3);
    //Just to set random coefficients
    xd3->evaluate(Xvar, ofs4);
    xdc3->evaluate(Xvard, ofsd4);

    res1  = ofs.evaluate(VAR)  + ofs4.evaluate(VAR)*ofs3.evaluate(VAR)*ofs2.evaluate(VAR);
    resd1 = ofsd.evaluate(VAR) + ofsd4.evaluate(VAR)*ofsd3.evaluate(VAR)*ofsd2.evaluate(VAR);

    //Take conjugate
    //--------------
    tic();
    xd->oftsh_smprod_t(*xd2, *xd3, ofs4, temp);
    cout << "Product in : " << toc() << "s." << endl;
    tic();
    xdc->oftsh_smprod_t(*xdc2, *xdc3, ofsd4, tempc);
    cout << "Product in : " << toc() << "s." << endl;


    //Set results
    //--------------
    xd->evaluate(Xvar, ofs);
    xdc->evaluate(Xvard, ofsd);

    res2  = ofs.evaluate(VAR);
    resd2 = ofsd.evaluate(VAR);

    //Comparison. double case
    //--------------
    double maxC = fabs(pow(Xvar[0], 1.0*OFTS_ORDER));
    for(int i = 1; i < REDUCED_NV ; i++) if(fabs(pow(Xvar[i], 1.0*OFTS_ORDER)) > maxC) maxC = fabs(pow(Xvar[i], 1.0*OFTS_ORDER));

    cout << "1. double case.         " << endl;
    cout << "Expected result: " << creal(res1) << cimag(res1) << "i" << endl;
    cout << "Obtained result: " << creal(res2) << cimag(res2) << "i" << endl;
    res2 -= res1;
    cout << "Delta: " << cabs(res2) << endl;
    //How do we estimate the error made on the product? Compute the number of single product in Oftsh*Oftsh at given order/nv = nmon(nv, order/2)^2
    //Then an estimation of this error is \pm Ofs*Ofs error times this number times xvar(0)^order.
    cout <<  setw(5) << setprecision(2) << std::showpos  <<  setiosflags(ios::scientific);
    cout << "Expected error : up to ";
    cout <<  maxC*cabs(smprod_expct_error(xd2->get_coef(0), xd3->get_coef(0), c, VAR))*pow(Manip::nmon(REDUCED_NV, OFTS_ORDER/2), 2.0) << endl;
    cout <<  setw(5) << setprecision(15) << std::showpos  <<  setiosflags(ios::scientific);


    //Comparison. cdouble case
    //--------------
    maxC = cabs(cpow(Xvard[0], 1.0*OFTS_ORDER+0.0*I));
    for(int i = 1; i < REDUCED_NV ; i++) if(cabs(cpow(Xvard[i], 1.0*OFTS_ORDER+0.0*I)) > maxC) maxC = cabs(cpow(Xvard[i], 1.0*OFTS_ORDER+0.0*I));

    cout << "2. cdouble case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;
    cout <<  setw(5) << setprecision(2) << std::showpos  <<  setiosflags(ios::scientific);
    cout << "Expected error : up to ";
    cout << maxC*cabs(smprod_expct_error(xdc2->get_coef(0), xdc3->get_coef(0), cd, VAR))*pow(Manip::nmon(REDUCED_NV, OFTS_ORDER/2), 2.0) << endl;
    cout <<  setw(5) << setprecision(15) << std::showpos  <<  setiosflags(ios::scientific);
}

/**
 * \fn void oftsh_test_sderh()
 * \brief Test of the routine: Oftsh<T>& Oftsh<T>::sderh(Oftsh< T > const& a, int ni) and its specializations.
 */
void oftsh_test_sderh()
{
    cout << "----------------------------------------------------------" << endl;
    cout << "    Test of the routine:                                  " << endl;
    cout << " Oftsh<T>& Oftsh<T>::sderh(Oftsh< T > const& a, int ni).  " << endl;
    cout << "----------------------------------------------------------" << endl;
    double epsilon = 1e-8;
    cout << " The partial derivatives are estimated with a simple Euler scheme with epsilon = " << epsilon << "." << endl;

    //Initialization of OFTS objects, used to initialize the OFTSH objects
    //--------------
    Oftsd xd_ofts(REDUCED_NV, OFTS_ORDER-1, 1, OFS_ORDER);
    Oftsc xdc_ofts(REDUCED_NV, OFTS_ORDER-1, 1, OFS_ORDER);
    Oftsd xd2_ofts(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc2_ofts(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);


    //Initialization
    //--------------
    Oftsh< Ofsd > *xd = xd_ofts.get_term(OFTS_ORDER-1);
    Oftsh< Ofsc > *xdc = xdc_ofts.get_term(OFTS_ORDER-1);
    Oftsh< Ofsd > *xd2 = xd2_ofts.get_term(OFTS_ORDER);
    Oftsh< Ofsc > *xdc2 = xdc2_ofts.get_term(OFTS_ORDER);

    double Xvar[REDUCED_NV];
    double Xvar2[REDUCED_NV];
    cdouble Xvard[REDUCED_NV];
    cdouble Xvard2[REDUCED_NV];

    Ofsd ofs(OFS_ORDER);
    Ofsc ofsd(OFS_ORDER);
    Ofsd ofs2(OFS_ORDER);
    Ofsc ofsd2(OFS_ORDER);
    Ofsd ofs3(OFS_ORDER);
    Ofsc ofsd3(OFS_ORDER);

    //result storage
    cdouble res1, res2;
    cdouble resd1, resd2;

    //Loop on all the variables
    for(int ni = 1; ni<= REDUCED_NV ; ni++)
    {
        cout << "---------------------------" << endl;
        cout << "    Variable " << ni << ": " << endl;
        cout << "---------------------------" << endl;
        //array of coordinates (randomly init)
        for(int i = 0; i < REDUCED_NV; i++)
        {
            Xvar[i] = 1.0*rand()/RAND_MAX;
            Xvar2[i] = Xvar[i];
        }
        Xvar2[ni-1] += epsilon;

        for(int i = 0; i < REDUCED_NV; i++)
        {
            Xvard[i] = 1.0*rand()/RAND_MAX + I*rand()/RAND_MAX ;
            Xvard2[i] = Xvard[i];
        }
        Xvard2[ni-1] += epsilon;

        //Set random coefs in xd and xdc
        xd->set_random_coefs();
        xdc->set_random_coefs();
        xd2->set_random_coefs();
        xdc2->set_random_coefs();

        //Set results
        //--------------
        xd2->evaluate(Xvar, ofs);
        xdc2->evaluate(Xvard, ofsd);
        xd2->evaluate(Xvar2, ofs2);
        xdc2->evaluate(Xvard2, ofsd2);

        xd->evaluate(Xvar2, ofs3);
        xdc->evaluate(Xvard2, ofsd3);

        res1  = ofs3.evaluate(VAR) + (ofs2.evaluate(VAR) - ofs.evaluate(VAR))/epsilon;
        resd1 = ofsd3.evaluate(VAR) + (ofsd2.evaluate(VAR) - ofsd.evaluate(VAR))/epsilon;

        //Take sderh
        //--------------
        xd->sderh(*xd2, ni);
        xdc->sderh(*xdc2, ni);

        //Set results
        //--------------
        xd->evaluate(Xvar, ofs);
        xdc->evaluate(Xvard, ofsd);

        res2  = ofs.evaluate(VAR);
        resd2 = ofsd.evaluate(VAR);

        //Comparison. double case
        //--------------
        cout << "1. double case.         " << endl;
        cout << "Expected result: " << creal(res1) << cimag(res1) << "i" << endl;
        cout << "Obtained result: " << creal(res2) << cimag(res2) << "i" << endl;
        res2 -= res1;
        cout << "Delta: " << cabs(res2) << endl;

        //Comparison. cdouble case
        //--------------
        cout << "2. cdouble case.         " << endl;
        cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
        cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
        resd2 -= resd1;
        cout << "Delta: " << cabs(resd2) << endl;
    }

}

/**
 * \fn void oftsh_test_derh()
 * \brief Test of the routine: Oftsh<T>& Oftsh<T>::derh(Oftsh< T > const& a, int ni) and its specializations.
 */
void oftsh_test_derh()
{
    cout << "----------------------------------------------------------" << endl;
    cout << "    Test of the routine:                                  " << endl;
    cout << " Oftsh<T>& Oftsh<T>::derh(Oftsh< T > const& a, int ni).   " << endl;
    cout << "----------------------------------------------------------" << endl;
    double epsilon = 1e-8;
    cout << " The partial derivatives are estimated with a simple Euler scheme with epsilon = " << epsilon << "." << endl;

    //Initialization of OFTS objects, used to initialize the OFTSH objects
    //--------------
    Oftsd xd_ofts(REDUCED_NV, OFTS_ORDER-1, 1, OFS_ORDER);
    Oftsc xdc_ofts(REDUCED_NV, OFTS_ORDER-1, 1, OFS_ORDER);
    Oftsd xd2_ofts(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc2_ofts(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);


    //Initialization
    //--------------
    Oftsh< Ofsd > *xd = xd_ofts.get_term(OFTS_ORDER-1);
    Oftsh< Ofsc > *xdc = xdc_ofts.get_term(OFTS_ORDER-1);
    Oftsh< Ofsd > *xd2 = xd2_ofts.get_term(OFTS_ORDER);
    Oftsh< Ofsc > *xdc2 = xdc2_ofts.get_term(OFTS_ORDER);

    double Xvar[REDUCED_NV];
    double Xvar2[REDUCED_NV];
    cdouble Xvard[REDUCED_NV];
    cdouble Xvard2[REDUCED_NV];

    Ofsd ofs(OFS_ORDER);
    Ofsc ofsd(OFS_ORDER);
    Ofsd ofs2(OFS_ORDER);
    Ofsc ofsd2(OFS_ORDER);

    //result storage
    cdouble res1, res2;
    cdouble resd1, resd2;

    //Loop on all the variables
    for(int ni = 1; ni<= REDUCED_NV ; ni++)
    {
        cout << "---------------------------" << endl;
        cout << "    Variable " << ni << ": " << endl;
        cout << "---------------------------" << endl;
        //array of coordinates (randomly init)
        for(int i = 0; i < REDUCED_NV; i++)
        {
            Xvar[i] = 1.0*rand()/RAND_MAX;
            Xvar2[i] = Xvar[i];
        }
        Xvar2[ni-1] += epsilon;

        for(int i = 0; i < REDUCED_NV; i++)
        {
            Xvard[i] = 1.0*rand()/RAND_MAX + I*rand()/RAND_MAX ;
            Xvard2[i] = Xvard[i];
        }
        Xvard2[ni-1] += epsilon;

        //Set random coefs in xd and xdc
        xd->set_random_coefs();
        xdc->set_random_coefs();
        xd2->set_random_coefs();
        xdc2->set_random_coefs();

        //Set results
        //--------------
        xd2->evaluate(Xvar, ofs);
        xdc2->evaluate(Xvard, ofsd);
        xd2->evaluate(Xvar2, ofs2);
        xdc2->evaluate(Xvard2, ofsd2);


        res1  = (ofs2.evaluate(VAR) - ofs.evaluate(VAR))/epsilon;
        resd1 = (ofsd2.evaluate(VAR) - ofsd.evaluate(VAR))/epsilon;

        //Take derh
        //--------------
        xd->derh(*xd2, ni);
        xdc->derh(*xdc2, ni);

        //Set results
        //--------------
        xd->evaluate(Xvar, ofs);
        xdc->evaluate(Xvard, ofsd);

        res2  = ofs.evaluate(VAR);
        resd2 = ofsd.evaluate(VAR);

        //Comparison. double case
        //--------------
        cout << "1. double case.         " << endl;
        cout << "Expected result: " << creal(res1) << cimag(res1) << "i" << endl;
        cout << "Obtained result: " << creal(res2) << cimag(res2) << "i" << endl;
        res2 -= res1;
        cout << "Delta: " << cabs(res2) << endl;

        //Comparison. cdouble case
        //--------------
        cout << "2. cdouble case.         " << endl;
        cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
        cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
        resd2 -= resd1;
        cout << "Delta: " << cabs(resd2) << endl;
    }

}

/**
 * \fn void oftsh_test_dot()
 * \brief Test of the routine: Oftsh<T>& Oftsh<T>::dot(Oftsh<T> const& a, double const&  n).
 * Note that the Ofs objects are evaluated @ n*var, and not var,
 * to take into account the presence of the pulsation n.
 */
void oftsh_test_dot()
{
    cout << "--------------------------------------------------------------" << endl;
    cout << "    Test of the routine:                                      " << endl;
    cout << " Oftsh<T>& Oftsh<T>::dot(Oftsh<T> const& a, double const&  n)." << endl;
    cout << "--------------------------------------------------------------" << endl;
    double epsilon = 1e-8;
    cout << " The time derivatives are estimated with a simple Euler scheme with epsilon = " << epsilon << "." << endl;

    //Initialization of OFTS objects, used to initialize the OFTSH objects
    //--------------
    Oftsc xdc_ofts(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc2_ofts(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);

    //Initialization
    //--------------
    double n  = 0.925195985520347;

    //Initialization
    //--------------
    Oftsh< Ofsc > *xdc = xdc_ofts.get_term(OFTS_ORDER);
    Oftsh< Ofsc > *xdc2 = xdc2_ofts.get_term(OFTS_ORDER);

    cdouble Xvard[REDUCED_NV];

    Ofsc ofsd(OFS_ORDER);
    Ofsc ofsd2(OFS_ORDER);

    //result storage
    cdouble resd1, resd2;

    //array of coordinates (randomly init)
    for(int i = 0; i < REDUCED_NV; i++) Xvard[i] = 1.0*rand()/RAND_MAX + I*rand()/RAND_MAX ;

    //Set random coefs in xd and xdc
    xdc->set_random_coefs();
    xdc2->set_random_coefs();

    //Set results
    //--------------
    xdc2->evaluate(Xvard, ofsd);
    resd1 = (ofsd.evaluate(n*(VAR+epsilon)) - ofsd.evaluate(n*VAR))/epsilon;

    //Take dot
    //--------------
    xdc->dot(*xdc2, n);

    //Set results
    //--------------
    xdc->evaluate(Xvard, ofsd);
    resd2 = ofsd.evaluate(n*VAR);

    //Comparison. cdouble case
    //--------------
    cout << "Double complex case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;

}

