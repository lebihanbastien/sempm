/**
 * \file ofts_test.cpp
 * \brief Test file for the Ofts (Fourier-Taylor series) class.
 * \author BLB
 */

#include "timec.h"
#include "ofts_test.h"
#define VAR 2.0

using namespace std;

/**
 * \fn void ofts_test()
 * \brief Main routine for Ofts class test.
 */
void ofts_test()
{
    char ch;
    //------------------------------------------------------------------------------------
    // Splash
    //------------------------------------------------------------------------------------
    cout << "---------------------------------------------------" << endl;
    cout << "                                                   " << endl;
    cout << "         Test routine for the Ofts class           " << endl;
    cout << "        specialized in the form Ofts<Ofs>          " << endl;
    cout << "             (Fourier-Taylor series)               " << endl;
    cout << "                                                   " << endl;
    cout << "---------------------------------------------------" << endl;
    cout << "The following public functions are tested:         " << endl;
    cout << "---------------------------------------------------" << endl;
    cout << "Ofts<T>::ofts_smult_t(Ofts<T> const& a, T const& m)" << endl;
    cout << "Ofts<T>::ofts_smult_u(Ofts< Ofs<U> > const& a, U const& c)" << endl;
    cout << "Ofts<T>::ofts_smult_tu(Ofts< Ofs<U> > const& a, T const& m, U const& c)" << endl;
    cout << "Ofts<T>::ofts_sfsum_t(Ofts<T> const& a, T const& ma, Ofts<T> const& b, T const& mb)" << endl;
    cout << "Ofts<T>::ofts_sprod(Ofts<T> const& a, Ofts<T> const& b)" << endl;
    cout << "Ofts<T>::ofts_smprod_t(Ofts<T> const& a, Ofts<T> const& b, T const& m)" << endl;
    cout << "Ofts<T>::ofts_smprod_u(Ofts<T> const& a, Ofts<T> const& b, U const& c)" << endl;
    cout << "Note that the other variations of ofts_sfsum (fsum_t, sfsum_tt, fsum_u) " << endl;
    cout << "are not tested since they present little difference with ofts_sfsum" << endl;
    cout << "See source code for details. " << endl;

    cout <<  setw(5) << setprecision(1) << setiosflags(ios::scientific);
    cout << "---------------------------------------------------" << endl;
    cout << " 1. Tests are made with series of order " << OFTS_ORDER << "." << endl;
    cout << " 2. Tests are made with number of variables: " << REDUCED_NV << "." << endl;
    cout << " 3. Each coefficient is a series of order " << OFS_ORDER << "." << endl;
    cout << " 4. The series are evaluated at arbitrary coordinates." << endl;
    cout << " 5. Random coefficients are set as input into every Ofs coefficient. " << endl;
    cout << "However, an arbitrary decreasing of the coefficients is set: " << endl;
    cout << " around " << (double)rand()/(10.0*(pow(0,7.0)+1)*RAND_MAX) << " @ order 0." << endl;
    cout << " around " << (double)rand()/(10.0*(pow(OFS_ORDER,7.0)+1)*RAND_MAX) << " @ order "<< OFS_ORDER << "." << endl;
    cout << " 7. WARNING: in order for the Ofts*Ofts products to make sense, " << endl;
    cout << " the original series are filled with random coefficients up to order " << OFTS_ORDER/2 << "(half the max order)" << endl;
    cout << " so that the results of the product contains all relevant terms. " << endl;
    cout << " Note that a more realistic process would be to ensure a decrease of the coefficients " << endl;
    cout << " with respect to the order of the homogeneous FT polynomials. " << endl;
    cout << "---------------------------------------------------" << endl;
    cout << " WARNING: orders > 10 coupled with a number of variable > 4 " << endl;
    cout << " can lead to a time-consuming test session. " << endl;
    cout <<  setw(5) << setprecision(15) << std::showpos  <<  setiosflags(ios::scientific);
    cout << "Press Enter to proceed with the tests." << endl;
    scanf("%c",&ch);


    //------------------------------------------------------------------------------------
    // Tests.
    //------------------------------------------------------------------------------------
    ofts_test_sprod();
    cout << "Press Enter to proceed with the tests." << endl;
    scanf("%c",&ch);
    cout << endl;

    ofts_test_smult_t();
    cout << "Press Enter to proceed with the tests." << endl;
    scanf("%c",&ch);
    cout << endl;


    ofts_test_pows();
    cout << "Press Enter to proceed with the tests." << endl;
    scanf("%c",&ch);
    cout << endl;


    ofts_test_smprod_t();
    cout << "Press Enter to proceed with the tests." << endl;
    scanf("%c",&ch);
    cout << endl;


    ofts_test_conjugate();
    cout << "Press Enter to proceed with the tests." << endl;
    scanf("%c",&ch);
    cout << endl;



    ofts_test_smult_u();
    cout << "Press Enter to proceed with the tests." << endl;
    scanf("%c",&ch);
    cout << endl;

    ofts_test_smult_tu();
    cout << "Press Enter to proceed with the tests." << endl;
    scanf("%c",&ch);
    cout << endl;

    ofts_test_sfsum_t();
    cout << "Press Enter to proceed with the tests." << endl;
    scanf("%c",&ch);
    cout << endl;

    ofts_test_smprod_u();
    cout << "Press Enter to proceed with the tests." << endl;
    scanf("%c",&ch);
    cout << endl;

    cout << "---------------------------------------------------" << endl;
    cout << "End of test session.                               " << endl;
    cout << "---------------------------------------------------" << endl;
}

/**
 * \fn void ofts_test_conjugate()
 * \brief Test of the routine: Ofts<T>::conjugate().
 */
void ofts_test_conjugate()
{
    cout << "--------------------------------------" << endl;
    cout << " Test of the routines:                " << endl;
    cout << " Ofts<T>::conjugate()                 " << endl;
    cout << "--------------------------------------" << endl;
    //Initialization of OFTS objects, used to initialize the OFTSH objects
    //--------------
    Oftsd xd(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsd xd2(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc2(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);

    double Xvar[REDUCED_NV];
    cdouble Xvard[REDUCED_NV];

    Ofsd ofs(OFS_ORDER);
    Ofsc ofsd(OFS_ORDER);
    Ofsd ofs2(OFS_ORDER);
    Ofsc ofsd2(OFS_ORDER);

    //result storage
    cdouble res1, res2;
    cdouble resd1, resd2;
    //array of coordinates (randomly init)
    for(int i = 0; i < REDUCED_NV; i++) Xvar[i] = 1.0*rand()/RAND_MAX;
    for(int i = 0; i < REDUCED_NV; i++) Xvard[i] = 1.0*rand()/RAND_MAX + I*rand()/RAND_MAX ;

    //Set random coefs in xd and xdc
    xd.set_random_coefs();
    xdc.set_random_coefs();
    xd2.set_random_coefs();
    xdc2.set_random_coefs();


    cout << "-----------------------" << endl;
    cout << " 1. At all terms       " << endl;
    cout << "-----------------------" << endl;

    //Set results
    //--------------
    xd.evaluate(Xvar, ofs);
    xdc.evaluate_conjugate(Xvard, ofsd);
    res1  = conj(ofs.evaluate(VAR));
    resd1 = conj(ofsd.evaluate(VAR));

    //Take operation
    //--------------
    xd.conjugate();
    xdc.conjugate();

    //Set results
    //--------------
    xd.evaluate(Xvar, ofs);
    xdc.evaluate(Xvard, ofsd);

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


    cout << "-----------------------" << endl;
    cout << " 2. At order n         " << endl;
    cout << "-----------------------" << endl;

    //Set results
    //--------------
    xd.evaluate(Xvar, ofs);
    xdc.evaluate_conjugate(Xvard, ofsd);
    res1  = conj(ofs.evaluate(VAR));
    resd1 = conj(ofsd.evaluate(VAR));

    //Take operation
    //--------------
    for(int k = 0; k<= xd.get_order(); k++) xd.conjugate(k);
    for(int k = 0; k<= xd.get_order(); k++) xdc.conjugate(k);

    //Set results
    //--------------
    xd.evaluate(Xvar, ofs);
    xdc.evaluate(Xvard, ofsd);

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
 * \fn void ofts_test_smult_t()
 * \brief Test of the routine: Ofts<T>::ofts_smult_t(Ofts<T> const& a, T const& m).
 */
void ofts_test_smult_t()
{
    cout << "-------------------------------------------------------------------------" << endl;
    cout << " Test of the routines:                                                   " << endl;
    cout << " Ofts<T>::ofts_smult_t(Ofts<T> const& a, T const& m)                     " << endl;
    cout << "-------------------------------------------------------------------------" << endl;
    //Initialization of OFTS objects, used to initialize the OFTSH objects
    //--------------
    Oftsd xd(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsd xd2(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc2(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsd xd3(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc3(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    //Coefficient array
    double Xvar[REDUCED_NV];
    cdouble Xvard[REDUCED_NV];
    //Ofs objects to receive the evaluations
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
    for(int i = 0; i < REDUCED_NV; i++) Xvar[i]  = 1.0*rand()/RAND_MAX;
    for(int i = 0; i < REDUCED_NV; i++) Xvard[i] = 1.0*rand()/RAND_MAX + I*rand()/RAND_MAX ;
    //Set random coefs
    xd.set_random_coefs();
    xdc.set_random_coefs();
    xd2.set_random_coefs();
    xdc2.set_random_coefs();
    xd3.set_random_coefs();
    xdc3.set_random_coefs();

    cout << "-----------------------" << endl;
    cout << " 1. At all terms       " << endl;
    cout << "-----------------------" << endl;

    //Set results
    //--------------
    xd.evaluate(Xvar, ofs);
    xdc.evaluate(Xvard, ofsd);

    xd2.evaluate(Xvar, ofs2);
    xdc2.evaluate(Xvard, ofsd2);

    xd3.evaluate(Xvar, ofs3);
    xdc3.evaluate(Xvard, ofsd3);

    res1  = ofs.evaluate(VAR)  + ofs3.evaluate(VAR)*ofs2.evaluate(VAR);
    resd1 = ofsd.evaluate(VAR) + ofsd3.evaluate(VAR)*ofsd2.evaluate(VAR);

    //Take operation
    //--------------
    xd.ofts_smult_t(xd2, ofs3);
    xdc.ofts_smult_t(xdc2, ofsd3);


    //Set results
    //--------------
    xd.evaluate(Xvar, ofs);
    xdc.evaluate(Xvard, ofsd);

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


    cout << "-----------------------" << endl;
    cout << " 2. At order n         " << endl;
    cout << "-----------------------" << endl;

    //Set results
    //--------------
    xd.evaluate(Xvar, ofs);
    xdc.evaluate(Xvard, ofsd);

    xd2.evaluate(Xvar, ofs2);
    xdc2.evaluate(Xvard, ofsd2);

    xd3.evaluate(Xvar, ofs3);
    xdc3.evaluate(Xvard, ofsd3);

    res1  = ofs.evaluate(VAR)  + ofs3.evaluate(VAR)*ofs2.evaluate(VAR);
    resd1 = ofsd.evaluate(VAR) + ofsd3.evaluate(VAR)*ofsd2.evaluate(VAR);

    //xdc in TFS format
    //--------------
    Oftsc xdc4(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    xdc4.tfs_from_ofs(xdc);         //xdc in TFS format

    //Take operation, in OFS format
    //--------------
    for(int k = 0; k<= xd.get_order(); k++) xd.ofts_smult_t(xd2, ofs3, k);
    for(int k = 0; k<= xd.get_order(); k++) xdc.ofts_smult_t(xdc2, ofsd3, k);

    //Take operation, in TFS format
    //--------------
    xdc3.tfs_from_ofs(xdc2);        //xdc2 in TFS format
    ofsd2.tfs_from_ofs(ofsd3);      //ofsd3 in TFS format
    for(int k = 0; k<= xd.get_order(); k++) xdc4.tfts_smult_t(xdc3, ofsd2, k);
    xdc.tfs_to_ofs(xdc4);           //Back to OFS format

    //Set results
    //--------------
    xd.evaluate(Xvar, ofs);
    xdc.evaluate(Xvard, ofsd);

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
 * \fn void ofts_test_smult_u()
 * \brief Test of the routine: Ofts<T>::ofts_smult_u(Ofts< Ofs<U> > const& a, U const& c).
 */
void ofts_test_smult_u()
{
    cout << "-------------------------------------------------------------------------" << endl;
    cout << " Test of the routines:                                                   " << endl;
    cout << " Ofts<T>::ofts_smult_u(Ofts< Ofs<U> > const& a, U const& c)                     " << endl;
    cout << "-------------------------------------------------------------------------" << endl;
    //Initialization of OFTS objects, used to initialize the OFTSH objects
    //--------------
    Oftsd xd(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsd xd2(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc2(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    //Coefficient array
    double Xvar[REDUCED_NV];
    cdouble Xvard[REDUCED_NV];
    //Ofs objects to receive the evaluations
    Ofsd ofs(OFS_ORDER);
    Ofsc ofsd(OFS_ORDER);
    Ofsd ofs2(OFS_ORDER);
    Ofsc ofsd2(OFS_ORDER);
    //result storage
    cdouble res1, res2;
    cdouble resd1, resd2;
    //coeff
    double c = 2.0;
    cdouble cd = 2.0+I;
    //array of coordinates (randomly init)
    for(int i = 0; i < REDUCED_NV; i++) Xvar[i]  = 1.0*rand()/RAND_MAX;
    for(int i = 0; i < REDUCED_NV; i++) Xvard[i] = 1.0*rand()/RAND_MAX + I*rand()/RAND_MAX ;
    //Set random coefs
    xd.set_random_coefs();
    xdc.set_random_coefs();
    xd2.set_random_coefs();
    xdc2.set_random_coefs();

    cout << "-----------------------" << endl;
    cout << " 1. At all terms       " << endl;
    cout << "-----------------------" << endl;

    //Set results
    //--------------
    xd.evaluate(Xvar, ofs);
    xdc.evaluate(Xvard, ofsd);

    xd2.evaluate(Xvar, ofs2);
    xdc2.evaluate(Xvard, ofsd2);

    res1  = ofs.evaluate(VAR)  + c*ofs2.evaluate(VAR);
    resd1 = ofsd.evaluate(VAR) + cd*ofsd2.evaluate(VAR);

    //Take operation
    //--------------
    xd.ofts_smult_u(xd2, c);
    xdc.ofts_smult_u(xdc2, cd);

    //Set results
    //--------------
    xd.evaluate(Xvar, ofs);
    xdc.evaluate(Xvard, ofsd);

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


    cout << "-----------------------" << endl;
    cout << " 2. At order n         " << endl;
    cout << "-----------------------" << endl;

    //Set results
    //--------------
    xd.evaluate(Xvar, ofs);
    xdc.evaluate(Xvard, ofsd);

    xd2.evaluate(Xvar, ofs2);
    xdc2.evaluate(Xvard, ofsd2);

    res1  = ofs.evaluate(VAR)  + c*ofs2.evaluate(VAR);
    resd1 = ofsd.evaluate(VAR) + cd*ofsd2.evaluate(VAR);

    //Take operation
    //--------------
    for(int k = 0; k<= xd.get_order(); k++) xd.ofts_smult_u(xd2, c, k);
    for(int k = 0; k<= xd.get_order(); k++) xdc.ofts_smult_u(xdc2, cd, k);

    //Set results
    //--------------
    xd.evaluate(Xvar, ofs);
    xdc.evaluate(Xvard, ofsd);

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
 * \fn void ofts_test_smult_tu()
 * \brief Test of the routine: Ofts<T>::ofts_mult_tu(Ofts< Ofs<U> > const& a, T const& m, U const& c).
 */
void ofts_test_smult_tu()
{
    cout << "-------------------------------------------------------------------------" << endl;
    cout << " Test of the routines:                                                   " << endl;
    cout << " Ofts<T>::ofts_mult_tu(Ofts< Ofs<U> > const& a, T const& m, U const& c)  " << endl;
    cout << "-------------------------------------------------------------------------" << endl;
    //Initialization of OFTS objects, used to initialize the OFTSH objects
    //--------------
    Oftsd xd(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsd xd2(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc2(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsd xd3(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc3(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    //Coefficient array
    double Xvar[REDUCED_NV];
    cdouble Xvard[REDUCED_NV];
    //Ofs objects to receive the evaluations
    Ofsd ofs(OFS_ORDER);
    Ofsc ofsd(OFS_ORDER);
    Ofsd ofs2(OFS_ORDER);
    Ofsc ofsd2(OFS_ORDER);
    Ofsd ofs3(OFS_ORDER);
    Ofsc ofsd3(OFS_ORDER);
    //result storage
    cdouble res1, res2;
    cdouble resd1, resd2;
    //coeff
    double c = 2.0;
    cdouble cd = 2.0+I;
    //array of coordinates (randomly init)
    for(int i = 0; i < REDUCED_NV; i++) Xvar[i]  = 1.0*rand()/RAND_MAX;
    for(int i = 0; i < REDUCED_NV; i++) Xvard[i] = 1.0*rand()/RAND_MAX + I*rand()/RAND_MAX ;
    //Set random coefs
    xd.set_random_coefs();
    xdc.set_random_coefs();
    xd2.set_random_coefs();
    xdc2.set_random_coefs();
    xd3.set_random_coefs();
    xdc3.set_random_coefs();

    cout << "-----------------------" << endl;
    cout << " 1. At all terms       " << endl;
    cout << "-----------------------" << endl;

    //Set results
    //--------------
    xd.evaluate(Xvar, ofs);
    xdc.evaluate(Xvard, ofsd);

    xd2.evaluate(Xvar, ofs2);
    xdc2.evaluate(Xvard, ofsd2);

    xd3.evaluate(Xvar, ofs3);
    xdc3.evaluate(Xvard, ofsd3);

    res1  = ofs.evaluate(VAR)  + c*ofs3.evaluate(VAR)*ofs2.evaluate(VAR);
    resd1 = ofsd.evaluate(VAR) + cd*ofsd3.evaluate(VAR)*ofsd2.evaluate(VAR);

    //Take operation
    //--------------
    xd.ofts_smult_tu(xd2, ofs3, c);
    xdc.ofts_smult_tu(xdc2, ofsd3, cd);

    //Set results
    //--------------
    xd.evaluate(Xvar, ofs);
    xdc.evaluate(Xvard, ofsd);

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


    cout << "-----------------------" << endl;
    cout << " 2. At order n         " << endl;
    cout << "-----------------------" << endl;

    //Set results
    //--------------
    xd.evaluate(Xvar, ofs);
    xdc.evaluate(Xvard, ofsd);

    xd2.evaluate(Xvar, ofs2);
    xdc2.evaluate(Xvard, ofsd2);

    xd3.evaluate(Xvar, ofs3);
    xdc3.evaluate(Xvard, ofsd3);

    res1  = ofs.evaluate(VAR)  + c*ofs3.evaluate(VAR)*ofs2.evaluate(VAR);
    resd1 = ofsd.evaluate(VAR) + cd*ofsd3.evaluate(VAR)*ofsd2.evaluate(VAR);

    //Take operation
    //--------------
    for(int k = 0; k<= xd.get_order(); k++)  xd.ofts_smult_tu(xd2, ofs3, c, k);
    for(int k = 0; k<= xd.get_order(); k++)  xdc.ofts_smult_tu(xdc2, ofsd3, cd, k);

    //Set results
    //--------------
    xd.evaluate(Xvar, ofs);
    xdc.evaluate(Xvard, ofsd);

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
 * \fn void ofts_test_sfsum_t()
 * \brief Test of the routine: Ofts<T>::ofts_sfsum_t(Ofts<T> const& a, T const& m).
 */
void ofts_test_sfsum_t()
{
    cout << "-------------------------------------------------------------------------" << endl;
    cout << " Test of the routines:                                                   " << endl;
    cout << " Ofts<T>::ofts_sfsum_t(Ofts<T> const& a, T const& m)                     " << endl;
    cout << "-------------------------------------------------------------------------" << endl;
    //Initialization of OFTS objects, used to initialize the OFTSH objects
    //--------------
    Oftsd xd(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsd xd2(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc2(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsd xd3(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc3(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsd xd4(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc4(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsd xd5(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc5(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);

    //Coefficient array
    double Xvar[REDUCED_NV];
    cdouble Xvard[REDUCED_NV];
    //Ofs objects to receive the evaluations
    Ofsd ofs(OFS_ORDER);
    Ofsc ofsd(OFS_ORDER);
    Ofsd ofs2(OFS_ORDER);
    Ofsc ofsd2(OFS_ORDER);
    Ofsd ofs3(OFS_ORDER);
    Ofsc ofsd3(OFS_ORDER);
    Ofsd ofs4(OFS_ORDER);
    Ofsc ofsd4(OFS_ORDER);
    Ofsd ofs5(OFS_ORDER);
    Ofsc ofsd5(OFS_ORDER);
    //result storage
    cdouble res1, res2;
    cdouble resd1, resd2;
    //array of coordinates (randomly init)
    for(int i = 0; i < REDUCED_NV; i++) Xvar[i]  = 1.0*rand()/RAND_MAX;
    for(int i = 0; i < REDUCED_NV; i++) Xvard[i] = 1.0*rand()/RAND_MAX + I*rand()/RAND_MAX ;
    //Set random coefs
    xd.set_random_coefs();
    xdc.set_random_coefs();
    xd2.set_random_coefs();
    xdc2.set_random_coefs();
    xd3.set_random_coefs();
    xdc3.set_random_coefs();
    xd4.set_random_coefs();
    xdc4.set_random_coefs();
    xd5.set_random_coefs();
    xdc5.set_random_coefs();

    cout << "-----------------------" << endl;
    cout << " 1. At all terms       " << endl;
    cout << "-----------------------" << endl;

    //Set results
    //--------------
    xd.evaluate(Xvar, ofs);
    xdc.evaluate(Xvard, ofsd);

    xd2.evaluate(Xvar, ofs2);
    xdc2.evaluate(Xvard, ofsd2);

    xd3.evaluate(Xvar, ofs3);
    xdc3.evaluate(Xvard, ofsd3);

    //Just to set random coeffs
    xd4.evaluate(Xvar, ofs4);
    xdc4.evaluate(Xvard, ofsd4);
    xd5.evaluate(Xvar, ofs5);
    xdc5.evaluate(Xvard, ofsd5);


    res1  = ofs.evaluate(VAR)  + ofs4.evaluate(VAR)*ofs2.evaluate(VAR)  + ofs5.evaluate(VAR)*ofs3.evaluate(VAR);
    resd1 = ofsd.evaluate(VAR)  + ofsd4.evaluate(VAR)*ofsd2.evaluate(VAR) + ofsd5.evaluate(VAR)*ofsd3.evaluate(VAR);

    //Take operation
    //--------------
    xd.ofts_sfsum_t(xd2, ofs4, xd3, ofs5);
    xdc.ofts_sfsum_t(xdc2, ofsd4, xdc3, ofsd5);


    //Set results
    //--------------
    xd.evaluate(Xvar, ofs);
    xdc.evaluate(Xvard, ofsd);

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


    cout << "-----------------------" << endl;
    cout << " 2. At order n         " << endl;
    cout << "-----------------------" << endl;

    //Set results
    //--------------
    xd.evaluate(Xvar, ofs);
    xdc.evaluate(Xvard, ofsd);

    xd2.evaluate(Xvar, ofs2);
    xdc2.evaluate(Xvard, ofsd2);

    xd3.evaluate(Xvar, ofs3);
    xdc3.evaluate(Xvard, ofsd3);

    //Just to set random coeffs
    xd4.evaluate(Xvar, ofs4);
    xdc4.evaluate(Xvard, ofsd4);
    xd5.evaluate(Xvar, ofs5);
    xdc5.evaluate(Xvard, ofsd5);

    res1  = ofs.evaluate(VAR)  + ofs4.evaluate(VAR)*ofs2.evaluate(VAR)  + ofs5.evaluate(VAR)*ofs3.evaluate(VAR);
    resd1 = ofsd.evaluate(VAR)  + ofsd4.evaluate(VAR)*ofsd2.evaluate(VAR) + ofsd5.evaluate(VAR)*ofsd3.evaluate(VAR);

    //Take operation
    //--------------
    for(int k = 0; k<= xd.get_order(); k++) xd.ofts_sfsum_t(xd2, ofs4, xd3, ofs5, k);
    for(int k = 0; k<= xd.get_order(); k++) xdc.ofts_sfsum_t(xdc2, ofsd4, xdc3, ofsd5, k);


    //Set results
    //--------------
    xd.evaluate(Xvar, ofs);
    xdc.evaluate(Xvard, ofsd);

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
 * \fn void ofts_test_sprod()
 * \brief Test of the routine: Ofts<T>::ofts_sprod(Ofts<T> const& a, Ofts<T> const& b).
 */
void ofts_test_sprod()
{
    cout << "-------------------------------------------------------------------------" << endl;
    cout << " Test of the routines:                                                   " << endl;
    cout << " Ofts<T>::ofts_sprod(Ofts<T> const& a, Ofts<T> const& b)                 " << endl;
    cout << "-------------------------------------------------------------------------" << endl;
    //Initialization of OFTS objects, used to initialize the OFTSH objects
    //--------------
    Oftsd xd(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsd xd2(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc2(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsd xd3(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc3(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    //Coefficient array
    double Xvar[REDUCED_NV];
    cdouble Xvard[REDUCED_NV];
    //Ofs objects to receive the evaluations
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
    for(int i = 0; i < REDUCED_NV; i++) Xvar[i]  = 1.0*rand()/RAND_MAX;
    for(int i = 0; i < REDUCED_NV; i++) Xvard[i] = 1.0*rand()/RAND_MAX + I*rand()/RAND_MAX ;


    cout << "-----------------------" << endl;
    cout << " 1. At all terms       " << endl;
    cout << "-----------------------" << endl;
    //Set random coefs
    xd.set_random_coefs();
    xdc.set_random_coefs();
    xd2.set_random_coefs();
    xdc2.set_random_coefs();
    xd3.set_random_coefs();
    xdc3.set_random_coefs();

    //Set results
    //--------------
    xd.evaluate(Xvar, ofs);
    xdc.evaluate(Xvard, ofsd);

    xd2.evaluate(Xvar, ofs2);
    xdc2.evaluate(Xvard, ofsd2);

    xd3.evaluate(Xvar, ofs3);
    xdc3.evaluate(Xvard, ofsd3);

    res1  = ofs.evaluate(VAR)  + ofs3.evaluate(VAR)*ofs2.evaluate(VAR);
    resd1 = ofsd.evaluate(VAR) + ofsd3.evaluate(VAR)*ofsd2.evaluate(VAR);

    //Take product
    //--------------
    tic();
    xd.ofts_sprod(xd2, xd3);
    cout << "Product in : " << toc() << "s." << endl;
    tic();
    xdc.ofts_sprod(xdc2, xdc3);
    cout << "Product in : " << toc() << "s." << endl;


    //Set results
    //--------------
    xd.evaluate(Xvar, ofs);
    xdc.evaluate(Xvard, ofsd);

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


    cout << "-----------------------" << endl;
    cout << " 2. At order n         " << endl;
    cout << "-----------------------" << endl;
    xd.zero();
    xdc.zero();

    //Set random coefs
    xd.set_random_coefs();
    xdc.set_random_coefs();
    xd2.set_random_coefs();
    xdc2.set_random_coefs();
    xd3.set_random_coefs();
    xdc3.set_random_coefs();

    //Set results
    //--------------
    xd.evaluate(Xvar, ofs);
    xdc.evaluate(Xvard, ofsd);

    xd2.evaluate(Xvar, ofs2);
    xdc2.evaluate(Xvard, ofsd2);

    xd3.evaluate(Xvar, ofs3);
    xdc3.evaluate(Xvard, ofsd3);

    res1  = ofs.evaluate(VAR)  + ofs3.evaluate(VAR)*ofs2.evaluate(VAR);
    resd1 = ofsd.evaluate(VAR) + ofsd3.evaluate(VAR)*ofsd2.evaluate(VAR);

    //xdc in TFS format
    //--------------
    xdc.tfs_from_ofs_inline();   //xdc in TFS format
    xdc2.tfs_from_ofs_inline();  //xdc in TFS format
    xdc3.tfs_from_ofs_inline();  //xdc in TFS format

    //Take product, in OFS format
    //--------------
    tic();
    for(int k = 0; k<= xd.get_order(); k++) xd.ofts_sprod(xd2, xd3, k);
    cout << "Product in : " << toc() << "s." << endl;

    //Take operation, in TFS format
    //--------------
    for(int k = 0; k<= xd.get_order(); k++) xdc.tfts_sprod(xdc2, xdc3, k);
    xdc.tfs_to_ofs_inline(); //back in OFS format

    //Set results
    //--------------
    xd.evaluate(Xvar, ofs);
    xdc.evaluate(Xvard, ofsd);

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
 * \fn void ofts_test_smprod_t()
 * \brief Test of the routine: Ofts<T>::ofts_smprod_t(Ofts<T> const& a, Ofts<T> const& b, T const& m).
 */
void ofts_test_smprod_t()
{
    cout << "-------------------------------------------------------------------------" << endl;
    cout << " Test of the routines:                                                   " << endl;
    cout << " Ofts<T>::ofts_smprod_t(Ofts<T> const& a, Ofts<T> const& b, T const& m)  " << endl;
    cout << "-------------------------------------------------------------------------" << endl;
    //Initialization of OFTS objects, used to initialize the OFTSH objects
    //--------------
    Oftsd xd(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsd xd2(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc2(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsd xd3(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc3(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);

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


    cout << "-----------------------" << endl;
    cout << " 1. At all terms       " << endl;
    cout << "-----------------------" << endl;
    //Set random coefs in xd and xdc
    xd.set_random_coefs();
    xdc.set_random_coefs();
    xd2.set_random_coefs();
    xdc2.set_random_coefs();
    xd3.set_random_coefs();
    xdc3.set_random_coefs();

    //Set results
    //--------------
    xd.evaluate(Xvar, ofs);
    xdc.evaluate(Xvard, ofsd);
    xd2.evaluate(Xvar, ofs2);
    xdc2.evaluate(Xvard, ofsd2);
    xd3.evaluate(Xvar, ofs3);
    xdc3.evaluate(Xvard, ofsd3);
    //Just to set random coefficients
    xd.evaluate(Xvar, ofs4);
    xdc.evaluate(Xvard, ofsd4);

    res1  = ofs.evaluate(VAR)  + ofs4.evaluate(VAR)*ofs3.evaluate(VAR)*ofs2.evaluate(VAR);
    resd1 = ofsd.evaluate(VAR) + ofsd4.evaluate(VAR)*ofsd3.evaluate(VAR)*ofsd2.evaluate(VAR);

    //Take operation
    //--------------
    tic();
    xd.ofts_smprod_t(xd2, xd3, ofs4, temp);
    cout << "Product in : " << toc() << "s." << endl;
    tic();
    xdc.ofts_smprod_t(xdc2, xdc3, ofsd4, tempc);
    cout << "Product in : " << toc() << "s." << endl;


    //Set results
    //--------------
    xd.evaluate(Xvar, ofs);
    xdc.evaluate(Xvard, ofsd);

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


    cout << "-----------------------" << endl;
    cout << " 2. At order n         " << endl;
    cout << "-----------------------" << endl;
    xd.zero();
    xdc.zero();

    //Set random coefs in xd and xdc
    xd.set_random_coefs();
    xdc.set_random_coefs();
    xd2.set_random_coefs();
    xdc2.set_random_coefs();
    xd3.set_random_coefs();
    xdc3.set_random_coefs();

    //Set results
    //--------------
    xd.evaluate(Xvar, ofs);
    xdc.evaluate(Xvard, ofsd);
    xd2.evaluate(Xvar, ofs2);
    xdc2.evaluate(Xvard, ofsd2);
    xd3.evaluate(Xvar, ofs3);
    xdc3.evaluate(Xvard, ofsd3);
    //Just to set random coefficients
    xd.evaluate(Xvar, ofs4);
    xdc.evaluate(Xvard, ofsd4);

    res1  = ofs.evaluate(VAR)  + ofs4.evaluate(VAR)*ofs3.evaluate(VAR)*ofs2.evaluate(VAR);
    resd1 = ofsd.evaluate(VAR) + ofsd4.evaluate(VAR)*ofsd3.evaluate(VAR)*ofsd2.evaluate(VAR);

    //Take operation
    //--------------
    tic();
    for(int k = 0; k<= xd.get_order(); k++) xd.ofts_smprod_t(xd2, xd3, ofs4, k, temp);
    cout << "Product in : " << toc() << "s." << endl;
    tic();
    for(int k = 0; k<= xd.get_order(); k++) xdc.ofts_smprod_t(xdc2, xdc3, ofsd4, k, tempc);
    cout << "Product in : " << toc() << "s." << endl;


    //Set results
    //--------------
    xd.evaluate(Xvar, ofs);
    xdc.evaluate(Xvard, ofsd);

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
 * \fn void ofts_test_smprod_u()
 * \brief Test of the routine: Ofts<T>::ofts_smprod_u(Ofts<T> const& a, Ofts<T> const& b, U const& c).
 */
void ofts_test_smprod_u()
{
    cout << "-------------------------------------------------------------------------" << endl;
    cout << " Test of the routines:                                                   " << endl;
    cout << " Ofts<T>::ofts_smprod_u(Ofts<T> const& a, Ofts<T> const& b, U const& c)  " << endl;
    cout << "-------------------------------------------------------------------------" << endl;
    //Initialization of OFTS objects, used to initialize the OFTSH objects
    //--------------
    Oftsd xd(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsd xd2(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc2(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsd xd3(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc3(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);

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

    //coeff
    double c = 2.0;
    cdouble cd = 2.0+I;

    //array of coordinates (randomly init)
    for(int i = 0; i < REDUCED_NV; i++) Xvar[i] = 1.0*rand()/RAND_MAX;
    for(int i = 0; i < REDUCED_NV; i++) Xvard[i] = 1.0*rand()/RAND_MAX + I*rand()/RAND_MAX ;

    cout << "-----------------------" << endl;
    cout << " 1. At all terms       " << endl;
    cout << "-----------------------" << endl;

    //Set random coefs in xd and xdc
    xd.set_random_coefs();
    xdc.set_random_coefs();
    xd2.set_random_coefs();
    xdc2.set_random_coefs();
    xd3.set_random_coefs();
    xdc3.set_random_coefs();

    //Set results
    //--------------
    xd.evaluate(Xvar, ofs);
    xdc.evaluate(Xvard, ofsd);
    xd2.evaluate(Xvar, ofs2);
    xdc2.evaluate(Xvard, ofsd2);
    xd3.evaluate(Xvar, ofs3);
    xdc3.evaluate(Xvard, ofsd3);

    res1  = ofs.evaluate(VAR)  + c*ofs3.evaluate(VAR)*ofs2.evaluate(VAR);
    resd1 = ofsd.evaluate(VAR) + cd*ofsd3.evaluate(VAR)*ofsd2.evaluate(VAR);

    //Take operation
    //--------------
    tic();
    xd.ofts_smprod_u(xd2, xd3, c);
    cout << "Product in : " << toc() << "s." << endl;
    tic();
    xdc.ofts_smprod_u(xdc2, xdc3, cd);
    cout << "Product in : " << toc() << "s." << endl;


    //Set results
    //--------------
    xd.evaluate(Xvar, ofs);
    xdc.evaluate(Xvard, ofsd);

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

    cout << "-----------------------" << endl;
    cout << " 2. At order n         " << endl;
    cout << "-----------------------" << endl;
    xd.zero();
    xdc.zero();

    //Set random coefs in xd and xdc
    xd.set_random_coefs();
    xdc.set_random_coefs();
    xd2.set_random_coefs();
    xdc2.set_random_coefs();
    xd3.set_random_coefs();
    xdc3.set_random_coefs();

    //Set results
    //--------------
    xd.evaluate(Xvar, ofs);
    xdc.evaluate(Xvard, ofsd);
    xd2.evaluate(Xvar, ofs2);
    xdc2.evaluate(Xvard, ofsd2);
    xd3.evaluate(Xvar, ofs3);
    xdc3.evaluate(Xvard, ofsd3);

    res1  = ofs.evaluate(VAR)  + c*ofs3.evaluate(VAR)*ofs2.evaluate(VAR);
    resd1 = ofsd.evaluate(VAR) + cd*ofsd3.evaluate(VAR)*ofsd2.evaluate(VAR);

    //Take operation
    //--------------
    tic();
    for(int k = 0; k<= xd.get_order(); k++) xd.ofts_smprod_u(xd2, xd3, c, k);
    cout << "Product in : " << toc() << "s." << endl;
    tic();
    for(int k = 0; k<= xd.get_order(); k++) xdc.ofts_smprod_u(xdc2, xdc3, cd, k);
    cout << "Product in : " << toc() << "s." << endl;


    //Set results
    //--------------
    xd.evaluate(Xvar, ofs);
    xdc.evaluate(Xvard, ofsd);

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
 * \fn void ofts_test_pows()
 * \brief Test of the routine: Ofts<T>& Ofts<T>::pows(Ofts<T> const& a,  U const& alpha) .
 */
void ofts_test_pows()
{
    cout << "---------------------------------------------------------------------------" << endl;
    cout << " Test of the routine with cdouble coefficients:                     " << endl;
    cout << " Ofts<T>& Ofts<T>::pows(Ofts<T> const& a,  U const& alpha)                 " << endl;
    cout << "---------------------------------------------------------------------------" << endl;
    cout << "WARNING: Works if and only if the (Fourier) order 0                        " << endl;
    cout << "of the (Taylor) order 0 is >> to everybody else. Otherwise: pb with inverse" << endl;
    cout << "---------------------------------------------------------------------------" << endl;
    //------------------------------------------------------------------------------
    //Initialization of OFTS objects, used to initialize the OFTSH objects
    //------------------------------------------------------------------------------
    Oftsc xdc(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc2(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc3(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);
    Oftsc xdc4(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);

    cdouble Xvard[REDUCED_NV];
    Ofsc ofsd(OFS_ORDER);
    Ofsc ofsd2(OFS_ORDER);

    //result storage
    cdouble resd1, resd2;
    //Exponent
    double alpha = -1.5;

    //array of coordinates (randomly init)
    for(int i = 0; i < REDUCED_NV; i++) Xvard[i] = 1.0*rand()/RAND_MAX + I*rand()/RAND_MAX ;

    cout << "-----------------------" << endl;
    cout << " 1. At all terms       " << endl;
    cout << "-----------------------" << endl;

    //Read coefs from file
    read_ofts_txt(xdc2, "test/W[0].txt", OFS_ORDER);
    //xdc2.set_random_coefs();

    //------------------------------------------------------------------------------
    //Take operation
    //------------------------------------------------------------------------------
    //--------------
    //1. Mixed power
    //--------------
    tic();
    //xdc.ofts_pows(xdc2, alpha);
    double tmix = toc();
    cout << "Mixed Power in : " << tmix << "s." << endl;

    //--------------
    //2. Pure TFS power
    //--------------
    tic();
    xdc3.tfs_from_ofs(xdc2);
    double t1 = toc();

    tic();
    xdc4.tfts_pows(xdc3, alpha);
    double tpure = toc();
    cout << "Pure TFS Power in : " << tpure << "s." << endl;

    tic();
    xdc.tfs_to_ofs(xdc4);
    double t2 = toc();

    cout << "OFS to TFS back to OFS in : " << t1+t2 << "s." << endl;
    cout << "Thus the ratio is: " << tmix/(tpure+t1+t2) << endl;
    cout << "Without the OFS <-> TFS: " << tmix/(tpure) << endl;
    cout << "To be compared with " << 1.0*OFTS_ORDER/log(OFTS_ORDER) << endl;

    //Set results
    //--------------
    xdc2.evaluate(Xvard, ofsd);
    resd1 = cpow(ofsd.evaluate(VAR), alpha);
    xdc.evaluate(Xvard, ofsd);
    resd2 = ofsd.evaluate(VAR);


    //Comparison. cdouble case
    //--------------
    cout << "2. cdouble case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;


    cout << "-----------------------" << endl;
    cout << " 1. At order n         " << endl;
    cout << "-----------------------" << endl;
    //------------------------------------------------------------------------------
    //Take operation
    //------------------------------------------------------------------------------

    //--------------
    //1. Mixed power
    //--------------
    tic();
    for(int k = 0; k<= xdc.get_order(); k++) xdc.ofts_pows(xdc2, alpha, k);
    cout << "Mixed Pows in : " << toc() << "s." << endl;

    //--------------
    //2. Pure TFS power
    //--------------
    tic();
    xdc3.tfs_from_ofs(xdc2);
    for(int k = 0; k<= xdc4.get_order(); k++) xdc4.tfts_pows(xdc3, alpha, k);
    xdc.tfs_to_ofs(xdc4);
    cout << "Pure Pows in : " << toc() << "s." << endl;



    //Set results
    //--------------
    xdc2.evaluate(Xvard, ofsd);
    resd1 = cpow(ofsd.evaluate(VAR), alpha);
    xdc.evaluate(Xvard, ofsd);
    resd2 = ofsd.evaluate(VAR);


    //Comparison. cdouble case
    //--------------
    cout << "2. cdouble case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;
}
