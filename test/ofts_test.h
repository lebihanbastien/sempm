#ifndef OFTS_TEST_H_INCLUDED
#define OFTS_TEST_H_INCLUDED

/**
 * \file ofts_test.h
 * \brief Header test file for the Ofts (Fourier-Taylor series) class.
 * \author BLB.
 * \date May 2015
 * \version 1.0
 */

//std
#include <iostream>
#include <fstream>
#include <time.h>
#include <stdio.h>
#include <complex.h>
#include <sstream>

//custom
#include "ofts.h"



/**
 * \fn void ofts_test()
 * \brief Main routine for Ofts class test.
 */
void ofts_test();

/**
 * \fn void ofts_test_conjugate()
 * \brief Test of the routine: Ofts<T>::conjugate().
 */
void ofts_test_conjugate();

/**
 * \fn void ofts_test_smult_t()
 * \brief Test of the routine: Ofts<T>::ofts_smult_t(Ofts<T> const& a, T const& m).
 */
void ofts_test_smult_t();

/**
 * \fn void ofts_test_smult_u()
 * \brief Test of the routine: Ofts<T>::ofts_smult_u(Ofts< Ofs<U> > const& a, U const& c).
 */
void ofts_test_smult_u();

/**
 * \fn void ofts_test_smult_tu()
 * \brief Test of the routine: Ofts<T>::ofts_mult_tu(Ofts< Ofs<U> > const& a, T const& m, U const& c).
 */
void ofts_test_smult_tu();

/**
 * \fn void ofts_test_sfsum_t()
 * \brief Test of the routine: Ofts<T>::ofts_sfsum_t(Ofts<T> const& a, T const& m).
 */
void ofts_test_sfsum_t();

/**
 * \fn void ofts_test_smprod_t()
 * \brief Test of the routine: Ofts<T>::ofts_smprod_t(Ofts<T> const& a, Ofts<T> const& b, T const& m).
 */
void ofts_test_smprod_t();

/**
 * \fn void ofts_test_smprod_u()
 * \brief Test of the routine: Ofts<T>::ofts_smprod_u(Ofts<T> const& a, Ofts<T> const& b, U const& c).
 */
void ofts_test_smprod_u();

/**
 * \fn void ofts_test_sprod()
 * \brief Test of the routine: Ofts<T>::ofts_sprod(Ofts<T> const& a, Ofts<T> const& b).
 */
void ofts_test_sprod();

/**
 * \fn void ofts_test_pows()
 * \brief Test of the routine: Ofts<T>& Ofts<T>::pows(Ofts<T> const& a,  U const& alpha) .
 */
void ofts_test_pows();
#endif // OFTS_TEST_H_INCLUDED
