#ifndef OFTSH_TEST_H_INCLUDED
#define OFTSH_TEST_H_INCLUDED

/**
 * \file oftsh_test.h
 * \brief Header test file for the Oftsh (Homogeneous Fourier-Taylor polynomial) class.
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
 * \fn void oftsh_test()
 * \brief Main routine for Oftsh class test.
 */
void oftsh_test();

/**
 * \fn void oftsh_test_conjugate()
 * \brief Test of the routine: Oftsh<T>::conjugate() and its specializations.
 */
void oftsh_test_conjugate();

/**
 * \fn void oftsh_test_conjugate()
 * \brief Test of the routine: Oftsh<T>& Oftsh<T>::oftsh_smult_t(Oftsh<T> const& a, T const& c) (and mult)
 */
void oftsh_test_smult_t();

/**
 * \fn void oftsh_test_smult_tt()
 * \brief Test of the routine: Oftsh<T>& oftsh_smult_tt(Oftsh<T> const& a, T const& c1, T const& c2).
 */
void oftsh_test_smult_tt();

/**
 * \fn void oftsh_test_smult_u()
 * \brief Test of the routine: Oftsh<T>& Oftsh<T>::oftsh_smult_t(Oftsh<T> const& a, T const& c) (and mult)
 */
void oftsh_test_smult_u();

/**
 * \fn void oftsh_test_sprod()
 * \brief Test of the routine: Oftsh<T>::oftsh_sprod(Oftsh<T> const& a, Oftsh<T> const& b) and its specializations.
 */
void oftsh_test_sprod();

/**
 * \fn void oftsh_test_smprod()
 * \brief Test of the routine: Oftsh<T>::oftsh_smprod_u(Oftsh< Ofs<U> > const& a, Oftsh< Ofs<U> > const& b, U const& m).
 */
void oftsh_test_smprod();

/**
 * \fn void oftsh_test_sderh()
 * \brief Test of the routine: Oftsh<T>& Oftsh<T>::sderh(Oftsh< T > const& a, int ni) and its specializations.
 */
void oftsh_test_sderh();

/**
 * \fn void oftsh_test_derh()
 * \brief Test of the routine: Oftsh<T>& Oftsh<T>::derh(Oftsh< T > const& a, int ni) and its specializations.
 */
void oftsh_test_derh();

/**
 * \fn void oftsh_test_dot()
 * \brief Test of the routine: Oftsh<T>& Oftsh<T>::dot(Oftsh<T> const& a, double const&  n).
 */
void oftsh_test_dot();

#endif // OFTSH_TEST_H_INCLUDED
