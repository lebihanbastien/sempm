#ifndef OFS_TEST_H_INCLUDED
#define OFS_TEST_H_INCLUDED

/**
 * \file ofs_test.h
 * \brief Header test file for the Ofs (Fourier series) class.
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
#include "ofs.h"
#include "timec.h"



/**
 * \fn void ofs_test()
 * \brief Main routine for Ofs class test.
 */
void ofs_test();

//---------------------

/**
 * \fn void ofs_test_conjugate()
 * \brief Test of the routine: void conjugate().
 */
void ofs_test_conjugate();
/**
 * \fn void ofs_test_operators()
 * \brief Test of the operator routines: +=, -=, *= and /=
 */
void ofs_test_operators();
/**
 * \fn void ofs_test_sprod()
 * \brief Test of the routine: Ofs<T>&  sprod(Ofs<T> const& a, Ofs<T> const& b)
 */
void ofs_test_sprod();
/**
 * \fn void ofs_test_prod()
 * \brief Test of the routine: Ofs<T>&  prod(Ofs<T> const& a, Ofs<T> const& b)
 */
void ofs_test_prod();
/**
 * \fn void ofs_test_smprod()
 * \brief Test of the routine: Ofs<T>&  smprod(Ofs<T> const& a, Ofs<T> const& b)
 */
void ofs_test_smprod();
/**
 * \fn void ofs_test_smprod_t()
 * \brief Test of the routine: Ofs<T>&  ofs_smprod_t(Ofs<T> const& a, Ofs<T> const& b, Ofs<T> const& c)
 */
void ofs_test_smprod_t();
/**
 * \fn void ofs_test_mprod()
 * \brief Test of the routine: Ofs<T>&  mprod(Ofs<T> const& a, Ofs<T> const& b)
 */
void ofs_test_mprod();
/**
 * \fn void ofs_smult()
 * \brief Test of the operator routine: Ofs<T>& smult(Ofs<T> const& a, T const& c)
 */
void ofs_test_smult();
/**
 * \fn void ofs_mult()
 * \brief Test of the operator routine: Ofs<T>& mult(Ofs<T> const& a, T const& c)
 */
void ofs_test_mult();
/**
 * \fn void ofs_fsum()
 * \brief Test of the operator routine: Ofs<T>& fsum(Ofs<T> const& a, T const& ma, Ofs<T> const& b, T const& mb)
 */
void ofs_test_fsum();
/**
 * \fn void ofs_test_dot()
 * \brief Test of the routine: void dot().
 */
void ofs_test_dot();
/**
 * \fn void ofs_test_epspow()
 * \brief Test of the routine: Ofs<T>& epspow(Ofs<T> const& a, T const& alpha).
 */
void ofs_test_epspow();

/**
 * \fn void ofs_test_pows()
 * \brief Test of the routine: Ofs<T>& ofs_pows(Ofs<T> const& a, T const& alpha, Ofs<T> const& temp).
 */
void ofs_test_pows();
//---------------------

#endif // OFS_TEST_H_INCLUDED
