#ifndef FTDA_H_INCLUDED
#define FTDA_H_INCLUDED

/**
 * \file manip.h
 * \brief Basic operations for manipulation of polynomials
 *       (exponents in reverse lexicographic order, number of coefficient in
 *       homogeneous polynomials...).
 *       Note that the algebraic operations on series (sum, product) are not done here
 *       but in ofts.h/ofts.tpp.
 * \author BLB using code by Angel Jorba.
 *
 *      Based on Jorba 1999 (http://www.maia.ub.es/~angel/soft.html).
 */

#include <iostream>
#include <iomanip>
#include <limits.h>
#include <cstdio>
#include <cstdlib>
#include <math.h>


/**
 *  \brief Basic operations for manipulation of polynomials.
 *         Manip::m_num_mon(i,j) for i=2..m_nvar, j=1...m_max_deg, is the number of monomials of
 *         degree j with i variables.
 *         Note that the algebraic operations on series (sum, product) are not done here
 *         but in ofts.h/ofts.tpp.
 **/
class Manip;
class Manip
{
private:
    ///static maximum degree
    static int m_max_deg;
    ///static number of variables
    static int m_nvar;
    ///contains m_num_mon(i,j) for i=2..m_nvar, which is the number of monomials of degree j with i variables.
    static long  int **m_num_mon;

public:

    /**
     *   \brief Initializes the table m_num_mon(i,j), which contains the number of
     *          homogeneous polynomials of order j with i variables.
     *
     *   parameters:
     *   t_max_deg: maximum degree we are going to work with. it can not be greater than 63.
     *   t_nvar: number of variables.
     *
     *   returned value: number of kbytes allocated by the internal tables.
     *   Based on a routine by Angel Jorba, 1999.
     **/
    static int init(int t_nvar, int t_max_deg);
    /**
     *  \brief Frees the space allocated by Manip::init.
     **/
    static void free();
    /**
     *   \brief Returns the number of monomials of degree t_max_deg with t_nvar variables,
     *          making use of the table Manip::m_num_mon.
     *
     *  parameters:
     *  t_nvar: number of variables
     *  t_max_deg: order we are interested in (input).
     *  returned value: number of monomials of order no.
     **/
    static long int nmon(int t_nvar, int t_max_deg);
    /**
     *  \brief given a multiindex k, this routine computes the next one
     *  according to the (reverse) lexicographic order.
     *
     *  parameters:
     *  k: array of t_nvar components containing the multiindex. It is overwritten on exit
     *   (input and output).
     **/
    static void prxkt(int k[], int t_nvar);
};


//----------------------------------------------------------------------------------------
//Binomial coefficients. These routines do not make use of Manip::m_num_mon.
//----------------------------------------------------------------------------------------
/**
 * \brief Computes the binomial coefficient (x y).
 **/
unsigned long binomial(unsigned long n, unsigned long k);

/**
 * \brief Computes the binomial coefficient (n k).
 **/
unsigned long gcd_ui(unsigned long x, unsigned long y);

#endif // FTDA_H_INCLUDED
