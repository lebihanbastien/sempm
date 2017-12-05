#ifndef OFS_H_INCLUDED
#define OFS_H_INCLUDED

/**
 * \file ofs.h
 * \brief Fourier series template class
 * \author BLB
 */

//std
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//GSL
#include <gsl/gsl_vector.h>
#include <gsl/gsl_fft_complex.h>
#include "gslc.h"

extern "C"
{
    void foun_(double *F, int *N, int *M, double *CSF, double *SIF);
}

using namespace std;

//----------------------------------------------------------------------------------------
// OFS class
//----------------------------------------------------------------------------------------
template<typename T>
class Ofs;

/**
 * \typedef Ofsc
 * \brief Ofs<double complex> class.
 **/
typedef Ofs<cdouble> Ofsc;
/**
 * \typedef Ofsd
 * \brief Ofs<double> class.
 **/
typedef Ofs<double> Ofsd;

template<typename T>
std::ostream& operator << (std::ostream& stream, Ofs<T> const& ofs);

/**
 *\class Ofs
 *  \brief Fourier series template class.
 *
 *  The Ofs class is an implementation of Fourier series with one pulsation omega:
 *
 *              w(t) = sum_j=-order^+order w_j exp(i j omega t). (1)
 *
 *  In practice, however, the Fourier series are used in two different formats:
 *              - In frequency domain (called ofs format), with the form (1).
 *  In ofs format, the array containing the coefficients w_j is indexed from 0 to 2*order
 *  of the Fourier series, with the positions 0 to order-1 containing the terms of order
 *  -order to -1, and the positions order to 2*order containing the terms of order 0 to
 *  order.
 *              - In time domain (called tfs format), where the 2*order+1 coefficients
 *  contain the evaluations of the series on a time grid over [0 2*pi].
 *
 *  Both format are linked with direct and inverse FFT using the GSL C library.
 *
 *  This class contains an algebra that makes use of both formats. In particular,
 *  the product of two Fourier series is made in tfs format.
 *  On the efficiency of this approach, see e.g.
 *  G. Huguet, R. De La Llave, and Y. Sire (2009). Fast numerical algorithms for the
 *  computation of invariant tori in hamiltonian systems. Prepublicacions del Centre
 *  de Recerca Matemàtica.
 *
 *  Ofs is a templated class. In practice, it is used in the forms
 *          - Ofs<double> (Fourier series with real coefficients, denoted Ofsd)
 *          - Ofs<double complex> (Fourier series with complex coefficients, denoted Ofsc)
 *  The necessary specializations of the templated routines are made with these two forms.
 *
 *
 *  Since Ofs is a templated class, a .tpp is used instead of a .cpp, and this source file
 *  is called directly in the header (the current file), via ' #include "ofs.tpp" '
 *
 */
template<typename T>
class Ofs
{

private:
    int m_ofs_order; // order of the expansion (from -m_ofs_order to +m_ofs_order)
    T  *m_ofs_coef;  // array of coefficients

public:
    //------------------------------------------------------------------------------------
    //Create
    //------------------------------------------------------------------------------------
    /**
     *  \brief Default constructor of the class Ofs<T>.
     */
    Ofs<T>();

    /**
     *  \brief Constructor with a given order.
     *  \param order_: order of the series
     */
    Ofs<T>(const int order_);

    /**
     *  \brief Constructor from a given Ofs object.
     *  \param ofs_:  a reference to the Ofs object to copy in the new object
     */
    Ofs<T>(Ofs const& ofs_);

    //------------------------------------------------------------------------------------
    //Delete
    //------------------------------------------------------------------------------------
    /**
     *  \brief Default destructor of the class Ofs<T>.
     */
    ~Ofs<T>();

    //------------------------------------------------------------------------------------
    //Copy
    //------------------------------------------------------------------------------------
    /**
     *  \brief  Copy from a given Ofs object (only the coefficients).
     *  \param  ofs_: a reference to the Ofs object to copy

     */
    Ofs<T>& ccopy(Ofs<T> const& ofs_);

    /**
     *  \brief  Linked copy from a given Ofs object (exact same object is obtained).
     *  \param  ofs_: a reference to the Ofs object to copy

     */
    Ofs<T>& lcopy(Ofs<T> const& ofs_);

    //------------------------------------------------------------------------------------
    //Setters
    //------------------------------------------------------------------------------------
    /**
     *  \brief Sets a coefficient at a given position in the series.
     *  \param value: the value to set
     *  \param pos: the position (indix in the series) to modify
     */
    void set_coef(T  const&  value, int const& pos);
    template <typename U> void set_coef(U const& value, int const& pos);

    /**
     *  \brief Adds a coefficient at a given position in the series.
     *  \param value: the value to add
     *  \param pos: the position (indix in the series) to modify
     */
    void add_coef(T  const&  value, int const& pos);

    /**
     *  \brief Sets a coefficient to all positions in the series.
     *  \param value: the value to set
     */
    void set_all_coefs(T  const& value);

    /**
     *  \brief Sets random coefficients to all positions in the series.
     */
    void set_random_coefs();

    //------------------------------------------------------------------------------------
    //Getters
    //------------------------------------------------------------------------------------
    /**
     *  \brief  Gets the order of the series.
     *  \return the order of the series as an \c int
     */
    int get_order() const;

    /**
     *  \brief  Gets the pointer address of the Ofs object
     *  \return A pointer to the Ofs object
     */
    Ofs* get_ptr() const;

    /**
     *  \brief  Gets the coefficient at a given position.
     *  \param  pos: the position to get
     *  \return the coefficient of type \c T at the position \c pos
     */
    T ofs_get_coef(int pos) const;

    /**
     *  \brief  Gets the coefficient at a given position.
     *  \param  pos: the position to get
     *  \return the coefficient of type \c T at the position \c pos
     */
    T tfs_get_coef(int pos) const;

    /**
     *  \brief  Computes the maximum coefficient in norm.
     *  \param  maxAbx  a pointer to a \c double array of size 2, so that, for the series
     *  \f$ F_s = \sum \limits_{k = -J}^{J} c_k e^{i k \theta }\f$,
     *   1. <tt>maxAbs[0]</tt> = \f$ \max \limits_{k \in [-J, J]} |c_k|\f$.
     *   2. <tt>maxAbs[1]</tt> = \f$ \arg\max \limits_{k \in [-J, J]} |c_k|\f$.
     */
    void get_coef_max_norm(double maxAbs[]) const; //Get the maximum coefficient (in norm) + its indix

    /**
     *  \brief  Computes the maximum coefficient in norm.
     *  \return a \c double containing <tt>maxAbs[0]</tt> = \f$ \max \limits_{k \in [-J, J]} |c_k|\f$,
     *   for the series \f$ F_s = \sum \limits_{k = -J}^{J} c_k e^{i k \theta }\f$.
     */
    double get_coef_max_norm() const; //get only the maximum in norm

    /**
    *  \brief  Is the Ofs object equal to zero, at order ofs_order?
    */
    bool is_null(const int ofs_order) const;

    //------------------------------------------------------------------------------------
    //Zeroing
    //------------------------------------------------------------------------------------
    /**
     *  \brief  Sets all coefficients to zero.
     */
    void zero();

    //------------------------------------------------------------------------------------
    // Frequency domain <--> Time domain
    //------------------------------------------------------------------------------------
    /**
     *  \brief  From Frequency domain to time domain.
     */
    void tfs_from_ofs(Ofs<T> const& a);
    /**
     *  \brief  Inline from Frequency domain to time domain.
     */
    void tfs_from_ofs_inline(Ofs<T>& temp);
    /**
     *  \brief  From Time domain to Frequency domain.
     */
    void tfs_to_ofs(Ofs<T>const& a);
    /**
     *  \brief  Inline from Time domain to Frequency domain.
     */
    void tfs_to_ofs_inline();
    /**
     *  \brief  From Time domain to Frequency domain. Alternative version with Fortran code, for real values only
     */
    void tfs_to_ofs_F(Ofs<T>& a);
    /**
     *  \brief  From Time domain to Frequency domain. Version with externalized GSL structures.
     */
    void tfs_to_ofs(Ofs<T>& a, gsl_vector_complex *data_fft,
                               gsl_fft_complex_wavetable *wavetable,
                               gsl_fft_complex_workspace *workspace);

    //------------------------------------------------------------------------------------
    // TFS operations
    //------------------------------------------------------------------------------------
    /**
     *  \brief  An operation. Performs the power this = a^alpha in time domain.
     */
    void tfs_pows(Ofs<T> const& a, T const& alpha);

    /**
     *  \brief  An operation. Performs the power this = this^alpha in time domain.
     */
    void tfs_pows(T const& alpha);

    /**
     *  \brief  An operation. Adds the product: \c this \f$ += a \times b \f$ in time domain.
     */
    void tfs_sprod(Ofs<T> const& a, Ofs<T> const& b);

    /**
     *  \brief  An operation. Adds the product: \c this \f$ += a \times b \times c \f$ in time domain.
     */
    void tfs_smprod_t(Ofs<T> const& a, Ofs<T> const& b, Ofs<T> const& c);

    /**
     *  \brief  An operation. Adds the product: \c this \f$ += a \times b \times c \f$ in time domain.
     */
    template<typename U> void tfs_smprod_tu(Ofs<T> const& a, Ofs<T> const& b, Ofs<T> const& c, U const& m);

    /**
     *  \brief  An operation. Adds the product: \c this \f$ += a \times b \times c \f$ in time domain.
     */
    void tfs_smprod(Ofs<T> const& a, Ofs<T> const& b,  T const& m);

    //------------------------------------------------------------------------------------
    // Mixed operations
    //------------------------------------------------------------------------------------
    /**
     *  \brief An operation. Performs the expansion this \f$ = a^{\alpha} \f$ using
     *         inverse FFT, the power function in time domain, and finally direct FFT.
     **/
    void ofs_pows(Ofs<T> const& a, T const& alpha);

    //------------------------------------------------------------------------------------
    //Operators
    //------------------------------------------------------------------------------------
    /**
     *  \brief  An operator. Constructor from a given Ofs object (only the coefficients).
     *  \param  ofs_: a reference to the Ofs object to copy
     */
    Ofs<T>& operator  = (Ofs<T> const& ofs_);

    /**
     *  \brief  An operator. Constructor from a given coefficient. The returned object is of order zero.
     *  \param  coef0: a reference to the coefficient to copy
     */
    Ofs<T>& operator  = (T const& coef0);

    /**
     *  \brief  An operator. Adds all coefficients term by term  from a given Ofs object.
     *  \param  ofs_: the Ofs object to add
     */
    Ofs<T>& operator += (Ofs<T> const& ofs_);

    /**
     *  \brief  An operator. Subtracts all coefficients term by term  from a given Ofs object.
     *  \param  ofs_: the Ofs object to subtract
     */
    Ofs<T>& operator -= (Ofs<T> const& ofs_);

    /**
     *  \brief  An operator. Multiplies all coefficients by a given \c T coefficient.
     *  \param  c: a reference to the multiplicating coefficient.
     */
    Ofs<T>& operator *= (T const& c);

    /**
     *  \brief  An operator. Divides all coefficients by a given \c T coefficient.
     *  \param  c : a reference to the dividing coefficient.
     */
    Ofs<T>& operator /= (T const& c);

    //------------------------------------------------------------------------------------
    //Operations
    //------------------------------------------------------------------------------------
    /**
     *  \brief Conjugates the Ofs object. The conjugate of the series
     *  \f$ F_s = \sum \limits_{k = -J}^{J} c_k e^{i k \theta }\f$ is
     *  \f$ \bar{F}_s = \sum \limits_{k = -J}^{J} \bar{c}_{-k} e^{i k \theta }\f$
     */
    void conjugate();

    /**
     *  \brief  An operation. Adds the product: \c this \f$ += a \times b \f$.
     *  \param  a: a reference to an Ofs object
     *  \param  b: reference to an Ofs object
     */
    void ofs_sprod(Ofs<T> const& a, Ofs<T> const& b);

    /**
     *  \brief  An operation. Adds the product: \c this \f$ += a \times b \times c \f$.
     *  \param  a: a reference to an Ofs object
     *  \param  b: reference to an Ofs object
     *  \param  c: reference to an Ofs object
     */
    void ofs_smprod_t(Ofs<T> const& a, Ofs<T> const& b, Ofs<T> const& c, Ofs<T>& temp);

    /**
     *  \brief  An operation. Sets the product: \c this \f$  = a \times b \f$.
     *  \param  a: a reference to an Ofs object
     *  \param  b: reference to an Ofs object
     */
    void ofs_prod(Ofs<T> const& a, Ofs<T> const& b);

    /**
     *  \brief  An operation. Adds the product: \c this \f$  += m a \times b \f$.
     *  \param  a: a reference to an Ofs object
     *  \param  b: reference to an Ofs object
     *  \param  m: reference to an T coefficient
     */
    void ofs_smprod(Ofs<T> const& a, Ofs<T> const& b, T const& m);

    /**
     *  \brief  An operation. Sets the product: \c this \f$  = m a \times b \f$.
     *  \param  a: a reference to an Ofs object
     *  \param  b: reference to an Ofs object
     *  \param  m: reference to an T coefficient
     */
    void ofs_mprod(Ofs<T> const& a, Ofs<T> const& b, T const& m);

    /**
     *  \brief  An operation. Adds the product: \c this \f$  += m a \f$.
     *  \param  a: a reference to an Ofs object
     *  \param  c: reference to an T coefficient
     */
    void ofs_smult(Ofs<T> const& a, T const& c);

    /**
     *  \brief  An operation. Adds the product: \c this \f$  += m a \f$ at order eff_order.
     *  \param  a: a reference to an Ofs object
     *  \param  c: reference to an T coefficient
     *  \param  eff_order: the effective order of the operation
     */
    void ofs_smult(Ofs<T> const& a, T const& c, int eff_order);

    /**
     *  \brief  An operation. Sets the product: \c this \f$  = m a \f$.
     *  \param  a: a reference to an Ofs object
     *  \param  c: reference to an T coefficient
     */
    void ofs_mult(Ofs<T> const& a, T const& c);

    /**
     *  \brief  An operation. Sets the sum-product: \c this \f$  = m_a a + m_b b \f$.
     *  \param  a: a reference to an Ofs object
     *  \param  ma: reference to an T coefficient
     *  \param  b: a reference to an Ofs object
     *  \param  mb: reference to an T coefficient
     */
    void ofs_fsum(Ofs<T> const& a, T const& ma, Ofs<T> const& b, T const& mb);

    /**
     *  \brief  An operation. Set the time derivative of object \c a with
     *          pulsation \f$ \omega = n \f$, so that \c this \f$ = \dot{a} \f$.
     *  \param  a: a reference to an Ofs object
     *  \param  n: reference to the pulsation
     *
     */
    void dot(Ofs<T> const& a, double const& n);

    /**
     *  \brief  An operation. Applies the time derivative with
     *          pulsation \f$ \omega = n \f$ directly to this.
     *  \param  n: reference to the pulsation
     *
     */
    void dot(double const& n);

    /**
     *  \brief  An operation. Performs the expansion this \f$ = a^{\alpha} \f$
     *          when \c a is of the form \c this \f$ a = 1 + d \f$ with \f$ |d|_1 << 1 \f$.
     *  \param  a: a reference to an Ofs object
     *  \param  alpha: a reference to a T coefficient

     */
    void ofs_epspow(Ofs<T> const& a, T const& alpha);

    /**
     *  \brief  Performs the comparison with Ofs object \c b.
     *  \param  b: a reference to an Ofs object
     *  \return \c true if \c this = \c b (same order and same coefficients
     *          at every position). \c false otherwise.
     */
    bool isEqual(Ofs<T> const& b) const;


    /**
     *  \brief  Evaluates the Ofs object at time t.
     *  \param  t: a reference to a \c double t
     *  \param  eff_order: the effective order of the operation
     *  \return the evaluation \f$ F_s(t) \f$ at time t.
     */
    cdouble evaluate(double const& t, int eff_order);
    cdouble fevaluate(double cR[], double sR[], int eff_order);

    /**
     *  \brief  Evaluates the Ofs object at time t.
     *  \param  t: a reference to a \c double t
     *  \return the evaluation \f$ F_s(t) \f$ at time t.
     */
    cdouble evaluate(double const& t);

    /**
     *  \brief  Evaluates the derivatives of the Ofs object at angle theta,
     *          with pulsation n (theta = nt).
     *
     *  The sum is given in the form:
     *      \f$ c_0 + n \sum \limits_{k = 1}^{J} - k \sin(k\theta) (c_k + c_{-k}) + i k \cos(k\theta) (c_k - c_{-k}) \f$
     *  in order to avoid any additionnal roundoff errors by using cos and sin
     *  functions as little as possible.
     */
    cdouble evaluatedot(double const& theta, double const& n);

    /**
     *  \brief  Evaluates the derivatives of the Ofs object at angle theta,
     *          with pulsation n (theta = nt), at a certain order eff_order.
     *
     *  The sum is given in the form:
     *      \f$ c_0 + n \sum \limits_{k = 1}^{J} - k \sin(k\theta) (c_k + c_{-k}) + i k \cos(k\theta) (c_k - c_{-k}) \f$
     *  in order to avoid any additionnal roundoff errors by using cos and sin
     *  functions as little as possible.
     */
    cdouble evaluatedot(double const& theta, double const& n, int eff_order);

    /**
     *  \brief  Evaluates the \f$ L_1 \f$ norm of the current Ofs object.
     *  \return a \c double containing \f$ \sum \limits_{k = -J}^{J} |c_k| \f$
     *          for the series \f$ F_s = \sum \limits_{k = -J}^{J} c_k e^{i k \theta }\f$.
     */
    double l1norm();

    /**
     *  \brief  Number of small divisors under a certain value
     */
    int nsd(int odmax, double sdmax);

    /**
     *  \brief  A stream operator
     */
    friend std::ostream& operator << <>(std::ostream& stream, Ofs<T> const& ofs);

    /**
     *  \brief  A stream operator. Print only the autonomous term (term of order zero).
     */
    void fprint_0(ofstream& stream);
};

//----------------------------------------------------------------------------------------
//Functions
//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------
//Operators
//----------------------------------------------------------------------------------------
/**
 * \fn template<typename T> bool   operator == (Ofs<T> const& a, Ofs<T> const& b)
 * \brief An operator. Compares two Ofs objects
 * \return \c true if the order are equals and all the coefficients,
 *          term by term. \c false otherwise.
 */
template<typename T> bool   operator == (Ofs<T> const& a, Ofs<T> const& b);

/**
 * \fn template<typename T> Ofs<T> operator - (Ofs<T> const& b)
 * \brief Returns -b at order b.order
 * \param b: a reference to an Ofs object
 * \return -b at order b.order.
 */
template<typename T> Ofs<T> operator - (Ofs<T> const& b);

/**
 * \fn template<typename T> Ofs<T> operator + (Ofs<T> const& a, Ofs<T> const& b)
 * \brief An operator. Makes the sum a+b
 * \param a: a reference to an Ofs object
 * \param b: a reference to an Ofs object
 * \return a+b at order max(a.order, b.order).
 */
template<typename T> Ofs<T> operator + (Ofs<T> const& a, Ofs<T> const& b);

/**
 * \fn template<typename T> Ofs<T> operator - (Ofs<T> const& a, Ofs<T> const& b)
 * \brief An operator. Makes the sum a-b
 * \param a: a reference to an Ofs object
 * \param b: a reference to an Ofs object
 * \return a-b at order max(a.order, b.order).
 */
template<typename T> Ofs<T> operator - (Ofs<T> const& a, Ofs<T> const& b);

/**
 * \fn template<typename T> Ofs<T> operator * (Ofs<T> const& a, T const& c)
 * \brief An operator. Makes the product c*a
 * \param a: a reference to an Ofs object
 * \param c: a reference to a T coefficient
 * \return c*a at order a.order.
 */
template<typename T> Ofs<T> operator * (Ofs<T> const& a, T const& c);

/**
 * \fn template<typename T> Ofs<T> operator / (Ofs<T> const& a, T const& c)
 * \brief An operator. Makes the division a/c
 * \param a: a reference to an Ofs object
 * \param c: a reference to a T coefficient
 * \return a/c at order a.order.
 */
template<typename T> Ofs<T> operator / (Ofs<T> const& a, T const& c);

/**
 * \fn template<typename T> Ofs<T> operator * (T const& c, Ofs<T> const& a)
 * \brief An operator. Makes the product c*a
 * \param a: a reference to an Ofs object
 * \param c: a reference to a T coefficient
 * \return c*a at order a.order.
 */
template<typename T> Ofs<T> operator * (T const& c, Ofs<T> const& a);

/**
 * \fn template<typename T> Ofs<T> operator *  (Ofs<T> const& a, Ofs<T> const& b)
 * \brief An operator. Makes the product a*b
 * \param a: a reference to an Ofs object
 * \param b: a reference to an Ofs object
 * \return a*b at order a.order.
 */
template<typename T> Ofs<T> operator *  (Ofs<T> const& a, Ofs<T> const& b);


//----------------------------------------------------------------------------------------
// I/O
//----------------------------------------------------------------------------------------
/**
 * \fn void inline read_ofs_txt(Ofsc& xFFT, string filename, int fftN)
 * \brief Reading an Ofsc object from a text file.
 * \param xFFT: a reference to the Ofs object to update.
 * \param filename: a \c string containing the path to the txt file
 *       (without the suffix ".txt").
 */
void inline read_ofs_txt(Ofsc& xFFT, string filename);

/**
 * \fn void inline ofs_double_to_complex(Ofsd& xr, Ofsc& xc)
 * \brief Copy from Ofsd to Ofsc.
 * \param xr: a reference to the Ofsd object to copy
 * \param xc: a reference to the Ofsc object to update
 */
void inline ofs_double_to_complex(Ofsd const& xr, Ofsc& xc);

/**
 * \fn void inline ofs_real_part(Ofsc& x, Ofsc& const xr)
 * \brief Takes the real part of an Ofsc object: xr = real(x).
 * \param x: a reference to the Ofsc to use
 * \param xr: a reference to the Ofsc to update
 */
void inline ofs_real_part(Ofsc const& x, Ofsc&  xr);

/**
 * \fn void inline ofs_real_part(Ofsc& x, Ofsc& const xr)
 * \brief Takes the real part of an Ofsc object: xi = imag(x).
 * \param x: a reference to the Ofsc to use
 * \param xi: a reference to the Ofsc to update
 */
void inline ofs_imag_part(Ofsc const& x, Ofsc& xi);

/**
 *  \fn     template<typename T>  cdouble sprod_expct_error(Ofs<T> const& a, Ofs<T> const& b, double const& t)
 *  \brief  Expected error on the product: \c this \f$ += a \times b \f$ at times t. Works only when a.order = b.order which is the default case.
 *  \param  a: a reference to an Ofs object
 *  \param  b: reference to an Ofs object
 *  \param  t: a reference to a \c double t
 *  \return the error in the product
 */
template<typename T>  cdouble sprod_expct_error(Ofs<T> const& a, Ofs<T> const& b, double const& t);

/**
 *  \brief  Expected error on the product: \c this \f$ += a \times b \f$ at times t.
 *          Works only when a.order = b.order which is the default case.
 *  \param  a: a reference to an Ofs object
 *  \param  b: reference to an Ofs object
 *  \param  c: reference to an T coefficient
 *  \param  t: a reference to a \c double t
 *  \return the error in the product
 */
template<typename T>  cdouble smprod_expct_error(Ofs<T> const& a, Ofs<T> const& b, T const& c,  double const& t);

/**
 *  \brief Initializes the arrays cR[] and sR[], containing the numbers
 *         cos(t), cos(nt), ... cos(ofs_order*n*t) and
 *         sin(t), sin(nt), ... sin(ofs_order*n*t), respectively.
 **/
void init_cR_sR(double t, double cR[], double sR[], int ofs_order);

//----------------------------------------------------------------------------------------
// Single storage
//----------------------------------------------------------------------------------------
/**
 *  \brief Single storage of a QBTBP Ofs object in txt files + its cosinus/sinus version.
 *  \param xOFS: a reference to the Ofs object to update.
 *  \param filename: a \c string containing the txt file to write, without the suffix ".txt".
 *  \param flag: if flag==1, the function to store is even, otherwise it is odd.
 */
void ofs_sst(Ofsc &xOFS, string filename, int flag, string suffix);

//Include the implementation .tpp
#include "ofs.tpp"

//----------------------------------------------------------------------------------------
// end of OFS class
//----------------------------------------------------------------------------------------

#endif // OFS_H_INCLUDED
