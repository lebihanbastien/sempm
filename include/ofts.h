#ifndef OFTS_H_INCLUDED
#define OFTS_H_INCLUDED

/**
 * \file ofts.h
 * \brief Fourier-Taylor series template class
 * \author BLB
 */

//std
#include <iostream>
#include <iomanip>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <vector>
//custom
#include "manip.h"
#include "oftsh.h"
#include "ofs.h"
#include "timec.h"
//gsl
#include <gsl/gsl_errno.h>

using namespace std;

template <typename T>
class Ofts;

typedef Ofts<Ofsc> Oftsc;
typedef Ofts<Ofsd> Oftsd;

template<typename T>
std::ostream& operator << (std::ostream& stream, Ofts<T> const& ofts);

/** \class Ofts
 *  \brief Fourier-Taylor series class.
 *
 *  The Ofts class implements Fourier-Taylor series where each homogeneous polynomial is
 *  an Oftsh object.
 *
 *   *  Ofts is a templated class. In practice, it is used in the forms
 *          - Ofts<Ofs<double>> i.e. Fourier-Taylor series with Fourier series of
 *            real coefficients (Ofsd) as coefficients. They are denoted Oftsd.
 *          - Ofts<Ofs<double complex>> i.e. Fourier-Taylor series with Fourier series of
 *            complex coefficients (Ofsc) as coefficients. They are denoted Oftsc.
 *
 *  Since the Ofs class implements Fourier series in two formats (frequency domain - ofs -
 *  and time domain - tfs), the Ofts also contain operation routines that works with these
 *  two formats. The routines that must be used with ofs format have the prefix "oftsh_",
 *  while the routines that must be used with tfs format have the prefix "tftsh_".
 *
 *  Since Ofts is a templated class, a .tpp is used instead of a .cpp, and this source
 *  file is called directly in the header (the current file), via ' #include "ofts.tpp" '
 *
 *
 *  Finally, an important note: the delete of Ofts objects is badly handled because
 *  There is no proper destructor for the Oftsh objects.
 *  As a consequence, the Ofts object must be handled with care (do not create too many!).
 */
template <typename T>
class Ofts
{

private:

    //number of variables
    int m_ofts_nvar;
    //order
    int m_ofts_order;
    //number of variables of the coefficients
    int m_ofs_nvar;
    //order of the coefficients
    int corder;
    //Array pointing to the homogeneous polynomials
    Oftsh<T> **m_ofts_term;
    //Array pointing to the coefficients
    T *m_ofts_coefs;

public:

    //------------------------------------------------------------------------------------
    //Create
    //------------------------------------------------------------------------------------
    /**
     *  \brief Default constructor of the class Ofts<T>.
     */
    Ofts<T>();

    /**
     *  \brief Constructor with given order/number of variables both for the Fourier-Taylor series and the coefficients.
     *  \param newNv: number of variables of the series
     *  \param newOrder: order of the series
     *  \param newCnv: number of variables of the coefficients
     *  \param newOrder: order of the coefficients
     */
    Ofts<T>(int newNv, int newOrder, int newCnv, int newCorder);

    /**
     *  \brief Constructor with given order/number of variables both for the Fourier-Taylor series and the coefficients.
     *   The coefficients to set are given in an array.
     *  \param newNv: number of variables of the series
     *  \param newOrder: order of the series
     *  \param newCnv: number of variables of the coefficients
     *  \param newOrder: order of the coefficients
     *  \param coef0: the array of coefficients to set
     */
    Ofts<T>(int newNv, int newOrder, int newCnv, int newCorder, T *coef0);

    /**
     *  \brief Constructor from a given Ofts object (without any link).
     *  \param b:  a reference to the Ofts object to copy in the new object
     */
    Ofts<T>(Ofts<T> const& b);

    //------------------------------------------------------------------------------------
    //Delete
    //------------------------------------------------------------------------------------
    /**
     *  \brief Default destructor of the class Ofts<T>. WARNING: potential memory leak here, through the terms of type Oftsh.
     */
    ~Ofts<T>();

    //------------------------------------------------------------------------------------
    //Copy
    //------------------------------------------------------------------------------------
    /**
     *  \brief  Linked copy from a given Ofts object (exact same object is obtained).
     *  \param  b: a reference to the Ofts object to copy
     *  \return a reference to the current object
     */
    Ofts<T>& lcopy(Ofts<T> const& b); //linked copy (probably useless at this step)
    /**
     *  \brief  Copy from a given Ofts object (only the coefficients).
     *  \param  b: a reference to the Ofts object to copy
     *  \return a reference to the current object
     */
    Ofts<T>& ccopy(Ofts<T> const& b); //coefficient copy (restricted to same order, same number of variables)
    /**
     *  \brief  Copy from a given Ofts object (only the coefficients) at order nrc.
     *  \param  b: a reference to the Ofts object to copy
     *  \param  nrc: a reference to the order to update
     *  \return a reference to the current object
     */
    Ofts<T>& ccopy(Ofts<T> const& b, int const& nrc); //coefficient copy (restricted to same order, same number of variables)
    /**
     *  \brief  An operator. Constructor from a given Ofts object (only the coefficients).
     *  \param  b: a reference to the Ofts object to copy
     *  \return a reference to the current object
     */
    Ofts<T>& operator = (Ofts<T> const& b); //unlinked copy

    //------------------------------------------------------------------------------------
    //Setters
    //------------------------------------------------------------------------------------
    /**
     *  \brief Sets a coefficient at a given order \c ord and a given position \c i at this order in the series.
     *  \param m: the value to set
     *  \param ord: the order to modify
     *  \param i: the position to modify
     */
    void set_coef(T const& m, int ord, int i);

    /**
     *  \brief Adds a coefficient at a given order \c ord and a given position \c i at this order in the series.
     *  \param m: the value to set
     *  \param ord: the order to modify
     *  \param i: the position to modify
     */
    void add_coef(T const& m, int ord, int i);

    /**
     *  \brief Sets of a given double/complex (typename U) subcoefficient everywhere (at each order and in each coefficient).
     *  \param m: the value to set
     */
    template <typename U> void set_all_coefs(U const& m);

    /**
     *  \brief Sets of a given U subcoefficient at order zero of the coefficient at position \c i of term of order \c n.
     *  \param m: the value to set
     *  \param ord: the order to modify
     *  \param i: the position to modify
     */
    template < typename U > void set_coef0(U const& m, int const& ord, int const& i);

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
     *  \brief  Gets the number of variables of series.
     *  \return the number of variables of the series as an \c int
     */
    int get_nvar() const;

    /**
     *  \brief  Gets the order of the coefficients.
     *  \return the order of the coefficients as an \c int
     */
    int get_coef_order() const;

    /**
     *  \brief  Gets the number of variables of the coefficients.
     *  \return the number of variables of the coefficients as an \c int
     */
    int get_coef_nvar() const;

    /**
     *  \brief  Gets the adress of the coefficient at order \c ord and position \c pos
     *  \return a pointer to the desired coefficient
     */
    T* get_coef(int const& ord, int const& pos) const;

    /**
     *  \brief  Gets the adress of the term at order \c ord
     *  \return a pointer to the desired term
     */
    Oftsh<T>* get_term(int const& ord) const;

    /**
     *  \brief  Gets the adress of the Ofts object
     *  \return a pointer to the Ofts object
     */
    Ofts<T>* get_ptr() const;


    //------------------------------------------------------------------------------------
    //Zeroing
    //------------------------------------------------------------------------------------
    /**
     *  \brief  Sets all coefficients to zero.
     */
    void zero();

    //------------------------------------------------------------------------------------
    //Operations
    //------------------------------------------------------------------------------------
    //------------------
    // Conjugate
    //------------------
    /**
     *  \brief Conjugates  all terms (Oftsh object), and only them! To be used with evaluate_conjugate to have the true conjugate.
     */
    Ofts<T>& conjugate();

    /**
     *  \brief Conjugates the order \c nrc. To be used with evaluate_conjugate to have the true conjugate.
     */
    Ofts<T>& conjugate(int const& nrc);

    //------------------
    // smult
    //------------------
    /**
     *  \brief  An operation. Adds the product: \c this \f$  += m a \f$ with m a coefficient.
     *  \param  a: a reference to an Ofts object
     *  \param  m: a reference to a coefficient
     */
    Ofts<T>& ofts_smult_t(Ofts<T> const& a, T const& m);

    //------------------

    /**
     *  \brief  An operation. Adds the product: \c this \f$  += m a \f$ with m a coefficient, at order n
     *  \param  a: a reference to an Ofts object
     *  \param  m: a reference to a coefficient
     *  \param  n: a reference to the order to update
     */
    Ofts<T>& ofts_smult_t(Ofts<T> const& a, T const& m, int const& n);

    /**
    *  \brief  An operation. Sets the product: \c this \f$  = m a \f$ with m a coefficient, at order n
    *  \param  a: a reference to an Oftsh object
    *  \param  m: a reference to a coefficient
    *  \param  n: a reference to the order to update
    */
    Ofts<T>& ofts_mult_t(Ofts<T> const& a, T const& m, int const& n);

    /**
     *  \brief  An operation. Adds the product: \c this \f$  += c a \f$ with c a subcoefficient.
     *  \param  a: a reference to an Ofts object
     *  \param  c: a reference to a subcoefficient
     */
    template<typename U> Ofts< Ofs<U> >& ofts_smult_u(Ofts< Ofs<U> > const& a, U const& c);

    /**
     *  \brief  An operation. Adds the product: \c this \f$  += c a \f$ with c a subcoefficient at order k.
     *  \param  a: a reference to an Ofts object
     *  \param  c: a reference to a subcoefficient
     *  \param  k: a reference to the order to update
     */
    template<typename U> Ofts< Ofs<U> >& ofts_smult_u(Ofts< Ofs<U> > const& a, U const& c, int const& k);

    /**
     *  \brief  An operation. Adds the product: \c this \f$  += c m a \f$ with m a coefficient and c a subcoefficient.
     *  \param  a: a reference to an Ofts object
     *  \param  c: a reference to a subcoefficient
     */
    template<typename U> Ofts< Ofs<U> >& ofts_smult_tu(Ofts< Ofs<U> > const& a, T const& m, U const& c);
    /**
     *  \brief  An operation. Adds the product: \c this \f$  += c m a \f$ with m a coefficient and c a subcoefficient, at order k.
     *  \param  a: a reference to an Ofts object
     *  \param  c: a reference to a subcoefficient
     *  \param  k: a reference to the order to update
     */
    template<typename U> Ofts< Ofs<U> >& ofts_smult_tu(Ofts< Ofs<U> > const& a, T const& m, U const& c, int const& k);


    //------------------
    // mult
    //------------------
    /**
     *  \brief  An operation. Sets the product: \c this \f$  = m a \f$ with m a coefficient.
     *  \param  a: a reference to an Ofts object
     *  \param  m: a reference to a coefficient
     */
    Ofts<T>& ofts_mult_t(Ofts<T> const& a, T const& m);

    //------------------

    /**
     *  \brief  An operation. Sets the product: \c this \f$  = c a \f$ with c a subcoefficient.
     *  \param  a: a reference to an Ofts object
     *  \param  c: a reference to a subcoefficient
     */
    template<typename U> Ofts< Ofs<U> >&  ofts_mult_u(Ofts< Ofs<U> > const& a, U const& c);

    /**
     *  \brief  An operation. Sets the product: \c this \f$  = c a \f$ with c a subcoefficient at order k.
     *  \param  a: a reference to an Ofts object
     *  \param  c: a reference to a subcoefficient
     *  \param  k: a reference to the order to update
     */
    template<typename U> Ofts< Ofs<U> >&  ofts_mult_u(Ofts< Ofs<U> > const& a, U const& c, int const& k);

    /**
     *  \brief  An operation. Sets the product: \c this \f$  = c m a \f$ with m a coefficient and c a subcoefficient.
     *  \param  a: a reference to an Ofts object
     *  \param  c: a reference to a subcoefficient
     */
    template<typename U> Ofts< Ofs<U> >&  ofts_mult_tu(Ofts< Ofs<U> > const& a, T const& m, U const& c);

    //------------------
    // sfsum
    //------------------
    /**
     *  \brief  An operation. Adds the product: \c this \f$  += ma*a + mb*b \f$ with ma and mb coefficients.
     *  \param  a: a reference to an Ofts object
     *  \param  ma: a reference to a coefficient
     *  \param  b: a reference to an Ofts object
     *  \param  mb: a reference to a coefficient
     */
    Ofts<T>& ofts_sfsum_t(Ofts<T> const& a,  T const& ma, Ofts<T> const& b, T const& mb);
    /**
     *  \brief  An operation. Adds the product: \c this \f$  += ma*a + mb*b \f$ with ma and mb coefficients at order n.
     *  \param  a: a reference to an Ofts object
     *  \param  ma: a reference to a coefficient
     *  \param  b: a reference to an Ofts object
     *  \param  mb: a reference to a coefficient
     *  \param  n: a reference to the order to update
     */
    Ofts<T>& ofts_sfsum_t(Ofts<T> const& a,  T const& ma, Ofts<T> const& b, T const& mb, int n);
    /**
     *  \brief  An operation. Adds the product: \c this \f$  += ma*a + mb*b + mc*c \f$ with ma and mb coefficients at order n.
     *  \param  a: a reference to an Ofts object
     *  \param  ma: a reference to a coefficient
     *  \param  b: a reference to an Ofts object
     *  \param  mb: a reference to a coefficient
     *  \param  c: a reference to an Ofts object
     *  \param  mc: a reference to a coefficient
     *  \param  n: a reference to the order to update
     */
    Ofts<T>& ofts_sfsum_tt(Ofts<T> const& a,  T const& ma, Ofts<T> const& b, T const& mb, Ofts<T> const& c, T const& mc, int n);

    //------------------
    // fsum
    //------------------
    /**
     *  \brief  An operation. Sets the product: \c this \f$  = ma*a + mb*b \f$ with ma and mb coefficients.
     *  \param  a: a reference to an Ofts object
     *  \param  ma: a reference to a coefficient
     *  \param  b: a reference to an Ofts object
     *  \param  mb: a reference to a coefficient
     */
    Ofts<T>&  ofts_fsum_t(Ofts<T> const& a,  T const& ma, Ofts<T> const& b, T const& mb);

    /**
     *  \brief  An operation. Sets the product: \c this \f$  = ma*a + mb*b \f$ with ma and mb coefficients at order n.
     *  \param  a: a reference to an Ofts object
     *  \param  ma: a reference to a coefficient
     *  \param  b: a reference to an Ofts object
     *  \param  mb: a reference to a coefficient
     *  \param  n: a reference to the order to update
     */
    Ofts<T>&  ofts_fsum_t(Ofts<T> const& a,  T const& ma, Ofts<T> const& b, T const& mb, int n);

    /**
     *  \brief  An operation. Sets the product: \c this \f$  = ca*a + cb*b \f$ with ca and cb subcoefficients.
     *  \param  a: a reference to an Ofts object
     *  \param  ca: a reference to a subcoefficient
     *  \param  b: a reference to an Ofts object
     *  \param  cb: a reference to a subcoefficient
     */
    template<typename U> Ofts<T>&  ofts_fsum_u(Ofts<T> const& a,  U const& ca, Ofts<T> const& b, U const& cb);

    /**
     *  \brief  An operation. Sets the product: \c this \f$  = ca*a + cb*b \f$ with ca and cb subcoefficients at order m
     *  \param  a: a reference to an Ofts object
     *  \param  ca: a reference to a subcoefficient
     *  \param  b: a reference to an Ofts object
     *  \param  cb: a reference to a subcoefficient
     *  \param  m: a reference to the order to update
     */
    template<typename U> Ofts<T>&  ofts_fsum_u(Ofts<T> const& a,  U const& ca, Ofts<T> const& b, U const& cb, int const& m);


    //------------------
    // sprod
    //------------------
    /**
     *  \brief  An operation. Adds the product: \c this \f$  += a*b \f$.
     *  \param  a: a reference to an Ofts object
     *  \param  b: a reference to an Ofts object
     */
    Ofts<T>& ofts_sprod(Ofts<T> const& a, Ofts<T> const& b);

    /**
     *  \brief  An operation. Adds the product: \c this \f$  += a*b \f$ at order n.
     *  \param  a: a reference to an Ofts object
     *  \param  b: a reference to an Ofts object
     *  \param  n: a reference to the order to update
     */
    Ofts<T>& ofts_sprod(Ofts<T> const& a, Ofts<T> const& b, int const& n);

    //------------------
    // smprod
    //------------------
    /**
     *  \brief  An operation. Adds the product: \c this \f$  += m*a*b \f$ with m a coefficient.
     *  \param  a: a reference to an Ofts object
     *  \param  b: a reference to an Ofts object
     *  \param  m: a reference to a coefficient
     */
    Ofts<T>& ofts_smprod_t(Ofts<T> const& a, Ofts<T> const& b, T const& m, T& temp);

    /**
     *  \brief  An operation. Adds the product: \c this \f$  += m*a*b \f$ with m a coefficient at order n.
     *  \param  a: a reference to an Ofts object
     *  \param  b: a reference to an Ofts object
     *  \param  m: a reference to a coefficient
     */
    Ofts<T>& ofts_smprod_t(Ofts<T> const& a, Ofts<T> const& b, T const& m, int const& n, T& temp);

    /**
     *  \brief  An operation. Adds the product: \c this \f$  += m*a*b \f$ with c a subcoefficient.
     *  \param  a: a reference to an Ofts object
     *  \param  b: a reference to an Ofts object
     *  \param  c: a reference to a subcoefficient
     */
    template<typename U> Ofts<T>& ofts_smprod_u(Ofts<T> const& a, Ofts<T> const& b, U const& c);

    /**
     *  \brief  An operation. Adds the product: \c this \f$  += m*a*b \f$ with c a subcoefficient at order n.
     *  \param  a: a reference to an Ofts object
     *  \param  b: a reference to an Ofts object
     *  \param  c: a reference to a subcoefficient
     *  \param  n: a reference to the order to update
     */
    template<typename U> Ofts<T>& ofts_smprod_u(Ofts<T> const& a, Ofts<T> const& b, U const& c, int const& n);


    //------------------
    // prod
    //------------------
    /**
     *  \brief  An operation. Sets the product: \c this \f$  = a*b \f$.
     *  \param  a: a reference to an Ofts object
     *  \param  b: a reference to an Ofts object
     */
    Ofts<T>& prod(Ofts<T> const& a, Ofts<T> const& b);
    /**
     *  \brief  An operation. Sets the product: \c this \f$  = m*a*b \f$ with m a coefficient.
     *  \param  a: a reference to an Ofts object
     *  \param  b: a reference to an Ofts object
     *  \param  m: a reference to a coefficient
     */
    Ofts<T>& mprod(Ofts<T> const& a, Ofts<T> const& b, T const& m);
    /**
     *  \brief  An operation. Sets the product: \c this \f$  = m*a*b \f$ with c a subcoefficient.
     *  \param  a: a reference to an Ofts object
     *  \param  b: a reference to an Ofts object
     *  \param  c: a reference to a subcoefficient
     */
    template<typename U> Ofts<T>&  mprod(Ofts<T> const& a, Ofts<T> const& b, U const& c);


    //------------------
    // pows
    //------------------
    /**
     *   \brief Power function: p = a^alpha at order n. This generic routine is not directly used, only specialized versions are (in particular for the Ofts<Ofs> form).
     **/
    template<typename U> Ofts<T>& ofts_pows(Ofts<T> const& a,  U const& alpha, int const& n);

    /**
     *  \brief Power function: p = a^alpha at order n, for Oftsd objects with Fourier coefficients in Frequency domain.
     *         Moreover, the order zero of the FT series a[0] must satisfy: a[0] >> a[i], for all i > 0
     *         In this routine, the coefficient a0inv = 1/(a[0]) and a0palpha = a[0]^alpha
     **/
    template<typename U> Ofts<T>& pows(Ofts< Ofsd > const& a,  //intiial Ofts
                                       Ofsd a0inv,             //inverse of order 0
                                       Ofsd a0palpha,          //order 0 ^alpha
                                       U const& alpha);        //power coef

    //------------------
    // Order 0 routines
    //------------------
    /**
     * \brief Returns the address of the first coefficient of order 0 of the taylor series s
     */
    T* coef0s(Ofts<T> const& a);

    /**
     * \brief Sets the coefficient of order 0 of the taylor series s equal to x0
     */
    void acoef0s(T const& x0);

    //------------------------------------------------------------------------------------
    //Derivation
    //------------------------------------------------------------------------------------
    /**
     *  \brief Partial derivative wrt to the variable z[ni], with \f$ n_i \in [[1, n]] \f$: \c this \f$ = \frac{\partial a}{\partial z_{n_i}} \f$.
     **/
    Ofts<T>& der(Ofts<T> const& a, int ni);

    /**
     *  \brief Partial derivative wrt to the variable z[ni], with \f$ n_i \in [[1, n]] \f$: \c this \f$ += \frac{\partial a}{\partial z_{n_i}} \f$.
     **/
    Ofts<T>& sder(Ofts< T > const& a, int ni);

    /**
     *  \brief Partial derivative wrt to the variable z[ni], with \f$ n_i \in [[1, n]] \f$, at order m of the expansions:
     *          \f$ this_{m-1}  = \left[ \frac{\partial a}{\partial z_{n_i}} \right]_m \f$.
     *
     *          Note: order m is derived and set in order m-1 of this. If m==0, nothing is done, but a warning is sent to the user.
     **/
    Ofts<T>& der(Ofts< T > const& a, int ni, int m);

    /**
     *  \brief Partial derivative wrt to the variable z[ni], with \f$ n_i \in [[1, n]] \f$, at order m of the expansions:
     *          \f$ this_{m-1}  += \left[ \frac{\partial a}{\partial z_{n_i}} \right]_m \f$.
     *
     *          Note: order m is derived and set in order m-1 of this. If m==0, nothing is done, but a warning is sent to the user.
     **/
    Ofts<T>& sder(Ofts< T > const& a, int ni, int m);

    /**
     *  \brief Partial derivative wrt to time:
     *          \f$ this  += \left[ \frac{\partial a}{\partial t} \right] \f$.
     **/
    Ofts<T>& dot(Ofts<T> const& a, double const&  n);

    /**
     *  \brief Partial derivative wrt to the time, at order k of the expansions:
     *          \f$ this_{k-1}  += \left[ \frac{\partial a}{\partial t} \right]_k \f$.
     **/
    Ofts<T>& dot(Ofts<T> const& a, double const&  n, int const& k);


    //---------------------------------------------------------------------------
    //Integral
    //---------------------------------------------------------------------------
    /**
     *  \brief Primitive wrt to the variable z[ni], with \f$ n_i \in [[1, n]] \f$: \c this \f$ = \int a(z) dz_{n_i} \f$.
     **/
    Ofts<T>& sprim(Ofts< T > const& a, int ni);


    //--------------------------------------------------------------------------------
    //Evaluate
    //--------------------------------------------------------------------------------
    /**
     *  \brief Generic routine for the evaluation of an Ofts object at the state X: z = this(X).
     **/
    template<typename U> void evaluate(U X[], T& z);

    /**
     *  \brief Generic routine for the conjugate evaluation of an Ofts object at the state X: z = conj(this(X)).
     **/
    template<typename U> void evaluate_conjugate(U X[], T& z);

    /**
     *  \brief Toutine for the evaluation of an Ofts object at the state X, at order m: z = [this(X)]_m.
     **/
    template<typename U> void evaluate(U X[], T& z, int const& m, int const& ofs_order);

    /**
     *  \brief Routine for the evaluation of an Ofts object at the state X, at order m, and at time t: z = [this(X, t)]_(<=m).
     **/
    template<typename U> cdouble fevaluate(U X[], double const& t, int const& m, int const& ofs_order);

    /**
     *  \brief Routine for the evaluation of an Ofts object at the state X, at order m, and at time t: z = [this(X, t)]_(<=m).
     *
     *         Contrary to the routine fevaluate(U X[], double const& t, int const& m, int const& ofs_order), the cosinus/sinus arrays are given as inputs:
     *              - cR[] = [cos(t), ..., cos(ofs_order*t)]
     *              - sR[] = [sin(t), ..., sin(ofs_order*t)]
     **/
    template<typename U> cdouble fevaluate(U X[], double cR[], double sR[], int const& m, int const& ofs_order);

    /**
     *  \brief Contribution of the order m of this to the evaluation: z = [this(X)]_=m).
     **/
    template<typename U> void contribution(U X[], T& z, int const& m);


    //---------------------------------------------------------------------------
    //Norms
    //---------------------------------------------------------------------------
    /**
     *  \brief L1 norm of the term of order m: returns \f$ L_1 \left( [this]_m \right) \f$
     **/
    double l1norm(int const& m);

    /**
     *  \brief Infinity norm of the term of order m: returns \f$ L_\infty \left( [this]_m \right) \f$
     **/
    double linfnorm(int const& m);

    /**
     *  \brief  Number of small divisors under a certain value sdmax in the term of order m
     */
    int nsd(int const& m, int odmax, double sdmax);


    //------------------------------------------------------------------------------------
    //Operations with TFS coefficients - pure operations
    //------------------------------------------------------------------------------------
    //----------------
    // Frequency to Time domain
    //---------------
    /**
     *  \brief Puts all Ofs coefficients in a into Tfs coefficients and puts it into this.
     **/
    Ofts<T>& tfs_from_ofs(Ofts<T> const& a);

    /**
     *  \brief Puts all Ofs coefficients in a into Tfs coefficients and puts it into this. Inline
     **/
    Ofts<T>& tfs_from_ofs_inline();

    /**
     *  \brief Puts all Tfs coefficients in a into ofs coefficients and puts it into this.
     *         Makes use of FFT routines from GSL on each Fourier coefficient.
     **/
    Ofts<T>& tfs_to_ofs(Ofts<T> const& a);

    /**
     *  \brief Puts all Tfs coefficients in a into ofs coefficients and puts it into this. Inline
     *         Makes use of FFT routines from GSL on each Fourier coefficient.
     **/
    Ofts<T>& tfs_to_ofs_inline();

    /**
     *  \brief Puts all Ofs coefficients in a into Tfs coefficients and puts it into this. Inline.
     **/
    Ofts<T>& tfs_from_ofs_inline(int nrc);

    /**
     *  \brief Puts all Tfs coefficients in a into ofs coefficients and puts it into this. Inline
     *         Makes use of FFT routines from GSL on each Fourier coefficient.
     **/
    Ofts<T>& tfs_to_ofs_inline(int nrc);

    //----------------
    // pows
    //---------------
    /**
     *   \brief Power function: p = a^alpha, with Tfs coefficients.
     *          WARNING: for alpha < 0, a good precision is achieved ONLY if
     *          the (Fourier) order 0 of the (Taylor) order 0 is >> wrt to the rest of the (Fourier) coefficients of the (Taylor) order 0.
     **/
    template<typename U> Ofts<T>& tfts_pows(Ofts<T> const& a,  U const& alpha);

    /**
     *   \brief Power function: p = a^alpha, with Tfs coefficients, at order n
     *          WARNING: for alpha < 0, a good precision is achieved ONLY if
     *          the (Fourier) order 0 of the (Taylor) order 0 is >> wrt to the rest of the (Fourier) coefficients of the (Taylor) order 0.
     **/
    template<typename U> Ofts<T>& tfts_pows(Ofts<T> const& a,  U const& alpha, int const& n);

    /**
     *   \brief Compute the order zero of the power function: p = a^alpha, with Tfs coefficients at order n
     *          WARNING: for alpha < 0, a good precision is achieved ONLY if
     *          the (Fourier) order 0 of the (Taylor) order 0 is >> wrt to the rest of the (Fourier) coefficients of the (Taylor) order 0.
     **/
    template<typename U> Ofts<T>& tfts_spows_zero(Ofts<T> const& a,  U const& alpha, int const& n);

    //----------------
    // smult
    //----------------
    /**
     *  \brief  An operation. Adds the product: \c this \f$  += m a \f$ with m a coefficient, at order n
     *  \param  a: a reference to an Oftsh object
     *  \param  m: a reference to a coefficient
     *  \param  n: a reference to the order to update
     */
    Ofts<T>& tfts_smult_t(Ofts<T> const& a, T const& m, int const& n);
    /**
     *  \brief  An operation. Adds the product: \c this \f$  += c a \f$ with c a subcoefficient at order n.
     *  \param  a: a reference to an Oftsh object
     *  \param  c: a reference to a subcoefficient
     */
    template<typename U> Ofts< Ofs<U> >& tfts_smult_u(Ofts< Ofs<U> > const& a, U const& c, int const& n);

    //----------------
    // sprod
    //----------------
    /**
     *  \brief  An operation. Adds the product: \c this \f$  += a*b \f$ at order n.
     *
     *  Handle the case for which n >= max(a.order, b.order)
     */
    Ofts<T>& tfts_sprod(Ofts<T> const& a, Ofts<T> const& b, int const& n);

    /**
     *  \brief  An operation. Adds the order zero of the product: \c this \f$  += a*b \f$ at order n.
     *
     *  Handle the case for which n >= max(a.order, b.order)
     */
    Ofts<T>& tfts_sprod_zero(Ofts<T> const& a, Ofts<T> const& b, int const& n);

    /**
     *  \brief  An operation. Adds the product: \c this \f$  += m*a*b \f$ with m a coefficient at order n
     *  \param  a: a reference to an Ofts object
     *  \param  b: a reference to an Ofts object
     *  \param  n: a reference to a coefficient
     */
    template<typename U> Ofts<T>& tfts_smprod_u(Ofts<T> const& a, Ofts<T> const& b, U const& c, int const& n);

    /**
     *  \brief  An operation. Adds the order zero of the product: \c this \f$  += m*a*b \f$ with m a coefficient at order n
     *  \param  a: a reference to an Ofts object
     *  \param  b: a reference to an Ofts object
     *  \param  n: a reference to a coefficient
     */
    template<typename U> Ofts<T>& tfts_smprod_u_zero(Ofts<T> const& a, Ofts<T> const& b, U const& c, int const& n);

    //----------------
    // sfsum
    //----------------
    /**
     *  \brief  An operation. Adds the product: \c this \f$  += ma*a + mb*b \f$ with ma and mb coefficients at order n.
     *  \param  a: a reference to an Ofts object
     *  \param  ma: a reference to a coefficient
     *  \param  b: a reference to an Ofts object
     *  \param  mn: a reference to a coefficient
     */
    Ofts<T>& tfts_sfsum_t(Ofts<T> const& a, T const& ma, Ofts<T> const& b, T const& mb, int n);

    /**
     *  \brief  An operation. Adds the product: \c this \f$  += ma*a + mb*b + mc*c \f$ with ma and mb coefficients at order n.
     *  \param  a: a reference to an Ofts object
     *  \param  ma: a reference to a coefficient
     *  \param  b: a reference to an Ofts object
     *  \param  mb: a reference to a coefficient
     *  \param  c: a reference to an Ofts object
     *  \param  mc: a reference to a coefficient
     *  \param  n: a reference to the order to update
     */
    Ofts<T>& tfts_sfsum_tt(Ofts<T> const& a, T const& ma, Ofts<T> const& b, T const& mb, Ofts<T> const& c, T const& mc, int n);

    /**
     *  \brief  An operation. Sets the product: \c this \f$  = ca*a + cb*b \f$ with ca and cb subcoefficients at order m
     *  \param  a: a reference to an Ofts object
     *  \param  ca: a reference to a subcoefficient
     *  \param  b: a reference to an Ofts object
     *  \param  cb: a reference to a subcoefficient
     *  \param  m: a reference to the order to update
     */
    template<typename U> Ofts<T>& tfts_fsum_u(Ofts<T> const& a, U const& ca, Ofts<T> const& b, U const& cb, int const& m);

    //----------------
    // der
    //----------------
    /**
     *   \brief Partial derivative at order m: works for order = a.order. TFS format.
     *          WARNING: order m is derived and set in order m-1 of this!
     *          If m==0, nothing is done.
     **/
    Ofts<T>& tfts_der(Ofts< T > const& a, int ni, int m);

    /**
     *   \brief Same as der but the result is added to the current Ofts instance.
     **/
    Ofts<T>& tfts_sder(Ofts< T > const& a, int ni, int m);

    //------------------------------------------------------------------------------------
    //Stream
    //------------------------------------------------------------------------------------
    /**
     *   \brief Stream operator << for Ofts objects.
     **/
    friend std::ostream& operator << <>(std::ostream& stream, Ofts<T> const& ofts);

    /**
     *   \brief Stream operator for Ofts objects. Print only the order zero of each coefficients.
     **/
    void fprint_0(ofstream& stream);

};


//----------------------------------------------------------------------------------------
// Ofts 2 Ofs
//----------------------------------------------------------------------------------------
/**
 *  \brief Transform a one-variable Ofts< Ofs<U> > fts_z into an Ofs<U> object fs object such that fts_z(epsilon) = fs.
 **/
template <typename U>  void fts2fs(Ofs<U> *fs, Ofts< Ofs<U> > const& fts_z, double epsilon);


/**
 *  \brief Transform a one-variable Ofts< Ofs<U> > fts_z into a complex number fts_z(epsilon, t).
 **/
template <typename U>  cdouble fts2scalar(Ofts< Ofs<U> > const& fts_z, double epsilon, double t);


//----------------------------------------------------------------------------------------
//          Reading & writing (I/O)
//----------------------------------------------------------------------------------------
//----------------------------------------------
// Text format, write
//----------------------------------------------
/**
 * \brief Writes a given vector W of type \c Ofts<Ofsc >  in a txt files of the form "filename+i.txt", with i = 0, length(W)-1.
 **/
inline void  write_vofts_txt(vector<Ofts<Ofsc > > &W, string filename);

//----------------------------------------------
// Text format, read
//----------------------------------------------
/**
 * \brief Reads a given \c Ofsc  object within a \c Ofts<Ofsc >, in txt format.
 **/
inline void read_ofs_txt(Ofsc &xFFT, ifstream &readStream, int fftN);
/**
 * \brief Reads a given \c Ofts<Ofsc >  object, in txt format.
 **/
inline int read_ofts_txt(Ofts<Ofsc > &x, string filename, int fftN);
/**
 * \brief Reads a given vector W of type \c Ofts<Ofsc >  in a txt files of the form "filename+i.txt", with i = 0, length(W)-1.
 **/
inline void read_vofts_txt(vector<Ofts<Ofsc > >  &W, string filename, int fftN);

//----------------------------------------------
// Binary format, write
//----------------------------------------------
/**
 * \brief Writes a given \c Ofsc  object within a \c Ofts<Ofsc >, in bin format.
 **/
inline void  write_ofs_bin(Ofsc  &xFFT, fstream &myfile);
/**
 * \brief Writes a given \c Ofts<Ofsc >  object, in bin format.
 **/
inline void  write_ofts_bin(Ofts<Ofsc > &W, string filename);
/**
 * \brief Writes a given vector W of type \c Ofts<Ofsc >  in a binary files of the form "filename+i.bin", with i = 0, length(W)-1.
 **/
inline void  write_vofts_bin(vector<Ofts<Ofsc > > &W, string filename);


//----------------------------------------------
// Binary format, read
//----------------------------------------------
/**
 * \brief Reads a given \c Ofsc  object within a \c Ofts<Ofsc >, in bin format.
 **/
inline void read_ofs_bin(Ofsc  &xFFT, fstream &myfile);

/**
 * \brief Reads a given \c Ofts<Ofsc >  object, in bin format.
 **/
inline int read_ofts_bin(Ofts<Ofsc > &W, string filename);

/**
 * \brief Reads a given vector W of type \c Ofts<Ofsc >  in binary files of the form "filename+i.bin", with i = 0, length(W)-1.
 **/
inline void read_vofts_bin(vector<Ofts<Ofsc > >  &W, string filename, int fftN);

//----------------------------------------------
// Text 2 binary
//----------------------------------------------
/**
 * \brief Reads a given vector W of type \c Ofts<Ofsc >  in a txt files of the form "filename+i.txt", with i = 0, length(W)-1.
 *        Writes it again in binary form, in binary files of the form "filename+i.bin", with i = 0, length(W)-1.
 **/
inline void txt2bin_vofts(vector<Ofts<Ofsc > > &W, string filename, int fftN);


//Include the implementation .tpp
#include "ofts.tpp"



#endif // OFTS_H_INCLUDED
