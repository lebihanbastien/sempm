#ifndef OFTSH_H_INCLUDED
#define OFTSH_H_INCLUDED

/**
 * \file oftsh.h
 * \brief Homogeneous Fourier-Taylor series template class
 * \author BLB, with insights from Alex Haro.
 */

//std
#include <iostream>
#include <iomanip>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//custom
#include "manip.h"
#include "ofs.h"
#include "constants.h"

using namespace std;

//----------------------------------------------------------------------------------------
// Oftsh class
//----------------------------------------------------------------------------------------
template <typename T>
class Oftsh;

template<typename T>
std::ostream& operator << (std::ostream& stream, Oftsh<T> const& oftsh);

/** \class Oftsh
 *  \brief Homogeneous Fourier-Taylor series class.
 *
 *  The Otsh class implements homogeneous Fourier-Taylor polynomials and an algebra
 *  on such objects.
 *  The Otsh class implementation is based on the "tree" structure described in
 *  A. Haro, Automatic differentiation tools in computational dynamical systems, 2008.
 *  (http://webmail.maia.ub.es/dsg/2008/0806haro.pdf).
 *  The coefficients are stored in reverse lexicographical order.
 *
 *  Oftsh is a templated class. In practice, it is used in the forms
 *          - Oftsh<Ofs<double>> i.e. Fourier-Taylor polynomials with Fourier series of
 *            real coefficients (Ofsd) as coefficients
 *          - Oftsh<Ofs<double complex>> i.e. Fourier-Taylor polynomials with Fourier
 *            series of complex coefficients (Ofsc) as coefficients
 *
 *  Since the Ofs class implements Fourier series in two formats (frequency domain - ofs -
 *  and time domain - tfs), the Oftsh also contain operation routines that works with
 *  these two formats. The routines that must be used with ofs format have the
 *  prefix "oftsh_", while the routines that must be used with tfs format have
 *  the prefix "tftsh_".
 *
 *  Since Oftsh is a templated class, a .tpp is used instead of a .cpp, and this source
 *  file is called directly in the header (the current file), via ' #include "oftsh.tpp" '
 *
 *  Finally, note that there is potential memory leak in the destructor ~Oftsh<T>.
 *  I haven't succeeded in implementing such a destructor given the tree (recursive)
 *  structure of the Oftsh objects. Hence, one must not use too many temporary Oftsh
 *  objects.
 *
 */
template <typename T>
class Oftsh
{

private:

    // number of variables
    int m_oftsh_nvar;
    // degree of the homogeneous polynomial
    int m_oftsh_order;
    // ptr to the first coefficient
    T *m_oftsh_coef;
    // ptr to the first child
    Oftsh<T> *m_oftsh_term;

public:

    //------------------------------------------------------------------------------------
    //Create
    //------------------------------------------------------------------------------------
    /**
     *  \brief Default constructor of the class Oftsh<T>.
     */
    Oftsh<T>();
    /**
     *  \brief Constructor with given order and number of variables.
     *  \param t_oftsh_nvar: number of variables of the series
     *  \param t_oftsh_order: order of the series
     */
    Oftsh<T>(const int t_oftsh_nvar, const int t_oftsh_order);
    /**
     *  \brief Constructor from a given Oftsh object (without any link).
     *  \param b:  a reference to the Oftsh object to copy in the new object
     */
    Oftsh<T>(Oftsh<T> const& b);

    //------------------------------------------------------------------------------------
    //Delete
    //------------------------------------------------------------------------------------
    /**
     *  \brief Default destructor of the class Oftsh<T>. WARNING: potential memory leak here.
     */
    ~Oftsh<T>();


    //------------------------------------------------------------------------------------
    //Copy
    //------------------------------------------------------------------------------------
    /**
     *  \brief  Linked copy from a given Oftsh object (exact same object is obtained).
     *  \param  b: a reference to the Oftsh object to copy
     *  \return a reference to the current object
     */
    Oftsh<T>& lcopy(Oftsh<T> const& b);         //linked copy (for recursivity)
    /**
     *  \brief  Copy from a given Oftsh object (only the coefficients).
     *  \param  b: a reference to the Oftsh object to copy
     *  \return a reference to the current object
     */
    Oftsh<T>& ccopy(Oftsh<T> const& b); //coefficient copy (restricted to same order, same number of variables)
    /**
     *  \brief  An operator. Constructor from a given Oftsh object (only the coefficients).
     *  \param  b: a reference to the Oftsh object to copy
     *  \return a reference to the current object
     */
    Oftsh<T>& operator = (Oftsh<T> const& b);

    //------------------------------------------------------------------------------------
    //Linking
    //------------------------------------------------------------------------------------
    /**
     *  \brief  Performs linking between the Oftsh<T> and a coefficient array.
     *  \param  coef0: the array of coefficients to link
     *  \return a reference to the current object
     */
    void link_coefs(T *coef0);

    //------------------------------------------------------------------------------------
    //Setters
    //------------------------------------------------------------------------------------
    /**
     *  \brief Sets a coefficient at a given position in the polynomial.
     *  \param value: the value to set
     *  \param pos: the position to modify
     */
    void set_coef(T  const&  value, int pos);

    /**
     *  \brief Adds a coefficient at a given position in the polynomial.
     *  \param value: the value to add
     *  \param pos: the position to modify
     */
    void add_coef(T const& value, int pos);

    /**
     *  \brief Sets all subcoefficient of the coefficient \c pos to \c value.
     *  \param value: the value to set
     *  \param pos: the position to modify
     */
    template<typename U> void set_sub_coef(U value, int pos);

    /**
     *  \brief Sets subcoefficient \c i of the coefficient \c pos to \c value.
     *  \param value: the value to set
     *  \param pos: the position to modify
     */
    template<typename U> void set_sub_coef(U value, int pos, int i);

    /**
     *  \brief Sets random coefficients to all positions in the polynomial.
     */
    void set_random_coefs();


    //------------------------------------------------------------------------------------
    //Getters
    //------------------------------------------------------------------------------------
    /**
     *  \brief  Gets the first child in the tree.
     *  \return and Oftsh object
     */
    Oftsh<T> get_term() const;

    /**
     *  \brief  Gets the child \c i in the tree.
     *  \return and Oftsh object
     */
    Oftsh<T> get_term(int i) const;

    /**
     *  \brief  Gets the coefficient at a given position.
     *  \param  i: the position to get
     *  \return the coefficient of type \c T at the position \c i
     */
    T get_coef(int i) const;

    /**
     *  \brief  Gets the order of the polynomial.
     *  \return the order of the polynomial as an \c int
     */
    int get_order() const;

    /**
     *  \brief  Gets the number of variables of the polynomial.
     *  \return the number of variables of the polynomial as an \c int
     */
    int get_nvar() const;

    /**
     *  \brief  Gets the address of the first coefficient
     *  \return the address of the first coefficient
     */
    T* get_ptr_first_coef() const;

    //------------------------------------------------------------------------------------
    //Zeroing
    //------------------------------------------------------------------------------------
    /**
     *  \brief  Sets all coefficients to zero.
     */
    void zero();

    //------------------
    //Operators
    //------------------
    //bool isEqual(Oftsh<T> const& b) const;
    //Oftsh<T>& operator += (Oftsh<T> const& b);
    //Oftsh<T>& operator -= (Oftsh<T> const& b);
    //Oftsh<T>& operator *= (T const& c);
    //Oftsh<T>& operator /= (T const& c);

    //------------------------------------------------------------------------------------
    //Operations
    //------------------------------------------------------------------------------------
    //------------------
    // Conjugate
    //------------------
    /**
     *  \brief Conjugates the coefficients the Oftsh object (and only them!).
     *         WARNING: To be used with evaluate_conjugate to have the true conjugate.
     *   The conjugate of the series \f$ T_n = \sum \limits_{|r| = n} c_r x^r \f$ is
     *  \f$ \bar{T}_n = \sum \sum \limits_{|r| = n} \bar{c}_{r} x^r\f$
     */
    Oftsh<T>& conjugate();

    //------------------
    // smult
    //------------------
    /**
     *  \brief  An operation. Adds the product: \c this \f$  += c a \f$ with c a coefficient.
     *  \param  a: a reference to an Oftsh object
     *  \param  c: reference to a coefficient
     */
    Oftsh<T>& oftsh_smult_t(Oftsh<T> const& a, T const& c);

    /**
     *  \brief  An operation. Adds the product: \c this \f$  += c1*c2 a \f$ with c1 and
     *          c2 coefficients. Warning: set for Ofs coefficient by default.
     */
    Oftsh<T>& oftsh_smult_tt(Oftsh<T> const& a, T const& c1, T const& c2, T& temp);

    /**
     *  \brief  An operation. Adds the product: \c this \f$  += m a \f$ with m a
     *          subcoefficient. Oftsh< Ofs<U> > case.
     *  \param  a: a reference to an Oftsh object
     *  \param  m: reference to a subcoefficient
     */
    template<typename U>
    Oftsh< Ofs<U> >& oftsh_smult_u(Oftsh< Ofs<U> > const& a, U const& m);

    /**
     *  \brief  An operation. Adds the product: \c this \f$  += m * r_a * a \f$ with m a
     *          subcoefficient and ra a coefficient. Oftsh< Ofs<U> > case.
     *  \param  a: a reference to an Oftsh object
     *  \param  ra: reference to a coefficient
     *  \param  m: reference to a subcoefficient
     */
    template<typename U>
    Oftsh< Ofs<U> >& oftsh_smult_tu(Oftsh< Ofs<U> > const& a, Ofs<U> const& ra, U const& m);

    //------------------

    //------------------
    // mult
    //------------------
    /**
     *  \brief  An operation. Sets the product: \c this \f$  += c a \f$ with c a
     *          coefficient.
     */
    Oftsh<T>&  oftsh_mult_t(Oftsh<T> const& a, T const& c);

    /**
     *  \brief  An operation. Sets the product: \c this \f$  = m a \f$ with m a
     *          subcoefficient. Oftsh< Ofs<U> > case.
     */
    template<typename U>
    Oftsh< Ofs<U> >& oftsh_mult_u(Oftsh< Ofs<U> > const& a, U const& m);

    //------------------
    // sprod
    //------------------
    /**
     *  \brief  An operation. Adds the product: \c this \f$ += a \times b \f$ with
     *          a and b Oftsh objects
     *  \param  a: a reference to an Oftsh object
     *  \param  b: reference to an Oftsh object
     *  \return a reference to the current object
     */
    Oftsh<T>&  oftsh_sprod(Oftsh<T> const& a, Oftsh<T> const& b);

    //------------------
    // smprod
    //------------------
    /**
     *  \brief  An operation. Adds the product: \c this \f$ += c a \times b \f$ with
     *          a and b Oftsh objects.
     *  \param  a: a reference to an Oftsh object
     *  \param  b: reference to an Oftsh object
     *  \param  c: reference to a coefficient
     *  \return a reference to the current object
     */
    Oftsh<T>& oftsh_smprod_t(Oftsh< T > const& a, Oftsh< T> const& b, T const& c, T& temp);


    /**
     *  \brief  An operation. Adds the product: \c this \f$ += m a \times b \f$ with
     *          a and b Oftsh objects. Oftsh< Ofs<U> > case.
     *  \param  a: a reference to an Oftsh object
     *  \param  b: reference to an Oftsh object
     *  \param  m: reference to a subcoefficient
     *  \return a reference to the current object
     */
    template<typename U>
    Oftsh< Ofs<U> >& oftsh_smprod_u(Oftsh< Ofs<U> > const& a, Oftsh< Ofs<U> > const& b, U const& m);

    //------------------
    // Derivation
    //------------------
    /**
     *  \brief  An operation. Applies the partial derivative with respect to the
     *          variable \c ni: this \f$ = \frac{\partial a}{\partial x_{ni}} \f$
     *
     *  \param  a: a reference to an Oftsh object
     *  \param ni: an \c int representing the derivative variable
     *  \return a reference to the current object
     */
    Oftsh<T>&  derh(Oftsh< T > const& a, int ni);

    /**
     *  \brief  An operation. Adds the partial derivative with respect to the
     *          variable \c ni: this \f$ += \frac{\partial a}{\partial x_{ni}} \f$
     *
     *  \param  a: a reference to an Oftsh object
     *  \param ni: an \c int representing the derivative variable
     *  \return a reference to the current object
     */
    Oftsh<T>& sderh(Oftsh< T > const& a, int ni);

    /**
     *  \brief  An operation. Set the time derivative of object \c a with
     *          pulsation \f$ \omega = n \f$, so that \c this \f$ = \dot{a} \f$.
     *
     *  \param  a: a reference to an Ofs object
     *  \param  n: reference to the pulsation
     *  \return a reference to the current object
     */
    Oftsh<T>& dot(Oftsh<T> const& a, double const&  n);


    //------------------------------------------------------------------------------------
    // TFS operations:
    // They allow to do the same operations as their OFS ("oftsh_") counterpart.
    // E.g: tftsh_sprod does the same thing as oftsh_sprod but for Oftsh<Ofs> objects with
    // the Ofs in time domain. Note that there is probably no need to template these
    // routines but the template has been kept in order to have the same form
    // as the OFS equivalents.
    //------------------------------------------------------------------------------------
    //------------------
    // sprod
    //------------------
    Oftsh<T>& tftsh_sprod(Oftsh<T> const& a, Oftsh<T> const& b);

    //------------------
    // smprod
    //------------------
    Oftsh<T>& tftsh_smprod_t(Oftsh< T > const& a, Oftsh< T> const& b, T const& c);
    template<typename U> Oftsh<T>& tftsh_smprod_tu(Oftsh< T > const& a, Oftsh< T> const& b, T const& c, U const& m);
    template<typename U> Oftsh< Ofs<U> >& tftsh_smprod_u(Oftsh< Ofs<U> > const& a, Oftsh< Ofs<U> > const& b, U const& m);

    //------------------
    // smult
    //------------------
    template<typename U> Oftsh< Ofs<U> >& tftsh_smult_u(Oftsh< Ofs<U> > const& a, U const& c);
    template<typename U> Oftsh< Ofs<U> >& tftsh_mult_u(Oftsh< Ofs<U> > const& a, U const& c);
    Oftsh<T>& tftsh_smult_t(Oftsh<T> const& a, T const& c);
    Oftsh<T>& tftsh_smult_tt(Oftsh<T> const& a, T const& c1, T const& c2);
    template<typename U> Oftsh<T>& tftsh_smult_ttu(Oftsh<T> const& a, T const& c1, T const& c2, U const& m);
    template<typename U> Oftsh< Ofs<U> >& tftsh_smult_tu(Oftsh< Ofs<U> > const& a, Ofs<U> const& ra, U const& c);

    //------------------
    // derh
    //------------------
    Oftsh<T>& tfts_derh(Oftsh<T> const& a, int ni);
    Oftsh<T>& tfts_sderh(Oftsh<T> const& a, int ni);

    //------------------
    // Integral
    //------------------
    /**
    *  \brief  An operation. Adds the partial primitive with respect to the variable \c ni
    *
    *  Notes:
    *  1. If a is of order n, this is of order n+1.
    *  2. Need an extension to autonomous case: WORKS ONLY FOR Ofs<double/cdouble> coefficients.
    */
    Oftsh<T>& sprimh(Oftsh< T > const& a, int ni);

    //------------------
    // Evaluate
    //------------------
    /**
     *  \brief  Evaluates the Ofs object at coordinates X and set it in \c z: \c z \f$ = T_n(X) \f$.
     *  \param  X: an array of coefficients
     *  \param  z: a reference to the coefficient to update.
     */
    template<typename U> void evaluate(U X[], T& z);
    void evaluate(double X[], T& z);//Double case

    /**
     *  \brief  Evaluates the Ofs object at coordinates X and set it in \c z: \c z \f$ += T_n(X) \f$.
     */
    template<typename U> void sevaluate(U X[], T& z);
    template<typename U> void sevaluate(U X[], T& z, int const& ofs_order);
    //template<typename U> void sevaluate(U X[], T* z);
    void sevaluate(double X[], T& z); //Double case

    /**
     *  \brief  Evaluates the Ofs object at coordinates X and time t and adds it in \c z: \c z \f$ += T_n(X) \f$. cdouble version
     */
    template<typename U> void sevaluatedc(U X[], U& z, double const& t, int const& ofs_order);
    template<typename U> void fevaluate(U X[], U& z, int kv[], double cR[], double sR[], int const& ofs_order);

    /**
    *  \brief  Evaluates the Ofs object at the conjugate of coordinates X and set it in \c z: \c z \f$ = T_n(\bar{X}) \f$.
    */
    template<typename U> void evaluate_conjugate(U X[], T& z);
    void evaluate_conjugate(double X[], T& z); //Double case

    /**
    *  \brief  Evaluates the Ofs object at the conjugate of coordinates X and set it in \c z: \c z \f$ += T_n(\bar{X}) \f$.
    */
    template<typename U> void sevaluate_conjugate(U X[], T& z);
    void sevaluate_conjugate(double X[], T& z); //Double case

    //------------------
    // Norms
    //------------------
    /**
     *  \brief  Evaluates the \f$ L_1 \f$ norm of the current Oftsh object.
     *  \return a \c double containing \f$ \sum \limits_{|r| = n} ||c_r||_{L_1} \f$ for the pol. \f$ T_n = \sum \limits_{|r| = n} c_r x^r\f$.
     */
    double l1norm();

    /**
     *  \brief  Evaluates the \f$ L_{\infty} \f$ norm of the current Oftsh object.
     */
    double linfnorm();

    //------------------
    // Small divisors
    //------------------
    /**
    *  \brief  Number of small divisors under a certain value
    */
    int nsd(int odmax, double sdmax);

    //------------------
    // I/O
    //------------------
    /**
     *  \brief  A stream operator
     */
    friend std::ostream& operator << <>(std::ostream& stream, Oftsh<T> const& oftsh);
};

//----------------------------------------------------------------------------------------
// Functions
//----------------------------------------------------------------------------------------
/**
 * \fn template<typename T> Oftsh<T> operator + (Oftsh<T> const& a, Oftsh<T> const& b)
 * \brief An operator. Makes the sum a+b
 * \param a: a reference to an Oftsh object
 * \param b: a reference to an Oftsh object
 * \return a+b at order max(a.order, b.order).
 */
template<typename T> Oftsh<T> operator + (Oftsh<T> const& a, Oftsh<T> const& b);

/**
 * \fn template<typename T> Oftsh<T> operator - (Oftsh<T> const& a, Oftsh<T> const& b)
 * \brief An operator. Makes the sub a-b
 * \param a: a reference to an Oftsh object
 * \param b: a reference to an Oftsh object
 * \return a-b at order max(a.order, b.order).
 */
template<typename T> Oftsh<T> operator - (Oftsh<T> const& a, Oftsh<T> const& b);

/**
 *  \brief  Expected error on the product: \c this \f$ += c \times a \f$ at times t.
 *          Works only when a.order = b.order which is the default case.
 */
template <typename U> cdouble smult_expct_error(Oftsh<Ofs<U> > const& a, Ofs<U> const& c, U X[], double const& t);

/**
 *  \brief  Expected error on the product: \c this \f$ += c \times a \f$ at times t.
 *          Works only when a.order = b.order which is the default case. Real case.
 */
cdouble smult_expct_error(Oftsh<Ofsd > const& a, Ofsd const& c, double X[], double const& t);

//Include the implementation .tpp
#include "oftsh.tpp"

#endif // OFTSH_H_INCLUDED
