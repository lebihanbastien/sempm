/**
 * \file ofs_test.cpp
 * \brief Test file for the Ofs (Fourier series) class.
 * \author BLB
 */

#include "ofs_test.h"
#define VAR 2.0

using namespace std;

/**
 * \fn void ofs_test()
 * \brief Main routine for Ofs class test.
 */
void ofs_test()
{
    char ch;
    //------------------------------------------------------------------------------------
    // Splash
    //------------------------------------------------------------------------------------
    cout << "---------------------------------------------------" << endl;
    cout << "                                                   " << endl;
    cout << "    Test routine for the Ofs class                 " << endl;
    cout << "           (Fourier series)                        " << endl;
    cout << "                                                   " << endl;
    cout << "---------------------------------------------------" << endl;
    cout << "The following public functions are NOT tested:     " << endl;
    cout << "---------------------------------------------------" << endl;
    cout << "- The constructors" << endl;
    cout << "- The destructor " << endl;
    cout << "- Ofs<T>& ccopy(Ofs<T> const& ofs_); " << endl;
    cout << "- Ofs<T>& lcopy(Ofs<T> const& ofs_); " << endl;
    cout << "- void set_coef(T const value, int const pos); " << endl;
    cout << "- void add_coef(T const value, int const pos); " << endl;
    cout << "- void set_all_coefs(T const value); " << endl;
    cout << "- int get_order() const; " << endl;
    cout << "- Ofs* get_ptr() const; " << endl;
    cout << "- T get_coef(int pos) const; " << endl;
    cout << "- void get_coef_max_norm(double maxAbs[]) const; " << endl;
    cout << "- double get_coef_max_norm() const; " << endl;
    cout << "- void zero(); " << endl;
    cout << "- Ofs<T>& operator  = (Ofs<T> const& ofs_); " << endl;
    cout << "- Ofs<T>& operator  = (T const& coef0); " << endl;
    cout << "- cdouble evaluate(double const& t); " << endl;
    cout << "- double l1norm(); " << endl;
    cout << "- bool operator == (Ofs<T> const& a, Ofs<T> const& b); " << endl;
    cout << "- bool isEqual(Ofs<T> const& b) const;" << endl;


    cout << "---------------------------------------------------" << endl;
    cout << "The following public functions are tested:                 " << endl;
    cout << "---------------------------------------------------" << endl;
    cout << "- void    conjugate();" << endl;
    cout << "- Ofs<T>& operator += (Ofs<T> const& ofs_);" << endl;
    cout << "- Ofs<T>& operator -= (Ofs<T> const& ofs_);" << endl;
    cout << "- Ofs<T>& operator *= (T const& c);" << endl;
    cout << "- Ofs<T>& operator /= (T const& c);" << endl;
    cout << "- Ofs<T>& sprod(Ofs<T> const& a, Ofs<T> const& b);" << endl;
    cout << "- Ofs<T>& prod(Ofs<T> const& a, Ofs<T> const& b);" << endl;
    cout << "- Ofs<T>& smprod(Ofs<T> const& a, Ofs<T> const& b, T const& m);" << endl;
    cout << "- Ofs<T>& mprod(Ofs<T> const& a, Ofs<T> const& b, T const& m);" << endl;
    cout << "- Ofs<T>& smult(Ofs<T> const& a, T const& c);" << endl;
    cout << "- Ofs<T>& mult(Ofs<T> const& a, T const& c);" << endl;
    cout << "- Ofs<T>& fsum(Ofs<T> const& a, T const& ma, Ofs<T> const& b, T const& mb);" << endl;
    cout << "- Ofs<T>& dot(Ofs<T> const& a, double const& n););" << endl;
    cout << "- Ofs<T>& dot(double const& n););" << endl;
    cout << "- Ofs<T>& epspow(Ofs<T> const& a, T const& alpha);" << endl;
    cout << "- Ofs<T>& ofs_pow(Ofs<T> const& a, T const& alpha);" << endl;
    cout << "- Ofs<T>  operator + (Ofs<T> const& a, Ofs<T> const& b)" << endl;
    cout << "- Ofs<T>  operator - (Ofs<T> const& a, Ofs<T> const& b)" << endl;
    cout << "- Ofs<T>  operator - (Ofs<T> const& b)" << endl;
    cout << "- Ofs<T>  operator * (Ofs<T> const& a, T const& c)" << endl;
    cout << "- Ofs<T>  operator * (T const& c, Ofs<T> const& a)" << endl;
    cout << "- Ofs<T>  operator / (Ofs<T> const& a, T const& c)" << endl;
    cout << "- Ofs<T>  operator * (Ofs<T> const& a, Ofs<T> const& b)" << endl;

    cout <<  setw(5) << setprecision(1) << setiosflags(ios::scientific);
    cout << "---------------------------------------------------" << endl;
    cout << " 1. Tests are made with series of order " << OFS_ORDER << "." << endl;
    cout << " 2. The series are evaluated at arbitrary time omega*t = " << VAR << "." << endl;
    cout << " 3. Random coefficients are set as input into every Ofs. " << endl;
    cout << "However, an arbitrary decreasing of the coefficients is set: " << endl;
    cout << " around " << (double)rand()/(10.0*(pow(0,7.0)+1)*RAND_MAX) << " @ order 0." << endl;
    cout << " around " << (double)rand()/(10.0*(pow(OFS_ORDER,7.0)+1)*RAND_MAX) << " @ order "<< OFS_ORDER << "." << endl;
    cout << " 4. Whenever it is possible, expected error is displayed. " << endl;
    cout << "---------------------------------------------------" << endl;
    cout <<  setw(5) << setprecision(15) << std::showpos  <<  setiosflags(ios::scientific);
    cout << "Press Enter to proceed with the tests." << endl;
    scanf("%c",&ch);


    //------------------------------------------------------------------------------------
    // Tests.
    //------------------------------------------------------------------------------------
    ofs_test_pows();;
    cout << "Press Enter to proceed with the tests." << endl;
    scanf("%c",&ch);

    ofs_test_conjugate();
    cout << "Press Enter to proceed with the tests." << endl;
    scanf("%c",&ch);

    ofs_test_operators();
    cout << "Press Enter to proceed with the tests." << endl;
    scanf("%c",&ch);

    ofs_test_sprod();
    cout << "Press Enter to proceed with the tests." << endl;
    scanf("%c",&ch);

    ofs_test_prod();
    cout << "Press Enter to proceed with the tests." << endl;
    scanf("%c",&ch);

    ofs_test_smprod();
    cout << "Press Enter to proceed with the tests." << endl;
    scanf("%c",&ch);

    ofs_test_mprod();
    cout << "Press Enter to proceed with the tests." << endl;
    scanf("%c",&ch);

    ofs_test_smult();
    cout << "Press Enter to proceed with the tests." << endl;
    scanf("%c",&ch);

    ofs_test_mult();
    cout << "Press Enter to proceed with the tests." << endl;
    scanf("%c",&ch);

    ofs_test_fsum();
    cout << "Press Enter to proceed with the tests." << endl;
    scanf("%c",&ch);

    ofs_test_dot();
    cout << "Press Enter to proceed with the tests." << endl;
    scanf("%c",&ch);

    ofs_test_epspow();
    cout << "Press Enter to proceed with the tests." << endl;
    scanf("%c",&ch);

    ofs_test_smprod_t();
    cout << "Press Enter to proceed with the tests." << endl;
    scanf("%c",&ch);

    cout << "---------------------------------------------------" << endl;
    cout << "End of test session.                               " << endl;
    cout << "---------------------------------------------------" << endl;
}

/**
 * \fn void ofs_test_conjugate()
 * \brief Test of the routine: void conjugate().
 */
void ofs_test_conjugate()
{

    cout << "---------------------------------------------------" << endl;
    cout << "    Test of the routine: void conjugate()          " << endl;
    cout << "---------------------------------------------------" << endl;


    //Initialization
    Ofsd xd(OFS_ORDER);
    Ofsc xdc(OFS_ORDER);
    cdouble res1, res2;
    cdouble resd1, resd2;

    //Set random coefs in xd and xdc
    xd.set_random_coefs();
    xdc.set_random_coefs();

    //Set results
    res1  = conj(xd.evaluate(VAR));
    resd1 = conj(xdc.evaluate(VAR));

    //Take conjugate
    xd.conjugate();
    xdc.conjugate();

    //Set results
    res2  = xd.evaluate(VAR);
    resd2 = xdc.evaluate(VAR);

    //Comparison. double case
    cout << "1. double case.         " << endl;
    cout << "Expected result: " << creal(res1) << cimag(res1) << "i" << endl;
    cout << "Obtained result: " << creal(res2) << cimag(res2) << "i" << endl;
    res2 -= res1;
    cout << "Delta: " << cabs(res2) << endl;

    //Comparison. double case
    cout << "2. cdouble case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;
}

/**
 * \fn void ofs_test_operators()
 * \brief Test of the operator routines: +=, -=, *= and /=
 */
void ofs_test_operators()
{

    cout << "---------------------------------------------------" << endl;
    cout << "  Test of the operator routines: +=, -=, *= and /= " << endl;
    cout << "---------------------------------------------------" << endl;


    //Initialization
    Ofsd xd(OFS_ORDER);
    Ofsd xd2(OFS_ORDER);
    Ofsc xdc(OFS_ORDER);
    Ofsc xdc2(OFS_ORDER);
    cdouble res1, res2;
    cdouble resd1, resd2;
    //Arbitrary coefficients
    double c = 2.0;
    cdouble cd = 2.0 + 1.0*I;

    //Set random coefs in xd and xdc
    xd.set_random_coefs();
    xd2.set_random_coefs();
    xdc.set_random_coefs();
    xdc2.set_random_coefs();

    //------------------------
    // Operator +=
    //------------------------
    cout << "---------------------" << endl;
    cout << "Operator +=          " << endl;
    cout << "---------------------" << endl;
    //Set results
    res1  = xd.evaluate(VAR) + xd2.evaluate(VAR);
    resd1  = xdc.evaluate(VAR) + xdc2.evaluate(VAR);

    //Take sum
    xd  += xd2;
    xdc += xdc2;

    //Set results
    res2  = xd.evaluate(VAR);
    resd2 = xdc.evaluate(VAR);

    //Comparison. double case
    cout << "1. double case.         " << endl;
    cout << "Expected result: " << creal(res1) << cimag(res1) << "i" << endl;
    cout << "Obtained result: " << creal(res2) << cimag(res2) << "i" << endl;
    res2 -= res1;
    cout << "Delta: " << cabs(res2) << endl;

    //Comparison. double case
    cout << "2. cdouble case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;


    //------------------------
    // Operator -=
    //------------------------
    cout << "---------------------" << endl;
    cout << "Operator -=          " << endl;
    cout << "---------------------" << endl;
    //Set results
    res1  = xd.evaluate(VAR) - xd2.evaluate(VAR);
    resd1  = xdc.evaluate(VAR) - xdc2.evaluate(VAR);

    //Take sum
    xd  -= xd2;
    xdc -= xdc2;

    //Set results
    res2  = xd.evaluate(VAR);
    resd2 = xdc.evaluate(VAR);

    //Comparison. double case
    cout << "1. double case.         " << endl;
    cout << "Expected result: " << creal(res1) << cimag(res1) << "i" << endl;
    cout << "Obtained result: " << creal(res2) << cimag(res2) << "i" << endl;
    res2 -= res1;
    cout << "Delta: " << cabs(res2) << endl;

    //Comparison. double case
    cout << "2. cdouble case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;

    //------------------------
    // Operator *=
    //------------------------
    cout << "---------------------" << endl;
    cout << "Operator *=          " << endl;
    cout << "---------------------" << endl;
    //Set results
    res1  =  c*xd.evaluate(VAR);
    resd1  = cd*xdc.evaluate(VAR);

    //Take sum
    xd  *= c;
    xdc *= cd;

    //Set results
    res2  = xd.evaluate(VAR);
    resd2 = xdc.evaluate(VAR);

    //Comparison. double case
    cout << "1. double case.         " << endl;
    cout << "Expected result: " << creal(res1) << cimag(res1) << "i" << endl;
    cout << "Obtained result: " << creal(res2) << cimag(res2) << "i" << endl;
    res2 -= res1;
    cout << "Delta: " << cabs(res2) << endl;

    //Comparison. double case
    cout << "2. cdouble case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;


    //------------------------
    // Operator /=
    //------------------------
    cout << "---------------------" << endl;
    cout << "Operator /=          " << endl;
    cout << "---------------------" << endl;
    //Set results
    res1  =  xd.evaluate(VAR)/c;
    resd1  = xdc.evaluate(VAR)/cd;

    //Take sum
    xd  /= c;
    xdc /= cd;

    //Set results
    res2  = xd.evaluate(VAR);
    resd2 = xdc.evaluate(VAR);

    //Comparison. double case
    cout << "1. double case.         " << endl;
    cout << "Expected result: " << creal(res1) << cimag(res1) << "i" << endl;
    cout << "Obtained result: " << creal(res2) << cimag(res2) << "i" << endl;
    res2 -= res1;
    cout << "Delta: " << cabs(res2) << endl;

    //Comparison. double case
    cout << "2. cdouble case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;


    cout << "---------------------------------------------------" << endl;
    cout << "  Test of the operator routines: +, -, * and / " << endl;
    cout << "---------------------------------------------------" << endl;

    //------------------------
    // Routine +
    //------------------------
    cout << "---------------------" << endl;
    cout << " Routine +           " << endl;
    cout << "---------------------" << endl;
    //Set results
    res1  =  xd.evaluate(VAR)+xd2.evaluate(VAR);
    resd1  = xdc.evaluate(VAR)+xdc2.evaluate(VAR);


    //Set results
    res2  = (xd + xd2).evaluate(VAR);
    resd2 = (xdc + xdc2).evaluate(VAR);

    //Comparison. double case
    cout << "1. double case.         " << endl;
    cout << "Expected result: " << creal(res1) << cimag(res1) << "i" << endl;
    cout << "Obtained result: " << creal(res2) << cimag(res2) << "i" << endl;
    res2 -= res1;
    cout << "Delta: " << cabs(res2) << endl;

    //Comparison. double case
    cout << "2. cdouble case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;


    //------------------------
    // Routine -
    //------------------------
    cout << "---------------------" << endl;
    cout << " Routine -           " << endl;
    cout << "---------------------" << endl;
    //Set results
    res1  =  xd.evaluate(VAR)-xd2.evaluate(VAR);
    resd1  = xdc.evaluate(VAR)-xdc2.evaluate(VAR);


    //Set results
    res2  = (xd - xd2).evaluate(VAR);
    resd2 = (xdc - xdc2).evaluate(VAR);

    //Comparison. double case
    cout << "1. double case.         " << endl;
    cout << "Expected result: " << creal(res1) << cimag(res1) << "i" << endl;
    cout << "Obtained result: " << creal(res2) << cimag(res2) << "i" << endl;
    res2 -= res1;
    cout << "Delta: " << cabs(res2) << endl;

    //Comparison. double case
    cout << "2. cdouble case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;


    //------------------------
    // Operator -
    //------------------------
    cout << "---------------------" << endl;
    cout << " Operator -          " << endl;
    cout << "---------------------" << endl;
    //Set results
    res1  =  0.0*I-xd2.evaluate(VAR);
    resd1  = 0.0*I-xdc2.evaluate(VAR);


    //Set results
    res2  = (-xd2).evaluate(VAR);
    resd2 = (-xdc2).evaluate(VAR);

    //Comparison. double case
    cout << "1. double case.         " << endl;
    cout << "Expected result: " << creal(res1) << cimag(res1) << "i" << endl;
    cout << "Obtained result: " << creal(res2) << cimag(res2) << "i" << endl;
    res2 -= res1;
    cout << "Delta: " << cabs(res2) << endl;

    //Comparison. double case
    cout << "2. cdouble case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;


    //------------------------
    // Routine c*a
    //------------------------
    cout << "---------------------" << endl;
    cout << " Routine c*a         " << endl;
    cout << "---------------------" << endl;
    //Set results
    res1  =  c*xd2.evaluate(VAR);
    resd1  = cd*xdc2.evaluate(VAR);


    //Set results
    res2  = (c*xd2).evaluate(VAR);
    resd2 = (cd*xdc2).evaluate(VAR);

    //Comparison. double case
    cout << "1. double case.         " << endl;
    cout << "Expected result: " << creal(res1) << cimag(res1) << "i" << endl;
    cout << "Obtained result: " << creal(res2) << cimag(res2) << "i" << endl;
    res2 -= res1;
    cout << "Delta: " << cabs(res2) << endl;

    //Comparison. double case
    cout << "2. cdouble case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;

    //------------------------
    // Routine a*c
    //------------------------
    cout << "---------------------" << endl;
    cout << " Routine a*c         " << endl;
    cout << "---------------------" << endl;
    //Set results
    res1  =  c*xd2.evaluate(VAR);
    resd1  = cd*xdc2.evaluate(VAR);


    //Set results
    res2  = (xd2*c).evaluate(VAR);
    resd2 = (xdc2*cd).evaluate(VAR);

    //Comparison. double case
    cout << "1. double case.         " << endl;
    cout << "Expected result: " << creal(res1) << cimag(res1) << "i" << endl;
    cout << "Obtained result: " << creal(res2) << cimag(res2) << "i" << endl;
    res2 -= res1;
    cout << "Delta: " << cabs(res2) << endl;

    //Comparison. double case
    cout << "2. cdouble case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;

    //------------------------
    // Routine a/c
    //------------------------
    cout << "---------------------" << endl;
    cout << " Routine a/c         " << endl;
    cout << "---------------------" << endl;
    //Set results
    res1  =  xd2.evaluate(VAR)/c;
    resd1  = xdc2.evaluate(VAR)/cd;


    //Set results
    res2  = (xd2/c).evaluate(VAR);
    resd2 = (xdc2/cd).evaluate(VAR);

    //Comparison. double case
    cout << "1. double case.         " << endl;
    cout << "Expected result: " << creal(res1) << cimag(res1) << "i" << endl;
    cout << "Obtained result: " << creal(res2) << cimag(res2) << "i" << endl;
    res2 -= res1;
    cout << "Delta: " << cabs(res2) << endl;

    //Comparison. double case
    cout << "2. cdouble case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;

    //------------------------
    // Routine a*b
    //------------------------
    cout << "---------------------" << endl;
    cout << " Routine a*b         " << endl;
    cout << "---------------------" << endl;
    //Set results
    res1  =  xd.evaluate(VAR)*xd2.evaluate(VAR);
    resd1  = xdc.evaluate(VAR)*xdc2.evaluate(VAR);


    //Set results
    res2  = (xd*xd2).evaluate(VAR);
    resd2 = (xdc*xdc2).evaluate(VAR);

    //Comparison. double case
    cout << "1. double case.         " << endl;
    cout << "Expected result: " << creal(res1) << cimag(res1) << "i" << endl;
    cout << "Obtained result: " << creal(res2) << cimag(res2) << "i" << endl;
    res2 -= res1;
    cout << "Delta: " << cabs(res2) << endl;
    cout << "Expected error : " << cabs(sprod_expct_error(xd, xd2, VAR)) << endl;
    cout << "Corrected delta: " << cabs(res2 + sprod_expct_error(xd, xd2, VAR)) << endl;

    //Comparison. double case
    cout << "2. cdouble case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;
    cout << "Expected error : " << cabs(sprod_expct_error(xdc, xdc2,  VAR)) << endl;
    cout << "Corrected delta: " << cabs(resd2 + sprod_expct_error(xdc, xdc2, VAR)) << endl;
}

/**
 * \fn void ofs_test_sprod()
 * \brief Test of the routine: Ofs<T>&  sprod(Ofs<T> const& a, Ofs<T> const& b)
 */
void ofs_test_sprod()
{

    cout << "---------------------------------------------------------------" << endl;
    cout << " Test of the routine: sprod(Ofs<T> const& a, Ofs<T> const& b)  " << endl;
    cout << "---------------------------------------------------------------" << endl;


    //Initialization
    Ofsd xd1(OFS_ORDER);
    Ofsd xd2(OFS_ORDER);
    Ofsd xd3(OFS_ORDER);
    Ofsc xdc1(OFS_ORDER);
    Ofsc xdc2(OFS_ORDER);
    Ofsc xdc3(OFS_ORDER);
    cdouble res1, res2;
    cdouble resd1, resd2;

    //Set random coefs in xd and xdc
    xd1.set_random_coefs();
    xd2.set_random_coefs();
    xd3.set_random_coefs();

    xdc1.set_random_coefs();
    xdc2.set_random_coefs();
    xdc3.set_random_coefs();

    //Set results
    res1  = xd1.evaluate(VAR) + xd2.evaluate(VAR)*xd3.evaluate(VAR);
    resd1 = xdc1.evaluate(VAR) + xdc2.evaluate(VAR)*xdc3.evaluate(VAR);

    //Take produc
    xd1.ofs_sprod(xd2, xd3);
    xdc1.ofs_sprod(xdc2, xdc3);

    //Set results
    res2  = xd1.evaluate(VAR);
    resd2 = xdc1.evaluate(VAR);

    //Comparison. double case
    cout << "1. double case.         " << endl;
    cout << "Expected result: " << creal(res1) << cimag(res1) << "i" << endl;
    cout << "Obtained result: " << creal(res2) << cimag(res2) << "i" << endl;
    res2 -= res1;
    cout << "Delta: " << cabs(res2) << endl;
    cout << "Expected error : " << cabs(sprod_expct_error(xd2, xd3, VAR)) << endl;
    cout << "Corrected delta: " << cabs(res2 + sprod_expct_error(xd2, xd3, VAR)) << endl;


    //Comparison. double case
    cout << "2. cdouble case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;
    cout << "Expected error : " << cabs(sprod_expct_error(xdc2, xdc3, VAR)) << endl;
    cout << "Corrected delta: " << cabs(resd2 + sprod_expct_error(xdc2, xdc3, VAR)) << endl;
}

/**
 * \fn void ofs_test_prod()
 * \brief Test of the routine: Ofs<T>&  prod(Ofs<T> const& a, Ofs<T> const& b)
 */
void ofs_test_prod()
{

    cout << "---------------------------------------------------------------" << endl;
    cout << " Test of the routine: prod(Ofs<T> const& a, Ofs<T> const& b)  " << endl;
    cout << "---------------------------------------------------------------" << endl;


    //Initialization
    Ofsd xd1(OFS_ORDER);
    Ofsd xd2(OFS_ORDER);
    Ofsd xd3(OFS_ORDER);
    Ofsc xdc1(OFS_ORDER);
    Ofsc xdc2(OFS_ORDER);
    Ofsc xdc3(OFS_ORDER);
    cdouble res1, res2;
    cdouble resd1, resd2;

    //Set random coefs in xd and xdc
    xd1.set_random_coefs();
    xd2.set_random_coefs();
    xd3.set_random_coefs();

    xdc1.set_random_coefs();
    xdc2.set_random_coefs();
    xdc3.set_random_coefs();

    //Set results
    res1  = xd2.evaluate(VAR)*xd3.evaluate(VAR);
    resd1 = xdc2.evaluate(VAR)*xdc3.evaluate(VAR);

    //Take produc
    xd1.ofs_prod(xd2, xd3);
    xdc1.ofs_prod(xdc2, xdc3);

    //Set results
    res2  = xd1.evaluate(VAR);
    resd2 = xdc1.evaluate(VAR);

    //Comparison. double case
    cout << "1. double case.         " << endl;
    cout << "Expected result: " << creal(res1) << cimag(res1) << "i" << endl;
    cout << "Obtained result: " << creal(res2) << cimag(res2) << "i" << endl;
    res2 -= res1;
    cout << "Delta: " << cabs(res2) << endl;
    cout << "Expected error : " << cabs(sprod_expct_error(xd2, xd3, VAR)) << endl;
    cout << "Corrected delta: " << cabs(res2 + sprod_expct_error(xd2, xd3, VAR)) << endl;


    //Comparison. double case
    cout << "2. cdouble case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;
    cout << "Expected error : " << cabs(sprod_expct_error(xdc2, xdc3, VAR)) << endl;
    cout << "Corrected delta: " << cabs(resd2 + sprod_expct_error(xdc2, xdc3, VAR)) << endl;
}

/**
 * \fn void ofs_test_smprod()
 * \brief Test of the routine: Ofs<T>&  smprod(Ofs<T> const& a, Ofs<T> const& b, T const& m)
 */
void ofs_test_smprod()
{

    cout << "---------------------------------------------------------------" << endl;
    cout << " Test of the routine: " << endl;
    cout << " Ofs<T>&  smprod(Ofs<T> const& a, Ofs<T> const& b, T const& m)  " << endl;
    cout << "---------------------------------------------------------------" << endl;


    //Initialization
    Ofsd xd1(OFS_ORDER);
    Ofsd xd2(OFS_ORDER);
    Ofsd xd3(OFS_ORDER);
    Ofsc xdc1(OFS_ORDER);
    Ofsc xdc2(OFS_ORDER);
    Ofsc xdc3(OFS_ORDER);
    cdouble res1, res2;
    cdouble resd1, resd2;
    //Arbitrary coefficients
    double c = 2.0;
    cdouble cd = 2.0 + 1.0*I;

    //Set random coefs in xd and xdc
    xd1.set_random_coefs();
    xd2.set_random_coefs();
    xd3.set_random_coefs();

    xdc1.set_random_coefs();
    xdc2.set_random_coefs();
    xdc3.set_random_coefs();

    //Set results
    res1  = xd1.evaluate(VAR) + c*xd2.evaluate(VAR)*xd3.evaluate(VAR);
    resd1 = xdc1.evaluate(VAR) + cd*xdc2.evaluate(VAR)*xdc3.evaluate(VAR);

    //Take produc
    xd1.ofs_smprod(xd2, xd3, c);
    xdc1.ofs_smprod(xdc2, xdc3, cd);

    //Set results
    res2  = xd1.evaluate(VAR);
    resd2 = xdc1.evaluate(VAR);

    //Comparison. double case
    cout << "1. double case.         " << endl;
    cout << "Expected result: " << creal(res1) << cimag(res1) << "i" << endl;
    cout << "Obtained result: " << creal(res2) << cimag(res2) << "i" << endl;
    res2 -= res1;
    cout << "Delta: " << cabs(res2) << endl;
    cout << "Expected error : " << cabs(smprod_expct_error(xd2, xd3, c, VAR)) << endl;
    cout << "Corrected delta: " << cabs(res2 + smprod_expct_error(xd2, xd3, c, VAR)) << endl;


    //Comparison. double case
    cout << "2. cdouble case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;
    cout << "Expected error : " << cabs(smprod_expct_error(xdc2, xdc3, cd, VAR)) << endl;
    cout << "Corrected delta: " << cabs(resd2 + smprod_expct_error(xdc2, xdc3, cd, VAR)) << endl;
}

/**
 * \fn void ofs_test_smprod_t()
 * \brief Test of the routine: Ofs<T>&  ofs_smprod_t(Ofs<T> const& a, Ofs<T> const& b, Ofs<T> const& c)
 */
void ofs_test_smprod_t()
{

    cout << "---------------------------------------------------------------" << endl;
    cout << " Test of the routine: " << endl;
    cout << " Ofs<T>& ofs_smprod_t(Ofs<T> const& a, Ofs<T> const& b, Ofs<T> const& c)  " << endl;
    cout << "---------------------------------------------------------------" << endl;


    //Initialization
    Ofsd xd1(OFS_ORDER);
    Ofsd xd2(OFS_ORDER);
    Ofsd xd3(OFS_ORDER);
    Ofsd xd4(OFS_ORDER);
    Ofsd temp(OFS_ORDER);
    Ofsc xdc1(OFS_ORDER);
    Ofsc xdc2(OFS_ORDER);
    Ofsc xdc3(OFS_ORDER);
    Ofsc xdc4(OFS_ORDER);
    Ofsc tempc(OFS_ORDER);
    cdouble res1, res2;
    cdouble resd1, resd2;

    //Set random coefs in xd and xdc
    xd1.set_random_coefs();
    xd2.set_random_coefs();
    xd3.set_random_coefs();
    xd4.set_random_coefs();

    xdc1.set_random_coefs();
    xdc2.set_random_coefs();
    xdc3.set_random_coefs();
    xdc4.set_random_coefs();

    //Set results
    res1  = xd1.evaluate(VAR) + xd4.evaluate(VAR)*xd2.evaluate(VAR)*xd3.evaluate(VAR);
    resd1 = xdc1.evaluate(VAR) + xdc4.evaluate(VAR)*xdc2.evaluate(VAR)*xdc3.evaluate(VAR);

    //Take produc
    xd1.ofs_smprod_t(xd2, xd3, xd4, temp);
    xdc1.ofs_smprod_t(xdc2, xdc3, xdc4, tempc);

    //Set results
    res2  = xd1.evaluate(VAR);
    resd2 = xdc1.evaluate(VAR);

    //Comparison. double case
    cout << "1. double case.         " << endl;
    cout << "Expected result: " << creal(res1) << cimag(res1) << "i" << endl;
    cout << "Obtained result: " << creal(res2) << cimag(res2) << "i" << endl;
    res2 -= res1;
    cout << "Delta: " << cabs(res2) << endl;
    cout << "Note: the expected error is trickier to obtain since we do two products in a row in ofs_smprod_t " << endl;

    //Comparison. double case
    cout << "2. cdouble case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;
    cout << "Note: same remark here." << endl;
}

/**
 * \fn void ofs_test_mprod()
 * \brief Test of the routine: Ofs<T>&  mprod(Ofs<T> const& a, Ofs<T> const& b, T const& m)
 */
void ofs_test_mprod()
{
    cout << "---------------------------------------------------------------" << endl;
    cout << " Test of the routine: mprod(Ofs<T> const& a, Ofs<T> const& b, T const& c)  " << endl;
    cout << "---------------------------------------------------------------" << endl;


    //Initialization
    Ofsd xd1(OFS_ORDER);
    Ofsd xd2(OFS_ORDER);
    Ofsd xd3(OFS_ORDER);
    Ofsc xdc1(OFS_ORDER);
    Ofsc xdc2(OFS_ORDER);
    Ofsc xdc3(OFS_ORDER);
    cdouble res1, res2;
    cdouble resd1, resd2;
    //Arbitrary coefficients
    double c = 2.0;
    cdouble cd = 2.0 + 1.0*I;

    //Set random coefs in xd and xdc
    xd1.set_random_coefs();
    xd2.set_random_coefs();
    xd3.set_random_coefs();

    xdc1.set_random_coefs();
    xdc2.set_random_coefs();
    xdc3.set_random_coefs();

    //Set results
    res1  = c*xd2.evaluate(VAR)*xd3.evaluate(VAR);
    resd1 = cd*xdc2.evaluate(VAR)*xdc3.evaluate(VAR);

    //Take produc
    xd1.ofs_mprod(xd2, xd3, c);
    xdc1.ofs_mprod(xdc2, xdc3, cd);

    //Set results
    res2  = xd1.evaluate(VAR);
    resd2 = xdc1.evaluate(VAR);

    //Comparison. double case
    cout << "1. double case.         " << endl;
    cout << "Expected result: " << creal(res1) << cimag(res1) << "i" << endl;
    cout << "Obtained result: " << creal(res2) << cimag(res2) << "i" << endl;
    res2 -= res1;
    cout << "Delta: " << cabs(res2) << endl;
    cout << "Expected error : " << cabs(smprod_expct_error(xd2, xd3, c, VAR)) << endl;
    cout << "Corrected delta: " << cabs(res2 + smprod_expct_error(xd2, xd3, c, VAR)) << endl;


    //Comparison. double case
    cout << "2. cdouble case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;
    cout << "Expected error : " << cabs(smprod_expct_error(xdc2, xdc3, cd, VAR)) << endl;
    cout << "Corrected delta: " << cabs(resd2 + smprod_expct_error(xdc2, xdc3, cd, VAR)) << endl;
}

/**
 * \fn void ofs_smult()
 * \brief Test of the operator routine: Ofs<T>& smult(Ofs<T> const& a, T const& c)
 */
void ofs_test_smult()
{
    cout << "---------------------------------------------------" << endl;
    cout << "  Test of the routine Ofs<T>& smult(Ofs<T> const& a, T const& c) " << endl;
    cout << "---------------------------------------------------" << endl;


    //Initialization
    Ofsd xd(OFS_ORDER);
    Ofsd xd2(OFS_ORDER);
    Ofsc xdc(OFS_ORDER);
    Ofsc xdc2(OFS_ORDER);
    cdouble res1, res2;
    cdouble resd1, resd2;
    //Arbitrary coefficients
    double c = 2.0;
    cdouble cd = 2.0 + 1.0*I;

    //Set random coefs in xd and xdc
    xd.set_random_coefs();
    xd2.set_random_coefs();
    xdc.set_random_coefs();
    xdc2.set_random_coefs();

    //Set results
    res1  =  xd.evaluate(VAR) + c*xd2.evaluate(VAR);
    resd1  = xdc.evaluate(VAR) + cd*xdc2.evaluate(VAR);

    //Take sum
    xd.ofs_smult(xd2, c);
    xdc.ofs_smult(xdc2, cd);

    //Set results
    res2  = xd.evaluate(VAR);
    resd2 = xdc.evaluate(VAR);

    //Comparison. double case
    cout << "1. double case.         " << endl;
    cout << "Expected result: " << creal(res1) << cimag(res1) << "i" << endl;
    cout << "Obtained result: " << creal(res2) << cimag(res2) << "i" << endl;
    res2 -= res1;
    cout << "Delta: " << cabs(res2) << endl;

    //Comparison. double case
    cout << "2. cdouble case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;
}

/**
 * \fn void ofs_mult()
 * \brief Test of the operator routine: Ofs<T>& mult(Ofs<T> const& a, T const& c)
 */
void ofs_test_mult()
{
    cout << "---------------------------------------------------" << endl;
    cout << "  Test of the routine Ofs<T>& mult(Ofs<T> const& a, T const& c) " << endl;
    cout << "---------------------------------------------------" << endl;


    //Initialization
    Ofsd xd(OFS_ORDER);
    Ofsd xd2(OFS_ORDER);
    Ofsc xdc(OFS_ORDER);
    Ofsc xdc2(OFS_ORDER);
    cdouble res1, res2;
    cdouble resd1, resd2;
    //Arbitrary coefficients
    double c = 2.0;
    cdouble cd = 2.0 + 1.0*I;

    //Set random coefs in xd and xdc
    xd.set_random_coefs();
    xd2.set_random_coefs();
    xdc.set_random_coefs();
    xdc2.set_random_coefs();

    //Set results
    res1  =  c*xd2.evaluate(VAR);
    resd1  = cd*xdc2.evaluate(VAR);

    //Take sum
    xd.ofs_mult(xd2, c);
    xdc.ofs_mult(xdc2, cd);

    //Set results
    res2  = xd.evaluate(VAR);
    resd2 = xdc.evaluate(VAR);

    //Comparison. double case
    cout << "1. double case.         " << endl;
    cout << "Expected result: " << creal(res1) << cimag(res1) << "i" << endl;
    cout << "Obtained result: " << creal(res2) << cimag(res2) << "i" << endl;
    res2 -= res1;
    cout << "Delta: " << cabs(res2) << endl;

    //Comparison. double case
    cout << "2. cdouble case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;
}


/**
 * \fn void ofs_fsum()
 * \brief Test of the operator routine: Ofs<T>& fsum(Ofs<T> const& a, T const& ma, Ofs<T> const& b, T const& mb)
 */
void ofs_test_fsum()
{

    cout << "---------------------------------------------------" << endl;
    cout << "  Test of the routine:" << endl;
    cout << "  Ofs<T>& fsum(Ofs<T> const& a, T const& ma, Ofs<T> const& b, T const& mb) " << endl;
    cout << "---------------------------------------------------" << endl;


    //Initialization
    Ofsd xd(OFS_ORDER);
    Ofsd xd2(OFS_ORDER);
    Ofsd xd3(OFS_ORDER);
    Ofsc xdc(OFS_ORDER);
    Ofsc xdc2(OFS_ORDER);
    Ofsc xdc3(OFS_ORDER);
    cdouble res1, res2;
    cdouble resd1, resd2;
    //Arbitrary coefficients
    double c = 2.0;
    double c2 = 3.0;
    cdouble cd = 2.0 + 1.0*I;
    cdouble cd2 = 1.0 + 3.0*I;

    //Set random coefs in xd and xdc
    xd.set_random_coefs();
    xd2.set_random_coefs();
    xd3.set_random_coefs();
    xdc.set_random_coefs();
    xdc2.set_random_coefs();
    xdc3.set_random_coefs();

    //Set results
    res1  =  c*xd2.evaluate(VAR) + c2*xd3.evaluate(VAR);
    resd1  = cd*xdc2.evaluate(VAR) + cd2*xdc3.evaluate(VAR);

    //Take sum
    xd.ofs_fsum(xd2, c, xd3, c2);
    xdc.ofs_fsum(xdc2, cd, xdc3, cd2);

    //Set results
    res2  = xd.evaluate(VAR);
    resd2 = xdc.evaluate(VAR);

    //Comparison. double case
    cout << "1. double case.         " << endl;
    cout << "Expected result: " << creal(res1) << cimag(res1) << "i" << endl;
    cout << "Obtained result: " << creal(res2) << cimag(res2) << "i" << endl;
    res2 -= res1;
    cout << "Delta: " << cabs(res2) << endl;

    //Comparison. double case
    cout << "2. cdouble case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;
}

/**
 * \fn void ofs_test_dot()
 * \brief Test of the routine: void dot().
 * Note that the Ofs objects are evaluated @ n*var, and not var,
 * to take into account the presence of the pulsation n.
 */
void ofs_test_dot()
{

    cout << "---------------------------------------------------" << endl;
    cout << "    Test of the routine:                           " << endl;
    cout << "    - Ofs<T>& dot(Ofs<T> const& a, double const& n)" << endl;
    cout << "    - Ofs<T>& dot(double const& n)                 " << endl;
    cout << "---------------------------------------------------" << endl;
    double epsilon = 1e-9;
    cout << " The time derivatives are estimated with a simple Euler scheme with epsilon = " << epsilon << "." << endl;

    //Initialization
    Ofsc xdc(OFS_ORDER);
    Ofsc xdc2(OFS_ORDER);
    cdouble resd1, resd2;
    double n  = 0.925195985520347;

    //Set random coefs in xd and xdc
    xdc.set_random_coefs();
    xdc2.set_random_coefs();

    //------------------------
    // Ofs<T>& dot(Ofs<T> const& a, double const& n)
    //------------------------
    cout << "---------------------" << endl;
    cout << "Ofs<T>& dot(Ofs<T> const& a, double const& n):         " << endl;
    cout << "---------------------" << endl;

    //Set results
    resd1 = (xdc2.evaluate(n*(VAR+epsilon))- xdc2.evaluate(n*(VAR)))/(epsilon);

    //Take dot
    xdc.dot(xdc2, n);

    //Set results
    resd2 = xdc.evaluate(n*VAR);

    //Comparison. double case
    cout << "Double complex case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;

    //------------------------
    // Ofs<T>& dot(Ofs<T> const& a)
    //------------------------
    cout << "---------------------" << endl;
    cout << "Ofs<T>& dot(Ofs<T> const& a):         " << endl;
    cout << "---------------------" << endl;


    //Set results
    resd1 = (xdc.evaluate(n*(VAR+epsilon))- xdc.evaluate(n*(VAR)))/(epsilon);

    //Take dot
    xdc.dot(n);

    //Set results
    resd2 = xdc.evaluate(n*VAR);

    //Comparison. double case
    cout << "Double complex case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;
}

/**
 * \fn void ofs_test_epspow()
 * \brief Test of the routine: Ofs<T>& epspow(Ofs<T> const& a, T const& alpha).
 */
void ofs_test_epspow()
{

    cout << "---------------------------------------------------" << endl;
    cout << "    Test of the routine:                           " << endl;
    cout << "    Ofs<T>& epspow(Ofs<T> const& a, T const& alpha)" << endl;
    cout << "---------------------------------------------------" << endl;
    cout << "This routine is only valid when" << endl;
    cout << "a is of the form \a = 1 + d  with  |d|_1 << 1." << endl;
    cout << "---------------------------------------------------" << endl;


    //Initialization
    Ofsc xdc(OFS_ORDER);
    Ofsc xdc2(OFS_ORDER);
    cdouble resd1, resd2;
    double n  = 0.925195985520347;
    double alpha = 2.0;

    //Set random coefs in xd and xdc
    xdc.set_random_coefs();
    xdc2.set_random_coefs();
    //Increase the condition |d| << 1
    xdc2.ofs_mult(xdc, 1e-2+0.0*I);
    //Put coef zero to 1.0
    xdc2.set_coef(1.0,0);

    double factor = cabs(xdc2.ofs_get_coef(1));
    cout <<  setw(5) << setprecision(1) << setiosflags(ios::scientific);
    cout << "Here, there is a factor ~" << factor << " between the order zero an one " << endl;
    cout << "and the order one." << endl;
    cout <<  setw(5) << setprecision(15) << std::showpos  <<  setiosflags(ios::scientific);

    cout << "1. Results for alpha = +2.0" << endl;
    //Set results
    resd1 = cpow(xdc2.evaluate(n*(VAR)), alpha+0.0*I);
    //Take dot
    xdc.ofs_epspow(xdc2, alpha+0.0*I);
    //Set results
    resd2 = xdc.evaluate(n*VAR);
    //Comparison. double case
    cout << "Double complex case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;

    alpha = -1.0;
    cout << "2. Results for alpha = -1.0" << endl;
    //Set results
    resd1 = cpow(xdc2.evaluate(n*(VAR)), alpha+0.0*I);
    //Take dot
    xdc.ofs_epspow(xdc2, alpha+0.0*I);
    //Set results
    resd2 = xdc.evaluate(n*VAR);
    //Comparison. double case
    cout << "Double complex case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;

}

/**
 * \fn void ofs_test_pows()
 * \brief Test of the routine: Ofs<T>& ofs_pows(Ofs<T> const& a, T const& alpha, Ofs<T> const& temp).
 */
void ofs_test_pows()
{
    cout << "----------------------------------------------------------------------" << endl;
    cout << "    Test of the power routine:                                        " << endl;
    cout << " Ofs<T>& ofs_pows(Ofs<T> const& a, T const& alpha, Ofs<T> const& temp)" << endl;
    cout << "----------------------------------------------------------------------" << endl;
    //Initialization
    Ofsc xdc(OFS_ORDER);
    Ofsc xdc2(OFS_ORDER);
    Ofsc temp(OFS_ORDER);
    cdouble resd1, resd2;
    double n  = 0.925195985520347;

    //Set random coefs in xd and xdc
    //xdc2.set_random_coefs();
    read_ofs_txt(xdc2, "test/test2_fft");

    double alpha = 2.0;
    cout << "1. Results for alpha = +2.0" << endl;
    //Set results
    resd1 = cpow(xdc2.evaluate(n*(VAR)), alpha+0.0*I);
    //Take power
    tic();
    xdc.ofs_pows(xdc2, alpha+0.0*I);
    double T = toc();
    cout << "Power via FFT in " << T << endl;
    //Set results
    resd2 = xdc.evaluate(n*VAR);
    //Comparison. double case
    cout << "Double complex case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;

    alpha = -1.0;
    cout << "2. Results for alpha = -1.0" << endl;
    //Set results
    resd1 = cpow(xdc2.evaluate(n*(VAR)), alpha+0.0*I);
    //Take power
    xdc.ofs_pows(xdc2, alpha+0.0*I);
    //Set results
    resd2 = xdc.evaluate(n*VAR);
    //Comparison. double case
    cout << "Double complex case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;

    cout << "-------------------------" << endl;
    cout << " Comparison with product:" << endl;
    cout << "-------------------------" << endl;
    tic();
    xdc.ofs_sprod(xdc, xdc2);
    double Tprod = toc();
    cout << "Product  in " << Tprod << endl;

}







