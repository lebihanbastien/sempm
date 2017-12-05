#ifndef GSLC_H_INCLUDED
#define GSLC_H_INCLUDED

/**
 * \file gslc.h
 * \brief Additional operations on GSL objects.
 * \author BLB
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
#include <math.h>

//GSL
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl_complex_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_eigen.h>

//Custom
#include "constants.h"
#include "config.h"



//----------------------------------------------------------------------------------------
// Matrix and vectors
//----------------------------------------------------------------------------------------
/**
 * \brief Transform a vector into a matrix with a given shift in the initial vector
 **/
void gslc_vec_to_mat(gsl_matrix* m, const double z[], int rows, int columns, int shift);

/**
 * \brief Transform a vector into a complex matrix with a given shift in the initial vector
 **/
void gslc_vec_to_mat_complex(gsl_matrix_complex* m, const double z[],
                             int rows, int columns, int shift);

/**
 * \brief Transform a matrix into a vector with a given shift in the final vector
 **/
void gslc_mat_to_vec(double z[], const gsl_matrix* m, int rows, int columns, int shift);

/**
 * \brief Transform a complex matrix into a vector with a given shift in the final vector
 **/
void gslc_mat_complex_to_vec(double z[], const gsl_matrix_complex*  m,
                             int rows, int columns, int shift);

//----------------------------------------------------------------------------------------
// Complex numbers
//----------------------------------------------------------------------------------------
/**
 * \brief Create a gsl_complex from (real, imag) in double format.
 **/
gsl_complex gslc_complex(double real, double imag);

/**
 * \brief Create a gsl_complex from a complex number
 **/
gsl_complex gslc_complex(cdouble x);

/**
 * \brief Create a complex number from a gsl_complex
 **/
cdouble gslc_complex(gsl_complex c);

//----------------------------------------------------------------------------------------
// Printing a real matrix
//----------------------------------------------------------------------------------------
/**
 * \brief Print a real matrix
 **/
void gslc_matrix_printf(const gsl_matrix* m);

/**
 * \brief Print an approximation of a real matrix: the coefficients under 1e-5 are displayed as null.
 **/
void gslc_matrix_printf_approx(const gsl_matrix* m);

/**
 * \brief Print a real matrix into a txt file.
 **/
void gslc_matrix_fprintf(const gsl_matrix* m, char* filename);

//----------------------------------------------------------------------------------------
// Printing a complex matrix
//----------------------------------------------------------------------------------------
/**
 * \brief Print a complex matrix.
 **/
void gslc_matrix_complex_printf(const gsl_matrix_complex* m);

/**
 * \brief Print the real part of a complex matrix into a txt file.
 **/
void gslc_matrix_complex_printf_real(const gsl_matrix_complex* m);

/**
 * \brief Print the imag part of a complex matrix into a txt file.
 **/
void gslc_matrix_complex_printf_imag(const gsl_matrix_complex* m);

/**
 * \brief Print an approximation of a complex matrix: the coefficients under 1e-5 are displayed as null.
 **/
void gslc_matrix_complex_printf_approx(const gsl_matrix_complex* m);

/**
 * \brief Print a complex matrix into a txt file.
 **/
void gslc_matrix_complex_fprintf(const gsl_matrix_complex* m, char* filename);

/**
 * \brief Print a complex matrix into a txt file, separating the real and imaginary part, in order to be easier to read.
 **/
void gslc_matrix_complex_fprintf_pretty(const gsl_matrix_complex* m, char* filename);

/**
 * \brief Print the real part of a complex matrix into a txt file.
 **/
void gslc_matrix_complex_fprintf_real(const gsl_matrix_complex* m, char* filename);

//----------------------------------------------------------------------------------------
// Printing a complex vector
//----------------------------------------------------------------------------------------
/**
 * \brief Print a complex vector into a txt file.
 **/
void gslc_vector_complex_fprintf(const gsl_vector_complex* z, char* filename);

/**
 * \brief Print a complex vector.
 **/
void gslc_vector_complex_printf(const gsl_vector_complex* z);

//----------------------------------------------------------------------------------------
// Printing a vector
//----------------------------------------------------------------------------------------
/**
 * \brief Print a real vector.
 **/
void gslc_vector_printf(const gsl_vector* z);

//----------------------------------------------------------------------------------------
// Printing an eigensystem
//----------------------------------------------------------------------------------------
/**
 * \brief Print an eigensystem with the eigenvalues given as a complex vector.
 **/
void gslc_eigensystem_printf(gsl_vector_complex* eval, gsl_matrix_complex* evec, int k);

/**
 * \brief Print an eigensystem with the eigenvalues given as a complex matrix (on the diagonal).
 **/
void gslc_eigensystem_printf(gsl_matrix_complex* eval, gsl_matrix_complex* evec, int k);

/**
 * \brief Print an eigensystem in a txt file with the eigenvalues given as a complex vector.
 **/
void gslc_eigensystem_fprintf(gsl_vector_complex* eval, gsl_matrix_complex* evec,
                              int k, char* filename);

/**
 * \brief Print an eigensystem in a txt file with the eigenvalues given as a complex matrix (on the diagonal).
 **/
void gslc_eigensystem_fprintf(gsl_matrix_complex* eval, gsl_matrix_complex* evec,
                              int k, char* filename);

//----------------------------------------------------------------------------------------
// Reading GSL objects
//--------------------------------------------------------------------------------------
/**
 * \brief Read a complex vector from a txt file, obtained with the routine
 *        gslc_matrix_complex_fprintf(const gsl_matrix_complex *m, char* filename).
 **/
void glsc_matrix_complex_read(gsl_matrix_complex* m, std::string filename);

//----------------------------------------------------------------------------------------
// Views of a matrix
//--------------------------------------------------------------------------------------
/**
 * \brief Set the kth column of the complex matrix m in the complex vector z using the GSL "views" structures.
 **/
void gslc_matrix_complex_column(gsl_vector_complex* z, gsl_matrix_complex* m, int k);

/**
 * \brief Set the kth row of the complex matrix m in the complex vector z using the GSL "views" structures.
 **/
void gslc_matrix_complex_row(gsl_vector_complex* z, gsl_matrix_complex* m, int k);

/**
 * \brief Set the complex vector z in the kth column of the complex matrix m using the GSL "views" structures.
 **/
void gslc_matrix_complex_column_V(gsl_matrix_complex* m, gsl_vector_complex* z, int k);

//----------------------------------------------------------------------------------------
// Misc manipulations
//----------------------------------------------------------------------------------------
/**
 * \brief Gives the infinity norm of a complex matrix m
 **/
double gslc_matrix_complex_infinity_norm(gsl_matrix_complex* m);

/**
 *  \brief Normalization: xm = xm/norm(xm)
 **/
void gslc_vector_complex_normalize(gsl_vector_complex* xm);

/**
 * \brief Isolate the real part of a complex matrix m: Mr = real(m)
 **/
void gslc_matrix_complex_real(gsl_matrix* Mr, gsl_matrix_complex* m);

/**
 * \brief Copy the conjugate of a complex vector xm into xc.
 **/
void gslc_vector_complex_conjugate_memcpy(gsl_vector_complex* xc, gsl_vector_complex* xm);


//----------------------------------------------------------------------------------------
// Specific routines for Monodromy and STM matrices manipulations
//----------------------------------------------------------------------------------------
/**
 * \brief Delete one raw and one column of a given square GSL matrix
 **/
gsl_matrix_complex* gslc_matrix_complex_deleteRC(gsl_matrix_complex* m, int k);

/**
 * \brief Inverse transformation of a vector during Wielandt deflation: w = 1/vw*(w + vx/(vw-vx)*(z.w)*x)
 **/
void gslc_wielandt_inv_trans(gsl_vector_complex* w, gsl_complex vw,
                             gsl_vector_complex const* x, gsl_complex vx, gsl_vector_complex* z);

/**
 *  \brief Inverse a symplectic complex matrix S0 into Sinv:
 *
 *  S0  =  S11 S12
 *         S21 S22
 *
 *  S-1 =  S22 -S12T
 *        -S21T S11T
 **/
void gslc_matrix_complex_symplectic_inverse(const gsl_matrix_complex* S0, gsl_matrix_complex* Sinv);

/**
 *  \brief Inverse a symplectic real matrix S0 into Sinv:
 *
 *  S0  =  S11 S12
 *         S21 S22
 *
 *  S-1 =  S22 -S12T
 *        -S21T S11T
 **/
void gslc_matrix_symplectic_inverse(const gsl_matrix* S0, gsl_matrix* Sinv);

/**
 *  \brief Transpose + conjugate the complex matrix S intro SH: SH = S^H != S^T
 **/
void gslc_matrix_complex_H_memcpy(gsl_matrix_complex* SH, const gsl_matrix_complex* S);

/**
 *  \brief Setting the 2n*2n complex matrix J equal to the fundamental symplectic matrix:
 *
 *  real(J) = | 0  In |   and imag(J) = 0
 *            |-In 0  |
 **/
void glsc_matrix_complex_set_J(gsl_matrix_complex* J);

/**
 * \brief Setting the 2n*2n complex matrix J equal to the fundamental symplectic matrix:
 *
 *  J = | 0  In |
 *      |-In 0  |
 **/
void glsc_matrix_set_J(gsl_matrix* J);

/**
 * \brief Use the symmetry of the QBCP to compute the stable (resp. unstable) from the
 *        unstable (resp. stabl) eigenvector
 **/
void gslc_vector_complex_symcpy(gsl_vector_complex* vep2, gsl_vector_complex* vep1);
/**
 * \brief Matrix-vector product when the matrix is given as a product of matrices:
 *        ym = ten[1]...ten[m] * xm, with
 *        ym  a 6x1 vector
 *        xm  a 6x1 vector
 *        ten a 6x6xM matrix
 *        Careful: the product of matrices ten is indexed from 1 to m, and not 0 to m-1.
 **/
void gslc_matrix_vector_product(gsl_matrix_complex** ten,
                                const gsl_vector_complex* xm, gsl_vector_complex* ym, int k);
/**
 * \brief Matrix-matrix product when the matrix is given as a product of matrices:
 *        ym = ten[1]...ten[m] * xm, with
 *        ym  a 6x6 matrix
 *        xm  a 6x6 matrix
 *        ten a 6x6xM matrix
 *        Careful: the product of matrices ten is indexed from 1 to m, and not 0 to m-1.
 **/
void gslc_matrix_matrix_product(gsl_matrix_complex** ten, const gsl_matrix_complex* xm, gsl_matrix_complex* ym, int k);

/**
 * \brief Matrix inverse-vector product when the matrix is given as a product of matrices:
 *        ym = ten[m]^(-1)...ten[1]^(-1) * xm, with
 *        ym  a 6x1 vector
 *        xm  a 6x1 vector
 *        ten a 6x6xM matrix
 *        Careful: the product of matrices ten is indexed from 1 to m, and not 0 to m-1.
 **/
void gslc_matrix_vector_invproduct(gsl_matrix_complex** ten,
                                   const gsl_vector_complex* xm, gsl_vector_complex* ym, int k);

/**
 * \brief Initialize and return a complex tensor (an array of GSL complex matrices)
 *        CAREFUL: the array of matrices ten is shifted of one & so that the storage
 *        is easier in other routines (e.g. vepro)
 *        ==> ten has to be used from ten[1] to ten[m] and ten[0] is USELESS.
 **/
gsl_matrix_complex** gslc_matrix_complex_product_alloc(int size1, int size2, int k);

/**
 * \brief Free a complex tensor (an array of GSL complex matrices)
 *        CAREFUL: the array of matrices ten is shifted of one & so that the storage
 *        is easier in other routines (e.g. vepro)
 *        ==> ten has to be used from ten[1] to ten[m] and ten[0] is USELESS.
 **/
void gslc_matrix_complex_product_free(gsl_matrix_complex** ten, int k);

//----------------------------------------------------------------------------------------
// Allocation
//----------------------------------------------------------------------------------------
/**
 * \brief Initialize and return a tensor (an array of GSL matrices)
 *        NO shift in this routine, contrary to gslc_matrix_complex_product_alloc
 **/
gsl_matrix** gslc_matrix_array_alloc(int size1, int size2, int k);

/**
 * \brief Free a tensor (an array of GSL matrices)
 *        NO shift in this routine, contrary to gslc_matrix_complex_product_alloc
 **/
void gslc_matrix_array_free(gsl_matrix** ten, int k);


/**
 * \brief Initialize and return an array of GSL matrices initialized at zero
 **/
gsl_matrix** gslc_matrix_array_calloc(int size1, int size2, int k);

//----------------------------------------------------------------------------------------
// Utilitary routines on C vectors are gathered here
//----------------------------------------------------------------------------------------
/**
 *  \brief Prints an array of double using cout.
 **/
void vector_printf(double* z, int k);

/**
 *  \brief Prints an array of complex double using cout.
 **/
void vector_complex_printf(cdouble* z, int k);

/**
 *  Euclidian norm computed on the first k components of a complex vector:
 *         \f$ ENorm(z_0, k) = \left( \sum_{p = 0}^{k-1}  z_0[p] ^2 \right)^{-1/2} \f$
 **/
double ENorm(cdouble z0[], int k);

/**
 *  Euclidian norm computed on the first k components of a double vector:
 *          \f$ ENorm(z_0, k) = \left( \sum_{p = 0}^{k-1}  z_0[p] ^2 \right)^{-1/2} \f$
 **/
double ENorm(double z0[], int k);

#endif // GSLC_H_INCLUDED
