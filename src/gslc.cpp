#include "gslc.h"
using namespace std;

/**
 * \file gslc.cpp
 * \brief Additional operations on GSL objects.
 * \author BLB
 */


//----------------------------------------------------------------------------------------
// Matrix and vectors
//----------------------------------------------------------------------------------------
/**
 * \brief Transform a vector into a matrix with a given shift in the initial vector
 **/
void gslc_vec_to_mat(gsl_matrix* m, const double z[], int rows, int columns, int shift)
{
    int i,j,k;
    for(i=1; i<=rows ; i++)
        for(j=1; j<=columns; j++)
        {
            k = rows*(i-1)+j;
            gsl_matrix_set(m, i-1, j-1, z[shift+k-1]);
        }
}

/**
 * \brief Transform a vector into a complex matrix with a given shift in the initial vector
 **/
void gslc_vec_to_mat_complex(gsl_matrix_complex* m, const double z[],
                             int rows, int columns, int shift)
{
    int i,j,k;
    gsl_complex zc;
    for(i=1; i<=rows ; i++)
        for(j=1; j<=columns; j++)
        {
            k = rows*(i-1)+j;
            GSL_SET_COMPLEX(&zc, z[shift+k-1], 0.0);
            gsl_matrix_complex_set(m, i-1, j-1, zc);
        }
}

/**
 * \brief Transform a matrix into a vector with a given shift in the final vector
 **/
void gslc_mat_to_vec(double z[], const gsl_matrix* m, int rows, int columns, int shift)
{
    int i,j,k;
    for(i=1; i<=rows ; i++)
        for(j=1; j<=columns; j++)
        {
            k = rows*(i-1)+j;
            z[shift+k-1] = gsl_matrix_get(m, i-1, j-1);
        }
}

/**
 * \brief Transform a complex matrix into a vector with a given shift in the final vector
 **/
void gslc_mat_complex_to_vec(double z[], const gsl_matrix_complex*  m,
                             int rows, int columns, int shift)
{
    int i,j,k;
    for(i=1; i<=rows ; i++)
        for(j=1; j<=columns; j++)
        {
            k = rows*(i-1)+j;
            z[shift+k-1] = GSL_REAL(gsl_matrix_complex_get(m, i-1, j-1));
        }
}

//----------------------------------------------------------------------------------------
// Complex numbers
//----------------------------------------------------------------------------------------
/**
 * \brief Create a gsl_complex from (real, imag) in double format.
 **/
gsl_complex gslc_complex(double real, double imag)
{
    gsl_complex one_c;
    GSL_SET_COMPLEX(&one_c, real, imag);
    return one_c;
}

/**
 * \brief Create a gsl_complex from a complex number
 **/
gsl_complex gslc_complex(cdouble x)
{
    gsl_complex one_c;
    GSL_SET_COMPLEX(&one_c, creal(x), cimag(x));
    return one_c;
}

/**
 * \brief Create a complex number from a gsl_complex
 **/
cdouble gslc_complex(gsl_complex c)
{
    return GSL_REAL(c) + I*GSL_IMAG(c);
}

//----------------------------------------------------------------------------------------
// Printing a real matrix
//----------------------------------------------------------------------------------------
/**
 * \brief Print a real matrix
 **/
void gslc_matrix_printf(const gsl_matrix* m)
{
    int im = m->size1;
    int jm = m->size2;

    for(int i = 0; i < im ; i++)
    {
        for(int j=0; j <jm; j++) printf("%+5.15e ", gsl_matrix_get(m, i, j));
        printf("\n");
    }
}

/**
 * \brief Print an approximation of a real matrix: the coefficients under 1e-5 are displayed as null.
 **/
void gslc_matrix_printf_approx(const gsl_matrix* m)
{
    int im = m->size1;
    int jm = m->size2;
    double c;
    for(int i = 0; i < im ; i++)
    {
        for(int j=0; j <jm; j++)
        {
            c = gsl_matrix_get(m, i, j);
            if(fabs(c) > 1e-5) printf("%+1.3e ", c);
            else printf("%+1.3e ", 0.0);
        }
        printf("\n");
    }
}

/**
 * \brief Print a real matrix into a txt file.
 **/
void gslc_matrix_fprintf(const gsl_matrix* m, char* filename)
{

    FILE* f;
    f = fopen(filename, "w");
    int im = m->size1;
    int jm = m->size2;

    for(int i = 0; i < im ; i++)
    {
        for(int j=0; j <jm; j++) fprintf(f, "%+5.16e ", gsl_matrix_get(m, i, j));
        fprintf(f, "\n");
    }

    fclose(f);
}

//----------------------------------------------------------------------------------------
// Printing a complex matrix
//----------------------------------------------------------------------------------------
/**
 * \brief Print a complex matrix.
 **/
void gslc_matrix_complex_printf(const gsl_matrix_complex* m)
{

    int im = m->size1;
    int jm = m->size2;

    for(int i = 0; i < im ; i++)
    {
        for(int j=0; j <jm; j++)
        {
            printf("%+1.0e%+1.0ei ", GSL_REAL(gsl_matrix_complex_get(m, i, j)),
                   GSL_IMAG(gsl_matrix_complex_get(m, i, j)));
        }
        printf("\n");
    }
}

/**
 * \brief Print the real part of a complex matrix into a txt file.
 **/
void gslc_matrix_complex_printf_real(const gsl_matrix_complex* m)
{

    int im = m->size1;
    int jm = m->size2;

    for(int i = 0; i < im ; i++)
    {
        for(int j=0; j <jm; j++) printf("%+5.5e ", GSL_REAL(gsl_matrix_complex_get(m, i, j)));
        printf("\n");
    }
}

/**
 * \brief Print the imag part of a complex matrix into a txt file.
 **/
void gslc_matrix_complex_printf_imag(const gsl_matrix_complex* m)
{

    int im = m->size1;
    int jm = m->size2;

    for(int i = 0; i < im ; i++)
    {
        for(int j=0; j <jm; j++) printf("%+5.5e ", GSL_IMAG(gsl_matrix_complex_get(m, i, j)));
        printf("\n");
    }
}

/**
 * \brief Print an approximation of a complex matrix: the coefficients under 1e-5 are displayed as null.
 **/
void gslc_matrix_complex_printf_approx(const gsl_matrix_complex* m)
{
    int im = m->size1;
    int jm = m->size2;
    double c;

    for(int i = 0; i < im ; i++)
    {
        for(int j=0; j <jm; j++)
        {
            c = GSL_REAL(gsl_matrix_complex_get(m, i, j));
            if(fabs(c) > 1e-5) printf("%+1.3e ", c);
            else printf("%+1.3e ", 0.0);
        }
        printf("\n");
    }
}

/**
 * \brief Print a complex matrix into a txt file.
 **/
void gslc_matrix_complex_fprintf(const gsl_matrix_complex* m, char* filename)
{
    FILE* f;

    f = fopen(filename, "w");
    int im = m->size1;
    int jm = m->size2;

    for(int i = 0; i < im ; i++)
    {
        for(int j=0; j <jm; j++)
            fprintf(f, "%+5.16e %+5.16e ", GSL_REAL(gsl_matrix_complex_get(m, i, j)),
                    GSL_IMAG(gsl_matrix_complex_get(m, i, j)));
        fprintf(f, "\n");
    }

    fclose(f);
}

/**
 * \brief Print a complex matrix into a txt file, separating the real and imaginary part, in order to be easier to read.
 **/
void gslc_matrix_complex_fprintf_pretty(const gsl_matrix_complex* m, char* filename)
{
    FILE* f;
    f = fopen(filename, "w");
    int im = m->size1;
    int jm = m->size2;

    fprintf(f,"real part:\n");
    for(int i = 0; i < im ; i++)
    {
        for(int j=0; j <jm; j++)
            fprintf(f, "%+5.16e ", GSL_REAL(gsl_matrix_complex_get(m, i, j)));
        fprintf(f, "\n");
    }

    fprintf(f,"imag part:\n");
    for(int i = 0; i < im ; i++)
    {
        for(int j=0; j <jm; j++)
            fprintf(f, "%+5.16e ", GSL_IMAG(gsl_matrix_complex_get(m, i, j)));
        fprintf(f, "\n");
    }
    fclose(f);
}

/**
 * \brief Print the real part of a complex matrix into a txt file.
 **/
void gslc_matrix_complex_fprintf_real(const gsl_matrix_complex* m, char* filename)
{
    FILE* f;

    f = fopen(filename, "w");
    int im = m->size1;
    int jm = m->size2;

    for(int i = 0; i < im ; i++)
    {
        for(int j=0; j <jm; j++)
            fprintf(f, "%+2.0e ", GSL_REAL(gsl_matrix_complex_get(m, i, j)));
        fprintf(f, "\n");
    }

    fclose(f);
}

//----------------------------------------------------------------------------------------
// Printing a complex vector
//----------------------------------------------------------------------------------------
/**
 * \brief Print a complex vector into a txt file.
 **/
void gslc_vector_complex_fprintf(const gsl_vector_complex* z, char* filename)
{
    FILE* f;

    f = fopen(filename, "w");
    int im = z->size;

    for(int i = 0; i < im ; i++)
        fprintf(f, "%+5.15e  %+5.15e i\n", GSL_REAL(gsl_vector_complex_get(z, i)),
                GSL_IMAG(gsl_vector_complex_get(z, i)));
    fclose(f);
}
/**
 * \brief Print a complex vector.
 **/
void gslc_vector_complex_printf(const gsl_vector_complex* z)
{
    int im = z->size;

    for(int i = 0; i < im ; i++)
        printf("%+5.15e  %+5.15e i\n", GSL_REAL(gsl_vector_complex_get(z, i)),
               GSL_IMAG(gsl_vector_complex_get(z, i)));
}

//----------------------------------------------------------------------------------------
// Printing a vector
//----------------------------------------------------------------------------------------
/**
 * \brief Print a real vector.
 **/
void gslc_vector_printf(const gsl_vector* z)
{
    int im = z->size;
    for(int i = 0; i < im ; i++)
        printf("%+5.15e \n", gsl_vector_get(z, i));
}

//----------------------------------------------------------------------------------------
// Printing an eigensystem
//----------------------------------------------------------------------------------------
/**
 * \brief Print an eigensystem with the eigenvalues given as a complex vector.
 **/
void gslc_eigensystem_printf(gsl_vector_complex* eval, gsl_matrix_complex* evec, int k)
{
    gsl_vector_complex_view evec_i;
    gsl_complex eval_i;
    for (int i = 0; i < k; i++)
    {
        eval_i = gsl_vector_complex_get (eval, i);
        evec_i = gsl_matrix_complex_column (evec, i);
        printf ("eigenvalue = %+5.15e + %+5.15ei\n",
                GSL_REAL(eval_i), GSL_IMAG(eval_i));
        printf ("eigenvector = \n");
        for (int j = 0; j < 6; ++j)
        {
            gsl_complex z =
                gsl_vector_complex_get(&evec_i.vector, j);
            printf("%+5.15e + %+5.15ei\n", GSL_REAL(z), GSL_IMAG(z));
        }
    }
}

/**
 * \brief Print an eigensystem with the eigenvalues given as a complex matrix (on the diagonal).
 **/
void gslc_eigensystem_printf(gsl_matrix_complex* eval, gsl_matrix_complex* evec, int k)
{
    gsl_vector_complex_view evec_i;
    gsl_complex eval_i;
    for (int i = 0; i < k; i++)
    {
        eval_i = gsl_matrix_complex_get (eval, i, i);
        evec_i = gsl_matrix_complex_column (evec, i);
        printf ("eigenvalue = %+5.15e + %+5.15ei\n",
                GSL_REAL(eval_i), GSL_IMAG(eval_i));
        printf ("eigenvector = \n");
        for (int j = 0; j < 6; ++j)
        {
            gsl_complex z =
                gsl_vector_complex_get(&evec_i.vector, j);
            printf("%+5.15e + %+5.15ei\n", GSL_REAL(z), GSL_IMAG(z));
        }
    }
}

/**
 * \brief Print an eigensystem in a txt file with the eigenvalues given as a complex vector.
 **/
void gslc_eigensystem_fprintf(gsl_vector_complex* eval, gsl_matrix_complex* evec,
                              int k, char* filename)
{
    FILE* f;
    f = fopen(filename, "w");
    gsl_vector_complex_view evec_i;
    gsl_complex eval_i;
    for (int i = 0; i < k; i++)
    {
        eval_i = gsl_vector_complex_get (eval, i);
        evec_i = gsl_matrix_complex_column (evec, i);
        fprintf (f,"eigenvalue = %+5.15e + %+5.15ei\n",
                 GSL_REAL(eval_i), GSL_IMAG(eval_i));
        fprintf (f,"||eigenvalue|| = %+5.15e\n",
                 sqrt(GSL_REAL(eval_i)*GSL_REAL(eval_i) + GSL_IMAG(eval_i)*GSL_IMAG(eval_i)));
        fprintf (f,"eigenvector = \n");
        for (int j = 0; j < k; ++j)
        {
            gsl_complex z =
                gsl_vector_complex_get(&evec_i.vector, j);
            fprintf(f,"%+5.15e + %+5.15ei\n", GSL_REAL(z), GSL_IMAG(z));
        }
    }

    fclose(f);
}

/**
 * \brief Print an eigensystem in a txt file with the eigenvalues given as a complex matrix (on the diagonal).
 **/
void gslc_eigensystem_fprintf(gsl_matrix_complex* eval, gsl_matrix_complex* evec,
                              int k, char* filename)
{
    FILE* f;
    f = fopen(filename, "w");
    gsl_vector_complex_view evec_i;
    gsl_complex eval_i;
    for (int i = 0; i < k; i++)
    {
        eval_i = gsl_matrix_complex_get (eval, i, i);
        evec_i = gsl_matrix_complex_column (evec, i);
        fprintf (f,"eigenvalue = %+5.15e + %+5.15ei\n",
                 GSL_REAL(eval_i), GSL_IMAG(eval_i));
        fprintf (f,"eigenvector = \n");
        for (int j = 0; j < 6; ++j)
        {
            gsl_complex z =
                gsl_vector_complex_get(&evec_i.vector, j);
            fprintf(f,"%+5.15e + %+5.15ei\n", GSL_REAL(z), GSL_IMAG(z));
        }
    }
    fclose(f);
}

//----------------------------------------------------------------------------------------
// Reading GSL objects
//--------------------------------------------------------------------------------------
/**
 * \brief Read a complex vector from a txt file, obtained with the routine
 *        gslc_matrix_complex_fprintf(const gsl_matrix_complex *m, char* filename).
 **/
void glsc_matrix_complex_read(gsl_matrix_complex* m, string filename)
{
    int k1 = m->size1;
    int k2 = m->size2;
    //Init
    ifstream readStream;
    double cr, ci;

    //Reading
    cout << "glsc_matrix_complex_read. cr = " << filename+".txt" << endl;
    readStream.open((filename+".txt").c_str());
    for(int i = 0; i < k1 ; i++)
    {
        for(int j = 0; j< k2 ; j++)
        {
            readStream >> cr;  //real part
            readStream >> ci;  //imag part

            gsl_matrix_complex_set(m, i, j, gslc_complex(cr, ci));
        }
    }
    readStream.close();
}


//----------------------------------------------------------------------------------------
// Views of a matrix
//--------------------------------------------------------------------------------------
/**
 * \brief Set the kth column of the complex matrix m in the complex vector z using the GSL "views" structures.
 **/
void gslc_matrix_complex_column(gsl_vector_complex* z, gsl_matrix_complex* m, int k)
{
    gsl_vector_complex_view vecview;               //Vector views
    vecview = gsl_matrix_complex_column(m, k);     //Select the kth column of m
    gsl_vector_complex_memcpy(z, &vecview.vector); //Copy in z
}

/**
 * \brief Set the kth row of the complex matrix m in the complex vector z using the GSL "views" structures.
 **/
void gslc_matrix_complex_row(gsl_vector_complex* z, gsl_matrix_complex* m, int k)
{
    gsl_vector_complex_view vecview;               //Vector views
    vecview = gsl_matrix_complex_row(m, k);        //Select the kth column of m
    gsl_vector_complex_memcpy(z, &vecview.vector); //Copy in z
}

/**
 * \brief Set the complex vector z in the kth column of the complex matrix m using the GSL "views" structures.
 **/
void gslc_matrix_complex_column_V(gsl_matrix_complex* m, gsl_vector_complex* z, int k)
{
    gsl_vector_complex_view vecview;               //Vector views
    vecview = gsl_matrix_complex_column(m, k);     //Select the kth column of m
    gsl_vector_complex_memcpy(&vecview.vector, z); //Copy of z in the kth column of m
}


//----------------------------------------------------------------------------------------
// Misc manipulations
//----------------------------------------------------------------------------------------
/**
 * \brief Gives the infinity norm of a complex matrix m
 **/
double gslc_matrix_complex_infinity_norm(gsl_matrix_complex* m)
{
    int irow = m->size1;
    int icol = m->size2;
    //Mr = |m|
    gsl_matrix* Mr = gsl_matrix_calloc(irow, icol);
    for(int i = 0; i < irow; i++)
        for(int j = 0; j < icol; j++)
            gsl_matrix_set(Mr, i, j, gsl_complex_abs(gsl_matrix_complex_get(m, i, j)));
    //maxr = max (Mr)
    double maxr = gsl_matrix_max(Mr);
    //Memory release
    gsl_matrix_free(Mr);
    //Return maxr
    return maxr;
}

/**
 *  \brief Normalization: xm = xm/norm(xm)
 **/
void gslc_vector_complex_normalize(gsl_vector_complex* xm)
{
    gsl_complex VN1invC;
    double VN1 = gsl_blas_dznrm2 (xm);
    GSL_SET_COMPLEX(&VN1invC, 1/VN1, 0.0);
    gsl_vector_complex_scale (xm, VN1invC);
}

/**
 * \brief Isolate the real part of a complex matrix m: Mr = real(m)
 **/
void gslc_matrix_complex_real(gsl_matrix* Mr, gsl_matrix_complex* m)
{
    int irow = m->size1;
    int icol = m->size2;
    for(int i=0; i<irow; i++)
        for(int j=0; j<icol; j++)
            gsl_matrix_set(Mr, i, j,  GSL_REAL(gsl_matrix_complex_get(m, i, j)));
}

/**
 * \brief Copy the conjugate of a complex vector xm into xc.
 **/
void gslc_vector_complex_conjugate_memcpy(gsl_vector_complex* xc, gsl_vector_complex* xm)
{
    for(int i=0; i< (int)xm->size; i++)
        gsl_vector_complex_set(xc, i, gsl_complex_conjugate(gsl_vector_complex_get(xm, i)));
}


//----------------------------------------------------------------------------------------
// Specific routines for Monodromy and STM matrices manipulations
//----------------------------------------------------------------------------------------
/**
 * \brief Delete one raw and one column of a given square GSL matrix
 **/
gsl_matrix_complex* gslc_matrix_complex_deleteRC(gsl_matrix_complex* m, int k)
{
    if(k >= (int)m->size1)
    {
        cout << "Error in gslc_matrix_complex_deleteRC: the provided indix is out of scope." << endl;
        return NULL;
    }
    else if((int)m->size1 != (int)m->size2)
    {
        cout << "Error in gslc_matrix_complex_deleteRC: the matrix is not square" << endl;
        return NULL;
    }
    else
    {
        int im = (int)m->size1;
        gsl_matrix_complex* result = gsl_matrix_complex_calloc(im-1, im-1);

        if(k == 0)
        {
            gsl_matrix_complex_view M11 = gsl_matrix_complex_submatrix (m , 1 , 1 , im-1 , im-1);
            gsl_matrix_complex_memcpy(result, &M11.matrix);
        }
        else if(k==im-1)
        {
            gsl_matrix_complex_view M11 = gsl_matrix_complex_submatrix (m , 0 , 0 , im-1 , im-1);
            gsl_matrix_complex_memcpy(result, &M11.matrix);

        }
        else
        {
            gsl_matrix_complex_view M11 = gsl_matrix_complex_submatrix (m , 0   , 0   , k   , k  );
            gsl_matrix_complex_view M12 = gsl_matrix_complex_submatrix (m , 0   , k+1 , k   , im-1-k);
            gsl_matrix_complex_view M21 = gsl_matrix_complex_submatrix (m , k+1 , 0   , im-1-k , k  );
            gsl_matrix_complex_view M22 = gsl_matrix_complex_submatrix (m , k+1 , k+1 , im-1-k , im-1-k);


            gsl_matrix_complex_view R11 = gsl_matrix_complex_submatrix (result , 0 , 0 , k    , k   );
            gsl_matrix_complex_view R12 = gsl_matrix_complex_submatrix (result , 0 , k , k    , im-1-k);
            gsl_matrix_complex_view R21 = gsl_matrix_complex_submatrix (result , k , 0 , im-1-k , k   );
            gsl_matrix_complex_view R22 = gsl_matrix_complex_submatrix (result , k , k , im-1-k , im-1-k);

            gsl_matrix_complex_memcpy(&R11.matrix, &M11.matrix);
            gsl_matrix_complex_memcpy(&R12.matrix, &M12.matrix);
            gsl_matrix_complex_memcpy(&R21.matrix, &M21.matrix);
            gsl_matrix_complex_memcpy(&R22.matrix, &M22.matrix);

        }

        return result;
    }

}

/**
 * \brief Inverse transformation of a vector during Wielandt deflation: w = 1/vw*(w + vx/(vw-vx)*(z.w)*x)
 **/
void gslc_wielandt_inv_trans(gsl_vector_complex* w, gsl_complex vw,
                             gsl_vector_complex const* x, gsl_complex vx, gsl_vector_complex* z)
{
    gsl_complex lm;
    gsl_complex lm2;
    gsl_vector_complex* xc = gsl_vector_complex_calloc(x->size);

    //lm2 = vx/(vw-vx)
    lm2 = gsl_complex_sub(vw, vx);
    lm2 = gsl_complex_div(vx, lm2);
    //lm = vx/(vw-vx)*(z^T*w) = vx/(vw-vx)*(z.w)
    gsl_blas_zdotu(z , w , &lm);  //zdotu or zdotc??
    lm = gsl_complex_mul(lm, lm2);
    //xc = lm*x
    gsl_vector_complex_memcpy(xc, x);
    gsl_vector_complex_scale(xc, lm);
    //w = w - xc;
    gsl_vector_complex_add(w, xc);
    //w = 1/vw*w
    gsl_vector_complex_scale(w, gsl_complex_inverse(vw));
    //Normalisation
    gsl_vector_complex_scale(w, gslc_complex(1/gsl_blas_dznrm2 (w),0.0));
    //memory release
    gsl_vector_complex_free(xc);
}

/**
 *  \brief Inverse a symplectic complex matrix S0 into Sinv:
 *
 *  S0  =  S11 S12
 *         S21 S22
 *
 *  S-1 =  S22 -S12T
 *        -S21T S11T
 **/
void gslc_matrix_complex_symplectic_inverse(const gsl_matrix_complex* S0,
        gsl_matrix_complex* Sinv)
{
    int k = S0->size1/2;

    gsl_complex minus_one_c = gslc_complex(-1.0, 0.0);

    gsl_matrix_complex* S = gsl_matrix_complex_calloc (2*k, 2*k);
    gsl_matrix_complex_memcpy(S, S0);

    //--------------------------------------------------------------------------------------------------
    //Views of S
    //--------------------------------------------------------------------------------------------------
    gsl_matrix_complex_view S11 = gsl_matrix_complex_submatrix (S , 0 , 0 , k , k );   //S11
    gsl_matrix_complex_view S12 = gsl_matrix_complex_submatrix (S , 0 , k , k , k );   //S12
    gsl_matrix_complex_view S21 = gsl_matrix_complex_submatrix (S , k , 0 , k , k );   //S21
    gsl_matrix_complex_view S22 = gsl_matrix_complex_submatrix (S , k , k , k , k );   //S22

    //--------------------------------------------------------------------------------------------------
    //Views of Sinv
    //--------------------------------------------------------------------------------------------------
    gsl_matrix_complex_view Sinv11 = gsl_matrix_complex_submatrix (Sinv , 0 , 0 , k , k ); //Sinv11
    gsl_matrix_complex_view Sinv12 = gsl_matrix_complex_submatrix (Sinv , 0 , k , k , k ); //Sinv12
    gsl_matrix_complex_view Sinv21 = gsl_matrix_complex_submatrix (Sinv , k , 0 , k , k ); //Sinv21
    gsl_matrix_complex_view Sinv22 = gsl_matrix_complex_submatrix (Sinv , k , k , k , k ); //Sinv22


    //--------------------------------------------------------------------------------------------------
    // Storage
    //--------------------------------------------------------------------------------------------------
    //Sij_m = -S12
    gsl_matrix_complex* Sij_m = gsl_matrix_complex_calloc (k, k);
    gsl_matrix_complex_memcpy(Sij_m, &S12.matrix);
    gsl_matrix_complex_scale(Sij_m, minus_one_c);
    //Sinv12 = -S12^T = Sij_m^T
    gsl_matrix_complex_transpose_memcpy(&Sinv12.matrix, Sij_m);
    //Sij_m = -S21
    gsl_matrix_complex_memcpy(Sij_m, &S21.matrix);
    gsl_matrix_complex_scale(Sij_m, minus_one_c);
    //Sinv12 = -S21^T = Sij_m^T
    gsl_matrix_complex_transpose_memcpy(&Sinv21.matrix, Sij_m);
    //Sinv11 = S22^T
    gsl_matrix_complex_transpose_memcpy(&Sinv11.matrix, &S22.matrix);
    //Sinv22 = S11^T
    gsl_matrix_complex_transpose_memcpy(&Sinv22.matrix, &S11.matrix);

    //Memory release
    gsl_matrix_complex_free(Sij_m);
}

/**
 *  \brief Inverse a symplectic real matrix S0 into Sinv:
 *
 *  S0  =  S11 S12
 *         S21 S22
 *
 *  S-1 =  S22 -S12T
 *        -S21T S11T
 **/
void gslc_matrix_symplectic_inverse(const gsl_matrix* S0, gsl_matrix* Sinv)
{
    int k = S0->size1/2;
    gsl_matrix* S = gsl_matrix_calloc (2*k, 2*k);
    gsl_matrix_memcpy(S, S0);

    //--------------------------------------------------------------------------------------------------
    //Views of S
    //--------------------------------------------------------------------------------------------------
    gsl_matrix_view S11 = gsl_matrix_submatrix (S , 0 , 0 , k , k );   //S11
    gsl_matrix_view S12 = gsl_matrix_submatrix (S , 0 , k , k , k );   //S12
    gsl_matrix_view S21 = gsl_matrix_submatrix (S , k , 0 , k , k );   //S21
    gsl_matrix_view S22 = gsl_matrix_submatrix (S , k , k , k , k );   //S22

    //--------------------------------------------------------------------------------------------------
    //Views of Sinv
    //--------------------------------------------------------------------------------------------------
    gsl_matrix_view Sinv11 = gsl_matrix_submatrix (Sinv , 0 , 0 , k , k ); //Sinv11
    gsl_matrix_view Sinv12 = gsl_matrix_submatrix (Sinv , 0 , k , k , k ); //Sinv12
    gsl_matrix_view Sinv21 = gsl_matrix_submatrix (Sinv , k , 0 , k , k ); //Sinv21
    gsl_matrix_view Sinv22 = gsl_matrix_submatrix (Sinv , k , k , k , k ); //Sinv22


    //--------------------------------------------------------------------------------------------------
    // Storage
    //--------------------------------------------------------------------------------------------------
    //Sij_m = -S12
    gsl_matrix* Sij_m = gsl_matrix_calloc (k, k);
    gsl_matrix_memcpy(Sij_m, &S12.matrix);
    gsl_matrix_scale(Sij_m, -1.0);
    //Sinv12 = -S12^T = Sij_m^T
    gsl_matrix_transpose_memcpy(&Sinv12.matrix, Sij_m);
    //Sij_m = -S21
    gsl_matrix_memcpy(Sij_m, &S21.matrix);
    gsl_matrix_scale(Sij_m, -1.0);
    //Sinv12 = -S21^T = Sij_m^T
    gsl_matrix_transpose_memcpy(&Sinv21.matrix, Sij_m);
    //Sinv11 = S22^T
    gsl_matrix_transpose_memcpy(&Sinv11.matrix, &S22.matrix);
    //Sinv22 = S11^T
    gsl_matrix_transpose_memcpy(&Sinv22.matrix, &S11.matrix);
    //Memory release
    gsl_matrix_free(Sij_m);
}

/**
 *  \brief Transpose + conjugate the complex matrix S intro SH: SH = S^H != S^T
 **/
void gslc_matrix_complex_H_memcpy(gsl_matrix_complex* SH, const gsl_matrix_complex* S)
{
    int k = S->size1;
    //gsl_matrix_complex_transpose_memcpy(SH, S);
    for(int i =0; i<k; i++)
        for(int j=0; j<k; j++)
            gsl_matrix_complex_set(SH, i, j, gsl_complex_conjugate(gsl_matrix_complex_get(S, j, i)));
}

/**
 *  \brief Setting the 2n*2n complex matrix J equal to the fundamental symplectic matrix:
 *
 *  real(J) = | 0  In |   and imag(J) = 0
 *            |-In 0  |
 **/
void glsc_matrix_complex_set_J(gsl_matrix_complex* J)
{
    if(J->size1%2!=0)
    {
        cout << "glsc_matrix_complex_set_J. J is not of the form 2n*2n." << endl;
        return;
    }
    else
    {
        gsl_complex one_c  = gslc_complex(1.0, 0.0);
        gsl_complex minus_one_c = gslc_complex(-1.0, 0.0);
        int k = J->size1/2;

        gsl_matrix_complex_set_zero(J);
        for(int j = 0; j< k; j++)
        {
            gsl_matrix_complex_set(J, j, j+k, one_c);
            gsl_matrix_complex_set(J, j+k, j, minus_one_c);
        }
    }

}

/**
 * \brief Setting the 2n*2n complex matrix J equal to the fundamental symplectic matrix:
 *
 *  J = | 0  In |
 *      |-In 0  |
 **/
void glsc_matrix_set_J(gsl_matrix* J)
{
    if(J->size1%2!=0)
    {
        cout << "glsc_matrix_complex_set_J. J is not of the form 2n*2n." << endl;
        return;
    }
    else
    {
        int k = J->size1/2;

        gsl_matrix_set_zero(J);
        for(int j = 0; j< k; j++)
        {
            gsl_matrix_set(J, j, j+k,  1.0);
            gsl_matrix_set(J, j+k, j, -1.0);
        }
    }

}

/**
 * \brief Use the symmetry of the QBCP to compute the stable (resp. unstable) from the
 *        unstable (resp. stabl) eigenvector
 **/
void gslc_vector_complex_symcpy(gsl_vector_complex* vep2, gsl_vector_complex* vep1)
{
    //Opposite for i = 0, 2, 4
    gsl_vector_complex_set(vep2, 0, gsl_complex_negative(gsl_vector_complex_get(vep1, 0)));
    gsl_vector_complex_set(vep2, 2, gsl_complex_negative(gsl_vector_complex_get(vep1, 2)));
    gsl_vector_complex_set(vep2, 4, gsl_complex_negative(gsl_vector_complex_get(vep1, 4)));
    //Equal for i = 1, 3, 5
    gsl_vector_complex_set(vep2, 1, gsl_vector_complex_get(vep1, 1));
    gsl_vector_complex_set(vep2, 3, gsl_vector_complex_get(vep1, 3));
    gsl_vector_complex_set(vep2, 5, gsl_vector_complex_get(vep1, 5));
}

/**
 * \brief Matrix-vector product when the matrix is given as a product of matrices:
 *        ym = ten[1]...ten[m] * xm, with
 *        ym  a 6x1 vector
 *        xm  a 6x1 vector
 *        ten a 6x6xM matrix
 *        Careful: the product of matrices ten is indexed from 1 to m, and not 0 to m-1.
 **/
void gslc_matrix_vector_product(gsl_matrix_complex** ten, const gsl_vector_complex* xm,
                                gsl_vector_complex* ym, int k)
{
    gsl_complex one_c  = gslc_complex(1.0, 0.0);
    gsl_complex zero_c = gslc_complex(0.0, 0.0);

    gsl_vector_complex* ym1 = gsl_vector_complex_calloc(6);
    gsl_vector_complex* ym2 = gsl_vector_complex_calloc(6);

    //ym1 = xm
    gsl_vector_complex_memcpy(ym1, xm);
    for(int i =1; i<=k; i++)
    {
        //ym2 = ten[i]*ym1
        gsl_blas_zgemv (CblasNoTrans, one_c , ten[i], ym1 , zero_c , ym2);
        //ym1 = ym2
        gsl_vector_complex_memcpy(ym1, ym2);
    }
    //ym = ym2
    gsl_vector_complex_memcpy(ym, ym2);

    gsl_vector_complex_free(ym1);
    gsl_vector_complex_free(ym2);
}

/**
 * \brief Matrix-matrix product when the matrix is given as a product of matrices:
 *        ym = ten[1]...ten[m] * xm, with
 *        ym  a 6x6 matrix
 *        xm  a 6x6 matrix
 *        ten a 6x6xM matrix
 *        Careful: the product of matrices ten is indexed from 1 to m, and not 0 to m-1.
 **/
void gslc_matrix_matrix_product(gsl_matrix_complex** ten, const gsl_matrix_complex* xm,
                                gsl_matrix_complex* ym, int k)
{
    gsl_complex one_c  = gslc_complex(1.0, 0.0);
    gsl_complex zero_c = gslc_complex(0.0, 0.0);

    gsl_matrix_complex* ym1 = gsl_matrix_complex_calloc(6,6);
    gsl_matrix_complex* ym2 = gsl_matrix_complex_calloc(6,6);

    //ym1 = xm
    gsl_matrix_complex_memcpy(ym1, xm);
    for(int i =1; i<=k; i++)
    {
        //ym2 = ten[i]*ym1
        gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, one_c , ten[i], ym1 , zero_c , ym2);
        //ym1 = ym2
        gsl_matrix_complex_memcpy(ym1, ym2);
    }
    //ym = ym2
    gsl_matrix_complex_memcpy(ym, ym2);

    gsl_matrix_complex_free(ym1);
    gsl_matrix_complex_free(ym2);
}

/**
 * \brief Matrix inverse-vector product when the matrix is given as a product of matrices:
 *        ym = ten[m]^(-1)...ten[1]^(-1) * xm, with
 *        ym  a 6x1 vector
 *        xm  a 6x1 vector
 *        ten a 6x6xM matrix
 *        Careful: the product of matrices ten is indexed from 1 to m, and not 0 to m-1.
 **/
void gslc_matrix_vector_invproduct(gsl_matrix_complex** ten,  const gsl_vector_complex* xm,
                                   gsl_vector_complex* ym, int k)
{
    int s;
    gsl_permutation* p = gsl_permutation_alloc (6);
    gsl_matrix_complex* AUX = gsl_matrix_complex_calloc(6,6);
    gsl_vector_complex* ym1 = gsl_vector_complex_calloc(6);
    gsl_vector_complex* ym2 = gsl_vector_complex_calloc(6);

    //ym1 = xm
    gsl_vector_complex_memcpy(ym1, xm);
    for(int i =k; i>=1; i--)
    {
        //AUX = ten(K). Here, memcpy is mandatory, otherwise the original ten would be modified by the LU decomposition
        gsl_matrix_complex_memcpy(AUX, ten[i]);
        //LU decomposition of AUX
        gsl_linalg_complex_LU_decomp (AUX, p, &s);
        //BUX = AUX^(-1)*VEP
        gsl_linalg_complex_LU_solve (AUX , p , ym1 , ym2);
        //ym1 = ym2;
        gsl_vector_complex_memcpy(ym1, ym2);
    }
    //ym = ym2
    gsl_vector_complex_memcpy(ym, ym2);

    gsl_vector_complex_free(ym1);
    gsl_vector_complex_free(ym2);
    gsl_matrix_complex_free(AUX);
    gsl_permutation_free(p);
}

/**
 * \brief Initialize and return a complex tensor (an array of GSL complex matrices)
 *        CAREFUL: the array of matrices ten is shifted of one & so that the storage
 *        is easier in other routines (e.g. vepro)
 *        ==> ten has to be used from ten[1] to ten[m] and ten[0] is USELESS.
 **/
gsl_matrix_complex** gslc_matrix_complex_product_alloc(int size1, int size2, int k)
{
    gsl_matrix_complex** ten    = (gsl_matrix_complex**) malloc((k+1)*sizeof(gsl_matrix_complex*));

    for(int i = 0; i<= k; i++)
    {
        ten[i] = gsl_matrix_complex_calloc(size1,size2);
        gsl_matrix_complex_set_zero(ten[i]);
    }

    return ten;
}

/**
 * \brief Free a complex tensor (an array of GSL complex matrices)
 *        CAREFUL: the array of matrices ten is shifted of one & so that the storage
 *        is easier in other routines (e.g. vepro)
 *        ==> ten has to be used from ten[1] to ten[m] and ten[0] is USELESS.
 **/
void gslc_matrix_complex_product_free(gsl_matrix_complex** ten, int k)
{
    for(int i = 0; i<= k; i++) gsl_matrix_complex_free(ten[i]);
}



//----------------------------------------------------------------------------------------
// Allocation
//----------------------------------------------------------------------------------------
/**
 * \brief Initialize and return a tensor (an array of GSL matrices)
 *        NO shift in this routine, contrary to gslc_matrix_complex_product_alloc
 **/
gsl_matrix** gslc_matrix_array_alloc(int size1, int size2, int k)
{
    gsl_matrix** ten = (gsl_matrix**) malloc((k)*sizeof(gsl_matrix*));

    for(int i = 0; i< k; i++)
    {
        ten[i] = gsl_matrix_calloc(size1,size2);
        gsl_matrix_set_zero(ten[i]);
    }

    return ten;
}

/**
 * \brief Free a tensor (an array of GSL matrices)
 *        NO shift in this routine, contrary to gslc_matrix_complex_product_alloc
 **/
void gslc_matrix_array_free(gsl_matrix** ten, int k)
{
    for(int i = 0; i< k; i++) gsl_matrix_free(ten[i]);
}


/**
 * \brief Initialize and return an array of GSL matrices initialized at zero
 **/
gsl_matrix** gslc_matrix_array_calloc(int size1, int size2, int k)
{
    gsl_matrix** ten = (gsl_matrix**) malloc((k)*sizeof(gsl_matrix*));
    for(int i = 0; i< k; i++) ten[i] = gsl_matrix_calloc(size1,size2);
    return ten;
}


//----------------------------------------------------------------------------------------
// Utilitary routines on C vectors are gathered here
//----------------------------------------------------------------------------------------
/**
 *  \brief Prints an array of double using cout.
 **/
void vector_printf(double* z, int k)
{
    Config::configManager().coutmp();
    for(int i = 0; i < k; i ++)
        cout << i << "  " << z[i] << endl;
}

/**
 *  \brief Prints an array of complex double using cout.
 **/
void vector_complex_printf(cdouble* z, int k)
{
    Config::configManager().coutmp();
    for(int i = 0; i < k; i ++)
        cout << i << "  " << creal(z[i]) << "  " << cimag(z[i]) << endl;
}

/**
 *  Euclidian norm computed on the first k components of a complex vector:
 *         \f$ ENorm(z_0, k) = \left( \sum_{p = 0}^{k-1}  z_0[p] ^2 \right)^{-1/2} \f$
 **/
double ENorm(cdouble z0[], int k)
{
    double res = cabs(z0[0])*cabs(z0[0]);
    for(int p = 1; p < k; p++) res += cabs(z0[p])*cabs(z0[p]);
    return(sqrt(res));
}

/**
 *  Euclidian norm computed on the first k components of a double vector:
 *          \f$ ENorm(z_0, k) = \left( \sum_{p = 0}^{k-1}  z_0[p] ^2 \right)^{-1/2} \f$
 **/
double ENorm(double z0[], int k)
{
    double res = cabs(z0[0])*cabs(z0[0]);
    for(int p = 1; p < k; p++) res += cabs(z0[p])*cabs(z0[p]);
    return(sqrt(res));
}

