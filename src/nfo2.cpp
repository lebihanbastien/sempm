/**
 * \file nfo2.cpp
 * \brief Computes the complete change of coordinates for
 *        the Normalized-Centered Hamiltonian of the QBCP (src). See the brief of the
 *        routine nfo2 for details.
 *
 *        WARNING: Note that the routine nfo2 has been tailored for the QBCP and surely
 *        requires a novel generic implementation so as to be able to tackle other
 *        dynamical models.
 *
 *        Moreover, from a pure coding p.o.v., and in that same perspective,
 *        the division of the code and naming of variables needs to be reviewed.
 *
 * \author BLB.
 */

#include "nfo2.h"


using namespace std;

#define CLEANED_FFT 1
#define ORIGINAL_FFT 0

// From C and Fortran
extern "C"
{
    #include "gnuplot_i.h"
    #include "nrutil.h"
    //void vepro_(double **DAT, int *M, double **VECS, int *NMX);
}


//----------------------------------------------------------------------------------------
// Main routine
//----------------------------------------------------------------------------------------
/**
 *  \brief Main routine. Stands for "normal form at order 2".
 *         Compute the complete change of coordinates for the QBCP about a
 *         given libration point, in order to get a diagonal form for the order 2 of the
 *         Hamiltonian.
 *         Results are stored in the folder data/COC.
 *
 *  We recall that the change of coordinates is of the form:
 *
 *  z = P(t)C zt + V(t) with zt the new set of coordinates and z the NC set of coordinates.
 *
 *  The following decompositions are computed as matrix/vector of Fourier series:
 *  - P
 *  - Q = inv(P)
 *  - FT11, FT12, FT21, FT22 that contains the c.o.c. for y_{theta} the (unused) momentum associated to theta.
 *  - G1, that contains V(t) = (g1, g2, 0, g3, g4, 0)^T (see code for details)
 *  - Xe, Xm, Xs, the translated positions of the Earth, Moon and Sun so that Xe[i] = real Xe[i] - V[i].
 *  See the modified potential expansion in BLB PhD manuscript for details. Note that:
 *  - Xe[0] contains the x position of the Earth
 *  - Xe[1] contains the y position of the Earth
 *  - Xe[2] does not contains the z position of the Earth (which is always zero)
 *    but contains sqrt(Xe[0]^2 + Xe[1]^2)
 *  - Xe[3:5] contains the UNtranslated positions of the Earth (xe, ye, ze).
 *
 * Note: G1 contains the gi coefficients with the form:
 *
 * G1 = | g1  g2 |
 *      | g3  g4 |
 *
 *  The idea behind this form is to be able to use the routines implemented for the matrices.
 *  This is a bit counterintuitive and should probably be changed in a future version of the code.
 *
 *  The results are stored in txt files.
 *  Each Fourer expansion is tested on a different grid than the one used for the FFT process
 *  For now, the test of the periodicity of P has been abandonned
 *  WARNING: to compute Q = inv(P), the following hypothesis is made
 *  (consistent with what mono_eigen produces):
 *
 *     | p11  p12   0   p14  p15    0  |
 *     | p21  p22   0   p24  p25    0  |
 * P = |  0    0   p33   0    0    p36 |
 *     | p41  p42   0   p44  p45    0  |
 *     | p51  p52   0   p54  p55    0  |
 *     |  0    0   p63   0    0    p66 |
 *
 *  WARNING: a lot of real matrices are manipulated as complex ones so that various
 *  scalar-matrix and vector-matrix operations can be made.
 */
void nfo2(FBPL &fbpl, int is_stored)
{
    cout << "nfo2 routine is called. " << endl;
    cout << "-------------------------------------------------------------" << endl;
    cout << "Compute the complete change of coordinates to:" << endl;
    cout << "- Get rid of order 1                          " << endl;
    cout << "- Get a normal form for the order 2           " << endl;
    cout << "of the Hamiltonian of the QBCP                " << endl;
    cout << "Results are stored in the folder data/COC/    " << endl;
    cout << "-------------------------------------------------------------" << endl;

    //------------------------------------------------------------------------------------
    //Model must be normalized.
    //------------------------------------------------------------------------------------
    if(!fbpl.is_norm)
    {
        cout << "nfo2. Warning: selected model is not normalized. Premature ending." << endl;
        return;
    }

    //------------------------------------------------------------------------------------
    //Plotting devices
    //------------------------------------------------------------------------------------
    gnuplot_ctrl  *h1;
    h1 = gnuplot_init();
    gnuplot_cmd(h1, "set grid");
    gnuplot_cmd(h1, "set title \" Dynamical equivalent of the libration point (synodical coordinates) \" ");

    //------------------------------------------------------------------------------------
    // Integration tools
    // Initial vector field include nonlinear variationnal equations
    // to get the periodic orbit
    //------------------------------------------------------------------------------------
    //System
    gsl_odeiv2_system sys;
    sys.function      = qbcp_vfn_varnonlin;
    sys.jacobian      = NULL;
    sys.dimension     = 42;
    sys.params        = &fbpl;
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    //Driver
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new (&sys, T,
                       Config::configManager().G_PREC_HSTART(),
                       Config::configManager().G_PREC_ABS(),
                       Config::configManager().G_PREC_REL());

    //====================================================================================
    // 1. Differential correction to get the dynamical equivalent of
    // the libration point (lpdyneq_single_shooting)
    //====================================================================================
    cout << "nfo2. Differential correction to get " << endl;
    cout << " the dynamical equivalent " << endl;
    cout << " of the libration point. " << endl;
    cout << "-------------------------------------------------------------" << endl;
    double y0[42];
    lpdyneq_single_shooting(d, y0);

    //Plotting the refined solution
    int Npoints = 5000;
    ode_plot_xy(y0, 42, +fbpl.us.T*0.5, d, h1, Npoints, 8, fbpl.is_norm, false, "", "");
    ode_plot_xy(y0, 42, -fbpl.us.T*0.5, d, h1, Npoints, 8, fbpl.is_norm, false, "", "");
    pressEnter(true, "Press ENTER to close the gnuplot window(s)");
    gnuplot_close(h1);

    //------------------------------------------------------------------------------------
    // Once the orbit is obtained, the 6th dimension vector field is kept
    // BUT the linearized translated variationnal equations are used
    // (qbcp_vfn_varlin_trans instead of qbcp_vfn_varnonlin).
    //------------------------------------------------------------------------------------
    //System
    gsl_odeiv2_system sys2;
    sys2.function      = qbcp_vfn_varlin_trans; //compare with "sys"
    sys2.jacobian      = NULL;
    sys2.dimension     = 42;
    sys2.params        = &fbpl;
    //Stepper
    const gsl_odeiv2_step_type *T2 = gsl_odeiv2_step_rk8pd;
    //Driver
    gsl_odeiv2_driver *d2 = gsl_odeiv2_driver_alloc_y_new (&sys2, T2,
                            Config::configManager().G_PREC_HSTART(),
                            Config::configManager().G_PREC_ABS(),
                            Config::configManager().G_PREC_REL());

    //------------------------------------------------------------------------------------
    // Initialization of various GSL objects
    // for next step (Monodromy matrix decomposition)
    //------------------------------------------------------------------------------------
    //MMc = DAT[M]*DAT[M-1]*...DAT[1]
    int M                       = 32;
    gsl_matrix_complex **DAT    = gslc_matrix_complex_product_alloc(6, 6, M);
    //MMc = MM = S*Dm*Sinv
    gsl_matrix_complex *S       = gsl_matrix_complex_calloc(6, 6);
    gsl_matrix_complex *MMc     = gsl_matrix_complex_calloc(6, 6);
    gsl_matrix_complex *Dm      = gsl_matrix_complex_calloc(6, 6);
    gsl_matrix* MM              = gsl_matrix_calloc(6, 6);
    //B = S*DB*Sinv
    gsl_matrix_complex *B       = gsl_matrix_complex_calloc(6, 6);
    //DB = log(Dm)/T
    gsl_matrix_complex *DB      = gsl_matrix_complex_calloc(6, 6);
    //Br = R*JB*Rinv
    gsl_matrix* JB              = gsl_matrix_calloc(6, 6);
    gsl_matrix* Br              = gsl_matrix_calloc(6, 6);
    gsl_matrix* R               = gsl_matrix_calloc(6, 6);

    //====================================================================================
    // 2. Decomposition of the Monodromy matrix through the use of custom routines
    // from a monodromy matrix given as a product of matrices
    //====================================================================================
    cout << "-------------------------------------------------------------" << endl;
    cout << "nfo2. Decomposition of the Monodromy matrix via mono_eigen.  " << endl;
    cout << "-------------------------------------------------------------" << endl;
    mono_eigen(d2, y0, fbpl.us.n, &fbpl, M, Csts::STABLE_DIR_POW, MMc, MM, Dm, DB, B, S, JB, Br, R, DAT, is_stored);

    //Change of Br to use B instead (usually better final precision).
    //Possible because we have ensured that B = real(B) inside mono_eigen.
    for(int i=0; i<6; i++) for(int j = 0; j<6; j++) gsl_matrix_set(Br, i, j, GSL_REAL(gsl_matrix_complex_get(B, i, j)));

    //Storage of the matrix Br (real form of B in params)
    gslc_mat_to_vec(fbpl.B, Br, 6, 6, 0);

    //------------------------------------------------------------------------------------
    // Change of derivative routines in the system sys2
    // in order to be able to compute P, given by
    //               P'(t) = Q(t)P(t) - P(t)*B
    //               P(0)  = Id
    //------------------------------------------------------------------------------------
    sys2.function = qbcp_vfn_varlin_trans_P;

    //------------------------------------------------------------------------------------
    // 3 Integration & storage of the COC
    //------------------------------------------------------------------------------------
    cout << "-------------------------------------------------------------" << endl;
    cout << "nfo2. Integration & storage of the COC                       " << endl;
    cout << "-------------------------------------------------------------" << endl;
    nfo2_coc(d2, y0, fbpl.us.n, &fbpl, Br, R, JB, (int) 8*OFS_ORDER, is_stored);

    if(is_stored)
    {
        cout << "-------------------------------------------------------------" << endl;
        cout << "nfo2. The COC has been stored in " << fbpl.cs.F_COC << endl;
        cout << "-------------------------------------------------------------" << endl;
    }
}



//----------------------------------------------------------------------------------------
// Monodromy and STM
//----------------------------------------------------------------------------------------
/**
 *  \brief Get the STM matrix at time \c t1, using GSL tools.
 *  \param y: the initial conditions in an array of \c double
 *  \param d: a pointer to a GSL driver for integration
 *  \param t1: the end time.
 *  \param M: a pointer to the gsl_matrix to update
 */
void stm_matrix_with_gsl(const double y[], gsl_odeiv2_driver *d, double t1, gsl_matrix* M)
{
    gsl_odeiv2_driver_reset(d);

    //Initial conditions
    double ys[42];
    for(int i=0; i<42; i++) ys[i] = y[i];

    //Final conditions
    double t = 0.0;

    //Integration
    gsl_odeiv2_driver_apply (d, &t, t1, ys);

    //Return the monodromy matrix
    gslc_vec_to_mat(M, ys, 6, 6 ,6);
}

/**
 *  \brief Get the STM matrix at time \c t1 in complex format, using GSL tools.
 *  \param y: the initial conditions in an array of \c double
 *  \param d: a pointer to a GSL driver for integration
 *  \param t1: the end time.
 *  \param M: a pointer to the gsl_matrix_complex to update
 */
void stm_complex(const double y[], gsl_odeiv2_driver *d, double t1, gsl_matrix_complex* M, double ys[])
{
    double t;
    //Initial conditions
    //double ys[42];
    for(int i=0; i<42; i++) ys[i]  = y[i];

    //Integration
    t = 0.0;
    gsl_odeiv2_driver_reset(d);
    gsl_odeiv2_driver_apply (d, &t, t1, ys);

    //Return the monodromy matrix
    gslc_vec_to_mat_complex(M, ys, 6, 6 ,6);
}

/**
 *  \brief Get the STM matrix on a \c M steps grid in complex format, using GSL tools.
 *  \param y: the initial conditions in an array of \c double
 *  \param d: a pointer to a GSL driver for integration
 *  \param t1: the end time.
 *  \param M: the number of steps
 *  \param DAT: a double pointer to the gsl_matrix_complex* array to update
 *
 *  CAREFUL: the matrix DAT is shifted of ONE & so that the storage is easier in vepro
 *  So, the DAT has to be used from DAT[1] to DAT[M] and DAT[0] is USELESS.
 *  Not a smart move, since vepro is obsolete for now (May 2015)...
 */
void stm_stepped_complex(const double y[], gsl_odeiv2_driver *d, double t1, int M, gsl_matrix_complex **DAT)
{
    gsl_odeiv2_driver_reset(d);
    gsl_matrix_complex *STM = gsl_matrix_complex_calloc(6,6);
    //CAREFUL: the matrix DAT is shifted of ONE & so that the storage is easier in vepro
    //SO the DAT has to be used from DAT[1] to DAT[M] and DAT[0] is USELESS.

    //Initial conditions
    double ys[42];
    for(int i=0; i < 42; i++) ys[i] = y[i];

    //Initial conditions
    double t = 0.0;
    double ti;
    //Loop on the grid [1...M]
    for(int i = 1; i<=M; i++)
    {
        ti = (double)i*t1/M;
        //gsl_odeiv2_driver_reset(d);
        //Integration
        gsl_odeiv2_driver_apply (d, &t, ti, ys);
        //Storage: DAT[i] = STM(ti)
        gslc_vec_to_mat_complex(STM, ys, 6, 6 ,6);
        gsl_matrix_complex_memcpy(DAT[i], STM);
        //The STM is restarded @ STM(ti) = eye(6);
        for(int j=6; j<42; j++) ys[j] = y[j];
    }
    return;
}

/**
 *  \brief Get the STM matrix on a \c M steps grid in complex format, using GSL tools. Alternative version that uses the symmetries of the native orbit.
 *  \param y: the initial conditions in an array of \c double
 *  \param d: a pointer to a GSL driver for integration
 *  \param t1: the end time.
 *  \param M: the number of steps
 *  \param DAT: a double pointer to the gsl_matrix_complex* array to update
 *
 *  CAREFUL: the matrix DAT is shifted of ONE & so that the storage is easier in vepro
 *  So, the DAT has to be used from DAT[1] to DAT[M] and DAT[0] is USELESS.
 *  Not a smart move, since vepro is obsolete for now (May 2015)...
 */
void stm_stepped_complex_alt(const double y[], gsl_odeiv2_driver *d, double t1, int M, gsl_matrix_complex **DAT)
{
    gsl_odeiv2_driver_reset(d);
    gsl_matrix_complex *STM = gsl_matrix_complex_calloc(6,6);
    //CAREFUL: the matrix DAT is shifted of ONE & so that the storage is easier in vepro
    //SO the DAT has to be used from DAT[1] to DAT[M] and DAT[0] is USELESS.

    //Initial conditions on the grid [1...M]
    double ic[M+1][42];
    //Error at a given point on the grid
    double dy;
    //Current state in the integration
    double ys[42];
    //Identity matrix set at all points
    for(int k = 0; k<=M; k++) for(int i=6; i<42; i++) ic[k][i] = y[i];

    //-----------------------------
    //Loop from 0 to T/2
    //-----------------------------
    //Initial conditions @t = 0
    for(int i=0; i<42; i++) ys[i]    = y[i];
    for(int i=0; i< 6; i++) ic[0][i] = y[i];
    double t = 0.0;
    double ti;
    for(int i = 1; i<=M/2; i++)
    {
        ti = (double)i*t1/M;
        //Integration
        gsl_odeiv2_driver_apply (d, &t, ti, ys);
        //Storage: DAT[i] = STM(ti)
        gslc_vec_to_mat_complex(STM, ys, 6, 6 ,6);
        gsl_matrix_complex_memcpy(DAT[i], STM);
        //Initial condition along the path
        for(int j=0; j< 6; j++) ic[i][j] = ys[j];
        //The STM is restarded @ STM(ti) = eye(6);
        for(int j=6; j<42; j++) ys[j] = y[j];
    }

    //-----------------------------
    //Loop from 0 to -T/2
    //-----------------------------
    //Initial conditions again in ys @ t= 0
    gsl_odeiv2_driver_reset(d);
    for(int i=0; i<42; i++) ys[i] = y[i];
    t = 0.0;
    d->h *= -1.0;
    for(int i = 1; i<=M/2; i++)
    {
        ti = -(double)i*t1/M;
        //Integration
        gsl_odeiv2_driver_apply (d, &t, ti, ys);
        //Initial condition along the path, in reverse order
        if(i < M/2) for(int j=0; j< 6; j++) ic[M-i][j] = ys[j];
        //no STM storage here
    }
    d->h *= -1.0;

    dy = 0.0;
    for(int j=0; j< 6; j++) dy+= (ic[M/2][j]-ys[j])*(ic[M/2][j]-ys[j]);


    //-----------------------------
    //Loop from T/2 to T
    //with IC from previous loop
    //-----------------------------
    //Initial conditions again in ys @ t = T/2
    gsl_odeiv2_driver_reset(d);
    for(int i=0; i<42; i++) ys[i] = ic[M/2][i];
    t = 0.5*t1;
    for(int i = 1; i<=M/2; i++)
    {
        ti = (double)(0.5*M+i)*t1/M;
        //Integration
        gsl_odeiv2_driver_apply (d, &t, ti, ys);
        //Storage: DAT[i] = STM(ti)
        gslc_vec_to_mat_complex(STM, ys, 6, 6 ,6);
        gsl_matrix_complex_memcpy(DAT[M/2+i], STM);
        dy = 0.0;
        for(int j=0; j< 6; j++) dy+= (ic[M/2+i][j]-ys[j])*(ic[M/2+i][j]-ys[j]);
        //Initial condition along the path
        for(int j=0; j< 6; j++) ys[j] = ic[M/2+i][j];
        //The STM is restarded @ STM(ti) = eye(6);
        for(int j=6; j<42; j++) ys[j] = y[j];
    }
    return;
}

/**
 *  \brief Computing the Monodromy matrix from the stepped STM stored in DAT[1, ..., M]: MMc = DAT[M]*DAT[M-1]*...*DAT[1]
 *  \param MMc: the matrix to update.
 *  \param DAT: a double pointer to the gsl_matrix_complex* array to update
 *  \param M: the number of steps
 **/
void mono_from_prod(gsl_matrix_complex* MMc, gsl_matrix_complex** DAT, int M)
{
    //GSL complex constants (1.0 and 0.0)
    gsl_complex one_c  = gslc_complex(1.0, 0.0);
    gsl_complex zero_c = gslc_complex(0.0, 0.0);
    //Temporary matrix
    gsl_matrix_complex *Mspare = gsl_matrix_complex_calloc(6, 6);
    //MMc = DAT[1]
    gsl_matrix_complex_memcpy(MMc, DAT[1]);
    //Loop: MMc = DAT[i]*DAT[i-1]*...*DAT[1]
    for(int i =2; i<=M; i++)
    {
        gsl_matrix_complex_memcpy(Mspare, MMc);
        gsl_blas_zgemm (CblasNoTrans,CblasNoTrans, one_c, DAT[i] , Mspare , zero_c , MMc);
    }
    //Memory release
    gsl_matrix_complex_free(Mspare);
}

/**
 *  \brief Direct diagonalization of the monodromy matrix
 *  \param MM    the monodromy matrix in real format
 *  \param evecd output: the eigenvectors in matrix format
 *  \param evald output: the eigenvalues  in vector format
 *  \param evmd  output: the eigenvalues  in matrix format (along the diagonal)
 **/
void mono_diag(gsl_matrix* MM, gsl_matrix_complex *evecd, gsl_vector_complex *evald, gsl_matrix_complex *evmd)
{
    //Init
    gsl_eigen_nonsymmv_workspace * wr = gsl_eigen_nonsymmv_alloc(6);  //workspace for eigenspaces determination
    gsl_matrix *AUX = gsl_matrix_calloc(6, 6);                         //auxiliary matrix
    //Diagonalization
    gsl_matrix_memcpy(AUX, MM);
    gsl_eigen_nonsymmv (AUX, evald, evecd, wr);
    gsl_eigen_nonsymmv_sort (evald, evecd, GSL_EIGEN_SORT_ABS_ASC);
    //Set evald in matrix format in evmd
    gsl_matrix_complex_set_zero(evmd);
    for(int i = 0; i <6; i++) gsl_matrix_complex_set(evmd, i, i, gsl_vector_complex_get(evald, i));
    //Memory release
    gsl_eigen_nonsymmv_free(wr);
    gsl_matrix_free(AUX);
}

//-----------------------------------------------------------------------------------------
// Manipulation of the matrix S, which contains the eigenvectors of the monodromy matrix,
// in columns
//----------------------------------------------------------------------------------------
/**
 *  \brief Normalisation of the matrix S in order to get a symplectic matrix.
 *
 *         Theory shows that
 *                    BUX = | B1  0  |   with B1 diagonal.
 *                          | 0   B1 |
 *         Normalization to obtain a true symplectic matrix S:
 *                  S = S x | B1  0 |^(-1)    then S^T x J x S = J
 *                          | 0   I |
 *
 *      See "Lectures on celestial mechanics", Siegel & Moser, p. 101
 **/
void norm_matrix_S(gsl_matrix_complex *S)
{
    //---------------------------------------------------------------------------------
    // Initialization of GSL objects
    //---------------------------------------------------------------------------------
    gsl_matrix_complex *ST   = gsl_matrix_complex_calloc(6, 6);  //Transpose (not hermitian inverse) of the matrix S
    gsl_matrix_complex *J    = gsl_matrix_complex_calloc(6, 6);  //J default symplectic matrix
    gsl_matrix_complex *AUX  = gsl_matrix_complex_calloc(6, 6);  //auxiliary matrix
    gsl_matrix_complex *BUX  = gsl_matrix_complex_calloc(6, 6);  //auxiliary matrix
    gsl_complex one_c  = gslc_complex(1.0, 0.0);                //1.0
    gsl_complex zero_c = gslc_complex(0.0, 0.0);                //0.0
    gsl_vector_complex_view vecview;                            //Vector view for vector selection in matrix

    //---------------------------------------------------------------------------------
    // BUX = S^T J S
    // Theory shows that
    //                  BUX = | B1  0  |   with B1 diagonal.
    //                        | 0   B1 |
    //---------------------------------------------------------------------------------
    //ST = transpose(S)
    gsl_matrix_complex_transpose_memcpy(ST, S);  //transpose only
    //Matrix J =  fundamental symplectic matrix
    glsc_matrix_complex_set_J(J);
    //AUX = J*S
    gsl_blas_zgemm (CblasNoTrans , CblasNoTrans , one_c , J , S , zero_c , AUX );
    //BUX = ST*AUX = ST*J*S
    gsl_blas_zgemm (CblasNoTrans , CblasNoTrans , one_c , ST , AUX , zero_c , BUX );

    //---------------------------------------------------------------------------------
    // Normalization to obtain a true symplectic matrix S:
    //
    //      S = S x | B1  0 |^(-1)    then S^T x J x S = J
    //              | 0   I |
    //---------------------------------------------------------------------------------
    gsl_complex coefc;
    for(int i = 0; i < 3 ; i++)
    {
        coefc = gsl_complex_pow_real(gsl_matrix_complex_get(BUX, i, i+3), -1.0/2);
        vecview = gsl_matrix_complex_column (S, i);
        gsl_vector_complex_scale(&vecview.vector, coefc);
        vecview = gsl_matrix_complex_column (S, i+3);
        gsl_vector_complex_scale(&vecview.vector, coefc);
    }

    gsl_matrix_complex_free(J);
    gsl_matrix_complex_free(AUX);
    gsl_matrix_complex_free(BUX);
    gsl_matrix_complex_free(ST);

}

/**
 *  \brief Prenormalisation of the matrix S. Used to enforce a decoupling between the real and imag parts
 *         of the vectors that bear the ELLIPTIC motion.
 *         Mainly here for reference, not used in current implementation. Use with care, hard coded values inside!
 **/
void prenorm_matrix_S(gsl_matrix_complex *S)
{
    gsl_vector_complex_view vecview;
    gsl_complex coefc;
    cout << "Moreover, you can force a decoupling between the real and imag parts" << endl;
    cout << "of the vectors that bear the ELLIPTIC motion, IF it's not already the case." << endl;
    cout << "This decoupling will lead to additional symmetries within the final matrix P(t)." << endl;
    cout << "If you want so, please enter 1" << endl;
    int ch;
    scanf("%d",&ch);

    if(ch == 1)
    {
        //--------------
        // Here can use hard-coded values because the user has guaranteed that the vertical motion is on the 3rd and sixth columns.
        // BUT be careful, not robust AT ALL to a change of system
        //--------------

        //Z motion
        coefc = gslc_complex(GSL_REAL(gsl_matrix_complex_get(S, 2,2)), -GSL_IMAG(gsl_matrix_complex_get(S, 2,2)));
        vecview = gsl_matrix_complex_column (S, 2);
        gsl_vector_complex_scale(&vecview.vector, coefc);

        coefc = gslc_complex(GSL_REAL(gsl_matrix_complex_get(S, 2,5)), -GSL_IMAG(gsl_matrix_complex_get(S, 2,5)));
        vecview = gsl_matrix_complex_column (S, 5);
        gsl_vector_complex_scale(&vecview.vector, coefc);


        //XY motion
        coefc = gslc_complex(GSL_REAL(gsl_matrix_complex_get(S, 0,0)), -GSL_IMAG(gsl_matrix_complex_get(S, 0,0)));
        vecview = gsl_matrix_complex_column (S, 0);
        gsl_vector_complex_scale(&vecview.vector, coefc);

        coefc = gslc_complex(GSL_REAL(gsl_matrix_complex_get(S, 0,3)), -GSL_IMAG(gsl_matrix_complex_get(S, 0,3)));
        vecview = gsl_matrix_complex_column (S, 3);
        gsl_vector_complex_scale(&vecview.vector, coefc);


        //--------------
        // Visual check
        //--------------
        cout << "real(S)" << endl;
        gslc_matrix_complex_printf_real(S);
        cout << "imag(S):" << endl;
        gslc_matrix_complex_printf_imag(S);
        cout << "Please press Enter to continue." << endl;
        char chc;
        scanf("%c",&chc);
    }
}

/**
 *  \brief Update the matrices S and Dm using the center directions in evecr, evalr, and the hyperbolic directions in eigenVu/Vs.
 *         Before this step:
 *              - evecr contains the center eigenvectors of the mon. matrix, at positions (1,2,3,4): vxy, conj(vxy), vz, conj(vz)
 *              - evalr contains the center eigenvalues  of the mon. matrix, at positions (1,2,3,4)
 *              - eigenVu/Lu contains the unstable couple of the mon. matrix.
 *              - eigenVs/Ls contains the stable couple of the mon. matrix.
 *         After this step:
 *              - S contains the eigenvector in the following order:
 *                                          (vxy, vu, vz, conj(vxy), vs, conj(vz))
 *              - Dm contains the corresponding eigenvalues.
 **/
void update_matrix_S(gsl_matrix_complex *S,
             gsl_matrix_complex *Dm,
             gsl_matrix_complex *evecr,
             gsl_vector_complex *evalr,
             gsl_vector_complex* eigenVu,
             gsl_vector_complex* eigenVs,
             gsl_complex eigenLu,
             gsl_complex eigenLs)
{
    //------------------------------------------------------------
    // At this step:
    // - evecr contains the center eigenvectors of MMc, at positions (1,2,3,4)
    // - evalr contains the center eigenvalues  of MMc, at positions (1,2,3,4)
    // - eigenVu/Lu contains the unstable couple of MMc.
    // - eigenVs/Ls contains the stable couple of MMc.
    // (As in the Power+Wielandt case).
    //------------------------------------------------------------
    //The permutation logic in S and Dm
    int keymap[6];
    keymap[0] = 5;
    keymap[1] = 2;
    keymap[2] = 1;
    keymap[3] = 4;
    keymap[4] = 3;
    keymap[5] = 0;

    //Perform permutation
    permutation_matrix_S(S, Dm, evecr, evalr, eigenVu, eigenVs, eigenLu, eigenLs, keymap);

    //Check that the vertical center is well placed.
    gsl_vector_complex_view V2;
    V2 = gsl_matrix_complex_column (S, 2);
    double v1, v2;

    v1 = gsl_complex_abs(gsl_vector_complex_get(&V2.vector, 0));  //first component of V(2)
    v2 = gsl_complex_abs(gsl_vector_complex_get(&V2.vector, 1));

    if(!(v1 < 1e-8 && v2 < 1e-8)) //if V2 does not bear the vertical motion
    {
        //New logic
        keymap[0] = 3;
        keymap[1] = 0;
        keymap[2] = 1;
        keymap[3] = 4;
        keymap[4] = 2;
        keymap[5] = 5;

        //Perform permutation
        permutation_matrix_S(S, Dm, evecr, evalr, eigenVu, eigenVs, eigenLu, eigenLs, keymap);
    }

    cout << "update_matrix_S. Check of the form of S:" << endl;
    cout << "real(S)" << endl;
    gslc_matrix_complex_printf_real(S);
    cout << "imag(S):" << endl;
    gslc_matrix_complex_printf_imag(S);
    cout << "The eigenvectors which bear the vertical motion" << endl;
    cout << "should appear in columns 3 and 6." << endl;
    cout << "Please press Enter to continue." << endl;
    char chc;
    scanf("%c",&chc);

    //------------------------------------------------------------
    // At this step:
    // - S contains the eigenvectors in the following order:
    //      (vxy, vu, vz, conj(vxy), vs, conj(vz))
    // - Dm contains the corresponding eigenvalues.
    //------------------------------------------------------------
}

/**
 *  \brief Perform a permutation within S and Dm.
 *
 *  Storage of:
 *  - The hyperbolic directions (input in eigenVu/Lu and eigenVs/Ls)
 *  - The central directions (input in evecr, evalr)
 *  Used in mono_eigen to ensure that the eigenvectors
 *  that bear the vertical motion are at the right position in S (and Dm)
 *
 */
void permutation_matrix_S(gsl_matrix_complex *S,
                  gsl_matrix_complex *Dm,
                  gsl_matrix_complex *evecr,
                  gsl_vector_complex *evalr,
                  gsl_vector_complex *eigenVu,
                  gsl_vector_complex *eigenVs,
                  gsl_complex eigenLu,
                  gsl_complex eigenLs,
                  int *keymap)
{
    int A0, B0;
    gsl_complex evalA, evalB;
    gsl_vector_complex_view evacA, evecB, S_V;

    //Center 1
    //-------------------------------------------
    A0 = keymap[0];
    B0 = keymap[1];
    evalA = gsl_vector_complex_get (evalr, 2);
    evalB = gsl_vector_complex_get (evalr, 3);
    evacA = gsl_matrix_complex_column (evecr, 2);
    evecB = gsl_matrix_complex_column (evecr, 3);
    //Set the eigenvalues
    gsl_matrix_complex_set(Dm , A0 , A0, evalA);
    gsl_matrix_complex_set(Dm , B0 , B0, evalB);
    //Set the eigenvectors
    S_V = gsl_matrix_complex_column (S, A0);
    gsl_vector_complex_memcpy(&S_V.vector, &evacA.vector);
    S_V = gsl_matrix_complex_column (S, B0);
    gsl_vector_complex_memcpy(&S_V.vector, &evecB.vector);

    //Saddle
    //-------------------------------------------
    A0 = keymap[2];
    B0 = keymap[3];
    //Unstable dir.
    gsl_matrix_complex_set(Dm , A0 , A0, eigenLu);
    gsl_matrix_complex_set_col (S , A0, eigenVu);

    //Stable dir.
    gsl_matrix_complex_set(Dm , B0 , B0, eigenLs);
    gsl_matrix_complex_set_col (S , B0, eigenVs);


    //Center 2
    //-------------------------------------------
    //WARNING: the center 2 should correspond to the vertical components (z, pz)
    //Decoupled from the rest of the dynamics
    A0 = keymap[4];
    B0 = keymap[5];
    evalA = gsl_vector_complex_get (evalr, 4);
    evalB = gsl_vector_complex_get (evalr, 5);
    evacA = gsl_matrix_complex_column (evecr, 4);
    evecB = gsl_matrix_complex_column (evecr, 5);
    //Set the eigenvalues
    gsl_matrix_complex_set(Dm , A0 , A0, evalA);
    gsl_matrix_complex_set(Dm , B0 , B0, evalB);
    //Set the eigenvectors
    S_V = gsl_matrix_complex_column (S, A0);
    gsl_vector_complex_memcpy(&S_V.vector, &evacA.vector);
    S_V = gsl_matrix_complex_column (S, B0);
    gsl_vector_complex_memcpy(&S_V.vector, &evecB.vector);
}

/**
 *  \brief Checks that R is symplectic (real version of S). If not, performs an additional permutation/normalization on the couple (S,Dm)
 **/
void symplecticity_test_for_R(gsl_matrix_complex *S, gsl_matrix_complex *Dm, gsl_matrix *R)
{
    gsl_matrix *J     = gsl_matrix_calloc(6, 6);
    gsl_matrix *RH    = gsl_matrix_calloc(6, 6);
    gsl_matrix *AUX   = gsl_matrix_calloc(6, 6);
    gsl_matrix *BUX   = gsl_matrix_calloc(6, 6);

    gsl_vector_complex_view evecA, evecB;
    gsl_complex evalA, evalB;
    gsl_vector_complex *AUXc  = gsl_vector_complex_calloc(6);  //auxiliary vector

    //---------------------------------------------------------------
    //  Update R from S columns
    //---------------------------------------------------------------
    for(int i=0; i<6; i++)
    {
        gsl_matrix_set(R, i, 0 , sqrt(2)*GSL_REAL(gsl_matrix_complex_get(S, i, 0))); //sqrt(2)*Re(c2)
        gsl_matrix_set(R, i, 1 , GSL_REAL(gsl_matrix_complex_get(S, i, 1)));         //c1
        gsl_matrix_set(R, i, 2 , sqrt(2)*GSL_REAL(gsl_matrix_complex_get(S, i, 2))); //sqrt(2)*Re(c3)
        gsl_matrix_set(R, i, 3 , sqrt(2)*GSL_IMAG(gsl_matrix_complex_get(S, i, 0))); //sqrt(2)*Im(c2)
        gsl_matrix_set(R, i, 4 , GSL_REAL(gsl_matrix_complex_get(S, i, 4)));         //c4
        gsl_matrix_set(R, i, 5 , sqrt(2)*GSL_IMAG(gsl_matrix_complex_get(S, i, 2))); //sqrt(2)*Im(c3)
    }

    //---------------------------------------------------------------
    //  Compute R^T J R
    //---------------------------------------------------------------
    //Transpose only
    gsl_matrix_transpose_memcpy(RH, R);  //transpose only
    //Matrix J
    glsc_matrix_set_J(J);
    //AUX = J*M
    gsl_blas_dgemm (CblasNoTrans , CblasNoTrans , 1.0, J , R , 0.0 , AUX );
    //BUX = MH*AUX = MH*J*M
    gsl_blas_dgemm (CblasNoTrans , CblasNoTrans , 1.0 , RH , AUX , 0.0 , BUX );

    //---------------------------------------------------------------
    //  Check the signs within BUX, which should be = J
    //  If the sign is wrong, execute a permutation in the couple (S, Dm)
    //---------------------------------------------------------------
    double sign;
    int flag = 0;
    for(int i = 0; i < 3; i ++)
    {
        sign = gsl_matrix_get(BUX, i, i +3);
        if(sign < 0)    //wrong sign, permutation!
        {
            evecA = gsl_matrix_complex_column (S, i);                //evecA = S(i)
            evecB = gsl_matrix_complex_column (S, i+3);              //evecB = S(i+3)

            evalA = gsl_matrix_complex_get(Dm, i, i);                //evalA = Dm(i,i)
            evalB = gsl_matrix_complex_get(Dm, i+3, i+3);            //evalB = Dm(i+3,i+3)

            gsl_vector_complex_memcpy(AUXc, &evecA.vector);          //AUXc = S(i)
            gsl_vector_complex_memcpy(&evecA.vector, &evecB.vector); //S(i) = S(i+3)
            gsl_vector_complex_memcpy(&evecB.vector, AUXc);          //S(i+3) = S(i)

            gsl_matrix_complex_set(Dm, i, i, evalB);                 //Dm(i,i) = Dm(i+3,i+3)
            gsl_matrix_complex_set(Dm, i+3, i+3, evalA);             //Dm(i+3,i+3) = Dm(i,i)
            flag = 1;
        }
    }

    //---------------------------------------------------------------
    //  Normalize S again, if necessary
    //---------------------------------------------------------------
    if(flag) norm_matrix_S(S);

    //---------------------------------------------------------------
    //  Update again R from S columns
    //---------------------------------------------------------------
    for(int i=0; i<6; i++)
    {
        gsl_matrix_set(R, i, 0 , sqrt(2)*GSL_REAL(gsl_matrix_complex_get(S, i, 0))); //sqrt(2)*Re(c2)
        gsl_matrix_set(R, i, 1 , GSL_REAL(gsl_matrix_complex_get(S, i, 1)));         //c1
        gsl_matrix_set(R, i, 2 , sqrt(2)*GSL_REAL(gsl_matrix_complex_get(S, i, 2))); //sqrt(2)*Re(c3)
        gsl_matrix_set(R, i, 3 , sqrt(2)*GSL_IMAG(gsl_matrix_complex_get(S, i, 0))); //sqrt(2)*Im(c2)
        gsl_matrix_set(R, i, 4 , GSL_REAL(gsl_matrix_complex_get(S, i, 4)));         //c4
        gsl_matrix_set(R, i, 5 , sqrt(2)*GSL_IMAG(gsl_matrix_complex_get(S, i, 2))); //sqrt(2)*Im(c3)
    }


    //---------------------------------------------------------------------------------
    //Visual check
    //---------------------------------------------------------------------------------
    cout << "---------------------------------------------" << endl;
    cout << "Symplectic test of S (after normalization):" << endl;
    cout << "---------------------------------------------" << endl;
    cout << "Please press Enter to continue." << endl;
    char chc;
    scanf("%c",&chc);
    symplecticity_test_complex(S, Csts::INVERSE_SYMP);
    cout << "Please press Enter to continue." << endl;
    scanf("%c",&chc);

    cout << "---------------------------------------------" << endl;
    cout << "Symplectic test of R (after normalization):" << endl;
    cout << "---------------------------------------------" << endl;
    cout << "Please press Enter to continue." << endl;
    scanf("%c",&chc);
    symplecticity_test_real(R, Csts::INVERSE_SYMP);
    cout << "Please press Enter to continue." << endl;
    scanf("%c",&chc);
}

//----------------------------------------------------------------------------------------
// Matrix decomposition
//----------------------------------------------------------------------------------------
/**
 *  \brief Obtain the monodromy matrix decomposition from a monodromy matrix given as a product of matrices, along with the matrix B of the Floquet theorem
 *   described in Dynamics and Mission Design near Libration Points IV, pp.67-90.
 *
 *  Steps of the computation:
 *  1. The monodromy matrix of the dynamical equivalent of EML1,2 is computed as a product of the State Transition Matrix (STM)
 *  computed on a M point grid in the variable **DAT: MMc = DAT[M]*DAT[M-1]*...*DAT[1].
 *  2. The resulting matrix is diagonalized with direct methods taken from GSL library.
 *  3. The hyperbolic directions are obtained via power methods (direct and inverse).
 *  4. A wielandt's deflation is applied on MMc to get rid of the unstable direction. The resulting matrix is diagonalized.
 *     - Dm contains the corresponding eigenvalues.
 *     - S  containes the corresponding eigenvectors in symplectic form.
 *  Note that a user-driven permutation within S and Dm may be done so that we ensure that:
 *
 *  S = | unstable-hyp. | center-xy | center-z | stable-hyp. | center-xy | center-z|
 *
 *  5. The routine monoDecomLog allows to build:
 *      - DB = log(Dm)/T
 *      - B = S*DB*Sinv
 *      - Br, the real form of B
 *      - JB, the Jordan form of B
 *      - the matrix R, so that Br = R*JB*Rinv
 *
 *   Note that the current code ensures that real(B) = Br and imag(B) = 0. In other words, the real form of B is already selected when computing B = S*DB*Sinv.
 */
void mono_eigen(gsl_odeiv2_driver *d,           //driver for odeRK78
                const double y0[],              //initial conditions
                const double n,                 //pulsation
                FBPL *fbpl,                     //Four-Body problem
                int M,                          //Number of steps for Monodromy matrix splitting in a product of matrices
                int STABLE_DIR_TYPE,            //Type of computation of the stable direction
                gsl_matrix_complex* MMc,        //Final monodromy matrix in complex form                                    |
                gsl_matrix* MM,                 //Final monodromy matrix in real form                                       |
                gsl_matrix_complex *Dm,         //The eigenvalues  are stored in Dm(i,i), i = 0,...,5                       |   Outputs
                gsl_matrix_complex *DB,         //DB = 1/T*log(Dm)                                                          |
                gsl_matrix_complex *B,          //B = S*DB*Sinv                                                             |
                gsl_matrix_complex *S,          //The eigenvectors are stored on  S(i,*), i = 0,...,5                       |
                gsl_matrix *JB,                 //Real Jordan form of B                                                     |
                gsl_matrix *Br,                 //Real form of B computed from the Jordan Form                              |
                gsl_matrix *R,                  //B = R*JB*Rinv                                                             |
                gsl_matrix_complex** DAT,       //Contains the STM at each of the M steps                                   |
                int is_stored)
{
    //------------------------------------------------------------------------------------
    // Initialization of GSL objects
    //------------------------------------------------------------------------------------
    //Hyperbolic directions
    gsl_vector_complex* eigenVu = gsl_vector_complex_calloc(6);  //Unstable eigenvector
    gsl_vector_complex* eigenVs = gsl_vector_complex_calloc(6);  //Stable eigenvector
    gsl_complex eigenLu, eigenLs;                               //(Un)stable eigenvalues

    //Misc complex numbers
    gsl_complex one_c  = gslc_complex(1.0, 0.0); //1.0

    //Auxiliary objects
    gsl_vector_complex *VEP  = gsl_vector_complex_calloc(6);     //auxiliary vector
    gsl_complex VUX, VAPL;
    int ch;

    //Constants
    double tend = 2*M_PI/n;  //period

    //------------------------------------------------------------------------------------
    // Initialization of outputs
    //------------------------------------------------------------------------------------
    //eigenvectors/eigenvalues of MMc obtained via Wieland deflation
    gsl_matrix_complex *evecr = gsl_matrix_complex_calloc(6, 6);
    gsl_vector_complex *evalr = gsl_vector_complex_calloc(6);
    gsl_matrix_complex *evmr  = gsl_matrix_complex_calloc(6, 6);
    //eigenvectors/eigenvalues of MMc obtained via direct diagonalization
    gsl_matrix_complex *evecd = gsl_matrix_complex_calloc(6, 6);
    gsl_matrix_complex *evmd  = gsl_matrix_complex_calloc(6, 6);
    gsl_vector_complex *evald = gsl_vector_complex_calloc(6);

    //------------------------------------------------------------------------------------
    // The Monodromy matrix obtained as a product of matrices in DAT.
    // CAREFUL: DAT has to be used from DAT[1] to DAT[M], and DAT[0] is USELESS.
    //------------------------------------------------------------------------------------
    cout << "mono_eigen. The monodromy matrix M is obtained as a product of matrices... " << endl;
    // DAT[1...M] is computed from the initial conditions y0
    // DAT[M] is computed at tend = orbital period of the system.
    // cout << "mono_eigen. Computation of DAT[1...M]... " << endl;
    stm_stepped_complex(y0, d, tend, M, DAT);

    // Alternative computation of DAT[1...M] in DAT2, making use of the symmetries of the orbit
    gsl_matrix_complex **DAT2 = gslc_matrix_complex_product_alloc(6, 6, M);
    gsl_matrix_complex **DDAT = gslc_matrix_complex_product_alloc(6, 6, M);
    // cout << "mono_eigen. Alternative computation of DAT[1...M]... " << endl;
    stm_stepped_complex_alt(y0, d, tend, M, DAT2);
    //DDAT = DAT-DAT2
    for(int k = 1; k <=M; k++)
    {
        gsl_matrix_complex_memcpy(DDAT[k], DAT2[k]);
        gsl_matrix_complex_sub(DDAT[k], DAT[k]);
    }

    //------------------------------------------------------
    // Comparison between DAT and DAT2 (uncomment to show)
    //------------------------------------------------------
    //    cout << Csts::SSEPR << endl;
    //    cout << "mono_eigen. The difference between the two computations is:" << endl;
    //    cout << setprecision(5);
    //    cout << gslc_matrix_complex_infinity_norm(DDAT[0])   << " at t = 0 "   << endl;
    //    cout << gslc_matrix_complex_infinity_norm(DDAT[M/2]) << " at t = T/2 "   << endl;
    //    cout << gslc_matrix_complex_infinity_norm(DDAT[M])   << " at t = T "   << endl;
    //    cout << setprecision(15);
    //    cout << Csts::SSEPR << endl;

    //DAT = DAT2 is performed by default, because numerical tests show that it usually
    //yields better precision. The interested user can comment the following line to see
    //the result with DAT.
    for(int k = 1; k <=M; k++) gsl_matrix_complex_memcpy(DAT[k], DAT2[k]);

    //Memory release of DAT2 and DDAT
    gslc_matrix_complex_product_free(DAT2, M);
    gslc_matrix_complex_product_free(DDAT, M);

    //------------------------------------------------------------------------------------
    // Monodromy matrix computed by product, stored in MMc. Real version stored in MM.
    //------------------------------------------------------------------------------------
    cout << "mono_eigen. The monodromy matrix M is computed by product... " << endl;
    mono_from_prod(MMc, DAT, M);
    //Real version stored in MM
    gslc_matrix_complex_real(MM, MMc);
    //Print on screen
    cout << "real(M) =  " << endl;
    gslc_matrix_complex_printf_real(MMc);
    cout << "--------------------------------" << endl;

    //------------------------------------------------------------------------------------
    // Diagonalization of MM: solving directly the system
    //------------------------------------------------------------------------------------
    cout << "mono_eigen. Direct diagonalization of M... " << endl;
    mono_diag(MM, evecd, evald, evmd);
    cout << "mono_eigen. Test of the eigensystem obtained via direct diagonalization:" << endl;
    eigensystem_test(evmd, evecd, DAT, M);
    //gslc_eigensystem_fprintf(evald, evecd, 6, (char*) "diagDirect.txt");
    cout << "--------------------------------" << endl;

    //------------------------------------------------------------------------------------
    //Shifted power method (direct on MMc, bad precision, uncomment to show)
    //------------------------------------------------------------------------------------
    //gsl_vector_complex* eigenV = gsl_vector_complex_calloc(6);
    //gsl_complex eigenL;
    //shifted_power_method(MMc, 1e-15, 6, gsl_vector_complex_get(DECSl, 0), eigenV, &eigenL);

    //====================================================================================
    // Use of Power+Wielandt?
    // For very unstable cases (eg EML1,2 in the QBCP), use Power+Wielandt (the user enter 1).
    // For other cases, enter 0 (eg SEML1,2).
    //====================================================================================
    printf("Please enter:\n 0 to use current eigendecomposition obtained via direct diagonalization \n 1 if you want the Power+Wielandt process to start\n");
    printf("(by default, use the Power+Wielandt process on the Earth-Moon system, and not for the Sun-Earth case):");
    scanf("%d",&ch);
    if(ch == 1)
    {
        //--------------------------------------------------------------------------------
        // 3. Hyp. unstable direction with product of matrices
        //--------------------------------------------------------------------------------
        cout << "mono_eigen. Power methods on hyperbolic directions." << endl;
        //Direct power method
        gsl_vector_complex_set_zero(VEP);
        gsl_vector_complex_set(VEP, 0, one_c);
        power_method_on_prod(DAT+1, M, 6, VEP, &VUX, &VAPL, 1);
        //Update of eigenVu and eigenLu
        gsl_vector_complex_memcpy(eigenVu, VEP);
        eigenLu = VUX;

        //--------------------------------------------------------------------------------
        // 4. Hyp. stable direction with product of matrices
        //--------------------------------------------------------------------------------
        if(STABLE_DIR_TYPE == Csts::STABLE_DIR_POW)
        {
            //Has to be used in priority: indeed, the computation using the symmetries
            //(through Csts::STABLE_DIR_SYM) does not take into account the numerical instabilities
            //Inverse power method
            gsl_vector_complex_set_zero(VEP);
            gsl_vector_complex_set(VEP, 1, one_c);
            power_method_on_prod(DAT+1, M, 6, VEP, &VUX, &VAPL, -1);

            //Update of eigenVs and eigenLs
            gsl_vector_complex_memcpy(eigenVs, VEP);
            eigenLs = gsl_complex_inverse(VUX);
        }
        else if(STABLE_DIR_TYPE == Csts::STABLE_DIR_SYM)
        {
            //Update of eigenVs and eigenLs
            gslc_vector_complex_symcpy(eigenVs, eigenVu);
            eigenLs = gsl_complex_inverse(eigenLu);
        }
        else
        {
            //stable direction with inverse power method directly on the mon. matrix
            //GENERALLY GIVES BAD RESULTS, JUST HERE FOR REFERENCE
            inverse_power_method(MMc, 1e-15, 6, 0, eigenVs, &eigenLs);
        }

        //--------------------------------------------------------------------------------
        // 5. Diagonalization of MMc via Wielandt deflation on the unstable direction
        // After this step:
        // - evecr contains the center eigenvectors of MMc.
        // - evalr contains the center eigenvalues  of MMc.
        // - eigenVu/Lu contains the unstable couple of MMc.
        // - eigenVs/Ls contains the stable couple of MMc.
        //--------------------------------------------------------------------------------
        if(gsl_complex_abs(eigenLu) > 1e3)
        {
            cout << "mono_eigen. Diagonalization of the mono. matrix via Wielandt deflation on the unstable direction." << endl;
            wielandt_deflation(MMc, eigenVu, eigenLu, evecr, evalr);
            //Uncomment to print in root folder
            //gslc_eigensystem_fprintf(evalr, evecr, 6, (char*) "diagWielandt.txt");
        }
        else //probably never used, needs a proper check if the case occurs.
        {
            cout << "mono_eigen. The unstable eigenvalue is small enough to avoid the use of Wielandt deflation." << endl;
            gsl_vector_complex_view viewd, viewl;

            viewd = gsl_matrix_complex_column(evecd, 0); //vs
            viewl = gsl_matrix_complex_column(evecr, 1);
            gsl_vector_complex_memcpy(&viewl.vector, &viewd.vector);
            gsl_vector_complex_set(evalr, 1, gsl_vector_complex_get(evald, 0));

            gsl_vector_complex_memcpy(eigenVs, &viewd.vector);
            eigenLs =  gsl_vector_complex_get(evald, 0);

            viewd = gsl_matrix_complex_column(evecd, 1); //vxy
            viewl = gsl_matrix_complex_column(evecr, 4);
            gsl_vector_complex_memcpy(&viewl.vector, &viewd.vector);
            gsl_vector_complex_set(evalr, 4, gsl_vector_complex_get(evald, 1));

            viewd = gsl_matrix_complex_column(evecd, 2); //vxy
            viewl = gsl_matrix_complex_column(evecr, 5);
            gsl_vector_complex_memcpy(&viewl.vector, &viewd.vector);
            gsl_vector_complex_set(evalr, 5, gsl_vector_complex_get(evald, 2));

            viewd = gsl_matrix_complex_column(evecd, 3); //vz
            viewl = gsl_matrix_complex_column(evecr, 2);
            gsl_vector_complex_memcpy(&viewl.vector, &viewd.vector);
            gsl_vector_complex_set(evalr, 2, gsl_vector_complex_get(evald, 3));

            viewd = gsl_matrix_complex_column(evecd, 4); //vz
            viewl = gsl_matrix_complex_column(evecr, 3);
            gsl_vector_complex_memcpy(&viewl.vector, &viewd.vector);
            gsl_vector_complex_set(evalr, 3, gsl_vector_complex_get(evald, 4));

            viewd = gsl_matrix_complex_column(evecd, 5); //vu
            viewl = gsl_matrix_complex_column(evecr, 0);
            gsl_vector_complex_memcpy(&viewl.vector, &viewd.vector);
            gsl_vector_complex_set(evalr, 0, gsl_vector_complex_get(evald, 5));

            gsl_vector_complex_memcpy(eigenVu, &viewd.vector);
            eigenLu =  gsl_vector_complex_get(evald, 5);
        }

    }
    else
    {
        //--------------------------------------------------------------------------------
        // Steps 3/4/5 in one step
        // After this step:
        // - evecr contains the center eigenvectors of MMc.
        // - evalr contains the center eigenvalues  of MMc.
        // - eigenVu/Lu contains the unstable couple of MMc.
        // - eigenVs/Ls contains the stable couple of MMc.
        //--------------------------------------------------------------------------------
        //Heuristic to ensure that the matrix R will be symplectic (see mono_eigen_log)
        int kmp[6];
        switch(fbpl->li)
        {
        case 1:
            kmp[0] = 0;
            kmp[1] = 2;
            kmp[2] = 1;
            kmp[3] = 3;
            kmp[4] = 4;
            kmp[5] = 5;
            cout << "L1 case in heuristic." << endl;
            break;
        case 2:
            kmp[0] = 0;
            kmp[1] = 2;
            kmp[2] = 1;
            kmp[3] = 4;
            kmp[4] = 3;
            kmp[5] = 5;
            cout << "L2 case in heuristic." << endl;
            break;
        default: //never here, L2 case by default
            kmp[0] = 0;
            kmp[1] = 2;
            kmp[2] = 1;
            kmp[3] = 4;
            kmp[4] = 3;
            kmp[5] = 5;
        }

        //--------------------------------------------------------------------------------
        // Copy of the direct eigendecomposition into evecr and evalr,
        // After this step:
        // - evecr contains the center eigenvectors of MMc.
        // - evalr contains the center eigenvalues  of MMc.
        // - eigenVu/Lu contains the unstable couple of MMc.
        // - eigenVs/Ls contains the stable couple of MMc.
        // (As in the Power+Wielandt case).
        //--------------------------------------------------------------------------------
        gsl_vector_complex_view viewd, viewl;
        viewd = gsl_matrix_complex_column(evecd, kmp[0]); //vs
        viewl = gsl_matrix_complex_column(evecr, 1);
        gsl_vector_complex_memcpy(&viewl.vector, &viewd.vector);
        gsl_vector_complex_set(evalr, 1, gsl_vector_complex_get(evald, kmp[0]));

        gsl_vector_complex_memcpy(eigenVs, &viewd.vector);
        eigenLs =  gsl_vector_complex_get(evald, kmp[0]);

        viewd = gsl_matrix_complex_column(evecd, kmp[1]); //vxy
        viewl = gsl_matrix_complex_column(evecr, 4);
        gsl_vector_complex_memcpy(&viewl.vector, &viewd.vector);
        gsl_vector_complex_set(evalr, 4, gsl_vector_complex_get(evald, kmp[1]));

        viewd = gsl_matrix_complex_column(evecd, kmp[2]); //vxy
        viewl = gsl_matrix_complex_column(evecr, 5);
        gsl_vector_complex_memcpy(&viewl.vector, &viewd.vector);
        gsl_vector_complex_set(evalr, 5, gsl_vector_complex_get(evald, kmp[2]));

        viewd = gsl_matrix_complex_column(evecd, kmp[3]); //vz
        viewl = gsl_matrix_complex_column(evecr, 2);
        gsl_vector_complex_memcpy(&viewl.vector, &viewd.vector);
        gsl_vector_complex_set(evalr, 2, gsl_vector_complex_get(evald, kmp[3]));

        viewd = gsl_matrix_complex_column(evecd, kmp[4]); //vz
        viewl = gsl_matrix_complex_column(evecr, 3);
        gsl_vector_complex_memcpy(&viewl.vector, &viewd.vector);
        gsl_vector_complex_set(evalr, 3, gsl_vector_complex_get(evald, kmp[4]));

        viewd = gsl_matrix_complex_column(evecd, kmp[5]); //vu
        viewl = gsl_matrix_complex_column(evecr, 0);
        gsl_vector_complex_memcpy(&viewl.vector, &viewd.vector);
        gsl_vector_complex_set(evalr, 0, gsl_vector_complex_get(evald, kmp[5]));

        gsl_vector_complex_memcpy(eigenVu, &viewd.vector);
        eigenLu =  gsl_vector_complex_get(evald, kmp[5]);
    }

    //Set evald in matrix format in evmr
    gsl_matrix_complex_set_zero(evmr);
    for(int i = 0; i <6; i++) gsl_matrix_complex_set(evmr, i, i, gsl_vector_complex_get(evalr, i));

    //------------------------------------------------------------------------------------
    // Force some coefficients to zero: uncomment to set on
    // Help decoupling the vertical and xy motion
    //------------------------------------------------------------------------------------
    for(int i = 0; i <6; i++)
    {
        for(int j = 0; j <6; j++)
        {
            if(gsl_complex_abs(gsl_matrix_complex_get(evecr, i, j)) < 1e-14) gsl_matrix_complex_set(evecr, i, j, gslc_complex(0.0, 0.0));
            if(gsl_complex_abs(gsl_vector_complex_get(eigenVu, j))  < 1e-14) gsl_vector_complex_set(eigenVu, j,  gslc_complex(0.0, 0.0));
            if(gsl_complex_abs(gsl_vector_complex_get(eigenVs, j))  < 1e-14) gsl_vector_complex_set(eigenVs, j,  gslc_complex(0.0, 0.0));
        }
    }

    //------------------------------------------------------------------------------------
    // At this step:
    // - evecr contains the center eigenvectors of MMc, at positions (1,2,3,4)
    // - evalr contains the center eigenvalues  of MMc, at positions (1,2,3,4)
    // - eigenVu/Lu contains the unstable couple of MMc.
    // - eigenVs/Ls contains the stable couple of MMc.
    //------------------------------------------------------------------------------------


    //====================================================================================
    // 6. Computation of stable, unstable and central directions through vepro
    // After this step:
    // - VECS(*, 2:5, 0) contains a 6*4 central subspace that contains no unstable
    //   component under DAT[M]*...*DAT[1]
    // - VECS(*, 2:5, M) contains the same central subspace as if it has been
    //   integrated through the variationnal equations until t = T
    // - DECSv = eigenvectors of the matrix DECS such that VECS[0]*DECS = VECS[M] in vepro
    // - DECSl = eigenvalues of the matrix DECS such that  VECS[0]*DECS = VECS[M] in vepro
    //
    // This part makes use of vepro, inspired from Fortran code by Masdemont.
    //
    //
    // The results are NOT conclusive, and the central and hyperbolic directions are kept
    // as is.
    //====================================================================================
    /*
     cout << "mono_eigen. Computation of stable, unstable and central directions through vepro." << endl;
     //VECS(*,*,0) contains a 4*4 central subspace that contains no unstable component under DAT[M]*...*DAT[1].
     gsl_matrix_complex **VECS = gslc_matrix_complex_product_alloc(6, 6, M);
     eigenvectors/eigenvalues of the matrix DECS such that VECS[0]*DECS = VECS[M] in vepro
     gsl_matrix_complex *DECSv = gsl_matrix_complex_calloc(6, 6);
     gsl_matrix_complex *DECS  = gsl_matrix_complex_calloc(6, 6);
     gsl_vector_complex *DECSl = gsl_vector_complex_calloc(6);
     vepro(DAT, M, VECS, DECS, DECSv, DECSl);
     */

    //====================================================================================
    // Annex: Solving the system by hand. NOT conclusive either
    //====================================================================================
    /*
    gsl_vector_complex* y  = gsl_vector_complex_calloc(6);
    cout << "Multimin for lambda = " << GSL_REAL(gsl_vector_complex_get(DECSl, 0)) << " " << GSL_IMAG(gsl_vector_complex_get(DECSl, 0)) << endl;
    cout << "The eigenvector approximation is taken from Wielandt's deflation procedure:" << endl;
    evec_i = gsl_matrix_complex_column (evecr, 5);
    gsl_vector_complex_set(&evec_i.vector, 2, zero_c);
    gsl_vector_complex_set(&evec_i.vector, 5, zero_c);
    gslc_vector_complex_printf(&evec_i.vector);


    //mm = MMc - lamba id
    cdouble **mm = dcmatrix(0, 5, 0, 5);
    cdouble *v1 = dcvector(0,5);
    cdouble lambda = GSL_REAL(gsl_vector_complex_get(DECSl, 0)) + GSL_IMAG(gsl_vector_complex_get(DECSl, 0))*I;
    for(int i = 0; i<6; i++)
    {
        for(int j = 0; j<6; j++)
        {
            mm[i][j] = GSL_REAL(gsl_matrix_complex_get(MMc, i, j)) + GSL_IMAG(gsl_matrix_complex_get(MMc, i, j))*I;
            if(i==j) mm[i][j] -= lambda;
        }
    }

    v1[0] = GSL_REAL(gsl_vector_complex_get(&evec_i.vector, 0)) + GSL_IMAG(gsl_vector_complex_get(&evec_i.vector, 0))*I;
    v1[1] = GSL_REAL(gsl_vector_complex_get(&evec_i.vector, 1)) + GSL_IMAG(gsl_vector_complex_get(&evec_i.vector, 1))*I;
    v1[2] = 0.0;
    v1[3] = GSL_REAL(gsl_vector_complex_get(&evec_i.vector, 3)) + GSL_IMAG(gsl_vector_complex_get(&evec_i.vector, 3))*I;
    v1[4] = GSL_REAL(gsl_vector_complex_get(&evec_i.vector, 4)) + GSL_IMAG(gsl_vector_complex_get(&evec_i.vector, 4))*I;
    v1[5] = 0.0;

    cout << "v1 = " << endl;
    for(int i = 0; i<6; i++)  cout << creal(v1[i]) << "  " << cimag(v1[i]) << endl;
    cout << "M*v1?" << endl;
    for(int i = 0; i < 6; i ++)
    {
        cdouble res = mm[i][0]*v1[0] + mm[i][1]*v1[1] + mm[i][2]*v1[2] + mm[i][3]*v1[3] + mm[i][4]*v1[4] + mm[i][5]*v1[5];
        cout << creal(res) << "  " << cimag(res) << endl;
    }

    cdouble x1     = 1.0+1.0*I;
    cdouble alpha1 = x1*(mm[1][0] - mm[1][1]*mm[0][0]/mm[0][1]);
    cdouble alpha4 = mm[1][3] - mm[0][3]*mm[1][1]/mm[0][1];
    cdouble alpha5 = mm[1][4] - mm[0][4]*mm[1][1]/mm[0][1];

    cdouble eps1   = x1*(mm[3][0] - mm[3][1]*mm[0][0]/mm[0][1]);
    cdouble eps4   = mm[3][3] - mm[3][1]*mm[0][3]/mm[0][1];
    cdouble eps5   = mm[3][4] - mm[3][1]*mm[0][4]/mm[0][1];

    cdouble det    = alpha4*eps5 - eps4*alpha5;
    cdouble **a    = dcmatrix(0, 2, 0, 2);

    a[0][0] = +eps5/det;
    a[0][1] = -alpha5/det;
    a[1][0] = -eps4/det;
    a[1][1] = +alpha4/det;

    cdouble x4 = -alpha1*a[0][0] - eps1*a[0][1];
    cdouble x5 = -alpha1*a[1][0] - eps1*a[1][1];
    cdouble x2 = -1.0/mm[0][1]*(mm[0][0]*x1 + mm[0][3]*x4 + mm[0][4]*x5);

    gsl_vector_complex_set(y, 0, gslc_complex(creal(x1), cimag(x1)));
    gsl_vector_complex_set(y, 1, gslc_complex(creal(x2), cimag(x2)));
    gsl_vector_complex_set(y, 3, gslc_complex(creal(x4), cimag(x4)));
    gsl_vector_complex_set(y, 4, gslc_complex(creal(x5), cimag(x5)));
    gsl_vector_complex_scale(y, gslc_complex(-1.807784977191888e-01, +3.006437368959317e-02));
    //gslc_vector_complex_normalize(y);

    v1[0] = GSL_REAL(gsl_vector_complex_get(y, 0)) + GSL_IMAG(gsl_vector_complex_get(y, 0))*I;
    v1[1] = GSL_REAL(gsl_vector_complex_get(y, 1)) + GSL_IMAG(gsl_vector_complex_get(y, 1))*I;
    v1[3] = GSL_REAL(gsl_vector_complex_get(y, 3)) + GSL_IMAG(gsl_vector_complex_get(y, 3))*I;
    v1[4] = GSL_REAL(gsl_vector_complex_get(y, 4)) + GSL_IMAG(gsl_vector_complex_get(y, 4))*I;

    cout << "final vector = " << endl;
    gslc_vector_complex_printf(y);



    cout << "M*v1?" << endl;
    for(int i = 0; i < 6; i ++)
    {
        long cdouble res = mm[i][0]*v1[0] + mm[i][1]*v1[1] + mm[i][2]*v1[2] + mm[i][3]*v1[3] + mm[i][4]*v1[4] + mm[i][5]*v1[5];
        cout << creal(res) << "  " << cimag(res) << endl;
    }

    cout << "-------------------------------------------------" << endl;
    cout << "-------------------------------------------------" << endl;
    cout << "-------------------------------------------------" << endl;

    //*/

    //----------------------------------------------------------------------------------
    // 7. Update S and Dm, with the form (center1 x saddle x center2),
    // with center2 corresponding to the vertical components (z, pz)
    // WARNING: the center 2 should really correspond to the vertical components (z, pz)
    // Decoupled from the rest of the dynamics.
    // As a consequence, the user is asked to manually check the former condition.
    // If it is not fulfilled, a permutation is made.
    //----------------------------------------------------------------------------------
    cout << "mono_eigen. Update the matrices S and Dm." << endl;
    cout << "After this step:                         " << endl;
    cout << "- S contains the eigenvectors in the following order:" << endl;
    cout << "      (vxy, vu, vz, conj(vxy), vs, conj(vz))" << endl;
    cout << "- Dm contains the corresponding eigenvalues." << endl;
    update_matrix_S(S, Dm, evecr, evalr, eigenVu, eigenVs, eigenLu, eigenLs);

    //----------------------------------------------------------------------------------
    // 6.2 Prenormalization (if necessary) of S. Uncomment to let the user decide.
    //----------------------------------------------------------------------------------
    prenorm_matrix_S(S);

    //----------------------------------------------------------------------------------
    // 7. Normalization of S to get a symplectic matrix
    //----------------------------------------------------------------------------------
    norm_matrix_S(S);

    //----------------------------------------------------------------------------------
    // 8. Second permutation within S to garantee that R is truly symplectic.
    //----------------------------------------------------------------------------------
    symplecticity_test_for_R(S, Dm, R);

    //----------------------------------------------------------------------------------
    // 9. Obtention of B, Br, DB and JB
    //----------------------------------------------------------------------------------
    //Obtain the real Jordan form of the matrix B
    mono_eigen_log(S, Dm, tend, DB, B, Br, JB, R);


    //----------------------------------------------------------------------------------
    // 10. Print in txt file
    //----------------------------------------------------------------------------------
    if(is_stored)
    {
        string F_COC = fbpl->cs.F_COC;
        string filename = F_COC+"complexCOC.txt";
        cout << "mono_eigen. For eye-checking: storage of the complex change of coord. in " << filename << endl;
        gslc_eigensystem_fprintf(Dm, S, 6, (char*) filename.c_str());

        filename = F_COC+"B.txt";
        cout << "mono_eigen. Storage of the complex matrix B  in  " << filename << endl;
        gslc_matrix_complex_fprintf(B, (char*) filename.c_str());

        filename = F_COC+"DB.txt";
        cout << "mono_eigen. Storage of the complex matrix DB in  " << filename << endl;
        gslc_matrix_complex_fprintf(DB, (char*) filename.c_str());

        filename = F_COC+"Dm.txt";
        cout << "mono_eigen. Storage of the complex matrix Dm in  " << filename << endl;
        gslc_matrix_complex_fprintf(Dm, (char*) filename.c_str());

        filename = F_COC+"S.txt";
        cout << "mono_eigen. Storage of the complex matrix S  in  " << filename << endl;
        gslc_matrix_complex_fprintf(S, (char*) filename.c_str());

        filename = F_COC+"Br.txt";
        cout << "mono_eigen. Storage of the real matrix    Br in  " << filename << endl;
        gslc_matrix_fprintf(Br, (char*) filename.c_str());

        filename = F_COC+"JB.txt";
        cout << "mono_eigen. Storage of the real matrix    JB in  " << filename << endl;
        gslc_matrix_fprintf(JB, (char*) filename.c_str());
    }


    //----------------------------------------------------------------------------------
    // 11. Tests
    //----------------------------------------------------------------------------------
    //Test of the eigenvectors
    //----------------------------------------------------------------------------------
    cout << "Test of Wielandt eigensystem:" << endl;
    eigensystem_test(Dm, S, DAT, M);

    //Test of the change of base - uncomment if needed
    //----------------------------------------------------------------------------------
    /*
      gsl_matrix_complex *MMc_estimate  = gsl_matrix_complex_calloc(6, 6);
      change_of_base(MMc_estimate, Dm, DAT, S, Csts::INVERSE_GSL, M);

      cout << "---------------------------------------------" << endl;
      cout << "M" << endl;
      cout << "---------------------------------------------" << endl;
      gslc_matrix_complex_printf_real(MMc);

      cout << "---------------------------------------------" << endl;
      cout << "M - M estimate (real part)"  << endl;
      cout << "---------------------------------------------" << endl;
      gsl_matrix_complex_sub(MMc_estimate, MMc);
      gslc_matrix_complex_printf_real(MMc_estimate);
      gsl_matrix_complex_free(MMc_estimate);
    //*/

    //----------------------------------------------------------------------------------
    // 12. Memory release
    //----------------------------------------------------------------------------------
    gsl_matrix_complex_free(evecr);
    gsl_vector_complex_free(VEP);
    gsl_vector_complex_free(eigenVu);
    gsl_vector_complex_free(eigenVs);
    gsl_vector_complex_free(evalr);
}


/**
 *  \brief Compute the following matrices (see nfo2 brief):
 *          - DB = 1/T*log(Dm)
 *          - B = S*DB*Sinv
 *          - the real Jordan form JB of the matrix B
 *          - Br = R*JB*Rinv
 *  With the current implementation, real(B) = Br and imag(B) = 0
 */
void mono_eigen_log(gsl_matrix_complex *S,       //The eigenvectors are stored on  S(*,i), i = 0,...,5
                   gsl_matrix_complex *Dm,      //The eigenvalues  are stored in Dm(i,i), i = 0,...,5
                   double T,                    //Full period
                   //double n,                    //Pulsation
                   gsl_matrix_complex *DB,      //DB = 1/T*log(Dm)
                   gsl_matrix_complex *B,       //B = S*DB*Sinv
                   gsl_matrix *Br,              //Real form of B computed from the Jordan Form
                   gsl_matrix *JB,              //Real Jordan form of B
                   gsl_matrix *R)               //B = R*JB*Rinv
{
    //Misc complex numbers
    gsl_complex one_c  = gslc_complex(1.0, 0.0);
    gsl_complex zero_c = gslc_complex(0.0, 0.0);
    //Complex matrices
    gsl_matrix_complex *Sinv = gsl_matrix_complex_calloc(6, 6);
    gsl_matrix_complex *AUX  = gsl_matrix_complex_calloc(6, 6);
    //Real matrices
    gsl_matrix *Rc   = gsl_matrix_calloc(6,6);
    gsl_matrix *Rinv = gsl_matrix_calloc(6,6);
    gsl_matrix *AUXr = gsl_matrix_calloc(6,6);
    //For GSL inversion
    gsl_permutation * p = gsl_permutation_alloc (6);

    //----------------------------------------------------------------------------------
    // DB = 1/T*log(Dm)
    //----------------------------------------------------------------------------------
    gsl_complex logL;
    for(int i =0; i < 6; i++)
    {
        logL = gsl_complex_log(gsl_matrix_complex_get(Dm, i, i));
        logL = gsl_complex_div_real(logL, T);

        //Forcing some coefficients to zero
        if(fabs(GSL_REAL(logL)) < 1e-10) logL = gslc_complex(0.0, GSL_IMAG(logL));
        if(fabs(GSL_IMAG(logL)) < 1e-10) logL = gslc_complex(GSL_REAL(logL), 0.0);

        gsl_matrix_complex_set(DB, i, i, logL);
    }
    cout << "mono_eigen_log. Small real and imag parts of coefs of DB have been set null." << endl;

    //----------------------------------------------------------------------------------
    // B = S*DB*Sinv
    //----------------------------------------------------------------------------------
    gslc_matrix_complex_symplectic_inverse(S, Sinv);
    //AUX = DB*Sinv
    gsl_blas_zgemm (CblasNoTrans , CblasNoTrans , one_c , DB , Sinv , zero_c , AUX );
    //B = S*AUX = S*DB*Sinv
    gsl_blas_zgemm (CblasNoTrans , CblasNoTrans , one_c , S , AUX , zero_c , B);

    //====================================================================================
    //Real normal form
    //====================================================================================
    double omega1 = GSL_IMAG(gsl_matrix_complex_get(DB, 0, 0));
    double omega2 = GSL_REAL(gsl_matrix_complex_get(DB, 1, 1));
    double omega3 = GSL_IMAG(gsl_matrix_complex_get(DB, 2, 2));

    cout << "---------------------------------------------" << endl;
    cout << "Frequencies in B:" << endl;
    cout << "---------------------------------------------" << endl;
    cout << "omega1 = " << omega1 << endl;
    cout << "omega2 = " << omega2 << endl;
    cout << "omega3 = " << omega3 << endl;

    gsl_matrix_set_zero(JB);
    gsl_matrix_set(JB, 1, 1,  omega2);
    gsl_matrix_set(JB, 4, 4, -omega2);

    gsl_matrix_set(JB, 3, 0, -omega1);
    gsl_matrix_set(JB, 0, 3,  omega1);

    gsl_matrix_set(JB, 5, 2, -omega3);
    gsl_matrix_set(JB, 2, 5,  omega3);

    //cout << "---------------------------------------------" << endl;
    //cout << "Real normal form of B:" << endl;
    //cout << "---------------------------------------------" << endl;
    //gslc_matrix_printf(JB);


    //====================================================================================
    //B = R*JB*Rinv, with R = [sqrt(2)*Re(c2), c1, sqrt(2)*Re(c3), sqrt(2)*Im(c2), c4, sqrt(2)*Im(c3)]
    //====================================================================================
    //Rinv = R^{-1} with Symplectic routine
    gslc_matrix_symplectic_inverse(R, Rinv);

    //Rinv = R^{-1} with GSL routines
    //    int s;
    //    gsl_matrix_memcpy(Rc, R);
    //    gsl_linalg_LU_decomp (Rc, p , &s );
    //    gsl_linalg_LU_invert (Rc, p , Rinv );

    //AUX = JB*Rinv
    gsl_blas_dgemm (CblasNoTrans , CblasNoTrans , 1.0 , JB , Rinv , 0.0 , AUXr );
    //Br = R*AUX = R*JB*Rinv
    gsl_blas_dgemm (CblasNoTrans , CblasNoTrans , 1.0 , R , AUXr , 0.0 , Br);

    cout << "---------------------------------------------" << endl;
    cout << "B (real):" << endl;
    cout << "---------------------------------------------" << endl;
    gslc_matrix_complex_printf_real(B);

    cout << "---------------------------------------------" << endl;
    cout << "B (imag):" << endl;
    cout << "---------------------------------------------" << endl;
    gslc_matrix_complex_printf_imag(B);

    cout << "---------------------------------------------" << endl;
    cout << "Br:" << endl;
    cout << "---------------------------------------------" << endl;
    gslc_matrix_printf(Br);


    //Memory release
    gsl_matrix_complex_free(Sinv);
    gsl_matrix_complex_free(AUX);
    gsl_matrix_free(Rc);
    gsl_matrix_free(Rinv);
    gsl_matrix_free(AUXr);
    gsl_permutation_free(p);

}

//----------------------------------------------------------------------------------------
//Resolution of eigensystem: power methods and gsl methods on single matrices
//----------------------------------------------------------------------------------------
/**
 *  \brief Power Method on a single matrix.
 */
void powerMethod(gsl_matrix_complex *Minit, double prec, int sizeM, int version, gsl_vector_complex *eigenV, gsl_complex *eigenL)
{
    //====================================================================================
    // Misc complex numbers
    //====================================================================================
    gsl_complex one_c = gslc_complex(1.0, 0.0);
    gsl_complex zero_c = gslc_complex(0.0, 0.0);

    //====================================================================================
    // Power method
    //====================================================================================
    cout << "-----------------------------------------" << endl;
    cout << "Power method" << endl;
    cout << "-----------------------------------------" << endl;
    gsl_matrix_complex* M = gsl_matrix_complex_calloc (sizeM, sizeM);
    gsl_matrix_complex_memcpy(M, Minit);

    //Initial guess of the eigenvector of the unstable direction
    gsl_vector_complex *xm = gsl_vector_complex_calloc(sizeM);
    gsl_complex ymnormc;


    if(version == 0) //Euclidian norm version, needs to normalize
    {
        //for(int i=0; i<sizeM; i++) gsl_vector_complex_set(xm, i, one_c);
        gsl_vector_complex_set(xm, 0, one_c);
        for(int i=1; i<sizeM; i++) gsl_vector_complex_set(xm, i, zero_c);
        GSL_SET_COMPLEX(&ymnormc, 1.0/gsl_blas_dznrm2(xm), 0.0);
        gsl_vector_complex_scale (xm, ymnormc);
    }
    else
    {
        gsl_vector_complex_set(xm, 0, one_c);
        for(int i=1; i<sizeM; i++) gsl_vector_complex_set(xm, i, zero_c);
    }

    //Loop
    gsl_vector_complex *ym = gsl_vector_complex_calloc(sizeM);
    gsl_vector_complex *dxm = gsl_vector_complex_calloc(sizeM);
    gsl_complex lm;
    double dx;
    gsl_complex lm_prec = gslc_complex(0.0, 0.0);
    double dl;
    int iter = 0;
    if(version == 0)
    {
        // Euclidian norm version
        //-------------------------------------
        do
        {
            //dx(m) = x(m-1);
            gsl_vector_complex_memcpy(dxm, xm);

            //y(m) = A*x(m-1)
            gsl_blas_zgemv (CblasNoTrans, one_c , M , xm , zero_c , ym);

            //l(m) = x(m-1)^H*y(m)
            gsl_blas_zdotc(xm , ym ,&lm);

            //norm(ym) = sqrt(y(m)^H*y(m))
            GSL_SET_COMPLEX(&ymnormc, 1.0/gsl_blas_dznrm2 (ym), 0.0);


            //----------
            //Error
            //----------
            //dxm = -ym
            gsl_vector_complex_memcpy(dxm, ym);
            for(int i=0; i<sizeM; i++) gsl_vector_complex_set(dxm, i, gsl_complex_negative(gsl_vector_complex_get(dxm, i)));
            //dxm+= lm*xm
            gsl_blas_zaxpy (lm , xm , dxm);
            //dxm/=ymnormc
            gsl_vector_complex_scale (dxm , ymnormc);

            //dx = norm(dx(m))
            dx = gsl_blas_dznrm2 (dxm);
            //----------

            //x(m) = y(m)/norm(y(m))
            gsl_vector_complex_memcpy(xm, ym);
            gsl_vector_complex_scale (xm , ymnormc);


            //cout << iter << " " << dx << " " << GSL_REAL(lm) << " + " << GSL_IMAG(lm) << "i " << endl;

            if(iter ==0)
            {
                dl = 1;
                lm_prec = lm;
            }
            else
            {
                dl = gsl_complex_abs(gsl_complex_sub(lm,lm_prec))/gsl_complex_abs(lm);
                lm_prec = lm;
            }

            iter++;

        }
        while((dx > prec || dl > prec) && iter < 50);
    }
    else
    {
        //-------------------------------------
        // Infinity norm version
        //-------------------------------------

        gsl_vector *ymabs = gsl_vector_calloc(sizeM);
        int maxIndex = 0;
        do
        {

            //dx(m) = x(m-1);
            gsl_vector_complex_memcpy(dxm, xm);

            //y(m) = A*x(m-1)
            gsl_blas_zgemv (CblasNoTrans, one_c , M , xm , zero_c , ym);

            //norm_inf(ym) = max(|y(m)|)
            for(int i = 0; i<sizeM; i++) gsl_vector_set(ymabs, i, fabs(gsl_complex_abs(gsl_vector_complex_get(ym, i))));
            maxIndex = gsl_vector_max_index (ymabs);

            //l(m) = x(m-1)^H*y(m)
            lm = gsl_vector_complex_get(ym, maxIndex);


            //x(m) = y(m)/norm(y(m))
            ymnormc =  gsl_complex_inverse(gsl_vector_complex_get(ym, maxIndex));
            gsl_vector_complex_memcpy(xm, ym);
            gsl_vector_complex_scale (xm , ymnormc);

            //dx(m) = dx(m-1) - dx(m)
            gsl_vector_complex_sub(dxm, xm);

            //dx = norm(dx(m))
            dx = gsl_blas_dznrm2 (dxm);

            cout << iter << " " << dx << " " << GSL_REAL(lm) << " + " << GSL_IMAG(lm) << "i " << endl;

            iter++;

        }
        while(dx > prec && iter < 50);

        gsl_vector_free(ymabs);
    }
    //-------------------------------------
    cout << "Eigenvalue: "  <<  GSL_REAL(lm) << " + " << GSL_IMAG(lm) << "i "  << endl;
    cout << "Eigenvector:" << endl;
    gsl_vector_complex_fprintf(stdout, xm, "%+5.5e");


    //Outputs
    gsl_vector_complex_memcpy(eigenV, xm);
    *eigenL =  lm;

    gsl_matrix_complex_free(M);
    gsl_vector_complex_free(xm);
    gsl_vector_complex_free(ym);
    gsl_vector_complex_free(dxm);
}

/**
 *  \brief Inverse Power Method on a single matrix.
 */
void inverse_power_method(gsl_matrix_complex *Minit, double prec, int sizeM, int version, gsl_vector_complex *eigenV, gsl_complex *eigenL)
{

    //====================================================================================
    // Misc complex numbers
    //====================================================================================
    gsl_complex one_c = gslc_complex(1.0, 0.0);
    gsl_complex zero_c = gslc_complex(0.0, 0.0);

    //====================================================================================
    // Inverse Power method
    //====================================================================================
    cout << "-----------------------------------------" << endl;
    cout << "Inverse Power method" << endl;
    cout << "-----------------------------------------" << endl;
    gsl_matrix_complex* MLU = gsl_matrix_complex_calloc (sizeM, sizeM);
    gsl_matrix_complex_memcpy(MLU, Minit);

    int s;
    gsl_permutation * p = gsl_permutation_alloc (sizeM);
    gsl_linalg_complex_LU_decomp (MLU, p, &s);

    //Initial guess of the eigenvector of the stable direction
    gsl_vector_complex *xm = gsl_vector_complex_calloc(sizeM);
    for(int i=0; i<sizeM; i++) gsl_vector_complex_set(xm, i, one_c);
    gsl_complex ymnormc;
    for(int i=0; i<sizeM; i++) gsl_vector_complex_set(xm, i, one_c);

    if(version == 0) //Euclidian norm version, needs to normalize
    {
        //for(int i=0; i<sizeM; i++) gsl_vector_complex_set(xm, i, one_c);

        gsl_vector_complex_set(xm, 0, one_c);
        for(int i=1; i<sizeM; i++) gsl_vector_complex_set(xm, i, zero_c);
        GSL_SET_COMPLEX(&ymnormc, 1.0/gsl_blas_dznrm2(xm), 0.0);
        gsl_vector_complex_scale (xm, ymnormc);
    }
    else
    {
        gsl_vector_complex_set(xm, 0, one_c);
        for(int i=1; i<sizeM; i++) gsl_vector_complex_set(xm, i, zero_c);
    }


    //Loop
    gsl_vector_complex *ym = gsl_vector_complex_calloc(sizeM);
    gsl_vector_complex *dxm = gsl_vector_complex_calloc(sizeM);
    gsl_complex lm;
    gsl_complex lm_prec = gslc_complex(0.0, 0.0);
    double dx;
    double dl;
    int iter = 0;
    if(version == 0)
    {
        do
        {
            //y(m) = A-1*x(m-1)
            gsl_linalg_complex_LU_solve (MLU , p , xm , ym);

            //l(m) = x(m-1)^H*y(m)
            gsl_blas_zdotc (xm , ym ,&lm);

            //norm(ym) = sqrt(y(m)^H*y(m))
            GSL_SET_COMPLEX(&ymnormc, 1.0/gsl_blas_dznrm2 (ym), 0.0);

            //----------
            //Error
            //----------
            //dxm = -ym
            gsl_vector_complex_memcpy(dxm, ym);
            for(int i=0; i<sizeM; i++) gsl_vector_complex_set(dxm, i, gsl_complex_negative(gsl_vector_complex_get(dxm, i)));
            //dxm+= lm*xm
            gsl_blas_zaxpy (lm , xm , dxm);
            //dxm/=ymnormc
            gsl_vector_complex_scale (dxm , ymnormc);

            //dx = norm(dx(m))
            dx = gsl_blas_dznrm2 (dxm);
            //----------

            //x(m) = y(m)/norm(y(m))
            gsl_vector_complex_memcpy(xm, ym);
            gsl_vector_complex_scale (xm , ymnormc);

            //cout << iter << " " << dx << " " << GSL_REAL(gsl_complex_inverse(lm)) << " + " << GSL_IMAG(gsl_complex_inverse(lm)) << "i " << endl;

            if(iter ==0)
            {
                dl = 1;
                lm_prec = lm;
            }
            else
            {
                dl = gsl_complex_abs(gsl_complex_sub(gsl_complex_inverse(lm),gsl_complex_inverse(lm_prec)))/gsl_complex_abs(gsl_complex_inverse(lm));
                lm_prec = lm;
            }
            iter++;
            //cout << dl << endl;
        }
        while((dx > prec || dl > prec) && iter < 50);
    }
    else
    {
        //-------------------------------------
        // Infinity norm version
        //-------------------------------------

        gsl_vector *ymabs = gsl_vector_calloc(sizeM);
        int maxIndex = 0;
        do
        {

            //dx(m) = x(m-1);
            gsl_vector_complex_memcpy(dxm, xm);

            //y(m) = A-1*x(m-1)
            gsl_linalg_complex_LU_solve (MLU , p , xm , ym);


            //norm_inf(ym) = max(|y(m)|)
            for(int i = 0; i<sizeM; i++) gsl_vector_set(ymabs, i, fabs(gsl_complex_abs(gsl_vector_complex_get(ym, i))));
            maxIndex = gsl_vector_max_index (ymabs);

            //l(m) = x(m-1)^H*y(m)
            lm = gsl_vector_complex_get(ym, maxIndex);

            //x(m) = y(m)/norm(y(m))
            ymnormc =  gsl_complex_inverse(gsl_vector_complex_get(ym, maxIndex));
            gsl_vector_complex_memcpy(xm, ym);
            gsl_vector_complex_scale (xm , ymnormc);

            //dx(m) = dx(m-1) - dx(m)
            gsl_vector_complex_sub(dxm, xm);

            //dx = norm(dx(m))
            dx = gsl_blas_dznrm2 (dxm);

            cout << iter << " " << dx << " " << GSL_REAL(gsl_complex_inverse(lm)) << " + " << GSL_IMAG(gsl_complex_inverse(lm)) << "i " << endl;

            iter++;

        }
        while(dx > prec && iter < 50);

        gsl_vector_free(ymabs);
    }
    cout << "Eigenvalue: "  <<  GSL_REAL(gsl_complex_inverse(lm)) << " + " << GSL_IMAG(gsl_complex_inverse(lm)) << "i "  << endl;
    cout << "Eigenvector:" << endl;
    gsl_vector_complex_fprintf(stdout, xm, "%+5.5e");


    //Outputs
    gsl_vector_complex_memcpy(eigenV, xm);
    *eigenL =  gsl_complex_inverse(lm);

    gsl_matrix_complex_free(MLU);
    gsl_vector_complex_free(xm);
    gsl_vector_complex_free(ym);
    gsl_vector_complex_free(dxm);
    gsl_permutation_free(p);

}

/**
 *  \brief Shifted Inverse Power Method on a single matrix.
 */
void shifted_power_method(gsl_matrix_complex *Minit, double prec, int sizeM, int version, gsl_complex shift, gsl_vector_complex *eigenV, gsl_complex *eigenL)
{
    //====================================================================================
    // Misc complex numbers
    //====================================================================================
    gsl_complex one_c = gslc_complex(1.0, 0.0);
    gsl_complex zero_c = gslc_complex(0.0, 0.0);

    //====================================================================================
    // Shifted Power method
    //====================================================================================
    cout << "-----------------------------------------" << endl;
    cout << "Shifted Power method" << endl;
    cout << "-----------------------------------------" << endl;
    //lI = lambda_i*Id
    gsl_matrix_complex* lI = gsl_matrix_complex_calloc (sizeM, sizeM);
    gsl_matrix_complex_set_zero(lI);
    for(int i = 0; i<sizeM; i++) gsl_matrix_complex_set(lI, i, i, shift);

    //MLU = M;
    gsl_matrix_complex* MLU = gsl_matrix_complex_calloc (sizeM, sizeM);
    gsl_matrix_complex_memcpy(MLU, Minit);

    //MLU-= lI
    gsl_matrix_complex_sub(MLU, lI);

    //LU decomposition
    int s;
    gsl_permutation * p = gsl_permutation_alloc (sizeM);
    gsl_linalg_complex_LU_decomp (MLU, p, &s);


    //Initial guess of the eigenvector of the stable direction
    gsl_vector_complex *xm = gsl_vector_complex_calloc(sizeM);
    gsl_complex ymnormc;

    gsl_vector_complex_set(xm, 0, one_c);
    for(int i=1; i<sizeM; i++) gsl_vector_complex_set(xm, i, zero_c);
    GSL_SET_COMPLEX(&ymnormc, 1.0/gsl_blas_dznrm2(xm), 0.0);
    gsl_vector_complex_scale (xm, ymnormc);

    //Loop
    gsl_vector_complex *ym = gsl_vector_complex_calloc(sizeM);
    gsl_vector_complex *dxm = gsl_vector_complex_calloc(sizeM);
    gsl_complex lm;
    double dx;
    int iter = 0;
    gsl_complex lmr;

    if(version == 0)
    {
        do
        {
            //dx(m) = x(m-1);
            gsl_vector_complex_memcpy(dxm, xm);

            //y(m) = A-1*x(m-1)
            gsl_linalg_complex_LU_solve (MLU , p , xm , ym);

            //l(m) = x(m-1)^H*y(m)
            gsl_blas_zdotc (xm , ym ,&lm);
            lmr =  gsl_complex_add(gsl_complex_inverse(lm), shift);

            //norm(ym) = sqrt(y(m)^H*y(m))
            GSL_SET_COMPLEX(&ymnormc, 1.0/gsl_blas_dznrm2 (ym), 0.0);

            //----------
            //Error
            //----------
            //dxm = -ym
            gsl_vector_complex_memcpy(dxm, ym);
            for(int i=0; i<sizeM; i++) gsl_vector_complex_set(dxm, i, gsl_complex_negative(gsl_vector_complex_get(dxm, i)));
            //dxm+= lm*xm
            gsl_blas_zaxpy (lm , xm , dxm);
            //dxm/=ymnormc
            gsl_vector_complex_scale (dxm , ymnormc);

            //dx = norm(dx(m))
            dx = gsl_blas_dznrm2 (dxm);
            //----------


            //x(m) = y(m)/norm(y(m))
            gsl_vector_complex_memcpy(xm, ym);
            gsl_vector_complex_scale (xm , ymnormc);


            //cout << iter << " " << dx << " " << GSL_REAL(lmr) << " + " << GSL_IMAG(lmr) << "i " << endl;

            iter++;
        }
        while(dx > prec && iter < 50);

    }

    else
    {
        //-------------------------------------
        // Infinity norm version
        //-------------------------------------

        gsl_vector *ymabs = gsl_vector_calloc(sizeM);
        int maxIndex = 0;
        do
        {

            //dx(m) = x(m-1);
            gsl_vector_complex_memcpy(dxm, xm);

            //y(m) = A-1*x(m-1)
            gsl_linalg_complex_LU_solve (MLU , p , xm , ym);


            //norm_inf(ym) = max(|y(m)|)
            for(int i = 0; i<sizeM; i++) gsl_vector_set(ymabs, i, fabs(gsl_complex_abs(gsl_vector_complex_get(ym, i))));
            maxIndex = gsl_vector_max_index (ymabs);

            //l(m) = x(m-1)^H*y(m)
            lm = gsl_vector_complex_get(ym, maxIndex);
            lmr =  gsl_complex_add(gsl_complex_inverse(lm), shift);

            //x(m) = y(m)/norm(y(m))
            ymnormc =  gsl_complex_inverse(gsl_vector_complex_get(ym, maxIndex));
            gsl_vector_complex_memcpy(xm, ym);
            gsl_vector_complex_scale (xm , ymnormc);

            //dx(m) = dx(m-1) - dx(m)
            gsl_vector_complex_sub(dxm, xm);

            //dx = norm(dx(m))
            dx = gsl_blas_dznrm2 (dxm);

            //cout << iter << " " << dx << " " << GSL_REAL(gsl_complex_inverse(lm)) << " + " << GSL_IMAG(gsl_complex_inverse(lm)) << "i " << endl;

            iter++;

        }
        while(dx > prec && iter < 50);

        gsl_vector_free(ymabs);


    }

    gslc_vector_complex_normalize(xm);
    cout << "Eigenvalue: "  << GSL_REAL(lmr) << " + " << GSL_IMAG(lmr) << "i " << endl;
    cout << "Eigenvector:" << endl;
    gsl_vector_complex_fprintf(stdout, xm, "%+5.8e");

    //Outputs
    gsl_vector_complex_memcpy(eigenV, xm);
    *eigenL =  lmr;

    gsl_matrix_complex_free(MLU);
    gsl_vector_complex_free(xm);
    gsl_vector_complex_free(ym);
    gsl_vector_complex_free(dxm);
    gsl_permutation_free(p);

}

/**
 *  \brief Shifted Inverse Power Method on a single matrix with alternative update of the eigenvalue
 */
void shifted_power_method_update(gsl_matrix_complex *Minit, double prec, int sizeM, gsl_complex shift, gsl_vector_complex *eigenV, gsl_complex *eigenL)
{
    //====================================================================================
    // Misc complex numbers
    //====================================================================================
    gsl_complex one_c = gslc_complex(1.0, 0.0);
    gsl_complex zero_c = gslc_complex(0.0, 0.0);

    //====================================================================================
    // Shifted Power method
    //====================================================================================
    cout << "-----------------------------------------" << endl;
    cout << "Shifted Power method" << endl;
    cout << "-----------------------------------------" << endl;
    //lI = lambda_i*Id
    gsl_matrix_complex* lI = gsl_matrix_complex_calloc (sizeM, sizeM);

    //MLU = M;
    gsl_matrix_complex* MLU = gsl_matrix_complex_calloc (sizeM, sizeM);

    //LU decomposition
    int s;
    gsl_permutation * p = gsl_permutation_alloc (sizeM);

    //Initial guess of the eigenvector
    gsl_vector_complex *xm = gsl_vector_complex_calloc(sizeM);
    gsl_vector_complex *zm = gsl_vector_complex_calloc(sizeM);
    gsl_complex ymnormc;

    gsl_vector_complex_set(xm, 0, one_c);
    for(int i=1; i<sizeM; i++) gsl_vector_complex_set(xm, i, zero_c);
    GSL_SET_COMPLEX(&ymnormc, 1.0/gsl_blas_dznrm2(xm), 0.0);
    gsl_vector_complex_scale (xm, ymnormc);
    //Initial guess for the eigenvalue
    gsl_complex lmr = shift;

    //Loop
    gsl_vector_complex *ym = gsl_vector_complex_calloc(sizeM);
    gsl_vector_complex *dxm = gsl_vector_complex_calloc(sizeM);
    gsl_complex lm;
    int iter = 0;
    double dx;
    do
    {
        //Li = lmr*Id
        gsl_matrix_complex_set_zero(lI);
        for(int i = 0; i<sizeM; i++) gsl_matrix_complex_set(lI, i, i, lmr);

        //y(m) = (M-lmr*Id)-1*x(m-1)
        gsl_matrix_complex_memcpy(MLU, Minit);
        gsl_matrix_complex_sub(MLU, lI);
        gsl_linalg_complex_LU_decomp (MLU, p, &s);
        gsl_linalg_complex_LU_solve (MLU , p , xm , ym);

        //norm(ym) = sqrt(y(m)^H*y(m))
        GSL_SET_COMPLEX(&ymnormc, 1.0/gsl_blas_dznrm2 (ym), 0.0);

        //zm = M*x(m-1)
        gsl_blas_zgemv (CblasNoTrans, one_c , Minit , xm , zero_c , zm);
        //lmr = x(m-1)^H*M*x(m-1)
        gsl_blas_zdotc (xm , zm , &lmr);
        //lm = norm(xm)^2
        gsl_blas_zdotc (xm , xm , &lm);
        //lmr = lmr/lm
        lmr = gsl_complex_div(lmr, lm);

        //----------
        //Error
        //----------
        //zm = Minit*xm
        gsl_blas_zgemv (CblasNoTrans, one_c , Minit, xm , zero_c , zm);
        //zm = Minit*xm/lmr
        gsl_vector_complex_scale(zm, gsl_complex_inverse(lmr));
        //zm = zm-xm
        gsl_vector_complex_sub(zm, xm);
        dx = gsl_blas_dznrm2(zm);

        //x(m) = y(m)/norm(y(m))
        gsl_vector_complex_memcpy(xm, ym);
        gsl_vector_complex_scale (xm , ymnormc);


        cout << iter << " " << dx << " " << GSL_REAL(lmr) << " + " << GSL_IMAG(lmr) << "i " << endl;

        iter++;
    }
    while(dx > prec && iter < 100);



    gslc_vector_complex_normalize(xm);
    cout << "Eigenvalue: "  << GSL_REAL(lmr) << " + " << GSL_IMAG(lmr) << "i " << endl;
    cout << "Eigenvector:" << endl;
    gsl_vector_complex_fprintf(stdout, xm, "%+5.15e");

    //Outputs
    gsl_vector_complex_memcpy(eigenV, xm);
    *eigenL =  lmr;

    gsl_matrix_complex_free(MLU);
    gsl_vector_complex_free(xm);
    gsl_vector_complex_free(ym);
    gsl_vector_complex_free(dxm);
    gsl_permutation_free(p);

}

/**
 *  \brief Decomposition of the Monodromy matrix through the use of GSL routines from a monodromy matrix given as a single matrix. Bad precision
 *
 *  Do not use in general: the precision achieved is not sufficient with a single matrix,
 *  the monodromy matrix should be given as a product, as in the routine mono_eigen.
 */
void gslc_mono_eigen(gsl_odeiv2_driver *d, const double y0[], double tend, gsl_matrix_complex* gsl_MMc, gsl_matrix* gsl_MM, gsl_vector_complex *eval, gsl_matrix_complex *evec)
{
    //====================================================================================
    // Monodromy matrix
    //====================================================================================
    double ys[42];
    stm_complex(y0, d, tend, gsl_MMc, ys);
    stm_matrix_with_gsl(y0, d, tend, gsl_MM);

    //====================================================================================
    //Eigenvectors & values from GSL routines
    //====================================================================================
    gsl_matrix* monodromy_copy = gsl_matrix_calloc (6, 6);
    gsl_matrix_memcpy(monodromy_copy, gsl_MM);
    gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (6);
    gsl_eigen_nonsymmv (monodromy_copy, eval, evec, w);
    gsl_eigen_nonsymmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);

    cout << "------------------------------------------" << endl;
    cout << "gslc_mono_eigen. Monodromy matrix" << endl;
    //gslc_matrix_printf(gsl_MM);
    gslc_matrix_complex_printf_real(gsl_MMc);

    cout << "------------------------------------------" << endl;
    cout << "gslc_mono_eigen. Eigenvalues from GSL routine" << endl;
    gslc_vector_complex_printf(eval);

    gsl_matrix_free(monodromy_copy);
    gsl_eigen_nonsymmv_free(w);
}

//----------------------------------------------------------------------------------------
//Resolution of eigensystem: inspired by Fortran code
//----------------------------------------------------------------------------------------
/**
 *  \brief Decomposition of a product of matrices DAT[M]*DAT[M-1]*...DAT[1] into unstable, stable and central part.
 *  \param DAT  the product of M matrices, given from indix 1 to indix M - and NOT 0 to M-1 as usually done in C.
 *  \param M    the number of matrices in DAT
 *  \param VECS output: hyperbolic and central directions at the beginning of each interval defined by DAT[K]:
 *         VECS[K] = DAT[K]*VECS[K] except for a multiplicative constant in each direction.
 *         VECS contains M+1 6*6 matrices and is of the of the form VECS(i,j,k) with
 *           - i = 0,5.
 *           - j = 0,5.
 *           - k = 0,M. Note that the indix M is included.
 *         The directions j = 0 and 1 refer to the unstable and stable directions respectively, and j = 2,5 to the central ones.
 *
 *  Inspired by Pr. Josep Masdemont's Fortran code and algorithm available at p. 86 of "Dynamics and mission design near libration points III", 2001.
 *  Recall that VECS(*,j,k)= DAT(*,*,k)*VECS(*,j,k-1) except for a multiplicative constant for each direction.
 */
void vepro(gsl_matrix_complex **DAT, int M, gsl_matrix_complex **VECS, gsl_matrix_complex *DECS, gsl_matrix_complex *DECSv, gsl_vector_complex *DECSl)
{
    //------------------------------------------------------------------------------------
    // Misc complex numbers
    //------------------------------------------------------------------------------------
    gsl_complex one_c  = gslc_complex(1.0, 0.0);
    gsl_complex zero_c = gslc_complex(0.0, 0.0);

    //------------------------------------------------------------------------------------
    // GSL objects
    //------------------------------------------------------------------------------------
    //Fortran
    gsl_matrix_complex **VECS_FORTRAN = gslc_matrix_complex_product_alloc(6, 6, M);
    gsl_matrix_complex **VECS_DELTA   = gslc_matrix_complex_product_alloc(6, 6, M);
    gsl_eigen_nonsymmv_workspace * wr = gsl_eigen_nonsymmv_alloc(6);
    //Hyperbolic eigenvalues
    gsl_complex VAP1, VAP2, VAPL1, VAPL2;
    //VEP1 = [1, 0, 0, ..., 0]
    gsl_vector_complex *VEP1 = gsl_vector_complex_calloc(6);
    gsl_vector_complex_set_zero(VEP1);
    gsl_vector_complex_set(VEP1, 0, one_c);
    //VEP2 = [0, 1, 0, ..., 0]
    gsl_vector_complex *VEP2 = gsl_vector_complex_calloc(6);
    gsl_vector_complex_set_zero(VEP2);
    gsl_vector_complex_set(VEP2, 1, one_c);
    //Matrices and vectors
    gsl_matrix *DECSr         = gsl_matrix_calloc(6, 6);
    gsl_matrix_complex *DECSm = gsl_matrix_complex_calloc(6, 6);
    gsl_matrix_complex *AUX   = gsl_matrix_complex_calloc(6, 6);
    gsl_vector_complex *BUX1  = gsl_vector_complex_calloc(6);
    gsl_vector_complex *BUX2  = gsl_vector_complex_calloc(6);
    gsl_vector_complex *CUX1  = gsl_vector_complex_calloc(6);
    gsl_vector_complex *CUX2  = gsl_vector_complex_calloc(6);
    //Temporary variables
    double **DELTA = dmatrix(0, 1, 0, M);
    double VN1, VN2, DEL;
    //Inverse and complexified version of VN1, VN2
    gsl_complex VN1invC, VN2invC, GAMMA;
    //For inverse power method
    int s;
    gsl_permutation * p = gsl_permutation_alloc(6);

    //====================================================================================
    // Computing stable and unstable direction
    //====================================================================================
    cout <<  Csts::SEPR << endl;
    cout << "vepro. Computing the initial and final hyperbolic directions." << endl;
    power_method_on_prod(DAT+1,M,6,VEP1,&VAP1,&VAPL1, 1);  //Unstable, direct power method, IS = +1
    power_method_on_prod(DAT+1,M,6,VEP2,&VAP2,&VAPL2,-1);  //Stable, inverse power method, IS = -1

    cout << "vepro. Hyperbolic eigenvalues: " << endl;
    cout <<  GSL_REAL(VAP1) << endl;
    cout <<  GSL_REAL(gsl_complex_inverse(VAP2)) << endl;
    cout <<  Csts::SSEPR << endl;

    cout << "vepro. Computig intermediate hyp. directions." << endl;
    gsl_vector_complex_view VECS_V;
    gsl_vector_complex_view VECS_V2;

    //Column 0 of VECS[0] = VEP1
    VECS_V = gsl_matrix_complex_column (VECS[0], 0);
    gsl_vector_complex_memcpy(&VECS_V.vector, VEP1);

    //Column 1 of VECS[M] = VEP2
    VECS_V = gsl_matrix_complex_column (VECS[M], 1);
    gsl_vector_complex_memcpy(&VECS_V.vector, VEP2);

    //Getting the hyperbolic directions along the path
    for(int K = 1; K <= M; K++)
    {
        //BUX1 = column 0 of VEC[K-1]
        gslc_matrix_complex_column(BUX1, VECS[K-1], 0);
        //VECS_V = gsl_matrix_complex_column (VECS[K-1], 0);
        //gsl_vector_complex_memcpy(BUX1, &VECS_V.vector);

        //BUX2 = column 1 of VEC[M-K+1]
        gslc_matrix_complex_column(BUX2, VECS[M-K+1], 1);
        //VECS_V = gsl_matrix_complex_column (VECS[M-K+1], 1);
        //gsl_vector_complex_memcpy(BUX2, &VECS_V.vector);

        //AUX = DAT(K)
        gsl_matrix_complex_memcpy(AUX, DAT[K]);
        //CUX1 = AUX*BUX1
        gsl_blas_zgemv (CblasNoTrans, one_c , AUX , BUX1, zero_c , CUX1);

        //AUX = DAT(M-K-1)
        gsl_matrix_complex_memcpy(AUX, DAT[M-K+1]);
        //LU decomposition of AUX
        gsl_linalg_complex_LU_decomp (AUX, p, &s);
        //CUX2 = AUX^(-1)*BUX2
        gsl_linalg_complex_LU_solve (AUX , p , BUX2 , CUX2);

        //VN1 = NORM(CUX1)
        VN1 = gsl_blas_dznrm2 (CUX1);
        //VN2 = NORM(CUX2)
        VN2 = gsl_blas_dznrm2 (CUX2);

        //Update the DELTAs
        DELTA[0][K]     = VN1;
        DELTA[1][M-K+1] = VN2;

        //VN1invC = Inverse and complexified version of VN1
        GSL_SET_COMPLEX(&VN1invC, 1.0/VN1, 0.0);
        //Colmun 0 of VECS[K] = CUX1/NORM(CUX1);
        VECS_V = gsl_matrix_complex_column (VECS[K], 0);
        gsl_vector_complex_memcpy(&VECS_V.vector, CUX1);
        gsl_vector_complex_scale (&VECS_V.vector, VN1invC);

        //VN2invC = Inverse and complexified version of VN2
        GSL_SET_COMPLEX(&VN2invC, 1.0/VN2, 0.0);
        //Colmun 1 of VECS[M-K] = CUX2/NORM(CUX2);
        VECS_V = gsl_matrix_complex_column (VECS[M-K], 1);
        gsl_vector_complex_memcpy(&VECS_V.vector, CUX2);
        gsl_vector_complex_scale (&VECS_V.vector, VN2invC);
    }

    //====================================================================================
    // Central components
    //====================================================================================
    for(int IJ = 2; IJ < 6; IJ++)
    {
        //Column IJ of VECS[0] = [0.0, 0.0, ..., 1.0 @ IJ, ..., 0.0];
        VECS_V = gsl_matrix_complex_column (VECS[0], IJ);
        gsl_vector_complex_set_zero(&VECS_V.vector);
        gsl_vector_complex_set(&VECS_V.vector, IJ, one_c);
        //BUX1 = [0.0, 0.0, ..., 1.0 @ IJ, ..., 0.0];
        gsl_vector_complex_set_zero(BUX1);
        gsl_vector_complex_set(BUX1, IJ, one_c);

        //------------------------------------------------------------------------------------
        // NOTE: the column IJ of VECS[K] is denoted VECS(*, IJ, K)
        // Example @ this point for IJ
        // VECS(*, IJ, 0) = [0.0, 0.0, ..., 1.0 @ IJ, ..., 0.0]^T = w_0^0;
        //
        // See Dynamics and mission design near libration points III p.86 for the associated formalism
        //------------------------------------------------------------------------------------

        //Loop
        for(int K= 1; K<=M; K++)
        {
            //AUX = DAT(K)
            gsl_matrix_complex_memcpy(AUX, DAT[K]);
            //CUX1 = AUX*BUX1
            gsl_blas_zgemv (CblasNoTrans, one_c , AUX , BUX1, zero_c , CUX1);
            //GAMMA = CUX1*(colum 0 of VECS[K])
            VECS_V = gsl_matrix_complex_column (VECS[K], 0);
            gsl_blas_zdotu(CUX1 , &VECS_V.vector , &GAMMA);
            //------------------------------------------------------------------------------------
            // Example @ this point for IJ, K = 1
            // CUX1  = A_1 w_0^0
            // gamma = A_1 w_0^0 u_1^T = lambda_1
            //
            // Example @ this point for IJ, K = 2
            // CUX1  = A_2 (wb_1^0 - lambda_1 u_1) = A_2 w_1^0 = wb_2^0
            // gamma = wb_2^0 u_2^T = lambda_2
            //
            // See Dynamics and mission design near libration points III p.86 for the associated formalism
            //------------------------------------------------------------------------------------


            //BUX1 = CUX1 - GAMMA*(colum 0 of VECS[K])
            gsl_vector_complex_memcpy(BUX1, CUX1);                                   //BUX1 = CUX1
            gsl_blas_zaxpy (gsl_complex_negative(GAMMA) , &VECS_V.vector , BUX1);    //BUX1 = -GAMMA*VECS_V + BUX1
            //(column IJ of VECS[K]) = BUX1
            VECS_V = gsl_matrix_complex_column (VECS[K], IJ);
            gsl_vector_complex_memcpy(&VECS_V.vector, BUX1);
            //------------------------------------------------------------------------------------
            // Example @ this point for IJ, K = 1
            // CUX1  = A_1 w_0^0
            // gamma = A_1 w_0^0 u_1^T = lambda_1
            // BUX1  = CUX1 - gamma*VECS(Id,0,1) = wb_1^0 - lambda_1 u_1 = w_1^0
            // VECS(*, IJ, 1) = BUX1 = w_1^0
            //
            // Example @ this point for IJ, K = 2
            // CUX1  = A_2 (wb_1^0 - lambda_1 u_1) = A_2 w_1^0 = wb_2^0
            // gamma = wb_2^0 u_2^T = lambda_2
            // BUX1  = CUX1 - gamma*VECS(Id,0,2) = wb_2^0 - lambda_2 u_2 = w_2^0
            // VECS(*, IJ, 2) = BUX1 = wb_2^0 - lambda_2 u_2 = w_2^0
            //
            // See Dynamics and mission design near libration points III p.86 for the associated formalism
            //------------------------------------------------------------------------------------


            DEL=1.0;
            for(int IK = K-1; IK>=0; IK--)
            {
                //Update of DEL
                DEL*=DELTA[0][IK+1];
                //(column IJ of VECS[IK]) += -GAMMA/DEL*(column 0 of VECS[IK])
                VECS_V  = gsl_matrix_complex_column (VECS[IK], IJ);
                VECS_V2 = gsl_matrix_complex_column (VECS[IK], 0);
                //VECS_V = -GAMMA/DEL*VECS_V2 + VECS_V
                gsl_blas_zaxpy (gsl_complex_negative(gsl_complex_div_real(GAMMA, DEL)), &VECS_V2.vector , &VECS_V.vector);
            }

            //------------------------------------------------------------------------------------
            // Example @ this point for IJ, K = 1
            // CUX1  = A_1 w_0^0
            // gamma = A_1 w_0^0 u_1^T = lambda_1
            // BUX1  = wb_1^0 - lambda_1 u_1 = w_1^0
            // VECS(*, IJ, 1) = w_1^0
            // VECS(*, IJ, 0) = w_0^0 - lambda_1/delta_1 u_0 = w_0^1
            //
            // Example @ this point for IJ, K = 2
            // CUX1  = A_2 (wb_1^0 - lambda_1 u_1) = A_2 w_1^0 = wb_2^0
            // gamma = wb_2^0 u_2^T = lambda_2
            // BUX1  = wb_2^0 - lambda_2 u_2 = w_2^0
            // VECS(*, IJ, 2) =  wb_2^0 - lambda_2 u_2 = w_2^0
            // VECS(*, IJ, 1) =  w_1^0  - lambda_2/delta_2 u_1 = w_1^1
            // VECS(*, IJ, 0) =  w_0^1  - lambda_2/delta_2 u_1 = w_0^2
            //
            //
            // See Dynamics and mission design near libration points III p.86 for the associated formalism
            //------------------------------------------------------------------------------------
        }

    }

    //====================================================================================
    // - VECS(*,j,0)*DECS = VECS(*,j,M), for j = 2,5.
    // - If the last input of diahip = 1, changing the intput VECS in such a way
    //   that the new vectors vecs(*,j,m) j=2,5 span the same 4th dimensional space
    //   that the vectors vecs(*,j,0) j=2,5
    //====================================================================================
    diahip(VECS, DELTA, DECS, M, 1);
    //print of DECS in txt file
    gslc_matrix_complex_fprintf_pretty(DECS, (char*) "fprint/DECSc.txt");

    //====================================================================================
    // (Ortho)normalisation
    // The line ortho(VECS, K, M) must be uncommented to get an orthogonal base.
    // Check also Fortran code for consistency between the 2 codes.
    //====================================================================================
    for(int K = 0; K <= M; K++)
    {
        for(int J = 2; J <6; J++)
        {
            VECS_V = gsl_matrix_complex_column (VECS[K], J);
            gslc_vector_complex_normalize(&VECS_V.vector);
        }
        //ortho(VECS, K, M); //for getting an orthonormalized base. Continuity is lost in this case.
    }

    //====================================================================================
    // Diag(DECS), output in DECSl, DECSv
    //====================================================================================
    gslc_matrix_complex_real(DECSr, DECS);
    mono_diag(DECSr, DECSv, DECSl, DECSm);

    //Print in txt file
    gslc_eigensystem_fprintf(DECSl, DECSv, 6,   (char*)"fprint/diagDECS.txt");
    gslc_matrix_complex_fprintf_pretty(VECS[M], (char*)"fprint/VECS[M]_after_diahip2.txt");
    gslc_matrix_complex_fprintf_pretty(VECS[0], (char*)"fprint/VECS[0]_after_diahip2.txt");


    //====================================================================================
    // Comparison with FORTRAN code
    //====================================================================================
    //Creating an array of DAT for fortran code. Note: all dimension are inversed to be able to use fortran array order.
    double datf[M][6][6];     //equivalent of DAT
    double vecsf[M+1][6][6];  //equivalent of VECS
    for(int k = 0; k < M; k++)
    {
        for(int i = 0; i <6; i ++)
        {
            for(int j = 0; j <6; j ++) datf[k][j][i] = GSL_REAL(gsl_matrix_complex_get(DAT[k+1], i, j));
        }
    }

    //Fortran vepro (does not work on server €€TODO)
    //vepro_((double**)datf, &M, (double**)vecsf, &M);
    cout << "vepro. datf[0][0][0] is displayed just to avoid compilation warning: " << datf[0][0][0] << endl;


    //Set vecsf in VECS_FORTRAN
    for(int k = 0; k <= M; k++)
    {
        for(int i = 0; i <6; i ++)
        {
            for(int j = 0; j <6; j ++) gsl_matrix_complex_set(VECS_FORTRAN[k], i, j, gslc_complex(vecsf[k][j][i]));
        }
    }

    //VECS_DELTA = VECS-VECS_FORTRAN
    for(int k = 0; k <=M; k++)
    {
        gsl_matrix_complex_memcpy(VECS_DELTA[k], VECS_FORTRAN[k]);
        gsl_matrix_complex_sub(VECS_DELTA[k], VECS[k]);
    }

    //Comparison
    cout << Csts::SSEPR << endl;
    cout << "vepro. The difference between C and Fortran computation is:" << endl;
    cout << setprecision(5);
    cout << gslc_matrix_complex_infinity_norm(VECS_DELTA[0])   << " at t = 0 "   << endl;
    cout << gslc_matrix_complex_infinity_norm(VECS_DELTA[M/2]) << " at t = T/2 "   << endl;
    cout << gslc_matrix_complex_infinity_norm(VECS_DELTA[M])   << " at t = T "   << endl;
    cout << setprecision(15);
    cout << Csts::SSEPR << endl;

    //------------------------------------------------------------------------------------------------------------------------
    // VECS = VECS_FORTRAN (uncomment if needed)
    //------------------------------------------------------------------------------------------------------------------------
    //    for(int k = 0; k <= M; k++)
    //    {
    //        for(int i = 0; i <6; i ++)
    //        {
    //            for(int j = 0; j <6; j ++) gsl_matrix_complex_set(VECS[k], i, j, gsl_matrix_complex_get(VECS_FORTRAN[k], i, j));
    //        }
    //    }

    //------------------------------------------------------------------------------------------------------------------------
    // Solving the eigensystem central part via the information on VECS[0]
    //------------------------------------------------------------------------------------------------------------------------
    center_diag(DAT, VECS[0], DECSl, DECSv, M);
    //center_diag(DAT, DECS, DECSl, DECSv, M);


    //------------------------------------------------------------------------------------------------------------------------
    // Memory release
    //------------------------------------------------------------------------------------------------------------------------
    free_dmatrix(DELTA, 0, 1, 0, M);
    gsl_permutation_free(p);
    gsl_matrix_free(DECSr);
    gsl_matrix_complex_free(AUX);
    gsl_vector_complex_free(BUX1);
    gsl_vector_complex_free(BUX2);
    gsl_vector_complex_free(CUX1);
    gsl_vector_complex_free(CUX2);
    gsl_vector_complex_free(VEP1);
    gsl_vector_complex_free(VEP2);
    gsl_eigen_nonsymmv_free(wr);
    gslc_matrix_complex_product_free(VECS_FORTRAN, M);
    gslc_matrix_complex_product_free(VECS_DELTA, M);

}

/**
 *  \brief Power (resp. Inverse Power) Method on a product of matrices to get the biggest (resp. smallest) eigenvalue in absolute value
 *         along with its associated eigenvectors. note that VAPL = log10(VAP)
 */
void power_method_on_prod(gsl_matrix_complex **DAT, int M, int N, gsl_vector_complex *VEP, gsl_complex *VAP, gsl_complex *VAPL,  int IS)
{
    //Number MAX of iterations allowed
    int ITM = 20;
    //Precision required
    double PREC = 1e-16;

    //------------------------------------------------------------------------------------
    // Misc complex numbers
    //------------------------------------------------------------------------------------
    gsl_complex one_c = gslc_complex(1.0, 0.0);
    gsl_complex zero_c = gslc_complex(0.0, 0.0);

    //------------------------------------------------------------------------------------
    // Specific gsl tools for the power_method_on_prod routine
    //------------------------------------------------------------------------------------
    gsl_matrix_complex *AUX = gsl_matrix_complex_calloc(N,N);
    gsl_vector_complex *BUX = gsl_vector_complex_calloc(N);
    gsl_vector_complex *CUX = gsl_vector_complex_calloc(N);
    double VN;
    gsl_complex VNinvC = gslc_complex(0.0, 0.0); //Inverse and complexified version of VN
    gsl_complex VAPV = gslc_complex(0.0, 0.0);
    double DVAP;

    //For inverse power method, may not be used
    int s;
    gsl_permutation * p = gsl_permutation_alloc (N);
    //Power method
    int IT = 0;
    do
    {
        //VAP = 1.0
        GSL_SET_COMPLEX(VAP, 1.0, 0.0);
        //VAPL = 0.0;
        GSL_SET_COMPLEX(VAPL, 0.0, 0.0);

        //Loop on the matrices
        int K ;
        for(int K1 = 0; K1 < M; K1++)
        {
            if(IS == 1) //direct  power method
            {
                //True indix
                K = K1;
                //AUX = DAT(K)
                gsl_matrix_complex_memcpy(AUX, DAT[K]);
                //BUX = AUX*VEP
                gsl_blas_zgemv (CblasNoTrans, one_c , AUX , VEP, zero_c , BUX);

            }
            else    //inverse power method
            {
                //True indix
                K = (M-1) - K1;
                //AUX = DAT(K). Here, memcpy is mandatory, otherwise the original DAT would be modified by the LU decomposition
                gsl_matrix_complex_memcpy(AUX, DAT[K]);
                //LU decomposition of AUX
                gsl_linalg_complex_LU_decomp (AUX, p, &s);
                //BUX = AUX^(-1)*VEP
                gsl_linalg_complex_LU_solve (AUX , p , VEP , BUX);
            }

            //VN = NORM(BUX)
            VN = gsl_blas_dznrm2 (BUX);
            //VNinvC = Inverse and complexified version of VN
            GSL_SET_COMPLEX(&VNinvC, 1.0/VN, 0.0);
            //VAP*=VN
            *VAP   = gsl_complex_mul_real(*VAP, VN);
            //VAP+=log10(VN);
            *VAPL  = gsl_complex_add_real(*VAPL, log10(VN));
            //VEP = BUX/NORM(BUX);
            gsl_vector_complex_memcpy(VEP, BUX);
            gsl_vector_complex_scale (VEP , VNinvC);
        }

        if(IT ==0) //first iteration, needs to update some parameters
        {
            //VAPV = VAPL
            VAPV = *VAPL;
            //CUX = VEP;
            gsl_vector_complex_memcpy(CUX, VEP);
        }
        else
        {
            //DVAP = |VAPL-VAPV|/|VAPL|
            DVAP = gsl_complex_abs(gsl_complex_sub(*VAPL, VAPV))/gsl_complex_abs(*VAPL);
            //CUX = CUX-VEP
            gsl_vector_complex_sub(CUX, VEP);
            //VN = NORM(CUX);
            VN = gsl_blas_dznrm2 (CUX);

            if(IT > ITM-5)
            {
                cout << "power_method_on_prod. REL. ERROR in VAPL: " << DVAP << endl;
                cout << "power_method_on_prod. ABS. ERROR in VECTOR: " << VN << endl;
            }

            //If desired precision is achieved, return!
            if(DVAP < PREC && VN < PREC)
            {
                //cout << "power_method_on_prod. DESIRED PRECISION IS ACHIEVED, END OF THE ITERATIVE PROCESS." << endl;
                gsl_matrix_complex_free(AUX);
                gsl_vector_complex_free(BUX);
                gsl_vector_complex_free(CUX);
                gsl_permutation_free(p);
                return;
            }

            //VAPV = VAPL
            VAPV = *VAPL;
            //CUX = VEP;
            gsl_vector_complex_memcpy(CUX, VEP);
        }
        IT++;

        if(IT > ITM)
        {
            cout << "power_method_on_prod. POWER METHOD SEEMS NOT CONVERGENT" << endl;
            cout << "        VALUE OF IS: " << IS << endl;
            gsl_matrix_complex_free(AUX);
            gsl_vector_complex_free(BUX);
            gsl_vector_complex_free(CUX);
            gsl_permutation_free(p);
            return;
        }
    }
    while(1);

    //Never here
    gsl_matrix_complex_free(AUX);
    gsl_vector_complex_free(BUX);
    gsl_vector_complex_free(CUX);
    gsl_permutation_free(p);
    return;
}

/**
 * \brief From J. Masdemont Fortran code. Change the intput VECS in such a way that
 *        the new vectors VECS(*,j,m) j=2,5 span the same 4th dimensional space that
 *        the vectors VECS(*,j,0) j=2,5.
 *
 *  1. The routine updates the matrix DECS so that VECS(*,j,0)*DECS = VECS(*,j,M), for j = 2,5.
 *
 *  2. If isSameBase == true: Given m+1 sets of basis of 6 vectors stored by columns in VECS(6,6,0:m)
 *  assuming that VECS(*,j,k+1)=a_k*VECS(*,j,k) where a_k is a certain
 *  6*6 matrix, and VECS(*,j,m)=di*VECS(*,j,0) for j=0,1, where
 *       - d1=delta(0,0)*delta(0,1)*..*delta(0,m-1) and
 *       - d2=1.d0/(delta(1,0)*delta(1,1)*...*delta(1,m-1))
 *  i.e. VECS(*,0,0) and VECS(*,1,0) are eigenvectors of the product
 *  of the a_k matrices, and delta must be given at the input,
 *  this routine modifies the vectors VECS(*,j,k) j=2,5 k=0,m, keeping
 *  always VECS(*,j,k+1)=a_k*VECS(*,j,k), in such a way that the new
 *  vectors VECS(*,j,m) j=2,5 span the same 4th dimensional space that
 *  the vectors VECS(*,j,0) j=2,5.
 *
 **/
void diahip(gsl_matrix_complex **VECS, double **DELTA, gsl_matrix_complex *DECS, int M, int isSameBase)
{
    gsl_matrix_complex *AUX   = gsl_matrix_complex_calloc(6, 6);
    gsl_matrix_complex *A4M   = gsl_matrix_complex_calloc(4, 4);
    gsl_matrix_complex *B4M   = gsl_matrix_complex_calloc(4, 4);
    gsl_matrix_complex *ACM   = gsl_matrix_complex_calloc(4, 4);

    gsl_vector_complex *BUX1  = gsl_vector_complex_calloc(6);
    gsl_vector_complex *ALF   = gsl_vector_complex_calloc(6);
    gsl_vector_complex *BET   = gsl_vector_complex_calloc(6);
    gsl_vector_complex *VA    = gsl_vector_complex_calloc(4);
    gsl_vector_complex *VB    = gsl_vector_complex_calloc(4);
    gsl_vector_complex *ALF4  = gsl_vector_complex_calloc(4);
    gsl_vector_complex *BET4  = gsl_vector_complex_calloc(4);


    gsl_vector_complex_view VECS_V;
    gsl_vector_complex_view VECS_V2;
    gsl_vector_complex_view DECS_V;
    //For inverse power method
    int s;
    gsl_permutation * p  = gsl_permutation_alloc(6);
    gsl_permutation * p14 = gsl_permutation_alloc(4);
    gsl_permutation * p24 = gsl_permutation_alloc(4);

    //Building the matrix DECS = VECS[0]^(-1)*VECS[M]
    for(int J = 0; J < 6; J++)
    {
        //AUX = VECS[0];
        gsl_matrix_complex_memcpy(AUX, VECS[0]);
        //LU decomposition of AUX
        gsl_linalg_complex_LU_decomp (AUX, p, &s);
        //VECS_V = view of VECS(*, J, M)
        VECS_V = gsl_matrix_complex_column (VECS[M], J);
        //BUX1 = AUX^(-1)*VECS_V
        gsl_linalg_complex_LU_solve (AUX , p , &VECS_V.vector , BUX1);
        //DECS_V = view of DECS(*, J)
        DECS_V = gsl_matrix_complex_column (DECS, J);
        //DECS(*,J) = BUX1 = VECS[0]^(-1)*VECS_V
        gsl_vector_complex_memcpy(&DECS_V.vector, BUX1);
    }

    gslc_matrix_complex_fprintf_pretty(VECS[M], (char*) "fprint/VECS[M]_before_diahip.txt");
    gslc_matrix_complex_fprintf_pretty(VECS[0], (char*) "fprint/VECS[0]_before_diahip.txt");

    double AQD = 1.0;
    double BQD = 1.0;
    for(int K = 1; K <= M; K++)
    {
        AQD /= DELTA[0][K];
        BQD /= DELTA[1][K];
    }

    //Building A4M, B4M and ACM
    gsl_complex decs_ij, decs_ji, temp;
    for(int J = 0; J < 4; J++)
    {
        for(int II = 0; II < 4; II++)
        {
            decs_ij = gsl_matrix_complex_get(DECS, II+2, J+2);
            decs_ji = gsl_matrix_complex_get(DECS, J+2, II+2);
            //A4M(II,J)=DESC(J+2,II+2)*AQD
            gsl_matrix_complex_set(A4M, II, J, gsl_complex_mul_real(decs_ji, AQD));
            //B4M(II,J)=DESC(J+2,II+2)
            gsl_matrix_complex_set(B4M, II, J, decs_ji);
            //ACM(II,J)=DESC(II+2,J+2)
            gsl_matrix_complex_set(ACM, II, J, decs_ij);
        }

        //A4M(J,J)=A4M(J,J)-1.D0
        temp = gsl_matrix_complex_get(A4M, J, J);
        gsl_matrix_complex_set(A4M, J, J, gsl_complex_sub_real(temp, 1.0));
        //B4M(J,J)=B4M(J,J)-BQD
        temp = gsl_matrix_complex_get(B4M, J, J);
        gsl_matrix_complex_set(B4M, J, J, gsl_complex_sub_real(temp, BQD));
        //VA(J)=DESC(0,J+2)
        temp = gsl_matrix_complex_get(DECS, 0, J+2);
        gsl_vector_complex_set(VA, J, temp);
        //VB(J)=DESC(1,J+2)
        temp = gsl_matrix_complex_get(DECS, 1, J+2);
        gsl_vector_complex_set(VB, J, temp);
    }

    gslc_matrix_complex_fprintf_pretty(A4M, (char*) "fprint/A4M.txt");
    gslc_matrix_complex_fprintf_pretty(B4M, (char*) "fprint/B4M.txt");

    //LU decomposition of A4M
    gsl_linalg_complex_LU_decomp (A4M, p14, &s);
    //ALF4 = A4M^(-1)*VA
    gsl_linalg_complex_LU_solve (A4M , p14 , VA , ALF4);
    //LU decomposition of B4M
    gsl_linalg_complex_LU_decomp (B4M, p24, &s);
    //BET4 = B4M^(-1)*VB
    gsl_linalg_complex_LU_solve (B4M , p24 , VB , BET4);

    //Shifted copy of ALF4 in ALF, and BET4 in BET
    for(int J = 0; J < 4; J++) gsl_vector_complex_set(ALF, J+2, gsl_vector_complex_get(ALF4, J));
    for(int J = 0; J < 4; J++) gsl_vector_complex_set(BET, J+2, gsl_vector_complex_get(BET4, J));

    gslc_vector_complex_fprintf(ALF, (char*) "fprint/ALF.txt");
    gslc_vector_complex_fprintf(BET, (char*) "fprint/BET.txt");

    //Update VECS if needed
    if(isSameBase)
    {
        for(int J = 2; J < 6; J++)
        {
            for(int K = M; K >= 0; K--)
            {
                //VECS(*,J,K)=VECS(*,J,K)+ALF(J)*VECS(*,0,K)
                temp    = gsl_vector_complex_get(ALF, J);
                VECS_V  = gsl_matrix_complex_column (VECS[K], J);
                VECS_V2 = gsl_matrix_complex_column (VECS[K], 0);
                gsl_blas_zaxpy (temp , &VECS_V2.vector , &VECS_V.vector);  //VECS_V = VECS_V + ALF(J)*VECS_V2

                //if K !=0,  ALF(J)=ALF(J)/DELTA(0,K)
                if(K!=0)
                {
                    gsl_vector_complex_set(ALF, J, gsl_complex_div_real(gsl_vector_complex_get(ALF, J), DELTA[0][K]));
                }
            }

            for(int K = 0; K <= M;  K++)
            {
                //VECS(*,J,K)=VECS(*,J,K)+BET(J)*VECS(*,1,K)
                temp    = gsl_vector_complex_get(BET, J);
                VECS_V  = gsl_matrix_complex_column (VECS[K], J);
                VECS_V2 = gsl_matrix_complex_column (VECS[K], 1);
                gsl_blas_zaxpy (temp , &VECS_V2.vector , &VECS_V.vector);  //VECS_V = VECS_V + BET(J)*VECS_V2


                //if K !=M,  BET(J)=BET(J)/DELTA(1,K+1)
                if(K!=M)
                {
                    gsl_vector_complex_set(BET, J, gsl_complex_div_real(gsl_vector_complex_get(BET, J), DELTA[1][K+1]));
                }

            }
        }
    }

    //Building (again) the matrix DECS = VECS[0]^(-1)*VECS[M]
    for(int J = 0; J < 6; J++)
    {
        //AUX = VECS[0];
        gsl_matrix_complex_memcpy(AUX, VECS[0]);
        //LU decomposition of AUX
        gsl_linalg_complex_LU_decomp (AUX, p, &s);
        //VECS_V = view of VECS(*, J, M)
        VECS_V = gsl_matrix_complex_column (VECS[M], J);
        //BUX1 = AUX^(-1)*VECS_V
        gsl_linalg_complex_LU_solve (AUX , p , &VECS_V.vector , BUX1);
        //DECS_V = view of DECS(*, J)
        DECS_V = gsl_matrix_complex_column (DECS, J);
        //DECS(*,J) = BUX1 = VECS[0]^(-1)*VECS_V
        gsl_vector_complex_memcpy(&DECS_V.vector, BUX1);
    }


    //Diag(ACM)
    gsl_matrix *ACMr = gsl_matrix_calloc(4, 4);
    gsl_eigen_nonsymmv_workspace * wr = gsl_eigen_nonsymmv_alloc(4);
    gsl_matrix_complex *ACMv   = gsl_matrix_complex_calloc(4, 4);
    gsl_vector_complex *ACMl   = gsl_vector_complex_calloc(4);
    for(int i=0; i<4; i++) for(int j=0; j<4; j++) gsl_matrix_set(ACMr, i, j,  GSL_REAL(gsl_matrix_complex_get(ACM, i, j)));
    gsl_matrix_complex_set_zero(ACMv);
    gsl_eigen_nonsymmv (ACMr, ACMl, ACMv, wr);
    gsl_eigen_nonsymmv_sort (ACMl, ACMv, GSL_EIGEN_SORT_ABS_ASC);
    //Print in txt file
    gslc_eigensystem_fprintf(ACMl, ACMv, 4, (char*)"fprint/diagACM.txt");
    gslc_matrix_complex_fprintf(ACM, (char*) "fprint/ACM.txt");
}

/**
 *  \brief Auxiliary routine. It orthonomalizes the bases VECS(*,*,K), for all K=0,M.
 **/
void ortho(gsl_matrix_complex **VECS, int IK, int M)
{
    gsl_vector_complex *B = gsl_vector_complex_calloc(6);
    gsl_vector_complex_view VECS_V;
    gsl_complex C = gslc_complex(0.0,0.0);
    gsl_complex temp1, temp2;

    for(int II = 2; II <6; II++)
    {
        //B = VECS(*, II, IK)
        VECS_V = gsl_matrix_complex_column (VECS[IK], II);
        gsl_vector_complex_memcpy(B, &VECS_V.vector);


        for(int K = 1; K <= II-1; K++)
        {
            //C = B*VECS(*, K, IK)
            C = gslc_complex(0.0,0.0);
            for(int L = 0; L <6; L++)
            {
                temp1 = gsl_vector_complex_get(B, L);
                temp2 = gsl_matrix_complex_get(VECS[IK], L, K);
                C = gsl_complex_add(C, gsl_complex_mul(temp1, temp2));
            }

            //B = B -C*A(*, K, IK)
            for(int L = 0; L <6; L++)
            {
                temp1 = gsl_vector_complex_get(B, L);
                temp2 = gsl_matrix_complex_get(VECS[IK], L, K);
                temp1 = gsl_complex_sub(temp1, gsl_complex_mul(C, temp2));
                gsl_vector_complex_set(B, L, temp1);
            }
        }

        //B = B/||B||
        gslc_vector_complex_normalize(B);

        //VECS(*,II, IJ)  = B
        for(int L = 0; L <6; L++)
        {
            temp1 = gsl_vector_complex_get(B, L);
            gsl_matrix_complex_set(VECS[IK], L, II, temp1);
        }


    }
}

/**
 *  \brief Diagonalization of the central part of VECS(*,*,0), wich at this point, has no unstable component under M = DAT[M]*DAT[M-1]*...*DAT[1].
 *
 *    This central part is denoted (v1 v2 v3 v4), with (v1 v2) associated to the xy center, (v3 v4) to the z center.
 *    Providing that (vi, vj) actually spans the 2th dimension subspace (either xy or z center), we have that:
 *
 *        vi = k1 u1 + k2 u1b
 *        vj = k3 u1 + k4 u1b
 *
 *        where (u1, u1b) is the couple of eigenvectors associated with either the xy or z center, with eigenvalues (l1, l1b).
 *
 *        Thus, the routines inverses following code inverses the following system:
 *
 *
 *             |  vi  |   |   I    I     0    0   |  |k1 u1  |       |k1 u1  |
 *             |  vj  |   |   0    0     I    I   |  |k2 u1b |       |k2 u1b |
 *             | M*vi | = | l1*I l1b*I   0    0   |  |k3 u1  | = M64 |k3 u1  |
 *             | M*vj |   |   0    0   l1*I l1b*I |  |k4 u1b |`      |k4 u1b |
 *
 *        with M the monodromy matrix, in order to get u1 and u1b.
 *
 **/
void center_diag(gsl_matrix_complex **DAT, gsl_matrix_complex *VECS_0, gsl_vector_complex *DECSl, gsl_matrix_complex *DECSv, int M)
{
    //Computing the matrix M64 and its inverse M64i that appears in the routine center_diag
    void setM64(gsl_matrix_complex *M64, gsl_matrix_complex *M64i, gsl_complex l1);
    //====================================================================================
    // GSL objects
    //====================================================================================
    gsl_complex one_c  = gslc_complex(1.0, 0.0);
    gsl_complex zero_c = gslc_complex(0.0, 0.0);
    gsl_vector_complex_view VECS_V;
    gsl_vector_complex *v1 = gsl_vector_complex_calloc(6);
    gsl_vector_complex *v2 = gsl_vector_complex_calloc(6);
    gsl_vector_complex *vu = gsl_vector_complex_calloc(6);
    gsl_vector_complex *vs = gsl_vector_complex_calloc(6);
    gsl_vector_complex *h1 = gsl_vector_complex_calloc(6);
    gsl_vector_complex *h2 = gsl_vector_complex_calloc(6);
    int d1, d2;
    gsl_complex l1 , l1b, factor;
    int isItZCenter;

    //Result vector
    gsl_vector_complex *V64     = gsl_vector_complex_calloc(4*6);
    gsl_vector_complex_view V1  = gsl_vector_complex_subvector (V64 , 0  , 6);
    gsl_vector_complex_view V2  = gsl_vector_complex_subvector (V64 , 6  , 6);
    gsl_vector_complex_view V3  = gsl_vector_complex_subvector (V64 , 12 , 6);
    gsl_vector_complex_view V4  = gsl_vector_complex_subvector (V64 , 18 , 6);
    gsl_vector_complex *H64     = gsl_vector_complex_calloc(4*6);
    gsl_vector_complex_view H1  = gsl_vector_complex_subvector (H64 , 0  , 6);
    gsl_vector_complex_view H2  = gsl_vector_complex_subvector (H64 , 6  , 6);
    gsl_vector_complex_view H3  = gsl_vector_complex_subvector (H64 , 12 , 6);
    gsl_vector_complex_view H4  = gsl_vector_complex_subvector (H64 , 18 , 6);
    //Result matrix
    gsl_matrix_complex *Mdiag   = gsl_matrix_complex_calloc(6, 6);
    gsl_matrix_complex *Ldiag   = gsl_matrix_complex_calloc(6, 6);
    //M64 and M64i
    gsl_matrix_complex *M64   = gsl_matrix_complex_calloc(4*6, 4*6);  //matrix M64
    gsl_matrix_complex *M64i  = gsl_matrix_complex_calloc(4*6, 4*6);  //its inverse
    gsl_matrix_complex *M64M64i = gsl_matrix_complex_calloc(4*6, 4*6);
    gsl_matrix_complex *Id64 = gsl_matrix_complex_calloc(4*6, 4*6);

    //------------------------------------------------------------------------------------------------------------------------
    // Set the hyperbolic directions in Mdiag/Ldiag
    //------------------------------------------------------------------------------------------------------------------------
    //vu = unstable direction, vs = stable direction
    gslc_matrix_complex_column(vu, VECS_0, 0);
    gslc_matrix_complex_column(vs, VECS_0, 1);
    gsl_complex dotu, dots;

    //------------------------------------------------------------------------------------------------------------------------
    //Selecting the right eigenvalues in DECSl
    //We know that there are two couples of central eigenvalues, one on the xy plane, the other on the z plane.
    //------------------------------------------------------------------------------------------------------------------------
    int i = 0;
    cout << Csts::SSEPR << endl;
    cout << "center_diag. Getting the central eigenvectors. " << endl;
    for(i = 0; i < 4; i+=2)
    {
        //Selecting the eigenvalue
        l1     = gsl_vector_complex_get(DECSl, i);  //DECSl[i]
        l1b    = gsl_complex_conjugate(l1);         //conjugate of l1
        factor = gsl_complex_sub(l1, l1b);          //factor = l1 - l1b
        //Checking wich direction was selected
        VECS_V = gsl_matrix_complex_column (DECSv, i);
        //If there is no z component, we are dealing with the xy center
        if(i == 0 /*gsl_complex_abs(gsl_vector_complex_get(&VECS_V.vector, 2)) == 0*/)
        {
            d1 = 3;
            d2 = 4;
            isItZCenter = 0;
        }
        else
        {
            d1 = 2;
            d2 = 5;
            isItZCenter = 1;
        }

        //-----------------------------
        // Selecting the central vectors v1, v2
        //-----------------------------
        //v1 = VECS(*,d1,0)
        gslc_matrix_complex_column(v1, VECS_0, d1);
        //v2 = VECS(*,d2,0)
        gslc_matrix_complex_column(v2, VECS_0, d2);

        //-----------------------------
        // Computing hi = A*vi
        //-----------------------------
        //h1 = DAT[N]*...*DAT[1]*v1
        gslc_matrix_vector_product(DAT, v1 , h1, M);
        // OR h1 = VECS(*,d1,M)
        //VECS_V = gsl_matrix_complex_column (VECS[M], d1);
        //gsl_vector_complex_memcpy(h1, &VECS_V .vector);

        //h2 = DAT[N]*...*DAT[1]*v2
        gslc_matrix_vector_product(DAT, v2 , h2, M);
        // OR h2 = VECS(*,d2,M)
        //VECS_V = gsl_matrix_complex_column (VECS[M], d2);
        //gsl_vector_complex_memcpy(h2, &VECS_V .vector);

        //cout << "diacp. h1 = " << endl;
        //gslc_vector_complex_printf(h1);
        //cout << "diacp. h2 = " << endl;
        //gslc_vector_complex_printf(h2);

        //get rid of vu projection: v1 = v1 - <v1,vu>vu
        gsl_blas_zdotu(v1, vu, &dotu);
        gsl_blas_zdotu(v1, vs, &dots);
        //cout << "center_diag. <v1,vu> = " << GSL_REAL(dotu) << "  " << GSL_IMAG(dotu) << endl;
        //cout << "center_diag. <v1,vs> = " << GSL_REAL(dots) << "  " << GSL_IMAG(dots) << endl;
        //gsl_blas_zaxpy(gsl_complex_negative(dots), vs, v1);
        //gsl_blas_zaxpy(gsl_complex_negative(dotu), vu, v1);

        //get rid of vu projection: v2 = v2 - <v2,vu>vu
        gsl_blas_zdotu(v2, vu, &dotu);
        gsl_blas_zdotu(v2, vs, &dots);
        //cout << "center_diag. <v2,vu> = " << GSL_REAL(dotu) << "  " << GSL_IMAG(dotu) << endl;
        //cout << "center_diag. <v2,vs> = " << GSL_REAL(dots) << "  " << GSL_IMAG(dots) << endl;
        //gsl_blas_zaxpy(gsl_complex_negative(dots), vs, v2);
        //gsl_blas_zaxpy(gsl_complex_negative(dotu), vu, v2);


        //------------------------------------------------------------------------------------------------------------------------
        // Providing that (v1, v2) actually spans the 2th dimension subspace (either xy or z center), we have that:
        //
        // v1 = k1 u1 + k2 u1b
        // v2 = k3 u1 + k4 u1b
        //
        // where (u1, u1b) is the couple of eigenvectors associated with either the xy or z center, with eigenvalues (l1, l1b).
        //
        // The following code inverses the following system:
        //
        //
        //      | v1 |   |   I    I     0    0   |  |k1 u1  |       |k1 u1  |
        //      | v2 |   |   0    0     I    I   |  |k2 u1b |       |k2 u1b |
        //      | h1 | = | l1*I l1b*I   0    0   |  |k3 u1  | = M64 |k3 u1  |
        //      | h2 |   |   0    0   l1*I l1b*I |  |k4 u1b |`      |k4 u1b |
        //
        //------------------------------------------------------------------------------------------------------------------------
        //V64 = [v1 v2 h1 h2]^T
        gsl_vector_complex_memcpy(&V1.vector, v1);
        gsl_vector_complex_memcpy(&V2.vector, v2);
        gsl_vector_complex_memcpy(&V3.vector, h1);
        gsl_vector_complex_memcpy(&V4.vector, h2);

        //M64 and M64i
        setM64(M64, M64i, l1);

        //M4M4i = M4*M4i
        gsl_matrix_complex_set_identity(Id64);
        gsl_blas_zgemm (CblasNoTrans , CblasNoTrans , one_c , M64 , M64i , zero_c , M64M64i );
        gsl_matrix_complex_sub(Id64, M64M64i);


        if(isItZCenter)
            cout << "Z  center. |Id - M64*M64i|_inf = " << gslc_matrix_complex_infinity_norm(Id64) << endl;
        else
            cout << "XY center. |Id - M64*M64i|_inf = " << gslc_matrix_complex_infinity_norm(Id64) << endl;


        //------------------------------------------------------------------------------------
        //Inversing the system
        //------------------------------------------------------------------------------------
        //H64 = M64i*V64
        gsl_blas_zgemv (CblasNoTrans, one_c , M64i, V64 , zero_c , H64);
        //V64 =  H64
        gsl_vector_complex_memcpy(&V1.vector, &H1.vector);
        gsl_vector_complex_memcpy(&V2.vector, &H2.vector);
        gsl_vector_complex_memcpy(&V3.vector, &H3.vector);
        gsl_vector_complex_memcpy(&V4.vector, &H4.vector);
        //Normalize the results
        gslc_vector_complex_normalize(&V1.vector);
        gslc_vector_complex_normalize(&V2.vector);
        gslc_vector_complex_normalize(&V3.vector);
        gslc_vector_complex_normalize(&V4.vector);
        //Set results in Mdiag
        gslc_matrix_complex_column_V(Mdiag, &V1.vector, d1);
        gslc_matrix_complex_column_V(Mdiag, &V2.vector, d2);
        //Set the results in Ldiag
        gsl_matrix_complex_set(Ldiag, d1, d1, l1);
        gsl_matrix_complex_set(Ldiag, d2, d2, l1b);
    }

    //Test of Mdiag/Ldiag
    cout << "center_diag. Test of the results. " << endl;
    eigensystem_test(Ldiag, Mdiag, DAT, M);
    gslc_matrix_complex_fprintf_pretty(Mdiag, (char*) "fprint/Mdiag.txt");
}

/**
 *  \brief Computing the matrix M64 and its inverse M64i that appears in the routine center_diag
 *
 * M64 is defined by:
 *
 *             |  vi  |   |   I    I     0    0   |  |k1 u1  |       |k1 u1  |
 *             |  vj  |   |   0    0     I    I   |  |k2 u1b |       |k2 u1b |
 *             | M*vi | = | l1*I l1b*I   0    0   |  |k3 u1  | = M64 |k3 u1  |
 *             | M*vj |   |   0    0   l1*I l1b*I |  |k4 u1b |`      |k4 u1b |
 **/
void setM64(gsl_matrix_complex *M64, gsl_matrix_complex *M64i, gsl_complex l1)
{
    //------------------------------------------------------------------------------------
    //Init
    //------------------------------------------------------------------------------------
    gsl_complex l1b    = gsl_complex_conjugate(l1);         //conjugate of l1
    gsl_complex factor = gsl_complex_sub(l1, l1b);          //factor = l1 - l1b
    gsl_complex one_c  = gslc_complex(1.0, 0.0);
    gsl_matrix_complex_set_zero(M64);
    gsl_matrix_complex_set_zero(M64i);


    //------------------------------------------------------------------------------------
    //Views of M64
    //------------------------------------------------------------------------------------
    gsl_matrix_complex_view S11 = gsl_matrix_complex_submatrix (M64 , 0 , 0  ,  6 , 6 );
    gsl_matrix_complex_view S12 = gsl_matrix_complex_submatrix (M64,  0 , 6  ,  6 , 6 );
    gsl_matrix_complex_view S13 = gsl_matrix_complex_submatrix (M64,  0 , 12 ,  6 , 6 );
    gsl_matrix_complex_view S14 = gsl_matrix_complex_submatrix (M64,  0 , 18 ,  6 , 6 );
    gsl_matrix_complex_view S21 = gsl_matrix_complex_submatrix (M64 , 6 , 0  ,  6 , 6 );
    gsl_matrix_complex_view S22 = gsl_matrix_complex_submatrix (M64,  6 , 6  ,  6 , 6 );
    gsl_matrix_complex_view S23 = gsl_matrix_complex_submatrix (M64,  6 , 12 ,  6 , 6 );
    gsl_matrix_complex_view S24 = gsl_matrix_complex_submatrix (M64,  6 , 18 ,  6 , 6 );
    gsl_matrix_complex_view S31 = gsl_matrix_complex_submatrix (M64 , 12 , 0  , 6 , 6 );
    gsl_matrix_complex_view S32 = gsl_matrix_complex_submatrix (M64,  12 , 6  , 6 , 6 );
    gsl_matrix_complex_view S33 = gsl_matrix_complex_submatrix (M64,  12 , 12 , 6 , 6 );
    gsl_matrix_complex_view S34 = gsl_matrix_complex_submatrix (M64,  12 , 18 , 6 , 6 );
    gsl_matrix_complex_view S41 = gsl_matrix_complex_submatrix (M64 , 18 , 0  , 6 , 6 );
    gsl_matrix_complex_view S42 = gsl_matrix_complex_submatrix (M64,  18 , 6  , 6 , 6 );
    gsl_matrix_complex_view S43 = gsl_matrix_complex_submatrix (M64,  18 , 12 , 6 , 6 );
    gsl_matrix_complex_view S44 = gsl_matrix_complex_submatrix (M64,  18 , 18 , 6 , 6 );


    gsl_matrix_complex_set_identity(&S11.matrix);  //S11 = Id
    gsl_matrix_complex_set_identity(&S12.matrix);  //S12 = Id
    gsl_matrix_complex_set_identity(&S23.matrix);  //S23 = Id
    gsl_matrix_complex_set_identity(&S24.matrix);  //S24 = Id


    gsl_matrix_complex_set_identity(&S31.matrix);  //S31 = Id
    gsl_matrix_complex_set_identity(&S32.matrix);  //S32 = Id
    gsl_matrix_complex_set_identity(&S43.matrix);  //S43 = Id
    gsl_matrix_complex_set_identity(&S44.matrix);  //S44 = Id

    gsl_matrix_complex_scale(&S31.matrix, l1);     //S31 = l1*Id
    gsl_matrix_complex_scale(&S32.matrix, l1b);    //S32 = l1b*Id
    gsl_matrix_complex_scale(&S43.matrix, l1);     //S43 = l1*Id
    gsl_matrix_complex_scale(&S44.matrix, l1b);    //S44 = l1b*Id

    //------------------------------------------------------------------------------------
    //Views of M64i
    //------------------------------------------------------------------------------------
    S11 = gsl_matrix_complex_submatrix (M64i , 0 , 0  , 6 , 6 );
    S12 = gsl_matrix_complex_submatrix (M64i,  0 , 6  , 6 , 6 );
    S13 = gsl_matrix_complex_submatrix (M64i,  0 , 12 , 6 , 6 );
    S14 = gsl_matrix_complex_submatrix (M64i,  0 , 18 , 6 , 6 );
    S21 = gsl_matrix_complex_submatrix (M64i , 6 , 0  , 6 , 6 );
    S22 = gsl_matrix_complex_submatrix (M64i,  6 , 6  , 6 , 6 );
    S23 = gsl_matrix_complex_submatrix (M64i,  6 , 12 , 6 , 6 );
    S24 = gsl_matrix_complex_submatrix (M64i,  6 , 18 , 6 , 6 );
    S31 = gsl_matrix_complex_submatrix (M64i , 12 , 0  , 6 , 6 );
    S32 = gsl_matrix_complex_submatrix (M64i,  12 , 6  , 6 , 6 );
    S33 = gsl_matrix_complex_submatrix (M64i,  12 , 12 , 6 , 6 );
    S34 = gsl_matrix_complex_submatrix (M64i,  12 , 18 , 6 , 6 );
    S41 = gsl_matrix_complex_submatrix (M64i , 18 , 0  , 6 , 6 );
    S42 = gsl_matrix_complex_submatrix (M64i,  18 , 6  , 6 , 6 );
    S43 = gsl_matrix_complex_submatrix (M64i,  18 , 12 , 6 , 6 );
    S44 = gsl_matrix_complex_submatrix (M64i,  18 , 18 , 6 , 6 );


    gsl_matrix_complex_set_identity(&S11.matrix);                        //S11 = Id
    gsl_matrix_complex_set_identity(&S13.matrix);                        //S13 = Id
    gsl_matrix_complex_set_identity(&S21.matrix);                        //S21 = Id
    gsl_matrix_complex_set_identity(&S23.matrix);                        //S23 = Id
    gsl_matrix_complex_set_identity(&S32.matrix);                        //S32 = Id
    gsl_matrix_complex_set_identity(&S34.matrix);                        //S34 = Id
    gsl_matrix_complex_set_identity(&S42.matrix);                        //S42 = Id
    gsl_matrix_complex_set_identity(&S44.matrix);                        //S44 = Id

    gsl_matrix_complex_scale(&S11.matrix, gsl_complex_negative(l1b));    //S11 = -l1b*Id
    gsl_matrix_complex_scale(&S21.matrix, l1);                           //S21 =  l1*Id
    gsl_matrix_complex_scale(&S23.matrix, gsl_complex_negative(one_c));  //S23 = -Id
    gsl_matrix_complex_scale(&S32.matrix, gsl_complex_negative(l1b));    //S32 = -l1b*Id
    gsl_matrix_complex_scale(&S42.matrix, l1);                           //S42 =  l1*Id
    gsl_matrix_complex_scale(&S44.matrix, gsl_complex_negative(one_c));  //S44 = -Id

    //M64i*=1/(l1-l1b)
    gsl_matrix_complex_scale(M64i, gsl_complex_inverse(factor));
}

//----------------------------------------------------------------------------------------
// Matrix integration
//----------------------------------------------------------------------------------------
/**
 *  \brief Obtain the FFT decomposition of various matrices P, Q = inv(P) FT11, FT12, FT21 and FT22.
 *
 *  We recall that the change of coordinates compute by nfo2 is of the form:
 *
 *  z = P(t)C zt + V(t)  <==> zt = Q(t) (z - V(t))  and  Q = inv(P)
 *
 *  The following decompositions are computed:
 *  - P
 *  - Q = inv(P)
 *  - FT11, FT12, FT21, FT22 that contains the c.o.c. for y_{theta}
 *  - G1, that contains V(t) = (g1, g2, 0, g3, g4, 0)^T (see code for details)
 *  - Xe, Xm, Xs, the translated positions of the Earth, Moon and Sun so that Xe[i] = real Xe[i] - V[i].
 *  See the modified potential expansion in pdf for details. Note that:
 *  - Xe[0] contains the x position of the Earth
 *  - Xe[1] contains the y position of the Earth
 *  - Xe[2] does not contains the z position of the Earth (which is always zero) but contains 1/sqrt(Xe[0]^2 + Xe[1]^2)
 * i.e. inverse of the euclidian distance from the origin, given as Fourier series.
 *  - Xe[3:5] contains the UNtranslated positions of the Earth (xe, ye, ze).
 *
 *
 * Note: G1 contains the gi coefficients with the form:
 *
 * G1 = | g1  g2 |
 *      | g3  g4 |
 *
 *  The idea behind this form is to be able to use the routines implemented for the matrices
 *
 *  The results are stored in txt files.
 *  Each Fourer expansion is tested on a different grid than the one used for the FFT process
 *  For now, the test of the periodicity of P has been abandonned
 *  WARNING: to compute Q = inv(P), the following hypothesis is made (consistent with what mono_eigen produces):
 *
 *     | p11  p12   0   p14  p15    0  |
 *     | p21  p22   0   p24  p25    0  |
 * P = |  0    0   p33   0    0    p36 |
 *     | p41  p42   0   p44  p45    0  |
 *     | p51  p52   0   p54  p55    0  |
 *     |  0    0   p63   0    0    p66 |
 *
 *
 *  Note that, as much as possible, the coefficients of P, Q, G1, Xe, Xs, and Xm
 *  are computed on half a period rather than a full period of the QBCP (see nfo2_int).
 *
 */
void nfo2_coc(gsl_odeiv2_driver *d,const double y0[], const double n,
             FBPL *fbpl, gsl_matrix* Br, gsl_matrix* R, gsl_matrix* JB,
             int N, int is_stored)
{
    cout << "-----------------------------------------------" << endl;
    cout << "               nfo2_coc                         " << endl;
    cout << "-----------------------------------------------" << endl;
    //------------------------------------------------------------------------------------
    // Initialization
    //------------------------------------------------------------------------------------
    //The validity of the FFTs is tested on a fftPlot points grid (different from the integer N)
    int fftPlot = 1000;
    //The FFT are computed up to order OFS_ORDER
    int fftN = max(OFS_ORDER,20);

    //Reset the driver
    gsl_odeiv2_driver_reset(d);

    //Allocation
    gsl_matrix** P    = gslc_matrix_array_alloc(6, 6, N);
    gsl_matrix** Pfb  = gslc_matrix_array_alloc(6, 6, N);
    gsl_matrix** Q    = gslc_matrix_array_alloc(6, 6, N);
    gsl_matrix** Qfb  = gslc_matrix_array_alloc(6, 6, N);
    gsl_matrix** FT11 = gslc_matrix_array_alloc(3, 3, N);
    gsl_matrix** FT12 = gslc_matrix_array_alloc(3, 3, N);
    gsl_matrix** FT21 = gslc_matrix_array_alloc(3, 3, N);
    gsl_matrix** FT22 = gslc_matrix_array_alloc(3, 3, N);
    gsl_matrix** G1   = gslc_matrix_array_alloc(2, 2, N);
    gsl_matrix** Xe   = gslc_matrix_array_alloc(6, 1, N);
    gsl_matrix** Xm   = gslc_matrix_array_alloc(6, 1, N);
    gsl_matrix** Xs   = gslc_matrix_array_alloc(6, 1, N);

    gsl_matrix** Pt    = gslc_matrix_array_alloc(6, 6, fftPlot);
    gsl_matrix** Pfbt  = gslc_matrix_array_alloc(6, 6, fftPlot);
    gsl_matrix** Qt    = gslc_matrix_array_alloc(6, 6, fftPlot);
    gsl_matrix** Qfbt  = gslc_matrix_array_alloc(6, 6, fftPlot);
    gsl_matrix** FT11b = gslc_matrix_array_alloc(3, 3, fftPlot);
    gsl_matrix** FT12b = gslc_matrix_array_alloc(3, 3, fftPlot);
    gsl_matrix** FT21b = gslc_matrix_array_alloc(3, 3, fftPlot);
    gsl_matrix** FT22b = gslc_matrix_array_alloc(3, 3, fftPlot);
    gsl_matrix** G1b   = gslc_matrix_array_alloc(2, 2, fftPlot);
    gsl_matrix** Xeb   = gslc_matrix_array_alloc(6, 1, fftPlot);
    gsl_matrix** Xmb   = gslc_matrix_array_alloc(6, 1, fftPlot);
    gsl_matrix** Xsb   = gslc_matrix_array_alloc(6, 1, fftPlot);

    //------------------------------------------------------------------------------------
    // Integration
    //------------------------------------------------------------------------------------
    //Integrate the matrices P, FT11, FT12, FT21 and FT22 on a N point grid
    nfo2_int(d, y0, n, fbpl, R, JB, N, P, Pfb, Q, Qfb, FT11, FT12, FT21, FT22, G1, Xe, Xm, Xs);
    //Integrate the matrices P, FT11, FT12, FT21 and FT22 on a fftPlot+1 point grid
    nfo2_int(d, y0, n, fbpl, R, JB, fftPlot, Pt, Pfbt, Qt, Qfbt, FT11b, FT12b, FT21b, FT22b, G1b, Xeb, Xmb, Xsb);


    //------------------------------------------------------------------------------------
    // Fourier Analysis of the coefficients of G1
    //------------------------------------------------------------------------------------
    string F_COC = fbpl->cs.F_COC;
    string filename = F_COC+"G1_";
    Ofsc xFFT(fftN);
    cout << "-----------------------------------------------" << endl;
    cout << "FFT comp & validity: G1                        " << endl;
    cout << "-----------------------------------------------" << endl;
    nfo2_fft(xFFT, G1, G1b, n,  N, fftN, fftPlot, filename, CLEANED_FFT, is_stored);

    //------------------------------------------------------------------------------------
    // Fourier Analysis of the coefficients of P
    //------------------------------------------------------------------------------------
    filename = F_COC+"P";
    cout << "-----------------------------------------------" << endl;
    cout << "FFT comp & validity: P                         " << endl;
    cout << "-----------------------------------------------" << endl;
    nfo2_fft(xFFT, P, Pt, n, N, fftN, fftPlot, filename, ORIGINAL_FFT, is_stored);

    //------------------------------------------------------------------------------------
    // Fourier Analysis of the coefficients of Pfb. Uncomment to replace P by Pfb
    //------------------------------------------------------------------------------------
    if(fbpl->model == Csts::QBCP) //does not work for BCP
    {
        cout << "-----------------------------------------------"   << endl;
        cout << "FFT comp & validity: Pfb                         " << endl;
        cout << "-----------------------------------------------"   << endl;
        nfo2_fft(xFFT, Pfb, Pfbt, n, N, fftN, fftPlot, filename, ORIGINAL_FFT, is_stored);
    }

    //------------------------------------------------------------------------------------
    // Fourier Analysis of the coefficients of Q
    //------------------------------------------------------------------------------------
    filename = F_COC+"Q";
    cout << "-----------------------------------------------" << endl;
    cout << "FFT comp & validity: Q                         " << endl;
    cout << "-----------------------------------------------" << endl;
    nfo2_fft(xFFT, Q, Qt, n, N, fftN, fftPlot, filename, ORIGINAL_FFT, is_stored);

    //------------------------------------------------------------------------------------
    // Fourier Analysis of the coefficients of Qfb. Uncomment to replace Q by Qfb
    //------------------------------------------------------------------------------------
    if(fbpl->model == Csts::QBCP) //does not work for BCP
    {
        cout << "-----------------------------------------------" << endl;
        cout << "FFT comp & validity: Qfb                         " << endl;
        cout << "-----------------------------------------------" << endl;
        nfo2_fft(xFFT, Qfb, Qfbt, n, N, fftN, fftPlot, filename, ORIGINAL_FFT, is_stored);
    }

    //------------------------------------------------------------------------------------
    // Fourier Analysis of the coefficients of Xe, Xm and Xs
    //------------------------------------------------------------------------------------
    filename = F_COC+"Xe";
    cout << "-----------------------------------------------" << endl;
    cout << "FFT comp & validity: Xe                        " << endl;
    cout << "-----------------------------------------------" << endl;
    //nfo2_fft(xFFT, Xe, Xeb,  N, fftN, fftPlot, filename, CLEANED_FFT, is_stored);
    nfo2_fft(xFFT, Xe, Xeb, n, N, fftN, fftPlot, filename, ORIGINAL_FFT, is_stored);

    filename = F_COC+"Xm";
    cout << "-----------------------------------------------" << endl;
    cout << "FFT comp & validity: Xm                        " << endl;
    cout << "-----------------------------------------------" << endl;
    //nfo2_fft(xFFT, Xm, Xmb, N, fftN, fftPlot, filename, CLEANED_FFT, is_stored);
    nfo2_fft(xFFT, Xm, Xmb, n, N, fftN, fftPlot, filename, ORIGINAL_FFT, is_stored);

    filename = F_COC+"Xs";
    cout << "-----------------------------------------------" << endl;
    cout << "FFT comp & validity: Xs                        " << endl;
    cout << "-----------------------------------------------" << endl;
    //nfo2_fft(xFFT, Xs, Xsb,  N, fftN, fftPlot, filename, CLEANED_FFT, is_stored);
    nfo2_fft(xFFT, Xs, Xsb, n, N, fftN, fftPlot, filename, ORIGINAL_FFT, is_stored);

    //------------------------------------------------------------------------------------
    // Optionnal tests
    //------------------------------------------------------------------------------------
    //Periodicity test on P
    //nfo2_periodicity_test(d, y0, R);
    //Symmetry test on P
    //nfo2_symmetry_test(d, y0, R, 10, 100);
    //Symplectic test on P
    symplecticity_test_for_P(d, y0, 1.0, R);

    //------------------------------------------------------------------------------------
    // Fourier Analysis of the coefficients of the Fijs
    //------------------------------------------------------------------------------------
    cout << "-----------------------------------------------" << endl;
    cout << "FFT comp & validity: FT11                      " << endl;
    cout << "-----------------------------------------------" << endl;
    filename = F_COC+"FT11_";
    nfo2_fft(xFFT, FT11, FT11b, n, N, fftN, fftPlot, filename, ORIGINAL_FFT, is_stored);
    cout << "-----------------------------------------------" << endl;
    cout << "FFT comp & validity: FT12                      " << endl;
    cout << "-----------------------------------------------" << endl;
    filename = F_COC+"FT12_";
    nfo2_fft(xFFT, FT12, FT12b, n, N, fftN, fftPlot, filename, ORIGINAL_FFT, is_stored);
    cout << "-----------------------------------------------" << endl;
    cout << "FFT comp & validity: FT21                      " << endl;
    cout << "-----------------------------------------------" << endl;
    filename = F_COC+"FT21_";
    nfo2_fft(xFFT, FT21, FT21b, n, N, fftN, fftPlot, filename, ORIGINAL_FFT, is_stored);
    cout << "-----------------------------------------------" << endl;
    cout << "FFT comp & validity: FT22                      " << endl;
    cout << "-----------------------------------------------" << endl;
    filename = F_COC+"FT22_";
    nfo2_fft(xFFT, FT22, FT22b, n, N, fftN, fftPlot, filename, ORIGINAL_FFT, is_stored);



    gslc_matrix_array_free(P, N);
    gslc_matrix_array_free(FT11, N);
    gslc_matrix_array_free(FT12, N);
    gslc_matrix_array_free(FT21, N);
    gslc_matrix_array_free(FT22, N);

    gslc_matrix_array_free(Pt, fftPlot);
    gslc_matrix_array_free(FT11b, fftPlot);
    gslc_matrix_array_free(FT12b, fftPlot);
    gslc_matrix_array_free(FT21b, fftPlot);
    gslc_matrix_array_free(FT22b, fftPlot);
}


/**
 *  \brief Integrate the various matrices referenced on a N point grid to seed FFT process.
 *
 *  Same remarks as the routine nfo2_coc.
 */
void nfo2_int(gsl_odeiv2_driver *d,
             const double y0[],
             const double n,
             FBPL *fbpl,
             gsl_matrix* R,
             gsl_matrix* JB,
             int N,
             gsl_matrix** P,                //6*6   |
             gsl_matrix** Pfb,              //6*6   |
             gsl_matrix** Q,                //6*6   |
             gsl_matrix** Qfb,              //6*6   |
             gsl_matrix** FT11,             //3*3   |
             gsl_matrix** FT12,             //3*3   | Outputs
             gsl_matrix** FT21,             //3*3   |
             gsl_matrix** FT22,             //3*3   |
             gsl_matrix** G1  ,             //2*2   |
             gsl_matrix** Xe  ,             //2*1   |
             gsl_matrix** Xm  ,             //2*1   |
             gsl_matrix** Xs  )             //2*1   |
{
    //------------------------------------------------------------------------------------
    // Initialization
    //------------------------------------------------------------------------------------
    double pe[3], ps[3], pm[3];         //Primaries position given from qbtbp
    double rc;                          //Euclidian distance from a given primary

    //Allocation
    gsl_matrix* Pb   = gsl_matrix_calloc (6, 6);
    gsl_matrix* Pc   = gsl_matrix_calloc (6, 6);
    gsl_permutation * p6 = gsl_permutation_alloc (6);

    gsl_matrix* Q1 = gsl_matrix_calloc (3, 3);
    gsl_matrix* Q2 = gsl_matrix_calloc (3, 3);
    gsl_matrix* Q3 = gsl_matrix_calloc (3, 3);

    gsl_matrix* AUX = gsl_matrix_calloc (3, 3);
    gsl_matrix* BUX = gsl_matrix_calloc (3, 3);
    gsl_matrix* CUX = gsl_matrix_calloc (3, 3);
    gsl_matrix* DUX = gsl_matrix_calloc (6, 6);

    //------------------------------------------------------------------------------------
    //Matrix views
    //------------------------------------------------------------------------------------
    gsl_matrix_view P11;
    gsl_matrix_view P12;
    gsl_matrix_view P21;
    gsl_matrix_view P22;

    gsl_matrix_view JB11 = gsl_matrix_submatrix (JB , 0 , 0 , 3 , 3 );   //JB11
    gsl_matrix_view JB12 = gsl_matrix_submatrix (JB , 0 , 3 , 3 , 3 );   //JB12
    gsl_matrix_view JB21 = gsl_matrix_submatrix (JB , 3 , 0 , 3 , 3 );   //JB21
    gsl_matrix_view JB22 = gsl_matrix_submatrix (JB , 3 , 3 , 3 , 3 );   //JB22

    //Update the starting point (including P(0) = Id)
    double y[42];
    for(int i=0; i<42; i++) y[i] = y0[i];

    //------------------------------------------------------------------------------------
    //Integration Loop for P except its 5th column (stable direction) and G1
    //------------------------------------------------------------------------------------
    double ti = 0.0;
    double t  = 0.0;
    double t1 = 2*M_PI/n;


    //Reset the driver
    gsl_odeiv2_driver_reset(d);
    int s;
    for(int i = 0; i < N; i++)
    {
        //------------------------------------
        //Integration step
        //------------------------------------
        if(i>0)
        {
            ti = i * t1 / N;
            gsl_odeiv2_driver_apply (d, &t , ti , y);
        }

        //------------------------------------
        //G1 update
        //------------------------------------
        gsl_matrix_set(G1[i], 0, 0, y[0]); //G1[0][0] = y[0] = x
        gsl_matrix_set(G1[i], 0, 1, y[1]); //G1[0][1] = y[1] = y
        gsl_matrix_set(G1[i], 1, 0, y[3]); //G1[1][0] = y[3] = px
        gsl_matrix_set(G1[i], 1, 1, y[4]); //G1[1][1] = y[4] = py

        //------------------------------------
        //Translated primaries positions update
        //------------------------------------
        eval_array_coef(ps, t, fbpl->us.n, fbpl->n_order_fourier, fbpl->cs.ps, 3);
        eval_array_coef(pe, t, fbpl->us.n, fbpl->n_order_fourier, fbpl->cs.pe, 3);
        eval_array_coef(pm, t, fbpl->us.n, fbpl->n_order_fourier, fbpl->cs.pm, 3);

        //Xe update
        gsl_matrix_set(Xe[i], 0, 0, pe[0] - y[0]); //tilde(xe) of the Earth
        gsl_matrix_set(Xe[i], 1, 0, pe[1] - y[1]); //tilde(ye) of the Earth
        //Last indix contains the term 1/le = 1/sqrt(tilde(xe)^2+tilde(ye)^2)
        rc = sqrt((pe[0] - y[0])*(pe[0] - y[0]) + (pe[1] - y[1])*(pe[1] - y[1]) );
        gsl_matrix_set(Xe[i], 2, 0, 1.0/rc );

        //Xm update
        gsl_matrix_set(Xm[i], 0, 0, pm[0] - y[0]); //tilde(xm) of the Moon
        gsl_matrix_set(Xm[i], 1, 0, pm[1] - y[1]); //tilde(ym) of the Moon
        //Last indix contains the term 1/lm = 1/sqrt(tilde(xm)^2+tilde(ym)^2)
        rc = sqrt((pm[0] - y[0])*(pm[0] - y[0]) + (pm[1] - y[1])*(pm[1] - y[1]) );
        gsl_matrix_set(Xm[i], 2, 0, 1.0/rc );

        //Xs update
        gsl_matrix_set(Xs[i], 0, 0, ps[0] - y[0]); //tilde(xs) of the Sun
        gsl_matrix_set(Xs[i], 1, 0, ps[1] - y[1]); //tilde(ys) of the Sun
        //Last indix contains the term 1/lm = 1/sqrt(tilde(xm)^2+tilde(ym)^2)
        rc = sqrt((ps[0] - y[0])*(ps[0] - y[0]) + (ps[1] - y[1])*(ps[1] - y[1]) );
        gsl_matrix_set(Xs[i], 2, 0, 1.0/rc );

        //------------------------------------
        //Non-translated primaries positions update
        //------------------------------------
        //Xe update
        gsl_matrix_set(Xe[i], 3, 0, pe[0]); //xe of the Earth
        gsl_matrix_set(Xe[i], 4, 0, pe[1]); //ye of the Earth
        //Last indix contains the term 1/le = 1/sqrt(xe^2+ye^2)
        rc = sqrt((pe[0])*(pe[0]) + (pe[1])*(pe[1]) );
        gsl_matrix_set(Xe[i], 5, 0, 1.0/rc );

        //Xm update
        gsl_matrix_set(Xm[i], 3, 0, pm[0]); //xm of the Moon
        gsl_matrix_set(Xm[i], 4, 0, pm[1]); //ym of the Moon
        //Last indix contains the term 1/lm = 1/sqrt(xm^2+ym^2)
        rc = sqrt((pm[0])*(pm[0]) + (pm[1])*(pm[1]) );
        gsl_matrix_set(Xm[i], 5, 0, 1.0/rc );

        //Xs update
        gsl_matrix_set(Xs[i], 3, 0, ps[0]); //xs of the Sun
        gsl_matrix_set(Xs[i], 4, 0, ps[1]); //ys of the Sun
        //Last indix contains the term 1/ls = 1/sqrt(xs^2+ys^2)
        rc = sqrt((ps[0])*(ps[0]) + (ps[1])*(ps[1]) );
        gsl_matrix_set(Xs[i], 5, 0, 1.0/rc );

        //------------------------------------
        //Pb update
        //------------------------------------
        gslc_vec_to_mat(Pb, y, 6, 6, 6);
        //P = Pb*R
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Pb, R, 0.0, P[i]);
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Pb, R, 0.0, Pfb[i]);

        //------------------------------------
        //Update Q1, Q2 and Q3
        //------------------------------------
        qbcp_vfn_Q(t, y, Q1, Q2, Q3, fbpl);
    }


    //------------------------------------------------------------------------------------
    // Integration Loop the 5th column of P (stable direction) and FT11, FT12, FT21, and FT22.
    // Note that the coefficients are computed backwards in time.
    // This is due to the specificity of the stable direction that has more precision when computing backwards.
    // The FTijs do not need backward computation but they need all the columns of P.
    //------------------------------------------------------------------------------------
    ti = 0.0;
    t  = 0.0;
    t1 = 2*M_PI/n;
    d->h = -1.0e-6;
    //Reset the driver
    gsl_odeiv2_driver_reset(d);
    //Reset the initial conditions
    for(int i=0; i<42; i++) y[i] = y0[i];
    for(int i = 1; i <= N; i++)
    {
        //Integration step
        ti = -i * t1 / N;
        gsl_odeiv2_driver_apply (d, &t , ti , y);

        //Pb update
        gslc_vec_to_mat(Pb, y, 6, 6, 6);
        //P = Pb*R for the 5th (indix 4) column only
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Pb, R, 0.0, DUX);
        for(int k = 0; k <6; k++) gsl_matrix_set(P[N-i], k, 4, gsl_matrix_get(DUX, k, 4));
        //for(int k = 0; k <6; k++) gsl_matrix_set(Pfb[N-i], k, 4, gsl_matrix_get(DUX, k, 4));
        //Note the other column have been computed in the previous loop,
        //so we can compute the FTijs without any danger

        //--------------------
        //Q = inv(P)
        //--------------------
        //Use of GSL library
        gsl_matrix_memcpy(Pc, P[N-i]);
        gsl_linalg_LU_decomp (Pc, p6, &s);
        gsl_linalg_LU_invert (Pc , p6 , Q[N-i]);

        //--------------------
        //Qfb = inv(Pfb)
        //--------------------
        //Use of GSL library
        gsl_matrix_memcpy(Pc, Pfb[N-i]);
        gsl_linalg_LU_decomp (Pc, p6, &s);
        gsl_linalg_LU_invert (Pc , p6 , Qfb[N-i]);

        //Update Q1, Q2 and Q3
        qbcp_vfn_Q(t, y, Q1, Q2, Q3, fbpl);

        //-------------------------------------------
        //Views
        //-------------------------------------------
        P11 = gsl_matrix_submatrix (P[N-i] , 0 , 0 , 3 , 3 );   //P11
        P12 = gsl_matrix_submatrix (P[N-i] , 0 , 3 , 3 , 3 );   //P12
        P21 = gsl_matrix_submatrix (P[N-i] , 3 , 0 , 3 , 3 );   //P21
        P22 = gsl_matrix_submatrix (P[N-i] , 3 , 3 , 3 , 3 );   //P22

        //-------------------------------------------
        //FT11[i]
        //-------------------------------------------
        //FT11[i] = -1/2*JB21
        gsl_matrix_memcpy(FT11[N-i], &JB21.matrix);
        gsl_matrix_scale(FT11[N-i], -0.5);

        //AUX = Q1*P11
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Q1, &P11.matrix, 0.0, AUX);
        //BUX = P11^T
        gsl_matrix_transpose_memcpy(BUX,  &P11.matrix);
        //CUX = BUX*AUX = P11^T*Q1*P11
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, BUX, AUX, 0.0, CUX);
        //FT11[i] = FT11[i] - CUX = -1/2*JB21-P11^T*Q1*P11
        gsl_matrix_sub(FT11[N-i], CUX);

        //AUX = Q2*P21
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Q2, &P21.matrix, 0.0, AUX);
        //BUX = P11^T
        gsl_matrix_transpose_memcpy(BUX,  &P11.matrix);
        //CUX = BUX*AUX = P11^T*Q2*P21
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, BUX, AUX, 0.0, CUX);
        //FT11[i] = FT11[i] - CUX = -1/2*JB21-P11^T*Q1*P11-P11^T*Q2*P21
        gsl_matrix_sub(FT11[N-i], CUX);

        //AUX = Q3*P21
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Q3, &P21.matrix, 0.0, AUX);
        //BUX = P21^T
        gsl_matrix_transpose_memcpy(BUX,  &P21.matrix);
        //CUX = BUX*AUX = P21^T*Q3*P21
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, BUX, AUX, 0.0, CUX);
        //FT11[i] = FT11[i] - CUX = -1/2*JB21-P11^T*Q1*P11-P11^T*Q2*P21-P21^T*Q3*P21
        gsl_matrix_sub(FT11[N-i], CUX);

        //FT11[i] = FT11[i]/n
        gsl_matrix_scale(FT11[N-i], 1.0/n);


        //-------------------------------------------
        //FT12[i]
        //-------------------------------------------
        //FT12[i] = 1/2*JB11
        gsl_matrix_memcpy(FT12[N-i], &JB11.matrix);
        gsl_matrix_scale(FT12[N-i], 0.5);

        //AUX = Q1*P11
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Q1, &P11.matrix, 0.0, AUX);
        //BUX = P12^T
        gsl_matrix_transpose_memcpy(BUX,  &P12.matrix);
        //CUX = BUX*AUX = P12^T*Q1*P11
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, BUX, AUX, 0.0, CUX);
        //FT12[i] = FT12[i] - CUX = 1/2*JB11-P12^T*Q1*P11
        gsl_matrix_sub(FT12[N-i], CUX);

        //AUX = Q2*P21
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Q2, &P21.matrix, 0.0, AUX);
        //BUX = P12^T
        gsl_matrix_transpose_memcpy(BUX,  &P12.matrix);
        //CUX = BUX*AUX = P12^T*Q2*P21
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, BUX, AUX, 0.0, CUX);
        //FT12[i] = FT12[i] - CUX = 1/2*JB11-P12^T*Q1*P11-P12^T*Q2*P21
        gsl_matrix_sub(FT12[N-i], CUX);

        //AUX = Q3*P21
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Q3, &P21.matrix, 0.0, AUX);
        //BUX = P22^T
        gsl_matrix_transpose_memcpy(BUX,  &P22.matrix);
        //CUX = BUX*AUX = P22^T*Q3*P21
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, BUX, AUX, 0.0, CUX);
        //FT12[i] = FT12[i] - CUX = 1/2*JB11-P12^T*Q1*P11-P22^T*Q3*P21
        gsl_matrix_sub(FT12[N-i], CUX);

        //FT12[i] = FT12[i]/n
        gsl_matrix_scale(FT12[N-i], 1.0/n);


        //-------------------------------------------
        //FT21[i]
        //-------------------------------------------
        //FT21[i] = -1/2*JB22
        gsl_matrix_memcpy(FT21[N-i], &JB22.matrix);
        gsl_matrix_scale(FT21[N-i], -0.5);

        //AUX = Q1*P12
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Q1, &P12.matrix, 0.0, AUX);
        //BUX = P11^T
        gsl_matrix_transpose_memcpy(BUX,  &P11.matrix);
        //CUX = BUX*AUX = P11^T*Q1*P12
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, BUX, AUX, 0.0, CUX);
        //FT21[i] = FT21[i] - CUX = -1/2*JB22-P11^T*Q1*P12
        gsl_matrix_sub(FT21[N-i], CUX);

        //AUX = Q2*P22
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Q2, &P22.matrix, 0.0, AUX);
        //BUX = P11^T
        gsl_matrix_transpose_memcpy(BUX,  &P11.matrix);
        //CUX = BUX*AUX = P11^T*Q2*P22
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, BUX, AUX, 0.0, CUX);
        //FT21[i] = FT21[i] - CUX = -1/2*JB22-P11^T*Q1*P12-P11^T*Q2*P22
        gsl_matrix_sub(FT21[N-i], CUX);

        //AUX = Q3*P22
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Q3, &P22.matrix, 0.0, AUX);
        //BUX = P21^T
        gsl_matrix_transpose_memcpy(BUX,  &P21.matrix);
        //CUX = BUX*AUX = P21^T*Q3*P22
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, BUX, AUX, 0.0, CUX);
        //FT21[i] = FT21[i] - CUX = -1/2*JB22-P11^T*Q1*P12-P11^T*Q2*P21-P21^T*Q3*P22
        gsl_matrix_sub(FT21[N-i], CUX);

        //FT21[i] = FT21[i]/n
        gsl_matrix_scale(FT21[N-i], 1.0/n);


        //-------------------------------------------
        //FT22[i]
        //-------------------------------------------
        //FT22[i] = 1/2*JB22
        gsl_matrix_memcpy(FT22[N-i], &JB12.matrix);
        gsl_matrix_scale(FT22[N-i], 0.5);

        //AUX = Q1*P12
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Q1, &P12.matrix, 0.0, AUX);
        //BUX = P12^T
        gsl_matrix_transpose_memcpy(BUX,  &P12.matrix);
        //CUX = BUX*AUX = P12^T*Q1*P12
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, BUX, AUX, 0.0, CUX);
        //FT22[i] = FT22[i] - CUX = 1/2*JB12-P12^T*Q1*P12
        gsl_matrix_sub(FT22[N-i], CUX);

        //AUX = Q2*P22
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Q2, &P22.matrix, 0.0, AUX);
        //BUX = P12^T
        gsl_matrix_transpose_memcpy(BUX,  &P12.matrix);
        //CUX = BUX*AUX = P12^T*Q2*P22
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, BUX, AUX, 0.0, CUX);
        //FT22[i] = FT22[i] - CUX = 1/2*JB12-P12^T*Q1*P12-P12^T*Q2*P22
        gsl_matrix_sub(FT22[N-i], CUX);

        //AUX = Q3*P22
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Q3, &P22.matrix, 0.0, AUX);
        //BUX = P22^T
        gsl_matrix_transpose_memcpy(BUX,  &P22.matrix);
        //CUX = BUX*AUX = P22^T*Q3*P22
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, BUX, AUX, 0.0, CUX);
        //FT22[i] = FT22[i] - CUX = 1/2*JB12-P12^T*Q1*P12-P12^T*Q2*P21-P22^T*Q3*P22
        gsl_matrix_sub(FT22[N-i], CUX);

        //FT22[i] = FT22[i]/n
        gsl_matrix_scale(FT22[N-i], 1.0/n);
    }
    d->h = 1.0e-6;


    //------------------------------------------------------------------------------------
    // Third loop, onward, on half period, to guarantee precision on G1 and Xs, Xe, Xm
    //------------------------------------------------------------------------------------
    ti = 0.0;
    t  = 0.0;
    t1 = 2*M_PI/n;
    d->h = 1.0e-6;
    //Reset the driver
    gsl_odeiv2_driver_reset(d);
    //Reset the initial conditions
    for(int i=0; i<42; i++) y[i] = y0[i];
    for(int i = 0; i <= N/2; i++)
    {
        //------------------------------------
        //Integration step
        //------------------------------------
        if(i>0)
        {
            ti = i * t1 / N;
            gsl_odeiv2_driver_apply (d, &t , ti , y);
        }

        //------------------------------------
        //G1 update
        //------------------------------------
        gsl_matrix_set(G1[i], 0, 0, y[0]); //G1[0][0] = y[0] = x
        gsl_matrix_set(G1[i], 0, 1, y[1]); //G1[0][1] = y[1] = y
        gsl_matrix_set(G1[i], 1, 0, y[3]); //G1[1][0] = y[3] = px
        gsl_matrix_set(G1[i], 1, 1, y[4]); //G1[1][1] = y[4] = py

        //------------------------------------
        //Translated primaries positions update
        //------------------------------------
        eval_array_coef(ps, t, fbpl->us.n, fbpl->n_order_fourier, fbpl->cs.ps, 3);
        eval_array_coef(pe, t, fbpl->us.n, fbpl->n_order_fourier, fbpl->cs.pe, 3);
        eval_array_coef(pm, t, fbpl->us.n, fbpl->n_order_fourier, fbpl->cs.pm, 3);

        //Xe update
        gsl_matrix_set(Xe[i], 0, 0, pe[0] - y[0]); //tilde(xe) of the Earth
        gsl_matrix_set(Xe[i], 1, 0, pe[1] - y[1]); //tilde(ye) of the Earth
        //Last indix contains the term 1/le = 1/sqrt(tilde(xe)^2+tilde(ye)^2)
        rc = sqrt((pe[0] - y[0])*(pe[0] - y[0]) + (pe[1] - y[1])*(pe[1] - y[1]) );
        gsl_matrix_set(Xe[i], 2, 0, 1.0/rc );

        //Xm update
        gsl_matrix_set(Xm[i], 0, 0, pm[0] - y[0]); //tilde(xe) of the Moon
        gsl_matrix_set(Xm[i], 1, 0, pm[1] - y[1]); //tilde(ye) of the Moon
        //Last indix contains the term 1/lm = 1/sqrt(tilde(xm)^2+tilde(ym)^2)
        rc = sqrt((pm[0] - y[0])*(pm[0] - y[0]) + (pm[1] - y[1])*(pm[1] - y[1]) );
        gsl_matrix_set(Xm[i], 2, 0, 1.0/rc );


        //Xs update
        gsl_matrix_set(Xs[i], 0, 0, ps[0] - y[0]); //tilde(xe) of the Moon
        gsl_matrix_set(Xs[i], 1, 0, ps[1] - y[1]); //tilde(ye) of the Moon
        //Last indix contains the term 1/lm = 1/sqrt(tilde(xm)^2+tilde(ym)^2)
        rc = sqrt((ps[0] - y[0])*(ps[0] - y[0]) + (ps[1] - y[1])*(ps[1] - y[1]) );
        gsl_matrix_set(Xs[i], 2, 0, 1.0/rc );


        //------------------------------------
        //Non-translated primaries positions update
        //------------------------------------
        //Xe update
        gsl_matrix_set(Xe[i], 3, 0, pe[0]); //xe of the Earth
        gsl_matrix_set(Xe[i], 4, 0, pe[1]); //ye of the Earth
        //Last indix contains the term 1/le = 1/sqrt(xe^2+tye^2)
        rc = sqrt((pe[0])*(pe[0]) + (pe[1])*(pe[1]) );
        gsl_matrix_set(Xe[i], 5, 0, 1.0/rc );

        //Xm update
        gsl_matrix_set(Xm[i], 3, 0, pm[0]); //tilde(xe) of the Moon
        gsl_matrix_set(Xm[i], 4, 0, pm[1]); //tilde(ye) of the Moon
        //Last indix contains the term 1/lm = 1/sqrt(tilde(xm)^2+tilde(ym)^2)
        rc = sqrt((pm[0])*(pm[0]) + (pm[1])*(pm[1]) );
        gsl_matrix_set(Xm[i], 5, 0, 1.0/rc );

        //Xs update
        gsl_matrix_set(Xs[i], 3, 0, ps[0]); //tilde(xe) of the Moon
        gsl_matrix_set(Xs[i], 4, 0, ps[1]); //tilde(ye) of the Moon
        //Last indix contains the term 1/lm = 1/sqrt(tilde(xm)^2+tilde(ym)^2)
        rc = sqrt((ps[0])*(ps[0]) + (ps[1])*(ps[1]) );
        gsl_matrix_set(Xs[i], 5, 0, 1.0/rc );

//        //------------------------------------
//        //Pfb update
//        //------------------------------------
//        //Pb update
//        gslc_vec_to_mat(Pb, y, 6, 6, 6);
//        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Pb, R, 0.0, Pfb[i]);
//
//        //--------------------
//        //Qfb = inv(Pfb)
//        //--------------------
//        //Use of GSL library
//        gsl_matrix_memcpy(Pc, Pfb[i]);
//        gsl_linalg_LU_decomp (Pc, p6, &s);
//        gsl_linalg_LU_invert (Pc , p6 , Qfb[i]);
    }

    //------------------------------------------------------------------------------------
    // Third loop, backward, on half period, to guarantee precision on G1 and Xs, Xe, Xm
    //------------------------------------------------------------------------------------
    ti = 2*M_PI/n;
    t  = 2*M_PI/n;
    t1 = 0.0;
    d->h = -1.0e-6;
    //Reset the driver
    gsl_odeiv2_driver_reset(d);
    //Reset the initial conditions
    for(int i=0; i<42; i++) y[i] = y0[i];
    for(int i=1; i < N/2; i++)
    {
        //------------------------------------
        //Integration step
        //------------------------------------
        if(i>0)
        {
            ti = 2*M_PI/n*(1.0-1.0*i/N);
            gsl_odeiv2_driver_apply (d, &t , ti , y);
        }

        //------------------------------------
        //G1 update
        //------------------------------------
        gsl_matrix_set(G1[N-i], 0, 0, y[0]); //G1[0][0] = y[0] = x
        gsl_matrix_set(G1[N-i], 0, 1, y[1]); //G1[0][1] = y[1] = y
        gsl_matrix_set(G1[N-i], 1, 0, y[3]); //G1[1][0] = y[3] = px
        gsl_matrix_set(G1[N-i], 1, 1, y[4]); //G1[1][1] = y[4] = py

        //------------------------------------
        //Translated primaries positions update
        //------------------------------------
        eval_array_coef(ps, t, fbpl->us.n, fbpl->n_order_fourier, fbpl->cs.ps, 3);
        eval_array_coef(pe, t, fbpl->us.n, fbpl->n_order_fourier, fbpl->cs.pe, 3);
        eval_array_coef(pm, t, fbpl->us.n, fbpl->n_order_fourier, fbpl->cs.pm, 3);

        //Xe update
        gsl_matrix_set(Xe[N-i], 0, 0, pe[0] - y[0]); //tilde(xe) of the Earth
        gsl_matrix_set(Xe[N-i], 1, 0, pe[1] - y[1]); //tilde(ye) of the Earth
        //Last indix contains the term 1/le = 1/sqrt(tilde(xe)^2+tilde(ye)^2)
        rc = sqrt((pe[0] - y[0])*(pe[0] - y[0]) + (pe[1] - y[1])*(pe[1] - y[1]) );
        gsl_matrix_set(Xe[N-i], 2, 0, 1.0/rc );

        //Xm update
        gsl_matrix_set(Xm[N-i], 0, 0, pm[0] - y[0]); //tilde(xe) of the Moon
        gsl_matrix_set(Xm[N-i], 1, 0, pm[1] - y[1]); //tilde(ye) of the Moon
        //Last indix contains the term 1/lm = 1/sqrt(tilde(xm)^2+tilde(ym)^2)
        rc = sqrt((pm[0] - y[0])*(pm[0] - y[0]) + (pm[1] - y[1])*(pm[1] - y[1]) );
        gsl_matrix_set(Xm[N-i], 2, 0, 1.0/rc );


        //Xs update
        gsl_matrix_set(Xs[N-i], 0, 0, ps[0] - y[0]); //tilde(xe) of the Moon
        gsl_matrix_set(Xs[N-i], 1, 0, ps[1] - y[1]); //tilde(ye) of the Moon
        //Last indix contains the term 1/lm = 1/sqrt(tilde(xm)^2+tilde(ym)^2)
        rc = sqrt((ps[0] - y[0])*(ps[0] - y[0]) + (ps[1] - y[1])*(ps[1] - y[1]) );
        gsl_matrix_set(Xs[N-i], 2, 0, 1.0/rc );


        //------------------------------------
        //Non-translated primaries positions update
        //------------------------------------
        //Xe update
        gsl_matrix_set(Xe[N-i], 3, 0, pe[0]); //xe of the Earth
        gsl_matrix_set(Xe[N-i], 4, 0, pe[1]); //ye of the Earth
        //Last indix contains the term 1/le = 1/sqrt(xe^2+tye^2)
        rc = sqrt((pe[0])*(pe[0]) + (pe[1])*(pe[1]) );
        gsl_matrix_set(Xe[N-i], 5, 0, 1.0/rc );

        //Xm update
        gsl_matrix_set(Xm[N-i], 3, 0, pm[0]); //tilde(xe) of the Moon
        gsl_matrix_set(Xm[N-i], 4, 0, pm[1]); //tilde(ye) of the Moon
        //Last indix contains the term 1/lm = 1/sqrt(tilde(xm)^2+tilde(ym)^2)
        rc = sqrt((pm[0])*(pm[0]) + (pm[1])*(pm[1]) );
        gsl_matrix_set(Xm[N-i], 5, 0, 1.0/rc );

        //Xs update
        gsl_matrix_set(Xs[N-i], 3, 0, ps[0]); //tilde(xe) of the Moon
        gsl_matrix_set(Xs[N-i], 4, 0, ps[1]); //tilde(ye) of the Moon
        //Last indix contains the term 1/lm = 1/sqrt(tilde(xm)^2+tilde(ym)^2)
        rc = sqrt((ps[0])*(ps[0]) + (ps[1])*(ps[1]) );
        gsl_matrix_set(Xs[N-i], 5, 0, 1.0/rc );

        //------------------------------------
        //Pb update
        //------------------------------------
        gslc_vec_to_mat(Pb, y, 6, 6, 6);
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Pb, R, 0.0, DUX);

        //------------------------------------
        //Pb update
        //------------------------------------
        gslc_vec_to_mat(Pb, y, 6, 6, 6);
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Pb, R, 0.0, Pfb[N-i]);

        //--------------------
        //Q = inv(P)
        //--------------------
        //Use of GSL library
        gsl_matrix_memcpy(Pc, Pfb[N-i]);
        gsl_linalg_LU_decomp (Pc, p6, &s);
        gsl_linalg_LU_invert (Pc , p6 , Qfb[N-i]);
    }


    d->h = 1.0e-6;

    //Memory release
    gsl_matrix_free(Pb);
    gsl_matrix_free(Pc);
    gsl_matrix_free(Q1);
    gsl_matrix_free(Q2);
    gsl_matrix_free(Q3);
    gsl_matrix_free(AUX);
    gsl_matrix_free(BUX);
    gsl_matrix_free(CUX);
    gsl_permutation_free(p6);
}





//----------------------------------------------------------------------------------------
// FFTs
//----------------------------------------------------------------------------------------
/**
 *  \brief GSL implementation of the FFT
 **/
void gslFFT_real(Ofsc &xFFT, int fftN, int N, gsl_vector *dEv, int parity)
{
    //--------------------------------------------------------------
    //FFT tools
    //--------------------------------------------------------------
    gsl_vector_complex *data_complex = gsl_vector_complex_calloc(N);
    gsl_vector *data = gsl_vector_calloc(N);
    gsl_fft_real_wavetable * wavetable = gsl_fft_real_wavetable_alloc (N);
    gsl_fft_real_workspace * workspace = gsl_fft_real_workspace_alloc (N);

    //--------------------------------------------------------------
    //Set data
    //--------------------------------------------------------------
    for(int i = 0; i< N; i++) gsl_vector_set(data, i, gsl_vector_get(dEv, i));

    //--------------------------------------------------------------
    //FFT transform
    //--------------------------------------------------------------
    gsl_fft_real_transform (data->data, 1, N, wavetable, workspace);
    gsl_fft_halfcomplex_unpack(data->data , data_complex->data ,  data->stride ,data->size);

    //--------------------------------------------------------------
    //Order 0
    //--------------------------------------------------------------
    switch(parity)
    {
    case 1: //even case (cosinus)
        xFFT.set_coef(+GSL_REAL(gsl_vector_complex_get(data_complex, 0))/(double)N,  0);
        break;
    case 0: //odd case (sinus)
        xFFT.set_coef(I*GSL_IMAG(gsl_vector_complex_get(data_complex, 0))/(double)N,  0);
        break;
    default: //unknown
        xFFT.set_coef(+GSL_REAL(gsl_vector_complex_get(data_complex, 0))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_complex, 0))/(double)N,  0);
        break;
    }

    //--------------------------------------------------------------
    //Order >< 0
    //--------------------------------------------------------------
    for(int i = 1; i<= fftN; i++)
    {
        switch(parity)
        {
        case 1: //even case (cosinus)
        {
            //Negative frequencies
            xFFT.set_coef(+GSL_REAL(gsl_vector_complex_get(data_complex, N-i))/(double)N, -i);
            //Positive frequencies
            xFFT.set_coef(+GSL_REAL(gsl_vector_complex_get(data_complex, i))/(double)N,  i);
            break;
        }
        case 0: //odd case (sinus)
        {
            //Negative frequencies
            xFFT.set_coef(I*GSL_IMAG(gsl_vector_complex_get(data_complex, N-i))/(double)N, -i);
            //Positive frequencies
            xFFT.set_coef(I*GSL_IMAG(gsl_vector_complex_get(data_complex, i))/(double)N,  i);
            break;
        }
        default: //unknown
        {
            //Negative frequecies
            xFFT.set_coef(+GSL_REAL(gsl_vector_complex_get(data_complex, N-i))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_complex, N-i))/(double)N, -i);
            //Positive frequencies
            xFFT.set_coef(+GSL_REAL(gsl_vector_complex_get(data_complex, i))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_complex, i))/(double)N,  i);
        }

        }
    }


    gsl_vector_complex_free(data_complex);
    gsl_vector_free(data);
    gsl_fft_real_wavetable_free(wavetable);
    gsl_fft_real_workspace_free(workspace);
}

/**
 *  \brief Fortran implementation of the FFT
 **/
void fortranFFT_real(Ofsc &xFFT, int fftN, int N, gsl_vector *dEv, int parity)
{

    //--------------------------------------------------------------
    //FFT tools
    //--------------------------------------------------------------
    int M = fftN;
    double CSF[M+1], SIF[M+1], F[N];
    for(int i = 0; i < N; i++) F[i]  = gsl_vector_get(dEv, i);

    //--------------------------------------------------------------
    //FFT (fortran, does not work on server, €€TODO)
    //--------------------------------------------------------------
    //foun_(F, &N, &M, CSF, SIF);
    cout << "fortranFFT_real. F[0] is displayed just to avoid compilation warning: " << F[0] << endl;

    //--------------------------------------------------------------
    //Order 0
    //--------------------------------------------------------------
    xFFT.set_coef(CSF[0],  0);


    //--------------------------------------------------------------
    //Order >< 0
    //--------------------------------------------------------------
    for(int i = 1; i<= fftN; i++)
    {
        switch(parity)
        {
        case 1: //even case (cosinus)
        {
            //Negative frequecies
            xFFT.set_coef(0.5*(CSF[i]), -i);
            //Positive frequencies
            xFFT.set_coef(0.5*(CSF[i]),  i);
            break;
        }
        case 0: //odd case (sinus)
        {
            //Negative frequecies
            xFFT.set_coef(0.5*(I*SIF[i]), -i);
            //Positive frequencies
            xFFT.set_coef(-0.5*(I*SIF[i]),  i);
            break;
        }
        default: //unknown
        {
            //Negative frequecies
            xFFT.set_coef(0.5*(CSF[i] + I*SIF[i]), -i);
            //Positive frequencies
            xFFT.set_coef(0.5*(CSF[i] - I*SIF[i]),  i);
        }

        }
    }
}

///**
// *  \brief GSL implementation of the FFT
// **/
//void gslFFT_complex(Ofsc &xFFT, int fftN, int N, gsl_vector *dEv, int parity)
//{
//    //TODO
//}

/**
 *  \brief FFT of the coefficients of a given matrix P obtained on a N points grid
 *
 *  Each FFt is tested on a fftPlot points grid with the use of the matrix Pt.
 *  Note on the flag (int): for the G1 matrix, symmetries allow us to clean the FFT by forcing some coefficients to zero
 */
void nfo2_fft(Ofsc& xFFT, gsl_matrix** P, gsl_matrix** Pt,  double n, int N, int fftN, int fftPlot, string filename, int flag, int is_stored)
{
    //-----------------------------------------
    // Init
    //-----------------------------------------
    //Storage vector
    gsl_vector *dEv = gsl_vector_calloc(N);

    //Storage tools (print in txt files)
    ofstream curentStream;
    string ss1, ss2;
    //Dimension of the matrix (either 6*6 for P, 3*3 for the Fij, 2*2 for G1 or 2*1 for the Xf)
    int size1 = P[0]->size1;
    int size2 = P[0]->size2;

    //-----------------------------------------
    //Loop on P
    //-----------------------------------------
    for(int i0 = 0; i0 < size1; i0++)
    {
        for(int j0 = 0; j0 < size2; j0++)
        {
            //-----------------------------------------
            // Store data in dEv
            //-----------------------------------------
            for(int i = 0; i< N; i++) gsl_vector_set(dEv, i, gsl_matrix_get(P[i], i0, j0));

            //-----------------------------------------
            // FFT with GSL tools
            //-----------------------------------------
            //For the G1, Xe, Xm and Xs matrices, symmetries allow us to clean the FFT by forcing some coefficients to zero:
            //G1[i][j] is even when (i+j)%2 == 0
            //            odd  when (i+j)%2 == 1
            if(flag)
            {
                if((i0+j0)%2 == 0) gslFFT_real(xFFT, fftN, N, dEv, 1);
                else gslFFT_real(xFFT, fftN, N, dEv, 0);
            }
            else gslFFT_real(xFFT, fftN, N, dEv, -1); //fortranFFT_real(xFFT, fftN, N, dEv, -1); /**/

            //-----------------------------------------
            //Storage in txt file, if necessary
            //-----------------------------------------
            if(is_stored)
            {
                cout << "nfo2_fft. The Fourier series is stored in txt file." << endl;
                ss1 = static_cast<ostringstream*>( &(ostringstream() << i0+1) )->str();
                ss2 = static_cast<ostringstream*>( &(ostringstream() << j0+1) )->str();
                curentStream.open((filename+ss1+ss2+".txt").c_str());
                curentStream << xFFT << endl;
                curentStream.close();
            }

            //-----------------------------------------
            //Each FFt is tested on a fftPlot points grid with the use of the matrix Pt
            //-----------------------------------------
            nfo2_test(xFFT, Pt, n, fftPlot, i0, j0);
        }
    }

    gsl_vector_free(dEv);
}

/**
 *  \brief Test of the validity of the FFT on a fftPlot points grid.
 */
void nfo2_test(Ofsc& xFFT, gsl_matrix** P, double n, int fftPlot, int i0, int j0)
{
    //Allocation
    gsl_vector * xxL2 = gsl_vector_calloc(fftPlot);
    double t1 = 2*M_PI/n;
    double ti = 0.0;

    //Test loop
    for(int i = 0; i < fftPlot; i++)
    {
        ti = i * t1 /fftPlot;
        //epsilon
        gsl_vector_set(xxL2, i, fabs(creal(xFFT.evaluate(n*ti))-gsl_matrix_get(P[i], i0, j0)));
    }
    cout << std::noshowpos << setiosflags(ios::scientific)  << setprecision(15);
    cout << "M" << i0 << j0 << ": ";
    cout << std::showpos << "||eps||_inf = " << gsl_vector_max(xxL2) << endl;
    gsl_vector_free(xxL2);
}


//----------------------------------------------------------------------------------------
// Tests
//----------------------------------------------------------------------------------------
/**
 *  \brief Symplectic test of P.
 */
void symplecticity_test_for_P(gsl_odeiv2_driver *d, const double y0[], double t1, gsl_matrix* R)
{
    //Allocation
    gsl_matrix* Pb = gsl_matrix_calloc (6, 6);
    gsl_matrix* Pi = gsl_matrix_calloc (6, 6);
    gsl_matrix* Pf = gsl_matrix_calloc (6, 6);

    //Reset the driver
    gsl_odeiv2_driver_reset(d);

    //Update the starting point (including P(0) = Id)
    double y[42];
    for(int i=0; i<42; i++) y[i] = y0[i];

    //Pb update
    gslc_vec_to_mat(Pb, y, 6, 6, 6);
    //Pi = Pb*R
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Pb, R, 0.0, Pi);

    //Time
    double t = 0.0;

    //Integration
    gsl_odeiv2_driver_apply (d, &t , t1 , y);

    //Pb update
    gslc_vec_to_mat(Pb, y, 6, 6, 6);
    //Pf = Pb*R
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Pb, R, 0.0, Pf);

    cout << "-----------------------------------------------" << endl;
    cout << "Symplectic test of the matrix Pb               " << endl;
    cout << "-----------------------------------------------" << endl;
    symplecticity_test_real(Pb, Csts::INVERSE_SYMP);
    cout << "-----------------------------------------------" << endl;
    cout << "Symplectic test of the matrix P                " << endl;
    cout << "-----------------------------------------------" << endl;
    symplecticity_test_real(Pf, Csts::INVERSE_SYMP);

    //Free memory
    gsl_matrix_free(Pb);
    gsl_matrix_free(Pi);
    gsl_matrix_free(Pf);
}

/**
 *  \brief Periodicity test of P.
 */
void nfo2_periodicity_test(gsl_odeiv2_driver *d, const double y0[], double n, gsl_matrix* R)
{
    //Allocation
    gsl_matrix* Pb = gsl_matrix_calloc (6, 6);
    gsl_matrix* Pi = gsl_matrix_calloc (6, 6);
    gsl_matrix* Pf = gsl_matrix_calloc (6, 6);

    //Reset the driver
    gsl_odeiv2_driver_reset(d);

    //Update the starting point (including P(0) = Id)
    double y[42];
    for(int i=0; i<42; i++) y[i] = y0[i];

    //Pb update
    gslc_vec_to_mat(Pb, y, 6, 6, 6);
    //Pi = Pb*R
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Pb, R, 0.0, Pi);

    //Time
    double t  = 0.0;
    double t1 = 2*M_PI/n;

    //Integration
    gsl_odeiv2_driver_apply (d, &t , t1 , y);

    //Pb update
    gslc_vec_to_mat(Pb, y, 6, 6, 6);
    //Pf = Pb*R
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Pb, R, 0.0, Pf);

    cout << "-----------------------------------------------" << endl;
    cout << "Periodicity test of the matrix P               " << endl;
    cout << "-----------------------------------------------" << endl;
    cout << "P[0]:" << endl;
    gslc_matrix_printf(Pi);
    cout << "P[end]:" << endl;
    gslc_matrix_printf(Pf);
    gsl_matrix_sub(Pi, Pf);
    cout << "P[0]-P[end]:" << endl;
    gslc_matrix_printf(Pi);

    //Free memory
    gsl_matrix_free(Pb);
    gsl_matrix_free(Pi);
    gsl_matrix_free(Pf);
}

/**
 *  \brief Symmetry test on P.
 */
void nfo2_symmetry_test(gsl_odeiv2_driver *d, const double y0[], double n, gsl_matrix* R, int p, int N)
{
    //Allocation
    gsl_matrix* Pb  = gsl_matrix_calloc (6, 6);
    gsl_matrix* Pi  = gsl_matrix_calloc (6, 6);
    gsl_matrix* Pf  = gsl_matrix_calloc (6, 6);
    gsl_matrix* Pfm = gsl_matrix_calloc (6, 6);

    //Reset the driver
    gsl_odeiv2_driver_reset(d);

    //Update the starting point (including P(0) = Id)
    double y[42];
    for(int i=0; i<42; i++) y[i] = y0[i];
    //Pb update
    gslc_vec_to_mat(Pb, y, 6, 6, 6);
    //Pi = Pb*R
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Pb, R, 0.0, Pi);
    //Time
    double t = 0.0;
    //Integration
    gsl_odeiv2_driver_apply (d, &t , (N-p)*2*M_PI/(n*N) , y);
    //Pb update
    gslc_vec_to_mat(Pb, y, 6, 6, 6);
    //Pf = Pb*R
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Pb, R, 0.0, Pf);

    //Re-init
    gsl_odeiv2_driver_reset(d);
    d->h = -d->h;
    for(int i=0; i<42; i++) y[i] = y0[i];
    t = 0.0;
    //Integration
    gsl_odeiv2_driver_apply (d, &t , (-p*2*M_PI/(n*N)) , y);
    //Pb update
    gslc_vec_to_mat(Pb, y, 6, 6, 6);
    //Pf = Pb*R
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Pb, R, 0.0, Pfm);
    d->h = -d->h;

    cout << "-----------------------------------------------" << endl;
    cout << "Symmetry test on P                             " << endl;
    cout << "-----------------------------------------------" << endl;
    cout << "Pf:" << endl;
    gslc_matrix_printf(Pf);
    cout << "Pf:" << endl;
    gslc_matrix_printf(Pfm);
    gsl_matrix_sub(Pf, Pfm);
    cout << "Pf-Pfm:" << endl;
    gslc_matrix_printf(Pf);

    //Free memory
    gsl_matrix_free(Pb);
    gsl_matrix_free(Pi);
    gsl_matrix_free(Pf);
    gsl_matrix_free(Pfm);
}

/**
 *  \brief Test of the symplectic character of a given complex matrix.
 */
void symplecticity_test_complex(const gsl_matrix_complex *M, int INVERSE_TYPE)
{
    cout << "---------------------------------------------" << endl;
    cout << "Test of the symplectic nature of a matrix.   " << endl;
    cout << "---------------------------------------------" << endl;
    int n = M->size1/2;
    gsl_complex one_c  = gslc_complex(1.0, 0.0);
    gsl_complex zero_c = gslc_complex(0.0, 0.0);

    gsl_matrix_complex *J     = gsl_matrix_complex_calloc(2*n, 2*n);
    gsl_matrix_complex *MH    = gsl_matrix_complex_calloc(2*n, 2*n);
    gsl_matrix_complex *AUX   = gsl_matrix_complex_calloc(2*n, 2*n);
    gsl_matrix_complex *BUX   = gsl_matrix_complex_calloc(2*n, 2*n);
    gsl_matrix_complex *Mcopy = gsl_matrix_complex_calloc (2*n, 2*n);


    //Transpose only
    gsl_matrix_complex_transpose_memcpy(MH, M);  //transpose only
    //gslc_matrix_complex_H_memcpy(MH, M);       //transpose + conjugate: probably not to be used!

    //Matrix J
    glsc_matrix_complex_set_J(J);

    cout << "---------------------------------------------" << endl;
    cout << "Initial matrix (real part) " << endl;
    cout << "---------------------------------------------" << endl;
    gslc_matrix_complex_printf_real(M);
    cout << "---------------------------------------------" << endl;
    cout << "Initial matrix (imag part) " << endl;
    cout << "---------------------------------------------" << endl;
    gslc_matrix_complex_printf_imag(M);


    //AUX = J*M
    gsl_blas_zgemm (CblasNoTrans , CblasNoTrans , one_c , J , M , zero_c , AUX );
    //BUX = MH*AUX = MH*J*M
    gsl_blas_zgemm (CblasNoTrans , CblasNoTrans , one_c , MH , AUX , zero_c , BUX );

//    cout << "---------------------------------------------" << endl;
//    cout << "BUX = MH*J*M (real.)" << endl;
//    cout << "---------------------------------------------" << endl;
//    gslc_matrix_complex_printf_real(BUX);
//    cout << "---------------------------------------------" << endl;
//    cout << "BUX = MH*J*M (imag.)" << endl;
//    cout << "---------------------------------------------" << endl;
//    gslc_matrix_complex_printf_imag(BUX);

    cout << "---------------------------------------------" << endl;
    cout << "BUX = MH*J*M (approx.)" << endl;
    cout << "---------------------------------------------" << endl;
    gslc_matrix_complex_printf_approx(BUX);

    //Inverse
    int s;
    gsl_permutation * p = gsl_permutation_alloc (2*n);
    if(INVERSE_TYPE == Csts::INVERSE_SYMP)
    {
        gslc_matrix_complex_symplectic_inverse(M, MH);

        //AUX = Minv*M
        gsl_blas_zgemm (CblasNoTrans , CblasNoTrans , one_c , MH , M , zero_c , AUX );

        //cout << "---------------------------------------------" << endl;
        //cout << "Minv*M (real part) using symplectic inverse." << endl;
        //cout << "---------------------------------------------" << endl;
        //gslc_matrix_complex_printf_real(AUX);
    }
    else
    {
        gsl_matrix_complex_memcpy(Mcopy, M);
        gsl_linalg_complex_LU_decomp (Mcopy, p , &s );
        gsl_linalg_complex_LU_invert (Mcopy, p , MH );

        //AUX = Minv*M
        gsl_blas_zgemm (CblasNoTrans , CblasNoTrans , one_c , MH , M , zero_c , AUX );

        //cout << "---------------------------------------------" << endl;
        //cout << "Minv*M (real part) using GSL inverse." << endl;
        //cout << "---------------------------------------------" << endl;
        //gslc_matrix_complex_printf_real(AUX);
    }


    cout << "---------------------------------------------" << endl;
    gsl_matrix_complex_free(J);
    gsl_matrix_complex_free(MH);
    gsl_matrix_complex_free(AUX);
    gsl_matrix_complex_free(BUX);
    gsl_matrix_complex_free(Mcopy);
    gsl_permutation_free(p);
}

/**
 *  \brief Test of the symplectic character of a given real matrix.
 */
void symplecticity_test_real(const gsl_matrix *M, int INVERSE_TYPE)
{
    cout << "---------------------------------------------" << endl;
    cout << "Test of the symplectic nature of a real matrix." << endl;
    cout << "---------------------------------------------" << endl;


    int n = M->size1/2;

    gsl_matrix *J     = gsl_matrix_calloc(2*n, 2*n);
    gsl_matrix *MH    = gsl_matrix_calloc(2*n, 2*n);
    gsl_matrix *AUX   = gsl_matrix_calloc(2*n, 2*n);
    gsl_matrix *BUX   = gsl_matrix_calloc(2*n, 2*n);
    gsl_matrix *Mcopy = gsl_matrix_calloc(2*n, 2*n);


    //Transpose only
    gsl_matrix_transpose_memcpy(MH, M);  //transpose only
    //gslc_matrix_complex_H_memcpy(MH, M);       //transpose + conjugate: probably not to be used!

    //Matrix J
    glsc_matrix_set_J(J);

    cout << "---------------------------------------------" << endl;
    cout << "Initial matrix" << endl;
    cout << "---------------------------------------------" << endl;
    gslc_matrix_printf(M);

    //AUX = J*M
    gsl_blas_dgemm (CblasNoTrans , CblasNoTrans , 1.0, J , M , 0.0 , AUX );
    //BUX = MH*AUX = MH*J*M
    gsl_blas_dgemm (CblasNoTrans , CblasNoTrans , 1.0 , MH , AUX , 0.0 , BUX );

//    cout << "---------------------------------------------" << endl;
//    cout << "M^T*J*M" << endl;
//    cout << "---------------------------------------------" << endl;
//    gslc_matrix_printf(BUX);

    cout << "---------------------------------------------" << endl;
    cout << "BUX = MH*J*M (approx.)" << endl;
    cout << "---------------------------------------------" << endl;
    gslc_matrix_printf_approx(BUX);


    //Inverse
    int s;
    gsl_permutation * p = gsl_permutation_alloc (2*n);
    if(INVERSE_TYPE == Csts::INVERSE_SYMP)
    {
        gslc_matrix_symplectic_inverse(M, MH);

        //AUX = Minv*M
        gsl_blas_dgemm(CblasNoTrans , CblasNoTrans , 1.0 , MH , M , 0.0 , AUX);

        //cout << "---------------------------------------------" << endl;
        //cout << "Minv*M  using symplectic inverse." << endl;
        //cout << "---------------------------------------------" << endl;
        //gslc_matrix_printf(AUX);
    }
    else
    {
        gsl_matrix_memcpy(Mcopy, M);
        gsl_linalg_LU_decomp (Mcopy, p , &s );
        gsl_linalg_LU_invert (Mcopy, p , MH );

        //AUX = Minv*M
        gsl_blas_dgemm (CblasNoTrans , CblasNoTrans , 1.0 , MH , M , 0.0 , AUX );

        //cout << "---------------------------------------------" << endl;
        //cout << "Minv*M using GSL inverse." << endl;
        //cout << "---------------------------------------------" << endl;
        //gslc_matrix_printf(AUX);
    }

    gsl_matrix_free(J);
    gsl_matrix_free(MH);
    gsl_matrix_free(AUX);
    gsl_matrix_free(BUX);
    gsl_matrix_free(Mcopy);
    gsl_permutation_free(p);
}

/**
 *  \brief Test of the monodromy eigensystem:
 *         For each eigencouple (y, ly) of the matrix Mc = DAT[1]*...*DAT[M],
 *         computation of the scalars |ly-Mc*y/ly| and  |y-Mc^-1*y*ly|.
 */
void eigensystem_test(gsl_matrix_complex *Dm, gsl_matrix_complex *S, gsl_matrix_complex **DAT, int M)
{
    //------------------------------------------------------------------------------------
    // Initialize the GSL objects
    //------------------------------------------------------------------------------------
    gsl_vector_complex_view evec_i;
    gsl_complex eval_i;
    gsl_vector_complex *x = gsl_vector_complex_calloc(6);

    //------------------------------------------------------------------------------------
    // Test of the eigenvectors/eigenvalues
    //------------------------------------------------------------------------------------
    cout << Csts::SSEPR << endl;
    cout << "Test of the eigensystem decomposition.         " << endl;
    cout << "The eigenvectors/eigenvalues are denoted y/ly.       " << endl;
    cout << Csts::SSEPR << endl;
    cout << "M*y is computed with the matrix M given as a product. "<< endl;
    cout << "    real(Eigenvalue) (ly)           |ly-M*y/ly| "   << endl;
    for(int eigN = 0; eigN < 6; eigN ++)
    {
        evec_i = gsl_matrix_complex_column (S, eigN);
        eval_i = gsl_matrix_complex_get (Dm, eigN, eigN);
        gslc_matrix_vector_product(DAT, &evec_i.vector , x, M);
        gsl_vector_complex_scale(x, gsl_complex_inverse(eval_i));
        //x = x-&evec_i.vector
        gsl_vector_complex_sub(x, &evec_i.vector);
        cout << eigN << " " << GSL_REAL(eval_i) << "   " << gsl_blas_dznrm2(x) << endl;
    }

    //------------------------------------------------------------------------------------
    // Test of the eigenvectors/eigenvalues with the inverse matrix
    //------------------------------------------------------------------------------------
    cout << Csts::SSEPR << endl;
    cout << "Test with matrix inverse, still with the matrix M given as a product. " << endl;
    cout << Csts::SSEPR << endl;
    cout << "    real(Eigenvalue) (ly)          |y-M^-1*y*ly| "   << endl;
    for(int eigN = 0; eigN < 6; eigN ++)
    {
        evec_i = gsl_matrix_complex_column (S, eigN);
        eval_i = gsl_matrix_complex_get (Dm, eigN, eigN);

        gslc_matrix_vector_invproduct(DAT, &evec_i.vector, x, M);
        gsl_vector_complex_scale(x, eval_i);
        gsl_vector_complex_sub(x, &evec_i.vector);
        cout << eigN << " " << GSL_REAL(eval_i) << "   " << gsl_blas_dznrm2(x) << endl;

    }

    cout << Csts::SSEPR << endl;
    cout << "Note: the eigenvalue leading to bad precision  " << endl;
    cout << "with the direct matrix corresponds to the stable direction " << endl;
    cout << "(Obtained by inverse power method).            " << endl;

    //------------------------------------------------------------------------------------
    // Comparison of the modulus of the eigenvectors with one
    //------------------------------------------------------------------------------------
    cout << Csts::SSEPR << endl;
    cout << "Modulus of the eigenvalues and comparison to 1" << endl;
    cout << Csts::SSEPR << endl;
    for(int eigN = 0; eigN < 6; eigN ++)
    {
        eval_i = gsl_matrix_complex_get (Dm, eigN, eigN);
        cout << eigN << " " << gsl_complex_abs(eval_i) << " comparison to 1: " << 1-gsl_complex_abs(eval_i) << endl;

    }

    //------------------------------------------------------------------------------------
    //Free
    //------------------------------------------------------------------------------------
    gsl_vector_complex_free(x);
}

//----------------------------------------------------------------------------------------
// Change of base (currently not used)
//----------------------------------------------------------------------------------------
/**
 *  \brief Change of base: Dm = Sinv*M*S with M given as a product of matrices
 */
void change_of_base(gsl_matrix_complex *MMc_estimate, gsl_matrix_complex *Dm, gsl_matrix_complex **DAT, const gsl_matrix_complex *S, int INVERSE_TYPE, int M)
{

    cout << "-------------------------------------------------------------------" << endl;
    cout << "Test of the change of base of the Monodromy matrix (Dm = Sinv*M*S) " << endl;
    cout << "The monodromy matrix is used as a product of matrices              " << endl;
    cout << "-------------------------------------------------------------------" << endl;
    int n = Dm->size1/2;
    gsl_complex one_c  = gslc_complex(1.0, 0.0);
    gsl_complex zero_c = gslc_complex(0.0, 0.0);

    gsl_matrix_complex *AUX = gsl_matrix_complex_calloc(2*n, 2*n);
    gsl_matrix_complex *BUX = gsl_matrix_complex_calloc(2*n, 2*n);
    gsl_matrix_complex *Sinv = gsl_matrix_complex_calloc(2*n, 2*n);
    gsl_matrix_complex *SLU = gsl_matrix_complex_calloc (2*n, 2*n);
    gsl_vector_complex_view evec_i;
    gsl_vector_complex *x = gsl_vector_complex_calloc(6);

    //Inverse
    int s;
    gsl_permutation * p = gsl_permutation_alloc (2*n);
    if(INVERSE_TYPE == Csts::INVERSE_SYMP)
    {
        gslc_matrix_complex_symplectic_inverse(S, Sinv);
    }
    else
    {
        gsl_matrix_complex_memcpy(SLU, S);
        gsl_linalg_complex_LU_decomp (SLU, p , &s );
        gsl_linalg_complex_LU_invert (SLU, p , Sinv );
    }

    //AUX = DAT[n]*...*DAT[1]*S
    gslc_matrix_matrix_product(DAT, S, AUX, M);

    //BUX = Sinv*AUX = Sinv*M*S
    if(INVERSE_TYPE == Csts::INVERSE_SYMP)
    {
        //Using  Sinv
        gsl_blas_zgemm (CblasNoTrans , CblasNoTrans , one_c , Sinv , AUX , zero_c , BUX );
    }
    else
    {
        //Through inversion of eigensystems on the column of AUX
        for(int i =0; i<6; i++)
        {
            evec_i = gsl_matrix_complex_column (AUX, i);
            gsl_linalg_complex_LU_solve (SLU , p , &evec_i.vector , x );
            evec_i = gsl_matrix_complex_column (BUX, i);
            gsl_vector_complex_memcpy(&evec_i.vector , x);
        }
    }

    cout << "---------------------------------------------" << endl;
    cout << "Sinv*M*S (real part)" << endl;
    cout << "---------------------------------------------" << endl;
    gslc_matrix_complex_printf_real(BUX);


    gsl_matrix_complex *SDm  = gsl_matrix_complex_calloc(6, 6);
    gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, one_c , S, Dm , zero_c , SDm);
    gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, one_c , SDm , Sinv, zero_c , MMc_estimate);
    cout << "---------------------------------------------" << endl;
    cout << "S*Dm*Sinv (real part)" << endl;
    cout << "---------------------------------------------" << endl;
    gslc_matrix_complex_printf_real(MMc_estimate);



    gsl_matrix_complex_free(AUX);
    gsl_matrix_complex_free(BUX);
    gsl_matrix_complex_free(Sinv);
    gsl_matrix_complex_free(SLU);
    gsl_permutation_free(p);
}

/**
 *  \brief Change of base: Dm = Sinv*M*S with M given as a single matrix
 */
void change_of_base_single_matrix(gsl_matrix_complex *Dm, gsl_matrix_complex *M, gsl_matrix_complex *S, int INVERSE_TYPE)
{
    int n = M->size1;
    gsl_complex one_c  = gslc_complex(1.0, 0.0);
    gsl_complex zero_c = gslc_complex(0.0, 0.0);


    gsl_matrix_complex *AUX = gsl_matrix_complex_calloc(n, n);
    gsl_matrix_complex *Sinv = gsl_matrix_complex_calloc(n, n);
    gsl_matrix_complex *SLU = gsl_matrix_complex_calloc (n, n);

    //Inverse
    int s;
    gsl_permutation * p = gsl_permutation_alloc (n);
    if(INVERSE_TYPE == Csts::INVERSE_SYMP)
    {
        gslc_matrix_complex_symplectic_inverse(S, Sinv);
    }
    else
    {
        gsl_matrix_complex_memcpy(SLU, S);
        gsl_linalg_complex_LU_decomp (SLU, p , &s );
        gsl_linalg_complex_LU_invert (SLU, p , Sinv );
    }

    //AUX = M*S
    gsl_blas_zgemm (CblasNoTrans , CblasNoTrans , one_c , M , S , zero_c , AUX );
    //Dm = Sinv*AUX = Sinv*M*S
    gsl_blas_zgemm (CblasNoTrans , CblasNoTrans , one_c , Sinv , AUX , zero_c , Dm );


    cout << "---------------------------------------------" << endl;
    cout << "Dm = Sinv*M*S (real)" << endl;
    cout << "---------------------------------------------" << endl;
    gslc_matrix_complex_printf_real(Dm);

    cout << "---------------------------------------------" << endl;
    cout << "Dm = Sinv*M*S (approx.)" << endl;
    cout << "---------------------------------------------" << endl;
    gslc_matrix_complex_printf_approx(Dm);


    gsl_matrix_complex_free(AUX);
    gsl_matrix_complex_free(Sinv);
    gsl_matrix_complex_free(SLU);
    gsl_permutation_free(p);
}

//----------------------------------------------------------------------------------------
// Wielandt's deflation
//----------------------------------------------------------------------------------------
/**
 *  \brief Application of the Wielandt deflation to the monodromy matrix. Algorithm taken from "Introduction to Numerical Analysis with C programs", A. Mate, 2004.
 *  \param MMc      the monodromy matrix to diagonalize.
 *  \param eigenVu  the unstable eigenvector, obtain via power method applied on MMc = DAT[M]*DAT[M-1]*...*DAT[1].
 *  \param eigenLu  the unstable eigenvalue,  obtain via power method.
 *  \param evecr    matrix of eigenvectors (in columns) of MMc (output).
 *  \param evalr    vectors of eigenvalues of MMc (output).
 *
 *  The algorithm is based on the following results:
 *  Given a eigenvector/value couple \f$ (\lambda,\mathbf{x}) \f$ of \f$ \mathbf{M} \f$, one can build the matrix
 *  \f$ \mathbf{B} = \mathbf{M} - \lambda \mathbf{x} \mathbf{z}^T \f$ where
 *  \f$ \mathbf{z} = \frac{1}{\lambda x_r} \mathbf{M}(r,*)^T \f$ where
 *  \f$ x_r \f$ is the component of \f$ \mathbf{x} \f$ of maximum absolute value, which guarantees \f$ x_r \neq 0 \f$.
 *
 *  Then, the if \f$ \rho \neq 0 \f$ is an eigenvalue of \f$ \mathbf{M} \f$ different from \f$ \lambda \f$, then \f$ \rho \f$ is also an eigenvalue of \f$ \mathbf{B} \f$.
 *  Moreover, if \f$ \mathbf{w} \f$ is the eigenvector of \f$ \mathbf{B} \f$ associated to \f$ \rho \f$ , then the eigenvector \f$ \mathbf{y} \f$ of \f$ \mathbf{M} \f$
 *  associated to \f$ \rho \f$ is given by:
 *
 *  \f$
 *  $$
 *  \mathbf{y} = \frac{1}{\rho} \left(\mathbf{w} + \frac{\lambda}{\rho - \lambda} (\mathbf{z}^T \mathbf{w}) \mathbf{x} \right).
 *  $$
 *  \f$
 **/
void wielandt_deflation(gsl_matrix_complex *MMc,
                        gsl_vector_complex *eigenVu,
                        gsl_complex eigenLu,
                        gsl_matrix_complex *evecr,
                        gsl_vector_complex *evalr)
{
    gsl_eigen_nonsymmv_workspace * wr = gsl_eigen_nonsymmv_alloc(6);     //workspace for eigenspaces determination
    gsl_matrix_complex *AUX   = gsl_matrix_complex_calloc(6, 6);          //auxiliary matrix
    gsl_matrix_complex *BUX   = gsl_matrix_complex_calloc(6, 6);          //auxiliary matrix
    gsl_matrix *Breal         = gsl_matrix_calloc(6, 6);                  //auxiliary matrix
    gsl_vector_complex *vux   = gsl_vector_complex_calloc(6);             //auxiliary vector
    gsl_vector* eigenVabs     = gsl_vector_calloc(6);                     //abs(eigenVu)
    int maxIndex;
    gsl_complex lambda;

    //----------------------------------------------------------------------------------
    // 1. AUX = MMc - eigenLu x*z^T  (see comments)
    //----------------------------------------------------------------------------------
    //AUX = MMc
    gsl_matrix_complex_memcpy(AUX, MMc);
    //eigenVabs = |eigenVu|
    for(int i = 0; i<6; i++) gsl_vector_set(eigenVabs, i, gsl_complex_abs(gsl_vector_complex_get(eigenVu, i)));
    //maxIndex = arsi_max_i |eigenVu(i)|
    maxIndex = gsl_vector_max_index (eigenVabs);
    //vux = maxIndex row of MMc
    gslc_matrix_complex_row(vux, AUX, maxIndex);
    //lambda = 1.0/(eigenLu*eigenVu[maxIndex])
    lambda = gsl_complex_inverse(gsl_complex_mul(eigenLu, gsl_vector_complex_get(eigenVu, maxIndex)));
    //vux *= lambda
    gsl_vector_complex_scale(vux, lambda);
    //BUX = eigenLu*eigenVu*vux^T
    gsl_blas_zgeru(eigenLu , eigenVu , vux , BUX );
    //AUX -= BUX
    gsl_matrix_complex_sub(AUX, BUX);

    //----------------------------------------------------------------------------------
    // 2. GSL routine applied again to get central directions
    //----------------------------------------------------------------------------------
    for(int i=0; i<6; i++) for(int j=0; j<6; j++) gsl_matrix_set(Breal, i, j,  GSL_REAL(gsl_matrix_complex_get(AUX, i, j)));
    gsl_eigen_nonsymmv (Breal, evalr, evecr, wr);
    gsl_eigen_nonsymmv_sort (evalr, evecr, GSL_EIGEN_SORT_ABS_ASC);

    //----------------------------------------------------------------------------------
    // 3. Wielandt inverse transformation to get the central directions in the good frame
    //----------------------------------------------------------------------------------
    gsl_vector_complex_view evec_i;
    gsl_complex eval_i;
    for(int i =0; i<6; i++)
    {
        eval_i = gsl_vector_complex_get (evalr, i);
        evec_i = gsl_matrix_complex_column (evecr, i);
        gslc_wielandt_inv_trans(&evec_i.vector, eval_i, eigenVu, eigenLu, vux);
    }
}
