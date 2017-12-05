/**
 * \file  diffcorr.cpp
 * \brief Contains all the routines that perform differential correction procedures, as
 *        well as some plotting routines (that should probably be set in a separate src
 *        file in the long run).
 * \author BLB.
 */

#include "diffcorr.h"


//----------------------------------------------------------------------------------------
//  Inner routines
//----------------------------------------------------------------------------------------
/**
 *  \brief Solve a definite-positive linear system:
 *         - Decompose a ncs x ncs matrix M that is definite-positive, using GSL routines.
 *         - Then inverse the system err_vec = M*corr_vec.
 **/
int ftc_inv_dfls(gsl_matrix* M, gsl_vector* err_vec, gsl_vector* corr_vec, int ncs)
{
    //Name of the routine
    string fname = "ftc_inv_dfls";

    //====================================================================================
    //We start by turning off the handler in the next chunk of code, so that we can test
    //if M is truly definite-positive.
    //We will restore it at the end of the decomposition + inversion process.
    //====================================================================================
    gsl_error_handler_t* error_handler = gsl_set_error_handler_off();

    //====================================================================================
    //Since M is in theory definite-positive,
    //a Cholesky decomposition can be used, instead of a LU decomposition.
    //====================================================================================
    int status = gsl_linalg_cholesky_decomp (M);

    //====================================================================================
    //If the decomposition went bad, status is different from GSL_SUCCESS. For example,
    //if M is not strictly definite-positive, the constant GSL_EDOM is returned by
    //gsl_linalg_cholesky_decomp. In any case, we try a LU decomposition.
    //====================================================================================
    if(status)
    {
        cerr << fname << ". Cholesky decomposition failed. A LU decomposition is tested. ref_errno = " << gsl_strerror(status) << endl;
        int s;
        gsl_permutation* p = gsl_permutation_alloc (ncs);
        status = gsl_linalg_LU_decomp (M, p, &s);
        if(!status) status = gsl_linalg_LU_solve(M, p, err_vec, corr_vec);
        gsl_permutation_free (p);

        //--------------------------------------------------------------------------------
        //If the decomposition or solving went bad, status is different from GSL_SUCCESS.
        //--------------------------------------------------------------------------------
        if(status)
        {
            cerr << fname << ". The LU decomposition failed. ref_errno = " << gsl_strerror(status) << endl;
            //Finally, we return an error to the user
            gsl_set_error_handler (error_handler);
            return GSL_FAILURE;
        }

    }
    else  //if not, the decomposition went fine, and we can go on with the inversion.
    {
        status = gsl_linalg_cholesky_solve(M, err_vec, corr_vec);
        //--------------------------------------------------------------------------------
        //If the solving went bad, status is different from GSL_SUCCESS.
        //--------------------------------------------------------------------------------
        if(status)
        {
            cerr << fname << ". Cholesky solving failed. ref_errno = " << gsl_strerror(status) << endl;
            //Finally, we return an error to the user
            gsl_set_error_handler (error_handler);
            return GSL_FAILURE;
        }
    }

    //We check that the result is not undefined
    if(gsl_isnan(gsl_vector_max(corr_vec)))
    {
        cerr << fname << ". The decomposition + solving failed, result contains NaN" << endl;
        //Finally, we return an error to the user
        gsl_set_error_handler (error_handler);
        return GSL_FAILURE;
    }


    //====================================================================================
    //We restore the error handler at the end
    //====================================================================================
    gsl_set_error_handler (error_handler);
    return GSL_SUCCESS;
}

/**
 *  \brief Computes the correction vector associated to the minimum norm solution.
 *         Given:
 *              - an ncs x 1   error vector err_vec
 *              - an nfv x ncs Jacobian err_jac,
 *         This routine computes the correction vector associated to
 *         the minimum norm solution:
 *
 *              corr_vec = err_jac^T x (err_jac x err_jac^T)^{-1} err_vec.
 **/
int ftc_corrvec_mn(gsl_vector* mn_vec, gsl_vector *err_vec, gsl_matrix* err_jac, int nfv, int ncs)
{
    //Name of the routine
    string fname = "ftc_corrvec_mn";

    //====================================================================================
    // Initialize the local GSL objects
    //====================================================================================
    gsl_matrix *temp_mat  = gsl_matrix_calloc(ncs, ncs);
    gsl_vector *temp_vec  = gsl_vector_calloc(ncs);

    //====================================================================================
    // Update temp_mat = err_jac x err_jac^T
    //====================================================================================
    int status = gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, err_jac , err_jac, 0.0, temp_mat);
    if(status)
    {
        cerr << fname << ". The computation of temp_mat failed. ref_errno = " << gsl_strerror(status) << endl;
        gsl_matrix_free(temp_mat);
        gsl_vector_free(temp_vec);
        return GSL_FAILURE;
    }

    //====================================================================================
    // Update correction vector corr_vec = err_jac^T x (err_jac x err_jac^T)^{-1} err_vec
    //====================================================================================
    //Inverse err_vec = temp_mat*temp_vec via Cholesky decomposition of temp_mat
    status = ftc_inv_dfls(temp_mat, err_vec, temp_vec, ncs);
    //We check that the inversing went well
    if(status)
    {
        cerr << fname << ". The decomposition + solving failed. ref_errno = " << gsl_strerror(status) << endl;
        gsl_matrix_free(temp_mat);
        gsl_vector_free(temp_vec);
        return GSL_FAILURE;
    }
    //Then, corr_vec = err_jac^T*temp_vec
    status = gsl_blas_dgemv(CblasTrans, -1.0, err_jac, temp_vec, 0.0, mn_vec);
    if(status)
    {
        cerr << fname << ". The computation of corr_vec failed. ref_errno = " << gsl_strerror(status) << endl;
        gsl_matrix_free(temp_mat);
        gsl_vector_free(temp_vec);
        return GSL_FAILURE;
    }

    //====================================================================================
    // Kill the local GSL objects and return GSL_SUCCESS
    //====================================================================================
    gsl_matrix_free(temp_mat);
    gsl_vector_free(temp_vec);
    return GSL_SUCCESS;
}

/**
 *  \brief Computes the correction vector for a square system.
 *         Given:
 *              - an ncs x 1   error vector err_vec
 *              - an ncs x ncs Jacobian err_jac,
 *         This routine computes the correction vector given by:
 *
 *              corr_vec = err_jac^{-1} x err_vec.
 **/
int ftc_corrvec_square(gsl_vector* corr_vec, gsl_vector *err_vec, gsl_matrix* err_jac, int ncs)
{
    //Name of the routine
    string fname = "ftc_corrvec_square";

    //====================================================================================
    // Initialize the local GSL objects
    //====================================================================================
    gsl_matrix *M   = gsl_matrix_calloc(ncs, ncs);
    gsl_permutation* p = gsl_permutation_alloc (ncs);
    int s;


    //M = err_jac
    gsl_matrix_memcpy(M, err_jac);

    //====================================================================================
    // Update correction vector corr_vec = err_jac^{-1} x err_vec
    //====================================================================================
    //Inverse err_vec = M*corr_vec via LU decomposition of M
    int status = gsl_linalg_LU_decomp (M, p, &s);
    if(!status) status = gsl_linalg_LU_solve(M, p, err_vec, corr_vec);

    //--------------------------------------------------------------------------------
    //If the decomposition or solving went bad, status is different from GSL_SUCCESS.
    //--------------------------------------------------------------------------------
    if(status)
    {
        cerr << fname << ". The LU decomposition failed. ref_errno = " << gsl_strerror(status) << endl;
        gsl_permutation_free (p);
        gsl_matrix_free(M);
        return GSL_FAILURE;
    }

    //====================================================================================
    // Kill the local GSL objects and return FTC_SUCCESS
    //====================================================================================
    gsl_permutation_free (p);
    gsl_matrix_free(M);
    return GSL_SUCCESS;
}



//----------------------------------------------------------------------------------------
//  Single shooting
//----------------------------------------------------------------------------------------
/**
 *  \brief Performs a differential correction procedure on zv0 in order to get a
 *         periodic orbit of period t_period.
 *         The algorithm assumes that the orbit is in the xy plane, is symmetric wrt to
 *         the x-axis, and has a period T = t_period.
 *
 *         The number of variable is 42: 6 for the state, 36 for the
 *         State Transition Matrix (STM).
 *
 *  \param zv0 the initial conditions that needs to be corrected.
 *  \param t_period the final integration time (period of the orbit).
 *  \param eps_diff: the desired precision (magnitude of the final error).
 *  \param d:  a gsl_odeiv2_driver object used to integrate the equations of motion.
 *  \param n_var: the dimension of zv0.
 *  \param is_plot: if true, the steps of the diffcorr procedure are plotted in a
 *         temporary gnuplot window.
 *
 *   In brief: the algorithm corrects  [ zv0[0] zv0[4] ] in order to get
 *   [y[1] y[3]](t1) = 0, i.e. when the trajectory crosses the line t = t1.
 **/
int single_shoot_sym(double zv0[], double t_period, double eps_diff,
                     gsl_odeiv2_driver *d, int is_plot)
{
    //====================================================================================
    // Init
    //====================================================================================
    int n_var = 42;
    int iter = 0;
    int status = 0;
    int itermax = 1000;

    double zv_int[n_var], zv_corr[n_var];
    double dzvf[2], dzv0[2];
    dzv0[0] = 0.0;
    dzv0[1] = 0.0;
    double t;

    double A[2][2];
    double invA[2][2];
    double deter;
    double eps_c = 10;
    double eps_m = 10;

    //====================================================================================
    // Optionnal plotting
    //====================================================================================
    gnuplot_ctrl  *hc = gnuplot_init();
    if(is_plot)
    {
        gnuplot_setstyle(hc, (char*)"lines");
        gnuplot_set_xlabel(hc, (char*)"X [-]");
        gnuplot_set_ylabel(hc, (char*)"Y [-]");
    }

    //====================================================================================
    //Update zv_corr which will be corrected in the loop (STM (0) = id_mat included)
    //====================================================================================
    for(int i=0; i< n_var; i++) zv_corr[i] = zv0[i];

    //====================================================================================
    //Correction loop
    //====================================================================================
    do
    {
        //================================================================================
        //Update zv_int (STM (0) = id_mat included)
        //================================================================================
        for(int i=0; i< n_var; i++) zv_int[i] = zv_corr[i];
        //Integration until t=t_period is reached
        t = 0.0;
        //Optionnal plotting
        if(is_plot)
        {
            ode_plot_xy(zv_int, n_var, t_period, d, hc, 500, 8, true, false, "Diff Corr n°"+iter, "");
        }

        //================================================================================
        //Integration
        //================================================================================
        gsl_odeiv2_driver_reset(d);
        set_dir_driver(d, t, t_period);
        status = gsl_odeiv2_driver_apply (d, &t , t_period , zv_int);

        if(status != GSL_SUCCESS)
        {
            printf("single_shoot_sym. WARNING: GSL driver failed to converge. Premature ending.\n");
            return GSL_FAILURE;
        }

        //================================================================================
        //Update the matrix A
        //================================================================================
        A[0][0] = zv_int[5+7];  //Phi21
        A[0][1] = zv_int[5+11]; //Phi25

        A[1][0] = zv_int[5+19]; //Phi41
        A[1][1] = zv_int[5+23]; //Phi45

        //Inverse A
        deter      = +A[0][0]*A[1][1] - A[0][1]*A[1][0];
        invA[0][0] = +A[1][1]/deter;
        invA[0][1] = -A[0][1]/deter;
        invA[1][0] = -A[1][0]/deter;
        invA[1][1] = +A[0][0]/deter;

        //================================================================================
        //Update the final error
        //================================================================================
        dzvf[0] = -zv_int[1];
        dzvf[1] = -zv_int[3];

        //================================================================================
        //Check that yc leads to a better precision, and update the output if so.
        //================================================================================
        //Current precision
        eps_c = sqrt(dzvf[0]*dzvf[0] + dzvf[1]*dzvf[1]);

        //Compare to the best precision yet
        if(eps_c < eps_m)
        {
            eps_m = eps_c;
            for(int i=0; i< n_var; i++) zv0[i] = zv_corr[i];
        }

        //================================================================================
        //Break the loop if precision is good
        //================================================================================
        if(eps_m < eps_diff) break;

        //================================================================================
        //Inversion of the error
        //================================================================================
        dzv0[0] = invA[0][0]*dzvf[0] + invA[0][1]*dzvf[1];
        dzv0[1] = invA[1][0]*dzvf[0] + invA[1][1]*dzvf[1];

        //================================================================================
        //Update the state
        //================================================================================
        zv_corr[0] += 0.1*dzv0[0];
        zv_corr[4] += 0.1*dzv0[1];
    }
    while((++iter) < itermax);

    //====================================================================================
    // Warning message if no convergence
    //====================================================================================
    if(iter>=itermax)
    {
        printf("WARNING: number of iter max exceeded during differential correction. Final precision is around %+5.2e. Premature ending\n", eps_m);
        return GSL_FAILURE;
    }

    //====================================================================================
    // Display info
    //====================================================================================
    //cout << "single_shoot_sym. Required eps = " << eps_diff << endl;
    //cout << "single_shoot_sym. Obtained eps = " << eps_m    << endl;

    gnuplot_close(hc);
    return GSL_SUCCESS;
}


//----------------------------------------------------------------------------------------
//  Multiple shooting
//----------------------------------------------------------------------------------------
/**
 * \brief Multiple shooting scheme with periodicity conditions, on n_patch+1 patch points.
 *
 *   \param (zmd[0:41][0:n_patch], tmd[0:n_patch]) is the first guess (input).
 *   \param (zmdn[0:41][0:n_patch], tmdn[0:n_patch]) is the corrected output.
 *   \param d:  a gsl_odeiv2_driver object used to integrate the equations of motion.
 *   \param n_patch: the number of patch points
 **/
int mult_shoot_period(double** zmd, double* tmd, double** zmdn, double* tmdn,
                      gsl_odeiv2_driver* d, int n_patch, double eps_diff,
                      int is_plot, gnuplot_ctrl* h1, int strong_conv)
{
    //====================================================================================
    // 1. Initialization
    //====================================================================================
    //Name of the routine
    string fname = "mult_shoot_period";

    //Cumulated norm of the error
    double eps_c = 0.0, eps_c_old = 1e5;

    //Status
    int status = GSL_SUCCESS;

    //Number of variables: 6 (state) + 36 (STM)
    int n_var = 42;

    //------------------------------------------------------------------------------------
    // GSL matrices and vectors
    //------------------------------------------------------------------------------------
    int nfv = 6*(n_patch+1); //free variables
    int ncs = 6*(n_patch);   //constraints

    // Correction vector at patch points
    gsl_vector *corr_vec = gsl_vector_calloc(nfv);
    // Error vector at patch points
    gsl_vector *err_vec  = gsl_vector_calloc(ncs);
    //Jacobian at patch points
    gsl_matrix **err_jac_patch = gslc_matrix_array_calloc(6, 6, n_patch);
    gsl_matrix *err_jac  = gsl_matrix_calloc(ncs, nfv);

    //Identity matrix eye(6)
    gsl_matrix *id_mat = gsl_matrix_calloc(6,6);
    gsl_matrix_set_identity (id_mat);

    //------------------------------------------------------------------------------------
    // Copy the departure state in zmdn
    //------------------------------------------------------------------------------------
    for(int k = 0; k <= n_patch; k++)
    {
        for(int i = 0; i < n_var; i++) zmdn[i][k] = zmd[i][k];
        tmdn[k] = tmd[k];
    }

    //====================================================================================
    // 2. Loop correction
    //====================================================================================
    //Maximum number of iterations is retrieved from config manager
    int itermax = Config::configManager().G_DC_ITERMAX();
    int iter = 0;

    //Reset driver
    gsl_odeiv2_driver_reset(d);
    set_dir_driver(d, tmd[0], tmd[1]);

    //Loop
    while(iter < itermax)
    {
        //--------------------------------------------------------------------------------
        // Build the Jacobian and other useful matrices
        //--------------------------------------------------------------------------------
        for(int k = 0; k <= n_patch-1; k++)
        {
            //----------------------------------------------------------------------------
            // Initialization for // computing
            //----------------------------------------------------------------------------
            double zve[n_var], tv = 0;

            //----------------------------------------------------------------------------
            // Integration
            //----------------------------------------------------------------------------
            for(int i = 0; i < n_var; i++) zve[i] = zmdn[i][k];
            tv = tmdn[k];

            //Storing eye(6) into the initial vector
            gslc_mat_to_vec(zve, id_mat, 6, 6, 6);

            //Integration in [tmdn[k], tmdn[k+1]]
            gsl_odeiv2_driver_apply (d, &tv, tmdn[k+1], zve);

            //----------------------------------------------------------------------------
            // Update the Jacobian
            //----------------------------------------------------------------------------
            gslc_vec_to_mat(err_jac_patch[k], zve, 6, 6, 6);

            //----------------------------------------------------------------------------
            // Update the error vector: F[k] = [zve[k] - zmdn[k+1]]
            //----------------------------------------------------------------------------
            for(int i = 0; i < 6; i++)
            {
                if(k != n_patch-1) gsl_vector_set(err_vec, i + 6*k, zve[i] - zmdn[i][k+1]);
                else gsl_vector_set(err_vec, i + 6*k, zve[i] - zmdn[i][0]);
            }

            //----------------------------------------------------------------------------
            // Update err_jac
            //----------------------------------------------------------------------------
            for(int i = 0; i < 6; i++)
            {
                if(k != n_patch-1)
                {
                    for(int j = 0; j < 6; j++)
                    {
                        gsl_matrix_set(err_jac, i + 6*k, j + 6*k,      gsl_matrix_get(err_jac_patch[k], i, j));
                        gsl_matrix_set(err_jac, i + 6*k, j + 6*(k+1), -gsl_matrix_get(id_mat, i, j));
                    }
                }
                else
                {
                    for(int j = 0; j < 6; j++)
                    {
                        gsl_matrix_set(err_jac, i + 6*k, j + 6*k,      gsl_matrix_get(err_jac_patch[k], i, j));
                        gsl_matrix_set(err_jac, i + 6*k, j      ,     -gsl_matrix_get(id_mat, i, j));
                    }


                }
            }
        }

        //================================================================================
        //Termination condition: if the desired precision is met,
        //the process is terminated.
        //================================================================================
        eps_c  = gsl_blas_dnrm2(err_vec);


        //--------------------------------------------------------------------------------
        // Check that all points are under a given threshold
        //--------------------------------------------------------------------------------
        cout << fname << ". n° " << iter+1 << "/" << itermax << ". nerror = " << eps_c << endl;
        if(eps_c < eps_diff)
        {
            cout << fname << ". Desired precision was reached. " << endl; cout << "nerror = " << eps_c << endl;
            break;
        }


        //--------------------------------------------------------------------------------
        // If a strong convergence is desired, check that we are always decreasing the error
        //--------------------------------------------------------------------------------
        if(strong_conv && eps_c_old < eps_c)
        {
            cout << fname << ". error[n+1] < error[n]. " << endl;
            return GSL_FAILURE;
        }
        eps_c_old = eps_c;


        //================================================================================
        //Compute the correction vector, with minimum norm
        //================================================================================
        status = ftc_corrvec_mn(corr_vec, err_vec, err_jac, nfv, ncs);
        if(status)
        {
            cerr << fname << ". The computation of the correction vector failed."  << endl;
            return GSL_FAILURE;
        }

        //================================================================================
        // Update the free variables
        //================================================================================
        for(int k = 0; k <= n_patch; k++)
        {
            for(int i = 0; i < 6; i++) zmdn[i][k] += gsl_vector_get(corr_vec, 6*k+i);
        }

        //--------------------------------------------------------------------------------
        // Update number of iterations
        //--------------------------------------------------------------------------------
        iter++;
    }

    //====================================================================================
    // Check that the desired condition was met
    //====================================================================================
    if(eps_c > eps_diff)
    {
        cerr << fname << ". The desired precision was not reached."  << endl;
        return GSL_FAILURE;
    }


    //====================================================================================
    // Free
    //====================================================================================
    gslc_matrix_array_free(err_jac_patch, n_patch);
    gsl_vector_free(corr_vec);
    gsl_vector_free(err_vec);
    gsl_matrix_free(err_jac);
    gsl_matrix_free(id_mat);


    return GSL_SUCCESS;
}




//----------------------------------------------------------------------------------------
//  Plot
//----------------------------------------------------------------------------------------
/**
 *  \brief Integrate the state zv[] up to t = t1 on a n_points grid, in either NC or SYS
 *         coordinates, and plot the corresponding result in SYS coordinates
 *         (e.g. EM or SEM coord.).
 *         Then the results are plotted in the xy-plane on a temporary gnuplot window
 *         via the handle *h1.
 *         Print in txt files is included via the integer is_stored.
 **/
int ode_plot_xy(const double zv[], int n_var, double t1, gsl_odeiv2_driver *d,
                gnuplot_ctrl  *h1, int n_points, int color, int is_norm, int is_stored,
                string legend, string filename)
{
    //------------------------------------------------------------------------------------
    // Init
    //------------------------------------------------------------------------------------
    gsl_odeiv2_driver_reset(d);
    set_dir_driver(d, 0.0, t1);

    //System state on a n_points grid
    double x_sys[n_points];
    double y_sys[n_points];

    //Retrieving the parameters
    FBPL* qbp = (FBPL *) d->sys->params;

    //Initial conditions
    double zv_nat[n_var], zv_sys[n_var];
    for(int i=0; i<n_var; i++) zv_nat[i] = zv[i];
    double ti = 0;

    //------------------------------------------------------------------------------------
    //First point
    //------------------------------------------------------------------------------------
    if(is_norm)
    {
        ncsys_m_to_sys_m(0.0, zv_nat, zv_sys, qbp);
        x_sys[0] = zv_sys[0];
        y_sys[0] = zv_sys[1];
    }
    else
    {
        x_sys[0] = zv_nat[0];
        y_sys[0] = zv_nat[1];
    }

    //------------------------------------------------------------------------------------
    //Loop and integration
    //------------------------------------------------------------------------------------
    double t = 0.0;
    set_dir_driver(d, t, t1);

    for(int i =1; i <= n_points; i++)
    {
        ti = i * t1 / n_points;
        gsl_odeiv2_driver_apply (d, &t, ti, zv_nat);

        //Update the sys state
        if(is_norm)
        {
            ncsys_m_to_sys_m(ti, zv_nat, zv_sys, qbp);
            x_sys[i] = zv_sys[0];
            y_sys[i] = zv_sys[1];
        }
        else
        {
            x_sys[i] = zv_nat[0];
            y_sys[i] = zv_nat[1];
        }
    }

    //------------------------------------------------------------------------------------
    //Plotting
    //------------------------------------------------------------------------------------
    gnuplot_setstyle(h1, (char*)"lines");
    gnuplot_set_xlabel(h1, (char*)"X [-]");
    gnuplot_set_ylabel(h1, (char*)"Y [-]");
    gnuplot_plot_xy(h1, x_sys, y_sys, n_points, (char*)legend.c_str(), "lines", "1", "4", color);


    //------------------------------------------------------------------------------------
    //In txt files
    //------------------------------------------------------------------------------------
    if(is_stored)
    gnuplot_fplot_xy(x_sys, y_sys, n_points, filename.c_str());   //orbit

    return GSL_SUCCESS;
}


/**
 *  \brief Integrate the states (zmdn**, tmdn*), in either NC or SYS coordinates,
 *         and plot the corresponding result in SYS coordinates (e.g. EM or SEM coord.).
 *         Then the results are plotted in the xy-plane  on a temporary gnuplot window
 *         via the handle *h1.
 **/
int ode_plot_vec(double **zmdn, double *tmdn, int n_var, int n_patch, gsl_odeiv2_driver *d,
               gnuplot_ctrl  *h1, int n_points, int color, int is_norm, string legend)
{
    //------------------------------------------------------------------------------------
    // Init
    //------------------------------------------------------------------------------------
    gsl_odeiv2_driver_reset(d);

    //EM state on a n_points grid
    double x_sys[n_points];
    double y_sys[n_points];

    //Retrieving the parameters
    FBPL* qbp = (FBPL *) d->sys->params;

    //Initial conditions
    double zv_nat[n_var], zv_sys[n_var];

    //------------------------------------------------------------------------------------
    //Loop and integration
    //------------------------------------------------------------------------------------
    double ts, ti;
    for(int k = 0; k < n_patch; k++)
    {
        for(int i= 0; i < n_var; i++) zv_nat[i] = zmdn[i][k];
        ts  = tmdn[k];

        //Update the SYS state
        if(is_norm)
        {
            ncsys_m_to_sys_m(ts, zv_nat, zv_sys, qbp);
            x_sys[0] = zv_sys[0];
            y_sys[0] = zv_sys[1];
        }
        else
        {
            x_sys[0] = zv_nat[0];
            y_sys[0] = zv_nat[1];
        }

        if(tmdn[n_patch] < 0) d->h = -d->h;
        for(int i = 1; i <= n_points; i++)
        {
            ti = tmdn[k] + 1.0*i/n_points*(tmdn[k+1] - tmdn[k]);

            //Integration
            gsl_odeiv2_driver_apply (d, &ts, ti, zv_nat);

            //Update the SYS state
            if(is_norm)
            {
                ncsys_m_to_sys_m(ts, zv_nat, zv_sys, qbp);
                x_sys[i] = zv_sys[0];
                y_sys[i] = zv_sys[1];
            }
            else
            {
                x_sys[i] = zv_nat[0];
                y_sys[i] = zv_nat[1];
            }
        }
        if(tmdn[n_patch] < 0) d->h = -d->h;

        //--------------------------------------------------------------------------------
        //Plotting
        //--------------------------------------------------------------------------------
        // Of each leg
        if(k == 0) gnuplot_plot_xy(h1, x_sys, y_sys, n_points, (char*)legend.c_str(), "lines", "2", "2", color);
        else gnuplot_plot_xy(h1, x_sys, y_sys, n_points, (char*)"", "lines", "2", "2", color);

        // Of each point
        gnuplot_plot_xy(h1, &x_sys[0], &y_sys[0], 1, (char*)"", "points", "2", "2", color);
    }

    //------------------------------------------------------------------------------------
    // Aesthetics
    //------------------------------------------------------------------------------------
    gnuplot_set_xlabel(h1, (char*)"X [-]");
    gnuplot_set_ylabel(h1, (char*)"Y [-]");

    return GSL_SUCCESS;
}

