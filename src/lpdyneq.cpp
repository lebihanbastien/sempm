#include "lpdyneq.h"

/**
 * \file lpdyneq_single_shooting.cpp
 * \brief Computation of the dynamical equivalent to the libration points for various models.
 * \author BLB.
 */


//========================================================================================
//         Dynamical equivalents to the Libration points
//========================================================================================
/**
 *  \brief Computation of the dynamical equivalents of the libration points. Test function.
 *
 *         - This function makes use of the subroutine lpdyneq_single_shooting.
 *           The latter contains the initialization of the first guess
 *           (the geometrical positions of the CRTBP libration points)
 *           and the differential corrector (single shooting).
 *         - After lpdyneq_single_shooting, the resulting initial conditions are
 *           integrated and plotted on a full orbit.
 *         - Then, a test of the periodicity of the orbit is performed, via the
 *           computation of the error |zv(0) - zv(T)|.
 *         - Finally, the dynamical equivalent is computed again via a multiple shooting
 *           approach, for consistency.
 *         - If is_stored is true, the results (x(t), y(t)) in synodical coordinates are
 *           stored in txt files of the form: "./plot/QBCP/DYNEQ_QBCP_EM_L1_NC.txt".
 **/
void compute_dyn_eq_lib_point(FBPL &fbpl, int is_stored)
{
    //====================================================================================
    // 1. Initialization
    //====================================================================================
    //------------------------------------------------------------------------------------
    //Plotting devices
    //------------------------------------------------------------------------------------
    gnuplot_ctrl  *h1 = gnuplot_init();

    //------------------------------------------------------------------------------------
    // Integration tools
    // Initial vector field include nonlinear variationnal equations
    // to get the periodic orbit
    //------------------------------------------------------------------------------------
    OdeStruct ode_s;
    init_ode_structure(&ode_s, gsl_odeiv2_step_rk8pd, gsl_root_fsolver_brent, 42,
                       fbpl.is_norm? qbcp_vfn_varnonlin:qbcp_vf, &fbpl);



    //====================================================================================
    // 2. Display
    //====================================================================================
    //Building the filename
    string f_model  = init_F_MODEL(fbpl.model);
    string f_csys   = init_F_COORDSYS(fbpl.coord_sys);
    string f_li     = init_F_LI(fbpl.li);
    string f_norm   = fbpl.is_norm?"_NC":"";
    string filename = "plot/"+f_model+"DYNEQ_"+f_model+"_"+f_csys+"_"+f_li+f_norm+".txt";

    cout << endl;
    cout << "===============================================" << endl;
    cout << "Computation of a dynamical equivalent of a     " << endl;
    cout << "libration point on a full period of the system." << endl;
    cout << "Characteristics:                               " << endl;
    cout << " - Model:  " << f_model << endl;
    cout << " - Number: " << f_li << endl;
    cout << " - Coordinates: " << f_norm+f_csys << endl;
    cout << " - Saved in: " << filename << endl;
    cout << "===============================================" << endl;

    //====================================================================================
    // 3. Single shooting to get the dynamical equivalent of the libration point (lpdyneq_single_shooting)
    //====================================================================================
    cout << "First, a single shooting procedure is used.    " << endl;
    cout << "===============================================" << endl;
    double z0v[42];
    lpdyneq_single_shooting(ode_s.d, z0v);

    //====================================================================================
    // 4. Plot & Print in file if necessary
    //====================================================================================
    int Npoints = 500;

    //Building the solution and printing in file
    ode_plot_xy(z0v, 42, fbpl.us.T, ode_s.d, h1, Npoints, 6, fbpl.is_norm, is_stored, "From single shooting", filename.c_str());

    //====================================================================================
    // 5. Periodicity condition for single shooting
    //====================================================================================
    periodicity_condition(z0v, 42, 6, fbpl.us.T, ode_s.d, fbpl.is_norm);

    //====================================================================================
    // 6. Multiple shooting
    //====================================================================================
    int M = 16;
    cout << "Second, a multiple shooting procedure is used,  " << endl;
    cout << "with " << M << " patch points.                  " << endl;
    cout << " Note that the multiple shooting procedure      " << endl;
    cout << " is used here only for the sake of completeness." << endl;
    cout << " In particular, it is not used in nfo2.cpp      " << endl;
    cout << " since it did not improve the results.          " << endl;
    cout << "=============================================== " << endl;
    double **zmd  = dmatrix(0, 41, 0, M);
    double *tmd   = dvector(0, M);

    //Building the patch points
    lpdyneq_patch_points(z0v, ode_s.d, fbpl.us.T, zmd, tmd, M);

    //Multiple shooting procedure
    mult_shoot_period(zmd, tmd, zmd, tmd, ode_s.d, M, 1e-14, false, h1, false);

    //Plot
    ode_plot_vec(zmd, tmd, 42, M, ode_s.d, h1, Npoints, 3, true, "From multiple shooting");

    //====================================================================================
    // 7. Close window
    //====================================================================================
    pressEnter(true, "Press ENTER to close the gnuplot window(s)\n");
    gnuplot_close(h1);
}


//========================================================================================
// Subroutines
//========================================================================================
/**
 *  \brief Main routine for the computation of the dynamical equivalent to the libration
 *         points via a single shooting procedure.
 *         The results are given in the form of the initial conditions (42 dimensions)
 *         updated in the array z0v[].
 **/
void lpdyneq_single_shooting(gsl_odeiv2_driver *d, double z0v[])
{
    //====================================================================================
    // 0. Misc init.
    //====================================================================================
    //Settings for iostd
    Config::configManager().coutlp();

    //Retrieving the parameters
    FBPL* qbp    = (FBPL *) d->sys->params;
    int is_norm  =  qbp->is_norm;
    int li       =  qbp->cs.li;
    double tend  =  qbp->us.T;

    //====================================================================================
    // 1. Initial guess
    //====================================================================================
    double z1v[6];

    //------------------------------------------------------------------------------------
    //Approximated position from CR3BP
    //------------------------------------------------------------------------------------
    if(is_norm)
    {
        //Approximated position from CR3BP: null vector
        for(int i =0; i<6; i++) z0v[i] = 0.0;
    }
    else
    {
        //Approximated position from CR3BP
        for(int i =0; i<6; i++) z1v[i] = 0.0;
        switch(li)
        {
        case 1:
            z1v[0] = -SEML.cs.cr3bp.l1.position[0];
            break;
        case 2:
            z1v[0] = -SEML.cs.cr3bp.l2.position[0];
            break;
        case 3:
            z1v[0] = -SEML.cs.cr3bp.l3.position[0];
            break;
        }

        switch(SEML.coord_sys)
        {
        case Csts::EM:
        {
            //To momenta
            em_v_to_em_m(0.0, z1v, z0v, qbp);
            break;
        }

        case Csts::SEM:
        {
            //To momenta
            se_v_to_se_m(0.0, z1v, z0v, qbp);
            break;
        }

        }
    }

    //------------------------------------------------------------------------------------
    // STM is appened to z0v (STM = eye(6) at t = t0)
    //------------------------------------------------------------------------------------
    //Identity matrix eye(6)
    gsl_matrix *Id = gsl_matrix_calloc(6,6);
    gsl_matrix_set_identity (Id);

    //Storing eye(6) into the initial vector
    gslc_mat_to_vec(z0v, Id, 6, 6, 6);

    //====================================================================================
    // 2. Differential correction scheme
    //====================================================================================
    //Initial IC
    cout << "lpdyneq_single_shooting. Initial state:" << endl;
    printf("(%3.5e, %3.5e, %3.5e, %3.5e, %3.5e, %3.5e)\n",
           z0v[0], z0v[1], z0v[2], z0v[3], z0v[4], z0v[5]);

    double prec = 2e-14;
    cout << "lpdyneq_single_shooting. Starting differential correction..." << endl;
    single_shoot_sym(z0v, 0.5*tend, prec, d, 0);
    cout << "lpdyneq_single_shooting. End of differential correction." << endl;

    //====================================================================================
    // 3. Print result
    //====================================================================================
    //Final IC
    cout << "lpdyneq_single_shooting. Final state:" << endl;
    printf("(%3.5f, %3.5e, %3.5e, %3.5e, %3.5e, %3.5e)\n",
           z0v[0], z0v[1], z0v[2], z0v[3], z0v[4], z0v[5]);
}



/**
 *   \brief Computes M+1 patch points, indexed between 0 and M,
 *          along the T-periodic orbit starting at z0v = [x0, y0, z0, vx0, vy0, vz0, eye(6)]
 *          Each STM is initialized with the identity matrix.
 **/
void lpdyneq_patch_points(const double *z0v, gsl_odeiv2_driver *d, double t_period,
                          double **zm, double *tm, int M)
{
    //====================================================================================
    //Initial conditions + eye(6)
    //====================================================================================
    double zv[42], t = 0.0;
    for(int i = 0; i < 42; i++) zv[i] = z0v[i];

    //First patch point is z0v
    for(int i = 0; i < 42; i++) zm[i][0] = z0v[i];
    tm[0] = 0.0;

    //====================================================================================
    //Loop on the grid [1...M]
    //====================================================================================
    for(int k = 1; k <= M; k++)
    {
        //Time
        tm[k] = 1.0*k*t_period/M;
        gsl_odeiv2_driver_reset(d);

        //Integration
        gsl_odeiv2_driver_apply (d, &t, tm[k], zv);

        //Storage in zm - each STM is initialized with the identity matrix
        for(int i= 0; i < 6;  i++)  zm[i][k] = zv[i];
        for(int i= 6; i < 42;  i++) zm[i][k] = z0v[i];

        //The STM is restarded: STM(tk) = eye(6);
        for(int j=6; j<42; j++) zv[j] = z0v[j];
    }
}



//========================================================================================
// Periodicity Condition
//========================================================================================
/**
 *  \brief Test the periodicity condition zv(0) = zv(t_period), with initial condition z0v,
 *         for the vector field contained in the driver d.
 *
 *  Note that the integer N is the number of variables associated to the driver d,
 *  and should be also the number of variables in zv.
 *  However, the periodicity condition is tested only on n_var_test variables (a usual
 *  example is n_var = 42 but n_var_test = 6).
 **/
int periodicity_condition(const double zv0[], int n_var, int n_var_test, double t_period,
                          gsl_odeiv2_driver *d, int is_norm)
{
    //------------------------------------------------------------------------------------
    // Init
    //------------------------------------------------------------------------------------
    gsl_odeiv2_driver_reset(d);

    //Retrieving the parameters
    FBPL* qbp = (FBPL *) d->sys->params;

    //Initial conditions
    double zv_nat[n_var], zv0_sys[n_var], zv1_sys[n_var], zv0_nat[n_var], zv1_nat[n_var];
    for(int i=0; i< n_var; i++) zv_nat[i] = zv0[i];

    //First point in system coordinates
    if(is_norm) ncsys_m_to_sys_m(0.0, zv_nat, zv0_sys, qbp);
    else for(int i=0; i< n_var; i++) zv0_sys[i] = zv_nat[i];

    //First point in native coordinates
    for(int i=0; i< n_var; i++) zv0_nat[i] = zv_nat[i];

    //------------------------------------------------------------------------------------
    //Integration
    //------------------------------------------------------------------------------------
    double t = 0.0;
    set_dir_driver(d, t, t_period);
    gsl_odeiv2_driver_apply (d, &t, t_period, zv_nat);

    //------------------------------------------------------------------------------------
    //Final state
    //------------------------------------------------------------------------------------
    // In system coordinates
    if(is_norm) ncsys_m_to_sys_m(t, zv_nat, zv1_sys, qbp);
    else for(int i=0; i< n_var; i++) zv1_sys[i] = zv_nat[i];

    //In native coordinates
    for(int i=0; i< n_var; i++) zv1_nat[i] = zv_nat[i];

    //------------------------------------------------------------------------------------
    //Comparison
    //------------------------------------------------------------------------------------
    // Norm |zv(0) - zv(T)|
    double zvn_sys = 0.0, zvn_nat = 0.0;
    for(int i = 0; i < n_var_test; i++) zvn_sys += (zv0_sys[i] - zv1_sys[i])*(zv0_sys[i] - zv1_sys[i]);
    for(int i = 0; i < n_var_test; i++) zvn_nat += (zv0_nat[i] - zv1_nat[i])*(zv0_nat[i] - zv1_nat[i]);
    zvn_sys = sqrt(zvn_sys);
    zvn_nat = sqrt(zvn_nat);

    cout << "===============================================" << endl;
    cout << "Periodicity condition after single shooting:   " << endl;
    cout << "- |zv(0) - zv(T)| in system coordinates = "      << zvn_sys << endl;
    cout << "- |zv(0) - zv(T)| in native coordinates = "      << zvn_nat << endl;
    cout << "===============================================" << endl;

    return GSL_SUCCESS;
}

