#include "pmode.h"

/**
 * \file pmode.cpp
 * \brief Integration and test of the outputs of the parameterization method.
 * \author BLB
 */


//----------------------------------------------------------------------------------------
//
//          Integration of the pm
//
//----------------------------------------------------------------------------------------
/**
 *  \brief Computes the reduce vector field dot(s) = f(s,t)
 *  \param t the current time
 *  \param y the complex state s (CCM), with real and imag parts separated (CCM8 form)
 *  \param f the complex reduced vector field f(s,t), with real and imag parts separated (CCM8 form)
 *  \param params_void a (void) pointer towards a RVF structure that contains:
 *           (1) the Fourier-Taylor expansions of the rvf
 *           (2) the order of the evaluation
 *           (3) the pulsation n of the system
 *           (4) an OFS temp variable.
 *
 *     Note that the use of void *params_void instead of the direct use of an
 *     RVF structure is to comply with GSL ODE integration tools,
 *     that require a vector field of the exact form :
 *          int vf(double t, const double y[], double f[], void *params_void)
 *
 **/
int qbcp_fh(double t, const double y[], double f[], void* params_void)
{
    //-------------------------------
    //Initialization
    //-------------------------------
    RVF* rvf = (RVF*) params_void;

    //-------------------------------
    //Evaluation of the reduced vector field at order rvf->order
    //-------------------------------
    rvf8_eval_cmm8(y, t, rvf->n, rvf->order, rvf->ofs_order, *rvf->fh, *rvf->ofs, f);

    return GSL_SUCCESS;
}


//----------------------------------------------------------------------------------------
//
//          Tests
//
//----------------------------------------------------------------------------------------
/**
 *  \brief Computations of the orbital and invariance errors along an orbit initialized by
 *         the parameterization of the center manifold W(s,t) = PC(t)*Wh(s,t) + V(t), at
 *         order ofts_order for the Taylor series and ofs_order for the Fourier series.
 *         The errors are computed on the interval tvec and plotted.
 **/
int orb_inv_error(const double st0[],      //RCM initial conditions
                  vector<Oftsc>& Wh,       //TFC manifold
                  matrix<Ofsc>& PC,        //COC matrix: z = PC*zh+V
                  vector<Ofsc>& V,         //COC vector: z = PC*zh+V
                  vector<Oftsc>& FW,       //NC vector field
                  vector<Oftsc>& DWf,      //Jacobian of FW
                  vector<Oftsc>& Wdot,     //Partial derivative of W wrt time
                  double tvec[2],          //Integration interval
                  OdeStruct* ode_s_nc,     //ode structure for NC integration
                  OdeStruct* ode_s_rvf,    //ode structure for RVF integration (reduced vector field)
                  int Npoints,             //Number of points on which the errors are estimated
                  FBPL& fbpl,              //current QBCP
                  int ofts_order,          //Order for the eval of the OFTS objects
                  int ofs_order,           //Order for the eval of the OFS objects
                  gnuplot_ctrl**  ht,      //Gnuplot handlers
                  int color)               //Color of the plots
{
    //------------------------------------------------------------------------------------
    // Strings
    //------------------------------------------------------------------------------------
    string F_GS    = fbpl.cs.F_GS;
    string F_PRINT  = fbpl.cs.F_PRINT;
    string F_COC   = fbpl.cs.F_COC;


    //------------------------------------------------------------------------------------
    //Reset integrator
    //------------------------------------------------------------------------------------
    reset_ode_structure(ode_s_nc);
    reset_ode_structure(ode_s_rvf);

    //------------------------------------------------------------------------------------
    //Variables for error computation & plotting
    //------------------------------------------------------------------------------------
    //Time
    double tv[Npoints+1];
    //Orbits (from NC or TFC integration)
    double z_nc[3][Npoints+1], z_nc_from_tfc[3][Npoints+1];
    //Errors and relative energy
    double eOc[Npoints+1], eIc[Npoints+1], rHc[Npoints+1];
    cdouble eId;
    double eIm, eOm, eI[6], eO[6];
    //Energy
    double h_sys, h_li_sys;

    //------------------------------------------------------------------------------------
    // Inner state (NC, TFC...)
    //------------------------------------------------------------------------------------
    double  s1_rcm[REDUCED_NV];    //RCM
    cdouble s1_ccm[REDUCED_NV];    //CCM
    double  s1_ccm8[2*REDUCED_NV]; //CCM8
    double  z1_nc[6], z0_nc[6];   //NC, integrated with NC vector field
    double  z1_nc_from_tfc[6];    //NC, but computed from z = W(s,t)
    double  z1_sys[6], z0_sys[6]; //SYS (either EM or SEM)
    double  f_nc[6];               //NC vector field
    double  n_z1_nc;              //Euclidian norm in NC coordinates
    double  n_z1_km;              //Euclidian norm in km
    Ofsc AUX;                     //OFS temp variable

    //------------------------------------------------------------------------------------
    //Retrieving the pulsation of the QBCP
    //------------------------------------------------------------------------------------
    double n = fbpl.us.n;

    //------------------------------------------------------------------------------------
    //Hamiltonian at the origin s = (0,0,0,0), in SYS coordinates
    //------------------------------------------------------------------------------------
    double st_li[REDUCED_NV];
    for(int i = 0; i < REDUCED_NV; i++) st_li[i] = 0.0;
    rcm_to_nc_by_tfc(st_li, tvec[0], n, ofts_order, ofs_order, Wh, PC, V, z0_nc, false);
    ncsys_m_to_sys_m(tvec[0], z0_nc, z0_sys, (FBPL*) ode_s_nc->d->sys->params);
    h_li_sys = qbcp_H(tvec[0], z0_sys, ode_s_nc->d->sys->params);

    //------------------------------------------------------------------------------------
    // RCM to NC for the initial conditions st0
    //------------------------------------------------------------------------------------
    rcm_to_nc_by_tfc(st0, tvec[0], n, ofts_order, ofs_order, Wh, PC, V, z1_nc, false);

    //------------------------------------------------------------------------------------
    // RCM to CCM8 (real reduced coordinates to complex reduced coordinates)
    //------------------------------------------------------------------------------------
    rcm_to_cmm8(st0, s1_ccm8);

    //------------------------------------------------------------------------------------
    // Euclidian Norm in NC and in km
    //------------------------------------------------------------------------------------
    n_z1_nc = ENorm(z1_nc, 3);
    n_z1_km = n_z1_nc*fbpl.cs.gamma*fbpl.cs.cr3bp.L;

    //------------------------------------------------------------------------------------
    // Print Initial Conditions
    //------------------------------------------------------------------------------------
    //Initial conditions
    printf("-----------------------------------------------------------\n");
    printf("orb_inv_error. With current reduced coordinates and order, \n");
    printf("the initial conditions in NC coordinates are: \n");
    printf("(%3.5f, %3.5f, %3.5f, %3.5f, %3.5f, %3.5f)\n",
           z1_nc[0], z1_nc[1], z1_nc[2], z1_nc[3], z1_nc[4], z1_nc[5]);
    //Euclidian norm in km
    printf("orb_inv_error. Approximated distance from Li: %3.5f km\n", n_z1_km);

    //------------------------------------------------------------------------------------
    // Initial relative Hamiltonian in SYS coordinates
    //------------------------------------------------------------------------------------
    ncsys_m_to_sys_m(tvec[0], z1_nc, z1_sys, (FBPL*) ode_s_nc->d->sys->params);
    h_sys = qbcp_H(tvec[0], z1_sys, ode_s_nc->d->sys->params) ;
    printf("orb_inv_error. Relative Hamiltonian in SYS coordinates: %3.5f\n", h_sys - h_li_sys);

    //------------------------------------------------------------------------------------
    // For plotting (first value)
    //------------------------------------------------------------------------------------
    //Time
    tv[0] = tvec[0];
    //NC
    for(int p = 0; p < 3; p++) z_nc[p][0] = z1_nc[p];
    //NC from manifold
    for(int p = 0; p < 3; p++) z_nc_from_tfc[p][0] = z1_nc[p];
    //Errors
    eIc[0]  = 0.0;
    eOc[0]  = 0.0;
    rHc[0]  = sqrt((h_sys - h_li_sys)*(h_sys - h_li_sys));

    //------------------------------------------------------------------------------------
    // Loop on time
    //------------------------------------------------------------------------------------
    double t  = tvec[0];
    double t2 = tvec[0];
    double ti = 0;
    //Change direction of integration if necessary
    set_dir_ode_structure(ode_s_nc, tvec[0], tvec[1]);
    set_dir_ode_structure(ode_s_rvf, tvec[0], tvec[1]);

    for(int i =0; i<= Npoints; i++)
    {
        //--------------------------------------------------------------------------------
        // Integrate until t = ti
        //--------------------------------------------------------------------------------
        ti = (double) i * (tvec[1] - tvec[0]) / Npoints + tvec[0];
        gsl_odeiv2_driver_apply (ode_s_nc->d, &t, ti, z1_nc);
        gsl_odeiv2_driver_apply(ode_s_rvf->d, &t2, ti, s1_ccm8);

        //--------------------------------------------------------------------------------
        // Semi-analytical computation
        //--------------------------------------------------------------------------------
        //CCM8 to RCM
        ccm8_to_rcm(s1_ccm8, s1_rcm);
        //CCM8 to CCM
        ccm8_to_cmm(s1_ccm8, s1_ccm);
        //z1_nc_from_tfc = W(s1_rcm, ti)
        rcm_to_nc_by_tfc(s1_rcm, ti, n, ofts_order, ofs_order, Wh, PC, V, z1_nc_from_tfc, false);

        //Evaluating the vector field at z1_nc_from_tfc = W(s1_rcm, ti)
        qbcp_vfn_novar(ti, z1_nc_from_tfc, f_nc, ode_s_nc->d->sys->params);

        //Hamiltonian, at the origin
        rcm_to_nc_by_tfc(st_li, ti, n, ofts_order, ofs_order, Wh, PC, V, z0_nc, false);
        ncsys_m_to_sys_m(ti, z0_nc, z0_sys, (FBPL*) ode_s_nc->d->sys->params);
        h_li_sys = qbcp_H(ti, z0_sys, ode_s_nc->d->sys->params);


        //Hamiltonian, taken at z1_nc_from_tfc = W(s1_rcm, ti)
        ncsys_m_to_sys_m(ti, z1_nc_from_tfc, z1_sys, (FBPL*) ode_s_nc->d->sys->params);
        h_sys  = qbcp_H(ti, z1_sys, ode_s_nc->d->sys->params);

        //--------------------------------------------------------------------------------
        //Error computation
        //--------------------------------------------------------------------------------
        for(int p = 0; p < Csts::NV; p++)
        {
            //eO = |z(t) - W(s(t), t)|
            eO[p] = cabs(z1_nc[p] - z1_nc_from_tfc[p]);

            //eI
            //-------------------------------
            //What is the definition of eI?
            // either :
            // - F(W(s,t)) - FW(s,t)
            // - F(W(s,t)) - DWf(s,t) - Wdot(s,t)
            //------------------------------
            eId = f_nc[p]+0.0*I;
            //FW[p].evaluate(s1_ccm, AUX, ofts_order);
            //eId = AUX.evaluate(n*ti);
            //OR
            DWf[p].evaluate(s1_ccm, AUX, ofts_order, ofs_order);
            eId -= AUX.evaluate(n*ti);
            Wdot[p].evaluate(s1_ccm, AUX, ofts_order, ofs_order);
            eId -= AUX.evaluate(n*ti);
            eI[p] = cabs(eId);
        }

        //Taking the maximum (infinity norm)
        eOm = eO[0];
        eIm = eI[0];
        for(int p = 1; p < Csts::NV; p++)
        {
            if(eO[p] > eOm) eOm = eO[p];
            if(eI[p] > eIm) eIm = eI[p];
        }

        //--------------------------------------------------------------------------------
        //Store for plotting
        //--------------------------------------------------------------------------------
        //Time
        tv[i]  = ti;
        //NC
        for(int p = 0; p < 3; p++) z_nc[p][i] = z1_nc[p];
        //NC from manifold
        for(int p = 0; p < 3; p++) z_nc_from_tfc[p][i] = z1_nc_from_tfc[p];
        //Errors
        eOc[i] = eOm;
        eIc[i] = eIm;
        rHc[i] = sqrt((h_sys - h_li_sys)*(h_sys - h_li_sys));
    }

    //------------------------------------------------------------------------------------
    //Define the output name
    //------------------------------------------------------------------------------------
    string str_eO    = orb_inv_error_output_name(F_PRINT, "eO", st0, ofts_order, ofs_order);
    string str_eI    = orb_inv_error_output_name(F_PRINT, "eI", st0, ofts_order, ofs_order);
    string str_rH    = orb_inv_error_output_name(F_PRINT, "rH", st0, ofts_order, ofs_order);
    string str_z_nc  = orb_inv_error_output_name(F_PRINT, "z_nc", st0, ofts_order, ofs_order);
    string str_z_tfc = orb_inv_error_output_name(F_PRINT, "z_nc_from_tfc", st0, ofts_order, ofs_order);


    cout << "orb_inv_error. Outputs are stored in " << F_PRINT+"orbits/" << endl;
    gnuplot_fplot_xy(tv, eOc, Npoints+1, str_eO.c_str()); //eO
    gnuplot_fplot_xy(tv, eIc, Npoints+1, str_eI.c_str()); //eI
    gnuplot_fplot_xy(tv, rHc, Npoints+1, str_rH.c_str()); //eH
    gnuplot_fplot_txyz(tv, z_nc[0], z_nc[1], z_nc[2],   Npoints, str_z_nc.c_str());   //orbit
    gnuplot_fplot_txyz(tv, z_nc_from_tfc[0], z_nc_from_tfc[1], z_nc_from_tfc[2],  Npoints, str_z_tfc.c_str());  //PM orbit


    //------------------------------------------------------------------------------------
    //Plotting
    //------------------------------------------------------------------------------------
    string str_otfs_order  = num_to_string(ofts_order);
    string str_ofs_order   = num_to_string(ofs_order);
    string str_legend      = "Order ("+str_otfs_order+", "+str_ofs_order+")";

    //Logscale if necessary (all errors)
    gnuplot_cmd(ht[0], "set logscale y");
    gnuplot_cmd(ht[0], "set format y \"1e\%%L\"");
    gnuplot_cmd(ht[4], "set logscale y");
    gnuplot_cmd(ht[4], "set format y \"1e\%%L\"");
    //Grid set by default
    for(int i = 0; i <6; i++) gnuplot_cmd(ht[i],  "set grid");

    //Orbital error
    //----------------
    gnuplot_set_xlabel(ht[0], (char*)"t [-]");
    gnuplot_set_ylabel(ht[0], (char*)"eO [-]");
    gnuplot_plot_xy(ht[0], tv, eOc, Npoints+1, (char*)str_legend.c_str(), "lines", "1", "2", color);

    //XY
    //----------------
    gnuplot_set_xlabel(ht[1], (char*)"x [-]");
    gnuplot_set_ylabel(ht[1], (char*)"y [-]");
    gnuplot_plot_xy(ht[1], z_nc[0], z_nc[1], Npoints+1, (char*)str_legend.c_str(), "lines", "1", "2", color);
    gnuplot_plot_xy(ht[1], z_nc_from_tfc[0], z_nc_from_tfc[1], Npoints+1, (char*)str_legend.c_str(), "lines", "dashed", "2", color);

    gnuplot_plot_xy(ht[1], &z_nc[0][0], &z_nc[1][0], 1, "", "points", "1", "2", color);

    //XZ
    //----------------
    gnuplot_set_xlabel(ht[2], (char*)"x [-]");
    gnuplot_set_ylabel(ht[2], (char*)"z [-]");
    gnuplot_plot_xy(ht[2], z_nc[0], z_nc[2], Npoints+1, (char*)str_legend.c_str(), "lines", "1", "2", color);
    gnuplot_plot_xy(ht[2], z_nc_from_tfc[0], z_nc_from_tfc[2], Npoints+1, (char*)str_legend.c_str(), "lines", "dashed", "2", color);

    //YZ
    //----------------
    gnuplot_set_xlabel(ht[3], (char*)"y [-]");
    gnuplot_set_ylabel(ht[3], (char*)"z [-]");
    gnuplot_plot_xy(ht[3], z_nc[1], z_nc[2], Npoints+1, (char*)str_legend.c_str(), "lines", "1", "2", color);
    gnuplot_plot_xy(ht[3], z_nc_from_tfc[1], z_nc_from_tfc[2], Npoints+1, (char*)str_legend.c_str(), "lines", "dashed", "2", color);


    //Invariance error
    //----------------
    gnuplot_set_xlabel(ht[4], (char*)"t [-]");
    gnuplot_set_ylabel(ht[4], (char*)"eI [-]");
    gnuplot_plot_xy(ht[4], tv, eIc, Npoints+1, (char*)str_legend.c_str(), "lines", "1", "2", color);


    //Relative Hamiltonian
    //----------------
    gnuplot_set_xlabel(ht[5], (char*)"t [-]");
    gnuplot_set_ylabel(ht[5], (char*)"rH [-]");
    gnuplot_plot_xy(ht[5], tv, rHc,  Npoints+1, (char*)str_legend.c_str(), "lines", "1", "3", color);


    //Titles
    //----------------
    gnuplot_cmd(ht[0], "set title \"Orbital error");
    gnuplot_cmd(ht[1], "set title \"Orbit in XY plane in NC coordinates");
    gnuplot_cmd(ht[2], "set title \"Orbit in XZ plane in NC coordinates");
    gnuplot_cmd(ht[3], "set title \"Orbit in YZ plane in NC coordinates");
    gnuplot_cmd(ht[4], "set title \"Invariance error");
    gnuplot_cmd(ht[5], "set title \"Relative Hamiltonian: H(t) - H_Li(t)");

    return GSL_SUCCESS;
}



/**
 *  \brief Test of the parameterization of the central manifold on a
 *         given orbit, through the computation of various errors (eO, eI) and the
 *         relative energy.
 *              - The orbit is initialized by the reduced state vector si[].
 *              - The errors are computed on the interval tvec.
 *              - The errors are computed for the n_ofts_order Taylor orders stored in v_ofts_order.
 *              - The errors are computed for the n_ofs_order Fourier orders stored in v_ofs_order.
 *
 *   Note that: the expansions FW, DWf and Wdot are taken from files,
 *   whereas the parameterization itself CM and CMh, and the reduced vector field Fh
 *   are taken from global objects, defined in init.cpp.
 *
 *   Requires init_inv_man().
 *
 **/
void pm_error_vs_orders_test(int n_ofts_order, int v_ofts_order[],
                             int n_ofs_order, int v_ofs_order[],
                             double si[], double tvec[])
{
    //------------------------------------------------------------------------------------
    // Strings
    //------------------------------------------------------------------------------------
    string F_GS    = SEML.cs.F_PMS;
    string F_PRINT  = SEML.cs.F_PRINT;
    string F_COC   = SEML.cs.F_COC;


    //------------------------------------------------------------------------------------
    // Initialization of additional expansions needed to compute the invariance error
    //------------------------------------------------------------------------------------
    vector<Oftsc> FW(6);
    vector<Oftsc> DWf(6);
    vector<Oftsc> Wdot(6);

    read_vofts_bin(FW,   F_GS+"FW/C_FW",       OFS_ORDER);
    read_vofts_bin(DWf,  F_GS+"DWf/C_DWf",     OFS_ORDER);
    read_vofts_bin(Wdot, F_GS+"Wdot/C_Wdot",   OFS_ORDER);


    //------------------------------------------------------------------------------------
    //Integration tools
    //------------------------------------------------------------------------------------
    // ODE structure for integration in NC coordinates
    OdeStruct ode_s_nc;
    init_ode_structure(&ode_s_nc, gsl_odeiv2_step_rk8pd, gsl_root_fsolver_brent, 6, qbcp_vfn_novar, &SEML);


    // RVF structure for integration in reduced coordinates
    RVF rvf;
    rvf.ofs_order = SEML.eff_nf;
    Ofsc AUX(rvf.ofs_order);
    rvf.fh        = &Fh;
    rvf.ofs       = &AUX;
    rvf.order     = OFTS_ORDER;
    rvf.n         = SEML.us.n;

    // ODE structure for integration in reduced coordinates
    OdeStruct ode_s_rvf;
    init_ode_structure(&ode_s_rvf, gsl_odeiv2_step_rk8pd, gsl_root_fsolver_brent,2*REDUCED_NV, qbcp_fh, &rvf);

    cout << "---------------------------------------------------" << endl;
    cout << "                                                   " << endl;
    cout << "               Test of the obtained PM             " << endl;
    cout << "                                                   " << endl;
    cout << "---------------------------------------------------" << endl;
    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(15);

    //------------------------------------------------------------------------------------
    //Misc variables
    //------------------------------------------------------------------------------------
    //Int variables
    double st0[REDUCED_NV];
    //temp
    Ofsc BUX(OFS_ORDER);


    //------------------------------------------------------------------------------------
    // Orders of the expansions
    //------------------------------------------------------------------------------------
    int ofts_order = OFTS_ORDER;
    int ofs_order  = OFS_ORDER;

    //------------------------------------------------------------------------------------
    // Initial conditions
    //------------------------------------------------------------------------------------
    string st0s = "";
    for(int p = 0; p < REDUCED_NV; p++) st0[p] = si[p];

    //------------------------------------------------------------------------------------
    //For plotting
    //------------------------------------------------------------------------------------
    //Gnuplot handlers
    gnuplot_ctrl**  ht = (gnuplot_ctrl**) calloc(6, sizeof(gnuplot_ctrl*));
    for(int i = 0; i <6; i++) ht[i] = gnuplot_init();

    //------------------------------------------------------------------------------------
    //Loop on orders
    //------------------------------------------------------------------------------------
    int color = 1;
    for(int ii = 0; ii< n_ofts_order; ii++)
    {
        //Update ofts_order
        ofts_order = min(v_ofts_order[ii], OFTS_ORDER);

        for(int jj = 0; jj< n_ofs_order; jj++)
        {
            //Update ofs_order
            ofs_order = min(v_ofs_order[jj], OFS_ORDER);

            // Update orders in rvf structure
            rvf.order     = ofts_order;
            rvf.ofs_order = ofs_order;

            // eO, eI, and orbit
            tic();
            orb_inv_error(st0, CMh, Mcoc, Vcoc, FW, DWf, Wdot, tvec, &ode_s_nc,
                          &ode_s_rvf, 1000, SEML, ofts_order, ofs_order, ht, color++);
            cout << "orb_inv_error in: " << toc() << endl;
        }
    }

    char ch;
    printf("Press ENTER to close the gnuplot window(s)\n");
    scanf("%c",&ch);


    //Closing handlers
    for(int i = 0; i <6; i++) gnuplot_close(ht[i]);
    free(ht);


    free_ode_structure(&ode_s_nc);
    free_ode_structure(&ode_s_rvf);

}

/**
 *  \brief Computations of the filenames used in the routine orb_inv_error (orbital error, invariance error, etc).
 **/
string orb_inv_error_output_name(string f_plot, string title, const double *st0, int ofts_order, int ofs_order)
{
    string str_p, str_otfs_order, str_ofs_order, str_st;

    str_otfs_order  = num_to_string(ofts_order);
    str_ofs_order   = num_to_string(ofs_order);

    string filename = (f_plot+"orbits/"+title+"_ofts_order_"+str_otfs_order+"_ofs_order_"+str_ofs_order);

    for(int p = 0; p < REDUCED_NV; p++)
    {
        str_st = num_to_string(st0[p]);
        str_p  = num_to_string(p+1);

        filename += "_s"+str_p+"_"+str_st;
    }
    filename += ".txt";

    return filename;
}

//----------------------------------------------------------------------------------------
//
//          DEPRECATED (to be cleaned up)
//
//----------------------------------------------------------------------------------------
/**
 *  \brief Evaluates the contributions of each order in W to the computation of W(s,t), with an arbitrary state (s,t)
 **/
void pmContributions()
{
    //------------------------------------------------------------------------------------
    // Strings
    //------------------------------------------------------------------------------------
    string F_GS    = SEML.cs.F_GS;
    string F_PRINT  = SEML.cs.F_PRINT;
    string F_COC   = SEML.cs.F_COC;

    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(5);

    //------------------------------------------------------------------------------------
    // Initialization of the NC manifolds
    //------------------------------------------------------------------------------------
    vector<Oftsc> W(6);
    vector<Oftsc> Wh(6);
    read_vofts_txt(Wh, F_GS+"W/Wh", OFS_ORDER);
    read_vofts_txt(W,  F_GS+"W/W", OFS_ORDER);


    //------------------------------------------------------------------------------------
    //Contribution of each order to the evaluation
    //------------------------------------------------------------------------------------
    // Initial conditions
    //------------------------------------------------------------------------------------
    string st0s;
    cdouble st0[4], s0[4];
    for(int p = 0; p < 4; p++) st0[p] = 5e-0+0.0*I;
    st0s = "5e-0";

    //------------------------------------------------------------------------------------
    //"Realification": s0 = REAL(st0)
    //------------------------------------------------------------------------------------
    s0[0] = 1.0/sqrt(2)*( st0[0]   - st0[2]*I);
    s0[2] = 1.0/sqrt(2)*( st0[2]   - st0[0]*I);
    s0[1] = 1.0/sqrt(2)*( st0[1]   - st0[3]*I);
    s0[3] = 1.0/sqrt(2)*( st0[3]   - st0[1]*I);

    //------------------------------------------------------------------------------------
    // Evaluate W
    //------------------------------------------------------------------------------------
    //z0t = W(s0), z0 = z0t(0.0)
    Ofsc AUX(OFS_ORDER);
    for(int k = 0; k <= OFTS_ORDER; k++)
    {
        cout << k << " ";
        for(int p = 0; p < 6; p++)
        {
            W[p].contribution(s0, AUX, k);
            cout << AUX.l1norm() << " ";
        }
        cout << endl;
    }
}

/**
 *  \brief Evaluates the l1/linf norm of each order in W.
 **/
void pmNorms()
{
    //------------------------------------------------------------------------------------
    // Strings
    //------------------------------------------------------------------------------------
    string F_GS    = SEML.cs.F_GS;
    string F_PRINT  = SEML.cs.F_PRINT;
    string F_COC   = SEML.cs.F_COC;


    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(5);
    //------------------------------------------------------------------------------------
    // Initialization of the NC manifolds
    //------------------------------------------------------------------------------------
    vector<Oftsc> W(6);
    vector<Oftsc> Wh(6);
    read_vofts_bin(Wh, F_GS+"W/Wh", OFS_ORDER);
    read_vofts_bin(W,  F_GS+"W/W", OFS_ORDER);

    //------------------------------------------------------------------------------------
    //L1-norm
    //------------------------------------------------------------------------------------
    gnuplot_ctrl*  h5;
    h5 = gnuplot_init();
    gnuplot_cmd(h5, "set logscale y");
    gnuplot_cmd(h5, "set format y \"1e\%%L\"");
    gnuplot_cmd(h5, "set grid");
    gnuplot_cmd(h5, "set title \"l_1(order)\" ");


    double l1n[OFTS_ORDER];
    double l1nMax[OFTS_ORDER+1];
    double kc[OFTS_ORDER+1];
    cout << "---------------------------"<< endl;
    cout << "Linf-norm of W:" << endl;

    for(int k = 0; k <= OFTS_ORDER; k++)
    {
        kc[k] = k;
        l1nMax[k] = 0.0;
        for(int p = 0; p < Csts::NV; p++)
        {
            //l1n[k] = Wh[p].linfnorm(k);
            //l1n[k] = W[p].l1norm(k);
            l1n[k] = W[p].linfnorm(k);
            if(l1n[k] > l1nMax[k]) l1nMax[k] = l1n[k];
        }
        cout << k << " " << l1nMax[k] << endl;
    }

    gnuplot_plot_xy(h5, kc, l1nMax, OFTS_ORDER+1, (char*)"", "points", "1", "2", 1);

    char ch;
    printf("Press ENTER to close the gnuplot window(s)\n");
    scanf("%c",&ch);
}

/**
 *  \brief Small divisors under a certain value
 **/
void pmSmallDivisors(double sdmax)
{
    //------------------------------------------------------------------------------------
    // Strings
    //------------------------------------------------------------------------------------
    string F_NF    = SEML.cs.F_NF;
    string F_PRINT  = SEML.cs.F_PRINT;
    string F_COC   = SEML.cs.F_COC;


    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(5);
    //------------------------------------------------------------------------------------
    // Initialization of the NC manifolds
    //------------------------------------------------------------------------------------
    vector<Oftsc> smallDiv(6);
    read_vofts_txt(smallDiv, F_NF+"W/smallDiv", OFS_ORDER);

    //------------------------------------------------------------------------------------
    //Small divisors
    //------------------------------------------------------------------------------------
    cout << "---------------------------"<< endl;
    cout << "Number of small divisors other than order 0:" << endl;
    for(int k = 2; k <= OFTS_ORDER; k++)
    {
        cout << k << " ";
        for(int p = 0; p < Csts::NV; p++)
        {
            cout <<  smallDiv[p].nsd(k, OFS_ORDER, sdmax) - smallDiv[p].nsd(k, 0, sdmax) << " ";

        }
        cout << endl;

    }



    cout << "---------------------------"<< endl;
    cout << "Number of small divisors at order 0:" << endl;
    for(int k = 2; k <= OFTS_ORDER; k++)
    {
        cout << k << " ";
        for(int p = 0; p < Csts::NV; p++)
        {
            cout <<  smallDiv[p].nsd(k, 0, sdmax) << " ";

        }
        cout << endl;

    }
}

/**
 *  \brief Evaluates the variations of the initial conditions (IC) wrt to the order in W(s,t), with an arbitrary state (s,t)
 **/
void pmTestIC()
{
    //------------------------------------------------------------------------------------
    // Initialization of the invariant manifold
    //------------------------------------------------------------------------------------
    init_inv_man(SEML);

    //------------------------------------------------------------------------------------
    // Strings
    //------------------------------------------------------------------------------------
    string F_GS    = SEML.cs.F_GS;
    string F_PRINT  = SEML.cs.F_PRINT;
    string F_COC   = SEML.cs.F_COC;


    //Initialization of the configuration
    double st0[4];
    for(int p = 0; p < 4; p++) st0[p] = 0.0;
    st0[0] = 4;

    //"Realification": Xeval = REAL(st0)
    cdouble Xeval[4];
    Xeval[0] = 1.0/sqrt(2)*( st0[0]   - st0[2]*I);
    Xeval[2] = 1.0/sqrt(2)*( st0[2]   - st0[0]*I);
    Xeval[1] = 1.0/sqrt(2)*( st0[1]   - st0[3]*I);
    Xeval[3] = 1.0/sqrt(2)*( st0[3]   - st0[1]*I);

    //------------------------------------------------------------------------------------
    // Evaluate Wc = |W(Xeval, 0)|
    //------------------------------------------------------------------------------------
    vector<Ofsc> z0t(6);
    double z0[4][6];

    //Keymap
    int keyMap[4];
    int m;
    keyMap[0] = 2;
    keyMap[1] = 5;
    keyMap[2] = 10;
    keyMap[3] = OFTS_ORDER;

    for(int i = 0; i < 4; i++)
    {
        m = keyMap[i];
        for(int p = 0; p < 6; p++)
        {
            CM[p].evaluate(Xeval, z0t[p], m, OFS_ORDER);
            z0[i][p] = creal(z0t[p].evaluate(0.0));
        }
    }

    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(5);
    cout << "IC " << m << ":" << endl;
    cout << "------------------------------" << endl;
    cout << " 2-5        5-10         10-15 " << endl;
    for(int p = 0; p <6; p++)
    {
        for(int i = 0; i < 3; i++)cout << z0[i][p]-z0[i+1][p] << "   ";
        cout << endl;
    }

}

