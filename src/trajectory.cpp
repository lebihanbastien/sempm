#include "pm.h"
#include "poincare.h"
#include "trajectory.h"
#include "config.h"

// €€ TODO:
// 0. Faire des fichiers de configuration propre pour TOUTES les variétés (RTBP/QBCP, SEM/EM, L1,L2, L3) + fichiers de test (E0, EI, EH + projection!!)
// 01. Verifier l'absence de petits diviseurs tant que le mixed style ne concerne que la séparation centre-hyperbolique. Pour cela, remettre en place la gestion des
//       petits diviseurs par un objet FT dédié dans pmt.
// 02. Verifier les COC EM --> SEM? Est-il normal de voir les variétés issus de EML2 "tourner" autour à ce point ? Oui probablement.
//      Dans Parker par exemple (fig 17, Shoot the Moon 3D), l'ensemble des branches de la variété issue de la Lune ont été démarrées au même instant.
//      Ce qui n'est PAS le cas pour nous, car tout dépend du temps ici.
// 1. Remettre au propre ce fichier. En particulier, commentaire et trajectory_CU.
// 2. Automatiser certaines tâches: plot des principaux points dans les deux repères...
// 3. Verifier le principe de projection en SEML2. Voir si on obtient bien Wh = Wh_proj pour les directions qui sont linéaires en les variables si
//      Peut-être que le gap observé actuellement est dû au fait que tout est mélangé dans W par rapport à Wh ?
//      Probablement: la différence en termes de vitesse (ou moments) est telle qu'elle est compensée par un écart énorme en position
// 4. Etape suivante: pour une orbite donnée (s0, t0) en EML2, on construit la variété instable + une carte des distance à la variété centrale de SEML2
// 5. En fait, ce qui nous intéresse, c'est la distance à la variété centrale-stable de SEML2. Peut-on envisager une projection sur elle ? Revoir les notes sur le passage
//      Center --> Center-Stable.


void NCEMmtoSEMm_vec(double **yNCEM, double *tNCEM, double **ySEM, int N, QBCP_L *qbcp_l)
{
    double yNCd[6], ySEMd[6];
    //Loop on all elements in yNCEM
    for(int p = 0; p <= N; p++)
    {
        //Copy step p in yNCd
        for(int k = 0; k <6; k++) yNCd[k] = yNCEM[k][p];
        //NC EM to SEM
        NCEMmtoSEMm(tNCEM[p], yNCd, ySEMd, qbcp_l);
        //Copy result in ySEM
        for(int k = 0; k <6; k++) ySEM[k][p] = ySEMd[k];
    }
}

void NCSEMmtoSEMm_vec(double **yNCSEM, double *tNCSEM, double **ySEM, int N, QBCP_L *qbcp_l)
{
    double yNCSd[6], ySEMd[6];

    //Loop on all elements in yNCEM
    for(int p = 0; p <= N; p++)
    {
        //Copy step p in yNCSd
        for(int k = 0; k <6; k++) yNCSd[k] = yNCSEM[k][p];
        //NC SEM to SEM
        NCtoSEM(tNCSEM[p], yNCSd, ySEMd, qbcp_l);
        //Copy result in ySEM
        for(int k = 0; k <6; k++) ySEM[k][p] = ySEMd[k];
    }
}

void NCprojCCMtoCUS(double *yv, double tv, double sti[5], Orbit &orbit, double epsilon)
{
    //Projection tools
    double omega1 = cimag(Fh[0].getCoef(1,0)->ofs_getCoef(0));
    double omega3 = cimag(Fh[1].getCoef(1,1)->ofs_getCoef(0));
    cdouble scp[REDUCED_NV];

    //Get the closest point on the center manifold, in scp[4]
    NCprojCCM(yv, tv, SEML.us.n, OFS_ORDER, MIcoc, Vcoc, omega1, omega3, scp, REDUCED_NV);
    //Get the correspondance in RCM coordinates
    CCMtoRCM(scp, sti, REDUCED_NV);
    //Add a given quantity on the hyperbolic direction
    sti[4] = epsilon;
}

/**
 *   \brief Integrates a given trajectory with IC in sti[]
 **/
void trajectory_CU(double sti[], Pmap &pmap, int label, int is3D)
{
    char ch;
    //-------------------------------------------------
    //Integration tools
    //-------------------------------------------------
    //------------------
    //For dot(z) = F(z)
    //------------------
    OdeStruct ode_s_6;
    init_ode_NC(ode_s_6, pmap);
    //For root finding
    OdeStruct ode_s_6_root;
    init_ode_NC(ode_s_6_root, pmap);

    //------------------
    //For dot(s) = fh(s)
    //------------------
    //Reduced vector field structure
    RVF rvf;
    OdeStruct ode_s_8;
    Ofsc rvf_ofs(OFS_ORDER);
    init_ode_CCM(ode_s_8, rvf, rvf_ofs, pmap);
    //Reduced vector field structure for root finding
    RVF rvf_root;
    OdeStruct ode_s_8_root;
    Ofsc rvf_ofs_root(OFS_ORDER);
    init_ode_CCM(ode_s_8_root, rvf_root, rvf_ofs_root, pmap);

    //------------------
    //Event structures
    //------------------
    value_params val_par;
    val_par.dim        = 2;                   //event on z component
    val_par.direction  = 0;                   //all zeros are detected
    val_par.max_events = pmap.max_events;     //maximum of events
    val_par.value      = 0.0;                 //event when z = 0
    value_function fvalue;
    fvalue.val_par     = &val_par;
    fvalue.value       = &linear_intersection;

    //------------------
    //Orbit structure
    //------------------
    Ofsc orbit_ofs(OFS_ORDER);
    Orbit orbit;
    init_orbit(&orbit, &CM, &CMh, &JCM, &Mcoc, &Vcoc, &val_par,
               &fvalue, &ode_s_6, &ode_s_8, &ode_s_6_root, &ode_s_8_root,
               &pmap, &SEML, &orbit_ofs, pmap.vdim, label);


    //------------------------------------------
    //Update the IC
    //------------------------------------------
    orbit_update_ic(orbit, sti, orbit.pmap->t0);
    orbit.int_method =  0;


    //------------------------------------------
    //Integration
    //------------------------------------------
    string filename;
    filename = SEML.cs.F_PLOT+"traj_IC_";
    filename += numTostring(sti[0])+"_";
    filename += numTostring(sti[1])+"_";
    filename += numTostring(sti[2])+"_";
    filename += numTostring(sti[3])+"_";
    filename += numTostring(sti[4])+"_";
    filename += "Order_"+numTostring(pmap.order);
    filename += ".dat";
    cout <<  "Trajectory save in " << filename << endl;

    //------------------------------------------
    //Plotting parameters
    //------------------------------------------
    int N = 10000;
    double **yNCE  = dmatrix(0, 5, 0, N);
    double **yNCEc = dmatrix(0, 5, 0, N);


    double *tNCE   = dvector(0, N);
    double *tNCEc  = dvector(0, N);
    double *tNCS   = dvector(0, N);

    double **ySEM   = dmatrix(0, 5, 0, N);
    double **yNCS   = dmatrix(0, 5, 0, N);

    //------------------------------------------
    //Integration
    //------------------------------------------
    trajectory_integration_grid(&orbit, orbit.pmap->t0, orbit.pmap->T, filename, yNCE, tNCE, N, 1);

    //------------------------------------------
    //In SEM coordinates
    //------------------------------------------
    NCEMmtoSEMm_vec(yNCE, tNCE,  ySEM, N, &SEML);

    //-------------------------------------------------
    //Plotting in EM
    //-------------------------------------------------
    // Init
    gnuplot_ctrl  *h1, *h13D;
    h1 = gnuplot_init();
    h13D = gnuplot_init();

    // Settings
    gnuplot_setstyle(h1,   (char*)"lines");
    gnuplot_set_xlabel(h1, (char*)"x [-]");
    gnuplot_set_ylabel(h1, (char*)"y [-]");

    //Plot
    gnuplot_plot_xy(h1, yNCE[0], yNCE[1], N+1, (char*)"NC coordinates", "lines", "1", "2", 1);
    gnuplot_fplot_xyz(yNCE[0], yNCE[1], yNCE[2], N+1, filename.c_str());
    //3D
    if(is3D) gnuplot_plot_xyz(h13D, yNCE[0], yNCE[1], yNCE[2], N+1, (char*)"NC coordinates", "lines", "1", "1", 1);

    //Pause
    //printf("Press ENTER to close the gnuplot window(s)\n");
    //scanf("%c",&ch);

    //-------------------------------------------------
    //Plotting in SEM
    //-------------------------------------------------
    //Init
    gnuplot_ctrl  *h2, *h23D;
    h2 = gnuplot_init();
    h23D = gnuplot_init();

    //Settings
    gnuplot_setstyle(h2, (char*)"lines");
    gnuplot_set_xlabel(h2, (char*)"xSEM [-]");
    gnuplot_set_ylabel(h2, (char*)"ySEM [-]");

    //Plot
    gnuplot_plot_xy(h2, ySEM[0], ySEM[1], N+1, (char*)"SEM coordinates", "lines", "1", "2", 3);
    //3D
    if(is3D) gnuplot_plot_xyz(h23D, ySEM[0], ySEM[1], ySEM[2], N+1, (char*)"SEM coordinates", "lines", "1", "2", 3);

    //printf("Press ENTER to close the gnuplot window(s)\n");
    //scanf("%c",&ch);

    //-------------------------------------------------
    //Earth, Sun and Moon in SEM
    //-------------------------------------------------
    double Pm[3], Pe[3], Ps[3];
    evaluateCoef(Pm, tNCE[0], SEML.us_sem.n, SEML.nf, SEML.cs_sem.Pm, 3);
    evaluateCoef(Pe, tNCE[0], SEML.us_sem.n, SEML.nf, SEML.cs_sem.Pe, 3);
    evaluateCoef(Ps, tNCE[0], SEML.us_sem.n, SEML.nf, SEML.cs_sem.Ps, 3);

    double L1[3], L2[3];
    L1[0] = -SEML.cs_sem_l1.cr3bp.l1.position[0];
    L1[1] = 0.0;
    L1[2] = 0.0;

    L2[0] = -SEML.cs_sem_l2.cr3bp.l2.position[0];
    L2[1] = 0.0;
    L2[2] = 0.0;

    //Plot
    gnuplot_plot_xy(h2, &Pm[0], &Pm[1], 1, (char*)"Moon",  "points", "4", "2", 4);
    gnuplot_plot_xy(h2, &Pe[0], &Pe[1], 1, (char*)"Earth", "points", "5", "2", 5);
    gnuplot_plot_xy(h2, &L1[0], &L1[1], 1, (char*)"L1",   "points", "7", "2", 6);
    gnuplot_plot_xy(h2, &L2[0], &L2[1], 1, (char*)"L2",   "points", "7", "2", 7);

    //3D
    if(is3D) gnuplot_plot_xyz(h23D, &Pm[0], &Pm[1], &Pm[2], 1, (char*)"Moon",  "points", "4", "2", 4);
    if(is3D) gnuplot_plot_xyz(h23D, &Pe[0], &Pe[1], &Pe[2], 1, (char*)"Earth", "points", "5", "2", 5);
    if(is3D) gnuplot_plot_xyz(h23D, &L1[0], &L1[1], &L1[2], 1, (char*)"L1",   "points", "7", "2", 6);
    if(is3D) gnuplot_plot_xyz(h23D, &L2[0], &L2[1], &L2[2], 1, (char*)"L2",   "points", "7", "2", 7);

    //-------------------------------------------------
    //Earth and Moon in NCEM
    //-------------------------------------------------
    double pm[3], pe[3];
    evaluateCoef(pm, tNCE[0], SEML.us_em.n, SEML.nf, SEML.cs_em.pm, 3);
    evaluateCoef(pe, tNCE[0], SEML.us_em.n, SEML.nf, SEML.cs_em.pe, 3);

    //Plot
    gnuplot_plot_xy(h1, &pm[0], &pm[1], 1, (char*)"Moon",  "points", "4", "2", 4);
    gnuplot_plot_xy(h1, &pe[0], &pe[1], 1, (char*)"Earth", "points", "5", "2", 5);
    //3D
    if(is3D) gnuplot_plot_xyz(h13D, &pm[0], &pm[1], &pm[2], 1, (char*)"Moon",  "points", "4", "2", 4);
    if(is3D) gnuplot_plot_xyz(h13D, &pe[0], &pe[1], &pe[2], 1, (char*)"Earth", "points", "5", "2", 5);

    printf("Press ENTER to close the gnuplot window(s)\n");
    scanf("%c",&ch);

    //------------------------------------------------------------------------------
    // Find closest point to SEML2
    //------------------------------------------------------------------------------
    double yvp[6], tvp;
    double minL2Dist, minL2Distc;
    minL2Dist = (ySEM[0][0]-L2[0])*(ySEM[0][0]-L2[0]);
    for(int i = 1; i < 3; i++) minL2Dist += ySEM[i][0]*ySEM[i][0];
    for(int i = 0; i < 6; i++) yvp[i] = ySEM[i][0];
    tvp = tNCE[0]*SEML.us_em.ns;

    //------------------------------------------------------------------------------
    // Second integration
    //------------------------------------------------------------------------------
    double yv[6], tv;
    double tman = 4*orbit.pmap->T;
    for(int p = 7000; p <= 10000; p+=200)
    {
        for(int i = 0; i < 6; i++) yv[i] = yNCE[i][p];
        tv = tNCE[p];

        //----------------------
        // Projection on the center-stable manifold
        //----------------------
        NCprojCCMtoCUS(yv, tv, sti, orbit, +1e-6);

        //----------------------
        // Update the state
        //----------------------
        orbit_update_ic(orbit, sti, tv);
        orbit.int_method =  0;

        //------------------------------------------
        //Integration
        //------------------------------------------
        trajectory_integration_grid(&orbit, tv, tv+tman, filename, yNCEc, tNCEc, N, 0);

        //------------------------------------------
        //In SEM coordinates
        //------------------------------------------
        NCEMmtoSEMm_vec(yNCEc, tNCEc,  ySEM, N, &SEML);

        //------------------------------------------
        //Plotting
        //------------------------------------------
        //EM
        gnuplot_plot_xy(h1, yNCEc[0], yNCEc[1], N+1, (char*)"", "lines", "dashed", "3", 2);
        gnuplot_plot_xy(h1, &yNCEc[0][0], &yNCEc[1][0], 1, (char*)"", "points", "1", "2", 2);
        //3D
        if(is3D) gnuplot_plot_xyz(h13D, yNCEc[0], yNCEc[1], yNCEc[2], N+1, (char*)"", "lines", "dashed", "2", 2);
        if(is3D) gnuplot_plot_xyz(h13D, &yNCEc[0][0], &yNCEc[1][0], &yNCEc[2][0], 1, (char*)"", "points", "1", "2", 2);

        //SEM
        gnuplot_plot_xy(h2, ySEM[0], ySEM[1], N+1, (char*)"", "lines", "dashed", "3", 2);
        gnuplot_plot_xy(h2, &ySEM[0][0], &ySEM[1][0], 1, (char*)"", "points", "1", "2", 2);
        //3D
        if(is3D) gnuplot_plot_xyz(h23D, ySEM[0], ySEM[1], ySEM[2], N+1, (char*)"", "lines", "dashed", "2", 2);
        if(is3D) gnuplot_plot_xyz(h23D, &ySEM[0][0], &ySEM[1][0], &ySEM[2][0], 1, (char*)"", "points", "1", "2", 2);


        //-------------------------------
        // Find closest point to SEML2
        //-------------------------------
        for(int k = 0; k<= N; k++)
        {
            minL2Distc = (ySEM[0][k]-L2[0])*(ySEM[0][k]-L2[0]);
            for(int i = 1; i < 3; i++) minL2Distc += ySEM[i][k]*ySEM[i][k];

            if(minL2Distc < minL2Dist)
            {
                for(int i = 0; i < 6; i++) yvp[i] = ySEM[i][k];
                tvp = tNCEc[k]*SEML.us_em.ns;
                minL2Dist = minL2Distc;
            }
        }


    }

    cout << "tv = " << tvp << endl;
    //Plot closest point
    gnuplot_plot_xy(h2, &yvp[0], &yvp[1], 1, (char*)"Closest to L2", "points", "1", "2", 2);
    printf("Press ENTER to close the gnuplot window(s)\n");
    scanf("%c",&ch);


    //Back to NC coordinates for this point
    double yvc[6];
    SEMtoNC(tvp, yvp, yvc, &SEML);


    //------------------------------------------------------------------------------
    // Integration in SEM
    //------------------------------------------------------------------------------

    //-------------------------------------------------
    //Change FOCUS
    //-------------------------------------------------
    changeCOORDSYS(SEML, F_SEM);
    //------------------------------------------
    // Initialisation of the central manifold
    //------------------------------------------
    updateCM(SEML);
    updateCOC(SEML);

    //----------------------
    // Projection on the center manifold
    //----------------------
    NCprojCCMtoCUS(yvc, tvp, sti, orbit, 0.0);

    //------------------------------------------
    //Update the IC
    //------------------------------------------
    orbit_update_ic(orbit, sti, tvp);
    orbit.int_method =  0;

    //Integration
    trajectory_integration_grid(&orbit, tvp, tvp+2*orbit.pmap->T, filename, yNCS, tNCS, N, 1);

    //NC-SEM to SEM
    NCSEMmtoSEMm_vec(yNCS, tNCS, ySEM, N, &SEML);

    //Plot
    gnuplot_plot_xy(h2, ySEM[0], ySEM[1], N+1, (char*)"SEM coordinates", "lines", "1", "2", 3);
    gnuplot_plot_xy(h2, &ySEM[0][0], &ySEM[1][0],  1, (char*)"SEM coordinates", "points", "1", "2", 3);

    //3D
    //if(is3D) gnuplot_plot_xyz(h23D, ySEM[0], ySEM[1], ySEM[2], N+1, (char*)"SEM coordinates", "lines", "1", "2", 3);


    //Pause
    printf("Press ENTER to close the gnuplot window(s)\n");
    scanf("%c",&ch);

    gnuplot_close(h1);
    gnuplot_close(h2);
    gnuplot_close(h13D);
    gnuplot_close(h23D);

    //-------------------------------------------------
    //Free
    //-------------------------------------------------
    free_dmatrix(yNCE, 0, 2, 0, N);
    free_dvector(tNCE, 0, N);
    free_dmatrix(yNCEc, 0, 2, 0, N);
    free_dvector(tNCEc, 0, N);
}

/**
 *   \brief Integrates a given trajectory with IC in sti[]
 **/
void trajectory(double sti[], Pmap &pmap, int label)
{
    //-------------------------------------------------
    //Integration tools
    //-------------------------------------------------
    //------------------
    //For dot(z) = F(z)
    //------------------
    OdeStruct ode_s_6;
    init_ode_NC(ode_s_6, pmap);
    //For root finding
    OdeStruct ode_s_6_root;
    init_ode_NC(ode_s_6_root, pmap);

    //------------------
    //For dot(s) = fh(s)
    //------------------
    //Reduced vector field structure
    RVF rvf;
    OdeStruct ode_s_8;
    Ofsc rvf_ofs(OFS_ORDER);
    init_ode_CCM(ode_s_8, rvf, rvf_ofs, pmap);
    //Reduced vector field structure for root finding
    RVF rvf_root;
    OdeStruct ode_s_8_root;
    Ofsc rvf_ofs_root(OFS_ORDER);
    init_ode_CCM(ode_s_8_root, rvf_root, rvf_ofs_root, pmap);


    //------------------
    //Event structures
    //------------------
    value_params val_par;
    val_par.dim        = 2;                   //event on z component
    val_par.direction  = 0;                   //all zeros are detected
    val_par.max_events = pmap.max_events;     //maximum of events
    val_par.value      = 0.0;                 //event when z = 0
    value_function fvalue;
    fvalue.val_par     = &val_par;
    fvalue.value       = &linear_intersection;

    //------------------
    //Orbit structure
    //------------------
    Ofsc orbit_ofs(OFS_ORDER);
    Orbit orbit;
    init_orbit(&orbit, &CM, &CMh, &JCM, &Mcoc, &Vcoc, &val_par,
               &fvalue, &ode_s_6, &ode_s_8, &ode_s_6_root, &ode_s_8_root,
               &pmap, &SEML, &orbit_ofs, pmap.vdim, label);


    //------------------------------------------
    //Update the IC
    //------------------------------------------
    orbit_update_ic(orbit, sti, orbit.pmap->t0);
    orbit.int_method =  0;


    //------------------------------------------
    //Integration
    //------------------------------------------
    string filename;
    filename = SEML.cs.F_PLOT+"traj_IC_";
    filename += numTostring(sti[0])+"_";
    filename += numTostring(sti[1])+"_";
    filename += numTostring(sti[2])+"_";
    filename += numTostring(sti[3])+"_";
    filename += "Order_"+numTostring(pmap.order);
    filename += ".dat";

    cout <<  "Trajectory save in " << filename << endl;

    //------------------------------------------
    //Plotting parameters
    //------------------------------------------
    int N = 10000;
    double **yNCE = dmatrix(0, 2, 0, N);
    double *tNCE = dvector(0, N);

    //------------------------------------------
    //Integration
    //------------------------------------------
    trajectory_integration_grid(&orbit, orbit.pmap->t0, 50*orbit.pmap->T, filename, yNCE, tNCE, N, 1);

    //-------------------------------------------------
    //Transformation
    //-------------------------------------------------
    //    double gamma = SEML.cs.gamma;
    //    double c1 = SEML.cs.c1;
    //    double L = SEML.cs.cr3bp.L;

    //    for(int i = 0; i <= N; i++)
    //    {
    //        //From NC to EM coefficients
    //        yNCE[0][i] = -gamma*(yNCE[0][i] - c1);
    //        yNCE[1][i] = -gamma*yNCE[1][i];
    //        yNCE[2][i] = +gamma*yNCE[2][i];
    //
    //        //From EM to Physical coefficients
    //        yNCE[0][i] = L*yNCE[0][i];
    //        yNCE[1][i] = L*yNCE[1][i];
    //        yNCE[2][i] = L*yNCE[2][i];
    //    }

    //-------------------------------------------------
    //Plotting
    //-------------------------------------------------
    char ch;            //Used to close the gnuplot windows at the very end of the program
    gnuplot_ctrl  *h1;
    h1 = gnuplot_init();

    gnuplot_setstyle(h1, (char*)"lines");
    gnuplot_set_xlabel(h1, (char*)"x [-]");
    gnuplot_set_ylabel(h1, (char*)"y [-]");
    gnuplot_plot_xyz(h1, yNCE[0], yNCE[1], yNCE[2], N+1, (char*)"NC coordinates", "lines", "1", "1", 1);
    gnuplot_fplot_xyz(yNCE[0], yNCE[1], yNCE[2], N+1, filename.c_str());

    printf("Press ENTER to close the gnuplot window(s)\n");
    scanf("%c",&ch);
    gnuplot_close(h1);

    //-------------------------------------------------
    //Free
    //-------------------------------------------------
    free_dmatrix(yNCE, 0, 2, 0, N);
}

/**
 *   \brief Integrates a given trajectory up to tf
 **/
int trajectory_integration(Orbit *orbit, double tf)
{
    //------------------------------------------
    //Initialization
    //------------------------------------------
    int status;            //current status
    double yv[6], t, t0;   //current state and time

    //Init the state
    for(int i = 0; i < 6; i++) yv[i] = orbit->z0[i];

    //Init the time;
    t0 = orbit->pmap->t0;
    t = t0;

    //Projection tools
    double ePm;
    int nreset = 1;

    //------------------------------------------
    //Evolving yv(t) up to tf
    //------------------------------------------
    if(tf <0) orbit->ode_s_6->h *= -1;
    status = gslc_proj_evolve(orbit, yv, &t, t0, tf, &ePm, &nreset, 1);


    if (status != GSL_SUCCESS)
    {
        cout << "error in trajectory_integration: integration of yv(t) has gone wrong. break." << endl;
        return GSL_FAILURE;
    }

    return GSL_SUCCESS;

}

/**
 *   \brief Integrates a given trajectory up to tf, on a given grid
 **/
int trajectory_integration_grid(Orbit *orbit, double t0, double tf, string filename, double **yNCE, double *tNCE, int N, int isResetOn)
{
    //------------------------------------------
    //Initialization
    //------------------------------------------
    int status;            //current status
    double yv[6], t;   //current state and time

    //Init the state & time
    for(int i = 0; i < 6; i++) yv[i] = orbit->z0[i];
    t = t0;

    //Projection tools
    double ePm;
    int nreset = 1;

    //Plot
    double ti;

    //------------------------------------------
    //Evolving yv(t) up to tf
    //------------------------------------------
    if(tf <0) orbit->ode_s_6->h *= -1;
    for(int i = 0; i <= N ; i ++)
    {
        ti = t0 + (double) i *tf/N;
        status = gslc_proj_evolve(orbit, yv, &t, t0, ti, &ePm, &nreset, isResetOn);
        for(int k = 0; k < 6; k++) yNCE[k][i] = yv[k];
        tNCE[i] = ti;

        if (status != GSL_SUCCESS)
        {
            cout << "error in trajectory_integration_grid: integration of yv(t) has gone wrong. break." << endl;
            return GSL_FAILURE;
        }
    }


    return GSL_SUCCESS;

}



//--------------------------------------------------------------------------------------------------------------
//
// Steppers
//
//--------------------------------------------------------------------------------------------------------------
/**
 *   \brief Integrates one step the current state yv using projection on the CM if necessary
 **/
int gslc_proj_step(Orbit *orbit,
                   double yv[],
                   double *t,
                   double t0,
                   double t1,
                   double *ePm,
                   int *nreset,
                   int isResetOn)
{
    int status;
    double yvp[6], yvi[6];
    cdouble scp[REDUCED_NV];

    //----------------------
    //Projection tools
    //----------------------
    double omega1 = cimag(Fh[0].getCoef(1,0)->ofs_getCoef(0));
    double omega3 = cimag(Fh[1].getCoef(1,1)->ofs_getCoef(0));

    //----------------------
    //Evolve one step of z(t)
    //----------------------
    status = gsl_odeiv2_evolve_apply (orbit->ode_s_6->e, orbit->ode_s_6->c, orbit->ode_s_6->s, &orbit->ode_s_6->sys, t, t1, &orbit->ode_s_6->h, yv);
    if (status != GSL_SUCCESS)
    {
        cout << "error in gslc_dual_step: integration of z(t) has gone wrong. break." << endl;
        return GSL_FAILURE;
    }

    //----------------------
    //Projection if necessary
    //----------------------
    if(isResetOn && fabs(*t-t0) > fabs(*nreset*orbit->pmap->tt))
    {
        //----------------------
        // Projection on the center manifold
        //----------------------
        //Get the closest point on the center manifold
        NCprojCCM(yv, *t, SEML.us_em.n, OFS_ORDER, MIcoc, Vcoc, omega1, omega3, scp, REDUCED_NV);
        //Update the state
        CCMtoNCbyTFC(scp, *t, orbit->qbcp_l->us.n, orbit->pmap->order,  orbit->pmap->ofs_order,  *orbit->Wh,  *orbit->ofs, Mcoc, Vcoc,  yvp,  orbit->pmap->isGS);
        //For comparison
        for(int i = 0; i <6; i++) yvi[i] = yv[i];
        // Copy of yvp in current state
        for(int i=0; i<6; i++) yv[i]  = yvp[i];

        //-----------------
        // Get the current projection error
        //-----------------
        //Get the current error
        *ePm = fabs(yvi[0] - yv[0]);
        for(int i = 1; i <6 ; i++)
        {
            if(fabs(yvi[i] - yv[i]) > *ePm) *ePm = fabs(yvi[i] - yv[i]);
        }

        if(*ePm > 1e-6)
        {
            cout << "Warning: Reset n° " << *nreset << ". ePm = " << *ePm << endl;
        }
        //-----------------
        //Reset ode structure for next step
        //-----------------
        reset_ode_structure(orbit->ode_s_6);

        //-----------------
        //One additional reset
        //-----------------
        *nreset = *nreset +1;
    }


    return GSL_SUCCESS;
}


/**
 *   \brief Integrates the current state yv up to t = t1, using projection on the CM if necessary
 **/
int gslc_proj_evolve(Orbit *orbit,
                     double yv[],
                     double *t,
                     double t0,
                     double t1,
                     double *ePm,
                     int *nreset,
                     int isResetOn)
{
    reset_ode_structure(orbit->ode_s_6);
    int status;
    do
    {
        status = gslc_proj_step(orbit, yv, t, t0, t1, ePm, nreset, isResetOn);

    }while(fabs(*t)<fabs(t1));

    return status;
}
