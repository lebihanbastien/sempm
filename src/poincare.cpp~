#include "poincare.h"
#include <omp.h>

/**
 * \file poincare.cpp
 * \brief Computation of Poincare, Error, and Period map for the QBCP.
 * \author BLB.
 * \date May 2015
 * \version 1.0
 */

//-----------------------------------------------------------------------------------------------------------------------------------------------------
//
//  Inner routines
//
//-----------------------------------------------------------------------------------------------------------------------------------------------------
/**
 *   \brief Initialize the grid on which the poincare map will be evaluated
 **/
void init_grid(double *grid, double gmin, double gmax, int gsize)
{
    double di, ds;
    for(int i = 0; i <= gsize; i++)
    {
        di = (double) i;
        ds = (double) gsize;
        grid[i] = gmin +  (gmax - gmin)*di/ds;
    }
}

/**
 *   \brief Initialize the ode structure for the NC vector field
 **/
void init_ode_NC(OdeStruct &ode_s_6, Pmap &pmap)
{
    //Root-finding
    const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent;
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    //Init ode structure
    init_ode_structure(&ode_s_6,
                       T,            //stepper: int
                       T_root,       //stepper: root
                       pmap.pabs,    //precision: int_abs
                       pmap.prel,    //precision: int_rel
                       pmap.proot,   //precision: root
                       1e-12,        //precision: diffcorr
                       6,            //dimension
                       1e-6,         //initial int step
                       qbfbp_vfn_novar,
                       NULL,
                       &SEML);
}

/**
 *   \brief Initialize the ode structure for the reduced vector field
 **/
void init_ode_CCM(OdeStruct &ode_s_8, RVF &rvf, Ofsc &rvf_ofs, Pmap &pmap)
{
    //Root-finding
    const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent;
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    //Reduced vector field structure
    rvf.fh        = &Fh;
    rvf.ofs       = &rvf_ofs;
    rvf.order     = pmap.order;
    rvf.ofs_order = pmap.ofs_order;
    rvf.n         = SEML.us.n;
    //Init ode structure
    init_ode_structure(&ode_s_8,
                       T,            //stepper:   int
                       T_root,       //stepper:   root
                       pmap.pabs,    //precision: int_abs
                       pmap.prel,    //precision: int_rel
                       pmap.proot,   //precision: root
                       1e-12,        //precision: diffcorr
                       8,            //dimension
                       1e-6,         //initial int step
                       qbfbp_fh,
                       NULL,
                       &rvf);
}


//-----------------------------------------------------------------------------------------------------------------------------------------------------
//
//  Poincare map
//
//-----------------------------------------------------------------------------------------------------------------------------------------------------
/**
 *   \brief Computes a Poincare map
 *   \param pmap a reference to the Poincare map parameters
 *   \param isPlot if true, the Poincare map is plotted during the computation
 *
 *    Requires initCM and initCOC
 **/
void pmap_build(Pmap &pmap, int append, int method, bool isPlot, bool isPar)
{
    cout << "---------------------------------------------------" << endl;
    cout << "                                                   " << endl;
    cout << "               Poincare map computation            " << endl;
    cout << "                                                   " << endl;
    cout << "---------------------------------------------------" << endl;
    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(15);


    //------------------------------------------
    // Strings
    //------------------------------------------
    string F_GS    = SEML.cs.F_GS;
    string F_PLOT  = SEML.cs.F_PLOT;
    string F_COC   = SEML.cs.F_COC;
    string F_PRINT = SEML.cs.F_PRINT;

    //------------------------------------------
    //Plot variables
    //------------------------------------------
    gnuplot_ctrl  *h1, *h2;
    h1  = gnuplot_init();
    h2  = gnuplot_init();
    char ch;

    //------------------------------------------
    //Energy of Li (st0 = 0)
    //------------------------------------------
    //Null initial conditions in TF coordinates
    double st0[4];
    for(int i =0; i<4; i++) st0[i] = 0.0;
    pmap.H0 = orbit_ham(pmap, st0);
    cout << "pmap.H0 =" << pmap.H0 << endl;

    return;

    //------------------------------------------
    //Set the Hamiltonian to H0+dHv
    //------------------------------------------
    pmap.Hv = pmap.H0+pmap.dHv;

    //------------------------------------------
    //Plot & Print
    //------------------------------------------
    //Filename to print
    string ssHv        = numTostring(pmap.dHv);
    string ssorder     = numTostring(pmap.order);
    string ssofs_order = numTostring(pmap.ofs_order);
    string type;

    //Get the type
    type = pmap.t0 == 0? "Serv_pm_t0_":"Serv_pm_";

    string smethod;
    switch(method)
    {
    case DUAL_INT:
        smethod = "_DUAL_INT";
        break;

    case DUAL_INT_NO_RESET:
        smethod = "_DUAL_INT_NOT_RESET";
        break;

    case DUAL_INT_STEPPED:
        smethod = "_DUAL_INT_STEPPED";
        break;

    case SINGLE_INT:
        smethod = "_SINGLE_INT";
        break;
    }

    string filename;

    //Get the PM style
    switch(SEML.pms)
    {
    case PMS_GRAPH:
        filename = F_PRINT+type; //default case, so no additionnal notations
        break;
    case PMS_NORMFORM:
        filename = F_PRINT+type+"NF_";//Normal form style
        break;
    case PMS_MIXED:
        filename = F_PRINT+type+"MX_";  //Normal form style
        break;
    }
    //Final name
    filename = filename+"Energy_"+ssHv+"_order_"+ssorder+"_ofs_"+ssofs_order+smethod;

    //------------------------------------------
    //If no appending, the right header is written in the txt file
    //------------------------------------------
    //if(!append) header_pmap_fprint(filename);
    if(!append) header_pmap_fprint_small(filename);
    cout << "pmap_build. data saved in " << filename << endl;

    //------------------------------------------
    //Building the working grid
    //------------------------------------------
    double *grid = dvector(0,  pmap.gsize);
    init_grid(grid, pmap.gmin, pmap.gmax, pmap.gsize);
    double numberOfOrbits = pow(1.0+pmap.gsize, 2.0);

    //------------------------------------------
    //Loop
    //------------------------------------------
    int label = 1;
    cout << std::noshowpos << resetiosflags(ios::scientific)  << setprecision(4);
    #pragma omp parallel for if(isPar) shared(label)
    for(int i = 0; i <= pmap.gsize; i++)
    {
        #pragma omp parallel for if(isPar) shared(label)
        for(int j = 0; j <= pmap.gsize; j++)
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
            double time;
            int status;
            Ofsc orbit_ofs(OFS_ORDER);
            Orbit orbit;
            init_orbit(&orbit, &CM, &CMh, &JCM, &Mcoc, &Vcoc, &val_par,
                       &fvalue, &ode_s_6, &ode_s_8, &ode_s_6_root, &ode_s_8_root,
                       &pmap, &SEML, &orbit_ofs, pmap.vdim, label);

            //------------------------------------------
            // Move on the grid
            //------------------------------------------
            double *sti = dvector(0,3);
            sti[0] = grid[i];
            sti[1] = 0.0;
            sti[2] = grid[j];
            sti[3] = 0.0;

            //------------------------------------------
            //Update the IC so that the energy matches the one the Pmap
            //------------------------------------------
            status = orbit_init_pmap(orbit, sti);

            //Between the lines: if we do not want to correct the energy
            //-------------------------//
            //status = GSL_SUCCESS;
            //orbit_update_ic(orbit, sti);
            //cout << "Energy = " << orbit_ham(pmap, sti) << endl;
            //-------------------------//

            //------------------------------------------
            //If the correction is succesful, we can go on with the computation of the return map
            //------------------------------------------
            if(status == GSL_SUCCESS)
            {
                //------------------------------------------
                // Computation
                //------------------------------------------
                tic();
                status = orbit_compute_pmap(&orbit, method);
                time = toc();

                //------------------------------------------
                // Postprocess
                //------------------------------------------
                if(status == GSL_SUCCESS)
                {
                    cout << "Return map " <<  label << "/" << numberOfOrbits << " computed in " << time << "s. " << endl;
                    //------------------------------------------
                    // Print in file
                    //------------------------------------------
                    #pragma omp critical
                    {
                        orbit.label = label++;
                        orbit.tf    = pmap.tf;
                        //orbit_pmap_fprint(&orbit, filename, 1);
                        orbit_pmap_fprint_small(&orbit, filename, 1);
                    }


                }
            }

            //Memory release
            free_orbit(&orbit);
            free_dvector(sti, 0, 3);

        }
    }

    //Memory release
    free_dvector(grid, 0,  pmap.gsize);

    //Plot handling at the end
    if(isPlot)
    {
        printf("Press ENTER to close the gnuplot window(s)\n");
        scanf("%c",&ch);

        //Save in EPS format
        gnuplot_cmd(h1, "set terminal postscript eps solid color enhanced");
        gnuplot_cmd(h1, "set output \"pmap.eps\"");
        gnuplot_cmd(h1, "replot");
    }
    gnuplot_close(h1);
    gnuplot_close(h2);
}

/**
 *   \brief Precision on a Poincare map. Parallelized version
 *   \param pmap a reference to the Poincare maps parameters
 *
 *    Requires initCM and initCOC
 *
 *   REMARK: may be good to "force" p36(t0= 0.0) (e.g. pmap.t0 =  +1.044814582930593 for L2) so that each IC begins on z = 0 plane
 *   If so, the way H(0) = cst is guaranteed must be changed because we need also to ensure that s4 = 0.0 (which is not the case, since it is the variable
 *   that ensures H(0) = cst
 **/
void pmap_precision(Pmap &pmap, int append, bool isPar)
{
    cout << "---------------------------------------------------" << endl;
    cout << "                                                   " << endl;
    cout << "               Poincare map computation            " << endl;
    cout << "                                                   " << endl;
    cout << "---------------------------------------------------" << endl;
    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(15);

    //------------------------------------------
    // Strings
    //------------------------------------------
    string F_GS    = SEML.cs.F_GS;
    string F_PLOT  = SEML.cs.F_PLOT;
    string F_COC   = SEML.cs.F_COC;
    string F_PRINT = SEML.cs.F_PRINT;

    cout << "F_GS    = " << F_GS << endl;
    cout << "F_PLOT  = " << F_PLOT << endl;
    cout << "F_COC   = " << F_COC << endl;
    cout << "F_PRINT = " << F_PRINT << endl;

    //------------------------------------------
    //Plot & Print
    //------------------------------------------
    //Filename to print
    string ssOrder, ssOfs;
    ssOrder = numTostring(pmap.order);
    ssOfs   = numTostring(pmap.ofs_order);
    //------------------------------------------
    //Energy of Li (sti = 0)
    //------------------------------------------
    double st0[4];   //initial conditions in real TFC coordinates
    for(int i =0; i<4; i++) st0[i] = 0.0;
    pmap.H0 = orbit_ham(pmap, st0);

    //------------------------------------------
    //Set the Hamiltonian to H0+dHv
    //------------------------------------------
    pmap.Hv = pmap.H0+pmap.dHv;

    //------------------------------------------
    //Building the working grid
    //------------------------------------------
    double *grid_s1 = dvector(0,  pmap.gsize);
    double *grid_s2 = dvector(0,  pmap.gsize);

    init_grid(grid_s1, pmap.gmin,  pmap.gmax, pmap.gsize);
    init_grid(grid_s2, pmap.gmin,  pmap.gmax, pmap.gsize);
    double numberOfOrbits = pow(1.0+pmap.gsize, 1.0);

    //------------------------------------------
    //Header in txt file if no appending (reset the file!)
    //
    // €€ TODO: better way to init the filename! earlier in the loop
    // €€ TODO: we need to take into account the various "standard" way of evaluating the precision (planar0, z1, 3d...)
    //
    //------------------------------------------
    string filename;

    switch(SEML.pms)
    {
    case PMS_GRAPH:
        filename = F_PRINT+"eOm_planar0_ofs_"+ssOfs+"_order_"+ssOrder; //default case, so no additionnal notations
        break;
    case PMS_NORMFORM:
        filename = F_PRINT+"eOm_NF_planar0_ofs_"+ssOfs+"_order_"+ssOrder; //Normal form style
        break;
    case PMS_MIXED:
        filename = F_PRINT+"eOm_MX_planar0_ofs_"+ssOfs+"_order_"+ssOrder; //Normal form style
        break;
    }
    if(!append) header_precision_fprint(filename);
    cout << "emap_build. data saved in " << filename << endl;
    //------------------------------------------
    //Loop
    //------------------------------------------
    int label = 1;
    #pragma omp parallel for if(isPar) shared(label)
    for(int i = 0; i <= pmap.gsize; i++)
    {
        /*#pragma omp parallel for if(isPar)
        for(int j = 0; j <= pmap.gsize; j++)
        {
            #pragma omp parallel  for if(isPar)
            for(int k = 0; k <= pmap.gsize; k++)
            {
                #pragma omp parallel  for if(isPar)
                for(int l = 0; l <= pmap.gsize; l++)
                {*/
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
        int status;
        Ofsc orbit_ofs(OFS_ORDER);
        Orbit orbit;
        init_orbit(&orbit, &CM, &CMh, &JCM, &Mcoc, &Vcoc, &val_par,
                   &fvalue, &ode_s_6, &ode_s_8, &ode_s_6_root, &ode_s_8_root,
                   &pmap, &SEML, &orbit_ofs, pmap.vdim, label);
        //------------------------------------------
        // Move on the grid
        //------------------------------------------
        double *sti = dvector(0, 3);
        sti[0] = grid_s1[i];
        sti[1] = 0.0;//grid_s1[i];
        sti[2] = 0.0;//grid_s1[j];
        sti[3] = 0.0;

        //------------------------------------------
        //Update the IC so that the energy matches the one of the Pmap
        //------------------------------------------
        //status = orbit_init_pmap(orbit, sti);

        //Between the lines: if we do not want to correct the energy
        //-------------------------//
        status = GSL_SUCCESS;
        orbit_update_ic(orbit, sti, orbit.pmap->t0);
        orbit.int_method =  0;
        //-------------------------//

        //If the correction is successful, we can go on with the computation of the return map
        if(status == GSL_SUCCESS)
        {
            //------------------------------------------
            // Computation
            //------------------------------------------
            dual_emap(&orbit);
            #pragma omp critical
            {
                cout << setprecision(3) <<  "Orbital error after " << (pmap.tt/SEML.us.T) << "T = " << orbit.eOm << endl;
                orbit.label = ++label;
                if(!isPar) cout << "Label "  << label << "/" << (int) numberOfOrbits << endl;  //display label only if no parallel computation
                cout << setprecision(15);

                //------------------------------------------
                // Print in file
                //------------------------------------------
                orbit_precision_fprint(&orbit, filename, 1);
            }
        }
        else cout << ""; //cout << "pmap_build. Energy mismatch prevent return map to be computed." << endl;

        //memory release
        free_orbit(&orbit);
        free_dvector(sti, 0, 3);
        /*}
        }*/
        /*}*/
    }

    //Memory release
    free_dvector(grid_s1, 0,  pmap.gsize);
    free_dvector(grid_s2, 0,  pmap.gsize);

}

/**
 *   \brief Invariance error computation
 *   \param pmap a reference to the Poincare maps parameters
 *   \param hzmax the maximum energy value allowed.
 *          Note that only positive dhz are selected
 *
 *    Requires initCM and initCOC
 *
 *   REMARK: may be good to "force" p36(t0= 0.0) (e.g. pmap.t0 =  +1.044814582930593 for L2) so that each IC begins on z = 0 plane
 *   If so, the way H(0) = cst is guaranteed must be changed because we need also to ensure that s4 = 0.0 (which is not the case, since it is the variable
 *   that ensures H(0) = cst
 **/
void pmap_invariance_error(Pmap &pmap, int append, bool isPar, double hzmax)
{
    cout << "---------------------------------------------------" << endl;
    cout << "                                                   " << endl;
    cout << "          Invariance error map computation         " << endl;
    cout << "                                                   " << endl;
    cout << "---------------------------------------------------" << endl;
    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(15);

    //------------------------------------------
    // Strings
    //------------------------------------------
    string F_GS    = SEML.cs.F_GS;
    string F_PLOT  = SEML.cs.F_PLOT;
    string F_COC   = SEML.cs.F_COC;
    string F_PRINT = SEML.cs.F_PRINT;


    //------------------------------------------
    //Plot & Print
    //------------------------------------------
    //Filename to print
    string ssOrder, ssOfs, ssEnergy;
    ssOrder = numTostring(pmap.order);
    ssOfs   = numTostring(pmap.ofs_order);
    ssEnergy   = numTostring(hzmax);

    //------------------------------------------
    //Energy of Li (sti = 0)
    //------------------------------------------
    double st0[4];   //initial conditions in real TFC coordinates
    for(int i =0; i<4; i++) st0[i] = 0.0;
    pmap.H0 = orbit_ham(pmap, st0);

    //------------------------------------------
    //Set the Hamiltonian to H0+dHv
    //------------------------------------------
    pmap.Hv = pmap.H0+pmap.dHv;
    cout << "pmap.Hv = " << pmap.Hv << endl;

    //------------------------------------------
    //Building the working grid
    //------------------------------------------
    double *grid_s1 = dvector(0,  pmap.gsize);
    double *grid_s2 = dvector(0,  pmap.gsize);

    init_grid(grid_s1, pmap.gmin,  pmap.gmax, pmap.gsize);
    init_grid(grid_s2, pmap.gmin,  pmap.gmax, pmap.gsize);
    double numberOfOrbits = pow(1.0+pmap.gsize, 2.0);

    //------------------------------------------
    //Header in txt file if no appending (reset the file!)
    //------------------------------------------
    string filename;
    string type = "s1s3";

    switch(SEML.pms)
    {
    case PMS_GRAPH:
        filename = F_PRINT+"eIm_"+type+"_ofs_"+ssOfs+"_order_"+ssOrder+"_hmax_"+ssEnergy; //default case, so no additionnal notations
        break;
    case PMS_NORMFORM:
        filename = F_PRINT+"eIm_NF_"+type+"_ofs_"+ssOfs+"_order_"+ssOrder+"_hmax_"+ssEnergy; //Normal form style
        break;
    case PMS_MIXED:
        filename = F_PRINT+"eIm_MX_"+type+"_ofs_"+ssOfs+"_order_"+ssOrder+"_hmax_"+ssEnergy; //Normal form style
        break;
    }
    if(!append) header_precision_fprint(filename);


    double percent;
    int label = 1;
    //------------------------------------------
    //Loop
    //------------------------------------------
    #pragma omp parallel for if(isPar)
    for(int i = 0; i <= pmap.gsize; i++)
    {
        #pragma omp parallel  for if(isPar)
        for(int j = 0; j <= pmap.gsize; j++)
        {
            /*#pragma omp parallel  for if(isPar)
            for(int k = 0; k <= pmap.gsize; k++)
            {*/
            /*#pragma omp parallel  for if(isPar)
            for(int l = 0; l <= pmap.gsize; l++)
            {*/
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
            int status;
            Ofsc orbit_ofs(OFS_ORDER);
            Orbit orbit;
            init_orbit(&orbit, &CM, &CMh, &JCM, &Mcoc, &Vcoc, &val_par,
                       &fvalue, &ode_s_6, &ode_s_8, &ode_s_6_root, &ode_s_8_root,
                       &pmap, &SEML, &orbit_ofs, pmap.vdim, label);

            //------------------------------------------
            // Move on the grid
            //------------------------------------------
            double *sti = dvector(0, 3);
            sti[0] = grid_s1[i];
            sti[1] = 0.0;
            sti[2] = grid_s1[j];
            sti[3] = 0.0;

            //Between the lines: if we do not want to correct the energy
            //-------------------------//
            status = GSL_SUCCESS;
            orbit_update_ic(orbit, sti, orbit.pmap->t0);
            orbit.int_method =  0;
            //-------------------------//

            //If the correction is successful, we can go on with the computation of the return map
            if(status == GSL_SUCCESS)
            {
                double eI[6], fnc[6], eIm;
                cdouble eId;
                //------------------------------------------
                // Computation
                //------------------------------------------
                //Evluating the vector field at z0 = W(s0, t0)
                qbfbp_vfn_novar(orbit.pmap->t0, orbit.z0, fnc, orbit.qbcp_l);

                //-------------------------------
                // eI:
                // - F(W(s,t)) - DWf(s,t) - Wdot(s,t)
                //------------------------------
                for(int p = 1; p < NV; p++)
                {
                    //The vector field was computed in fnc
                    eId = fnc[p]+0.0*I;
                    //Computing DW*Fh(s0, t0)
                    DWFh[p].evaluate(orbit.s0, *orbit.ofs, orbit.pmap->order, orbit.pmap->ofs_order);
                    eId -= orbit.ofs->evaluate(orbit.n*orbit.pmap->t0, orbit.pmap->ofs_order);

                    //Computing Wdot(s0, t0)
                    CMdot[p].evaluate(orbit.s0, *orbit.ofs, orbit.pmap->order, orbit.pmap->ofs_order);
                    eId -= orbit.ofs->evaluate(orbit.n*orbit.pmap->t0, orbit.pmap->ofs_order);
                    //Taking the error as the absolute value (norm)
                    eI[p] = cabs(eId);
                }

                //-------------------------------
                //Taking the maximum (infinity norm)
                //-------------------------------
                eIm = eI[0];
                for(int p = 1; p < NV; p++) if(eI[p] > eIm) eIm = eI[p];

                //-------------------------------
                // Set the error in orbit.eOm
                //-------------------------------
                orbit.eOm = eIm;

                #pragma omp critical
                {
                    //-------------------------------
                    // Displays the current state of the computation
                    //-------------------------------
                    cout << std::noshowpos << resetiosflags(ios::scientific)  << setprecision(4);
                    if(label%1000 == 0)
                    {
                        percent = (double) label/numberOfOrbits*100;
                        std::cout << "\r" << percent << "% completed: ";
                        std::cout << std::string(floor(0.1*percent), '|') << endl;
                        std::cout.flush();
                    }
                    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(15);


                    //------------------------------------------
                    // Print in file
                    //------------------------------------------
                    orbit.label = ++label;
                    orbit_energy_fprint(&orbit, filename, hzmax, 1);
                }
            }
            else cout << "";

            //memory release
            free_orbit(&orbit);
            free_dvector(sti, 0, 3);
        }
        /*}*/
    }


    //Memory release
    free_dvector(grid_s1, 0,  pmap.gsize);
    free_dvector(grid_s2, 0,  pmap.gsize);
}

/**
 *   \brief Test error computation
 *   \param pmap a reference to the Poincare maps parameters
 *   \param hzmax the maximum energy value allowed.
 *          Note that only positive dhz are selected
 *
 *    Requires initCM and initCOC
 *
 *   REMARK: may be good to "force" p36(t0= 0.0) (e.g. pmap.t0 =  +1.044814582930593 for L2) so that each IC begins on z = 0 plane
 *   If so, the way H(0) = cst is guaranteed must be changed because we need also to ensure that s4 = 0.0 (which is not the case, since it is the variable
 *   that ensures H(0) = cst
 **/
void pmap_test_error(Pmap &pmap, int append, bool isPar, double hzmax)
{
    cout << "---------------------------------------------------" << endl;
    cout << "                                                   " << endl;
    cout << "          Invariance error map computation         " << endl;
    cout << "                                                   " << endl;
    cout << "---------------------------------------------------" << endl;
    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(15);

    //------------------------------------------
    // Strings
    //------------------------------------------
    string F_GS    = SEML.cs.F_GS;
    string F_PLOT  = SEML.cs.F_PLOT;
    string F_COC   = SEML.cs.F_COC;
    string F_PRINT = SEML.cs.F_PRINT;


    //------------------------------------------
    //Plot & Print
    //------------------------------------------
    //Filename to print
    string ssOrder, ssOfs, ssEnergy;
    ssOrder = numTostring(pmap.order);
    ssOfs   = numTostring(pmap.ofs_order);
    ssEnergy   = numTostring(hzmax);

    //------------------------------------------
    //Energy of Li (sti = 0)
    //------------------------------------------
    double st0[4];   //initial conditions in real TFC coordinates
    for(int i =0; i<4; i++) st0[i] = 0.0;
    pmap.H0 = orbit_ham(pmap, st0);

    //------------------------------------------
    //Set the Hamiltonian to H0+dHv
    //------------------------------------------
    pmap.Hv = pmap.H0+pmap.dHv;
    cout << "pmap.Hv = " << pmap.Hv << endl;

    //------------------------------------------
    //Building the working grid
    //------------------------------------------
    double *grid_s1 = dvector(0,  pmap.gsize);
    double *grid_s2 = dvector(0,  pmap.gsize);

    init_grid(grid_s1, pmap.gmin,  pmap.gmax, pmap.gsize);
    init_grid(grid_s2, pmap.gmin,  pmap.gmax, pmap.gsize);
    double numberOfOrbits = pow(1.0+pmap.gsize, 2.0);

    //------------------------------------------
    //Header in txt file if no appending (reset the file!)
    //------------------------------------------
    string filename;
    string type = "s1s3";

    switch(SEML.pms)
    {
    case PMS_GRAPH:
        filename = F_PRINT+"eIm_"+type+"_ofs_"+ssOfs+"_order_"+ssOrder+"_hmax_"+ssEnergy; //default case, so no additionnal notations
        break;
    case PMS_NORMFORM:
        filename = F_PRINT+"eIm_NF_"+type+"_ofs_"+ssOfs+"_order_"+ssOrder+"_hmax_"+ssEnergy; //Normal form style
        break;
    case PMS_MIXED:
        filename = F_PRINT+"eIm_MX_"+type+"_ofs_"+ssOfs+"_order_"+ssOrder+"_hmax_"+ssEnergy; //Normal form style
        break;
    }
    if(!append) header_precision_fprint(filename);


    double percent;
    int label = 1;

    //------------------------------------------
    // Init cos/sin arrays
    //------------------------------------------
    double cR[pmap.ofs_order];
    double sR[pmap.ofs_order];
    initcRsR(SEML.us_em.n*pmap.t0, cR, sR, pmap.ofs_order);

    //------------------------------------------
    //Loop
    //------------------------------------------
    tic();
    #pragma omp parallel for if(isPar)
    for(int i = 0; i <= pmap.gsize; i++)
    {
        #pragma omp parallel  for if(isPar)
        for(int j = 0; j <= pmap.gsize; j++)
        {
            /*#pragma omp parallel  for if(isPar)
            for(int k = 0; k <= pmap.gsize; k++)
            {*/
            /*#pragma omp parallel  for if(isPar)
            for(int l = 0; l <= pmap.gsize; l++)
            {*/
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
            int status;
            Ofsc orbit_ofs(OFS_ORDER);
            Orbit orbit;
            init_orbit(&orbit, &CM, &CMh, &JCM, &Mcoc, &Vcoc, &val_par,
                       &fvalue, &ode_s_6, &ode_s_8, &ode_s_6_root, &ode_s_8_root,
                       &pmap, &SEML, &orbit_ofs, pmap.vdim, label);

            //------------------------------------------
            // Move on the grid
            //------------------------------------------
            double *sti = dvector(0, 3);
            sti[0] = grid_s1[i];
            sti[1] = 0.0;
            sti[2] = grid_s1[j];
            sti[3] = 0.0;

            //Between the lines: if we do not want to correct the energy
            //-------------------------//
            status = GSL_SUCCESS;
            orbit_update_ic(orbit, sti, orbit.pmap->t0);
            orbit.int_method =  0;
            //-------------------------//

            //If the correction is successful, we can go on with the computation of the return map
            if(status == GSL_SUCCESS)
            {
                double eI[6], fnc[6], eIm;
                cdouble eId;
                //------------------------------------------
                // Computation
                //------------------------------------------
                //Evluating the vector field at z0 = W(s0, t0)
                qbfbp_vfn_novar(orbit.pmap->t0, orbit.z0, fnc, orbit.qbcp_l);

                //-------------------------------
                // eI:
                // - F(W(s,t)) - DWf(s,t) - Wdot(s,t)
                //------------------------------
                for(int p = 1; p < NV; p++)
                {
                    //Computing DW*Fh(s0, t0)
                    DWFh[p].evaluate(orbit.s0, *orbit.ofs, orbit.pmap->order, orbit.pmap->ofs_order);
                    eId = orbit.ofs->evaluate(orbit.n*orbit.pmap->t0, orbit.pmap->ofs_order);
                    //Computing it in another manner
                    //eId = DWFh[p].fevaluate(orbit.s0, orbit.n*orbit.pmap->t0, orbit.pmap->order, orbit.pmap->ofs_order);
                    eId = DWFh[p].fevaluate(orbit.s0, cR, sR, orbit.pmap->order, orbit.pmap->ofs_order);
                    //Taking the error as the absolute value (norm)
                    eI[p] = cabs(eId);
                }

                //-------------------------------
                //Taking the maximum (infinity norm)
                //-------------------------------
                eIm = eI[0];
                for(int p = 1; p < NV; p++) if(eI[p] > eIm) eIm = eI[p];

                //-------------------------------
                // Set the error in orbit.eOm
                //-------------------------------
                orbit.eOm = eIm;
                cout << std::noshowpos << resetiosflags(ios::scientific)  << setprecision(4);
                if(label%1000 == 0)
                {
                    percent = (double) label/numberOfOrbits*100;
                    std::cout << "\r" << percent << "% completed: ";
                    std::cout << std::string(floor(0.1*percent), '|') << endl;
                    std::cout.flush();
                }
                cout << std::showpos << setiosflags(ios::scientific)  << setprecision(15);


                //------------------------------------------
                // Print in file
                //------------------------------------------
                orbit.label = ++label;
                orbit_energy_fprint(&orbit, filename, hzmax, 1);
            }
            else cout << "";

            //memory release
            free_orbit(&orbit);
            free_dvector(sti, 0, 3);
        }
        /*}*/
    }
    cout << "imap computed in " << toc() << " s. " << endl;

    //Memory release
    free_dvector(grid_s1, 0,  pmap.gsize);
    free_dvector(grid_s2, 0,  pmap.gsize);
}

/**
 *   \brief Energy of the initial conditions on a Poincare map. Parallelized version
 *   \param pmap a reference to the Poincare maps parameters
 *   \param hzmax the maximum energy value allowed.
 *          Note that only positive dhz are selected
 *
 *    Requires initCM and initCOC
 *
 *   REMARK: may be good to "force" p36(t0= 0.0) (e.g. pmap.t0 =  +1.044814582930593 for L2) so that each IC begins on z = 0 plane
 *   If so, the way H(0) = cst is guaranteed must be changed because we need also to ensure that s4 = 0.0 (which is not the case, since it is the variable
 *   that ensures H(0) = cst
 **/
void pmap_energy(Pmap &pmap, int append, bool isPar, double hzmax)
{
    cout << "---------------------------------------------------" << endl;
    cout << "                                                   " << endl;
    cout << "               Poincare map computation            " << endl;
    cout << "                                                   " << endl;
    cout << "---------------------------------------------------" << endl;
    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(15);

    //------------------------------------------
    // Strings
    //------------------------------------------
    string F_GS    = SEML.cs.F_GS;
    string F_PLOT  = SEML.cs.F_PLOT;
    string F_COC   = SEML.cs.F_COC;
    string F_PRINT = SEML.cs.F_PRINT;

    //------------------------------------------
    //Energy of Li (st0 = 0)
    //------------------------------------------
    //Null initial conditions in TF coordinates
    double st0[4];
    for(int i =0; i<4; i++) st0[i] = 0.0;
    pmap.H0 = orbit_ham(pmap, st0);

    //------------------------------------------
    //Set the Hamiltonian to H0+dHv
    //------------------------------------------
    pmap.Hv = pmap.H0+pmap.dHv;
    cout << "pmap.Hv = " << pmap.Hv << endl;

    //------------------------------------------
    //Plot & Print
    //------------------------------------------
    //Filename to print
    string ssHv        = numTostring(pmap.dHv);
    string ssorder     = numTostring(pmap.order);
    string ssofs_order = numTostring(pmap.ofs_order);
    string type;

    //Get the type
    type = "Serv_hm_";
    string filename;

    //Get the PM style
    switch(SEML.pms)
    {
    case PMS_GRAPH:
        filename = F_PRINT+type; //default case, so no additionnal notations
        break;
    case PMS_NORMFORM:
        filename = F_PRINT+type+"NF_";//Normal form style
        break;
    case PMS_MIXED:
        filename = F_PRINT+type+"MX_";  //Normal form style
        break;
    }
    //Final name
    filename = filename+"Energy_"+ssHv+"_order_"+ssorder+"_ofs_"+ssofs_order;

    //------------------------------------------
    //If no appending, the right header is written in the txt file
    //------------------------------------------
    if(!append) header_energy_fprint(filename);
    cout << "pmap_build. data saved in " << filename << endl;

    //------------------------------------------
    //Building the working grid
    //------------------------------------------
    double *grid = dvector(0,  pmap.gsize);
    init_grid(grid, pmap.gmin, pmap.gmax, pmap.gsize);
    double numberOfOrbits = pow(1.0+pmap.gsize, 3.0);

    //------------------------------------------
    //Loop
    //------------------------------------------
    int label = 1;
    double percent;
    #pragma omp parallel for if(isPar) shared(label)
    for(int i = 0; i <= pmap.gsize; i++)
    {
        #pragma omp parallel  for if(isPar) shared(label)
        for(int j = 0; j <= pmap.gsize; j++)
        {
            #pragma omp parallel  for if(isPar) shared(label)
            for(int k = 0; k <= pmap.gsize; k++)
            {
                /*#pragma omp parallel  for if(isPar)
                for(int l = 0; l <= pmap.gsize; l++)
                {*/
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
                int status;
                Ofsc orbit_ofs(OFS_ORDER);
                Orbit orbit;
                init_orbit(&orbit, &CM, &CMh, &JCM, &Mcoc, &Vcoc, &val_par,
                           &fvalue, &ode_s_6, &ode_s_8, &ode_s_6_root, &ode_s_8_root,
                           &pmap, &SEML, &orbit_ofs, pmap.vdim, label);

                //------------------------------------------
                // Move on the grid
                //------------------------------------------
                double *sti = dvector(0, 3);
                sti[0] = grid[i];
                sti[1] = grid[j]/4;  //for maxvalue = 60, s2 = 15
                sti[2] = grid[k];
                sti[3] = 0.0;

                //Between the lines: if we do not want to correct the energy
                //-------------------------//
                status = GSL_SUCCESS;
                orbit_update_ic(orbit, sti, orbit.pmap->t0);
                orbit.int_method =  0;
                //-------------------------//

                //If the correction is successful, we can go on with the computation of the return map
                if(status == GSL_SUCCESS)
                {
                    #pragma omp critical
                    {
                        //------------------------------------------
                        // Computation
                        //------------------------------------------
                        cout << std::noshowpos << resetiosflags(ios::scientific)  << setprecision(4);
                        if(label%1000 == 0)
                        {
                            percent = (double) label/numberOfOrbits*100;
                            std::cout << "\r" << percent << "% completed: ";
                            std::cout << std::string(floor(0.1*percent), '|') << endl;
                            std::cout.flush();
                        }
                        cout << std::showpos << setiosflags(ios::scientific)  << setprecision(15);


                        //------------------------------------------
                        // Print in file
                        //------------------------------------------
                        orbit.label = ++label;
                        orbit_energy_fprint(&orbit, filename, hzmax, 1);
                    }
                }
                else cout << "";

                //memory release
                free_orbit(&orbit);
                free_dvector(sti, 0, 3);
            }
        }
    }


    //Memory release
    free_dvector(grid, 0,  pmap.gsize);
}

/**
 *   \brief Computes a Stroboscopic map
 *   \param pmap a reference to the Stroboscopic map parameters
 *   \param isPlot if true, the Stroboscopic map is plotted during the computation
 *
 *    Requires initCM and initCOC
 **/
void tmap_build(Pmap &pmap, int append, int method, bool isPlot, bool isPar)
{
    cout << "---------------------------------------------------" << endl;
    cout << "                                                   " << endl;
    cout << "               Stroboscopic map computation        " << endl;
    cout << "                                                   " << endl;
    cout << "---------------------------------------------------" << endl;
    //cout << std::showpos << setiosflags(ios::scientific)  << setprecision(15);


    //------------------------------------------
    // Strings
    //------------------------------------------
    string F_GS    = SEML.cs.F_GS;
    string F_PLOT  = SEML.cs.F_PLOT;
    string F_COC   = SEML.cs.F_COC;
    string F_PRINT = SEML.cs.F_PRINT;

    //------------------------------------------
    //Plot variables
    //------------------------------------------
    gnuplot_ctrl  *h1, *h2;
    h1  = gnuplot_init();
    h2  = gnuplot_init();
    char ch;

    //------------------------------------------
    //Energy of Li (st0 = 0)
    //------------------------------------------
    //Null initial conditions in TF coordinates
    double st0[4];
    for(int i =0; i<4; i++) st0[i] = 0.0;
    pmap.H0 = orbit_ham(pmap, st0);

    //------------------------------------------
    //Set the Hamiltonian to H0+dHv
    //------------------------------------------
    pmap.Hv = pmap.H0+pmap.dHv;

    //------------------------------------------
    //Plot & Print
    //------------------------------------------
    //Filename to print
    string ssHv = numTostring(pmap.dHv);
    string ssorder     = numTostring(pmap.order);
    string ssofs_order = numTostring(pmap.ofs_order);

    //Get the type
    string type = "Serv_tm_";

    //Get the method
    string smethod;
    switch(method)
    {
    case DUAL_INT:
        smethod = "_DUAL_INT";
        break;

    case DUAL_INT_NO_RESET:
        smethod = "_DUAL_INT_NOT_RESET";
        break;

    case DUAL_INT_STEPPED:
        smethod = "_DUAL_INT_STEPPED";
        break;

    case SINGLE_INT:
        smethod = "_SINGLE_INT";
        break;
    }

    //Get the PM style
    string filename;
    switch(SEML.pms)
    {
    case PMS_GRAPH:
        filename = F_PRINT+type; //default case, so no additionnal notations
        break;
    case PMS_NORMFORM:
        filename = F_PRINT+type+"NF_";//Normal form style
        break;
    case PMS_MIXED:
        filename = F_PRINT+type+"MX_";  //Normal form style
        break;
    }
    //Final name
    filename = filename+"Energy_"+ssHv+"_order_"+ssorder+"_ofs_"+ssofs_order+smethod;

    //------------------------------------------
    //If no appending, the right header is written in the txt file
    //------------------------------------------
    //if(!append) header_pmap_fprint(filename);
    if(!append) header_pmap_fprint_small(filename);
    cout << "pmap_build. data saved in " << filename << endl;

    //------------------------------------------
    //Building the working grid
    //------------------------------------------
    double *gridx = dvector(0,  pmap.gsize);
    double *gridy = dvector(0,  pmap.gsize);
    init_grid(gridx, pmap.gmin, pmap.gmax, pmap.gsize);
    init_grid(gridy, pmap.gmin, pmap.gmax, pmap.gsize);

    double numberOfOrbits = pow(1.0+pmap.gsize, 2.0);

    //------------------------------------------
    //Loop
    //------------------------------------------
    int label = 1;
    #pragma omp parallel for if(isPar)
    for(int i = 0; i <= pmap.gsize; i++)
    {
        #pragma omp parallel for if(isPar)
        for(int j = 0; j <= pmap.gsize; j++)
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
            int status;
            double time;
            Ofsc orbit_ofs(OFS_ORDER);
            Orbit orbit;
            init_orbit(&orbit, &CM, &CMh, &JCM, &Mcoc, &Vcoc, &val_par,
                       &fvalue, &ode_s_6, &ode_s_8, &ode_s_6_root, &ode_s_8_root,
                       &pmap, &SEML, &orbit_ofs, pmap.vdim, label);


            //------------------------------------------
            // Move on the grid
            //------------------------------------------
            double *sti = dvector(0,3);
            sti[0] = gridx[i];
            sti[1] = 0.0;
            sti[2] = gridy[j];
            sti[3] = 0.0;

            //------------------------------------------
            //Update the IC so that the energy matches the one the Pmap
            //------------------------------------------
            //status = orbit_init_pmap(orbit, sti);

            //------------------------------------------
            //Between the lines: if we do not want to correct the energy
            //------------------------------------------
            //-------------------------//
            status = GSL_SUCCESS;
            orbit_update_ic(orbit, sti, orbit.pmap->t0);
            //-------------------------//

            //------------------------------------------
            //If the correction is succesful, we can go on with the computation of the return map
            //------------------------------------------
            //Additional condition on the initial state, if desired
            bool additionalCondition = true; //fabs(SEML.cs_em.gamma*orbit.z0[0]) > 0.019 && fabs(SEML.cs_em.gamma*orbit.z0[1]) < 0.08;
            if(status == GSL_SUCCESS && additionalCondition)
            {
                //------------------------------------------
                // Computation
                //------------------------------------------
                tic();
                //dual_tmap(&orbit);
                status = single_tmap(&orbit);
                //single_tmap_tested(&orbit);
                time = toc();

                //------------------------------------------
                // Postprocess
                //------------------------------------------
                if(status == GSL_SUCCESS)
                {
                    cout << "tmap. Orbit " << label << "/" << numberOfOrbits << " completed. in " << time << "s. " << endl;

                    //------------------------------------------
                    // Print in file
                    //------------------------------------------
                    orbit.label = label++;
                    orbit.tf = pmap.tf;
                    orbit_pmap_fprint_small(&orbit, filename, 1);

                }
            }

            //Memory release
            free_orbit(&orbit);
            free_dvector(sti, 0, 3);

        }
    }

    //Memory release
    free_dvector(gridx, 0,  pmap.gsize);
    free_dvector(gridy, 0,  pmap.gsize);

    //Plot handling at the end
    if(isPlot)
    {
        printf("Press ENTER to close the gnuplot window(s)\n");
        scanf("%c",&ch);

        //Save in EPS format
        gnuplot_cmd(h1, "set terminal postscript eps solid color enhanced");
        gnuplot_cmd(h1, "set output \"pmap.eps\"");
        gnuplot_cmd(h1, "replot");
    }
    gnuplot_close(h1);
    gnuplot_close(h2);
}


//-----------------------------------------------------------------------------------------------------------------------------------------------------
//
//  Print & Read
//
//-----------------------------------------------------------------------------------------------------------------------------------------------------
/**
 *   \brief Print the poincare map of and orbit in a txt file
 **/
void orbit_pmap_fprint(Orbit *orbit, string filename, int append)
{
    if(orbit->int_method == -1)
    {
        cout << "orbit_print. The orbit was not previously computed. No print." << endl;
    }
    else
    {
        //ifstream readStream;
        ofstream myfile;
        string ss1;
        double zEM[6], hz;


        //Open stream
        if(append) myfile.open ((filename+".txt").c_str(), std::ios::app); //for appending at the end of the file
        else myfile.open ((filename+".txt").c_str());
        myfile << std::showpos << setiosflags(ios::scientific)  << setprecision(15);

        //Line labels
        if(!append)
        {
            myfile << "  label,"
                   << "  order,"
                   << "  ofs_order,"
                   << "  x,   y,   z,   px,   py,   pz,"
                   << "  xt,   yt,   zt,   pxt,   pyt,   pzt,"
                   << "  xh,   yh,   zh,   pxh,   pyh,   pzh,"
                   << "  s1,   s2,   s3,   s4,"
                   << "  t,"
                   << "  dHz,"
                   << "  dHw,"
                   << "  dHz-dHw,"
                   << "  pmap.H0,"
                   << "  pmap.dHv,"
                   << "  ePm,"
                   << "  number,"
                   << "  reset" << endl;
        }

        //First line
        NCtoSYS(orbit->pmap->t0, orbit->z0, zEM, (QBCP_L*) orbit->ode_s_6->d->sys->params);
        hz = qbfbp_H(orbit->pmap->t0, zEM, orbit->ode_s_6->d->sys->params);

        myfile << orbit->label << ",  ";
        myfile << orbit->pmap->order << ",  ";
        myfile << orbit->pmap->ofs_order << ",  ";
        //NC
        for(int k = 0; k < 6; k++) myfile << orbit->z0[k] << ",  ";
        //Translated: NC - V
        for(int k = 0; k < 6; k++) myfile << orbit->z0[k]-creal(orbit->V->at(k).evaluate(orbit->pmap->t0*orbit->n, orbit->pmap->ofs_order)) << ",  ";
        //TFR
        for(int k = 0; k < 6; k++) myfile << orbit->zh0[k] << ",  ";
        //RCM
        for(int k = 0; k < 4; k++) myfile << orbit->si[k] << ",  ";
        myfile << orbit->pmap->t0
               << ",  " <<  hz - orbit->pmap->H0
               << ",  " <<  hz - orbit->pmap->H0
               << ",  " <<  0.0
               << ",  " <<  orbit->pmap->H0
               << ",  " <<  orbit->pmap->dHv
               << ",  " <<  0.0
               << ",  " <<  0
               << ",  " <<  orbit->reset_number << endl;

        //Loop on all events
        for(int i = 0; i <= orbit->last_indix; i++)
        {
            myfile << orbit->label << ",  ";
            myfile << orbit->pmap->order << ",  ";
            myfile << orbit->pmap->ofs_order << ",  ";
            //NC
            for(int k = 0; k <6; k++) myfile  << orbit->z0_mat[k][i]  << ",  ";
            //Translated: NC - V
            for(int k = 0; k < 6; k++) myfile << orbit->z0_mat[k][i]-creal(orbit->V->at(k).evaluate(orbit->te_mat[i]*orbit->n, orbit->pmap->ofs_order)) << ",  ";
            //TFR
            for(int k = 0; k <6; k++) myfile  << orbit->zh0_mat[k][i] << ",  ";
            //RCM
            for(int k = 0; k <4; k++) myfile  << orbit->s0_mat[k][i]  << ",  ";
            myfile << orbit->te_mat[i]
                   << ",  " <<  orbit->hz[i] - orbit->pmap->H0
                   << ",  " <<  orbit->hw[i] - orbit->pmap->H0
                   << ",  " <<  orbit->hz[i] - orbit->hw[i]
                   << ",  " <<  orbit->pmap->H0
                   << ",  " <<  orbit->pmap->dHv
                   << ",  " <<  orbit->ePm[i]
                   << ",  " <<  orbit->nevent[i]
                   << ",  " <<  orbit->reset_number << endl;
        }
        myfile.close();
    }
}

/**
 *   \brief Print the poincare map of and orbit in a txt file. Lighter version
 **/
void orbit_pmap_fprint_small(Orbit *orbit, string filename, int append)
{
    if(orbit->int_method == -1)
    {
        cout << "orbit_print. The orbit was not previously computed. No print." << endl;
    }
    else
    {
        //ifstream readStream;
        ofstream myfile;
        string ss1;
        double zEM[6], hz;


        //Open stream
        if(append) myfile.open ((filename+".txt").c_str(), std::ios::app); //for appending at the end of the file
        else myfile.open ((filename+".txt").c_str());
        myfile << std::showpos << setiosflags(ios::scientific)  << setprecision(15);

        //Line labels
        if(!append)
        {
            myfile << "  label,"
                   << "  order,"
                   << "  ofs_order,"
                   << "  x,   y,   z,   px,   py,   pz,"
                   << "  s1,   s2,   s3,   s4,"
                   << "  t,"
                   << "  dHz,"
                   << "  dHw,"
                   << "  dHz-dHw,"
                   << "  pmap.H0,"
                   << "  pmap.dHv,"
                   << "  ePm,"
                   << "  number,"
                   << "  reset" << endl;
        }

        //First line
        NCtoSYS(orbit->pmap->t0, orbit->z0, zEM, (QBCP_L*) orbit->ode_s_6->d->sys->params);
        hz = qbfbp_H(orbit->pmap->t0, zEM, orbit->ode_s_6->d->sys->params);

        myfile << orbit->label << ",  ";
        myfile << orbit->pmap->order << ",  ";
        myfile << orbit->pmap->ofs_order << ",  ";
        //NC
        for(int k = 0; k < 6; k++) myfile << orbit->z0[k] << ",  ";
        //RCM
        for(int k = 0; k < 4; k++) myfile << orbit->si[k] << ",  ";
        myfile << orbit->pmap->t0
               << ",  " <<  hz - orbit->pmap->H0
               << ",  " <<  hz - orbit->pmap->H0
               << ",  " <<  0.0
               << ",  " <<  orbit->pmap->H0
               << ",  " <<  orbit->pmap->dHv
               << ",  " <<  0.0
               << ",  " <<  0
               << ",  " <<  orbit->reset_number << endl;

        //Loop on all events
        for(int i = 0; i <= orbit->last_indix; i++)
        {
            myfile << orbit->label << ",  ";
            myfile << orbit->pmap->order << ",  ";
            myfile << orbit->pmap->ofs_order << ",  ";
            //NC
            for(int k = 0; k <6; k++) myfile  << orbit->z0_mat[k][i]  << ",  ";
            //RCM
            for(int k = 0; k <4; k++) myfile  << orbit->s0_mat[k][i]  << ",  ";
            myfile << orbit->te_mat[i]
                   << ",  " <<  orbit->hz[i] - orbit->pmap->H0
                   << ",  " <<  orbit->hw[i] - orbit->pmap->H0
                   << ",  " <<  orbit->hz[i] - orbit->hw[i]
                   << ",  " <<  orbit->pmap->H0
                   << ",  " <<  orbit->pmap->dHv
                   << ",  " <<  orbit->ePm[i]
                   << ",  " <<  orbit->nevent[i]
                   << ",  " <<  orbit->reset_number << endl;
        }
        myfile.close();
    }
}

/**
 *   \brief Read the poincare map of and orbit in a txt file. DEPRECATED
 *   TO BE ENHANCED WITH R? (More robust to change of columns, etc)
 **/
void orbit_pmap_fread(Orbit *orbit, string filename, int label)
{
    ifstream read;
    string ss1;
    double x1;      //garbage string
    int i1, events;

    //-----------------
    //Open stream
    //-----------------
    read.open((filename+".txt").c_str());
    //Line labels (first line)
    getline(read, ss1);

    //-----------------
    //Search for the right label
    //-----------------
    do
    {
        for(int k = 0; k <6; k++) read >> orbit->z0[k]   >> ss1;
        for(int k = 0; k <6; k++) read >> orbit->zh0[k]  >> ss1;
        for(int k = 0; k <4; k++) read >> orbit->si[k]   >> ss1;
        read >> orbit->pmap->t0     >> ss1;
        read >> x1                  >> ss1;
        read >> x1                  >> ss1;
        read >> x1                  >> ss1;
        read >> orbit->pmap->H0     >> ss1;
        read >> orbit->pmap->dHv    >> ss1;
        read >> x1                  >> ss1;
        read >> i1                  >> ss1;
        read >> x1                  >> ss1;
        read >> orbit->reset_number;
    }
    while(i1 != label && !read.eof());

    //if the label is found
    if(i1 == label)
    {
        orbit->label = i1;
        events = 0;
        while(!read.eof() && events < orbit->pmap->max_events && i1 == label)
        {
            for(int k = 0; k <6; k++) read >> orbit->z0_mat[k][events]  >> ss1;
            for(int k = 0; k <6; k++) read >> orbit->zh0_mat[k][events] >> ss1;
            for(int k = 0; k <4; k++) read >> orbit->s0_mat[k][events]  >> ss1;
            read >> orbit->te_mat[events]                >> ss1;
            read >> x1                                   >> ss1;
            orbit->hz[events] = x1 + orbit->pmap->H0;
            read >> x1                                   >> ss1;
            orbit->hw[events] = x1 + orbit->pmap->H0;
            read >> x1                                   >> ss1;
            read >> orbit->pmap->H0                      >> ss1;
            read >> orbit->pmap->dHv                     >> ss1;
            read >> orbit->ePm[events]                   >> ss1;
            read >> i1                                   >> ss1;
            read >> orbit->nevent[events]                >> ss1;
            read >> orbit->reset_number;
            events++;
        }

        read.close();
        orbit->int_method = DUAL_INT;
        orbit->last_indix = events-1;
        orbit_update_ic(*orbit, orbit->si, orbit->pmap->t0);
        orbit_pmap_fprint(orbit, "test" , 0); //printf for test
    }
    else cout << "orbit_pmap_fread. The desired label was not found in file." << endl;
}

/**
    \brief Print the porecision map of and orbit in a txt file
**/
void orbit_precision_fprint(Orbit *orbit, string filename, int append)
{
    if(orbit->int_method == -1)
    {
        cout << "orbit_print. The orbit was not previously computed. No print." << endl;
    }
    else
    {
        //ifstream readStream;
        ofstream myfile;
        string ss1;
        double zEM[6], hz;

        //First line
        NCtoSYS(orbit->pmap->t0, orbit->z0, zEM, (QBCP_L*) orbit->ode_s_6->d->sys->params);
        hz = qbfbp_H(orbit->pmap->t0, zEM, orbit->ode_s_6->d->sys->params);

        //Open stream
        if(append) myfile.open ((filename+".txt").c_str(), std::ios::app); //for appending at the end of the file
        else myfile.open ((filename+".txt").c_str());
        myfile << std::showpos << setiosflags(ios::scientific)  << setprecision(15);

        //Line labels
        if(!append)
        {
            myfile << "  label,"
                   << "  order,"
                   << "  x,   y,   z,   px,   py,   pz,"
                   << "  xh,   yh,   zh,   pxh,   pyh,   pzh,"
                   << "  s1,   s2,   s3,   s4,"
                   << "  t,"
                   << "  dHz,"
                   << "  dHw,"
                   << "  dHz-dHw,"
                   << "  pmap.H0,"
                   << "  pmap.dHv,"
                   << "  ePm,"
                   << "  eOm" << endl;
        }

        myfile <<  orbit->label << ",  " << orbit->pmap->order << ",  ";
        for(int k = 0; k <6; k++) myfile << orbit->z0[k]  << ",  ";
        for(int k = 0; k <6; k++) myfile << orbit->zh0[k] << ",  ";
        for(int k = 0; k <4; k++) myfile << orbit->si[k]  << ",  ";
        myfile << orbit->pmap->t0
               << ",  " <<  hz - orbit->pmap->H0
               << ",  " <<  hz - orbit->pmap->H0
               << ",  " <<  0.0
               << ",  " <<  orbit->pmap->H0
               << ",  " <<  orbit->pmap->dHv
               << ",  " <<  0.0
               << ",  " <<  orbit->eOm << endl;

        myfile.close();
    }
}

/**
 *   \brief Writing the precision map header in the txt file
 **/
void header_precision_fprint(string filename)
{
    //ifstream readStream;
    ofstream myfile;
    string ss1;
    //Open stream
    myfile.open ((filename+".txt").c_str());
    myfile << "  label,"
           << "  order,"
           << "  x,   y,   z,   px,   py,   pz,"
           << "  xh,   yh,   zh,   pxh,   pyh,   pzh,"
           << "  s1,   s2,   s3,   s4,"
           << "  t,"
           << "  dHz,"
           << "  dHw,"
           << "  dHz-dHw,"
           << "  pmap.H0,"
           << "  pmap.dHv,"
           << "  cmdist.proj,"
           << "  eOm" << endl;
    myfile.close();
}

/**
 *   \brief Writing the poincare map header in the txt file
 **/
void header_pmap_fprint(string filename)
{
    //ifstream readStream;
    ofstream myfile;
    string ss1;
    //Open stream
    myfile.open ((filename+".txt").c_str());
    myfile << "  label,"
           << "  order,"
           << "  ofs_order,"
           << "  x,   y,   z,   px,   py,   pz,"
           << "  xt,   yt,   zt,   pxt,   pyt,   pzt,"
           << "  xh,   yh,   zh,   pxh,   pyh,   pzh,"
           << "  s1,   s2,   s3,   s4,"
           << "  t,"
           << "  dHz,"
           << "  dHw,"
           << "  dHz-dHw,"
           << "  pmap.H0,"
           << "  pmap.dHv,"
           << "  ePm,"
           << "  number,"
           << "  reset" << endl;
    myfile.close();
}

/**
 *   \brief Writing the poincare map header in the txt file. Lighter version
 **/
void header_pmap_fprint_small(string filename)
{
    //ifstream readStream;
    ofstream myfile;
    string ss1;
    //Open stream
    myfile.open ((filename+".txt").c_str());
    myfile << "  label,"
           << "  order,"
           << "  ofs_order,"
           << "  x,   y,   z,   px,   py,   pz,"
           << "  s1,   s2,   s3,   s4,"
           << "  t,"
           << "  dHz,"
           << "  dHw,"
           << "  dHz-dHw,"
           << "  pmap.H0,"
           << "  pmap.dHv,"
           << "  ePm,"
           << "  number,"
           << "  reset" << endl;
    myfile.close();
}

/**
    \brief Print the energy map of and orbit in a txt file
**/
void orbit_energy_fprint(Orbit *orbit, string filename, double hzmax, int append)
{
    if(orbit->int_method == -1)
    {
        cout << "orbit_print. The orbit was not previously computed. No print." << endl;
    }
    else
    {
        //ifstream readStream;
        ofstream myfile;
        string ss1;
        double zEM[6], hz;

        //First line
        NCtoSYS(orbit->pmap->t0, orbit->z0, zEM, (QBCP_L*) orbit->ode_s_6->d->sys->params);
        hz = qbfbp_H(orbit->pmap->t0, zEM, orbit->ode_s_6->d->sys->params);

        if( fabs(hz - orbit->pmap->H0) <= hzmax) //do we need the condition >= 0.0?
        {
            //Open stream
            if(append) myfile.open ((filename+".txt").c_str(), std::ios::app); //for appending at the end of the file
            else myfile.open ((filename+".txt").c_str());
            myfile << std::showpos << setiosflags(ios::scientific)  << setprecision(15);

            //Line labels
            if(!append)
            {
                myfile << "  label,"
                       << "  order,"
                       << "  x,   y,   z,   px,   py,   pz,"
                       << "  xh,   yh,   zh,   pxh,   pyh,   pzh,"
                       << "  s1,   s2,   s3,   s4,"
                       << "  t,"
                       << "  dHz,"
                       << "  dHw,"
                       << "  dHz-dHw,"
                       << "  pmap.H0,"
                       << "  pmap.dHv,"
                       << "  ePm,"
                       << "  eOm" << endl;
            }

            myfile <<  orbit->label << ",  " << orbit->pmap->order << ",  ";
            for(int k = 0; k <6; k++) myfile << orbit->z0[k]  << ",  ";
            for(int k = 0; k <6; k++) myfile << orbit->zh0[k] << ",  ";
            for(int k = 0; k <4; k++) myfile << orbit->si[k]  << ",  ";
            myfile << orbit->pmap->t0
                   << ",  " <<  hz - orbit->pmap->H0
                   << ",  " <<  hz - orbit->pmap->H0
                   << ",  " <<  0.0
                   << ",  " <<  orbit->pmap->H0
                   << ",  " <<  orbit->pmap->dHv
                   << ",  " <<  0.0
                   << ",  " <<  orbit->eOm << endl;

            myfile.close();
        }
    }
}

/**
 *   \brief Writing the energy map header in the txt file
 **/
void header_energy_fprint(string filename)
{
    //ifstream readStream;
    ofstream myfile;
    string ss1;
    //Open stream
    myfile.open ((filename+".txt").c_str());
    myfile << "  label,"
           << "  order,"
           << "  x,   y,   z,   px,   py,   pz,"
           << "  xh,   yh,   zh,   pxh,   pyh,   pzh,"
           << "  s1,   s2,   s3,   s4,"
           << "  t,"
           << "  dHz,"
           << "  dHw,"
           << "  dHz-dHw,"
           << "  pmap.H0,"
           << "  pmap.dHv,"
           << "  ePm,"
           << "  eOm" << endl;
    myfile.close();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------
//
//  Init of one orbit
//
//-----------------------------------------------------------------------------------------------------------------------------------------------------
/**
 *   \brief Update the array st0 with values from orbit.s0 and the value sr. Used in orbit_init_pmap and deltham.
 *
 *          The following patterns are followed (examples):
 *              - if SEML.model = CRTBP and orbit.dim = 2: s2 = sr, s4 = 0
 *                   st0[0] = orbit->si[0];;
 *                   st0[1] = sr;
 *                   st0[2] = orbit->si[2];
 *                   st0[3] = 0.0;
 *
 *             - if SEML.model = QBCP and orbit.dim = 2: s2 = s4 = sr
 *                   st0[0] = sr;
 *                   st0[1] = orbit->si[1];
 *                   st0[2] = sr;
 *                   st0[3] = orbit->si[3];
 **/
void update_s0(Orbit *orbit, double st0[], double sr)
{
    //Set the current state
    if(SEML.model == M_RTBP) //if model = RTBP
    {
        switch(orbit->vdim)
        {
        case 1:
        {
            st0[0] = sr;
            st0[1] = orbit->si[1];
            st0[2] = orbit->si[2];
            st0[3] = 0.0;

            break;
        }

        case 2:
        {
            st0[0] = orbit->si[0];
            st0[1] = sr;
            st0[2] = orbit->si[2];
            st0[3] = 0.0;

            break;
        }

        case 3:
        {
            st0[0] = orbit->si[0];
            st0[1] = orbit->si[1];
            st0[2] = sr;
            st0[3] = 0.0;

            break;
        }

        default:  //as 3
        {
            st0[0] = orbit->si[0];
            st0[1] = orbit->si[1];
            st0[2] = sr;
            st0[3] = 0.0;
            break;
        }
        }
    }
    else //if model = QBCP
    {
        //If fwrk = F_SEM
        //The starting condition z(t0) = 0 is s2 == s4
        //AND there is not root of P36(nt)
        switch(orbit->vdim)
        {
        case 1:
        {
            st0[0] = sr;
            st0[1] = orbit->si[1];
            st0[2] = sr;
            st0[3] = orbit->si[3];
            break;
        }
        case 2:
        {
            st0[0] = orbit->si[0];
            st0[1] = sr;
            st0[2] = orbit->si[2];
            st0[3] = sr;
            break;
        }
        case 3:
        {
            st0[0] = sr;
            st0[1] = orbit->si[1];
            st0[2] = sr;
            st0[3] = orbit->si[3];
            break;
        }
        case 4:
        {
            st0[0] = orbit->si[0];
            st0[1] = sr;
            st0[2] = orbit->si[2];
            st0[3] = sr;
            break;
        }
        default: //as case 3
        {
            st0[0] = sr;
            st0[1] = orbit->si[1];
            st0[2] = sr;
            st0[3] = orbit->si[3];
            break;
        }
        }
    }
}

/**
    \brief Initialize an orbit wrt a Poincare map so that H(orbit.s0) = H(Pmap)
 **/
int orbit_init_pmap(Orbit &orbit, double st0[])
{
    //------------------------------------------
    //Starting by updating the IC
    //------------------------------------------
    orbit_update_ic(orbit, st0, orbit.pmap->t0);

    //------------------------------------------
    //Root finding
    //------------------------------------------
    gsl_function F;                  //called function in the root finding routine
    double s_low, s_high;            //variables for root bracketing
    double Hup, Hdown, r, fy;        //Energy values and root
    int status;
    int itermax = 100;

    //Initialization of F
    F.function = &deltham; //Energy delta with respect to target
    F.params   = &orbit;   //orbit structure contains the parameters

    //Bracketing the root
    s_low  = 0.0;
    s_high = 1e-6*orbit.pmap->gmax;  //arbitrary small value

    //Increasing s_high until a root is found (Hup*Hdown < 0)
    int iter = 0;
    do
    {
        //Evaluating the bracket
        Hup   = deltham (s_high , &orbit);
        Hdown = deltham (s_low, &orbit);
        s_high *= 1.1;
    }
    while(Hup*Hdown > 0 && iter < itermax);

    //Swap the two variables if needed
    if(s_low > s_high)
    {
        s_low  = s_low + s_high;
        s_high = s_low - s_high;
        s_low  = s_low - s_high;
    }

    // TURN OFF GSL ERROR ERROR HANDLER
    // (GSL JUST PRINT ERROR MESSAGE AND KILL THE PROGRAM IF FLAG IS ON)
    gsl_set_error_handler_off();

    //If a root is found, refine root
    if(Hup*Hdown < 0)
    {
        //Setting the solver
        status = gsl_root_fsolver_set (orbit.ode_s_6->s_root, &F, s_low, s_high);
        //Loop
        iter = 0;
        do
        {
            status = gsl_root_fsolver_iterate (orbit.ode_s_6->s_root);         //updating the solver
            r = gsl_root_fsolver_root (orbit.ode_s_6->s_root);                 //updating the root
            fy = deltham(r, &orbit);                                           //Checking convergence
            status = gsl_root_test_residual (fy , orbit.ode_s_6->eps_root);    //Checking convergence
        }
        while (status == GSL_CONTINUE && (++iter)<50);

        if(status == GSL_SUCCESS)
        {
            //Update st0
            update_s0(&orbit, st0, r);
            //Update the orbit once st0 is good
            orbit_update_ic(orbit, st0, orbit.pmap->t0);
            return GSL_SUCCESS;
        }
        else
        {
            //cout <<  "orbit_init_pmap: No refined root was found (1). The orbit is unchanged" << endl;
            return GSL_FAILURE;
        }
    }
    else
    {
        //cout <<  "orbit_init_pmap: No root was found (2). The orbit is unchanged" << endl;
        return GSL_FAILURE;
    }
}

/**
    \brief Update the initial conditions (si, s0, z0 and s0d) of the orbit given an array of initial TFC conditions si
 **/
void orbit_update_ic(Orbit &orbit, const double si[], double t0)
{
    //------------------------------------------
    // 1. Update si
    //------------------------------------------
    for(int p = 0; p < REDUCED_NV; p++) orbit.si[p] = si[p];

    //------------------------------------------
    // 2. Update s0
    //------------------------------------------
    RCMtoCCM(si, orbit.s0, REDUCED_NV);

    //------------------------------------------
    // 2. Update s0d
    //------------------------------------------
    RCMtoCCM8(si, orbit.s0d);

    //------------------------------------------
    // 4. Update z0
    //------------------------------------------
    //z0 = W(si, 0.0)
    RCMtoNCbyTFC(si,
                 t0,
                 orbit.qbcp_l->us.n,
                 orbit.pmap->order,
                 orbit.pmap->ofs_order,
                 *orbit.Wh,
                 *orbit.ofs,
                 *orbit.PC,
                 *orbit.V,
                 orbit.z0,
                 orbit.pmap->isGS);

    //------------------------------------------
    // 5. Update zh0
    //------------------------------------------
    //zh0 = wh(si, 0.0)
    RCMtoTF(si,
            t0,
            orbit.qbcp_l->us.n,
            orbit.pmap->order,
            orbit.pmap->ofs_order,
            *orbit.Wh,
            *orbit.ofs,
            orbit.zh0,
            orbit.pmap->isGS);
}

/**
 *   \brief Computes the hamiltonian at the position st0, in system coordinates and units.
 *   \param pmap a reference to the pmap that carries a set of useful parameters
 *   \param st0 the input state
 *
 *   WARNING: Direct use of CM inside this
 **/
double orbit_ham(Pmap &pmap, double st0[])
{
    //------------------------------------------
    // Inner variables (NC, TFC)
    //------------------------------------------
    double zEM[6];
    double z0d[6];
    Ofsc AUX(OFS_ORDER);

    //------------------------------------------
    // RCM to NC
    //------------------------------------------
    RCMtoNCbyTFC(st0, pmap.t0, SEML.us.n, pmap.order, pmap.ofs_order, CMh, AUX, Mcoc, Vcoc, z0d, pmap.isGS);

    //------------------------------------------
    // Computation
    //------------------------------------------
    NCtoSYS(pmap.t0, z0d, zEM, &SEML);
    return qbfbp_H(pmap.t0, zEM, &SEML);
}

/**
    \brief Computes the difference between an given Ham value and the state configuration defined in the routine update_s0 (see comments therein).
    \param  sr a double to complete the current tested configuration
    \param  params a pointer to the orbit with a given Ham value (why void? the idea is to a have a generic function, but might be useless at this point)
    \return the difference between the two hamiltonians
 **/
double deltham(double sr, void *params)
{
    Orbit *orbit = (Orbit *) params;
    double st0[4];

    //Update st0
    update_s0(orbit, st0, sr);

    //Return the hamiltonian value
    return (orbit_ham(*orbit->pmap, st0) - orbit->pmap->Hv);
}


//--------------------------------------------------------------------------------------------------------------
//
// Steppers with projection
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
                   double omega1,
                   double omega3,
                   int isResetOn)
{
    int status;
    double yvp[6], yvi[6];
    cdouble scp[REDUCED_NV];

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
                     double omega1,
                     double omega3,
                     int isResetOn)
{
    reset_ode_structure(orbit->ode_s_6);
    int status;
    do
    {
        status = gslc_proj_step(orbit, yv, t, t0, t1, ePm, nreset, omega1, omega3, isResetOn);

    }while(fabs(*t)<fabs(t1));

    return status;
}


//-----------------------------------------------------------------------------------------------------------------------------------------------------
//
// Computation of one orbit
//
//-----------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------
// Poincaré maps (pmaps)
//------------------------------------
/**
 *   \brief Computes the poincare map of one given orbit, with a given method
 *   \param orbit a pointer to the orbit
 **/
int orbit_compute_pmap(Orbit *orbit, int method)
{
    int status = GSL_FAILURE;
    //------------------------------------------------
    //Init
    //------------------------------------------------
    double tf = orbit->pmap->tf;
    double t  = orbit->pmap->t0;

    //Starting in the right direction
    orbit->ode_s_6->d->h = (tf>t) ? fabs(orbit->ode_s_6->d->h) : -fabs(orbit->ode_s_6->d->h);
    orbit->ode_s_6->h    = (tf>t) ? fabs(orbit->ode_s_6->h)    : -fabs(orbit->ode_s_6->h);
    orbit->ode_s_8->d->h = (tf>t) ? fabs(orbit->ode_s_8->d->h) : -fabs(orbit->ode_s_8->d->h);
    orbit->ode_s_8->h    = (tf>t) ? fabs(orbit->ode_s_8->h)    : -fabs(orbit->ode_s_8->h);

    //Integration
    switch(method)
    {
    case DUAL_INT:
    {
        //Dual integration
        //--------------------------------
        status = dual_pmap(orbit);
        break;
    }

    case DUAL_INT_STEPPED:
    {
        //Dual integration, with roots refinement during integration
        //--------------------------------
        status = dual_pmap_stepped(orbit);
        //status = dual_pmap_stepped_plot(orbit);
        break;
    }

    case SINGLE_INT:
    {
        //Reset @crossing (needs improvement)
        //--------------------------------
        //status = single_pmap(orbit);
        status = single_pmap_proj(orbit);
        //status = single_pmap_plot(orbit);
        break;
    }

    default:
    {
        //Dual integration
        //--------------------------------
        status = dual_pmap(orbit);
        break;
    }
    }

    //Reset
    reset_ode_structure(orbit->ode_s_6);
    reset_ode_structure(orbit->ode_s_8);

    return status;
}

/**
 *   \brief Computes the poincare map of one given orbit, with a dual method: both vector fields (NC and reduced) are integrated in parallel.
 *          The precise roots of z = 0 are computed at the end of the integration.
 *   \param orbit a pointer to the orbit
 **/
int dual_pmap(Orbit *orbit)
{
    //------------------------------------------------------------------------------
    //Initialization
    //------------------------------------------------------------------------------
    int i, k;                               //loop parameter
    int status;                             //status for GSL function integer output

    //-----------------
    //int & events
    //-----------------
    int events = 0;                         //number of events during integration
    double previous_yv[6], yv[6];           //for event detection purposes, configuration at step n and n-1
    double sv[8];                           //for integration purposes, CCM coordinates (dim = 8)
    double ys[6];                           //for integration purposes, NC coordinates  (dim = 6)
    double yhs[6];                          //for integration purposes, TF coordinates (dim = 6)
    double previous_t, t, tr, tf;           //times (event and integration)
    cdouble zh1[6];                         //for integration purposes, current state in complex form (TFC)
    double z1[6];                           //for integration purposes, current state in real form (NC)
    double previous_si[4], si[4];           //for event detection purposes, configuration at step n and n-1 (RCM)
    double eO[6];                                      //orbit error between z1 and W(s1, t)
    double yv_mat[6][orbit->pmap->max_events+1];         //buffer to store the events (configuration right before the event)
    double sv_mat[4][orbit->pmap->max_events+1];         //buffer to store the events (configuration right before the event)
    double  t_mat[2][orbit->pmap->max_events+1];         //buffer to store the time of the events by pair of 2 (time before and after the event)
    value_output previous_s;                           //value output structure before the event
    value_output new_s;                                //value output structure after  the event

    //-----------------
    //Root finding
    //-----------------
    gsl_function F;                         //called function in the root finding routine
    struct OdeParams params;               //params for F
    double t_low, t_high;                   //times for root bracketing
    double r;                               //root for GSL routine
    double fy;                              //error on the zero for GSL routine
    int iter;                               //cumulated iterations during root finding
    //Initialization of objects for the root finding (at least what can be init at this step).
    params.d      = orbit->ode_s_6->d;
    params.fvalue = orbit->fvalue;
    F.function    = &odezero_event;
    F.params      = &params;

    //------------------------------------------------------------------------------
    //Initialization of the state & time
    //------------------------------------------------------------------------------
    for(i=0; i<6; i++) yv[i] = orbit->z0[i];     //NC coordinates
    for(i=0; i<8; i++) sv[i] = orbit->s0d[i];    //CCM coordinates
    orbit->reset_number      = 0;
    t  = orbit->pmap->t0;
    tr = orbit->pmap->t0;
    tf = orbit->pmap->tf;

    //------------------------------------------------------------------------------
    //Update of the method
    //------------------------------------------------------------------------------
    orbit->int_method = DUAL_INT;

    //------------------------------------------------------------------------------
    //Integration and event detection up to t = tf or the number of events is maximal
    //------------------------------------------------------------------------------
    do
    {
        //Keeping track of the n-1 step
        previous_t = t;
        for(i=0; i<6; i++) previous_yv[i] = yv[i];  //NC
        for(i=0; i<4; i++) previous_si[i] = si[i];  //RCM

        //Evolve one step (dual integration)
        status = gslc_dual_step(orbit, yv, sv, eO, z1, &t, &tr, tf, 1);
        CCM8toRCM(sv, si); //From CCM to RCM

        //Event detection
        previous_s = orbit->fvalue->value(previous_t,previous_yv, orbit->fvalue->val_par);
        new_s      = orbit->fvalue->value(t,yv, orbit->fvalue->val_par);

        if(previous_s.val*new_s.val < 0 && fabs(t -orbit->pmap->t0) > 1e-6)  //a zero of the value function has been crossed
        {
            cout << "dual_pmap. Cross n°" << events+1 << " has been detected at time t ~  " << t << endl;
            //Storage of values
            for(i=0; i<6; i++) yv_mat[i][events] = previous_yv[i]; //NC
            for(i=0; i<4; i++) sv_mat[i][events] = previous_si[i]; //RCM
            t_mat[0][events] = previous_t;
            t_mat[1][events] = t;
            orbit->nevent[events] = events+1;
            //Add one event
            events++;
        }
    }
    while(fabs(t)<fabs(tf) && events < new_s.max_events);

    //Termination detection
    if(events>=new_s.max_events)  printf("dual_pmap. Maximum number of events reached. break.\n");


    //Detect if the final time was reached
    if(fabs(t)>=fabs(tf))
    {
        printf("dual_pmap. Final time was reached, last state is returned at the end of orbit->z0_mat.\n");
        for(i=0; i<6; i++) orbit->z0_mat[i][events] = yv[i];      //Store position (NC)
        for(i=0; i<4; i++) orbit->s0_mat[i][events] = si[i];      //Store position (RCM)
        orbit->te_mat[events] = t;                                //Store time
    }

    //Store the last indix
    if(events >0) orbit->last_indix = events-1;
    else orbit->last_indix          = 0;

    //------------------------------------------------------------------------------
    // Root finding for all events stored in z0_mat
    // For this step, only the NC vector field is integrated for speed
    // The orbit error is not checked since it has been bounded during the whole integration
    //------------------------------------------------------------------------------
    if(events==0)
    {
        printf("dual_pmap. No events was found, last state is returned. \n");
        for(i=0; i<6; i++) orbit->z0_mat[i][0] = yv[i];      //Store position (NC)
        for(i=0; i<4; i++) orbit->s0_mat[i][0] = si[i];      //Store position (RCM)
        orbit->te_mat[0]  = t;
        orbit->last_indix = 0;
    }
    else
    {
        for(k = 0; k < events; k++)
        {
            //Copy of initial values for storage in z0_mat after the root finding
            for(i=0; i<6; i++) yv[i] = yv_mat[i][k];
            for(i=0; i<4; i++) si[i] = sv_mat[i][k];

            //New start for root finding
            for(i=0; i<6; i++) ys[i] = yv_mat[i][k];

            //Initialization of the parameters for F
            params.t0 = t_mat[0][k];       //new initial time is previous_t
            params.y0 = ys;

            //Bracketing the root
            t_low  = (t_mat[0][k] < t_mat[1][k])? t_mat[0][k] : t_mat[1][k];
            t_high = (t_mat[0][k] < t_mat[1][k])? t_mat[1][k] : t_mat[0][k];

            //Setting the solver
            status = gsl_root_fsolver_set (orbit->ode_s_6->s_root, &F, t_low, t_high);

            //Loop
            fy = ys[1];
            iter = 0;
            do
            {
                status = gsl_root_fsolver_iterate (orbit->ode_s_6->s_root);         //updating the solver
                r = gsl_root_fsolver_root (orbit->ode_s_6->s_root);                 //updating the root
                previous_t = gsl_root_fsolver_x_lower (orbit->ode_s_6->s_root);     //updating t_low
                t = gsl_root_fsolver_x_upper (orbit->ode_s_6->s_root);              //updating t_high

                //Checking convergence
                fy = odezero_event (r, &params);
                status = gsl_root_test_residual (fy , orbit->ode_s_6->eps_root);
            }
            while (status == GSL_CONTINUE && (++iter)<100);

            if(iter>=100) cout << "WARNING: number of iter max exceeded in refine_root_dual. Returned precision = " << fy << endl;

            //-----------------------------------
            //Reset to go on with integration
            //-----------------------------------
            reset_ode_structure(orbit->ode_s_6);
            reset_ode_structure(orbit->ode_s_8);

            //Updating the outputs @ crossing time t = r in NC
            t = t_mat[0][k];
            status = gsl_odeiv2_driver_apply(orbit->ode_s_6->d, &t, r, yv);
            //Updating the outputs @ crossing time t = r in CCM/RCM
            tr = t_mat[0][k];
            RCMtoCCM8(si, sv); //from RCM to CCM8
            status = gsl_odeiv2_driver_apply(orbit->ode_s_8->d, &tr, r, sv);

            //-----------------
            //Storage in orbit
            //-----------------
            //NC
            for(i=0; i<6; i++) orbit->z0_mat[i][k] = yv[i];
            //RCM
            CCM8toRCM(sv, si); //from CCM8 to RCM
            for(i=0; i<4; i++) orbit->s0_mat[i][k] = si[i];
            //TFC
            RCMtoTFC(si, t, orbit->qbcp_l->us.n, orbit->pmap->order,  orbit->pmap->ofs_order,  *orbit->Wh,  *orbit->ofs, zh1,  orbit->pmap->isGS);
            TFCtoTF(zh1, yhs); //from TFC to TF
            for(i=0; i<6; i++) orbit->zh0_mat[i][k] = yhs[i];
            //Time (common to NC, RCM and TFC)
            orbit->te_mat[k] = t;

            //-----------------
            //Energy
            //-----------------
            NCtoSYS(t, yv, ys, (QBCP_L*) orbit->ode_s_6->d->sys->params);
            orbit->hz[k] = qbfbp_H(t, ys, orbit->ode_s_6->d->sys->params);
            //Evaluate pm at s1 and store it in yv
            RCMtoNCbyTFC(si, t, orbit->qbcp_l->us.n, orbit->pmap->order,  orbit->pmap->ofs_order,  *orbit->Wh,  *orbit->ofs, Mcoc, Vcoc,  yv,  orbit->pmap->isGS);
            //NC to EM coordinates
            NCtoSYS(t, yv, ys, (QBCP_L*) orbit->ode_s_6->d->sys->params);
            orbit->hw[k] = qbfbp_H(t, ys, orbit->ode_s_6->d->sys->params);

            //------------------------------------------------------------------------------
            //Reset to go on with next crossing
            //------------------------------------------------------------------------------
            reset_ode_structure(orbit->ode_s_6);
            reset_ode_structure(orbit->ode_s_8);
        }

    }

    return GSL_SUCCESS;
}

/**
 *   \brief Computes the poincare map of one given orbit, with a dual method: both vector fields (NC and reduced) are integrated in parallel.
 *          The precise roots of z = 0 are computed during the integration.
 *   \param orbit a pointer to the orbit
 **/
int dual_pmap_stepped(Orbit *orbit)
{
    //------------------------------------------------------------------------------
    //Initialization
    //------------------------------------------------------------------------------
    //-----------------
    //int & events
    //-----------------
    int events = 0;                         //number of events during integration
    double previous_yv[6], yv[6];           //for event detection purposes, configuration at step n and n-1
    double sv[8];                           //for integration purposes, CCM coordinates (dim = 8)
    double yhs[6];                          //for integration purposes, TF coordinates (dim = 6)
    double previous_t, t, tr, tf;           //times (event and integration)
    double z1[6];                           //for integration purposes, current state in real form (NC)
    double previous_si[4], si[4];           //for event detection purposes, configuration at step n and n-1 (RCM)
    double eO[6];                           //orbit error between z1 and W(s1, t)
    double yv_mat[6];                       //buffer to store the events (configuration right before the event)
    double sv_mat[4];                       //buffer to store the events (configuration right before the event)
    double  t_mat[2];                       //buffer to store the time of the events by pair of 2 (time before and after the event)
    value_output previous_s;                //value output structure before the event
    value_output new_s;                     //value output structure after  the event

    double yv_c[6];
    double si_c[4];
    double t_c;
    cdouble s1_c[4];
    double ys_c[6];

    //------------------------------------------------------------------------------
    //Initialization of the state & time
    //------------------------------------------------------------------------------
    for(int i=0; i<6; i++) yv[i] = orbit->z0[i];     //NC coordinates
    for(int i=0; i<8; i++) sv[i] = orbit->s0d[i];    //CCM coordinates
    for(int i=0; i<4; i++) si[i] = orbit->si[i];     //RCM coordinates
    orbit->reset_number      = 0;
    t  = orbit->pmap->t0;
    tr = orbit->pmap->t0;
    tf = orbit->pmap->tf;

    //------------------------------------------------------------------------------
    //Update of the method
    //------------------------------------------------------------------------------
    orbit->int_method = DUAL_INT_STEPPED;

    //------------------------------------------------------------------------------
    //Integration and event detection up to t = tf or the number of events is maximal
    //------------------------------------------------------------------------------
    do
    {
        //Keeping track of the n-1 step
        previous_t = t;
        for(int i=0; i<6; i++) previous_yv[i] = yv[i];  //NC
        for(int i=0; i<4; i++) previous_si[i] = si[i];  //RCM

        //Evolve one step (dual integration)
        gslc_dual_step(orbit, yv, sv, eO, z1, &t, &tr, tf, 1);
        CCM8toRCM(sv, si); //From CCM8 to RCM

        //Event detection
        previous_s = orbit->fvalue->value(previous_t,previous_yv, orbit->fvalue->val_par);
        new_s      = orbit->fvalue->value(t,yv, orbit->fvalue->val_par);

        if(previous_s.val*new_s.val < 0  && fabs(t -orbit->pmap->t0) > 1e-4)  //a zero of the value function has been crossed
        {
            cout << "dual_pmap. Cross n°" << events+1 << " has been detected at time t ~  " << t << endl;

            //-----------------
            //Storage of values
            //-----------------
            for(int i=0; i<6; i++) yv_mat[i] = previous_yv[i]; //NC
            for(int i=0; i<4; i++) sv_mat[i] = previous_si[i]; //RCM
            t_mat[0] = previous_t;
            t_mat[1] = t;
            orbit->nevent[events] = events+1;

            //-----------------
            //Refine root so that yv is really on the z=0 plane
            //-----------------
            refine_root_dual(orbit, yv_c, si_c, s1_c, &t_c, yv_mat, sv_mat, t_mat, events);

            //-----------------
            //Storage in orbit
            //-----------------
            //NC
            for(int i=0; i<6; i++) orbit->z0_mat[i][events] = yv_c[i];
            //RCM
            for(int i=0; i<4; i++) orbit->s0_mat[i][events] = si_c[i];
            //TF
            RCMtoTF(si_c, t_c, orbit->qbcp_l->us.n, orbit->pmap->order, orbit->pmap->ofs_order, *orbit->Wh, *orbit->ofs, yhs, orbit->pmap->isGS);
            for(int i=0; i<6; i++) orbit->zh0_mat[i][events] = yhs[i];
            //Time (common to NC, RCM and TFC)
            orbit->te_mat[events] = t_c;

            //-----------------
            //Energy
            //-----------------
            NCtoSYS(t_c, yv_c, ys_c, (QBCP_L*) orbit->ode_s_6->d->sys->params);
            orbit->hz[events] = qbfbp_H(t_c, ys_c, orbit->ode_s_6->d->sys->params);

            //-----------------
            //Energy (2)
            //-----------------
            RCMtoNCbyTFC(si_c, t, orbit->qbcp_l->us.n, orbit->pmap->order,  orbit->pmap->ofs_order,  *orbit->Wh,  *orbit->ofs, Mcoc, Vcoc,  yv_c,  orbit->pmap->isGS);
            NCtoSYS(t, yv_c, ys_c, (QBCP_L*) orbit->ode_s_6->d->sys->params);
            orbit->hw[events] = qbfbp_H(t, ys_c, orbit->ode_s_6->d->sys->params);
            orbit->nevent[events] = events+1; //number of the event

            //Add one event
            events++;
        }
    }
    while(fabs(t)<fabs(tf) && events < new_s.max_events);

    //------------------------------------------------------------------------------
    //Termination detection
    //------------------------------------------------------------------------------
    //Detect if the maximum number of events was reached
    if(events>=new_s.max_events)  printf("dual_pmap. Maximum number of events reached. break.\n");


    //Detect if the final time was reached
    if(fabs(t)>=fabs(tf))
    {
        printf("dual_pmap. Final time was reached, last state is returned at the end of orbit->z0_mat.\n");
        for(int i=0; i<6; i++) orbit->z0_mat[i][events] = yv[i];      //Store position (NC)
        for(int i=0; i<4; i++) orbit->s0_mat[i][events] = si[i];      //Store position (RCM)
        orbit->te_mat[events] = t;                                    //Store time
    }

    //------------------------------------------------------------------------------
    //Store the last indix
    //------------------------------------------------------------------------------
    if(events >0) orbit->last_indix = events-1;
    else orbit->last_indix          = 0;

    //------------------------------------------------------------------------------
    // Sucess!
    //------------------------------------------------------------------------------
    return GSL_SUCCESS;
}

/**
 *   \brief Computes the poincare map of one given orbit, with a dual method: both vector fields (NC and reduced) are integrated in parallel.
 *          The precise roots of z = 0 are computed during the integration.
 *   \param orbit a pointer to the orbit
 **/
int dual_pmap_stepped_plot(Orbit *orbit)
{
    int Npoints = 10000;
    int iterPlot = 0;
    double yPlot[3][Npoints];
    //------------------------------------------------------------------------------
    //Initialization
    //------------------------------------------------------------------------------
    //-----------------
    //int & events
    //-----------------
    int events = 0;                         //number of events during integration
    double previous_yv[6], yv[6];           //for event detection purposes, configuration at step n and n-1
    double sv[8];                           //for integration purposes, CCM coordinates (dim = 8)
    double yhs[6];                          //for integration purposes, TF coordinates (dim = 6)
    double previous_t, t, tr, tf;           //times (event and integration)
    double z1[6];                           //for integration purposes, current state in real form (NC)
    double previous_si[4], si[4];           //for event detection purposes, configuration at step n and n-1 (RCM)
    double eO[6];                           //orbit error between z1 and W(s1, t)
    double yv_mat[6];                       //buffer to store the events (configuration right before the event)
    double sv_mat[4];                       //buffer to store the events (configuration right before the event)
    double  t_mat[2];                       //buffer to store the time of the events by pair of 2 (time before and after the event)
    value_output previous_s;                //value output structure before the event
    value_output new_s;                     //value output structure after  the event

    double yv_c[6];
    double si_c[4];
    double t_c;
    cdouble s1_c[4];
    double ys_c[6];

    //------------------------------------------------------------------------------
    //Initialization of the state & time
    //------------------------------------------------------------------------------
    for(int i=0; i<6; i++) yv[i] = orbit->z0[i];     //NC coordinates
    for(int i=0; i<8; i++) sv[i] = orbit->s0d[i];    //CCM coordinates
    for(int i=0; i<4; i++) si[i] = orbit->si[i];     //RCM coordinates
    orbit->reset_number      = 0;
    t  = orbit->pmap->t0;
    tr = orbit->pmap->t0;
    tf = orbit->pmap->tf;

    //------------------------------------------------------------------------------
    //Update of the method
    //------------------------------------------------------------------------------
    orbit->int_method = DUAL_INT_STEPPED;

    //------------------------------------------------------------------------------
    //Integration and event detection up to t = tf or the number of events is maximal
    //------------------------------------------------------------------------------
    do
    {
        //Keeping track of the n-1 step
        previous_t = t;
        for(int i=0; i<6; i++) previous_yv[i] = yv[i];  //NC
        for(int i=0; i<4; i++) previous_si[i] = si[i];  //RCM

        //Evolve one step (dual integration)
        gslc_dual_step(orbit, yv, sv, eO, z1, &t, &tr, tf, 1);
        CCM8toRCM(sv, si); //From CCM8 to RCM

        //Event detection
        previous_s = orbit->fvalue->value(previous_t,previous_yv, orbit->fvalue->val_par);
        new_s      = orbit->fvalue->value(t,yv, orbit->fvalue->val_par);

        if(previous_s.val*new_s.val < 0  && fabs(t -orbit->pmap->t0) > 1e-4)  //a zero of the value function has been crossed
        {
            cout << "dual_pmap. Cross n°" << events+1 << " has been detected at time t ~  " << t << endl;

            //-----------------
            //Storage of values
            //-----------------
            for(int i=0; i<6; i++) yv_mat[i] = previous_yv[i]; //NC
            for(int i=0; i<4; i++) sv_mat[i] = previous_si[i]; //RCM
            t_mat[0] = previous_t;
            t_mat[1] = t;
            orbit->nevent[events] = events+1;

            //-----------------
            //Refine root so that yv is really on the z=0 plane
            //-----------------
            refine_root_dual(orbit, yv_c, si_c, s1_c, &t_c, yv_mat, sv_mat, t_mat, events);

            //-----------------
            //Storage in orbit
            //-----------------
            //NC
            for(int i=0; i<6; i++) orbit->z0_mat[i][events] = yv_c[i];
            //RCM
            for(int i=0; i<4; i++) orbit->s0_mat[i][events] = si_c[i];
            //TF
            RCMtoTF(si_c, t_c, orbit->qbcp_l->us.n, orbit->pmap->order, orbit->pmap->ofs_order, *orbit->Wh, *orbit->ofs, yhs, orbit->pmap->isGS);
            for(int i=0; i<6; i++) orbit->zh0_mat[i][events] = yhs[i];
            //Time (common to NC, RCM and TFC)
            orbit->te_mat[events] = t_c;

            //-----------------
            //Energy
            //-----------------
            NCtoSYS(t_c, yv_c, ys_c, (QBCP_L*) orbit->ode_s_6->d->sys->params);
            orbit->hz[events] = qbfbp_H(t_c, ys_c, orbit->ode_s_6->d->sys->params);

            //-----------------
            //Energy (2)
            //-----------------
            RCMtoNCbyTFC(si_c, t, orbit->qbcp_l->us.n, orbit->pmap->order,  orbit->pmap->ofs_order,  *orbit->Wh,  *orbit->ofs, Mcoc, Vcoc,  yv_c,  orbit->pmap->isGS);
            NCtoSYS(t, yv_c, ys_c, (QBCP_L*) orbit->ode_s_6->d->sys->params);
            orbit->hw[events] = qbfbp_H(t, ys_c, orbit->ode_s_6->d->sys->params);
            orbit->nevent[events] = events+1; //number of the event

            //Add one event
            events++;
        }


        //For plotting
        if(iterPlot < Npoints)
        {
            for(int i = 0; i < 3; i++) yPlot[i][iterPlot] = yv[i];
            iterPlot++;
        }
    }
    while(fabs(t)<fabs(tf) && events < new_s.max_events);

    //------------------------------------------------------------------------------
    //Termination detection
    //------------------------------------------------------------------------------
    //Detect if the maximum number of events was reached
    if(events>=new_s.max_events)  printf("dual_pmap. Maximum number of events reached. break.\n");


    //Detect if the final time was reached
    if(fabs(t)>=fabs(tf))
    {
        printf("dual_pmap. Final time was reached, last state is returned at the end of orbit->z0_mat.\n");
        for(int i=0; i<6; i++) orbit->z0_mat[i][events] = yv[i];      //Store position (NC)
        for(int i=0; i<4; i++) orbit->s0_mat[i][events] = si[i];      //Store position (RCM)
        orbit->te_mat[events] = t;                                    //Store time
    }

    //------------------------------------------------------------------------------
    //Store the last indix
    //------------------------------------------------------------------------------
    if(events >0) orbit->last_indix = events-1;
    else orbit->last_indix          = 0;


    //------------------------------------------------------------------------------
    //Plotting devices
    //------------------------------------------------------------------------------
    char ch;            //Used to close the gnuplot windows at the very end of the program
    gnuplot_ctrl  *h1;
    h1 = gnuplot_init();


    gnuplot_setstyle(h1, (char*)"lines");
    gnuplot_set_xlabel(h1, (char*)"x [-]");
    gnuplot_set_ylabel(h1, (char*)"y [-]");
    gnuplot_plot_xyz(h1, yPlot[0], yPlot[1], yPlot[2], iterPlot-1, (char*)"NC coordinates", "lines", "1", "1", 1);


    printf("Press ENTER to close the gnuplot window(s)\n");
    scanf("%c",&ch);
    gnuplot_close(h1);

    //------------------------------------------------------------------------------
    // Sucess!
    //------------------------------------------------------------------------------
    return GSL_SUCCESS;
}

/**
 *   \brief Computes the poincare map of one given orbit, with a single method: only the Nc vector field is computed. The distance with respect to the
 *          central manifold is evaluated each time that z = 0.
 *   \param orbit a pointer to the orbit
 **/
int single_pmap(Orbit *orbit)
{
    //------------------------------------------------------------------------------
    //Initialization
    //------------------------------------------------------------------------------
    //-----------------
    //int & events
    //-----------------
    int events = 0;                         //number of events during integration
    double previous_yv[6], yv[6];           //for event detection purposes, configuration at step n and n-1
    double previous_t, t, tf;               //times (event and integration)
    double s1[4];                           //for integration purposes, real TFC current configuration
    double yv_mat[6];                       //buffer to store the events (configuration right before the event)
    double t_mat[2];                        //buffer to store the time of the events by pair of 2 (time before and after the event)
    value_output previous_s;                //value output structure before the event
    value_output new_s;                     //value output structure after  the event
    double radius;

    //------------------------------------------------------------------------------
    //Initialization of the state & time
    //------------------------------------------------------------------------------
    for(int i = 0; i < 6; i++) yv[i] = orbit->z0[i];
    for(int i = 0; i < 4; i++) s1[i] = orbit->si[i];
    orbit->reset_number = 0;
    tf = orbit->pmap->tf;
    t  = orbit->pmap->t0;

    //------------------------------------------------------------------------------
    //Update of the method
    //------------------------------------------------------------------------------
    orbit->int_method = SINGLE_INT;

    //------------------------------------------------------------------------------
    //Find the projection parameters (€€TODO: make it more general)
    //------------------------------------------------------------------------------
    double omega1 = cimag(Fh[0].getCoef(1,0)->ofs_getCoef(0));
    double omega3 = cimag(Fh[1].getCoef(1,1)->ofs_getCoef(0));

    //------------------------------------------------------------------------------
    //Integration and event detection up to t = tf or the number of events is maximal
    //------------------------------------------------------------------------------
    double status;
    do
    {
        //Keep track of step n-1
        previous_t = t;
        for(int i=0; i<6; i++) previous_yv[i] = yv[i];

        //----------------------
        //Evolve one step of z(t)
        //----------------------
        status = gsl_odeiv2_evolve_apply (orbit->ode_s_6->e, orbit->ode_s_6->c, orbit->ode_s_6->s, &orbit->ode_s_6->sys, &t, tf, &orbit->ode_s_6->h, yv);
        if(status != GSL_SUCCESS) return GSL_FAILURE;

        //event detection
        previous_s = orbit->fvalue->value(previous_t,previous_yv,orbit->fvalue->val_par);
        new_s      = orbit->fvalue->value(t,yv,orbit->fvalue->val_par);

        //a zero of the value function has been crossed
        if(previous_s.val*new_s.val < 0)
        {
            //Storage of values
            for(int i=0; i<6; i++) yv_mat[i] = previous_yv[i];
            t_mat[0] = previous_t;
            t_mat[1] = t;

            //Refine root so that yv is really on the z=0 plane
            status = refine_root(orbit, yv, &t, s1, yv_mat, t_mat, omega1, omega3, events);
            if(status != GSL_SUCCESS) return GSL_FAILURE;


            //Storage in orbit
            for(int i=0; i<6; i++) orbit->z0_mat[i][events] = yv[i];                 //Store position (NC)
            for(int i=0; i<4; i++) orbit->s0_mat[i][events] = s1[i];                 //Store position (TFC)
            orbit->te_mat[events] = t;                                               //Store time
            orbit->nevent[events] = events+1;                                        //Store the number of the event

            //Reset ode structure for next step
            reset_ode_structure(orbit->ode_s_6);

            //Add one event
            events++;
        }

        //Check if the system is diverging
        radius = sqrt(yv[0]*yv[0]+ yv[1]*yv[1] + yv[2]*yv[2]);
        if(radius > orbit->pmap->maxRad)
        {
            cout << "single_pmap. the system is divergent: radius = " << radius << ". break." << endl;
            return GSL_FAILURE;
        }
    }
    while(status == GSL_SUCCESS && fabs(t)<fabs(tf) && events < orbit->pmap->max_events);


    if(fabs(t)>=fabs(tf))
    {
        printf("single_pmap. Final time was reached, last state is returned at the end of ye.\n");
        for(int i=0; i<6; i++)
        {
            orbit->z0_mat[i][events] = yv[i];
            orbit->te_mat[events]    = t;
            //WARNING: no TFC storage at this step!
        }
    }

    //------------------------
    //Last update of the orbit
    //------------------------
    //Last indix in the event storage
    if(events > 0) orbit->last_indix = events-1;
    else orbit->last_indix = 0;

    return GSL_SUCCESS;
}

/**
 *   \brief Computes the poincare map of one given orbit, with a single method: only the Nc vector field is computed. The distance with respect to the
 *          central manifold is evaluated each time that z = 0.
 *   \param orbit a pointer to the orbit
 **/
int single_pmap_proj(Orbit *orbit)
{
    //------------------------------------------------------------------------------
    //Initialization
    //------------------------------------------------------------------------------
    //-----------------
    //int & events
    //-----------------
    int events = 0;                         //number of events during integration
    double previous_yv[6], yv[6];           //for event detection purposes, configuration at step n and n-1
    double previous_t, t, tf;               //times (event and integration)
    double s1[4];                           //for integration purposes, real TFC current configuration
    double yv_mat[6];                       //buffer to store the events (configuration right before the event)
    double t_mat[2];                        //buffer to store the time of the events by pair of 2 (time before and after the event)
    value_output previous_s;                //value output structure before the event
    value_output new_s;                     //value output structure after  the event
    double radius;

    //------------------------------------------------------------------------------
    //Initialization of the state & time
    //------------------------------------------------------------------------------
    for(int i = 0; i < 6; i++) yv[i] = orbit->z0[i];
    for(int i = 0; i < 4; i++) s1[i] = orbit->si[i];
    orbit->reset_number = 0;
    tf = orbit->pmap->tf;
    t  = orbit->pmap->t0;

    //------------------------------------------------------------------------------
    //Update of the method
    //------------------------------------------------------------------------------
    orbit->int_method = SINGLE_INT;

    //------------------------------------------------------------------------------
    //Find the projection parameters (€€TODO: make it more general)
    //------------------------------------------------------------------------------
    double omega1 = cimag(Fh[0].getCoef(1,0)->ofs_getCoef(0));
    double omega3 = cimag(Fh[1].getCoef(1,1)->ofs_getCoef(0));

    //------------------------------------------------------------------------------
    //Integration and event detection up to t = tf or the number of events is maximal
    //------------------------------------------------------------------------------
    //Projection tools
    double ePm;
    int nreset = 1;
    double status;
    do
    {
        //Keep track of step n-1
        previous_t = t;
        for(int i=0; i<6; i++) previous_yv[i] = yv[i];

        //----------------------
        //Evolve one step of z(t)
        //----------------------
        status = gslc_proj_step(orbit, yv, &t, orbit->pmap->t0, tf, &ePm, &nreset, omega1, omega3, 1);
        if(status != GSL_SUCCESS) return GSL_FAILURE;

        //event detection
        previous_s = orbit->fvalue->value(previous_t,previous_yv,orbit->fvalue->val_par);
        new_s      = orbit->fvalue->value(t,yv,orbit->fvalue->val_par);

        //a zero of the value function has been crossed
        if(previous_s.val*new_s.val < 0)
        {
            //Storage of values
            for(int i=0; i<6; i++) yv_mat[i] = previous_yv[i];
            t_mat[0] = previous_t;
            t_mat[1] = t;

            //Refine root so that yv is really on the z=0 plane
            status = refine_root(orbit, yv, &t, s1, yv_mat, t_mat, omega1, omega3, events);
            if(status != GSL_SUCCESS) return GSL_FAILURE;

            //Storage in orbit
            for(int i=0; i<6; i++) orbit->z0_mat[i][events] = yv[i];                 //Store position (NC)
            for(int i=0; i<4; i++) orbit->s0_mat[i][events] = s1[i];                 //Store position (TFC)
            orbit->te_mat[events] = t;                                               //Store time
            orbit->nevent[events] = events+1;                                        //Store the number of the event

            //Reset ode structure for next step
            reset_ode_structure(orbit->ode_s_6);

            //Add one event
            events++;
        }

        //Check if the system is diverging
        radius = sqrt(yv[0]*yv[0]+ yv[1]*yv[1] + yv[2]*yv[2]);
        if(radius > orbit->pmap->maxRad)
        {
            cout << "single_pmap. the system is divergent: radius = " << radius << ". break." << endl;
            return GSL_FAILURE;
        }
    }
    while(status == GSL_SUCCESS && fabs(t)<fabs(tf) && events < orbit->pmap->max_events);


    if(fabs(t)>=fabs(tf))
    {
        printf("single_pmap. Final time was reached, last state is returned at the end of ye.\n");
        for(int i=0; i<6; i++)
        {
            orbit->z0_mat[i][events] = yv[i];
            orbit->te_mat[events]    = t;
            //WARNING: no TFC storage at this step!
        }
    }

    //------------------------
    //Last update of the orbit
    //------------------------
    //Last indix in the event storage
    if(events > 0) orbit->last_indix = events-1;
    else orbit->last_indix = 0;

    return GSL_SUCCESS;
}

/**
 *   \brief Computes the poincare map of one given orbit, with a single method: only the Nc vector field is computed. The distance with respect to the
 *          central manifold is evaluated each time that z = 0.
 *   \param orbit a pointer to the orbit
 **/
int single_pmap_plot(Orbit *orbit)
{
    int Npoints = 10000;
    int iterPlot = 0;
    double yPlot[3][Npoints];

    //------------------------------------------------------------------------------
    //Initialization
    //------------------------------------------------------------------------------
    //-----------------
    //int & events
    //-----------------
    int events = 0;                         //number of events during integration
    double previous_yv[6], yv[6];           //for event detection purposes, configuration at step n and n-1
    double previous_t, t, tf;               //times (event and integration)
    double s1[4];                           //for integration purposes, real TFC current configuration
    double yv_mat[6];                       //buffer to store the events (configuration right before the event)
    double t_mat[2];                        //buffer to store the time of the events by pair of 2 (time before and after the event)
    value_output previous_s;                //value output structure before the event
    value_output new_s;                     //value output structure after  the event
    double radius;

    //------------------------------------------------------------------------------
    //Initialization of the state & time
    //------------------------------------------------------------------------------
    for(int i = 0; i < 6; i++) yv[i] = orbit->z0[i];
    for(int i = 0; i < 4; i++) s1[i] = orbit->si[i];
    orbit->reset_number = 0;
    tf = orbit->pmap->tf;
    t  = orbit->pmap->t0;

    //------------------------------------------------------------------------------
    //Update of the method
    //------------------------------------------------------------------------------
    orbit->int_method = SINGLE_INT;

    //------------------------------------------------------------------------------
    //Find the projection parameters (€€TODO: make it more general)
    //------------------------------------------------------------------------------
    double omega1 = cimag(Fh[0].getCoef(1,0)->ofs_getCoef(0));
    double omega3 = cimag(Fh[1].getCoef(1,1)->ofs_getCoef(0));

    //------------------------------------------------------------------------------
    //Integration and event detection up to t = tf or the number of events is maximal
    //------------------------------------------------------------------------------
    do
    {
        //Keep track of step n-1
        previous_t = t;
        for(int i=0; i<6; i++) previous_yv[i] = yv[i];

        //----------------------
        //Evolve one step of z(t)
        //----------------------
        gsl_odeiv2_evolve_apply (orbit->ode_s_6->e, orbit->ode_s_6->c, orbit->ode_s_6->s, &orbit->ode_s_6->sys, &t, tf, &orbit->ode_s_6->h, yv);

        //event detection
        previous_s = orbit->fvalue->value(previous_t,previous_yv,orbit->fvalue->val_par);
        new_s      = orbit->fvalue->value(t,yv,orbit->fvalue->val_par);

        //a zero of the value function has been crossed
        if(previous_s.val*new_s.val < 0)
        {
            //Storage of values
            for(int i=0; i<6; i++) yv_mat[i] = previous_yv[i];
            t_mat[0] = previous_t;
            t_mat[1] = t;

            //Refine root so that yv is really on the z=0 plane
            refine_root(orbit, yv, &t, s1, yv_mat, t_mat, omega1, omega3, events);

            //Storage in orbit
            for(int i=0; i<6; i++) orbit->z0_mat[i][events] = yv[i];                 //Store position (NC)
            for(int i=0; i<4; i++) orbit->s0_mat[i][events] = s1[i];                 //Store position (TFC)
            orbit->te_mat[events] = t;                                               //Store time
            orbit->nevent[events] = events+1;                                        //Store the number of the event

            //Reset ode structure for next step
            reset_ode_structure(orbit->ode_s_6);

            //Add one event
            events++;
        }

        //Check if the system is diverging
        radius = sqrt(yv[0]*yv[0]+ yv[1]*yv[1] + yv[2]*yv[2]);
        if(radius > orbit->pmap->maxRad)
        {
            cout << "single_pmap. the system is divergent: radius = " << radius << ". break." << endl;
            return GSL_FAILURE;
        }

        //For plotting
        if(iterPlot < Npoints)
        {
            for(int i = 0; i < 3; i++) yPlot[i][iterPlot] = yv[i];
            iterPlot++;
        }
    }
    while(fabs(t)<fabs(tf) && events < orbit->pmap->max_events);


    if(fabs(t)>=fabs(tf))
    {
        printf("single_pmap. Final time was reached, last state is returned at the end of ye.\n");
        for(int i=0; i<6; i++)
        {
            orbit->z0_mat[i][events] = yv[i];
            orbit->te_mat[events]    = t;
            //WARNING: no TFC storage at this step!
        }
    }

    //-------------------------------------------------
    //Plotting devices
    //-------------------------------------------------
    char ch;            //Used to close the gnuplot windows at the very end of the program
    gnuplot_ctrl  *h1;
    h1 = gnuplot_init();


    gnuplot_setstyle(h1, (char*)"lines");
    gnuplot_set_xlabel(h1, (char*)"x [-]");
    gnuplot_set_ylabel(h1, (char*)"y [-]");
    gnuplot_plot_xyz(h1, yPlot[0], yPlot[1], yPlot[2], iterPlot-1, (char*)"NC coordinates", "lines", "1", "1", 1);


    printf("Press ENTER to close the gnuplot window(s)\n");
    scanf("%c",&ch);
    gnuplot_close(h1);


    //------------------------
    //Last update of the orbit
    //------------------------
    //Last indix in the event storage
    if(events > 0) orbit->last_indix = events-1;
    else orbit->last_indix = 0;

    return GSL_SUCCESS;
}

/**
 *   \brief Refine the root z = 0 for a dual integration NC+rvf. Used in dual_pmap_stepped.
 **/
int refine_root_dual(Orbit *orbit,
                     double *yv,
                     double *si,
                     cdouble *s1,
                     double *t,
                     double *yv_mat,
                     double *sv_mat,
                     double *t_mat,
                     int events)
{
    gsl_function F;             //called function in the root finding routine
    struct OdeParams params;    //params for F
    double ys[6];               //new start for root finding
    double sv[8];               //for integration purposes, CCM coordinates (dim = 8)
    double t_low, t_high;       //bracketing the root
    int status;                 //status of the search
    double r;                   //root for GSL routine
    double fy;                  //error on the zero for GSL routine
    int iter;                   //cumulated iterations during root finding
    double ti, tr;              //inner routine times

    //Copy of initial values for storage in ye after the root finding
    for(int i=0; i<6; i++) yv[i] = yv_mat[i];
    for(int i=0; i<4; i++) si[i] = sv_mat[i];

    //New start for root finding
    for(int i=0; i<6; i++) ys[i] = yv_mat[i];

    //Initialization of objects for the root finding
    params.d      = orbit->ode_s_6_root->d;
    params.fvalue = orbit->fvalue;
    F.function    = &odezero_event;
    F.params      = &params;
    params.t0     = t_mat[0];       //new initial time is previous_t - epsilon
    params.y0     = ys;

    //Bracketing the root
    t_low  = (t_mat[0] < t_mat[1])? t_mat[0] : t_mat[1];
    t_high = (t_mat[0] < t_mat[1])? t_mat[1] : t_mat[0];

    //Setting the solver
    gsl_root_fsolver_set (orbit->ode_s_6_root->s_root, &F, t_low, t_high);

    //Loop
    fy = ys[1];
    iter = 0;
    do
    {
        status = gsl_root_fsolver_iterate (orbit->ode_s_6_root->s_root);              //updating the solver
        r      = gsl_root_fsolver_root (orbit->ode_s_6_root->s_root);                 //updating the root
        fy = odezero_event (r, &params);                                              //Checking convergence (1)
        status = gsl_root_test_residual (fy , orbit->ode_s_6_root->eps_root);         //Checking convergence (2)
    }
    while (status == GSL_CONTINUE && (++iter)<100);

    if(iter>=100) cout << "WARNING: number of iter max exceeded in refine_root_dual. Returned precision = " << fy << endl;

    //-----------------
    //Reset ode structure for next step
    //-----------------
    reset_ode_structure(orbit->ode_s_6_root);
    reset_ode_structure(orbit->ode_s_8_root);

    //Updating the outputs @time t = r
    ti = t_mat[0];      //set the time to previous step in order to integrate from ti to r
    status = gsl_odeiv2_driver_apply(orbit->ode_s_6_root->d, &ti, r, yv);  //updating y
    *t = r;                                                                //updating t
    //Updating the outputs @ crossing time t = r in CCM/RCM
    tr = t_mat[0];
    RCMtoCCM8(si, sv); //from RCM to CCM8
    status = gsl_odeiv2_driver_apply(orbit->ode_s_8_root->d, &tr, r, sv);
    CCM8toRCM(sv, si);  //from CCM8 to RCM
    CCM8toCCM(sv, s1);  //from CCM8 to CCM4

    //-----------------
    //Reset ode structure for next step
    //-----------------
    reset_ode_structure(orbit->ode_s_6_root);
    reset_ode_structure(orbit->ode_s_8_root);


    return GSL_SUCCESS;
}


/**
 *   \brief Refine the root z = 0 for a single integration NC+rvf. Used in single_pmap.
 **/
int refine_root(Orbit *orbit,
                double *yv,
                double *t,
                double *s1,
                double *yv_mat,
                double *t_mat,
                double omega1,
                double omega3,
                int events)
{
    gsl_function F;             //called function in the root finding routine
    struct OdeParams params;   //params for F
    double ys[6];               //new start for root finding
    double t_low, t_high;       //bracketing the root
    int status;                 //status of the search
    double r;                   //root for GSL routine
    double fy;                  //error on the zero for GSL routine
    int iter;                   //cumulated iterations during root finding
    double ti;                  //inner routine time

    //Copy of initial values for storage in ye after the root finding
    for(int i=0; i<6; i++) yv[i] = yv_mat[i];

    //New start for root finding
    for(int i=0; i<6; i++) ys[i] = yv_mat[i];

    //Initialization of objects for the root finding
    params.d      = orbit->ode_s_6->d;
    params.fvalue = orbit->fvalue;
    F.function    = &odezero_event;
    F.params      = &params;
    params.t0     = t_mat[0]-1e-15;       //new initial time is previous_t - epsilon
    params.y0     = ys;

    //Bracketing the root
    t_low  = (t_mat[0] < t_mat[1])? t_mat[0] : t_mat[1];
    t_high = (t_mat[0] < t_mat[1])? t_mat[1] : t_mat[0];

    //Setting the solver
    gsl_root_fsolver_set (orbit->ode_s_6->s_root, &F, t_low, t_high);

    //Loop
    fy = ys[1];
    iter = 0;
    do
    {
        status = gsl_root_fsolver_iterate (orbit->ode_s_6->s_root);              //updating the solver
        r      = gsl_root_fsolver_root (orbit->ode_s_6->s_root);                 //updating the root
        //Checking convergence
        fy = odezero_event (r, &params);
        status = gsl_root_test_residual (fy , orbit->ode_s_6->eps_root);
    }
    while (status == GSL_CONTINUE && (++iter)<500);

    if(iter>=500)
    {
        cout << "WARNING: number of iter max exceeded in custom_odezero, with precision: " <<  fy  << "Premature ending." << endl;
        return GSL_FAILURE;
    }

    //Updating the outputs @time t = r
    ti = t_mat[0];                                                    //set the time to previous step in order to integrate from ti to r
    status = gsl_odeiv2_driver_apply(orbit->ode_s_6->d, &ti, r, yv);  //updating y
    *t = r;                                                           //updating t


    //Energy @ yv(t) before projection
    NCtoSYS(*t, yv, ys, (QBCP_L*) orbit->ode_s_6->d->sys->params);
    orbit->hz[events] = qbfbp_H(*t, ys, orbit->ode_s_6->d->sys->params);

    //--------------------------------------------------------
    // Project the state on the center manifold
    //--------------------------------------------------------
    //For comparison, the state before projection is stored in yvi
    double yvi[6];
    for(int i = 0; i <6; i++) yvi[i] = yv[i];

    //Project on the center manifold
    cdouble scp[4];
    NCprojCCM(yv, r, SEML.us_em.n, OFS_ORDER, MIcoc, Vcoc, omega1, omega3, scp, 4);
    CCMtoNCbyTFC(scp, r, orbit->qbcp_l->us.n, orbit->pmap->order,  orbit->pmap->ofs_order,  *orbit->Wh,  *orbit->ofs, Mcoc, Vcoc,  yv,  orbit->pmap->isGS);

    //Compute projected state in RCM coordinates
    CCMtoRCM(scp, s1, REDUCED_NV);

    //Get the current error
    double ePm = fabs(yvi[0] - yv[0]);
    for(int i = 1; i <6 ; i++)
    {
        if(fabs(yvi[i] - yv[i]) > ePm) ePm = fabs(yvi[i] - yv[i]);
    }

    //Update the error in ePm
    orbit->ePm[events] = ePm;

    //Important: forcing exactly z = 0, otherwise a false event may be triggered @next step
    yv[2] = 0.0;

    //Energy @ yv(t) after projection
    NCtoSYS(*t, yv, ys, (QBCP_L*) orbit->ode_s_6->d->sys->params);
    orbit->hw[events] = qbfbp_H(*t, ys, orbit->ode_s_6->d->sys->params);

    return GSL_SUCCESS;
}


//------------------------------------
// Stroboscopic maps (tmaps)
//------------------------------------
/**
 *   \brief Computes the period/stroboscopic/poincare map of one given orbit, with a dual method: both vector fields (NC and reduced) are integrated in parallel.
 *   \param orbit a pointer to the orbit
 **/
int dual_tmap(Orbit *orbit)
{
    //------------------------------------------------------------------------------
    //Initialization
    //------------------------------------------------------------------------------
    int i;                               //loop parameter
    //-----------------
    //int & events
    //-----------------
    int events = 0;         //number of events during integration
    double yv[6], ys[6];    //for event detection purposes, configuration at step n and n-1
    double sv[8];           //for integration purposes, CCM coordinates (dim = 8)
    double yhs[6];          //for integration purposes, TF coordinates (dim = 6)
    double t, tr, ti;       //times (event and integration)
    cdouble zh1[6];         //for integration purposes, current state in complex form (TFC)
    double z1[6];           //for integration purposes, current state in real form (NC)
    double si[4];           //for event detection purposes, configuration at step n and n-1 (RCM)
    double eO[6];           //orbit error between z1 and W(s1, t)


    //------------------------------------------------------------------------------
    //Initialization of the state & time
    //------------------------------------------------------------------------------
    for(i=0; i<6; i++) yv[i] = orbit->z0[i];     //NC coordinates
    for(i=0; i<8; i++) sv[i] = orbit->s0d[i];    //CCM coordinates
    orbit->reset_number      = 0;
    t  = orbit->pmap->t0;
    tr = orbit->pmap->t0;

    //------------------------------------------------------------------------------
    //Update of the method
    //------------------------------------------------------------------------------
    orbit->int_method = DUAL_INT;
    //------------------------------------------------------------------------------
    //Integration and event detection up to t = tf or the number of events is maximal
    //------------------------------------------------------------------------------
    for(int p = 0; p <= orbit->pmap->max_events; p++)
    {
        ti = orbit->pmap->t0 + 1.0*p*orbit->pmap->T;
        do
        {
            //Evolve one step (dual integration)
            gslc_dual_step(orbit, yv, sv, eO, z1, &t, &tr, ti, 1);
        }
        while(fabs(t)<fabs(ti));

        cout << "dual_tmap. Cross n°" << events+1 << " has been detected at time t ~  " << t/orbit->pmap->T << "T" << endl;

        //-----------------
        //Storage in orbit
        //-----------------
        //NC
        for(i=0; i<6; i++) orbit->z0_mat[i][events] = yv[i];
        //RCM
        CCM8toRCM(sv, si); //from CCM8 to RCM
        for(i=0; i<4; i++) orbit->s0_mat[i][events] = si[i];
        //TFC
        RCMtoTFC(si, t, orbit->qbcp_l->us.n, orbit->pmap->order,  orbit->pmap->ofs_order,  *orbit->Wh,  *orbit->ofs, zh1,  orbit->pmap->isGS);
        TFCtoTF(zh1, yhs); //from TFC to TF
        for(i=0; i<6; i++) orbit->zh0_mat[i][events] = yhs[i];
        //Time (common to NC, RCM and TFC)
        orbit->te_mat[events] = t;

        //-----------------
        //Energy
        //-----------------
        NCtoSYS(t, yv, ys, (QBCP_L*) orbit->ode_s_6->d->sys->params);
        orbit->hz[events] = qbfbp_H(t, ys, orbit->ode_s_6->d->sys->params);

        //-----------------
        //Energy (2)
        //-----------------
        RCMtoNCbyTFC(si, t, orbit->qbcp_l->us.n, orbit->pmap->order,  orbit->pmap->ofs_order,  *orbit->Wh,  *orbit->ofs, Mcoc, Vcoc,  yv,  orbit->pmap->isGS);
        NCtoSYS(t, yv, ys, (QBCP_L*) orbit->ode_s_6->d->sys->params);
        orbit->hw[events] = qbfbp_H(t, ys, orbit->ode_s_6->d->sys->params);
        orbit->nevent[events] = events+1; //number of the event

        //-----------------
        // If we are going in the wrong direction
        //-----------------
        //if(events > 0 && orbit->z0_mat[1][events-1] < orbit->z0_mat[1][events])
        //{
        //    break;
        //}

        //Add one event
        events++;


    }

    //Store the last indix
    if(events >0) orbit->last_indix = events-1;
    else orbit->last_indix          = 0;

    return GSL_SUCCESS;
}

/**
 *   \brief Computes the period/stroboscopic/poincare map of one given orbit, with a single method: only the NC vector field is integrated.
 *          The current state is projected on the center manifold after each half-period.
 *   \param orbit a pointer to the orbit
 **/
int single_tmap(Orbit *orbit)
{
    //------------------------------------------------------------------------------
    //Initialization
    //------------------------------------------------------------------------------
    int i;                               //loop parameter
    //-----------------
    //int & events
    //-----------------
    int events = 0;         //number of events during integration
    double yv[6], ys[6];    //for event detection purposes, configuration at step n and n-1
    double yhs[6];          //for integration purposes, TF coordinates (dim = 6)
    double ti, tp;          //times (event and integration)
    cdouble zh1[6];         //for integration purposes, current state in complex form (TFC)
    double si[4];           //for event detection purposes, configuration at step n and n-1 (RCM)


    double yvp[6],  yvi[6];
    cdouble scp[4];
    double radius;


    double omega1 = cimag(Fh[0].getCoef(1,0)->ofs_getCoef(0));
    double omega3 = cimag(Fh[1].getCoef(1,1)->ofs_getCoef(0));

    //------------------------------------------------------------------------------
    //Initialization of the state & time
    //------------------------------------------------------------------------------
    for(i=0; i<6; i++) yv[i]  = orbit->z0[i];     //NC coordinates
    for(i=0; i<6; i++) yvp[i] = orbit->z0[i];     //NC coordinates
    for(i=0; i<4; i++) scp[i] = orbit->s0[i];     //RCM coordinates
    tp = orbit->pmap->t0;
    orbit->reset_number = 0;


    //------------------------------------------------------------------------------
    //Update of the method
    //------------------------------------------------------------------------------
    orbit->int_method = DUAL_INT;
    //------------------------------------------------------------------------------
    //Integration and event detection up to t = tf or the number of events is maximal
    //------------------------------------------------------------------------------
    int steps = 5;
    double stepper = 1.0/steps;
    for(int p = 1; p <= orbit->pmap->max_events; p++)
    {
        if(p > 0)
        {
            //Inner loop for projection procedure
            for(int k = 1; k <= steps; k++)
            {
                ti = orbit->pmap->t0 + (p-1+stepper*k)*orbit->pmap->T;
                //------------------------------------------------------------------------------
                //Evolve the yv integration
                //------------------------------------------------------------------------------
                gsl_odeiv2_driver_apply (orbit->ode_s_6->d, &tp, ti, yv);
                //Project on the center manifold
                NCprojCCM(yv, ti, SEML.us_em.n, OFS_ORDER, MIcoc, Vcoc, omega1, omega3, scp, 4);
                CCMtoNCbyTFC(scp, ti, orbit->qbcp_l->us.n, orbit->pmap->order,  orbit->pmap->ofs_order,  *orbit->Wh,  *orbit->ofs, Mcoc, Vcoc,  yvp,  orbit->pmap->isGS);
                //Reset ode structure for next step
                reset_ode_structure(orbit->ode_s_6);
                //------------------------------------------------------------------------------
                //For comparison
                //------------------------------------------------------------------------------
                for(int i = 0; i <6; i++) yvi[i] = yv[i];
                //------------------------------------------------------------------------------
                // Copy of yvp in current state
                //------------------------------------------------------------------------------
                for(int i=0; i<6; i++) yv[i]  = yvp[i];

                //------------------------------------------------------------------------------
                //Check if the system is diverging
                //------------------------------------------------------------------------------
                radius = sqrt(yv[0]*yv[0]+ yv[1]*yv[1] + yv[2]*yv[2]);
                if(radius > orbit->pmap->maxRad)
                {
                    cout << "single_tmap. the system is divergent: radius = " << radius << ". break." << endl;
                    return GSL_FAILURE;
                }
            }
        }

        //-----------------
        //Storage in orbit
        //-----------------
        //NC
        for(i=0; i<6; i++) orbit->z0_mat[i][events] = yv[i];
        //RCM
        CCMtoRCM(scp, si, REDUCED_NV); //from CCM to RCM
        for(i=0; i<4; i++) orbit->s0_mat[i][events] = si[i];
        //TFC
        RCMtoTFC(si, tp, orbit->qbcp_l->us.n, orbit->pmap->order,  orbit->pmap->ofs_order,  *orbit->Wh,  *orbit->ofs, zh1,  orbit->pmap->isGS);
        TFCtoTF(zh1, yhs); //from TFC to TF
        for(i=0; i<6; i++) orbit->zh0_mat[i][events] = yhs[i];
        //Time (common to NC, RCM and TFC)
        orbit->te_mat[events] = tp;

        //-----------------
        //Energy (1)
        //-----------------
        NCtoSYS(tp, yvi, ys, (QBCP_L*) orbit->ode_s_6->d->sys->params);
        orbit->hz[events] = qbfbp_H(tp, ys, orbit->ode_s_6->d->sys->params);

        //-----------------
        // Get the current projection error
        //-----------------
        //Get the current error
        double ePm = fabs(yvi[0] - yv[0]);
        for(int i = 1; i <6 ; i++)
        {
            if(fabs(yvi[i] - yv[i]) > ePm) ePm = fabs(yvi[i] - yv[i]);
        }

        //Update the error in ePm
        orbit->ePm[events] = ePm;

        //-----------------
        //Energy (2)
        //-----------------
        NCtoSYS(tp, yv, ys, (QBCP_L*) orbit->ode_s_6->d->sys->params);
        orbit->hw[events] = qbfbp_H(tp, ys, orbit->ode_s_6->d->sys->params);
        orbit->nevent[events] = events+1; //number of the event

        //Add one event
        events++;
    }

    //Store the last indix
    if(events >0) orbit->last_indix = events-1;
    else orbit->last_indix          = 0;


    return GSL_SUCCESS;
}

/**
 *   \brief Computes the period/stroboscopic/poincare map of one given orbit, with a single method: only the NC vector field is integrated.
 *          The current state is projected on the center manifold after each half-period.
 *   \param orbit a pointer to the orbit
 **/
int single_tmap_tested(Orbit *orbit)
{
    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(15);
    //------------------------------------------------------------------------------
    //Initialization
    //------------------------------------------------------------------------------
    int i;                               //loop parameter
    //-----------------
    //int & events
    //-----------------
    int events = 0;         //number of events during integration
    double yv[6], ys[6];    //for event detection purposes, configuration at step n and n-1
    double sv[8];           //for integration purposes, CCM coordinates (dim = 8)
    double yhs[6];          //for integration purposes, TF coordinates (dim = 6)
    double tr, ti;       //times (event and integration)
    cdouble zh1[6];         //for integration purposes, current state in complex form (TFC)
    double si[4];           //for event detection purposes, configuration at step n and n-1 (RCM)


    double yvp[6],  yve[6], yvc[6];
    double tp;
    cdouble scp[4];


    double omega1 = cimag(Fh[0].getCoef(1,0)->ofs_getCoef(0));
    double omega3 = cimag(Fh[1].getCoef(1,1)->ofs_getCoef(0));

    //------------------------------------------------------------------------------
    //Initialization of the state & time
    //------------------------------------------------------------------------------
    for(i=0; i<6; i++) yv[i]  = orbit->z0[i];     //NC coordinates
    for(i=0; i<6; i++) yvp[i] = orbit->z0[i];     //NC coordinates
    for(i=0; i<8; i++) sv[i]  = orbit->s0d[i];    //CCM coordinates
    tr = orbit->pmap->t0;
    tp = orbit->pmap->t0;

    orbit->reset_number      = 0;


    //------------------------------------------------------------------------------
    //Update of the method
    //------------------------------------------------------------------------------
    orbit->int_method = DUAL_INT;
    //------------------------------------------------------------------------------
    //Integration and event detection up to t = tf or the number of events is maximal
    //------------------------------------------------------------------------------
    int steps = 20;
    double stepper = 1.0/steps;
    for(int p = 0; p <= orbit->pmap->max_events; p++)
    {
        if(p > 0)
        {
            //Inner loop for projection procedure
            for(int k = 1; k <= steps; k++)
            {
                ti = orbit->pmap->t0 + (p-1+stepper*k)*orbit->pmap->T;
                //------------------------------------------------------------------------------
                //Evolve the sv integration
                //------------------------------------------------------------------------------
                gsl_odeiv2_driver_apply (orbit->ode_s_8->d, &tr, ti, sv);
                //Copy in yvc
                CCM8toNCbyTFC(sv, ti, orbit->qbcp_l->us.n, orbit->pmap->order,  orbit->pmap->ofs_order,  *orbit->Wh,  *orbit->ofs, Mcoc, Vcoc,  yvc,  orbit->pmap->isGS);
                //------------------------------------------------------------------------------
                //Evolve the yv integration
                //------------------------------------------------------------------------------
                gsl_odeiv2_driver_apply (orbit->ode_s_6->d, &tp, ti, yv);
                //Project on the center manifold
                NCprojCCM(yv, ti, SEML.us_em.n, OFS_ORDER, MIcoc, Vcoc, omega1, omega3, scp, 4);
                CCMtoNCbyTFC(scp, ti, orbit->qbcp_l->us.n, orbit->pmap->order,  orbit->pmap->ofs_order,  *orbit->Wh,  *orbit->ofs, Mcoc, Vcoc,  yvp,  orbit->pmap->isGS);
                //Reset ode structure for next step
                reset_ode_structure(orbit->ode_s_6);

                //------------------------------------------------------------------------------
                //Comparison
                //------------------------------------------------------------------------------
                if(k%10 == 0)
                {
                    cout << "-----------------------------------" << endl;
                    cout << "Comparison single/dual integration:" << endl;
                    cout << "-----------------------------------" << endl;
                    cout << "yv[i]     yvp[i]     yvc[i]     (yvp - yvc)  (yv - yvc)" << endl;
                    for(i=0; i<6; i++)
                    {
                        cout << setprecision(5) << yv[i] << "    " << yvp[i]  << "    "
                             << yvc[i] << "    " << setprecision(2) << yvp[i]-yvc[i] << "    "
                             << yv[i]-yvc[i]<< endl;
                    }
                }
                //------------------------------------------------------------------------------
                // Copy of yvp in current state
                //------------------------------------------------------------------------------
                for(i=0; i<6; i++) yv[i]  = yvp[i];
            }
        }
        cout << "dual_tmap. Cross n°" << events+1 << " has been detected at time t ~  " << tr/orbit->pmap->T << "T" << endl;

        //-----------------
        //Storage in orbit
        //-----------------
        //NC
        for(i=0; i<6; i++) orbit->z0_mat[i][events] = yv[i];
        //RCM
        CCMtoRCM(scp, si, REDUCED_NV); //from CCM to RCM
        for(i=0; i<4; i++) orbit->s0_mat[i][events] = si[i];
        //TFC
        RCMtoTFC(si, tr, orbit->qbcp_l->us.n, orbit->pmap->order,  orbit->pmap->ofs_order,  *orbit->Wh,  *orbit->ofs, zh1,  orbit->pmap->isGS);
        TFCtoTF(zh1, yhs); //from TFC to TF
        for(i=0; i<6; i++) orbit->zh0_mat[i][events] = yhs[i];
        //Time (common to NC, RCM and TFC)
        orbit->te_mat[events] = tr;

        //-----------------
        //Energy
        //-----------------
        NCtoSYS(tr, yv, ys, (QBCP_L*) orbit->ode_s_6->d->sys->params);
        orbit->hz[events] = qbfbp_H(tr, ys, orbit->ode_s_6->d->sys->params);

        //-----------------
        //Energy (2)
        //-----------------
        RCMtoNCbyTFC(si, tr, orbit->qbcp_l->us.n, orbit->pmap->order,  orbit->pmap->ofs_order,  *orbit->Wh,  *orbit->ofs, Mcoc, Vcoc,  yve,  orbit->pmap->isGS);
        NCtoSYS(tr, yv, ys, (QBCP_L*) orbit->ode_s_6->d->sys->params);
        orbit->hw[events] = qbfbp_H(tr, ys, orbit->ode_s_6->d->sys->params);
        orbit->nevent[events] = events+1; //number of the event

        //-----------------
        // If we are going in the wrong direction
        //-----------------
        //if(events > 0 && orbit->z0_mat[1][events-1] < orbit->z0_mat[1][events])
        //{
        //    break;
        //}

        //Add one event
        events++;
    }

    //Store the last indix
    if(events >0) orbit->last_indix = events-1;
    else orbit->last_indix          = 0;


    return GSL_SUCCESS;
}

//------------------------------------
// Error maps (emaps)
//------------------------------------
/**
 *   \brief Computes the error/precision map of one given orbit, with a dual method: both vector fields (NC and reduced) are integrated in parallel.
 *   \param orbit a pointer to the orbit
 **/
int dual_emap(Orbit *orbit)
{
    //------------------------------------------------------------------------------
    //Initialization
    //------------------------------------------------------------------------------
    int i;                  //loop parameter
    //-----------------
    //int & events
    //-----------------
    double yv[6];           //for event detection purposes, configuration at step n and n-1
    double sv[8];           //for integration purposes, CCM coordinates (dim = 8)
    double t, tr, tf;       //times (event and integration)
    double z1[6];          //for integration purposes, current state in real form (NC)
    double eO[6];           //orbit error between z1 and W(s1, t)

    //------------------------------------------------------------------------------
    //Initialization of the state & time
    //------------------------------------------------------------------------------
    for(i=0; i<6; i++) yv[i] = orbit->z0[i];     //NC coordinates
    for(i=0; i<8; i++) sv[i] = orbit->s0d[i];    //CCM coordinates
    orbit->reset_number      = 0;
    t  = orbit->pmap->t0;
    tr = orbit->pmap->t0;
    tf = orbit->pmap->t0+orbit->pmap->tt;

    //------------------------------------------------------------------------------
    //Update of the method
    //------------------------------------------------------------------------------
    orbit->int_method = DUAL_INT_NO_RESET;
    //------------------------------------------------------------------------------
    //Integration and event detection up to t = tf or the number of events is maximal
    //------------------------------------------------------------------------------
    do
    {
        //Evolve one step (dual integration)
        orbit->eOm = gslc_dual_step(orbit, yv, sv, eO, z1, &t, &tr, tf, 0); //no reset at this step!
    }
    while(fabs(t)<fabs(tf));

    return GSL_SUCCESS;
}


//--------------------------------------------------------------------------------------------------------------
//
// Steppers
//
//--------------------------------------------------------------------------------------------------------------
/**
 *   \brief Integrates the current state yv/sv one step, using both NC and rvf vector fields.
 **/
double gslc_dual_step( Orbit *orbit,
                       double yv[],
                       double sv[],
                       double eO[],
                       double zr1[],
                       double *t,
                       double *tr,
                       double t1,
                       int resetOn)
{
    int status;
    double eOm;
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
    //Evolving s(t) up to t
    //----------------------
    status = gsl_odeiv2_driver_apply (orbit->ode_s_8->d, tr, *t, sv);
    if (status != GSL_SUCCESS)
    {
        cout << "error in gslc_dual_step: integration of s(t) has gone wrong. break." << endl;
        return GSL_FAILURE;
    }

    //----------------------
    //z1 = W(s1);
    //----------------------
    CCM8toNCbyTFC(sv, *t, orbit->qbcp_l->us.n, orbit->pmap->order, orbit->pmap->ofs_order, *orbit->Wh, *orbit->ofs, Mcoc, Vcoc, zr1, orbit->pmap->isGS);

    //----------------------
    //Orbit error (infty norm)
    //----------------------
    eOm = fabs(zr1[0] - yv[0]);
    for(int p = 1; p < NV; p++) if(fabs(zr1[p] - yv[p]) > eOm) eOm = fabs(zr1[p] - yv[p]);
    orbit->eOm = eOm;

    //----------------------
    //Reset if necessary
    //----------------------
    if(resetOn && eOm > orbit->pmap->threshold)
    {
        //cout << "Reset is detected @t  = " << *t << ", z = " << creal(z1[2]) << endl;
        for(int p = 0; p < NV; p++) yv[p] = zr1[p];
        (orbit->reset_number)++;
        //Reset ode structure
        reset_ode_structure(orbit->ode_s_6);
    }

    return eOm;
}


/**
 *   \brief Integrates the current state yv/sv up to t = t1, using both NC and rvf vector fields.
 **/
int gslc_dual_evolve(Orbit *orbit,
                     double yv[],
                     double sv[],
                     double eO[],
                     double z1[],
                     double *t,
                     double t1,
                     double threshold)
{
    int status;
    double tr = *t;
    do
    {
        status = gslc_dual_step(orbit, yv, sv, eO, z1, t, &tr, t1, 1);
    }
    while(fabs(*t)<fabs(t1));

    return status;
}



//-----------------------------------------------------------------------------------------------------------------------------------------------------
//
// Plotting (deprecated)
//
//-----------------------------------------------------------------------------------------------------------------------------------------------------
//Plot one orbit
void orbit_plot(Orbit *orbit, gnuplot_ctrl *h1, int type, int points, OdeStruct *ode_s_6, OdeStruct *ode_s_8)
{
    //Checking initialization of h1
    if(h1 == NULL) h1 = gnuplot_init();

    //Init
    double t, t1, ti;
    double x1[points+1];
    double x2[points+1];
    double xr1[points+1];
    double xr2[points+1];
    double xe1[orbit->last_indix+1];
    double xe2[orbit->last_indix+1];

    double y[6];
    double s[8];
    double eO[6];
    int i;

    //Initial conditions in NC and TFC
    for(i=0; i<6; i++) y[i]  = orbit->z0[i];
    for(i=0; i<8; i++) s[i]  = orbit->s0d[i];


    double z1[6];
    //------------------------------------------
    //First value
    //------------------------------------------
    CCM8toNCbyTFC(s, 0.0, orbit->qbcp_l->us.n, orbit->pmap->order, orbit->pmap->ofs_order, *orbit->Wh, *orbit->ofs, Mcoc, Vcoc, z1, orbit->pmap->isGS);

    switch(type)
    {
    case XY:
        x1[0] = y[0];
        x2[0] = y[1];

        for(int i = 0; i <= orbit->last_indix; i++)
        {
            xe1[i] = orbit->z0_mat[0][i];
            xe2[i] = orbit->z0_mat[1][i];
        }

        xr1[0] = creal(z1[0]);
        xr2[0] = creal(z1[1]);


        break;
    case XZ:
        x1[0] = y[0];
        x2[0] = y[2];

        for(int i = 0; i <= orbit->last_indix; i++)
        {
            xe1[i] = orbit->z0_mat[0][i];
            xe2[i] = orbit->z0_mat[2][i];
        }

        xr1[0] = creal(z1[0]);
        xr2[0] = creal(z1[2]);


        break;
    case YZ:
        x1[0] = y[1];
        x2[0] = y[2];

        for(int i = 0; i <= orbit->last_indix; i++)
        {
            xe1[i] = orbit->z0_mat[1][i];
            xe2[i] = orbit->z0_mat[2][i];
        }

        xr1[0] = creal(z1[1]);
        xr2[0] = creal(z1[2]);

        break;
    }

    //------------------------------------------
    //Store for plotting
    //------------------------------------------
    //Integration over a period
    t = 0;
    t1 = orbit->tf;
    cout << "orbit_plot. t1 = " << t1 << endl;
    //Starting in the right direction
    ode_s_6->d->h = (t1>t) ? fabs(ode_s_6->d->h) : -fabs(ode_s_6->d->h);
    for(i =1; i<=points; i++)
    {
        ti = i * t1 / points;
        //gsl_odeiv2_driver_apply (ode_s_6->d, &t, ti, y);
        //gsl_odeiv2_driver_apply (ode_s_8->d, &tr, ti, s);
        gslc_dual_evolve(orbit, y, s, eO, z1, &t, ti, orbit->pmap->threshold);


        //current condition updated
        CCM8toNCbyTFC(s, ti, orbit->qbcp_l->us.n, orbit->pmap->order, orbit->pmap->ofs_order, *orbit->Wh, *orbit->ofs, Mcoc, Vcoc, z1, orbit->pmap->isGS);

        switch(type)
        {
        case XY:
            x1[i] = y[0];
            x2[i] = y[1];
            xr1[i] = creal(z1[0]);
            xr2[i] = creal(z1[1]);
            break;
        case XZ:
            x1[i] = y[0];
            x2[i] = y[2];
            xr1[i] = creal(z1[0]);
            xr2[i] = creal(z1[2]);
            break;
        case YZ:
            x1[i] = y[1];
            x2[i] = y[2];
            xr1[i] = creal(z1[1]);
            xr2[i] = creal(z1[2]);
            break;
        }
    }

    //Plotting
    switch(type)
    {
    case XY:
        gnuplot_set_xlabel (h1, (char*) "X");
        gnuplot_set_ylabel (h1, (char*) "Y");
        break;
    case XZ:
        gnuplot_set_xlabel (h1, (char*) "X");
        gnuplot_set_ylabel (h1, (char*) "Z");
        break;
    case YZ:
        gnuplot_set_xlabel (h1, (char*) "Y");
        gnuplot_set_ylabel (h1, (char*) "Z");
        break;
    }
    gnuplot_plot_xy(h1, x1, x2, points+1, (char*)"Orbit", "points", "1", "1", 2);
    gnuplot_plot_xy(h1, xr1, xr2, points+1, (char*)"Orbit", "lines", "dashed", "4", 3);
    gnuplot_plot_xy(h1, xe1, xe2, orbit->last_indix+1, (char*)"Events", "points", "1", "2", 2);

    //Reset
    reset_ode_structure(ode_s_6);
}

//Plot one orbit (3D)
void orbit_plot_3d(Orbit *orbit, gnuplot_ctrl *h1, int points, OdeStruct *ode_s_6, OdeStruct *ode_s_8)
{
    //Checking initialization of h1
    if(h1 == NULL) h1 = gnuplot_init();

    //Init
    double t, t1, ti;
    double x1[points+1];
    double x2[points+1];
    double x3[points+1];
    double xe1[orbit->last_indix+1];
    double xe2[orbit->last_indix+1];
    double xe3[orbit->last_indix+1];
    double y[6];
    double s[8];
    double eO[6];
    int i;

    //Initial conditions in NC and TFC
    for(i=0; i<6; i++) y[i]  = orbit->z0[i];
    for(i=0; i<8; i++) s[i]  = orbit->s0d[i];
    double z1[6];


    //------------------------------------------
    //First value
    //------------------------------------------
    //Initial condition updated
    //------------------------------------------
    x1[0] = y[0];
    x2[0] = y[1];
    x3[0] = y[2];

    for(int i = 0; i <= orbit->last_indix; i++)
    {
        xe1[i] = orbit->z0_mat[0][i];
        xe2[i] = orbit->z0_mat[1][i];
        xe3[i] = orbit->z0_mat[2][i];
    }

    //------------------------------------------
    //Store for plotting
    //------------------------------------------
    //Integration over a period
    t = 0;
    double tr = 0;
    t1 = orbit->tf;
    //Starting in the right direction
    ode_s_6->d->h = (t1>t) ? fabs(ode_s_6->d->h) : -fabs(ode_s_6->d->h);
    for(i =1; i<=points; i++)
    {
        ti = (double) i * t1 / points;
        cout << "ti/tf = " << ti/t1 << endl;
        do
        {
            gslc_dual_step(orbit, y, s, eO, z1, &t, &tr, ti, 1);
        }
        while(fabs(t) < fabs(ti));

        //current condition updated
        CCM8toNCbyTFC(s, ti, orbit->qbcp_l->us.n, orbit->pmap->order, orbit->pmap->ofs_order, *orbit->Wh, *orbit->ofs, Mcoc, Vcoc, z1, orbit->pmap->isGS);

        x1[i] = y[0];
        x2[i] = y[1];
        x3[i] = y[2];
    }

    gnuplot_set_xlabel (h1, (char*) "X");
    gnuplot_set_ylabel (h1, (char*) "Y");
    gnuplot_set_zlabel (h1, (char*) "Z");
    gnuplot_setstyle(h1, (char*)  "lines");
    gnuplot_setcolor(h1, 2);
    gnuplot_plot_xyz(h1, x1, x2, x3, points+1, (char*)"Orbit", "lines", "1", "1", 2);
    gnuplot_plot_xyz(h1, xe1, xe2, xe3, orbit->last_indix+1, (char*)"Events", "points", "1", "2", 3);
    gnuplot_plot_xyz(h1, orbit->z0, orbit->z0+1, orbit->z0+2, 1, (char*)"IC", "points", "1", "2", 4);


    gnuplot_fplot_xyz(x1, x2, x3, points+1, (char*)"gnuplot/orbit.txt");
    gnuplot_fplot_xyz(xe1, xe2, xe3, orbit->last_indix+1, (char*)"gnuplot/events.txt");

    //Reset
    reset_ode_structure(ode_s_6);
}

//Plot poincare map for one orbit
void orbit_poincare_plot(Orbit *orbit, gnuplot_ctrl *h1, gnuplot_ctrl *h2, int color)
{
    //Checking initialization of h1 & h2
    if(h1 == NULL) h1 = gnuplot_init();
    if(h2 == NULL) h2 = gnuplot_init();

    //Init
    double xe1[orbit->last_indix+1];
    double xe2[orbit->last_indix+1];

    double xs1[orbit->last_indix+1];
    double xs2[orbit->last_indix+1];

    for(int i = 0; i <= orbit->last_indix; i++)
    {
        xe1[i] = orbit->z0_mat[0][i];
        xe2[i] = orbit->z0_mat[1][i];

        xs1[i] = orbit->zh0_mat[0][i];
        xs2[i] = orbit->zh0_mat[3][i];
    }

    gnuplot_set_xlabel (h1, (char*) "X");
    gnuplot_set_ylabel (h1, (char*) "Y");
    gnuplot_cmd(h1,  "set title \"Semi-Poincare map (NC) \"");
    gnuplot_plot_xy(h1, xe1, xe2, orbit->last_indix+1, (char*)"", "points", "1", "2", color);

    gnuplot_set_xlabel (h2, (char*) "x1");
    gnuplot_set_ylabel (h2, (char*) "y1");
    gnuplot_cmd(h2,  "set title \"Semi-Poincare map (TF) \"");
    gnuplot_plot_xy(h2, xs1, xs2, orbit->last_indix+1, (char*)"", "points", "1", "2", color);
}

//Plot T map for one orbit
void orbit_Tmap_plot(Orbit *orbit, gnuplot_ctrl *h1, gnuplot_ctrl *h2, int color)
{
    //Checking initialization of h1 & h2
    if(h1 == NULL) h1 = gnuplot_init();
    if(h2 == NULL) h2 = gnuplot_init();

    //Init
    double xe1[orbit->last_indix+1];
    double xe2[orbit->last_indix+1];
    double xe3[orbit->last_indix+1];

    double xs1[orbit->last_indix+1];
    double xs2[orbit->last_indix+1];
    double xs3[orbit->last_indix+1];


    for(int i = 0; i <= orbit->last_indix; i++)
    {
        xe1[i] = orbit->z0_mat[0][i];
        xe2[i] = orbit->z0_mat[1][i];
        xe3[i] = orbit->z0_mat[2][i];

        xs1[i] = orbit->zh0_mat[0][i];
        xs2[i] = orbit->zh0_mat[1][i];
        xs3[i] = orbit->zh0_mat[2][i];
    }


    gnuplot_set_xlabel (h1, (char*) "X");
    gnuplot_set_ylabel (h1, (char*) "Y");
    gnuplot_set_zlabel (h1, (char*) "Z");
    gnuplot_cmd(h1,  "set title \"Tmap (NC) \"");
    gnuplot_plot_xyz(h1, xe1, xe2, xe3, orbit->last_indix+1, (char*)"", "points", "1", "2", color);

    gnuplot_set_xlabel (h2, (char*) "Xh");
    gnuplot_set_ylabel (h2, (char*) "Yh");
    gnuplot_set_zlabel (h2, (char*) "Zh");
    gnuplot_cmd(h2,  "set title \"Tmap (TF) \"");
    gnuplot_plot_xyz(h2, xs1, xs2, xs3, orbit->last_indix+1, (char*)"", "points", "1", "2", color);
}



//-----------------------------------------------------------------------------------------------------------------------------------------------------
//
//  Square distance & minimization of it (deprecated)
//
//-----------------------------------------------------------------------------------------------------------------------------------------------------
/**
    \brief Given an initial guess st0, computes the min argument of square distance between a given configuraton z1 and z = W(g(st), t) (see sqdist).
    i.e. st1 = argmin sqdist(st, z1, t, &orbit)

    \param Orbit a reference to the current orbit
    \param z1 the target NC configuration (dim 6)
    \param t the current time
    \param st0 the initial guess in real TFC configuration (dim 4)
    \param st1 the min argument to update real TFC configuration (dim 4)
    \param multimin_method an integer to select a multidimensionnal minimization method
**/
void argmin_sqdist(Orbit &orbit, double z1[], double t, double st0[], double st1[], int multimin_method)
{
    double ftol   = 1e-6;  //Common precision for the multimin routines
    int niter = 0;         //Number of iterations

    //Switch on the multidimensionnal minimization method to use
    switch(multimin_method)
    {

        //Simplex method
    case MULTIMIN_SIMPLEX:
    {
        double **p    = dmatrix(0, 4, 0, 3); //5*4 matrices to store the vertices of the initial search simplex
        double *y     = dvector(0, 4);       //5*1 vector so that y[i] = y[i] = sqdist(p[i])
        double **id   = dmatrix(0, 3, 0, 3); //4*4 identity matrix
        double lambda = 0.1;                 //size of the initial simplex (arbitrary, might be neater to set it higher in the code)

        //Compute the identity matrix (ej)
        for(int i= 0; i < 4; i++) for(int j = 0; j < 4; j++) id[i][j] = (i==j)? 1:0;
        //First vertix is the initial guess st0
        for(int j = 0; j < 4; j++) p[0][j] = st0[j];
        //4 others are st0+lambda*ej
        for(int i= 1; i <= 4; i++)  for(int j = 0; j < 4; j++) p[i][j] = st0[j] + lambda*id[i-1][j];
        //Evaluate the 4+1 vertices and store them in y
        for(int i = 0; i <= 4; i++) y[i] = sqdist(p[i], z1, t, &orbit);

        //Simplex method
        tic();
        amoeba(z1, t, &orbit, p, y, 4, ftol, sqdist, &niter);
        cout << "Simplex method. niter = " << niter << "in " << toc() << "s. " << endl;

        //Update st1
        for(int i = 0; i<4; i++) st1[i] = p[0][i]; //arbitrary, we can choose p[k], for k= 0...4

        //Memory release
        free_dmatrix(p,  0, 4, 0, 3);
        free_dmatrix(id, 0, 3, 0, 3);
        free_dvector(y, 0, 4);
        break;
    }



    //Fletcher-Reeves-Polak-Ribiere minimization
    case MULTIMIN_FRPRMN ... MULTIMIN_DFRPRMN:
    {

        double *pv    = dvector(0, 3);  //4*1 vector that contains the final argmin
        double fret;                    //final minimum found

        //Initialize pv to st0
        for(int j = 0; j < 4; j++) pv[j] = st0[j];

        //frprmn method
        tic();
        if(multimin_method == MULTIMIN_FRPRMN)
        {
            frprmn(z1, t, &orbit, pv, 3, ftol, &niter, &fret, sqdist, dsqdist);
            //cout << "frprmn method. niter = " << niter << " in " << toc() << "s. " << endl;
        }
        else
        {
            dfrprmn(z1, t, &orbit, pv, 3, ftol, &niter, &fret, sqdist, dsqdist);
            //cout << "dfrprmn method. niter = " << niter << "in " << toc() << "s. " << endl;
        }

        //Update st1
        for(int i = 0; i<4; i++) st1[i] = pv[i];

        //Free memory
        free_dvector(pv, 0, 3);
        break;
    }

    }

}

/**
    \brief Computes the square distance between a given configuraton z1 and z0 = W(g(st0), t)

    \param st0 an array of 4 double which gives the configuration to study in real TFC coordinates
    \param z1 the target NC configuration (dim 6)
    \param t the current time
    \param params a set of parameters needed to perform evaluations of W and g
    \return a double, the square distance between a given configuraton z1 and z0 = W(g(st0), t)
**/
double sqdist(double st0[], double z1[], double t, void *params)
{
    //Pointers
    Orbit *orbit = (Orbit*) params;

    //------------------------------------------
    // Inner variables (NC, TFC)
    //------------------------------------------
    cdouble s0[4];
    cdouble z0[6];
    double z0d[6];
//    double ys[6];

    //------------------------------------------
    //"Realification": s0 = REAL(st0)
    //------------------------------------------
    s0[0] = 1.0/sqrt(2)*(st0[0] - st0[2]*I);
    s0[2] = 1.0/sqrt(2)*(st0[2] - st0[0]*I);
    s0[1] = 1.0/sqrt(2)*(st0[1] - st0[3]*I);
    s0[3] = 1.0/sqrt(2)*(st0[3] - st0[1]*I);

    //------------------------------------------
    // 2. Update z0
    //------------------------------------------
    double res = 0;
    for(int p = 0; p < 6; p++)
    {
        orbit->W->at(p).evaluate(s0, *orbit->ofs, orbit->pmap->order, orbit->pmap->ofs_order);
        z0[p] = orbit->ofs->evaluate(orbit->qbcp_l->us.n*t);
        z0d[p] = creal(z0[p]);
        res += (z0d[p] - z1[p])*(z0d[p] - z1[p]);
        //Add the delta in energy?
        //NCtoSYS(t, z0d, ys, orbit->ode_s_6->d->sys->params);
        //res += (qbfbp_H(t, ys, orbit->ode_s_6->d->sys->params) - orbit->pmap->Hvalue)*(qbfbp_H(t, ys, orbit->ode_s_6->d->sys->params) - orbit->pmap->Hvalue);
    }

    return res;
}

/**
    \brief Computes the gradient square distance between a given configuraton z1 and z0 = W(g(st0), t) (see sqdist)

    \param st0 an array of 4 double which gives the configuration to study in real TFC coordinates
    \param df the gradient to update
    \param z1 the target NC configuration (dim 6)
    \param t the current time
    \param params a set of parameters needed to perform evaluations of W and g
**/
void dsqdist(double st0[], double df[], double z1[], double t, void *params)
{
    //Pointers
    Orbit *orbit = (Orbit*) params;

    //------------------------------------------
    // Inner variables (NC, TFC)
    //------------------------------------------
    cdouble s0[4];
    cdouble z0;
    double w[6];
    cdouble **dw    = dcmatrix(0,5,0,3);
    cdouble **dg    = dcmatrix(0,3,0,3);
    cdouble **dwdg  = dcmatrix(0,5,0,3);

    //------------------------------------------
    //"Realification": s0 = REAL(st0)
    //------------------------------------------
    s0[0] = 1.0/sqrt(2)*(st0[0] - st0[2]*I);
    s0[2] = 1.0/sqrt(2)*(st0[2] - st0[0]*I);
    s0[1] = 1.0/sqrt(2)*(st0[1] - st0[3]*I);
    s0[3] = 1.0/sqrt(2)*(st0[3] - st0[1]*I);

    //------------------------------------------
    // Evaluating W and the Jacobian of W
    //------------------------------------------
    for(int i = 0; i < 6; i++)
    {
        //AUX = Wi(s0)
        orbit->W->at(i).evaluate(s0, *orbit->ofs, orbit->pmap->order, orbit->pmap->ofs_order);
        z0 = orbit->ofs->evaluate(orbit->qbcp_l->us.n*t);
        //w[i] = Aux(nt)
        w[i] = creal(z0);

        for(int j = 0; j < 4; j++)
        {
            //AUX= DW_ij(s0)
            orbit->DW->getCA(i,j)->evaluate(s0, *orbit->ofs, orbit->pmap->order, orbit->pmap->ofs_order);
            //dw[i][j] = AUX(nt)
            dw[i][j] = orbit->ofs->evaluate(orbit->qbcp_l->us.n*t);
        }
    }

    //------------------------------------------
    // Building the matrix dg
    //------------------------------------------
    dg[0][0] = 1.0/sqrt(2)+0.0*I;
    dg[0][1] = 0.0+0.0*I;
    dg[0][2] = -1.0/sqrt(2)*I;
    dg[0][3] = 0.0+0.0*I;

    dg[1][0] = 0.0+0.0*I;
    dg[1][1] = 1.0/sqrt(2)+0.0*I;
    dg[1][2] = 0.0+0.0*I;
    dg[1][3] = -1.0/sqrt(2)*I;

    dg[2][0] = -1.0/sqrt(2)*I;
    dg[2][1] = 0.0+0.0*I;
    dg[2][2] = 1.0/sqrt(2)+0.0*I;
    dg[2][3] = 0.0+0.0*I;

    dg[3][0] = 0.0+0.0*I;
    dg[3][1] = -1.0/sqrt(2)*I;
    dg[3][2] = 0.0+0.0*I;
    dg[3][3] = 1.0/sqrt(2)+0.0*I;

    //------------------------------------------
    // Evaluating the product dw*dg
    //------------------------------------------
    for(int i = 0; i < 6; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            dwdg[i][j] = 0.0+0.0*I;
            for(int k = 0; k <4; k++) dwdg[i][j] += dw[i][k]*dg[k][j];
        }
    }

    //------------------------------------------
    // Evaluating the gradient
    //------------------------------------------
    for(int p = 0; p < 4; p++)
    {
        df[p] = 0.0;
        for(int k = 0; k < 6; k++)
        {
            df[p] += creal(2*(w[k] - z1[k])*dwdg[k][p]);
            //cout << "Real? " << creal(2*(w[k] - z1[k])*dwdg[k][p]) << " + " << cimag(2*(w[k] - z1[k])*dwdg[k][p]) << endl;
        }
    }


    //------------------------------------------
    // Free memory
    //------------------------------------------
    free_dcmatrix(dg, 0, 3, 0, 3);
    free_dcmatrix(dw, 0, 5, 0, 3);
    free_dcmatrix(dwdg, 0, 5, 0, 3);
}


//-----------------------------------------------------------------------------------------------------------------------------------------------------
//
//  Orbit C structure handling
//
//-----------------------------------------------------------------------------------------------------------------------------------------------------
/**
    \brief Initialize one orbit structure
 **/
void init_orbit(Orbit *orbit,
                vector<Oftsc>*  W,
                vector<Oftsc>*  Wh,
                matrix<Oftsc>* DW,
                matrix<Ofsc>*  PC,
                vector<Ofsc>*  V,
                value_params   *val_par,
                value_function *fvalue,
                OdeStruct *ode_s_6,
                OdeStruct *ode_s_8,
                OdeStruct *ode_s_6_root,
                OdeStruct *ode_s_8_root,
                Pmap *pmap,
                QBCP_L *qbcp_l,
                Ofsc* orbit_ofs,
                int vdim,
                int label)
{
    //-----------
    //Parent
    //-----------
    orbit->pmap    = pmap;     //Poincare map (parent)
    orbit->qbcp_l = qbcp_l;  //QBCP around a given Li point (parent)

    //-----------
    //Parameterization (common to all orbits)
    //-----------
    orbit->W    =  W;            //z(t) = W(s(t), t)
    orbit->Wh   =  Wh;           //zh(t) = Wh(s(t), t)
    orbit->DW   =  DW;           //Jacobian of W
    orbit->ofs  =  orbit_ofs;    //Auxiliary Ofs object

    //-----------
    //COC (common to all orbits)
    //-----------
    orbit->PC  = PC;            //COC matrix
    orbit->V   = V;             //COC vector

    //-----------
    //Pulsation
    //-----------
    orbit->n = qbcp_l->us.n;


    //-----------
    //For event detection
    //-----------
    orbit->val_par    = val_par;           //Event parameters
    orbit->fvalue     = fvalue;            //fvalue for event detection

    orbit->z0_mat  = dmatrix(0, 5, 0, val_par->max_events); //Pointer towards the stored position of events
    orbit->zh0_mat = dmatrix(0, 5, 0, val_par->max_events); //Pointer towards the stored position of events
    orbit->s0_mat  = dmatrix(0, REDUCED_NV-1, 0, val_par->max_events); //Pointer towards the stored position of events
    orbit->te_mat  = dvector(0, val_par->max_events);       //Pointer towards the stored time of events
    orbit->hz      = dvector(0, val_par->max_events);       //Energy gap between z(t) and the pm W(s(t),t) at each event
    orbit->hw      = dvector(0, val_par->max_events);       //Energy gap between W(s(t),t) and the initial energy at each event
    orbit->ePm     = dvector(0, val_par->max_events);       //Projection error at each event
    orbit->nevent  = ivector(0, val_par->max_events);       //Number of the event

    //-----------
    //Characteristics
    //-----------
    orbit->z0  = dvector(0, 5);                 //Initial position in NC coordinates dim = 6
    orbit->zh0 = dvector(0, 5);                 //Initial position in TFC coordinates dim = 6
    orbit->si  = dvector(0, REDUCED_NV-1);      //Initial RCM configuration dim = 4
    orbit->s0d = dvector(0, 2*REDUCED_NV-1);    //Initial position in CCM8 coordinates (real+imag part) dim = 8
    orbit->xf  = dvector(0, 5);                 //Final position NC dim = 6
    orbit->s0  = dcvector(0,REDUCED_NV-1);      //Initial position in CCM4 coordinates (real+imag part) dim = 4
    orbit->tf  = pmap->tf;                      //Final time after computation
    orbit->int_method  = -1;                    //Integration method. -1 if not computed
    orbit->label = label;                       //Label of the orbit in the map
    orbit->vdim = vdim;                         //Modified dimension during energy adjustment

    //-----------
    //ODE integration
    //-----------
    orbit->ode_s_6    = ode_s_6;              //NC ode struct
    orbit->ode_s_8    = ode_s_8;              //TFC ode struct
    orbit->ode_s_6_root    = ode_s_6_root;    //NC ode struct (root finding)
    orbit->ode_s_8_root    = ode_s_8_root;    //TFC ode struct (root finding)
}



/**
    \brief Free one orbit
 **/
void free_orbit(Orbit *orbit)
{
    //-----------
    //Characteristics
    //-----------
    free_dvector(orbit->z0,  0, 5);
    free_dvector(orbit->zh0, 0, 5);
    free_dvector(orbit->si,  0, REDUCED_NV-1);
    free_dvector(orbit->s0d, 0, 2*REDUCED_NV-1);
    free_dvector(orbit->xf,  0, 5);
    free_dcvector(orbit->s0, 0, REDUCED_NV-1);

    //-----------
    //For event detection
    //-----------
    free_dmatrix(orbit->z0_mat,  0, 5, 0, orbit->val_par->max_events);
    free_dmatrix(orbit->zh0_mat, 0, 5, 0, orbit->val_par->max_events);
    free_dmatrix(orbit->s0_mat,  0, REDUCED_NV-1, 0, orbit->val_par->max_events);

    free_dvector(orbit->te_mat,  0, orbit->val_par->max_events);
    free_dvector(orbit->ePm,     0, orbit->val_par->max_events);
    free_dvector(orbit->hz,      0, orbit->val_par->max_events);
    free_dvector(orbit->hw,      0, orbit->val_par->max_events);
    free_ivector(orbit->nevent,  0, orbit->val_par->max_events);

    free_ode_structure(orbit->ode_s_6);
    free_ode_structure(orbit->ode_s_6_root);
    free_ode_structure(orbit->ode_s_8);
    free_ode_structure(orbit->ode_s_8_root);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------
//Backup
//-----------------------------------------------------------------------------------------------------------------------------------------------------
int gslc_pmin_step(Orbit *orbit,
                   OdeStruct *ode_s,
                   double yv[],
                   double s1[],
                   int *reset_number,
                   double *t,
                   double t1,
                   double threshold)
{
    int status;
    double eOm;
    //----------------------
    //Evolve one step of z(t)
    //----------------------
    status = gsl_odeiv2_evolve_apply (ode_s->e, ode_s->c, ode_s->s, &ode_s->sys, t, t1, &ode_s->h, yv);
    if (status != GSL_SUCCESS)
    {
        cout << "error in gslc_dual_step: integration of z(t) has gone wrong. break." << endl;
        return GSL_FAILURE;
    }

    tic();
    argmin_sqdist(*orbit, yv, *t, s1, s1, MULTIMIN_FRPRMN);
    eOm = sqdist(s1, yv, *t, orbit);
    cout << "Minimization done in " << toc() << "s." << endl;

    cout << "gslc_pmin_step. t  = " << *t << ", z = " << yv[2] << ", pmin = " << eOm << endl;
    cout << "Hamiltonian = " << orbit_ham(*orbit->pmap, s1) << endl;

    //----------------------
    //Reset if necessary
    //----------------------
    if(eOm > threshold)
    {
        cout << "Reset is detected @t  = " << *t << ", z = " << yv[2] << ", pmin = " << eOm << endl;
        RCMtoNC(s1, *t, orbit->n, orbit->pmap->order, orbit->pmap->ofs_order, *orbit->W, *orbit->ofs, yv, true);
        (*reset_number)++;
    }

    return status;
}



