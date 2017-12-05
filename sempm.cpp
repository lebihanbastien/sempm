/**
 * \file main.cpp
 * \brief Main file (exec) for the computation of the dynamics about the neighborhood of
 *        the points EML1,2 and SEL1,2
 *
 * The following models of the Sun-Earth-Moon system are available:
 *   - Quasi-Bicircular Four-Body Problem (QBCP).
 *   - Coupled Circular Restricted Three-Body Problem (CRTBP) (EM and SE subsystems).
 * The following models of the Sun-Earth-Moon system may be available in the future:
 *   - Bircicular Four-Body Problem (BCP) (work in progress. the real issue being: no
 *     clear dynamical equivalent for EML2
 *   - Coupled Elliptic Restricted Three-Body Problem (CRTBP) (far from being completed).
 * \author BLB.
 */

//----------------------------------------------------------------------------------------
// Include
//----------------------------------------------------------------------------------------
//Parallel computing
#include <omp.h>
//Custom
#include "config.h"
#include "ofs.h"
#include "ofts.h"
#include "poincare.h"
#include "pmt.h"
#include "pmode.h"
//Tests
#include "ofs_test.h"
#include "oftsh_test.h"
#include "ofts_test.h"
//C routines and files
extern "C" {
    #include "gnuplot_i.h"
    void fortfunc_(int *ii, float *ff);
    void saxpy_(int *n, double *alpha, double *x, double *y);
    void matvec_(int *n, double **A, double *x, double *y, int *lda);
}

//----------------------------------------------------------------------------------------
// Namespaces
//----------------------------------------------------------------------------------------
using namespace std;


//----------------------------------------------------------------------------------------
// Main
//----------------------------------------------------------------------------------------
/**
 *  \fn int main()
 *  \brief Main routine. It makes of several arguments via argv (see code for details).
 */
int main(int argc, char *argv[])
{
    cout << "---------------------------------------------------" << endl;
    cout << "                                                   " << endl;
    cout << "                    SEMPM                          " << endl;
    cout << "                                                   " << endl;
    cout << "---------------------------------------------------" << endl;
    cout << argc << " arguments have been passed to sempm"        << endl;

    //====================================================================================
    // Initialization of Manip object, which allows for algebraic manipulation
    // of polynomials (exponents in reverse lexicographic order, number of coefficient in
    // homogeneous polynomials...).
    // Note that, for now, the default order of the OFS and OFTS objects are set at
    // compilation (OFS_ORDER_0OFTS_ORDER_0).
    //====================================================================================
    int OFS_ORDER_0  = 30;  //The default OFS_ORDER (Fourier series) is 30
    int OFTS_ORDER_0 = 30;  //The default OFTS_ORDER (Taylor series) is 30
    int nr = max(2*OFS_ORDER_0, OFTS_ORDER_0);
    int nv = Csts::NV;
    Manip::init(nv, nr);


    //====================================================================================
    // Get the user-defined arguments
    //====================================================================================
    // The variable index contains the index of the current argument argv[index]
    // of the main routine. The first index is 1 because index 0 is the application path.
    int index    = 1;

    //------------------------------------------------------------------------------------
    // Get the current orders
    //------------------------------------------------------------------------------------
    OFTS_ORDER   = atoi(argv[index++]); //Order of the Taylor series
    OFS_ORDER    = atoi(argv[index++]); //Order of the Fourier series
    REDUCED_NV   = atoi(argv[index++]); //Number of reduced variables (4, 5 or 6).

    //To check consistency with the desired values
    cout << "Current orders:" << endl;
    cout << "OFTS_ORDER = "   << OFTS_ORDER << endl;
    cout << "OFS_ORDER  = "   << OFS_ORDER  << endl;
    cout << "REDUCED_NV = "   << REDUCED_NV << endl;
    cout << "---------------------------------------------------" << endl;

    //------------------------------------------------------------------------------------
    // Retrieving the parameters, in this order:
    // 1. comp_type: Type of computation (QBTBP, NFO2, PM...)
    // 2. model: Type of model (QBCP, CRTBP)
    // 3. coord_sys: Default coordinate system (EM, SEM)
    // 4. dli: Default libration point
    // 5. param_style: Parameterization (PM) style (only used in some computations)
    // 6. man_type: Type of manifold (center, center-stable, center-unstable)
    // 7. storage: Boolean for storage (in txt/bin files) of the results
    //------------------------------------------------------------------------------------
    int comp_type   = atoi(argv[index++]);
    int model       = atoi(argv[index++]);
    int dli         = atoi(argv[index++]);
    int param_style = atoi(argv[index++]);
    int man_type    = atoi(argv[index++]);
    int storage     = atoi(argv[index++]);

    //Check
    cout << "Current parameters: " << endl;
    cout << "--------------------" << endl;
    cout << "comp_type   = "   << comp_type << endl;
    cout << "model       = "   << model << endl;
    cout << "li          = "   << dli << endl;
    cout << "param_style = "   << param_style << endl;
    cout << "man_type    = "   << man_type << endl;
    cout << "storage     = "   << storage << endl;
    cout << "---------------------------------------------------" << endl;

    //------------------------------------------------------------------------------------
    // Set the global variable MODEL_TYPE
    //------------------------------------------------------------------------------------
    MODEL_TYPE = model;

    //====================================================================================
    // Initialization of the environnement (the Four-Body Problem at hand).
    // Mandatory to perform any computation, except qbtbp(int)
    //====================================================================================
    init_env(model, dli, man_type, param_style);


    //====================================================================================
    // Master switch on the type of computation required by the user
    //====================================================================================
    switch(comp_type)
    {

    //------------------------------------------------------------------------------------
    // Unitary tests of routines for the Fourier-Taylor algebra
    //
    // In bash, launch with e.g. sempm.sh config/test_fourier_taylor.sh
    //------------------------------------------------------------------------------------
    case Csts::FT_TEST:
    {
        //Test of the Ofs class (Fourier series)
        ofs_test();
        //Test of the Oftsh class (Fourier-Taylor homogeneous polynomials)
        oftsh_test();
        //Test of the Ofts class (Fourier-Taylor series)
        ofts_test();
        break;
    }


    //------------------------------------------------------------------------------------
    // Sun-Earth-Moon three-body motion resolution up to OFS_ORDER.
    //
    // In bash, launch with e.g. sempm.sh config/qbcp.sh
    //------------------------------------------------------------------------------------
    case Csts::QBTBP:
    {
        switch(model)
        {
            //----------------------------------------------------------------------------
            // QBTBP & QBCP resolution up to OFS_ORDER.
            // If storage==true, results are stored in the folders
            // data/qbtbp and data/VF/QBCP.
            //----------------------------------------------------------------------------
        case Csts::QBCP:
            qbtbp_and_qbcp(true, storage);
            break;

            //----------------------------------------------------------------------------
            // BCP resolution in OFS format (in order to match QBCP format).
            // Results are stored in the folder data/VF/BCP.
            //----------------------------------------------------------------------------
        case Csts::BCP:
            bcp();
            break;
        }
        break;
    }

    //------------------------------------------------------------------------------------
    // Compute the complete change of coordinates to:
    // - Get rid of order 1,
    // - Get a normal form for the order 2 of the Hamiltonian of the QBCP.
    // If storage==true, results are stored in the folder data/COC
    //
    // In bash, launch with e.g. sempm.sh config/nfo2/nfo2_eml2.sh
    //------------------------------------------------------------------------------------
    case Csts::NFO2:
    {
        switch(model)
        {
        case Csts::QBCP:
            nfo2(SEML, storage);
            break;
        default:
            cout << "Error: nfo2 routine only implemented for the QBCP." << endl;
            break;
        }
        break;
    }

    //------------------------------------------------------------------------------------
    // Parameterization method up to OFTS_ORDER.
    //
    // In bash, launch with e.g. sempm.sh config/pm/pm_cm_eml2.sh
    //------------------------------------------------------------------------------------
    case Csts::PM:
    {
        //--------------------------------------------------------------------------------
        //threshold for small divisors = 1e-2 (could be propagated to the user)
        //--------------------------------------------------------------------------------
        double small_div_threshold = 1e-2;

        //--------------------------------------------------------------------------------
        //param method
        //--------------------------------------------------------------------------------
        pmt(SEML.cs.man_type, param_style, small_div_threshold, storage);
        break;
    }

    //------------------------------------------------------------------------------------
    // Test of the parameterization method
    //
    // In bash, launch with e.g. sempm.sh config/test_pm/test_cm_eml2.sh
    //------------------------------------------------------------------------------------
    case Csts::PM_TEST:
    {
        //--------------------------------------------------------------------------------
        // Additionnal parameters in bash
        //--------------------------------------------------------------------------------
        // 1. Initial conditions
        double si[REDUCED_NV];
        for(int i = 0; i < REDUCED_NV; i++) si[i] = atof(argv[index++]);

        //2. Array of ofts_orders to test
        int n_ofts_order = atoi(argv[index++]);
        int v_ofts_order[n_ofts_order];
        for(int i = 0; i<n_ofts_order; i++)
        {
            v_ofts_order[i] = atoi(argv[index++]);
            if(v_ofts_order[i] > OFTS_ORDER)
            {
                cout << "Warning: a required order is greater than OFTS_ORDER. " << endl;
                cout << "The computation will be cut at OFTS_ORDER." << endl;
            }
        }

        //3. Array of ofs_orders to test
        int n_ofs_order = atoi(argv[index++]);
        int v_ofs_order[n_ofs_order];
        for(int i = 0; i<n_ofs_order; i++)
        {
            v_ofs_order[i] = atoi(argv[index++]);
            if(v_ofs_order[i] > OFS_ORDER)
            {
                cout << "Warning: a required order is greater than OFS_ORDER. " << endl;
                cout << "The computation will be cut at OFS_ORDER." << endl;
            }
        }

        //4. Time interval
        double tvec[2];
        tvec[0] = atof(argv[index++])*SEML.us.T;
        tvec[1] = atof(argv[index++])*SEML.us.T;


        //--------------------------------------------------------------------------------
        // Initialization of the invariant manifold
        //--------------------------------------------------------------------------------
        init_inv_man(SEML);
        init_coc(SEML);

        //--------------------------------------------------------------------------------
        // Test
        //--------------------------------------------------------------------------------
        pm_error_vs_orders_test(n_ofts_order, v_ofts_order, n_ofs_order, v_ofs_order, si, tvec);

        break;
    }

    //------------------------------------------------------------------------------------
    // Poincare maps (minimum working example for quick testing)
    //
    // In bash, launch with e.g. sempm.sh config/pmap/pmap_example_qbcp.sh
    //------------------------------------------------------------------------------------
    case Csts::COMPMAP:
    {
        //--------------------------------------------------------------------------------
        // Initialization of the invariant manifold
        //--------------------------------------------------------------------------------
        init_inv_man(SEML);
        init_coc(SEML);


        //--------------------------------------------------------------------------------
        // Pmap init
        //--------------------------------------------------------------------------------
        Pmap pmap;
        if(MODEL_TYPE == Csts::CRTBP)
        {
            //Quite stable parameters (hard coded here)
            //------------------------
            pmap.pabs           =  1e-14;       //absolute precision in numerical integrator
            pmap.prel           =  1e-14;       //relative precision in numerical integrator
            pmap.proot          =  1e-10;       //root precision in root-finding routines
            pmap.eproj          =  1e-6;        //error limit above which a projection of the current state to the semi-analytical parameterization is required
            pmap.var_dim_h      =  2;           //the dimension of the reduced vector that varies in order to guarantee a certain initial energy value (see orbit_update_s0)
            pmap.max_rad_div    =  5.0;         //maximum radius about the origin (in NC units) beyond which the trajectory is supposed to be divergent
            pmap.T              =  SEML.us.T;;  //SEM period
        }
        else
        {
            //Quite stable parameters (hard coded here)
            //------------------------
            pmap.pabs        =  1e-14;
            pmap.prel        =  1e-14;
            pmap.proot       =  1e-10;
            pmap.eproj       =  1e-6;
            pmap.var_dim_h   =  2;
            pmap.max_rad_div =  5.0;
            pmap.T           =  SEML.us.T;
        }


        //--------------------------------------------------------------------------------
        // Parameters that may evolve
        //--------------------------------------------------------------------------------
        pmap.type       =  atoi(argv[index++]);  //the type of map: only PMAP is available in this version of the code
        pmap.ofts_order =  atoi(argv[index++]);  //the effective order of the Taylor expansions on the map.
        pmap.ofs_order  =  atoi(argv[index++]);  //the effective order of the Fourier expansions on the map.
        pmap.dh_nc_t0   =  atof(argv[index++]);  //the initial relative energy at t = t0, in NC coordinates
        pmap.f_proj_T   =  atof(argv[index++]);         //the frequency of projection inside the map (given as a fraction of the SEM period T)
        pmap.t_proj_nc  =  1.0/pmap.f_proj_T*SEML.us.T; //the corresponding time, in normalized time units
        pmap.t0_nc      =  atof(argv[index++])*SEML.us.T;  //the initial time for all trajectories on the map, in normalized time units
        pmap.tmax_nc    =  atof(argv[index++])*SEML.us.T;  //the maximum time for all trajectories on the map, in normalized time units (usually set to a very large value)
        pmap.max_events =  atoi(argv[index++]);  //the maximum number of events on the map, i.e. the number of crossings of the xy-plane, with dot(z) > 0
        pmap.n_sol      =  atoi(argv[index++]);  //number of computed trajectories
        pmap.si_min     =  atof(argv[index++]);  //minimum value for the initial conditions (s1, s3) on the map
        pmap.si_max     =  atof(argv[index++]);  //maximum value for the initial conditions (s1, s3) on the map

        //boolean to true of the graph style is used (change the evalutation of the parameterization)
        if(SEML.param_style == Csts::GRAPH) pmap.graph_style = 1;
        else pmap.graph_style = 0;


        pmap.is_plot    = atoi(argv[index++]);      //plot boolean
        pmap.is_par     = atoi(argv[index++]);      //parallel computation of not
        pmap.int_method = atoi(argv[index++]);      //method of integration


        cout << "pmap.type        = "  << pmap.type << endl;
        cout << "pmap.ofts_order  = "  << pmap.ofts_order << endl;
        cout << "pmap.ofs_order   = "  << pmap.ofs_order  << endl;
        cout << "pmap.dh_nc_t0    = "  << pmap.dh_nc_t0   << endl;
        cout << "pmap.f_proj_T    = "  << pmap.f_proj_T << endl;


        cout << "pmap.t0_nc       = "  << pmap.t0_nc  << endl;
        cout << "pmap.tmax_nc     = "  << pmap.tmax_nc << endl;

        cout << "pmap.n_sol       = "  << pmap.n_sol << endl;
        cout << "pmap.si_min      = "  << pmap.si_min << endl;
        cout << "pmap.si_max      = "  << pmap.si_max << endl;

        cout << "pmap.graph_style = "  << pmap.graph_style << endl;
        cout << "pmap.max_events  = "  << pmap.max_events << endl;

        cout << "pmap.is_plot     = "  << pmap.is_plot << endl;
        cout << "pmap.is_par      = "  << pmap.is_par << endl;
        cout << "pmap.int_method  = "  << pmap.int_method << endl;
        cout << "---------------------------------------------------" << endl;

        //--------------------------------------------------------------------------------
        // openMP settings
        //--------------------------------------------------------------------------------
        int num_threads = atoi(argv[index++]);
        omp_set_num_threads(num_threads);
        cout << "num_threads in openMP   = "  << num_threads << endl;
        cout << "---------------------------------------------------" << endl;


        //--------------------------------------------------------------------------------
        // Pmap computation
        //--------------------------------------------------------------------------------
        pmap_build_random(pmap);

        break; //and of case 4: pmaps
    }

    //------------------------------------------------------------------------------------
    // Compute the dynamical equivalent of the libration point
    // (minimum working example for quick testing)
    //
    // In bash, launch with e.g. sempm.sh config/dyneq/dyneq_eml2.sh
    //------------------------------------------------------------------------------------
    case Csts::DYNEQ:
    {
        //--------------------------------------------------------------------------------
        // Direction computation of the lpdyneq_single_shooting. Works for:
        // - EML1, EML2 in QBCP
        // - SEML1, SEML2 in QBCP
        //--------------------------------------------------------------------------------
        compute_dyn_eq_lib_point(SEML, storage);
        break;
    }


    }



    return 0;
}



