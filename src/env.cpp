/**
 * \file env.cpp
 * \brief Routines for the definition of the working environment via several C structure
 *        that represent Three or Four-Body models of the Sun-Earth-Moon system.
 *        The main routine is init_fbp_lib, which initializes a FBPL structure, i.e. a
 *        structure that contains all the necessary constants and coefficients
 *        (equations of motion, Hamiltonian) for a given Four-Body problem
 *        around a given libration point. The coefficients are retrieved from txt files
 *        stored in subfolders of ./data/
 *        that must have been computed with the routine qbtbp().
 *
 *        The physical constants hard-coded in this file have been taken from:
 *          - Andreu 1998 for the Sun-Earth-Moon system,
 *          - JPL/Goddard Space Flight Center websites for the rest of the solar system.
 *
 *        Note that most of the routines initialize C structures with C routines,
 *        because they were taken from original C code made in 2014.
 *        Most of these C structures could be implemented again as C++ classes with
 *        proper constructors. This would requires quite a bit of work because
 *        most of these structures are used on the fly in the current code.
 *        This file could also be merged with init.cpp in the long run, in order
 *        to have one single file for the initialization routines.
 *
 * \author BLB
 */

#include "env.h"

//----------------------------------------------------------------------------------------
//            Initialization routines
//----------------------------------------------------------------------------------------
/**
 * \brief Initialize a FBPL structure, i.e. a FBP (Four-Body Problem) focused on two
 *        libration points: one the EM system and one of the SEM system.
 *        As for now, the libration point must be L1 or L2 for both systems.
 *        The FBPL structure fbpl contains several variables, including:
 *
 *          - fbpl.li, the current libration point at hand (1 or 2).
 *
 *          - fbpl.us_em, a USYS (unit system) structure. It contains several
 *          constants (distances, frequencies) in Earth-Moon  normalized units.
 *          - fbpl.us_sem, a USYS (unit system) structure. It contains several
 *          constants (distances, frequencies) in Sun-Earth normalized units.
 *
 *          - fbpl.cs_em_l1, a CSYS (coordinate system) structure. It contains several
 *          constants in Earth-Moon normalized units that are associated with the EML1
 *          point (gamma, c1), as well as the coefficients of the vector field of
 *          the Four-Body Problem, both in EM and EMNC coordinates.
 *          - fbpl.cs_em_l2, cs_sem_l1, cs_sem_l2, equivalents structures
 *          for the other points.
 *
 *          Moreover, the USYS and CSYS structures are duplicated for easy access:
 *          - fbpl.us, equal one of the previous USYS structures:
 *                  * us_em if coord_sys == Csts::EM,
 *                  * us_Sem if coord_sys == Csts::SEM.
 *          - fbpl.cs_em, equal to cs_em_l1,2 if li_EM == 1,2.
 *          - fbpl.cs_sem, equal to cs_sem_l1,2 if li_SEM == 1,2.
 *          - fbpl.cs, equal to
 *                  * cs_em_l1,2 if li_EM == 1,2 & coord_sys == Csts::EM,
 *                  * cs_sem_l1,2 if li_SEM == 1,2 & coord_sys == Csts::SEM.
 *
 *        In the end, we mainly use fbpl.li (libration point), fbpl.us (constants in the
 *        suitable normalized units), and fbpl.cs (constants and coefficients about the
 *        proper libration point).
 *        See the definition of the FBPL structure to see the other variables contained
 *        in it.
 *
 * Note that the FBPL is focused on one libration point but contains the constants and
 * coefficients for all four points (EML1, EML2, SEL1, SEL2).
 * The focus of the FBPL structure can be changed from one point to another via
 * the routines change_coord and change_li_coord. Moreover, the coefficients for each
 * point can be accessed via calls of the form fbpl->cs_em_l1, fbpl->cs_em_l2, etc.
 *
 * The inputs of this routines are:
 * \param fbpl a pointer on the FBPL structure to initialize.
 * \param fbp a pointer on the FBP parent structure.
 * \param li_EM number of the libration point for the EM system.
 * \param li_SEM number of the libration point for the SEM system.
 * \param model: either Csts::QBCP, Csts::BCP, Csts::CRTBP, etc.
 * \param coord_sys: default coordinate system for this structure:
 *           - if coord_sys == Csts::EM,  the fbpl is focused on the li_EM  point of the
 *             EM system.
 *           - if coord_sys == Csts::SEM, the fbpl is focused on the li_SEM point of the
 *             SEM system.
 *        The focus can be change dynamically during computation, via the routines
 *        change_coord and change_li_coord.
 * \param param_style: type of parameterization of the manifolds (Csts::GRAPH, etc). Note that
 *        the param_style influences the number of coefficients taken into account in the
 *        Fourier series. Indeed, for graph method, the reduced vector field is non
 *        autonomous, and full Fourier series are used. For normal form, the reduced
 *        vector field is quasi autonomous and we can safely reduce the order of the
 *        series to 5 (11 coefficients taken into account).
 * \param man_type_EM: type of manifold about li_EM (Csts::MAN_CENTER, etc).
 * \param man_type_SEM: type of manifold about li_SEM.
 * \param is_new an integer: equal to 1 if no solution has been previously
 *        computed with the routine qbtbp(), 0 otherwise.
 * \param is_norm: are the equations of motion normalized?  Probably deprecated,
 *        should be always true. Kept for consistency with older code.
 *
 *
 * Note that the FBP structure fbp is used only for the initialization of the coordinate
 * systems. More precisely, it contains some parameters specific to each libration point
 * (gamma, c1, etc), via its CR3BP structures
 * (see init_fbp, the routine that initializes the QBCP structures).
 *
 * Note also that there is no speficic need to use two libration points at the
 * same time (one of the EM and one of the SEM systems). This choice was made with the
 * study of li_EM-li_SEM connections in mind, which are now computed in another package.
 *
 **/
void init_fbp_lib(FBPL *fbpl, FBP *fbp, int li_EM, int li_SEM,
               int model, int coord_sys, int param_style, int man_type_EM, int man_type_SEM,
               int is_new, int is_norm)
{
    //====================================================================================
    //      Common to all models
    //      These settings may be needed to initialize the CSYS and USYS structures
    //      For this reason, they are initialized in priority
    //====================================================================================

    //------------------------------------------------------------------------------------
    // - Hard-coded value of the order of the Fourier series. Equals zero for CRTBP
    // - Note that for the non-autonomous models, the order of the Fourier series is still
    // equals to OFS_ORDER.
    // This may need to be changed when some dynamic order will be used, to be able
    // to manipulate Fourier series with an order smaller than the value 30, hard coded
    // in the data files.
    // - Used in the CSYS init functions.
    //------------------------------------------------------------------------------------
    switch(model)
    {
    case Csts::QBCP:
    case Csts::BCP:
        fbpl->n_order_fourier = OFS_ORDER;
        break;
    case Csts::CRTBP:
        fbpl->n_order_fourier = 0;
        break;
    default:
        cout << "init_fbp_lib. Warning: unknown model." << endl;
    }

    //------------------------------------------------------------------------------------
    //Normalization or not. True is always preferable - may be deprecated
    //Kept for compatibility with old code.
    //------------------------------------------------------------------------------------
    fbpl->is_norm = is_norm;

    //------------------------------------------------------------------------------------
    // Libration point
    // Note: passing the complete libration point structure as argument may be done
    // in the near future. Not necessary for now.
    //------------------------------------------------------------------------------------
    fbpl->li_EM  = li_EM;
    fbpl->li_SEM = li_SEM;

    //------------------------------------------------------------------------------------
    // Model - used in the CSYS init functions.
    //------------------------------------------------------------------------------------
    fbpl->model  = model;

    //------------------------------------------------------------------------------------
    // Parameterization style - used in the CSYS init functions.
    //------------------------------------------------------------------------------------
    fbpl->param_style   = param_style;

    //------------------------------------------------------------------------------------
    // Effective order of the Fourier series in RVF (reduced vector field)
    //------------------------------------------------------------------------------------
    switch(param_style)
    {
        case Csts::GRAPH:
        case Csts::MIXED:
            fbpl->eff_nf = fbpl->n_order_fourier;
        break;
        //For normal form, the reduced vector field is quasi autonomous
        //and we can safely reduce the value to 5
        case Csts::NORMFORM:
            fbpl->eff_nf = min(fbpl->n_order_fourier, 5);
        break;
        default:
        cout << "init_fbp_lib. Warning: unknown param_style." << endl;
    }

    //====================================================================================
    // Unit systems
    //====================================================================================
    init_usys(&fbpl->us_em,  Csts::EM,  model);
    init_usys(&fbpl->us_sem, Csts::SEM, model);

    //====================================================================================
    // Coordinate systems
    //====================================================================================
    // Max number of coefficients in the equations of motion
    fbpl->numberOfCoefs = 15;

    // Inform the user if the QBCP is considered 'new'
    // (no data exist for the coefficients of the equation of motion)
    if(is_new) //the QBCP is new, no coefficients exist
    {
        cout << "init_fbp_lib. The considered Sun-Earth-Moon FBP is new:" << endl;
        cout << "coefficients of the vector field are not retrieved from existing files." << endl;
    }


    //------------------------------------------------------------------------------------
    // Coordinate systems for the Earth-Moon libration points
    //------------------------------------------------------------------------------------
    init_csys(&fbpl->cs_em_l1, fbpl, fbp, 1, Csts::EM, man_type_EM, is_new); //L1
    init_csys(&fbpl->cs_em_l2, fbpl, fbp, 2, Csts::EM, man_type_EM, is_new); //L2
    init_csys(&fbpl->cs_em_l3, fbpl, fbp, 3, Csts::EM, man_type_EM, is_new); //L3

    //------------------------------------------------------------------------------------
    // Coordinate systems for the Sun-Earth libration points
    //------------------------------------------------------------------------------------
    init_csys(&fbpl->cs_sem_l1, fbpl, fbp, 1, Csts::SEM, man_type_SEM, is_new); //L1
    init_csys(&fbpl->cs_sem_l2, fbpl, fbp, 2, Csts::SEM, man_type_SEM, is_new); //L2
    init_csys(&fbpl->cs_sem_l3, fbpl, fbp, 3, Csts::SEM, man_type_SEM, is_new); //L3


    //====================================================================================
    // DEFAULT SETTINGS
    // - The int li_EM determines the default Earth-Moon libration point.
    // - The int li_SEM determines the default Sun-(Earth+Moon) libration point.
    // - The int coord_sys determines the default focus; either on the Earth-Moon or
    //   Sun-Earth+Moon framework
    //====================================================================================
    //Coord. syst. for the EM point
    //------------------------------------------------------------------------------------
    switch(li_EM)
    {
    case 1:
        fbpl->cs_em  = fbpl->cs_em_l1;
        break;
    case 2:
        fbpl->cs_em  = fbpl->cs_em_l2;
        break;
    case 3:
        fbpl->cs_em  = fbpl->cs_em_l3;
        break;
    }

    //------------------------------------------------------------------------------------
    //Coord. syst. for the SE point
    //------------------------------------------------------------------------------------
    switch(li_SEM)
    {
    case 1:
        fbpl->cs_sem  = fbpl->cs_sem_l1;
        break;
    case 2:
        fbpl->cs_sem  = fbpl->cs_sem_l2;
        break;
    case 3:
        fbpl->cs_sem  = fbpl->cs_sem_l3;
        break;
    }

    //------------------------------------------------------------------------------------
    //Default Li, Unit & Coord. system
    //------------------------------------------------------------------------------------
    fbpl->coord_sys = coord_sys;
    switch(coord_sys)
    {
    case Csts::EM:
        fbpl->us = fbpl->us_em;
        fbpl->cs = fbpl->cs_em;
        fbpl->li = fbpl->li_EM;
        break;
    case Csts::SEM:
        fbpl->us = fbpl->us_sem;
        fbpl->cs = fbpl->cs_sem;
        fbpl->li = fbpl->li_SEM;
        break;
    default:
        cout << "init_fbp_lib. Warning: unknown coord_sys." << endl;
    }

    //------------------------------------------------------------------------------------
    // 6*6 matrix B used for Floquet transformation
    // Do we need it here?
    //------------------------------------------------------------------------------------
    fbpl->B  = (double*) calloc(36, sizeof(double));
}

/**
 * \brief Initializes a unit system in the form of a USYS structure,
 *        such as Earth-Moon or Sun-(Earth+Moon) unit systems.
 *        Consistent with appendices of the PhD manuscript (October 2017).
 *        The Earth-Moon constants are taken from Andreu 1998.
 * \param usys pointer on the USYS structure to initialize.
 * \param label the type of unit system (Csts:EM, Csts:SEM, etc).
 * \param model the type of model (Csts::QBCP, etc).
 **/
void init_usys(USYS *usys, int label, int model)
{
    //------------------------------------------------------------------------------------
    // These values are used as reference values for all other constants.
    //------------------------------------------------------------------------------------
    //Earth-Moon mass ratio
    usys->mu_EM  = +1.215058162343360e-02;
    //Lunar eccentricity
    usys->lecc   = +0.054900489;

    //------------------------------------------------------------------------------------
    // Physical params in EM units
    // These values are also used
    // as reference values
    // for all other constants
    //------------------------------------------------------------------------------------
    //Relative pulsation of the SEM system
    double n_EM  = 0.925195985520347;
    //Outer pulsation
    double ns_EM = 1.0-n_EM;
    //Sun radius
    double as_EM = 388.81114;
    //Sun mass
    double ms_EM = ns_EM*ns_EM*pow(as_EM,3.0)-1;
    //Earth mass
    double me_EM = 1.0-usys->mu_EM;
    //Moon mass
    double mm_EM = usys->mu_EM;

    //Rest of the parameters
    switch(model)
    {
    case Csts::QBCP:
    case  Csts::BCP:
    {
        switch(label)
        {
        case Csts::EM:
        {
            //Physical params in EM units
            //--------------------------
            usys->n  = n_EM;
            usys->ns = ns_EM;
            usys->ni = 1.0;
            usys->as = as_EM;
            usys->ai = 1.0;
            usys->T  = 2*M_PI/usys->n;
            usys->ms = usys->ns*usys->ns*pow(usys->as,3.0)-1;
            usys->me = 1.0-usys->mu_EM;
            usys->mm = usys->mu_EM;
            break;
        }

        case Csts::SEM:
        {
            //Physical params in SEM units
            //--------------------------
            usys->ns = 1.0;
            usys->ni = 1.0/ns_EM;
            usys->n  = usys->ni-usys->ns;
            usys->as = 1.0;
            usys->ai = 1.0/as_EM;
            usys->T  = 2*M_PI/usys->n;
            usys->ms = ms_EM/(1.0+ms_EM);
            usys->me = me_EM/(1.0+ms_EM);
            usys->mm = mm_EM/(1.0+ms_EM);
            break;
        }

        default :  //Earth-Moon
        {
            cout << "init_usys. Unrecognised unit system. Earth-Moon units are used by default. " << endl;
            //Physical params in EM units
            //--------------------------
            usys->n  = n_EM;
            usys->ns = ns_EM;
            usys->ni = 1.0;
            usys->as = as_EM;
            usys->ai = 1.0;
            usys->T  = 2*M_PI/usys->n;
            usys->ms = usys->ns*usys->ns*pow(usys->as,3.0)-1;
            usys->me = 1.0-usys->mu_EM;
            usys->mm = usys->mu_EM;
        }
        }

        //SEM mass ratio
        usys->mu_SEM = (usys->me+usys->mm)/(usys->ms+usys->me+usys->mm);
        //SE mass ratio
        usys->mu_SE  = (usys->me)/(usys->ms+usys->me);
        break;

    }
    case Csts::CRTBP:
    {
        //SEM mass ratio
        usys->mu_SEM = (me_EM + mm_EM)/(ms_EM + me_EM + mm_EM);
        //SE mass ratio
        usys->mu_SE  = usys->mu_EM*usys->mu_SEM/(1.0 + usys->mu_EM);

        switch(label)
        {
        case Csts::EM:

            //Physical params in EM units
            //--------------------------
            usys->n  = 1.0;
            usys->ns = 0.0;
            usys->ni = 1.0;
            usys->as = as_EM;
            usys->ai = 1.0;
            usys->T  = 2*M_PI/usys->n;
            usys->ms = 0.0;
            usys->me = 1.0-usys->mu_EM;
            usys->mm = usys->mu_EM;
            break;

        case Csts::SEM:

            //Physical params in SEM units
            //--------------------------
            usys->n  = 1.0;
            usys->ns = -1.0;
            usys->ni = 0.0;
            usys->as = 1.0;
            usys->ai = 1.0/as_EM;
            usys->T  = 2*M_PI/usys->n;
            usys->ms = 1.0-usys->mu_SEM;
            usys->me = usys->mu_SEM; //the Earth contains both the Earth's and the Moon's masses.
            usys->mm = 0.0;//mm_EM/(1.0+ms_EM);
            break;

        default:  //Earth-Moon
            cout << "init_usys. Unrecognised unit system. Earth-Moon units are used by default. " << endl;
            //Physical params in EM units
            //--------------------------
            usys->n  = n_EM;
            usys->ns = ns_EM;
            usys->ni = 1.0;
            usys->as = as_EM;
            usys->ai = 1.0;
            usys->T  = 2*M_PI/usys->n;

            usys->ms = 0.0;
            usys->me = 1.0-usys->mu_EM;
            usys->mm = usys->mu_EM;
        }
    }
    }


}

/**
 * \brief Initializes a coordinate systems (CSYS structure), with associated
 *        vector field coefficients, data folder names, and unit system.
 * \param csys pointer on the CSYS structure to initialize.
 * \param fbpl pointer on the FBPL structure that contains csys.
 * \param fbp pointer on the FBP structure that contains parameters specific to
 *        each libration points (namely, the value of gamma)
 * \param li number of the libration point to focus on (L1, L2).
 * \param coord_sys indix of the coordinate system to use (Csts::EM, Csts::SEM).
 * \param man_type: type of manifold about li (Csts::MAN_CENTER, Csts::MAN_CENTER_S...).
 * \param is_new boolean. if true, the qbtbp has not been computed via the qbtbp() routine,
 *        so the vector field coefficients is not initialized.
 *
 *   Note that the FBP structure is used only for the initialization of the coordinate
 *   systems. More precisely, it contains some parameters specific to each libration point
 *  (gamma), via its CR3BP structures.
 **/
void init_csys(CSYS *csys, FBPL *fbpl, FBP *fbp, int li,
               int coord_sys, int man_type, int is_new)
{
    //------------------------------------------------------------------------------------
    // Retrieve variables from fbpl
    //------------------------------------------------------------------------------------
    int n_order_fourier = fbpl->n_order_fourier;
    int model           = fbpl->model;
    int coef_number     = fbpl->numberOfCoefs;
    int param_style     = fbpl->param_style;

    //------------------------------------------------------------------------------------
    //Complete folders, all of the forms "data/..."
    //------------------------------------------------------------------------------------
    csys->F_GS     = init_F_FOLDER("data/CM",    model, coord_sys, li);     //Graph style (PM)
    csys->F_NF     = init_F_FOLDER("data/NF",    model, coord_sys, li);     //Normal form style(PM)
    csys->F_MS     = init_F_FOLDER("data/MS",    model, coord_sys, li);     //Mixed style (PM)

    csys->F_CS     = init_F_FOLDER("data/CS",    model, coord_sys, li);     //Center-stable (PM)
    csys->F_CU     = init_F_FOLDER("data/CU",    model, coord_sys, li);     //Center-unstable (PM)
    csys->F_CUS    = init_F_FOLDER("data/CUS",   model, coord_sys, li);     //Center-hyperbolic (PM)

    csys->F_COEF   = init_F_FOLDER("data/VF",    model, coord_sys, li);     //For integration in a given coord. system
    csys->F_COC    = init_F_FOLDER("data/COC",   model, coord_sys, li);     //For the change of coordinates of the PM
    csys->F_PRINT  = init_F_FOLDER("print",     model, coord_sys, li);      //For printing results output

    //------------------------------------------------------------------------------------
    //Parameterization folders
    //------------------------------------------------------------------------------------
    csys->man_type = man_type;
    switch(man_type)
    {
    //If the Center Manifold is selected,
    //we can choose between the styles
    case Csts::MAN_CENTER:

        switch(param_style)
        {
        case Csts::GRAPH:
            csys->F_PMS = csys->F_GS;
            break;
        case Csts::NORMFORM:
            csys->F_PMS = csys->F_NF;
            break;
        case Csts::MIXED:
            csys->F_PMS = csys->F_MS;
            break;

        }
        break;
    //Else, the style is fixed
    case Csts::MAN_CENTER_S:
        csys->F_PMS = csys->F_CS;
        break;

    case Csts::MAN_CENTER_U:
        csys->F_PMS = csys->F_CU;
        break;

    case Csts::MAN_CENTER_US:
        csys->F_PMS = csys->F_CUS;
        break;
    }

    //------------------------------------------------------------------------------------
    // Unit system associated with csys
    //------------------------------------------------------------------------------------
    switch(coord_sys)
    {
    case Csts::EM:
        csys->us = fbpl->us_em;
        break;
    case Csts::SEM:
        csys->us = fbpl->us_sem;
        break;
    default:
        cout << "init_csys. Warning: unknown model." << endl;
    }

    //------------------------------------------------------------------------------------
    // c1 and gamma, for normalized computation
    //------------------------------------------------------------------------------------
    // First, select the right framework
    // Gives the assosicate CR3BP and mu
    //------------------------------
    CR3BP cr3bp_root;
    switch(coord_sys)
    {
    case Csts::EM:
        cr3bp_root  = fbp->cr3bp1;
        csys->cr3bp = fbp->cr3bp1;
        csys->mu    = fbpl->us_em.mu_EM;
        break;
    case Csts::SEM:
        cr3bp_root  = fbp->cr3bp2;
        csys->cr3bp = fbp->cr3bp2;
        csys->mu    = fbpl->us_sem.mu_SEM;
        break;
    default:
        cout << "init_csys. Warning: unknown model." << endl;
        cr3bp_root  = fbp->cr3bp1; //cr3bp_root is set here to avoid compilation warning, but the default case should NOT be used!
        csys->cr3bp = fbp->cr3bp1;
        csys->mu    = fbpl->us_em.mu_EM;
    }

    //------------------------------
    //Then set the value of c1 and gamma
    //according to the selected libration point
    //------------------------------
    csys->li = li;
    switch(li)
    {
    case 1:
    {
        csys->gamma = cr3bp_root.l1.gamma_i;
        csys->c1    = (csys->mu-1+csys->gamma)/csys->gamma;
        break;
    }

    case 2:
    {
        csys->gamma =  cr3bp_root.l2.gamma_i;
        csys->c1    = (csys->mu-1-csys->gamma)/csys->gamma;
        break;
    }

    case 3: //same as L1/L2. Different from convention by Jorba & Masdemont 1999
    {
        csys->gamma =  cr3bp_root.l3.gamma_i;
        csys->c1    = (csys->mu + csys->gamma)/csys->gamma;
        break;
    }
    }

    //------------------------------------------------------------------------------------
    //c2 coefficient
    //------------------------------------------------------------------------------------
    csys->c2 = cn_coeff(csys->li, csys->gamma, csys->mu, 2);

    //------------------------------------------------------------------------------------
    // Creating the arrays of coefficients
    //------------------------------------------------------------------------------------
    csys->coeffs = (double*) calloc(coef_number*(fbpl->n_order_fourier+1), sizeof(double)); //Default set of vector field coefficients
    csys->Ps  = (double*) calloc(3*(n_order_fourier+1), sizeof(double));  //Sun   position in EM coordinates
    csys->Pe  = (double*) calloc(3*(n_order_fourier+1), sizeof(double));  //Earth position in EM coordinates
    csys->Pm  = (double*) calloc(3*(n_order_fourier+1), sizeof(double));  //Moon  position in EM coordinates
    csys->ps  = (double*) calloc(3*(n_order_fourier+1), sizeof(double));  //Sun   position in NC coordinates
    csys->pe  = (double*) calloc(3*(n_order_fourier+1), sizeof(double));  //Earth position in NC coordinates
    csys->pm  = (double*) calloc(3*(n_order_fourier+1), sizeof(double));  //Moon  position in NC coordinates

    //------------------------------------------------------------------------------------
    // Retrieving the arrays of coefficients
    //------------------------------------------------------------------------------------
    if(~is_new)    //the coefficients can be retrieved
    {
        //The flag comp_type is here to tell if we want the coefficients
        //computed from the FFT or from direct computation
        double comp_type;

        //Switch between the models
        switch(model)
        {
        case Csts::QBCP:
        case Csts::BCP:
        {
            comp_type = 1; //from FFTs
            //----------------------------------------------------------------------------
            // Default set of vector field coefficients
            //----------------------------------------------------------------------------
            read_fourier_coef(csys->F_COEF+"alpha", csys->coeffs, n_order_fourier, 0, comp_type, coef_number);
            //----------------------------------------------------------------------------
            // primaries position
            //----------------------------------------------------------------------------
            read_fourier_coef(csys->F_COEF+"Ps", csys->Ps, n_order_fourier, 0, comp_type, 3);
            read_fourier_coef(csys->F_COEF+"Pe", csys->Pe, n_order_fourier, 0, comp_type, 3);
            read_fourier_coef(csys->F_COEF+"Pm", csys->Pm, n_order_fourier, 0, comp_type, 3);
            read_fourier_coef(csys->F_COEF+"ps", csys->ps, n_order_fourier, 0, comp_type, 3);
            read_fourier_coef(csys->F_COEF+"pe", csys->pe, n_order_fourier, 0, comp_type, 3);
            read_fourier_coef(csys->F_COEF+"pm", csys->pm, n_order_fourier, 0, comp_type, 3);
            break;

        }

        case Csts::CRTBP: //CRTBP case: no need to read from file
        {
            //cout << "init_csys. The use of the CRTBP has been detected." << endl;
            //----------------------------------------------------------------------------
            // Default set of vector field coefficients
            //----------------------------------------------------------------------------
            csys->coeffs[0] = 1.0;
            csys->coeffs[1] = 0.0;
            csys->coeffs[2] = 1.0;
            csys->coeffs[3] = 0.0;
            csys->coeffs[4] = 0.0;
            csys->coeffs[5] = 1.0;

            switch(coord_sys)
            {
                case Csts::EM:
                    //Sun position
                    csys->coeffs[6] = 0.0;
                    csys->coeffs[7] = 0.0;
                    //Earth position
                    csys->coeffs[8] = csys->mu;
                    csys->coeffs[9] = 0.0;
                    //Moon position
                    csys->coeffs[10] = csys->mu-1.0;
                    csys->coeffs[11] = 0.0;
                break;

                case Csts::SEM:
                    //Sun position
                    csys->coeffs[6]  = csys->mu;
                    csys->coeffs[7]  = 0.0;
                    //Earth position
                    csys->coeffs[8]  = csys->mu-1.0;
                    csys->coeffs[9]  = 0.0;
                    //Moon position
                    csys->coeffs[10] = 0.0;
                    csys->coeffs[11] = 0.0;
                break;
            }
            //NC coeffs
            csys->coeffs[12] = -csys->c1;     //alpha13 = alpha[12] = -c1
            csys->coeffs[13] = 0.0;           //alpha14 = alpha[13] = 0


            //----------------------------------------------------------------------------
            // Earth, Moon, and Sun.
            //----------------------------------------------------------------------------
            switch(coord_sys)
            {
                case Csts::EM:
                    csys->Pe[0] = csys->mu;
                    csys->Pe[1] = 0.0;
                    csys->Pe[2] = 0.0;

                    csys->Pm[0] = csys->mu-1.0;
                    csys->Pm[1] = 0.0;
                    csys->Pm[2] = 0.0;

                    csys->Ps[0] = 0.0;
                    csys->Ps[1] = 0.0;
                    csys->Ps[2] = 0.0;
                    break;

                case Csts::SEM:
                    csys->Pe[0] = csys->mu-1.0;
                    csys->Pe[1] = 0.0;
                    csys->Pe[2] = 0.0;

                    csys->Pm[0] = 0.0;
                    csys->Pm[1] = 0.0;
                    csys->Pm[2] = 0.0;

                    csys->Ps[0] = csys->mu;
                    csys->Ps[1] = 0.0;
                    csys->Ps[2] = 0.0;
                    break;
            }


            //From SYS (Earth-Moon or Sun-Earth) to NC coordinates for the primaries
            sys_to_nc_prim(csys->Pe, csys->pe, csys->c1, csys->gamma);
            sys_to_nc_prim(csys->Pm, csys->pm, csys->c1, csys->gamma);
            sys_to_nc_prim(csys->Ps, csys->ps, csys->c1, csys->gamma);

            break;
        }

        default:
            cout << "init_fbp. Warning: unknown model." << endl;
        }

    }


    //------------------------------------------------------------------------------------
    // Storing the solutions of the QBTBP
    // This part could be made more generic, but works for now.
    //------------------------------------------------------------------------------------
    csys->zt = Ofsc(n_order_fourier);
    csys->Zt = Ofsc(n_order_fourier);
    csys->ztdot = Ofsc(n_order_fourier);
    csys->Ztdot = Ofsc(n_order_fourier);
    switch(model)
    {
    case Csts::CRTBP:
    case Csts::BCP:
    {
        //Perfect circles
        csys->zt.set_coef(1.0+0.0*I, 0);
        csys->Zt.set_coef(1.0+0.0*I, 0);
        break;
    }

    case Csts::QBCP:
    default:
    {
        //Taken from files
        string filename = "data/qbtbp/";
        read_ofs_txt(csys->zt, filename+"bjc");
        read_ofs_txt(csys->Zt, filename+"cjc");
        //Derivatives
        csys->ztdot.dot(csys->zt, csys->us.n);
        csys->Ztdot.dot(csys->Zt, csys->us.n);
        break;
    }
    }
}

/**
 * \brief Initialize a Four-Body Problem in the form of a FBP structure.
 * \param fbp pointer on the FBP structure to init.
 * \param first name of the first primary (e.g Sun)
 * \param second name of the second primary (e.g. Earth)
 * \param third name of the third primary  (e.g. Moon)
 *
 * At the end of this routine:
 *    - fbp->cr3bp1 is initialized as the CRTBP of (second, third).
 *    - fbp->cr3bp2 is initialized as the CRTBP of (first, second).
 *
 * WARNING: in the case of the Sun-Earth-Moon system (in fact, the only interesting case
 * for us...), we need to set (SUN, EARTH_AND_MOON), instead of (SUN, EARTH), in
 * fbp->cr3bp2. So keep in mind that if the configuration is SUN/EARTH/MOON,
 * the result will NOT be (first, second) in fbp->cr3bp2, but (first, second+third).
 *
 **/
void init_fbp(FBP* fbp, int first, int second, int third)
{
    //E.g. Earth-Moon init
    init_cr3bp(&fbp->cr3bp1, second, third);
    //E.g. Sun-Earth+Moon init
    //WARNING: in the case of the Sun-Earth-Moon system,
    //we need to set EARTH_AND_MOON, instead of Csts::EARTH alone.
    if(second ==  Csts::EARTH && third == Csts::MOON)
    {
        init_cr3bp(&fbp->cr3bp2, first, Csts::EARTH_AND_MOON);
    }
    else
    {
        init_cr3bp(&fbp->cr3bp2, first, second);
    }
}

/**
 * \brief Initializes a certain Circular Restricted 3-Body Problem as a CR3BP structure.
 * \param cr3bp pointer on the CR3BP structure
 * \param name_1 name of the first primary
 * \param name_2 name of the second primary
 **/
void init_cr3bp(CR3BP *cr3bp, int name_1, int name_2)
{
    //Body initialization
    init_body(&(*cr3bp).m1, name_1);
    init_body(&(*cr3bp).m2, name_2);

    // Set constants
    cr3bp->mu = (*cr3bp).m2.M/( (*cr3bp).m1.M + (*cr3bp).m2.M );  // µ = m2/(m1 + m2)
    cr3bp->L  = (*cr3bp).m2.a;                                    // Distance parameter = semi major axis of m2
    cr3bp->T  = (*cr3bp).m2.T;                                    // Time parameter = sidereal period of m2
    cr3bp->R1 = (*cr3bp).m1.Req;                                  // Radius of m1
    cr3bp->R2 = (*cr3bp).m2.Req;                                  // Radius of m2
    cr3bp->rh = pow((*cr3bp).mu/3,1/3.0);                         // Hill's radius adim formula

    //Li initialization
    init_lib_point(&cr3bp->l1, *cr3bp, 1);
    init_lib_point(&cr3bp->l2, *cr3bp, 2);
    init_lib_point(&cr3bp->l3, *cr3bp, 3);
    init_lib_point(&cr3bp->l4, *cr3bp, 4);
    init_lib_point(&cr3bp->l5, *cr3bp, 5);

    //Name
    strcpy(cr3bp->name, cr3bp->m1.name);
    strcat(cr3bp->name, "-");
    strcat(cr3bp->name, cr3bp->m2.name);
}

/**
 * \brief Initializes a libration point.
 * \param libp a pointer towards the LibPoint structure to initialize.
 * \param cr3bp a CR3BP structure that contains useful coefficients.
 * \param number the indix of the libration point to init.
 **/
void init_lib_point(LibPoint *libp, CR3BP cr3bp, int number)
{
    //Value of gamma
    double gamma_i;

    //Number
    libp->number = number;

    switch(number)
    {
    case 1:
        //Gamma
        gamma_i = cr3bp.rh - 1.0/3.0*pow(cr3bp.rh,2.0)- 1/9*pow(cr3bp.rh,3); //initial guess
        gamma_i = rtnewt(quintic_eq, gamma_i, Config::configManager().G_PREC_LIB(), cr3bp.mu, number);   //newton-raphson method
        libp->gamma_i = gamma_i;

        //Position
        libp->position[0] = 1 - cr3bp.mu - gamma_i;
        libp->position[1] = 0;
        libp->position[2] = 0;
        break;

    case 2:
        //Gamma
        gamma_i = cr3bp.rh + 1.0/3.0*pow(cr3bp.rh,2.0)- 1/9*pow(cr3bp.rh,3);
        gamma_i = rtnewt(quintic_eq, gamma_i, Config::configManager().G_PREC_LIB(), cr3bp.mu, number);
        libp->gamma_i = gamma_i;

        //Position
        libp->position[0] = 1 - cr3bp.mu + gamma_i;
        libp->position[1] = 0;
        libp->position[2] = 0;
        break;

    case 3:
        //Gamma
        gamma_i = 7/12.0*cr3bp.mu + pow(237,2.0)/pow(12,4.0)*pow(cr3bp.mu,3.0);
        gamma_i = rtnewt(quintic_eq, gamma_i, Config::configManager().G_PREC_LIB(), cr3bp.mu, number);
        libp->gamma_i = 1-gamma_i;  //BEWARE: for L3, gamma3 = L3-M1 distance != L3-M2


        //Position
        libp->position[0] = - cr3bp.mu - libp->gamma_i;
        libp->position[1] = 0;
        libp->position[2] = 0;
        break;

    case 4:
        //Gamma
        libp->gamma_i = 1;

        //Position
        libp->position[0] = -cr3bp.mu + 0.5;
        libp->position[1] = sqrt(3)/2.0;
        libp->position[2] = 0;
        break;

    case 5:
        //Gamma
        libp->gamma_i = 1;

        //Position
        libp->position[0] = -cr3bp.mu + 0.5;
        libp->position[1] = -sqrt(3)/2.0;
        libp->position[2] = 0;
        break;
    }

    //Energy & Jacobi constant
    libp->Ei = crtbp_energy(libp->position, cr3bp.mu);
    libp->Ci = -2*libp->Ei;
}

/**
* \brief Initialize one celestial body
* \param body a pointer on the Body structure to init.
* \param name the name of the body in integer format (consistent with SPICE numerotation)
**/
void init_body(Body *body, int name)
{

    double days = 86400; //days to seconds

    switch(name)
    {

    case Csts::MERCURY:

        //Physical parameters
        body->Req = 2439.7;        //[km]
        body->Rm = 2439.7;         //[km]
        body->M = 0.330104e24;     //[kg]
        body->GM = 22032;          //[km^3.s^-2]

        //Orbital parameters
        body->a = 57.91e6;         //[kg]
        body->T = 87.9691*days;    //[s]

        strcpy(body->name, "Mercury");
        break;

    case Csts::VENUS:

        //Physical parameters
        body->Req = 6051.8;        //[km]
        body->Rm = 6501.8;         //[km]
        body->M = 4.86732e24;      //[kg]
        body->GM = 324858.63;      //[km^3.s^-2]

        //Orbital parameters
        body->a = 108.21e6;        //[km]
        body->T = 224.701*days;    //[s]

        strcpy(body->name, "Venus");
        break;


    case Csts::EARTH:

        //Physical parameters
        body->Req = 6378.14;        //[km]
        body->Rm  = 6371.00;        //[km]
        body->M   = 5.97219e24;     //[kg]
        body->GM  = 398600.440;     //[km^3.s^-2]

        //Orbital parameters
        body->a = 149.60e6;          //[km]
        body->T = 365.25636*days;    //[s]

        strcpy(body->name, "Earth");
        break;

    case Csts::MOON: //CONSISTENT WITH HARD CODED VALUE OF mu(Earth-Moon) IN USYS.

        //Physical parameters
        body->Req = 1737.5;       //[km]
        body->Rm = 1737.5;        //[km]
        body->M =  0.07345814120628661e24; //[kg]
        body->GM = 4902.801;      //[km^3.s^-2]

        //Orbital parameters
        body->a = 384400;           //[km]
        body->T = 27.321582*days;    //[s]

        strcpy(body->name, "Moon");
        break;

    case Csts::MARS:

        //Physical parameters
        body->Req = 3396.19;       //[km]
        body->Rm = 3389.50;        //[km]
        body->M = 0.641693e24;     //[kg]
        body->GM = 42828.3;        //[km^3.s^-2]

        //Orbital parameters
        body->a = 227.92e6;       //[kg]
        body->T = 686.98*days;     //[s]

        strcpy(body->name, "Mars");
        break;


    case Csts::JUPITER:

        //Physical parameters
        body->Req =  71492;      //[km]
        body->Rm = 69911;        //[km]
        body->M = 1898.13e24;     //[kg]
        body->GM = 126686511;       //[km^3.s^-2]

        //Orbital parameters
        body->a = 778.57e6;       //[kg]
        body->T = 4332.589*days;     //[s]

        strcpy(body->name, "Jupiter");
        break;

    case Csts::SATURN:

        //Physical parameters
        body->Req =  60268;    //[km]
        body->Rm = 58232;       //[km]
        body->M = 568.319e24;     //[kg]
        body->GM = 37931207.8;      //[km^3.s^-2]

        //Orbital parameters
        body->a = 1433.53e6;       //[kg]
        body->T = 10759.22*days;     //[s]

        strcpy(body->name, "Saturn");
        break;

    case Csts::URANUS:

        //Physical parameters
        body->Req = 25559;      //[km]
        body->Rm = 25362;        //[km]
        body->M =  86.8103e24;    //[kg]
        body->GM =  5793966;      //[km^3.s^-2]

        //Orbital parameters
        body->a =  2872.46e6;      //[kg]
        body->T =  30685.4*days;   //[s]

        strcpy(body->name, "Uranus");
        break;

    case Csts::NEPTUNE:

        //Physical parameters
        body->Req = 24764;      //[km]
        body->Rm = 24622;        //[km]
        body->M = 102.410e24;     //[kg]
        body->GM =  6835107;      //[km^3.s^-2]

        //Orbital parameters
        body->a =  4495.06e6;      //[kg]
        body->T =  60189*days;    //[s]

        strcpy(body->name, "Neptune");
        break;

    case Csts::PLUTO:

        //Physical parameters
        body->Req =  1195;     //[km]
        body->Rm =  1195;       //[km]
        body->M = .01309e24;     //[kg]
        body->GM =  872.4;      //[km^3.s^-2]

        //Orbital parameters
        body->a =  5906.38e6;      //[kg]
        body->T =  90465*days;    //[s]

        strcpy(body->name, "Pluto");
        break;


    case Csts::SUN:

        //Physical parameters
        body->Req = 696342;                //[km]
        body->Rm =  696342;                //[km]
        body->M  = 1988500e24;             //[kg]
        body->GM = 1.3271244004193938e11;  //[km^3.s^-2]

        //Orbital parameters
        body->a = 0;    //[kg]
        body->T = 0;    //[s]

        strcpy(body->name, "Sun");
        break;

    case Csts::EARTH_AND_MOON: //CONSISTENT WITH HARD CODED VALUES OF IN USYS.
        //Equivalent mass of the Earth+Moon system based at the center of mass
        //additionnal physical properties are those of the Earth for consistency)
        //Physical parameters
        body->Req = 6378.14;        //[km]
        body->Rm = 6371.00;         //[km]
        body->M = 6.04590064229622e+24; //[kg]
        body->GM = 398600.440+4902.801; //[km^3.s^-2]

        //Orbital parameters
        body->a = 149.60e6;          //[km]
        body->T = 365.25636*days;     //[s]

        strcpy(body->name, "Earth+Moon");
        break;
    }
}


/**
 *  \brief Initializes two FBPL (model1 and model2) with two different models for
 *         continuation process from one model (epsilon = 0.0) to the other
 *        (epsilon = 1.0).
 **/
void init_fbp_cont(QBCP_I *model, FBPL *model1, FBPL *model2,
                int name_1, int name_2, int name_3, int is_norm, int li_EM, int li_SEM,
                int is_new, int mod1, int mod2, int coord_sys, int param_style)
{
    //Initialize the models
    FBP qbp1, qbp2;
    init_fbp(&qbp1, name_1, name_2, name_3);
    init_fbp(&qbp2, name_1, name_2, name_3);

    //Initialize the models around the given li point
    init_fbp_lib(model1, &qbp1, li_EM, li_SEM, mod1, coord_sys, param_style, Csts::MAN_CENTER, Csts::MAN_CENTER, is_new, is_norm);
    init_fbp_lib(model2, &qbp2, li_EM, li_SEM, mod2, coord_sys, param_style, Csts::MAN_CENTER, Csts::MAN_CENTER, is_new, is_norm);

    //Store in model
    model->model1  = *model1;
    model->model2  = *model2;
    model->epsilon = 0.0;
}


//----------------------------------------------------------------------------------------
//            Change of coordinate systems
//----------------------------------------------------------------------------------------
/**
 *  \brief Change the default coordinate system of the FBPL structure to coord_sys.
 **/
void change_coord(FBPL &fbpl, int coord_sys)
{
    fbpl.coord_sys = coord_sys;
    switch(coord_sys)
    {
    case Csts::EM:
        fbpl.us = fbpl.us_em;
        fbpl.cs = fbpl.cs_em;
        fbpl.li = fbpl.li_EM;
        break;
    case Csts::SEM:
        fbpl.us = fbpl.us_sem;
        fbpl.cs = fbpl.cs_sem;
        fbpl.li = fbpl.li_SEM;
        break;
    default:
        cout << "change_coord. Warning: unknown coord_sys." << endl;
    }
}

/**
 *  \brief Change the default coordinate system to coord_sys and
 *         the libration point to li in the FBPL structure.
 **/
void change_li_coord(FBPL &fbpl, int coord_sys, int li)
{
    //Default settings
    fbpl.coord_sys = coord_sys;  //new default cs
    fbpl.li = li;      //new default libration point

    //Change the coord. system approprietly: "coord_sys" around "li"
    switch(coord_sys)
    {
    case Csts::EM:
        switch(li)
        {
        case 1:
            fbpl.cs_em  = fbpl.cs_em_l1;
            break;
        case 2:
            fbpl.cs_em  = fbpl.cs_em_l2;
            break;
        case 3:
            fbpl.cs_em  = fbpl.cs_em_l3;
            break;
        }
        fbpl.us = fbpl.us_em;
        fbpl.cs = fbpl.cs_em;
        break;
    case Csts::SEM:
        switch(li)
        {
        case 1:
            fbpl.cs_sem  = fbpl.cs_sem_l1;
            break;
        case 2:
            fbpl.cs_sem  = fbpl.cs_sem_l2;
            break;
        case 3:
            fbpl.cs_sem  = fbpl.cs_sem_l3;
            break;
        }
        fbpl.us = fbpl.us_sem;
        fbpl.cs = fbpl.cs_sem;
        break;
    default:
        cout << "change_li_coord. Warning: unknown coord_sys." << endl;
    }
}

//----------------------------------------------------------------------------------------
//            Subroutines - I/O
//----------------------------------------------------------------------------------------
/**
 *  \brief Return the string corresponding to the libration point number provided
 *         (e.g. "L1" if li == 1).
 **/
string init_F_LI(int li)
{
    switch(li)
    {
    case 1:
        return "L1";
    case 2:
        return "L2";
    case 3:
        return "L3";
    case 4:
        return "L4";
    case 5:
        return "L5";
    default:
        cout << "init_F_LI. Warning: supplied libration number is greater than 5." << endl;
    }
    return "L1"; //never here
}

/**
 *  \brief Return the string corresponding to the model indix provided
 *         (e.g. "QBCP" if model == Csts::QBCP).
 **/
string init_F_MODEL(int model)
{
    switch(model)
    {
    case Csts::QBCP:
        return "QBCP";
    case Csts::BCP:
        return "BCP";
    case Csts::CRTBP:
        return "CRTBP";
    case Csts::ERTBP:
        return "ERTBP";
    default:
        cout << "init_F_MODEL. Warning: unknown model." << endl;
    }
    return "QBCP"; //never here
}

/**
 *  \brief Return the string corresponding to the coordinate system
 *         provided (e.g. "EM" if coord_sys == Csts::EM).
 **/
string init_F_COORDSYS(int coord_sys)
{
    switch(coord_sys)
    {
    case Csts::EM:
        return  "EM";
    case Csts::SEM:
        return  "SEM";
    case Csts::SE:
        return  "SE";
    default:
        cout << "init_F_COORDSYS. Warning: unknown model." << endl;
    }
    return "EM"; //never here
}

/**
 *  \brief Return the folder name corresponding to the prefix/model/framework/libration
 *         point number combination provided (e.g. "prefix/QBCP/EM/L1").
 **/
string init_F_FOLDER(string prefix, int model, int coord_sys, int li)
{
    return prefix+"/"+init_F_MODEL(model)+"/"+init_F_COORDSYS(coord_sys)+"/"+init_F_LI(li)+"/";
}

/**
 * \brief Retrieve a set of coefficients, given as Fourier series from a txt file.
 * \param filename the name of the txt file.
 * \param params a pointer toward the array to update.
 * \param n_order_fourier the order of the Fourier series.
 * \param shift the indix from which to start the storage of the coefficients in params.
 * \param flag: if flag == 1, the coefficients computed via Fast Fourier Transform (FFT)
 *        are used. Otherwise, the expansions obtained through Fourier series algebraic
 *        manipulations are used.
 **/
void read_fourier_coef(string filename, double *params, int n_order_fourier, int shift, int flag, int number)
{
    //Reading tools
    ifstream readStream;
    double cDouble1;
    int alphaNumber = 1;
    string ss = static_cast<ostringstream*>( &(ostringstream() << alphaNumber) )->str();
    for(int header = shift; header <= (n_order_fourier+1)*(number-1)+shift; header+=(n_order_fourier+1))
    {
        if(flag) readStream.open((filename+ss+"c_fft.txt").c_str());
        else readStream.open((filename+ss+"c.txt").c_str());
        for(int i=0; i<= n_order_fourier; i++)
        {
            readStream >> cDouble1;  //current order
            readStream >> params[i+header];
        }
        readStream.close();
        alphaNumber++;
        ss = static_cast<ostringstream*>( &(ostringstream() << alphaNumber) )->str();
    }

    if(!flag) cout << "read_fourier_coef: the FFT coefficients have not been used." << endl;
}

//----------------------------------------------------------------------------------------
//            Subroutines - Computation
//----------------------------------------------------------------------------------------
/**
 * \brief Compute the CR3BP potential energy for the given state and mu.
 * \param y  the state array.
 * \param mu the mass ratio of the current CR3BP.
 **/
double crtbp_energy(double y[], double mu)
{
    double r1 = sqrt( pow(y[0]+mu,2.0) + pow(y[1],2.0) + pow(y[2],2.0) );
    double r2 = sqrt( pow(y[0]- 1 + mu,2.0) + pow(y[1],2.0) + pow(y[2],2.0) );
    return - ( 1.0/2.0*(pow(y[0],2) + pow(y[1],2)) + (1-mu)/r1 + mu/r2 + 1.0/2.0*mu*(1-mu) );
}

/**
 * \brief Compute the coefficient cn for a given libration point (L1, L2, and L3 for now)
 * \param fbpl a reference to the FBPL initialized around the selected libration point.
 * \param n the indix of the coefficient to compute.
 *  See double cn_coeff(int li, double gamma, double mu, int n).
 **/
double cn_coeff(FBPL& fbpl, int n)
{
    double mu;
    switch(fbpl.coord_sys)
    {
        case Csts::EM:
            mu   = fbpl.us.mu_EM;
        break;
        case Csts::SEM:
            mu   = fbpl.us.mu_SEM;
        break;
        case Csts::SE:
            mu   = fbpl.us.mu_SE;
        break;
        default: //EM by default
            cout << "WARNING in cn_coeff(): unknown framework. EM by default." << endl;
            mu   = fbpl.us.mu_EM;
    }

    return cn_coeff(fbpl.cs.li, fbpl.cs.gamma, mu, n);
}

/**
 * \brief Compute the coefficient cn for a given libration point (L1, L2, and L3 for now)
 * \param li the number of the current libration point (1,2,3)
 * \param gamma the gamma parameter associated to the current libration point
 * \param mu the mass ratio of the current TBP system
 * \param n the indix of the coefficient to compute.
 *
 * We recall that the cn coefficient are given by the following equation:
 * \f$  c_n = \frac{1}{\gamma_j^3} \left( (\pm1)^n \mu + (-1)^n \frac{(1-\mu) \gamma_j^{(n+1)}}{(1 -\mp \gamma_j)^{(n+1)}} \right) \f$ for \f$ L_j, j = 1,2 \f$, and where:
 * - The upper (resp. lower) signs stand for \f$ L_1 \f$ (resp. \f$ L_2 \f$).
 * - \f$ \gamma_j \f$ is the distance from \f$ L_j \f$ to the smallest primary.
 *
 * The L3 value is taken from Richardson, 1980:
 *
 * \f$  c_n = \frac{1}{\gamma_j^3} \left( 1 - \mu  + \frac{\mu \gamma_j^{n+1}}{(1 + \gamma_j)^{n+1}} \right) \f$.
 **/
double cn_coeff(int li, double gamma, double mu, int n)
{
    double res = 0.0;
    switch(li)
    {
    case 1:
        res =  pow(gamma,-3.0)*(pow(+1.0, n)*mu + pow(-1.0, n)*(1-mu)*pow(gamma/(1.0-gamma), n+1));
        break;
    case 2:
        res =  pow(gamma,-3.0)*(pow(-1.0, n)*mu + pow(-1.0, n)*(1-mu)*pow(gamma/(1.0+gamma), n+1));
        break;
    case 3:
        res =  pow(gamma,-3.0)*(1.0 - mu + mu*pow(gamma/1.0+gamma, n+1)); //convention by Richardson 1980
        break;
    default:
        cout << "cn_coeff. Warning: supplied Li number is out of scope. 0.0 is returned." << endl;
    }
    return res;
}

/**
 * \brief Using the Newton-Raphson method, find the root of a function known to lie close
 *        to the first guess x1. The root will be refined until its accuracy is known
 *        within ± xacc. funcd is a user-supplied routine that returns both the
 *        function value and the first derivative of the function at the point x.
 **/
double rtnewt(void (*funcd)(double, int, double, double *, double *), double x1, double xacc, double mu, int number)
{
    void nrerror(char error_text[]);
    int j;
    double df=0.0, f = 0.0;
    double dx,rtn;

    rtn=x1;   //initial guess

    for (j=1; j<=50; j++)
    {
        (*funcd)(mu, number, rtn,&f,&df);
        dx=f/df;
        rtn -= dx;

        if (fabs(dx) < xacc) return rtn;  //Convergence
    }
    printf("WARNING: Maximum number of iterations exceeded in rtnewt");
    return 0.0;   //Never get here.
}

/**
 * \brief Provides the function value and its first derivative for Newton's method.
 *        f corresponds to the quintic equation satisfied by the Li-m2 distance for
 *        the L1/L2 cases and by 1-(Li-m1 distance) for the L3 case.
 **/
void quintic_eq(double mu, int number, double y, double *f, double *df)
{
    switch(number)
    {

    case 1:
        *f =  pow(y,5.0)   - (3.0-mu)*pow(y,4.0) + (3-2*mu)*pow(y,3.0) - mu*pow(y,2.0) +  2*mu*y - mu;
        *df = 5*pow(y,4.0) - 4*(3.0-mu)*pow(y,3.0) + 3*(3-2*mu)*pow(y,2.0) - 2*mu*y    +  2*mu;
        break;

    case 2:
        *f =  pow(y,5.0)   + (3.0-mu)*pow(y,4.0) + (3-2*mu)*pow(y,3.0) - mu*pow(y,2.0) -  2*mu*y - mu;
        *df = 5*pow(y,4.0) + 4*(3.0-mu)*pow(y,3.0) + 3*(3-2*mu)*pow(y,2.0) - 2*mu*y    -  2*mu;
        break;

    case 3:
        //*f =  pow(y,5.0) + (2.0+mu)*pow(y,4.0) + (1+2*mu)*pow(y,3.0) + (1+mu)*pow(y,2.0) +  2*(1-mu)*y + 1-mu;
        //*df = 5*pow(y,4.0) + 4*(2.0+mu)*pow(y,3.0) + 3*(1+2*mu)*pow(y,2.0) + 2*(1+mu)*y +  2*(1-mu);
        *f= pow(y,5.0) + (7+mu)*pow(y,4.0) + (19+6*mu)*pow(y,3.0) -(24+13*mu)*pow(y,2.0) +  (12+14*mu)*y -7*mu;
        *df= 5*pow(y,4.0) + 4*(7+mu)*pow(y,3.0) + 3*(19+6*mu)*pow(y,2.0) -2*(24+13*mu)*pow(y,1.0) +  (12+14*mu);
        break;
    }
}

/**
 *  \brief This routine performs a change of coordinates:
 *         From SYSTEM (EM or SEM) to NORMALIZED-CENTERED (NC) coordinates.
 *         It is just used here for the primaries.
 */
void sys_to_nc_prim(double Zc[3], double zc[3], double c1, double gamma)
{
    zc[0] = c1 - Zc[0]/gamma;
    zc[1] =    - Zc[1]/gamma;
    zc[2] =    + Zc[2]/gamma;
}

