/**
 * \file init.cpp
 * \brief Allows to initialize the environment (QBCP, CRTBP)
 *        in the form of several global C structures.
 *        The main ones are:
 *          - SEM, a QBCP structure that describes the Sun-Earth-Moon QBCP.
 *          - SEML, a FBPL structure that describes the Sun-Earth-Moon QBCP
 *            about a given libration point.
 * \author BLB.
 */

#include "init.h"
#include <gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>

//========================================================================================
// Global structures
//========================================================================================
FBP SEM;    //global structure that describes the Sun-Earth-Moon system
FBPL SEML;  //global structure that describes the Sun-Earth-Moon system around Li

/**
 *   \brief Initialization of the environnement in the form of global C structures:
 *          SEM, and SEML are updated in order to describe the Sun-(Earth+Moon)
 *          Four-Body Problem about a given libration point.
 *          More precisely:
 *
 *      - SEM is a FBP (Four-Body Problem) structure that contains the numerical
 *        constants of the three primaries + the numerical constants of the associated
 *        Earth-Moon and Sun-(Earth+Moon) CRTBPs.
 *
 *      - SEML is a FBPL structure that contains the numerical constants and parameters
 *        (coefficients of the equations of motion, etc) of the Four-Body model t_model
 *        (QBCP, PACRTBP, BCP) about the libration point t_li (EML1, EML2, SEL1 or SEL2,
 *        as defined in the class Csts).
 *
 *        The structure is also initialized in terms of manifolds.
 *        The integer t_man_type define the type of manifold (stable, center-stable, etc,
 *        as defined in Csts), while the integer t_pms defines the style of
 *        parameterization that is used (graph, normal form... as defined in Csts).
 *
 *    Note that the routine qbtbp must have been run at least once for the
 *    coefficients of the equations of motion to be valid.
 **/
void init_env(int t_model, int t_li, int t_man_type, int t_pms)
{
    //Init the Sun-Earth-Moon problem
    init_fbp(&SEM, Csts::SUN, Csts::EARTH, Csts::MOON);
    //Init the Sun-Earth-Moon problem focused on one libration point
    switch(t_li)
    {
        case Csts::EML1:
            init_fbp_lib(&SEML, &SEM, 1, 1, t_model, Csts::EM, t_pms, t_man_type, t_man_type, false, true);
        break;

        case Csts::EML2:
            init_fbp_lib(&SEML, &SEM, 2, 1, t_model, Csts::EM, t_pms, t_man_type, t_man_type, false, true);
        break;

        case Csts::EML3:
            init_fbp_lib(&SEML, &SEM, 3, 1, t_model, Csts::EM, t_pms, t_man_type, t_man_type, false, true);
        break;

        case Csts::SEL1:
             init_fbp_lib(&SEML, &SEM, 1, 1, t_model, Csts::SEM, t_pms, t_man_type, t_man_type, false, true);
        break;

        case Csts::SEL2:
            init_fbp_lib(&SEML, &SEM,  1, 2, t_model, Csts::SEM, t_pms, t_man_type, t_man_type, false, true);
        break;

        case Csts::SEL3:
            init_fbp_lib(&SEML, &SEM,  1, 3, t_model, Csts::SEM, t_pms, t_man_type, t_man_type, false, true);
        break;

        default:
            cout << "init_env. Error: unknown libration point (t_li)." << endl;
        break;
    }
}


//========================================================================================
// Global variables for parameterization of the center manifold
//========================================================================================
vector<Oftsc>  CM;     //center manifold in NC coordinates
vector<Oftsc>  CMdot;  //time derivative of the center manifold in NC coordinates
vector<Oftsc> CMh;     //center manifold in TFC coordinates
matrix<Oftsc> JCM;     //jacobian of CM
vector<Oftsc>  Fh;     //reduced vector field
vector<Oftsc>  DWFh;   //DWFh =JCM * Fh

/**
 *   \brief Initialization of the parameterization of an invariant manifold around
 *          a given libration point, encoded in the fbpl structure.
 *
 *    The global variables initialized by this routine are:
 *    - vector<Oftsc>  CM, center manifold in NC coordinates
 *    - vector<Oftsc>  CMdot, time derivative of CM
 *    - vector<Oftsc> CMh, center manifold in TFC coordinates
 *    - matrix<Oftsc> JCM, jacobian of CM
 *    - vector<Oftsc>  Fh, reduced vector field
 *    - vector<Oftsc>  DWFh, equal to JCM * Fh
 *
 **/
void init_inv_man(FBPL &fbpl)
{
    //Memory allocation
    CM     = vector<Oftsc>(Csts::NV);
    CMdot  = vector<Oftsc>(Csts::NV);
    CMh    = vector<Oftsc>(Csts::NV);
    DWFh   = vector<Oftsc>(Csts::NV);
    Fh     = vector<Oftsc>(REDUCED_NV);
    JCM    = matrix<Oftsc>(Csts::NV,REDUCED_NV);

    //Update
    update_inv_man(fbpl);
}

/**
 *   \brief Update of the parameterization of an invariant manifold around
 *          a given libration point, encoded in the fbpl structure.
 *
 *    The parameterization is retrieved from text file given in the folder F_PMS,
 *    F_PMS being a string defined in fbpl.
 *    These data files must have been previously computed via the routine pmt.
 *    Moreover, as always, the Fourier-Taylor algebra must have been initialized
 *    by the routine Manip::init in the main program.
 *
 *    The global variables update by this routine are:
 *    - vector<Oftsc>  CM, center manifold in NC coordinates
 *    - vector<Oftsc>  CMdot, time derivative of CM
 *    - vector<Oftsc> CMh, center manifold in TFC coordinates
 *    - matrix<Oftsc> JCM, jacobian of CM
 *    - vector<Oftsc>  Fh, reduced vector field
 *    - vector<Oftsc>  DWFh, equal to JCM * Fh
 *
 **/
void update_inv_man(FBPL &qbcp)
{
    //Read from bin files
    read_vofts_bin(CM,    qbcp.cs.F_PMS+"W/W",         OFS_ORDER);
    read_vofts_bin(CMh,   qbcp.cs.F_PMS+"W/Wh",        OFS_ORDER);
    read_vofts_bin(Fh,    qbcp.cs.F_PMS+"rvf/fh",      OFS_ORDER);
    read_vofts_bin(DWFh,  qbcp.cs.F_PMS+"DWf/C_DWf",   OFS_ORDER);
    read_vofts_bin(CMdot, qbcp.cs.F_PMS+"Wdot/C_Wdot", OFS_ORDER);
    //Building the Jacobian via partial derivation of CM
    for(int i = 0 ; i < Csts::NV ; i++) for(int j = 0 ; j < REDUCED_NV; j++) JCM.der(CM[i], j+1, i, j);   //in NC
}


//========================================================================================
// COC: change of coordinates to achieve a diagonal form at order 2 of the Hamiltonian
//========================================================================================
matrix<Ofsc>  Mcoc;    //COC matrix
matrix<Ofsc>  Pcoc;    //COC matrix (Mcoc = Pcoc*Complex matrix)
matrix<Ofsc>  MIcoc;   //COC matrix = inv(Mcoc)
matrix<Ofsc>  Qcoc;    //COC matrix = inv(Pcoc)
vector<Ofsc>  Vcoc;    //COC vector

/**
 *  Main routine for the initialization of the COC, from TFC to NC coordinates.
 *  This COC is defined as follows:
 *
 *      z = Pcoc * C * zh + Vcoc =  Mcoc * zh + Vcoc,
 *
 *  where z is the NC state, zh is the TFC state, and C is a complex matrix.
 *  The global variable initialized by this routine are then:
 *      - matrix<Ofsc>  Mcoc;    //COC matrix
 *      - matrix<Ofsc>  Pcoc;    //COC matrix (Mcoc = Pcoc*C)
 *      - matrix<Ofsc>  MIcoc;   //COC matrix = inv(Mcoc)
 *      - matrix<Ofsc>  Qcoc;    //COC matrix = inv(Pcoc)
 *      - vector<Ofsc>  Vcoc;    //COC vector
 **/
void init_coc(FBPL &qbcp)
{
    //Memory allocation
    Pcoc  = matrix<Ofsc>(Csts::NV,Csts::NV);
    Mcoc  = matrix<Ofsc>(Csts::NV,Csts::NV);
    Qcoc  = matrix<Ofsc>(Csts::NV,Csts::NV);
    MIcoc = matrix<Ofsc>(Csts::NV,Csts::NV);
    Vcoc  = vector<Ofsc>(Csts::NV);

    //Read from files
    init_coc(Pcoc, Mcoc, Qcoc, MIcoc, Vcoc, qbcp);
}


/**
 *  \brief The several variables of the COC are retrieved from txt files, stored in the
 *         folder F_COC defined in fbpl.
 **/
void init_coc(matrix<Ofsc> &t_Pcoc, matrix<Ofsc> &t_Mcoc, matrix<Ofsc> &t_Qcoc,
              matrix<Ofsc> &t_MIcoc, vector<Ofsc> &t_Vcoc, FBPL& fbpl)
{
    //------------------------------------------------------------------------------------
    //Retrieve folder
    //------------------------------------------------------------------------------------
    string F_COC = fbpl.cs.F_COC;

    //------------------------------------------------------------------------------------
    //Switch CRTBP/QBCP/BCP
    //------------------------------------------------------------------------------------
    if(fbpl.model == Csts::QBCP || fbpl.model == Csts::BCP)
    {
        //Recovering the data: matrix t_Pcoc
        for(int i = 0; i < Csts::NV ; i++) for(int j = 0; j < Csts::NV ; j++) read_coc(*t_Pcoc.get_ptr_first_coef(i,j), F_COC+"P",  i+1, j+1);
        //Recovering the data: matrix t_Qcoc
        for(int i = 0; i < Csts::NV ; i++) for(int j = 0; j < Csts::NV ; j++) read_coc(*t_Qcoc.get_ptr_first_coef(i,j), F_COC+"Q",  i+1, j+1);

        //Recovering the data: vector t_Vcoc
        //t_Vcoc = [G1_11 G1_12 0 G1_21 G1_22 0]^T
        //Note:t_Vcoc[2] and t_Vcoc[5] are kept null
        read_coc(t_Vcoc[0], F_COC+"G1_",  1, 1);
        read_coc(t_Vcoc[1], F_COC+"G1_",  1, 2);
        read_coc(t_Vcoc[3], F_COC+"G1_",  2, 1);
        read_coc(t_Vcoc[4], F_COC+"G1_",  2, 2);
    }
    else //EM RTPB
    {
        //note that t_Vcoc is left untouched (set to zero by default)
        //--------------------------------------------------------------------------------
        //Init double variables
        //--------------------------------------------------------------------------------
        double eta1, eta2, la1, om1, om2, dl1, do1, s1, s2, c2;
        c2 = fbpl.cs.c2;
        eta1 = (c2 - 2.0 - sqrt(9*c2*c2 - 8*c2))/2.0;
        eta2 = (c2 - 2.0 + sqrt(9*c2*c2 - 8*c2))/2.0;
        om1 = sqrt(-eta1);
        la1 = sqrt(+eta2);
        om2 = sqrt(c2);
        dl1 = 2*la1*( (4.0+3*c2)*la1*la1 + 4 + 5*c2 - 6*c2*c2);
        do1 =   om1*( (4.0+3*c2)*om1*om1 - 4 - 5*c2 + 6*c2*c2);
        s1  = sqrt(dl1);
        s2  = sqrt(do1);

        //--------------------------------------------------------------------------------
        //Init t_Pcoc
        //--------------------------------------------------------------------------------
        t_Pcoc.set_coef(+2*la1/s1,                          0, 1);
        t_Pcoc.set_coef(+(la1*la1  - 2*c2 - 1)/s1,          1, 1);
        t_Pcoc.set_coef(+(la1*la1  + 2*c2 + 1)/s1,          3, 1);
        t_Pcoc.set_coef(+(la1*la1*la1 + (1 - 2*c2)*la1)/s1, 4, 1);

        t_Pcoc.set_coef(-(om1*om1 + 2*c2 + 1)/s2,           1, 0);
        t_Pcoc.set_coef(-(om1*om1 - 2*c2 - 1)/s2,           3, 0);

        t_Pcoc.set_coef(-2*la1/s1,                          0, 4);
        t_Pcoc.set_coef(+(la1*la1  - 2*c2 - 1)/s1,          1, 4);
        t_Pcoc.set_coef(+(la1*la1  + 2*c2 + 1)/s1,          3, 4);
        t_Pcoc.set_coef(-(la1*la1*la1 + (1 - 2*c2)*la1)/s1, 4, 4);

        t_Pcoc.set_coef(+2*om1/s2,                          0, 3);
        t_Pcoc.set_coef(-(om1*om1*om1 - (1 - 2*c2)*om1)/s2, 4, 3);

        t_Pcoc.set_coef(+1.0/sqrt(om2),                     2, 2);
        t_Pcoc.set_coef(+sqrt(om2),                         5, 5);

        //--------------------------------------------------------------------------------
        //t_Qcoc = inv(t_Pcoc) (GSL object)
        //--------------------------------------------------------------------------------
        int s;
        gsl_matrix* Pc   = gsl_matrix_calloc (Csts::NV, Csts::NV);
        gsl_matrix* Qc   = gsl_matrix_calloc (Csts::NV, Csts::NV);
        gsl_permutation * p6 = gsl_permutation_alloc (Csts::NV);

        //Init Pc
        for(int i =0; i < Csts::NV; i++) for(int j =0; j < Csts::NV; j++) gsl_matrix_set(Pc, i, j, creal(t_Pcoc.get_ptr_first_coef(i,j)->ofs_get_coef(0)));
        //Use of GSL library
        gsl_linalg_LU_decomp (Pc, p6, &s);
        gsl_linalg_LU_invert (Pc, p6, Qc);

        //--------------------------------------------------------------------------------
        // Init t_Qcoc
        //--------------------------------------------------------------------------------
        for(int i =0; i < Csts::NV; i++) for(int j =0; j < Csts::NV; j++) t_Qcoc.set_coef(gsl_matrix_get(Qc, i, j), i, j);

        //--------------------------------------------------------------------------------
        // Free GSL objects
        //--------------------------------------------------------------------------------
        gsl_matrix_free(Pc);
        gsl_matrix_free(Qc);
        gsl_permutation_free(p6);
    }

    //------------------------------------------------------------------------------------
    // Building t_Mcoc
    //------------------------------------------------------------------------------------
    Ofsc BUX(OFS_ORDER);
    //Keymap used to initialize t_Mcoc and t_MIcoc
    vector<int> keyMap(4);
    keyMap[0] = 0;
    keyMap[1] = 1;
    keyMap[2] = 3;
    keyMap[3] = 4;
    keyMap[4] = 2;
    keyMap[5] = 5;


    int ii;
    //Init t_Mcoc by rows
    for(int i = 0; i <= 3; i++)
    {
        ii = keyMap[i];
        BUX.ofs_fsum(t_Pcoc(ii,0),   1.0/sqrt(2)+0.0*I, t_Pcoc(ii,3), I*1.0/sqrt(2));
        t_Mcoc.set_coef(BUX, ii, 0);
        BUX.ofs_fsum(t_Pcoc(ii,0), I*1.0/sqrt(2), t_Pcoc(ii,3),   1.0/sqrt(2)+0.0*I);
        t_Mcoc.set_coef(BUX, ii, 3);
        t_Mcoc.set_coef(t_Pcoc(ii,1), ii, 1);
        t_Mcoc.set_coef(t_Pcoc(ii,4), ii, 4);
    }

    for(int i = 4; i <= 5; i++)
    {
        ii = keyMap[i];
        BUX.ofs_fsum(t_Pcoc(ii,2),   1.0/sqrt(2)+0.0*I, t_Pcoc(ii,5), I*1.0/sqrt(2));
        t_Mcoc.set_coef(BUX, ii, 2);
        BUX.ofs_fsum(t_Pcoc(ii,2), I*1.0/sqrt(2), t_Pcoc(ii,5),   1.0/sqrt(2)+0.0*I);
        t_Mcoc.set_coef(BUX, ii, 5);
    }

    //Init t_MIcoc by columns
    for(int i = 0; i <= 3; i++)
    {
        ii = keyMap[i];
        BUX.ofs_fsum(t_Qcoc(0,ii),    1.0/sqrt(2)+0.0*I, t_Qcoc(3,ii), -1.0/sqrt(2)*I);
        t_MIcoc.set_coef(BUX, 0, ii);
        BUX.ofs_fsum(t_Qcoc(0,ii), -1.0/sqrt(2)*I, t_Qcoc(3,ii),    1.0/sqrt(2)+0.0*I);
        t_MIcoc.set_coef(BUX, 3, ii);
        t_MIcoc.set_coef(t_Qcoc(1,ii), 1, ii);
        t_MIcoc.set_coef(t_Qcoc(4,ii), 4, ii);
    }

    for(int i = 4; i <= 5; i++)
    {
        ii = keyMap[i];
        BUX.ofs_fsum(t_Qcoc(2,ii),    1.0/sqrt(2)+0.0*I, t_Qcoc(5,ii), -1.0/sqrt(2)*I);
        t_MIcoc.set_coef(BUX, 2, ii);
        BUX.ofs_fsum(t_Qcoc(2,ii), -1.0/sqrt(2)*I, t_Qcoc(5,ii),    1.0/sqrt(2)+0.0*I);
        t_MIcoc.set_coef(BUX, 5, ii);
    }

}


/**
 *  \brief Update a complex Fourier series given as the component (k,p) of
 *         a certain matrix, from a given txt file.
 *  \param xFFT: the Ofsc object to update.
 *  \param prefix: the beginning of the name of the source txt file. e.g. "alpha"
 *  \param k the line indix of the desired component prefix(k,p).
 *  \param p the column indix of the desired component prefix(k,p).
 *
 *  As an example, the call read_coc(xFFT, "alpha", 2, 1) will update xFFT
 *  with the file "alpha21.txt".
 **/
void read_coc(Ofsc& xFFT, string prefix, int k, int p)
{
    //Init
    ifstream readStream;
    string ss1, ss2;
    ss1 = static_cast<ostringstream*>( &(ostringstream() << k) )->str();
    ss2 = static_cast<ostringstream*>( &(ostringstream() << p) )->str();
    //Reading an OFS from a text file
    read_ofs_txt(xFFT, (prefix+ss1+ss2));
}


/**
 *   \brief Number to string inner routine, using static_cast.
 *          Example: num_to_string(10) returns "10".
 **/
string num_to_string(double t_num)
{
    string res =  static_cast<ostringstream*>( &(ostringstream() << t_num) )->str();
    return res;
}
