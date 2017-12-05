#include "pmt.h"
#include <iostream>
#include <fstream>
extern "C" {
#include "nrutil.h"
}

/**
 * \file pmt.cpp
 * \brief Implements the parameterization method for both autonomous and non-autonomous
 *        Hamiltonian vector fields of the Sun-Earth-Moon problem.
 *        Uses time form of the Fourier-Taylor series to perform algebraic operations.
 * \author BLB.
 *
 *  More precisely, the parameterization method is applied to compute a parameterization
 *  of the central manifold of the Earth-Moon libration point L1,2:
 *
 *   - The parameterization is given as Fourier-Taylor expansions in the case
 *     of non-autonomous Hamiltonian vector fields (QBCP, BCP).
 *   - The parameterization is given as pure Taylor expansions in the case
 *     of autonomous Hamiltonian vector fields (CRTBP).
 *
 *  An adaptation to the current implementation to other dynamical models would require:
 *
 *      (i) A computation of the change of coordinates (COC) that brings
 *        the Hamiltonian in normal form at order two - Translated-Floquet-Complexified
 *        (TFC) coordinates (currently done in nfo2.cpp for the QBCP).
 *
 *      (ii) An adaptation of the following routines:
 *
 *      - update_potential/update_potential_m: routines for the computation of the
 *        potential associated with the gravitational influence of the Sun, Earth and Moon
 *      - apply_vector_field: routine for the computation of the vector field.
 */

//----------------------------------------------------------------------------------------
//          Main routine
//----------------------------------------------------------------------------------------
/**
 *  \brief Compute the parameterization of the central manifold of the dynamical
 *         equivalents to the Earth-Moon & Sun-Earth libration points up to a given order.
 *  \param man_type: type of manifold (center, center-stable, center-unstable or center-hyperbolic).
 *  \param param_style: parameterization style (graph, normal form or mixed)
 *  \param small_div_threshold: the limit above which small divisors are discarded.
 *         Only used with the NORMFORM and MIXED styles.
 *  \param is_stored: if true, the results are stored at the end of the computation.
 *
 *   Require inputs:
 *      - The FBPL structure SEML, which contains the addresses of the data folders (F_COC, F_PMS...),
 *        the order of the expansions (both Fourier and Fourier-Taylor),
 *        and some useful parameters (gamma, c1)
 *      - In the case of the QBCP (and BCP), the change of coordinates (COC) that brings
 *        the Hamiltonian in normal form at order two - Translated-Floquet-Complexified
 *        (TFC) coordinates - given in data txt files, computed with the routines
 *        contained in nfo2.cpp.
 *      - In the case of autonomous vector fields, such as the CRTBP, the COC is
 *        computed within the routine pmt.
 *
 *  Outputs (stored in the folder SEML.cs.F_PMS, in binary format):
 *      - W(s,t),  the parameterization of the central manifold in Normalized-Centered
 *       (NC) coordinates.
 *      - Wh(s,t), the parameterization of the central manifold in TFC coordinates.
 *      - FW(s,t), the complete vector field dot(W).
 *      - fh(s,t), the reduced vector field: dot(s) = fh(s,t).
 *      - DWf(s,t), the product DW * fh
 *      - Wdotc(s,t), the partial derivative of W(s,t) with respect to time.
 *        (not equal to dot(W) = dW/dt)
 *
 *  The majority of the operations inside the Fourier-Taylor algebra are performed with
 *  the Fourier series being in the time domain (TFS format, as introduced in ofs.h).
 *  Only the homological equations are computed in OFS (frequency domain) format.
 *
 **/
void pmt(int man_type, int param_style, double small_div_threshold, int is_stored)
{
    cout << "---------------------------------------------------" << endl;
    cout << "                                                   " << endl;
    cout << "    Application of the parameterization method     " << endl;
    cout << "  on the non-autonomous vector field of the QBCP   " << endl;
    cout << "                                                   " << endl;
    cout << "---------------------------------------------------" << endl;

    //------------------------------------------------------------------------------------
    // 0. PM Folder
    //------------------------------------------------------------------------------------
    string F_PMS    = SEML.cs.F_PMS;
    cout << "pmt. The PM will be saved in F_PMS = "   << F_PMS << endl;

    //------------------------------------------------------------------------------------
    // 1. Initialization of all necessary OFTS objects
    //------------------------------------------------------------------------------------

    //---------------------
    // Parameterizations
    //---------------------
    vector<Oftsc> Wh(6);    //Invariant manifold in TFC coordinates
    vector<Oftsc> W(6);     //W   = P(t)*Wh + V(t) (in NC coordinates)
    vector<Oftsc> E(6);     //E   = [DWhc<m fh<m]m-[Fh(Wh<m)]m
    vector<Oftsc> eta(6);   //eta = Hinv*E
    vector<Oftsc> xi(6);    //xi  = Hinv*Wh at order >=2
    //---------------------
    // DWh, DW (diff)
    //---------------------
    matrix<Oftsc> DWhc(6,REDUCED_NV);  //TFC, contribution of order < k to k
    matrix<Oftsc> DW(6,REDUCED_NV);    //NC,  contribution of order <=k to k

    //---------------------
    // Complete VF
    //---------------------
    vector<Oftsc> FW(6);      //NC,   contribution of order <=k to k
    vector<Oftsc> FWhc(6);    //TFC,  contribution of order < k to k
    vector<Oftsc> FWc(6);     //NC,   contribution of order < k to k
    vector<Oftsc> FW2c(6);    //Temp, for contribution of order < k to k
    vector<Oftsc> FW3c(6);    //Temp, for contribution of order < k to k
    //---------------------
    // Reduced VF
    //---------------------
    vector<Oftsc> fh(REDUCED_NV); //TFC & NC (the same in both coord. systems)
    //---------------------
    // DW x fh
    //---------------------
    vector<Oftsc> DWfhc(6);   //TFC, contribution of order < k to k
    vector<Oftsc> DWf(6);     //NC, contribution of order <=k to k
    //---------------------
    // dot(W)
    //---------------------
    vector<Oftsc> Wdot(6);    //NC, contribution of order <=k to k
    //---------------------
    // Potential (of the three primaries combined)
    //---------------------
    vector<Oftsc> Un(6);      //NC, contribution to order <=k to k
    vector<Oftsc> PrimPt(6);  //Intermediate steps for the computation of the gravitational potential
    vector<Oftsc> PrimPt2(6); //Intermediate steps for the computation of the gravitational potential

    //------------------------------------------------------------------------------------
    // 2. Initialization of the COC between TFC and NC coordinates
    //------------------------------------------------------------------------------------
    matrix<Ofsc> P(6,6);      //The matrix P of the c.o.c. (Floquet part)
    matrix<Ofsc> Q(6,6);      //The matrix Q of the c.o.c. (Floquet part). Q = inv(P)
    vector<Ofsc> V(6);        //The vector V of the c.o.c. (Translation part)
    matrix<Ofsc> PC(6,6);     //PC=P(theta)*C
    matrix<Ofsc> PCdot(6,6);  //PCdot=dot(P*C)=dot(P)C
    matrix<Ofsc> CQ(6,6);     //CQ=Cinv*Pinv=inv(PC)
    vector<Ofsc> Vdot(6);     //Vdot=dot(V)
    vector<Ofsc> Xe(3);       //The vectors Xe, Xm, Xs, true position
    vector<Ofsc> Xm(3);       //of the primaries used
    vector<Ofsc> Xs(3);       //to compute the potential
    Ofsc ILe(OFS_ORDER);      //The Ofs ILe, ILm and ILs, true inverse orbit radii
    Ofsc ILm(OFS_ORDER);      //of the primaries used to compute the potential
    Ofsc ILs(OFS_ORDER);      //ILf = 1.0/||Xf||, for f = e, m, s

    //Init routine
    init_coc_ofs(P, Q, PC, PCdot, CQ, Xe, Xm, Xs, V, Vdot, ILe, ILm, ILs, SEML);
    //TFS format of the COC
    tfs_from_ofs(P, Q, PC, PCdot, CQ, Xe, Xm, Xs, V, Vdot, ILe, ILm, ILs);

    //------------------------------------------------------------------------------------
    // 3. Initialization of the alphas (coefficients of the vector field)
    //------------------------------------------------------------------------------------
    int noc = SEML.numberOfCoefs;
    vector<Ofsc> alpha(noc);
    //Init routine
    tfts_init_vf(alpha);
    //TFS format
    tfs_from_ofs_inline(alpha);


    //------------------------------------------------------------------------------------
    // 4. Initialization of the order zero and one:
    //------------------------------------------------------------------------------------
    //------------------------------------------
    // Initialization of the matrices DB, H, L...
    // that encode the order one of the vector field in TFC coordinates.
    // See folded comments below for details
    //------------------------------------------
    /*
    //     | iw1 0    0    0    0   0   |
    //     | 0   w2   0    0    0   0   |
    //     | 0   0    iw3  0    0   0   |
    //DB = | 0   0    0   -iw1  0   0   |
    //     | 0   0    0    0   -w2  0   |
    //     | 0   0    0    0    0  -iw3 |
    //
    //DF0 = DB in matrix<Ofsc> format
    //
    //    |  iw1 0    0    0    0   0  |
    //    |  0   0    0    0    w2  0  |
    //    |  0   iw3  0    0    0   0  |
    //H = |  0   0   -iw1  0    0   0  | = (L N)
    //    |  0   0    0    0    0  -w2 |
    //    |  0   0    0   -iw3  0   0  |
    //
    //H2 = H in matrix<Ofsc> format
    //
    //    |  iw1 0    0    0    |
    //    |  0   0    0    0    |
    //    |  0   iw3  0    0    |
    //L = |  0   0   -iw1  0    |
    //    |  0   0    0    0    |
    //    |  0   0    0   -iw3  |
    //
    //L2 = L in matrix<Ofsc> format
    //
    //    |  0   0  |
    //    |  w2  0  |
    //    |  0   0  |
    //N = |  0   0  |
    //    |  0  -w2 |
    //    |  0   0  |
    //
    //       |1/iw1  0    0      0     0     0    |
    //       | 0     0   1/iw3   0     0     0    |
    //       | 0     0    0    -1/iw1  0     0    |
    //Hinv = | 0     0    0      0     0  -1/iw3  |
    //       | 0    1/w2  0      0     0     0    |
    //       | 0     0    0      0   -1/w2   0    |
    //
    //Hinv2 = Hinv in matrix<Ofsc> format
    */
    //------------------------------------------
    gsl_matrix_complex* DB    = gsl_matrix_complex_calloc(6, 6);
    gsl_matrix_complex* H     = gsl_matrix_complex_calloc(6, 6);
    gsl_matrix_complex* Hinv  = gsl_matrix_complex_calloc(6, 6);
    gsl_matrix_complex* La    = gsl_matrix_complex_calloc(6, 6);
    matrix<Ofsc> DF0(6, 6);
    matrix<Ofsc> H2(6, 6);
    matrix<Ofsc> Hinv2(6, 6);

    //Init routine
    tfts_init_order_one(DB, H, Hinv, La, DF0, H2, Hinv2, man_type);

    //------------------------------------------
    // Initialization of the parameterization
    // TFS formatting is included
    //------------------------------------------
    tfts_init_pm(W, Wh, fh, PC, V, H, La, man_type);

    //------------------------------------------------------------------------------------
    // 5. Initialization of various objects used to solve the pm equations
    //------------------------------------------------------------------------------------
    double Ke = SEML.us.me/pow(SEML.cs.gamma, 3.0);  //Earth potential factor
    double Km = SEML.us.mm/pow(SEML.cs.gamma, 3.0);  //Moon  potential factor
    double Ks = SEML.us.ms/pow(SEML.cs.gamma, 3.0);  //Sun   potential factor

    //------------------------------------------------------------------------------------
    // 6.  Indices for mixed style. The maximum number of invariant submanifold is 50 by default.
    //    (far enough to comprise any desire).
    //     In practice, there is only ims submanifolds, and for j = 1,..., ims, VIs[j]
    //     contains the set of VIn[j] directions that must be null for the submanifold j
    //     to be invariant.
    //------------------------------------------------------------------------------------
    int** VIs, *VIn, ims;
    VIs = (int**) calloc(50, sizeof(int*));
    VIn = (int*)  calloc(50, sizeof(int));
    init_mixed_style(VIs, VIn, &ims, man_type);

    //------------------------------------------------------------------------------------
    // 6. Computing PM at order 0 and 1
    //------------------------------------------------------------------------------------
    for(int m = 0; m < 2 ; m++)
    {
        //-----------------
        //Update potential at order m
        //-----------------
        update_potential(W, Xe, Xm, Xs, PrimPt, Ke, Km, Ks, Un, m);

        //-----------------
        //vector field at order m in NC coordinates
        //-----------------
        apply_vector_field(alpha, W, FWc, Un, m);  //contributions to orders < m to order m
        apply_vector_field(alpha, W, FW,  Un, m);  //contributions to orders <=m to order m

        //-----------------
        //vector field at order m in TFC coordinates
        //FWhc = COCinv(FWc)
        //-----------------
        apply_inv_coc_dot_tfts(CQ, PCdot, Vdot, Wh, FWc, FWhc, FW2c, FW3c, m);

        //-----------------
        //Jacobian DW/DWhc
        //-----------------
        if(m > 0) for(int i = 0 ; i < Csts::NV ; i++) for(int j = 0 ; j < REDUCED_NV; j++) DWhc.tfts_der(Wh[i], j+1, i, j, m);   //in TFC
        if(m > 0) for(int i = 0 ; i < Csts::NV ; i++) for(int j = 0 ; j < REDUCED_NV; j++) DW.tfts_der(W[i], j+1, i, j, m);      //in NC

        //-----------------
        // DWhc x f at order m
        //-----------------
        tfts_smvprod_t(DWhc, fh, DWfhc, m);       //in TFC
        tfts_smvprod_t(DW,  fh, DWf,  m);         //in NC
    }

    //------------------------------------------
    // Update the PrimPt2 at order 0
    //------------------------------------------
    for(int i = 0; i < 6; i++) PrimPt2[i].ccopy(PrimPt[i], 0);

    //------------------------------------------------------------------------------------
    // 7. Computing at orders >= 2
    //------------------------------------------------------------------------------------
    tic();
    for(int m = 2; m <= OFTS_ORDER; m++)
    {
        cout << "pmt. start of order " << m << endl;
        //------------------------------------------
        //The potential Un is updated.
        //What is updated: contributions of order < m to order m
        //------------------------------------------
        update_potential(W, Xe, Xm, Xs, PrimPt, Ke, Km, Ks, Un, m);

        //------------------------------------------
        // vector field at order m in NC coordinates: [F(W<m)]m
        // What is updated: contributions of order < m to order m
        //------------------------------------------
        //W is used in the vector field
        apply_vector_field(alpha, W, FWc, Un, m);

        //------------------------------------------
        //vector field at order m in TFC coordinates
        //FWhc = COCinv(FWc)
        //------------------------------------------
        apply_inv_coc_dot_tfts(CQ, PCdot, Vdot, Wh, FWc, FWhc, FW2c, FW3c, m);

        //------------------------------------------
        // [DW<m x f<m]m
        // What is updated: contributions of order < m to order m
        //------------------------------------------
        tfts_smvprod_t(DWhc, fh, DWfhc, m);

        //------------------------------------------
        // We have [Fh(Wh<m)]m and [DWhc<m fh<m]m
        // Update of E at order m
        //------------------------------------------
        // [E]m = [DWhc<m fh<m]m - [Fh(Wh<m)]m
        for(int i = 0 ; i< 6; i++) E[i].tfts_fsum_u(FWhc[i], -1.0+0.0*I, DWfhc[i], 1.0+0.0*I, m);

        //------------------------------------------
        // E to OFS format
        //------------------------------------------
        tfs_to_ofs_inline(E, m);

        //--------------------------------------------------------------------------------
        // START OF OFS FORMAT
        //--------------------------------------------------------------------------------

        //------------------------------------------
        //[eta]m = [Hinv*E]m;
        //------------------------------------------
        smvprod_u(Hinv2, E, eta, m);

        //------------------------------------------
        // Solving the cohomological equations
        //------------------------------------------
        solv_hom_eq(eta, xi, fh, La, m, param_style, small_div_threshold, VIs, VIn, ims);

        //------------------------------------------
        // Wh = H2*xi @ order m
        // Recall that H2 = H in matrix format
        //------------------------------------------
        smvprod_u(H2, xi, Wh, m);

        //------------------------------------------
        // Wh to TFS format
        //------------------------------------------
        tfs_from_ofs_inline(Wh, m);
        tfs_from_ofs_inline(fh, m);

        //--------------------------------------------------------------------------------
        // END OF OFS FORMAT
        //--------------------------------------------------------------------------------

        //------------------------------------------
        // Update W @ order m
        //------------------------------------------
        apply_coc_tfts(PC,  V,  Wh, W,  m, 1);     //order zero is added


        //------------------------------------------
        // Update the differential to take into account the new terms @ order m
        //------------------------------------------
        for(int i = 0 ; i < Csts::NV ; i++) for(int j = 0 ; j < REDUCED_NV; j++) DWhc.tfts_der(Wh[i], j+1, i, j, m);

        //------------------------------------------
        //The potential is updated: Un
        //What is updated: contributions of order >= m to order m
        //------------------------------------------
        update_potential_m(W, Xe, Xm, Xs, PrimPt, PrimPt2, Ke, Km, Ks, Un, m);

        //--------------------------------------------------------------------------------
        //
        // Note: after this point, only necessary for testing, can be discarded for speed
        //
        //---------------------------------------------------------------------------------
        //------------------------------------------
        // Differential of Wh
        //------------------------------------------
        for(int i = 0 ; i < Csts::NV ; i++) for(int j = 0 ; j < REDUCED_NV; j++) DW.tfts_der(W[i], j+1, i, j, m);

        //------------------------------------------
        // DWhc times f @order m
        //------------------------------------------
        tfts_smvprod_t(DW,  fh, DWf,  m);

        //------------------------------------------
        // Applying the VF: [F(W<=m)]m
        // What is updated: contributions of order <= m to order m
        //------------------------------------------
        apply_vector_field(alpha, W, FW, Un, m);

    }
    double tt = toc();
    cout << "  pmt. End of computation" <<  " in " << tt << " s (" << tt/60 << " mn)." << endl;

    //------------------------------------------------------------------------------------
    // 8. Back to TFS format at the very end.
    //------------------------------------------------------------------------------------
    tic();
    tfs_to_ofs_inline(W);
    tfs_to_ofs_inline(Wh);
    tfs_to_ofs_inline(fh);
    tfs_to_ofs_inline(DWf);
    tfs_to_ofs_inline(DWhc);
    cout << "  pmt. End of TFS formating" <<  " in " << toc() << " s. " << endl;

    //------------------------------------------
    // dot(W)
    //------------------------------------------
    tic();
    for(int m = 0; m <= OFTS_ORDER; m++) for(int i = 0; i < 6; i++) Wdot[i].dot(W[i], SEML.us.n, m);
    cout << "  pmt. End of dot" <<  " in " << toc() << " s. " << endl;

    //------------------------------------------------------------------------------------
    // 9. Printing
    //------------------------------------------------------------------------------------
    if(is_stored)
    {
        tic();
        //--------------------------------------------------------------------------------
        // Binary
        //--------------------------------------------------------------------------------

        //------------------------------------------
        //Vectors
        //------------------------------------------
        write_vofts_bin(W,    F_PMS+"W/W");
        write_vofts_bin(Wh,   F_PMS+"W/Wh");
        write_vofts_bin(fh,   F_PMS+"rvf/fh");
        write_vofts_bin(Wdot, F_PMS+"Wdot/C_Wdot");
        write_vofts_bin(DWf,  F_PMS+"DWf/C_DWf");
        write_vofts_bin(FW,   F_PMS+"FW/C_FW");

        //------------------------------------------
        // Jacobian matrix (TFC coordinates only)
        //------------------------------------------
        write_mofts_bin(DWhc, F_PMS+"DWf/DWhc");

        //------------------------------------------
        //If the manifold is CS or CU, the graph style used by default in the center part
        //and certain directions of the parameterization (0 and 3) can be stored in
        //one-dimensional series. Moreover, the Jacobian is computed and stored.
        //------------------------------------------
        if(man_type == Csts::MAN_CENTER_S || man_type == Csts::MAN_CENTER_U)
        {
            //Init
            Oftsc W1(1, OFTS_ORDER, Csts::OFS_NV, OFS_ORDER);
            Oftsc DW1(1, OFTS_ORDER, Csts::OFS_NV, OFS_ORDER);

            //Transform and store
            from_vofts_to_vofts_bin(Wh, W1, DW1, F_PMS+"W/F");
        }

        //--------------------------------------------------------------------------------
        // Txt: not used for now
        //--------------------------------------------------------------------------------
        //write_vofts_txt(Wh,   F_PMS+"W/Wh");
        //print W & FWc
        //write_vofts_txt(W,     F_PMS+"W/W");
        //write_vofts_txt(E,    F_PMS+"W/E");
        //vector_fprinf(fh, F_PMS+"rvf/fh");

        //--------------------------------------------------------------------------------
        // Info: stored in F_PMS+"INFO.txt"
        //--------------------------------------------------------------------------------
        cout << "  pmt. Info on the outputs are stored in " << F_PMS+"INFO.txt" << endl;
        ofstream myfile;
        string nfile = F_PMS+"INFO.txt";
        myfile.open (nfile.c_str());
        myfile << "OFTS_ORDER = " << OFTS_ORDER << endl;
        myfile << "OFS_ORDER = " << OFS_ORDER << endl;
        myfile.close();

        cout << "  pmt. End of printing" <<  " in " << toc() << " s. " << endl;
    }

    //------------------------------------------------------------------------------------
    // 10. Free
    //------------------------------------------------------------------------------------
    free_mixed_style(VIs, VIn, ims);
}

//----------------------------------------------------------------------------------------
//         Init routines
//----------------------------------------------------------------------------------------
/**
 *  \brief Initialization of the vector fields (alpha coefficients) used in pmt.
 **/
void tfts_init_vf(vector<Ofsc>& alpha)
{
    string F_COEF  = SEML.cs.F_COEF;    //VF folder
    int noc = SEML.numberOfCoefs;       //number of coefficients to retrieve

    //Switch on the model to initialize the alphas
    switch(SEML.model)
    {
    case Csts::QBCP:
    case Csts::BCP:
    {
        //Read from txt files
        ifstream readStream;
        string ss1;
        for(int i = 0; i < noc; i++)
        {
            ss1 = static_cast<ostringstream*>( &(ostringstream() << i+1) )->str();
            read_ofs_txt(alpha[i], F_COEF+"alpha"+ss1+"_fft");
        }

        break;
    }

    case Csts::CRTBP:
    {
        cout << "pmt. The use of the CRTBP has been detected when initializing the alphas." << endl;
        //All coeffs to zero
        for(int j = 0; j < noc;  j++) alpha[j].set_coef(0.0, 0);
        //Non-null coeffs
        alpha[0].set_coef(1.0, 0);               //alpha1  = 1.0
        alpha[2].set_coef(1.0, 0);               //alpha3  = 1.0
        alpha[5].set_coef(1.0, 0);               //alpha6  = 1.0
        alpha[12].set_coef(-SEML.cs.c1, 0);      //alpha13 = -c1

        //Not useful here:
        //        alpha[8].set_coef(SEML.cs.mu, 0);        //alpha9  = Xe
        //        alpha[10].set_coef(SEML.cs.mu, 0);       //alpha11 = Xm

        break;
    }

    default:
    {
        cout << "error in pmt. Unknown model." << endl;
        return;
    }

    }
}

/**
 *  \brief Initialization of various matrices that helps build the order one of the parameterization. Namely:
 *
 *      | iw1 0    0    0    0   0   |
 *      | 0   w2   0    0    0   0   |
 *      | 0   0    iw3  0    0   0   |
 * DB = | 0   0    0   -iw1  0   0   |
 *      | 0   0    0    0   -w2  0   |
 *      | 0   0    0    0    0  -iw3 |
 *
 * DF0 = DB in matrix<Ofsc> format
 *
 *     |  iw1 0    0    0    0   0  |
 *     |  0   0    0    0    w2  0  |
 *     |  0   iw3  0    0    0   0  |
 * H = |  0   0   -iw1  0    0   0  | = (L N)
 *     |  0   0    0    0    0  -w2 |
 *     |  0   0    0   -iw3  0   0  |
 *
 * H2 = H in matrix<Ofsc> format
 *
 *     |  iw1 0    0    0    |
 *     |  0   0    0    0    |
 *     |  0   iw3  0    0    |
 * L = |  0   0   -iw1  0    |
 *     |  0   0    0    0    |
 *     |  0   0    0   -iw3  |
 *
 * L2 = L in matrix<Ofsc> format
 *
 *     |  0   0  |
 *     |  w2  0  |
 *     |  0   0  |
 * N = |  0   0  |
 *     |  0  -w2 |
 *     |  0   0  |
 *
 *        |1/iw1  0    0      0     0     0    |
 *        | 0     0   1/iw3   0     0     0    |
 *        | 0     0    0    -1/iw1  0     0    |
 * Hinv = | 0     0    0      0     0  -1/iw3  |
 *        | 0    1/w2  0      0     0     0    |
 *        | 0     0    0      0   -1/w2   0    |
 *
 * with (w1, w2, w3) the frequency associated to the xy elliptic motion, hyperbolic motion and z elliptic motion, respectively.
 *
 *  WARNING: all these objects are initialized in OFS format, and NOT in TFS format, since they are only used in pure OFS computations.
 **/
void tfts_init_order_one(gsl_matrix_complex* DB,
                         gsl_matrix_complex* H,
                         gsl_matrix_complex* Hinv,
                         gsl_matrix_complex* La,
                         matrix<Ofsc>& DF0,
                         matrix<Ofsc>& H2,
                         matrix<Ofsc>& Hinv2,
                         int man_type)
{
    //------------------------------------------
    // Init
    //------------------------------------------
    string F_COC   = SEML.cs.F_COC;              //COC folder
    gsl_complex one_c  = gslc_complex(1.0, 0.0); //one
    gsl_complex zero_c = gslc_complex(0.0, 0.0); //zero

    //------------------------------------------
    // Getting DB from files
    //------------------------------------------
    //     | iw1 0    0    0    0   0   |
    //     | 0   w2   0    0    0   0   |
    //     | 0   0    iw3  0    0   0   |
    //DB = | 0   0    0   -iw1  0   0   |
    //     | 0   0    0    0   -w2  0   |
    //     | 0   0    0    0    0  -iw3 |
    //------------------------------------------
    if(SEML.model == Csts::QBCP || SEML.model == Csts::BCP)
    {
        //Reading an OFS from a text file
        string filename = F_COC+"DB";
        glsc_matrix_complex_read(DB, filename);

    }
    else //CRTBP
    {
        cout << "pmt. The use of the CRTBP has been detected when initializing DB." << endl;
        //--------------------
        //Init double variables
        //--------------------
        double eta1, eta2, la1, om1, om2, c2;
        c2 = SEML.cs.c2;
        eta1 = (c2 - 2 - sqrt(9*c2*c2 - 8*c2))/2.0;
        eta2 = (c2 - 2 + sqrt(9*c2*c2 - 8*c2))/2.0;
        om1 = sqrt(-eta1);
        la1 = sqrt(+eta2);
        om2 = sqrt(c2);
        //--------------------
        //Init DB
        //--------------------
        gsl_matrix_complex_set(DB, 0, 0, gslc_complex(+om1*I));
        gsl_matrix_complex_set(DB, 1, 1, gslc_complex(+la1+0.0*I));
        gsl_matrix_complex_set(DB, 2, 2, gslc_complex(+om2*I));
        gsl_matrix_complex_set(DB, 3, 3, gslc_complex(-om1*I));
        gsl_matrix_complex_set(DB, 4, 4, gslc_complex(-la1+0.0*I));
        gsl_matrix_complex_set(DB, 5, 5, gslc_complex(-om2*I));
    }

//    double omt = creal(gslc_complex(gsl_matrix_complex_get(DB, 1, 1)));
//    cout << "om2 = " << omt << endl;
//    cout << "la2 = " << exp(omt*SEML.us.T);

    //------------------------------------------
    //DF0 = DB in matrix<Ofsc> format
    //------------------------------------------
    Ofsc BUX(OFS_ORDER);
    BUX.set_coef(gslc_complex(gsl_matrix_complex_get(DB, 0, 0)), 0);
    DF0.set_coef(BUX, 0, 0);
    BUX.set_coef(gslc_complex(gsl_matrix_complex_get(DB, 1, 1)), 0);
    DF0.set_coef(BUX, 1, 1);
    BUX.set_coef(gslc_complex(gsl_matrix_complex_get(DB, 2, 2)), 0);
    DF0.set_coef(BUX, 2, 2);
    BUX.set_coef(gslc_complex(gsl_matrix_complex_get(DB, 3, 3)), 0);
    DF0.set_coef(BUX, 3, 3);
    BUX.set_coef(gslc_complex(gsl_matrix_complex_get(DB, 4, 4)), 0);
    DF0.set_coef(BUX, 4, 4);
    BUX.set_coef(gslc_complex(gsl_matrix_complex_get(DB, 5, 5)), 0);
    DF0.set_coef(BUX, 5, 5);


    //------------------------------------------
    // Setting H, Hinv, Lambda, LamdaL, LambdaN from DB
    // A this step, the central manifold is selected.
    //
    // TO-DO: something more generic to be able
    // to compute the central-stable and central-unstable manifold
    //
    //------------------------------------------

    //------------------------------------------
    //    |  iw1 0    0    0    0   0  |
    //    |  0   0    0    0    w2  0  |
    //    |  0   iw3  0    0    0   0  |
    //H = |  0   0   -iw1  0    0   0  | = (L N)
    //    |  0   0    0    0    0  -w2 |
    //    |  0   0    0   -iw3  0   0  |
    //------------------------------------------
    switch(man_type)
    {
    case Csts::MAN_CENTER:
    case Csts::MAN_CENTER_U:
    case Csts::MAN_CENTER_US:
    {
        gsl_matrix_complex_set(H, 0, 0, gsl_matrix_complex_get(DB, 0, 0));
        gsl_matrix_complex_set(H, 1, 4, gsl_matrix_complex_get(DB, 1, 1));
        gsl_matrix_complex_set(H, 2, 1, gsl_matrix_complex_get(DB, 2, 2));
        gsl_matrix_complex_set(H, 3, 2, gsl_matrix_complex_get(DB, 3, 3));
        gsl_matrix_complex_set(H, 4, 5, gsl_matrix_complex_get(DB, 4, 4));
        gsl_matrix_complex_set(H, 5, 3, gsl_matrix_complex_get(DB, 5, 5));

        //------------------------------------------
        //H2 = H in matrix<Ofsc> format
        //------------------------------------------
        BUX.zero();
        BUX.set_coef(gslc_complex(gsl_matrix_complex_get(DB, 0, 0)), 0);
        H2.set_coef(BUX, 0, 0);
        BUX.set_coef(gslc_complex(gsl_matrix_complex_get(DB, 1, 1)), 0);
        H2.set_coef(BUX, 1, 4);
        BUX.set_coef(gslc_complex(gsl_matrix_complex_get(DB, 2, 2)), 0);
        H2.set_coef(BUX, 2, 1);
        BUX.set_coef(gslc_complex(gsl_matrix_complex_get(DB, 3, 3)), 0);
        H2.set_coef(BUX, 3, 2);
        BUX.set_coef(gslc_complex(gsl_matrix_complex_get(DB, 4, 4)), 0);
        H2.set_coef(BUX, 4, 5);
        BUX.set_coef(gslc_complex(gsl_matrix_complex_get(DB, 5, 5)), 0);
        H2.set_coef(BUX, 5, 3);


        //------------------------------------------
        //       |1/iw1  0    0      0     0     0    |
        //       | 0     0   1/iw3   0     0     0    |
        //       | 0     0    0    -1/iw1  0     0    |
        //Hinv = | 0     0    0      0     0  -1/iw3  |
        //       | 0    1/w2  0      0     0     0    |
        //       | 0     0    0      0   -1/w2   0    |
        //------------------------------------------
        gsl_matrix_complex_set(Hinv, 0, 0, gsl_complex_inverse(gsl_matrix_complex_get(DB, 0, 0)));
        gsl_matrix_complex_set(Hinv, 4, 1, gsl_complex_inverse(gsl_matrix_complex_get(DB, 1, 1)));
        gsl_matrix_complex_set(Hinv, 1, 2, gsl_complex_inverse(gsl_matrix_complex_get(DB, 2, 2)));
        gsl_matrix_complex_set(Hinv, 2, 3, gsl_complex_inverse(gsl_matrix_complex_get(DB, 3, 3)));
        gsl_matrix_complex_set(Hinv, 5, 4, gsl_complex_inverse(gsl_matrix_complex_get(DB, 4, 4)));
        gsl_matrix_complex_set(Hinv, 3, 5, gsl_complex_inverse(gsl_matrix_complex_get(DB, 5, 5)));

        //------------------------------------------
        //Hinv2 = Hinv in matrix<Ofsc> format
        //------------------------------------------
        BUX.zero();
        BUX.set_coef(1.0/gslc_complex(gsl_matrix_complex_get(DB, 0, 0)), 0);
        Hinv2.set_coef(BUX, 0, 0);
        BUX.set_coef(1.0/gslc_complex(gsl_matrix_complex_get(DB, 1, 1)), 0);
        Hinv2.set_coef(BUX, 4, 1);
        BUX.set_coef(1.0/gslc_complex(gsl_matrix_complex_get(DB, 2, 2)), 0);
        Hinv2.set_coef(BUX, 1, 2);
        BUX.set_coef(1.0/gslc_complex(gsl_matrix_complex_get(DB, 3, 3)), 0);
        Hinv2.set_coef(BUX, 2, 3);
        BUX.set_coef(1.0/gslc_complex(gsl_matrix_complex_get(DB, 4, 4)), 0);
        Hinv2.set_coef(BUX, 5, 4);
        BUX.set_coef(1.0/gslc_complex(gsl_matrix_complex_get(DB, 5, 5)), 0);
        Hinv2.set_coef(BUX, 3, 5);

        break;
    }
    case Csts::MAN_CENTER_S:
    {
        //------------------------------------------
        // In the case of the Center-Stable manifold,
        // we need to swap two columns, so that we have:
        //    |  iw1 0    0    0    0   0  |
        //    |  0   0    0    0    0  w2  |
        //    |  0   iw3  0    0    0   0  |
        //H = |  0   0   -iw1  0    0   0  | = (L N)
        //    |  0   0    0    0   -w2  0  |
        //    |  0   0    0   -iw3  0   0  |
        //
        // In this way, we can still use the first 5 columns to
        // go on with the computation.
        //------------------------------------------
        gsl_matrix_complex_set(H, 0, 0, gsl_matrix_complex_get(DB, 0, 0));
        gsl_matrix_complex_set(H, 1, 5, gsl_matrix_complex_get(DB, 1, 1));
        gsl_matrix_complex_set(H, 2, 1, gsl_matrix_complex_get(DB, 2, 2));
        gsl_matrix_complex_set(H, 3, 2, gsl_matrix_complex_get(DB, 3, 3));
        gsl_matrix_complex_set(H, 4, 4, gsl_matrix_complex_get(DB, 4, 4));
        gsl_matrix_complex_set(H, 5, 3, gsl_matrix_complex_get(DB, 5, 5));

        //------------------------------------------
        //H2 = H in matrix<Ofsc> format
        //------------------------------------------
        BUX.zero();
        BUX.set_coef(gslc_complex(gsl_matrix_complex_get(DB, 0, 0)), 0);
        H2.set_coef(BUX, 0, 0);
        BUX.set_coef(gslc_complex(gsl_matrix_complex_get(DB, 1, 1)), 0);
        H2.set_coef(BUX, 1, 5);
        BUX.set_coef(gslc_complex(gsl_matrix_complex_get(DB, 2, 2)), 0);
        H2.set_coef(BUX, 2, 1);
        BUX.set_coef(gslc_complex(gsl_matrix_complex_get(DB, 3, 3)), 0);
        H2.set_coef(BUX, 3, 2);
        BUX.set_coef(gslc_complex(gsl_matrix_complex_get(DB, 4, 4)), 0);
        H2.set_coef(BUX, 4, 4);
        BUX.set_coef(gslc_complex(gsl_matrix_complex_get(DB, 5, 5)), 0);
        H2.set_coef(BUX, 5, 3);

        //------------------------------------------
        //       |1/iw1  0    0      0     0     0    |
        //       | 0     0   1/iw3   0     0     0    |
        //       | 0     0    0    -1/iw1  0     0    |
        //Hinv = | 0     0    0      0     0  -1/iw3  |
        //       | 0     0    0      0    -1/w2  0    |
        //       | 0    1/w2  0      0     0     0    |
        //------------------------------------------
        gsl_matrix_complex_set(Hinv, 0, 0, gsl_complex_inverse(gsl_matrix_complex_get(DB, 0, 0)));
        gsl_matrix_complex_set(Hinv, 5, 1, gsl_complex_inverse(gsl_matrix_complex_get(DB, 1, 1)));
        gsl_matrix_complex_set(Hinv, 1, 2, gsl_complex_inverse(gsl_matrix_complex_get(DB, 2, 2)));
        gsl_matrix_complex_set(Hinv, 2, 3, gsl_complex_inverse(gsl_matrix_complex_get(DB, 3, 3)));
        gsl_matrix_complex_set(Hinv, 4, 4, gsl_complex_inverse(gsl_matrix_complex_get(DB, 4, 4)));
        gsl_matrix_complex_set(Hinv, 3, 5, gsl_complex_inverse(gsl_matrix_complex_get(DB, 5, 5)));

        //------------------------------------------
        //Hinv2 = Hinv in matrix<Ofsc> format
        //------------------------------------------
        BUX.zero();
        BUX.set_coef(1.0/gslc_complex(gsl_matrix_complex_get(DB, 0, 0)), 0);
        Hinv2.set_coef(BUX, 0, 0);
        BUX.set_coef(1.0/gslc_complex(gsl_matrix_complex_get(DB, 1, 1)), 0);
        Hinv2.set_coef(BUX, 5, 1);
        BUX.set_coef(1.0/gslc_complex(gsl_matrix_complex_get(DB, 2, 2)), 0);
        Hinv2.set_coef(BUX, 1, 2);
        BUX.set_coef(1.0/gslc_complex(gsl_matrix_complex_get(DB, 3, 3)), 0);
        Hinv2.set_coef(BUX, 2, 3);
        BUX.set_coef(1.0/gslc_complex(gsl_matrix_complex_get(DB, 4, 4)), 0);
        Hinv2.set_coef(BUX, 4, 4);
        BUX.set_coef(1.0/gslc_complex(gsl_matrix_complex_get(DB, 5, 5)), 0);
        Hinv2.set_coef(BUX, 3, 5);

        break;
    }

    }


    //------------------------------------------
    //                     | LaL  T  |
    //Lambda = Hinv*DB*H = | 0  LaN  | with T = 0 in this context
    //
    //------------------------------------------
    gsl_matrix_complex* AUX = gsl_matrix_complex_calloc(6, 6);
    //AUX = DB*H
    gsl_blas_zgemm (CblasNoTrans , CblasNoTrans , one_c , DB , H , zero_c , AUX );
    //B = Hinv*AUX
    gsl_blas_zgemm (CblasNoTrans , CblasNoTrans , one_c , Hinv , AUX , zero_c , La);
    //------------------------------------------
}

/**
 *  \brief Initialization of the PM objects (W, Wh, and fh). TFS formatting is included.
 **/
void tfts_init_pm(vector<Oftsc>& W,
                  vector<Oftsc>& Wh,
                  vector<Oftsc>& fh,
                  matrix<Ofsc>& PC,
                  vector<Ofsc>& V,
                  gsl_matrix_complex* L,
                  gsl_matrix_complex* La,
                  int man_type)
{
    //------------------------------------------
    // Initialization of the TFC manifold ("W hat")
    //------------------------------------------
    cdouble db;
    gsl_complex dbc;

    //The order 0 is equal to zero, the order 1 is equal to L, which is diagonal
    //----------------------
    //     |   iw1*s1 |
    //     |     0  or w2*s5  |
    //     |   iw3*s2 |
    //Wh = |  -iw1*s3 |
    //     |     0  or -w2*s6 |
    //     |  -iw3*s4 |
    //----------------------
    // Wh[0] = L11*s1 = L[0][0]*s[0]
    dbc = gsl_matrix_complex_get(L, 0, 0);
    db = gslc_complex(dbc);
    Wh[0].set_coef0(db, 1, 0);
    // Wh[2] = L32*s2 = L[2][1]*s[1]
    dbc = gsl_matrix_complex_get(L, 2, 1);
    db = gslc_complex(dbc);
    Wh[2].set_coef0(db, 1, 1);
    // Wh[3] = L43*s3 = L[3][2]*s[2]
    dbc = gsl_matrix_complex_get(L, 3, 2);
    db = gslc_complex(dbc);
    Wh[3].set_coef0(db, 1, 2);
    // Wh[5] = L64*s4 = L[5][3]*s[3]
    dbc = gsl_matrix_complex_get(L, 5, 3);
    db = gslc_complex(dbc);
    Wh[5].set_coef0(db, 1, 3);

    //Special case: Wh[1] and Wh[4]
    //------------------------------------------
    switch(man_type)
    {
    case Csts::MAN_CENTER:
        //Nothing is done, since:
        //----------------------
        //Wh[1] = 0
        //Wh[4] = 0
        break;

    case Csts::MAN_CENTER_S:
        //Only the stable direction is added
        //----------------------
        // Wh[4] = L55*s5 = L[4][4]*s[4]
        dbc = gsl_matrix_complex_get(L, 4, 4);
        db  = gslc_complex(dbc);
        Wh[4].set_coef0(db, 1, 4);
        break;

    case Csts::MAN_CENTER_U:
        //Only the unstable direction is added
        //----------------------
        // Wh[1] = L25*s5 = L[1][4]*s[4]
        dbc = gsl_matrix_complex_get(L, 1, 4);
        db  = gslc_complex(dbc);
        Wh[1].set_coef0(db, 1, 4);
        break;

    case Csts::MAN_CENTER_US:
        //Both hyperbolic directions are added
        //----------------------
        // Wh[1] = L25*s5 = L[1][4]*s[4]
        dbc = gsl_matrix_complex_get(L, 1, 4);
        db  = gslc_complex(dbc);
        Wh[1].set_coef0(db, 1, 4);
        // Wh[4] = L56*s6 = L[4][5]*s[5]
        dbc = gsl_matrix_complex_get(L, 4, 5);
        db  = gslc_complex(dbc);
        Wh[4].set_coef0(db, 1, 5);
        break;
    }


    //------------------------------------------
    // TFS format for Wh
    //------------------------------------------
    tfs_from_ofs_inline(Wh, 0);
    tfs_from_ofs_inline(Wh, 1);

    //------------------------------------------
    // Initialization of the NC manifolds
    //------------------------------------------
    //Applying the COC (with and without zero order) @order 0 and 1
    for(int i = 0; i < 2; i++)
    {
        apply_coc_tfts(PC,  V,  Wh, W,  i, 1);   //order zero is added
    }

    //------------------------------------------
    // Initialization of the reduced VF
    //------------------------------------------
    for(int i = 0; i < REDUCED_NV; i++)
    {
        // fh[i] =  LaL[i][i]*s[i]
        dbc = gsl_matrix_complex_get(La, i, i);
        db = gslc_complex(dbc);
        fh[i].set_coef0(db, 1, i);
    }

    //------------------------------------------
    // TFS format FOR fh
    //------------------------------------------
    tfs_from_ofs_inline(fh, 0);
    tfs_from_ofs_inline(fh, 1);
}

//----------------------------------------------------------------------------------------
//         Solving the homological equations
//----------------------------------------------------------------------------------------
/**
 *  \brief Resolution of the cohomological equations at order m.
 *
 *         The right-hand side of the cohomological equation are in the expansion xi. The expansions eta and fh are updated.
 **/
void solv_hom_eq(vector<Oftsc>& eta,
             vector<Oftsc>& xi,
             vector<Oftsc>& fh,
             gsl_matrix_complex* La,
             int m,
             int param_style,
             double small_div_threshold,
             int** VIs,
             int*  VIn,
             int ims)
{
    //------------------------------------------
    // Initialization
    //------------------------------------------
    cdouble lap;        //contains LaL
    cdouble lar;        //contains LaL*r for each monomial
    cdouble aux;        //temp obj
    cdouble bux;        //temp obj
    cdouble div;
    int kv[REDUCED_NV]; //contains the indices for each monomial

    //------------------------------------------
    // Equation in the tangential space:
    // First d components of eta
    //------------------------------------------
    //Sum on the tangential components of eta (first d components)
    for(int p = 0; p < REDUCED_NV; p++)
    {
        //------------------------------------------
        //Update lap
        //------------------------------------------
        lap = gslc_complex(gsl_matrix_complex_get(La, p, p));

        //------------------------------------------
        //Reset kv
        //------------------------------------------
        kv[0] = m;
        for(int i = 1; i< REDUCED_NV; i++) kv[i] = 0;

        //------------------------------------------
        //Sum on all the coefficients of eta^p_m
        //------------------------------------------
        for(int k = 0; k< Manip::nmon(REDUCED_NV, m); k++)
        {
            //Reset aux
            aux = 0.0+0.0*I;
            //Update aux = laL.kv
            for(int ii = 0; ii < REDUCED_NV; ii++)
            {
                lar = gslc_complex(gsl_matrix_complex_get(La, ii, ii));
                aux += lar*(kv[ii]+0.0*I);
            }

            //Update fh
            for(int j = -OFS_ORDER; j<= OFS_ORDER; j++)
            {
                //divisor = jnI - (lap - laL.kv)
                div = SEML.us.n*I*j - (lap - aux);

                switch(param_style)
                {
                case Csts::GRAPH:
                {
                    //--------------------------------------------------------------------
                    // In GRAPH case: fh is updated, xi is kept
                    // as simple as possible
                    //--------------------------------------------------------------------
                    bux = eta[p].get_coef(m, k)->ofs_get_coef(j);
                    fh[p].get_coef(m, k)->set_coef(0.0*I-bux, j);
                    break;
                }

                case Csts::NORMFORM:
                {
                    //--------------------------------------------------------------------
                    // In NORMFORM case: fh is kept as simple as
                    // possible, unless a small divisor appears
                    //--------------------------------------------------------------------
                    if(cabs(div) < small_div_threshold) //if small divisor, we set the rvf
                    {
                        bux = eta[p].get_coef(m, k)->ofs_get_coef(j);
                        fh[p].get_coef(m, k)->set_coef(0.0*I-bux, j);

                        if(true)
                        {
                            cout <<  setprecision(1);
                            cout << "Small divisor: p = " << p << ", kv = ("<< kv[0] << ", " << kv[1] << ", " << kv[2]  << ", " << kv[3]  << "), j = " << j << "." << endl;
                            cout <<  setprecision(4);
                            cout << "Divisor = " << creal(div) << "  " << cimag(div)  << ", -eta    = " << creal(0.0*I-bux) << "  " << cimag(0.0*I-bux) << endl;
                            cout << "-------------" << endl;
                        }

                    }
                    else //normal form
                    {
                        bux = 0.0*I-eta[p].get_coef(m, k)->ofs_get_coef(j)/div;
                        xi[p].get_coef(m, k)->set_coef(bux, j);
                    }

                    break;
                }

                case Csts::MIXED:
                {
                    //--------------------------------------------------------------------
                    // Isolation of some submanifolds defined
                    // in VIs/VIn
                    // The idea is to simplify as much as possible
                    // the RVF fh for some subsets of indices.
                    //--------------------------------------------------------------------
                    //if i is in VIs[vi] and kv[j] == 0 for all j in VIs[vi]
                    if(is_mixed_style_on(p, kv, VIs, VIn, ims))
                    {
                        //Normal form
                        bux = 0.0*I-eta[p].get_coef(m, k)->ofs_get_coef(j)/div;
                        xi[p].get_coef(m, k)->set_coef(bux, j);

                        //if small divisor, send a warning. Should not happend if the set of indices is good
                        if(cabs(div) < small_div_threshold)
                        {
                            cout <<  setprecision(1);
                            cout << "Small divisor: p = " << p << ", kv = (";
                            for(int l = 0; l < REDUCED_NV-1; l++) cout << kv[l] << ", ";
                            cout << kv[REDUCED_NV-1] << "), j = " << j << "." << endl;
                            cout <<  setprecision(4);
                            cout << "Divisor = " << creal(div) << "  " << cimag(div)  << ", -eta    = " << creal(0.0*I-bux) << "  " << cimag(0.0*I-bux) << endl;
                            cout << "-------------" << endl;
                        }
                    }
                    else
                    {
                        //Graph style
                        bux = eta[p].get_coef(m, k)->ofs_get_coef(j);
                        fh[p].get_coef(m, k)->set_coef(0.0*I-bux, j);
                    }
                    break;
                }

                }
            }
            //Update kv
            if(k < Manip::nmon(REDUCED_NV, m)-1)  Manip::prxkt(kv, REDUCED_NV);
        }
    }

    //------------------------------------------
    // Equation in the normal space:
    // Last n-d components of eta
    //------------------------------------------
    //Sum on the normal components of eta (last n-d components)
    for(int p = REDUCED_NV; p<= Csts::NV-1; p++)
    {
        //Update lap
        lap = gslc_complex(gsl_matrix_complex_get(La, p, p));
        //Reset kv
        kv[0] = m;
        for(int i = 1; i< REDUCED_NV; i++) kv[i] = 0;

        //Sum on all the coefficients of eta^p_m
        for(int k = 0; k< Manip::nmon(REDUCED_NV, m); k++)
        {
            //Reset aux
            aux = 0.0+0.0*I;
            //Update aux = laL.kv
            for(int ii = 0; ii < REDUCED_NV; ii++)
            {
                lar = gslc_complex(gsl_matrix_complex_get(La, ii, ii));
                aux += lar*(kv[ii]+0.0*I);
            }

            //Update xi
            for(int j = -OFS_ORDER; j<= OFS_ORDER; j++)
            {
                //divisor = jnI - (lap - laL.kv)
                div = SEML.us.n*I*j - (lap - aux);

                bux = 0.0*I-eta[p].get_coef(m, k)->ofs_get_coef(j)/div;
                xi[p].get_coef(m, k)->set_coef(bux, j);
            }
            //Update kv
            if(k < Manip::nmon(REDUCED_NV, m)-1)  Manip::prxkt(kv, REDUCED_NV);
        }
    }
}

//----------------------------------------------------------------------------------------
//         Recurrence
//----------------------------------------------------------------------------------------
/**
 *  \brief Adds the contribution of the terms of order m to the order m of Un.
 *         These terms are part of the potential of the primary with position Xe and factor Ke.
 **/
void update_potential_prim_m(vector<Oftsc>& W,        //parameterization in NC coordinates
                              vector<Ofsc>& Xe,        //position of the primary
                              Oftsc& E1, Oftsc& E2,    //intermediate steps
                              Oftsc& Et1, Oftsc& Et2,  //intermediate steps
                              double Ke,               //primary factor
                              vector<Oftsc>& Un,       //potential to update
                              int m)                   //order
{
    //-------------------------
    // Initial operations
    //-------------------------
    //W = W - Xe @ order 0
    W[0].add_coef(-Xe[0], 0, 0);
    W[1].add_coef(-Xe[1], 0, 0);
    //In Et1: (X-Xe)^2+Ye^2+Ze^2 at order m
    Et1.tfts_sprod_zero(W[0], W[0], m);
    Et1.tfts_sprod_zero(W[1], W[1], m);
    Et1.tfts_sprod_zero(W[2], W[2], m);
    //In E1: (X-Xe)^2+Ye^2+Ze^2 at order m
    E1.tfts_sprod_zero(W[0], W[0], m);
    E1.tfts_sprod_zero(W[1], W[1], m);
    E1.tfts_sprod_zero(W[2], W[2], m);
    //In E2: E1^(-3/2): bottleneck
    Et2.tfts_spows_zero(Et1, -1.5, m);
    E2.tfts_spows_zero(Et1, -1.5, m);

    //-------------------------
    // Un
    //-------------------------
    Un[0].tfts_smprod_u_zero(Et2, W[0], Ke+0.0*I, m);
    Un[1].tfts_smprod_u_zero(Et2, W[1], Ke+0.0*I, m);
    Un[2].tfts_smprod_u_zero(Et2, W[2], Ke+0.0*I, m);

    //-------------------------
    // Back to real coordinates
    //-------------------------
    //W = W + Xe @ order 0
    W[0].add_coef(Xe[0], 0, 0);
    W[1].add_coef(Xe[1], 0, 0);
}

/**
 *  \brief Adds the contribution of the terms of order < m to the order m of Un.
 *         These terms are part of the the potential of the primary with position Xe
 *         and factor Ke.
 **/
void update_potential_prim(vector<Oftsc>& W,     //parameterization in NC coordinates
                         vector<Ofsc>& Xe,     //position of the primary
                         Oftsc& E1, Oftsc& E2, //intermediate steps
                         double Ke,            //primary factor
                         vector<Oftsc>& Un,    //potential to update
                         int m)                //order
{
    //-------------------------
    // Initial operations
    //-------------------------
    //W = W - Xe @ order 0
    W[0].add_coef(-Xe[0], 0, 0);
    W[1].add_coef(-Xe[1], 0, 0);
    //In E1: (X-Xe)^2+Ye^2+Ze^2 at order m
    E1.tfts_sprod(W[0], W[0], m);
    E1.tfts_sprod(W[1], W[1], m);
    E1.tfts_sprod(W[2], W[2], m);
    //In E2: E1^(-3/2): bottleneck
    E2.tfts_pows(E1, -1.5, m);

    //-------------------------
    // Un
    //-------------------------
    Un[0].tfts_smprod_u(E2, W[0], Ke+0.0*I, m);
    Un[1].tfts_smprod_u(E2, W[1], Ke+0.0*I, m);
    Un[2].tfts_smprod_u(E2, W[2], Ke+0.0*I, m);

    //-------------------------
    // Back to real coordinates
    //-------------------------
    //W = W + Xe @ order 0
    W[0].add_coef(Xe[0], 0, 0);
    W[1].add_coef(Xe[1], 0, 0);
}

/**
 *  \brief Adds the contribution of the terms of order m to the order m of the potential Un.
 **/
void update_potential_m(vector<Oftsc>& W,        //parameterization in NC coordinates
                         vector<Ofsc>& Xe,        //Earth position
                         vector<Ofsc>& Xm,        //Moon position
                         vector<Ofsc>& Xs,        //Sun position
                         vector<Oftsc>& PrimPt,   //Intermediate steps: contr. of order < m to order m
                         vector<Oftsc>& PrimPt2,  //Intermediate steps: contr. of order = m to order m
                         double Ke,               //Earth factor
                         double Km,               //Moon factor
                         double Ks,               //Sun factor
                         vector<Oftsc>& Un,       //potential to update
                         int m)                   //order
{
    //-------------------------
    // Earth, Moon and Sun potential
    //-------------------------
    if(Ke != 0) update_potential_prim_m(W, Xe, PrimPt[0], PrimPt[1], PrimPt2[0], PrimPt2[1], Ke, Un, m);
    if(Km != 0) update_potential_prim_m(W, Xm, PrimPt[2], PrimPt[3], PrimPt2[2], PrimPt2[3], Km, Un, m);
    if(Ks != 0) update_potential_prim_m(W, Xs, PrimPt[4], PrimPt[5], PrimPt2[4], PrimPt2[5], Ks, Un, m);
}

/**
 *  \brief Adds the contribution of the terms of order < m to the order m of the
 *         potential Un of the three primaries combined.
 **/
void update_potential(vector<Oftsc>& W,        //parameterization in NC coordinates
                    vector<Ofsc>& Xe,        //Earth position
                    vector<Ofsc>& Xm,        //Moon position
                    vector<Ofsc>& Xs,        //Sun position
                    vector<Oftsc>& PrimPt,   //Intermediate steps: contr. of order < m to order m
                    double Ke,               //Earth factor
                    double Km,               //Moon factor
                    double Ks,               //Sun factor
                    vector<Oftsc>& Un,       //potential to update
                    int m)                   //order
{
    //-------------------------
    // Earth, Moon and Sun potential
    //-------------------------
    if(Ke != 0) update_potential_prim(W, Xe, PrimPt[0], PrimPt[1], Ke, Un, m);
    if(Km != 0) update_potential_prim(W, Xm, PrimPt[2], PrimPt[3], Km, Un, m);
    if(Ks != 0) update_potential_prim(W, Xs, PrimPt[4], PrimPt[5], Ks, Un, m);
}

/**
 *  \brief Update the vector field FWc at order m, knowing the parameterization W and the
 *         potential Un of the three primaries combined.
 **/
void apply_vector_field(vector<Ofsc>& alpha,
                        vector<Oftsc>& W,
                        vector<Oftsc>& FWc,
                        vector<Oftsc>& Un,
                        int m)
{
    //-----------------------------------------
    //Configuration variables
    //-----------------------------------------
    //FWc[0] = alpha1*px + alpha2*x + alpha3*y
    FWc[0].tfts_sfsum_tt(W[3], alpha[0], W[0], alpha[1], W[1],  alpha[2], m);
    //FWc[1] = alpha1*py + alpha2*y - alpha3*x
    FWc[1].tfts_sfsum_tt(W[4], alpha[0], W[1], alpha[1], W[0], -alpha[2], m);
    //FWc[2] = alpha1*pz + alpha2*z
    FWc[2].tfts_sfsum_t(W[5], alpha[0], W[2], alpha[1], m);

    //-----------------------------------------
    //Momentum variables
    //-----------------------------------------
    //FWc[3] = -alpha2*px + alpha3*py + alpha24 + order m of the derivative of the potential
    FWc[3].tfts_sfsum_t(W[3], -alpha[1], W[4], alpha[2], m);
    //At order 0, need to add alpha13
    if(m == 0) FWc[3].add_coef(alpha[12], 0, 0);
    //Sum of the Earth+Moon+Sun potential up to degree m+1 (order m)
    FWc[3].tfts_smult_t(Un[0], -alpha[5], m);

    //FWc[4] = -alpha2*py - alpha3*px + alpha26 + order m of the derivative of the potential
    FWc[4].tfts_sfsum_t(W[4], -alpha[1], W[3], -alpha[2], m);
    //At order 0, need to add alpha14
    if(m == 0) FWc[4].add_coef(alpha[13], 0, 0);
    //Sum of the Earth+Moon+Sun potential up to degree m+1 (order m)
    FWc[4].tfts_smult_t(Un[1], -alpha[5], m);

    //FWc[5] = -alpha2*pz  + order m of the derivative of the potential
    FWc[5].tfts_smult_t(W[5], -alpha[1], m);
    //Sum of the Earth+Moon+Sun potential up to degree m+1 (order m)
    FWc[5].tfts_smult_t(Un[2], -alpha[5], m);
}


//----------------------------------------------------------------------------------------
//         Mixed style: indices handling
//----------------------------------------------------------------------------------------
/**
 *  \brief Initialize the arrays that contains the segregated families of solutions
 *         in the mixed style:
 *         For j = 1,..., ims, VIs[j] contains the set of VIn[j] directions that must
 *         be null for a certain submanifold to be invariant.
 **/
void init_mixed_style(int** VIs, int*  VIn, int*  ims, int man_type)
{
    switch(man_type)
    {
    case Csts::MAN_CENTER:
    {
        //--------------------------
        // Init ims
        //--------------------------
        *ims = 2;

        //--------------------------
        // Planar family
        //--------------------------
        VIn[0] = 2;
        VIs[0] = (int*) calloc(2, sizeof(int));
        VIs[0][0] = 2;
        VIs[0][1] = 4;

        //--------------------------
        // Vertical family
        //--------------------------
        VIn[1] = 2;
        VIs[1] = (int*) calloc(2, sizeof(int));
        VIs[1][0] = 1;
        VIs[1][1] = 3;
        break;
    }

    case Csts::MAN_CENTER_S:
    case Csts::MAN_CENTER_U:
    {
        //--------------------------
        // Init ims
        //--------------------------
        *ims = 2; //6

        //--------------------------
        // Hyperbolic normal part
        //--------------------------
        VIn[0] = 4;
        VIs[0] = (int*) calloc(4, sizeof(int));
        VIs[0][0] = 1;
        VIs[0][1] = 2;
        VIs[0][2] = 3;
        VIs[0][3] = 4;

        //--------------------------
        // Center manifold
        //--------------------------
        VIn[1] = 1;
        VIs[1] = (int*) calloc(1, sizeof(int));
        VIs[1][0] = 5;

        //--------------------------
        // Additional submanifolds
        // if uncommented, ims = 6
        //--------------------------
        /*

        //--------------------------
        // Planar family
        //--------------------------
        VIn[2] = 3;
        VIs[2] = (int *) calloc(3, sizeof(int));
        VIs[2][0] = 2;
        VIs[2][1] = 4;
        VIs[2][2] = 5;

        //--------------------------
        // Hyperbolic normal part of the planar family
        //--------------------------
        VIn[3] = 2;
        VIs[3] = (int *) calloc(2, sizeof(int));
        VIs[3][0] = 2;
        VIs[3][1] = 4;

        //--------------------------
        // Vertical family
        //--------------------------
        VIn[4] = 3;
        VIs[4] = (int *) calloc(3, sizeof(int));
        VIs[4][0] = 1;
        VIs[4][1] = 3;
        VIs[4][2] = 5;

        //--------------------------
        // Hyperbolic normal part of the vertical family
        //--------------------------
        VIn[5] = 2;
        VIs[5] = (int *) calloc(2, sizeof(int));
        VIs[5][0] = 1;
        VIs[5][1] = 3;

        */

        break;
    }

    case Csts::MAN_CENTER_US:
    {
        //--------------------------
        // Init ims
        //--------------------------
        *ims = 14; //6

        //--------------------------
        // Unstable manifold
        //--------------------------
        VIn[0] = 5;
        VIs[0] = (int*) calloc(5, sizeof(int));
        VIs[0][0] = 1;
        VIs[0][1] = 2;
        VIs[0][2] = 3;
        VIs[0][3] = 4;
        VIs[0][4] = 6;

        //--------------------------
        // Stable manifold
        //--------------------------
        VIn[1] = 5;
        VIs[1] = (int*) calloc(5, sizeof(int));
        VIs[1][0] = 1;
        VIs[1][1] = 2;
        VIs[1][2] = 3;
        VIs[1][3] = 4;
        VIs[1][4] = 5;

        //--------------------------
        // Hyperbolic normal part (stable + unstable)
        //--------------------------
        VIn[2] = 4;
        VIs[2] = (int*) calloc(4, sizeof(int));
        VIs[2][0] = 1;
        VIs[2][1] = 2;
        VIs[2][2] = 3;
        VIs[2][3] = 4;

        //--------------------------
        // Planar family
        //--------------------------
        VIn[3] = 4;
        VIs[3] = (int*) calloc(4, sizeof(int));
        VIs[3][0] = 2;
        VIs[3][1] = 4;
        VIs[3][2] = 5;
        VIs[3][3] = 6;

        //--------------------------
        // Unstable manifold of the planar family
        //--------------------------
        VIn[4] = 3;
        VIs[4] = (int*) calloc(3, sizeof(int));
        VIs[4][0] = 2;
        VIs[4][1] = 4;
        VIs[4][2] = 6;

        //--------------------------
        // Stable manifold of the planar family
        //--------------------------
        VIn[5] = 3;
        VIs[5] = (int*) calloc(3, sizeof(int));
        VIs[5][0] = 2;
        VIs[5][1] = 4;
        VIs[5][2] = 5;

        //--------------------------
        // Hyperbolic normal part of the planar family
        //--------------------------
        VIn[6] = 2;
        VIs[6] = (int*) calloc(2, sizeof(int));
        VIs[6][0] = 2;
        VIs[6][1] = 4;

        //--------------------------
        // Vertical family
        //--------------------------
        VIn[7] = 4;
        VIs[7] = (int*) calloc(4, sizeof(int));
        VIs[7][0] = 1;
        VIs[7][1] = 3;
        VIs[7][2] = 5;
        VIs[7][3] = 6;

        //--------------------------
        // Unstable manifold of the vertical family
        //--------------------------
        VIn[8] = 3;
        VIs[8] = (int*) calloc(3, sizeof(int));
        VIs[8][0] = 1;
        VIs[8][1] = 3;
        VIs[8][2] = 6;

        //--------------------------
        // Stable manifold of the vertical family
        //--------------------------
        VIn[9] = 3;
        VIs[9] = (int*) calloc(3, sizeof(int));
        VIs[9][0] = 1;
        VIs[9][1] = 3;
        VIs[9][2] = 5;

        //--------------------------
        // Hyperbolic normal part of the vertical family
        //--------------------------
        VIn[10] = 2;
        VIs[10] = (int*) calloc(2, sizeof(int));
        VIs[10][0] = 1;
        VIs[10][1] = 3;

        //--------------------------
        // Center manifold
        //--------------------------
        VIn[11] = 2;
        VIs[11] = (int*) calloc(2, sizeof(int));
        VIs[11][0] = 5;
        VIs[11][1] = 6;

        //--------------------------
        // Center-unstable manifold
        //--------------------------
        VIn[12] = 1;
        VIs[12] = (int*) calloc(1, sizeof(int));
        VIs[12][0] = 6;

        //--------------------------
        // Center-stable manifold
        //--------------------------
        VIn[13] = 1;
        VIs[13] = (int*) calloc(1, sizeof(int));
        VIs[13][0] = 5;

        break;
    }
    }
}

/**
 *  \brief Free the memory of VIs, VIn
 **/
void free_mixed_style(int** VIs, int*  VIn, int ims)
{
    for(int i = 0; i < ims; i++) free(VIs[i]);
    free(VIs);
    free(VIn);
}

/**
 * \brief Inner routine. Checks that the indix \c ind is in the array \c VI, of size \c sVI.
 *
 *  Note that there is a shif of indices: VI contains indices in the form 1, 2, 3, ...,
 *  whereas ind is an indix of the form 0, 1, 2, ...
 *
 **/
bool is_in_VI(int ind, int* VI, int sVI)
{
    for(int i = 0; i < sVI; i++)
    {
        if(VI[i]-1 == ind) return true;
    }
    return false;
}

/**
 * \brief Inner routine. Checks that all indices of the form \c kv[VI[i]-1] are null,
 *        where \c VI is an array of integers of size \c sVI.
 *
 *  Note that there is a shif of indices: VI contains indices in the form 1, 2, 3, ...,
 *  whereas ind is an indix of the form 0, 1, 2, ...
 *
 **/
bool null_ind_VI(int* kv, int* VI, int sVI)
{
    for(int i = 0; i < sVI; i++)
    {
        if(kv[VI[i]-1] != 0) return false;
    }
    return true;
}

/**
 * \brief Inner routine. For all arrays VIs[vi] of size VIn[vi], checks the boolean:
 *                                      is_in_VI(p, VIs[vi], VIn[vi]) && null_ind_VI(kv, VIs[vi], VIn[vi])
 *        This the standard test for mixed style parameterization
 *       (see Haro 2014, section 2.2.3, or BLB PhD manuscript, Section 5.2.2.3).
 **/
bool is_mixed_style_on(int p, int* kv, int** VIs, int* VIn, int numberOfSets)
{
    //Loop on all set of indices
    for(int vi = 0; vi < numberOfSets; vi ++)
    {
        if(is_in_VI(p, VIs[vi], VIn[vi]) && null_ind_VI(kv, VIs[vi], VIn[vi])) return true;
    }

    return false;
}
