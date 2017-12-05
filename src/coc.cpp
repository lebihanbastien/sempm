#include "coc.h"

/**
 * \file coc.cpp
 * \brief Implements the complete change of coordinates (coc) from the
 *        Translated-Floquet-Complexified (TFC) to the Normalized-Centered (NC)
 *        coordinates of the QBCP.
 * \author BLB
 *
 *  This coc is of the form:
 *                              zv_nc = P(t) * C *  zv_tfc(t) + V(t)          (1)
 *  where:
 *
 *
 *  - zv_tfc is the 6-dim state vector in TFC coordinates
 *  - zv_nc is the 6-dim state vector in NC coordinates
 *  - P(t) is a 6 x 6 periodic matrix, whose coefficients are Fourier series
 *  - V(t) is a 6 x 1 periodic vector, whose coefficients are Fourier series
 *  - C is a constant 6 x 6 complex matrix.
 *
 *  In practice, we define and use the matrix PC(t) = P(t)*C.
 *
 *  The Fourier series are manipulated as Ofsc == Ofs <double complex> objects, either in
 *  time (tfs) or frequency (ofs) format (see ofs.h and ofs.tpp for details).
 *
 *
 *
 *  The inverse of (1) is denoted:
 *                                zv_tfc = CQ(t) * (zv_nc(t) - V(t)           (2)
 *
 *  where CQ(t) = inv(PC(t)).
 *
 *  The operation (1) or its inverse is performed with zv_tfc and z as
 *      - vector<Oftsc = Otfs < Ofsc > > objects (i.e. vector of Fourier-Taylor series) or
 *      - vector< Ofsc > objects (i.e. vector of Fourier series).
 *
 */


//----------------------------------------------------------------------------------------
//
//          OFTS version of the coc
//
//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------
//          Initialization
//----------------------------------------------------------------------------------------
/**
 *  \brief Retrieves the Fouriers series coefficients of the matrices P and Q, the vectors
 *         V, Xe, Xm and Xs and the scalars ILe, ILs and ILm.
 *
 *  - We recall that the change of coordinates between the NC coordinates and the TFC
 *  coordinates is of the form:
 *                              zv_nc = P(t) * C *  zv_tfc(t) + V(t)          (1)
 *  where:
 *
 *
 *  - zv_tfc is the 6-dim state vector in TFC coordinates
 *  - zv_nc is the 6-dim state vector in NC coordinates
 *  - P(t) is a 6 x 6 periodic matrix, whose coefficients are Fourier series
 *  - V(t) is a 6 x 1 periodic vector, whose coefficients are Fourier series
 *  - C is a constant 6 x 6 complex matrix.
 *
 *  In practice, we define and use the matrix PC(t) = P(t)*C.
 *  The inverse of (1) is denoted:
 *                                zv_tfc = CQ(t) * (zv_nc(t) - V(t)           (2)
 *
 *  In this routine,
 *
 *  1. First step
 *      * P, Q and V are retrieved from txt files in the folder \c fbpl.F_COC.
 *      * dot(V) is computed and stored in Vdot.
 *      * The matrices PC = P(t)*C and CQ = inv(PC) are computed and stored in PC and CQ.
 *      * dot(PC) is computed and stored in PCdot.
 *
 *  2. Second step
 *     * The vectors Xe[0:2], Xm[0:2], and Xs[0:2] contain the true
 *       (non-shifted) time-dependent positions of the primaries in the xy-plane. They are
 *       retrieved from txt files in the folder \c fbpl.F_COC.
 *     * The Fourier series ILc = 1/sqrt(Xc[0]^2 + Xc[1]^2) are retrieved from txt files,
 *       for c = e, m, s.
 *
 *
 *  Note that, in the CRTBP case, the Fourier series become constant variables, and are
 *  computed from CRTBP constants rather than retrieved from txt files.
 *
 **/
void init_coc_ofs(matrix<Ofsc>& P, matrix<Ofsc>& Q, matrix<Ofsc>& PC, matrix<Ofsc>& PCdot,
                  matrix<Ofsc>& CQ, vector<Ofsc>& Xe, vector<Ofsc>& Xm, vector<Ofsc>& Xs,
                  vector<Ofsc>& V, vector<Ofsc>& Vdot, Ofsc& ILe, Ofsc& ILm, Ofsc& ILs,
                  FBPL& fbpl)
{
    //------------------------------------------------------------------------------------
    //Retrieve folder
    //------------------------------------------------------------------------------------
    string F_COC = fbpl.cs.F_COC;

    //------------------------------------------------------------------------------------
    //Switch CRTBP/QBCP
    //------------------------------------------------------------------------------------
    if(fbpl.model == Csts::QBCP || fbpl.model == Csts::BCP)
    {
        //Recovering the data: matrix P
        for(int i = 0; i < Csts::NV ; i++)
            for(int j = 0; j < Csts::NV ; j++)
                read_coc(*P.get_ptr_first_coef(i,j), F_COC+"P",  i+1, j+1);
        //Recovering the data: matrix Q
        for(int i = 0; i < Csts::NV ; i++)
            for(int j = 0; j < Csts::NV ; j++)
                read_coc(*Q.get_ptr_first_coef(i,j), F_COC+"Q",  i+1, j+1);

        //Recovering the data: vector V
        //V = [G1_11 G1_12 0 G1_21 G1_22 0]^T
        //Note:V[2] and V[5] are kept null
        read_coc(V[0], F_COC+"G1_",  1, 1);
        read_coc(V[1], F_COC+"G1_",  1, 2);
        read_coc(V[3], F_COC+"G1_",  2, 1);
        read_coc(V[4], F_COC+"G1_",  2, 2);

        //Recovering the data: vector Xm
        read_coc(Xm[0], F_COC+"Xm",  4, 1);
        read_coc(Xm[1], F_COC+"Xm",  5, 1);
        //Note: Xm[2] is kept null

        //Recovering the data: vector Xe
        read_coc(Xe[0], F_COC+"Xe",  4, 1);
        read_coc(Xe[1], F_COC+"Xe",  5, 1);
        //Note: Xe[2] is kept null

        //Recovering the data: vector Xs
        read_coc(Xs[0], F_COC+"Xs",  4, 1);
        read_coc(Xs[1], F_COC+"Xs",  5, 1);
        //Note: Xs[2] is kept null


        //Recovering the data: scalars ILe, ILm, ILs
        read_coc(ILe, F_COC+"Xe",  6, 1);
        read_coc(ILm, F_COC+"Xm",  6, 1);
        read_coc(ILs, F_COC+"Xs",  6, 1);

    }
    else //RTPB
    {
        //note that V is left untouched (set to zero by default)
        //--------------------------------------------------------------------------------

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
        //Init P
        //--------------------------------------------------------------------------------
        P.set_coef(+2*la1/s1,                          0, 1);
        P.set_coef(+(la1*la1  - 2*c2 - 1)/s1,          1, 1);
        P.set_coef(+(la1*la1  + 2*c2 + 1)/s1,          3, 1);
        P.set_coef(+(la1*la1*la1 + (1 - 2*c2)*la1)/s1, 4, 1);

        P.set_coef(-(om1*om1 + 2*c2 + 1)/s2,           1, 0);
        P.set_coef(-(om1*om1 - 2*c2 - 1)/s2,           3, 0);

        P.set_coef(-2*la1/s1,                          0, 4);
        P.set_coef(+(la1*la1  - 2*c2 - 1)/s1,          1, 4);
        P.set_coef(+(la1*la1  + 2*c2 + 1)/s1,          3, 4);
        P.set_coef(-(la1*la1*la1 + (1 - 2*c2)*la1)/s1, 4, 4);

        P.set_coef(+2*om1/s2,                          0, 3);
        P.set_coef(-(om1*om1*om1 - (1 - 2*c2)*om1)/s2, 4, 3);

        P.set_coef(+1.0/sqrt(om2),                     2, 2);
        P.set_coef(+sqrt(om2),                         5, 5);

        //--------------------------------------------------------------------------------
        //Check symplectic nature of P (uncomment if desired)
        //--------------------------------------------------------------------------------
        /*
        gsl_matrix *Ps = gsl_matrix_calloc(6,6);
        for(int i = 0; i <6; i++)
        {
            for(int j = 0; j <6; j++)
            {
                gsl_matrix_set(Ps, i, j, creal(P.get_coef(i,j).ofs_get_coef(0)));
            }
        }
        symplecticity_test_real(Ps, Csts::INVERSE_GSL);
        gsl_matrix_free(Ps);
        //*/

        //--------------------------------------------------------------------------------
        //Q = inv(P) (GSL is used)
        //--------------------------------------------------------------------------------
        int s;
        gsl_matrix* Pc   = gsl_matrix_calloc (Csts::NV, Csts::NV);
        gsl_matrix* Qc   = gsl_matrix_calloc (Csts::NV, Csts::NV);
        gsl_permutation* p6 = gsl_permutation_alloc (Csts::NV);

        //Init Pc
        for(int i =0; i < Csts::NV; i++) for(int j =0; j < Csts::NV; j++) gsl_matrix_set(Pc, i, j, creal(P.get_ptr_first_coef(i,j)->ofs_get_coef(0)));
        //Use of GSL library
        gsl_linalg_LU_decomp (Pc, p6, &s);
        gsl_linalg_LU_invert (Pc, p6, Qc);


        //--------------------------------------------------------------------------------
        // Init Q
        //--------------------------------------------------------------------------------
        for(int i =0; i < Csts::NV; i++) for(int j =0; j < Csts::NV; j++) Q.set_coef(gsl_matrix_get(Qc, i, j), i, j);


        //--------------------------------------------------------------------------------
        // Free GSL objects
        //--------------------------------------------------------------------------------
        gsl_matrix_free(Pc);
        gsl_matrix_free(Qc);
        gsl_permutation_free(p6);


        //--------------------------------------------------------------------------------
        // Xe, Xm, Xs, ILe, ILm, ILs
        //--------------------------------------------------------------------------------
        double gamma = fbpl.cs.gamma;
        switch(fbpl.coord_sys)
        {
        case Csts::EM:
            //Sun does not exist
            Xs[0].set_coef(  0.0+0.0*I, 0);
            Xs[1].set_coef(  0.0+0.0*I, 0);
            ILs.set_coef(    0.0+0.0*I, 0);

            //Earth and Moon
            if(fbpl.li_EM == 1)
            {
                Xe[0].set_coef((-1.0+gamma)/gamma+0.0*I, 0);
                Xe[1].set_coef( 0.0+0.0*I, 0);
                ILe.set_coef( 1.0/sqrt( pow((-1.0+gamma)/gamma, 2.0))+0.0*I, 0);

                Xm[0].set_coef( +1.0+0.0*I, 0);
                Xm[1].set_coef(  0.0+0.0*I, 0);
                ILm.set_coef(  1.0+0.0*I, 0);
            }
            else
            {
                Xe[0].set_coef((-1.0-gamma)/gamma+0.0*I, 0);
                Xe[1].set_coef( 0.0+0.0*I, 0);
                ILe.set_coef( gamma/(1.0+gamma)+0.0*I, 0);

                Xm[0].set_coef( -1.0+0.0*I, 0);
                Xm[1].set_coef(  0.0+0.0*I, 0);
                ILm.set_coef(    1.0+0.0*I, 0);
            }
            break;

        case Csts::SEM:
            //Moon does not exist
            Xm[0].set_coef(  0.0+0.0*I, 0);
            Xm[1].set_coef(  0.0+0.0*I, 0);
            ILm.set_coef(    0.0+0.0*I, 0);

            //Sun and Moon
            if(fbpl.li_SEM == 1)
            {
                Xs[0].set_coef((-1.0+gamma)/gamma+0.0*I, 0);
                Xs[1].set_coef( 0.0+0.0*I, 0);
                ILs.set_coef( 1.0/sqrt( pow((-1.0+gamma)/gamma, 2.0))+0.0*I, 0);

                Xe[0].set_coef( +1.0+0.0*I, 0);
                Xe[1].set_coef(  0.0+0.0*I, 0);
                ILe.set_coef(  1.0+0.0*I, 0);
            }
            else
            {
                Xs[0].set_coef((-1.0-gamma)/gamma+0.0*I, 0);
                Xs[1].set_coef( 0.0+0.0*I, 0);
                ILs.set_coef( gamma/(1.0+gamma)+0.0*I, 0);

                Xe[0].set_coef( -1.0+0.0*I, 0);
                Xe[1].set_coef(  0.0+0.0*I, 0);
                ILe.set_coef(    1.0+0.0*I, 0);
            }
            break;
        }
    }

    //------------------------------------------------------------------------------------
    // Building PC, CQ
    //------------------------------------------------------------------------------------
    Ofsc BUX(OFS_ORDER);
    //Keymap used to initialize PC and CQ
    vector<int> keyMap(6);
    keyMap[0] = 0;
    keyMap[1] = 1;
    keyMap[2] = 3;
    keyMap[3] = 4;
    keyMap[4] = 2;
    keyMap[5] = 5;

    int ii;
    //Init PC by rows
    for(int i = 0; i <= 3; i++)
    {
        ii = keyMap[i];
        BUX.ofs_fsum(P(ii,0),   1.0/sqrt(2)+0.0*I, P(ii,3), I*1.0/sqrt(2));
        PC.set_coef(BUX, ii, 0);
        BUX.ofs_fsum(P(ii,0), I*1.0/sqrt(2), P(ii,3),   1.0/sqrt(2)+0.0*I);
        PC.set_coef(BUX, ii, 3);
        PC.set_coef(P(ii,1), ii, 1);
        PC.set_coef(P(ii,4), ii, 4);
    }

    for(int i = 4; i <= 5; i++)
    {
        ii = keyMap[i];
        BUX.ofs_fsum(P(ii,2),   1.0/sqrt(2)+0.0*I, P(ii,5), I*1.0/sqrt(2));
        PC.set_coef(BUX, ii, 2);
        BUX.ofs_fsum(P(ii,2), I*1.0/sqrt(2), P(ii,5),   1.0/sqrt(2)+0.0*I);
        PC.set_coef(BUX, ii, 5);
    }

    //Init CQ by columns
    for(int i = 0; i <= 3; i++)
    {
        ii = keyMap[i];
        BUX.ofs_fsum(Q(0,ii),    1.0/sqrt(2)+0.0*I, Q(3,ii), -1.0/sqrt(2)*I);
        CQ.set_coef(BUX, 0, ii);
        BUX.ofs_fsum(Q(0,ii), -1.0/sqrt(2)*I, Q(3,ii),    1.0/sqrt(2)+0.0*I);
        CQ.set_coef(BUX, 3, ii);
        CQ.set_coef(Q(1,ii), 1, ii);
        CQ.set_coef(Q(4,ii), 4, ii);
    }

    for(int i = 4; i <= 5; i++)
    {
        ii = keyMap[i];
        BUX.ofs_fsum(Q(2,ii),    1.0/sqrt(2)+0.0*I, Q(5,ii), -1.0/sqrt(2)*I);
        CQ.set_coef(BUX, 2, ii);
        BUX.ofs_fsum(Q(2,ii), -1.0/sqrt(2)*I, Q(5,ii),    1.0/sqrt(2)+0.0*I);
        CQ.set_coef(BUX, 5, ii);
    }

    //------------------------------------------------------------------------------------
    // Building Vdot
    //------------------------------------------------------------------------------------
    for(int i = 0; i < Csts::NV; i++)
    {

        //------------------------------------------------------------------------------------
        Vdot[i].dot(V[i], fbpl.us.n);
    }
    //Init PCdot
    //------------------------------------------------------------------------------------
    PCdot.dot(PC, fbpl.us.n);
}


//----------------------------------------------------------------------------------------
//          Applying
//----------------------------------------------------------------------------------------
/**
 *  \brief Apply the direct change of variables (1) at every order.
 *         The Fourier series are supposed to be in ofs format.
 **/
void apply_coc_ofts(matrix<Ofsc>& PC,
                    vector<Ofsc>& V,
                    vector<Oftsc>& zv_tfc_in,
                    vector<Oftsc>& zv_nc_out)
{
    //zeroing the target
    for(unsigned int i = 0; i < zv_nc_out.size(); i++) zv_nc_out[i].zero();
    //zv_nc_out = PC*zv_tfc_in
    smvprod_u(PC, zv_tfc_in, zv_nc_out);
    //zv_nc_out+=V(theta)
    add_coef(V, zv_nc_out);
}


/**
 *  \brief Apply the change of variables (1) at order m, with or without the zero order
 *         V(t), depending on the boolean is_zero_order_shift.
 *         The Fourier series are supposed to be in ofs format.
 **/
void apply_coc_ofts(matrix<Ofsc>& PC,
                    vector<Ofsc>& V,
                    vector<Oftsc>& zv_tfc_in,
                    vector<Oftsc>& zv_nc_out,
                    int m,
                    int is_zero_order_shift)
{
    //zv_nc_out = PC*zv_tfc_in
    smvprod_u(PC, zv_tfc_in, zv_nc_out, m);
    //zv_nc_out += V(t)
    if(m == 0 && is_zero_order_shift == 1) add_coef(V, zv_nc_out);
}


/**
 *  \brief Apply the change of variables (1) at order m, with or without the zero order
 *         V(t), depending on the boolean is_zero_order_shift.
 *         The Fourier series are supposed to be in tfs format.
 **/
void apply_coc_tfts(matrix<Ofsc>& PC,
                    vector<Ofsc>& V,
                    vector<Oftsc>& zv_tfc_in,
                    vector<Oftsc>& zv_nc_out,
                    int m,
                    int flag)
{
    //zv_nc_out = PC*zv_tfc_in
    tfts_smvprod_u(PC, zv_tfc_in, zv_nc_out, m);

    //zv_nc_out+=V(theta)
    if(m == 0 && flag == 1) tfts_add_coef(V, zv_nc_out);
}


/**
 *  \brief Apply the inverse (2) of the change of variables (1) at every order.
 *         The Fourier series are supposed to be in ofs format.
 **/
void apply_inc_coc_ofts(matrix<Ofsc>& CQ,
                        vector<Ofsc>& V,
                        vector<Oftsc>& zv_nc_in,
                        vector<Oftsc>& zv_tfc_out,
                        vector<Oftsc>& zv_tfr_tp)
{
    //zv_tfr_tp = zv_nc_in - V
    sub_coef(V, zv_nc_in, zv_tfr_tp);

    //zv_tfc_out = CQ*zv_tfr_tp
    smvprod_u(CQ, zv_tfr_tp, zv_tfc_out);
}

/**
 *  \brief Apply the derivative of the direct change of variables (1) at every order.
 *         The Fourier series are supposed to be in ofs format.
 *
 *         The change of variables is of the form:
 *
 *                  zv_nc_dot = PCdot * zv_tfc + PC * zv_tfc_dot + Vdot
 **/
void apply_coc_dot_ofts(matrix<Ofsc>& PC,
                        matrix<Ofsc>& PCdot,
                        vector<Ofsc>& Vdot,
                        vector<Oftsc>& zv_tfc,
                        vector<Oftsc>& zv_tfc_dot,
                        vector<Oftsc>& zv_nc_dot)
{
    //Zeroing the target
    for(unsigned int i = 0; i < zv_nc_dot.size(); i++) zv_nc_dot[i].zero();
    //zv_nc_dot += PCdot*zv_tfc
    smvprod_u(PCdot, zv_tfc, zv_nc_dot);
    //zv_nc_dot += PC*zv_tfc_dot
    smvprod_u(PC, zv_tfc_dot, zv_nc_dot);
    //zv_nc_dot +=Vdot(theta)
    add_coef(Vdot, zv_nc_dot);
}


/**
 *  \brief Apply the derivative of the direct change of variables (1) at order m.
 *         The Fourier series are supposed to be in ofs format.
 *
 *         The change of variables is of the form:
 *
 *                  zv_nc_dot = PCdot * zv_tfc + PC * zv_tfc_dot + Vdot
 **/
void apply_coc_dot_ofts(matrix<Ofsc>& PC,
                        matrix<Ofsc>& PCdot,
                        vector<Ofsc>& Vdot,
                        vector<Oftsc>& zv_tfc,
                        vector<Oftsc>& zv_tfc_dot,
                        vector<Oftsc>& zv_nc_dot,
                        int m)
{
    //zv_nc_dot += PCdot*zv_tfc
    smvprod_u(PCdot, zv_tfc, zv_nc_dot, m);
    //zv_nc_dot += PC*zv_tfc_dot
    smvprod_u(PC, zv_tfc_dot, zv_nc_dot, m);
    //zv_nc_dot +=Vdot(theta)
    if(m == 0) add_coef(Vdot, zv_nc_dot);
}


/**
 *  \brief Apply the derivative of the inverse of the change of variables at every orders.
 *         The Fourier series are supposed to be in ofs format.
 *         The change of variables is of the form:
 *
 *                  zv_tfc_dot = CQ * (zv_nc_dot - PCdot * zv_tfc - Vdot)
 **/
void apply_inv_coc_dot_ofts(matrix<Ofsc>&  CQ,
                            matrix<Ofsc>&  PCdot,
                            vector<Ofsc>&  Vdot,
                            vector<Oftsc>& zv_tfc,
                            vector<Oftsc>& zv_nc_dot,
                            vector<Oftsc>& zv_tfc_dot,
                            vector<Oftsc>& ztp_1,
                            vector<Oftsc>& ztp_2)
{
    //Zeroing the target
    for(unsigned int i = 0; i < zv_tfc_dot.size(); i++) zv_tfc_dot[i].zero();
    for(unsigned int i = 0; i < zv_tfc_dot.size(); i++) ztp_1[i].zero();
    for(unsigned int i = 0; i < zv_tfc_dot.size(); i++) ztp_2[i].zero();

    //ztp_1 += PCdot*zv_tfc @order m
    smvprod_u(PCdot, zv_tfc, ztp_1);
    //ztp_2 += zv_nc_dot - PCdot*zv_tfc - Vdot @order m
    for(unsigned int i = 0; i < zv_tfc.size(); i++)
    {
        ztp_2[i].ofts_smult_u(zv_nc_dot[i],  +1.0+0.0*I);
        ztp_2[i].ofts_smult_u(ztp_1[i], -1.0+0.0*I);
    }
    sub_coef(Vdot, ztp_2);
    //zv_tfc_dot = CQ*ztp_2 = CQ*(zv_nc_dot - PCdot*zv_tfc - Vdot) @order m
    smvprod_u(CQ, ztp_2, zv_tfc_dot);

}


/**
 *  \brief Apply the derivative of the inverse of the change of variables at order m.
 *         The Fourier series are supposed to be in ofs format.
 *         The change of variables is of the form:
 *
 *                  zv_tfc_dot = CQ * (zv_nc_dot - PCdot * zv_tfc - Vdot)
 **/
void apply_inv_coc_dot_ofts(matrix<Ofsc>& CQ,
                            matrix<Ofsc>& PCdot,
                            vector<Ofsc>& Vdot,
                            vector<Oftsc>& zv_tfc,
                            vector<Oftsc>& zv_nc_dot,
                            vector<Oftsc>& zv_tfc_dot,
                            vector<Oftsc>& ztp_1,
                            vector<Oftsc>& ztp_2,
                            int m)
{
    //ztp_1 += PCdot*zv_tfc @order m
    smvprod_u(PCdot, zv_tfc, ztp_1, m);
    //ztp_2 += zv_nc_dot - PCdot*zv_tfc - Vdot @order m
    for(unsigned int i = 0; i < zv_tfc.size(); i++)
    {
        ztp_2[i].ofts_smult_u(zv_nc_dot[i],  +1.0+0.0*I, m);
        ztp_2[i].ofts_smult_u(ztp_1[i], -1.0+0.0*I, m);
    }
    if(m == 0) sub_coef(Vdot, ztp_2);
    //zv_tfc_dot = CQ*ztp_2 = CQ*(zv_nc_dot - PCdot*zv_tfc - Vdot) @order m
    smvprod_u(CQ, ztp_2, zv_tfc_dot, m);
}


/**
 *  \brief Apply the derivative of the inverse of the change of variables at order m.
 *         The Fourier series are supposed to be in tfs format.
 *         The change of variables is of the form:
 *
 *                  zv_tfc_dot = CQ * (zv_nc_dot - PCdot * zv_tfc - Vdot)
 **/
void apply_inv_coc_dot_tfts(matrix<Ofsc>&  CQ,
                            matrix<Ofsc>&  PCdot,
                            vector<Ofsc>&  Vdot,
                            vector<Oftsc>& zv_tfc,
                            vector<Oftsc>& zv_nc_dot,
                            vector<Oftsc>& zv_tfc_dot,
                            vector<Oftsc>& ztp_1,
                            vector<Oftsc>& ztp_2,
                            int m)
{
    //ztp_1 += PCdot*zv_tfc @order m
    tfts_smvprod_u(PCdot, zv_tfc, ztp_1, m);
    //ztp_2 += zv_nc_dot - PCdot*zv_tfc - Vdot @order m
    for(unsigned int i = 0; i < zv_tfc.size(); i++)
    {
        ztp_2[i].tfts_smult_u(zv_nc_dot[i],  +1.0+0.0*I, m);
        ztp_2[i].tfts_smult_u(ztp_1[i], -1.0+0.0*I, m);
    }
    if(m == 0) tfts_subCoef(Vdot, ztp_2);
    //zv_tfc_dot = CQ*ztp_2 = CQ*(zv_nc_dot - PCdot*zv_tfc - Vdot) @order m
    tfts_smvprod_u(CQ, ztp_2, zv_tfc_dot, m);
}

//----------------------------------------------------------------------------------------
//
//          Switch tfs/ofs format
//
//----------------------------------------------------------------------------------------
/**
 *  \brief From ofs to tfs format for the whole coc.
 **/
void tfs_from_ofs(matrix<Ofsc>& P,
                  matrix<Ofsc>& Q,
                  matrix<Ofsc>& PC,
                  matrix<Ofsc>& PCdot,
                  matrix<Ofsc>& CQ,
                  vector<Ofsc>& Xe,
                  vector<Ofsc>& Xm,
                  vector<Ofsc>& Xs,
                  vector<Ofsc>& V,
                  vector<Ofsc>& Vdot,
                  Ofsc& ILe,
                  Ofsc& ILm,
                  Ofsc& ILs)
{
    tfs_from_ofs_inline(P);
    tfs_from_ofs_inline(Q);
    tfs_from_ofs_inline(PC);
    tfs_from_ofs_inline(PCdot);
    tfs_from_ofs_inline(CQ);
    tfs_from_ofs_inline(V);
    tfs_from_ofs_inline(Xe);
    tfs_from_ofs_inline(Xm);
    tfs_from_ofs_inline(Xs);
    tfs_from_ofs_inline(Vdot);

    Ofsc temp(ILe.get_order());
    ILe.tfs_from_ofs_inline(temp);
    ILm.tfs_from_ofs_inline(temp);
    ILs.tfs_from_ofs_inline(temp);
}



//----------------------------------------------------------------------------------------
//
//          OFS version of the coc
//
//----------------------------------------------------------------------------------------
/**
 *  \brief Apply the direct change of variables (1) for zv_tfc_in/zv_nc_out in
 *         vector<Ofsc> format.
 *         The Fourier series are supposed to be in ofs format.
 **/
void apply_coc_ofs(matrix<Ofsc>& PC,
                   vector<Ofsc>& V,
                   vector<Ofsc>& zv_tfc_in,
                   vector<Ofsc>& zv_nc_out)
{
    //zeroing the target
    for(unsigned int i = 0; i < zv_nc_out.size(); i++) zv_nc_out[i].zero();
    //zv_nc_out = PC*zv_tfc_in
    smvprod_ofs(PC, zv_tfc_in, zv_nc_out);
    //zv_nc_out+=V(theta)
    for(int i = 0; i < (int) zv_nc_out.size(); i++) zv_nc_out[i] += V[i];
}

