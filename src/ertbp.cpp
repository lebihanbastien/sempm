#include "ertbp.h"

/**
 * \file ertbp.cpp
 * \brief Routines of the Elliptic Restricted Three-Body Problem (ERTBP) (src)
 * \author BLB.
 * \date November 2015
 * \version 1.0
 */

//-----------------------------------------------------------------------------
// Secondary routines
//-----------------------------------------------------------------------------
/**
* \brief Using the Newton-Raphson method, find the root of a function known to lie close to x1 The root
*  /rtnewt will be refined until its accuracy is known within ± xacc.
*  funcd is a user-supplied routine that returns both the function value and the first derivative of the
*  function at the point x.
**/
double rtnewtecc(void (*funcd)(double, double, double, double *, double *), double x1, double xacc, double e, double  M)
{
    void nrerror(char error_text[]);
    int j;
    double df,dx,f,rtn;

    rtn=x1;   //initial guess

    for (j=1; j<=50; j++)
    {
        (*funcd)(e, M, rtn,&f,&df);
        dx=f/df;
        rtn -= dx;

        //printf("rtnewtecc. %d: %5.15f\n", j, rtn);
        if (fabs(dx) < xacc) return rtn;  //Convergence
    }
    printf("WARNING: Maximum number of iterations exceeded in rtnewt");
    return rtn;   //Never get here.
}

/**
 * \brief Provides the function value and its first derivative for the newton-raphson method.
 *        The evaluated function is the Kepler equations: M = E - e*sin(E).
 **/
void eccDer(double e, double M, double E, double *f, double *df)
{
    *f  = E - e*sin(E) - M;
    *df = 1.0 -  e*cos(E);
}

/**
 *   \brief Get the FFT of the data stored in dEv. Note that, ideally:
 *          - xFFT is of order fftN, with 2*N+1 coefficients.
 *          - dEv contains 2*fftN+1 coefficients.
 **/
void eccFFT(Ofsc &xFFT, int fftN, int N, gsl_vector *dEv, int parity)
{
    if(N < 2*fftN+1)
    {
        cout << "eccFFT. Error: we must have N >= 2*fftN+1" << endl;
        return;
    }

    //--------------------------------------------------------------
    //FFT tools
    //--------------------------------------------------------------
    gsl_vector_complex *data_complex = gsl_vector_complex_calloc(N);
    gsl_vector *data = gsl_vector_calloc(N);
    gsl_fft_real_wavetable * wavetable = gsl_fft_real_wavetable_alloc (N);
    gsl_fft_real_workspace * workspace = gsl_fft_real_workspace_alloc (N);

    //--------------------------------------------------------------
    //Set data
    //--------------------------------------------------------------
    for(int i = 0; i< N; i++) gsl_vector_set(data, i, gsl_vector_get(dEv, i));

    //--------------------------------------------------------------
    //FFT transform
    //--------------------------------------------------------------
    gsl_fft_real_transform (data->data, 1, N, wavetable, workspace);
    gsl_fft_halfcomplex_unpack(data->data , data_complex->data ,  data->stride ,data->size);

    //--------------------------------------------------------------
    //Order 0
    //--------------------------------------------------------------
    switch(parity)
    {
    case 1: //even case (cosinus)
        xFFT.setCoef(+GSL_REAL(gsl_vector_complex_get(data_complex, 0))/(double)N,  0);
        break;
    case 0: //odd case (sinus)
        xFFT.setCoef(I*GSL_IMAG(gsl_vector_complex_get(data_complex, 0))/(double)N,  0);
        break;
    default: //unknown
        xFFT.setCoef(+GSL_REAL(gsl_vector_complex_get(data_complex, 0))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_complex, 0))/(double)N,  0);
        break;
    }

    //--------------------------------------------------------------
    //Order >< 0
    //--------------------------------------------------------------
    for(int i = 1; i<= fftN; i++)
    {
        switch(parity)
        {
        case 1: //even case (cosinus)
        {
            //Negative frequencies
            xFFT.setCoef(+GSL_REAL(gsl_vector_complex_get(data_complex, N-i))/(double)N, -i);
            //Positive frequencies
            xFFT.setCoef(+GSL_REAL(gsl_vector_complex_get(data_complex, i))/(double)N,  i);
            break;
        }
        case 0: //odd case (sinus)
        {
            //Negative frequencies
            xFFT.setCoef(I*GSL_IMAG(gsl_vector_complex_get(data_complex, N-i))/(double)N, -i);
            //Positive frequencies
            xFFT.setCoef(I*GSL_IMAG(gsl_vector_complex_get(data_complex, i))/(double)N,  i);
            break;
        }
        default: //unknown
        {
            //Negative frequecies
            xFFT.setCoef(+GSL_REAL(gsl_vector_complex_get(data_complex, N-i))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_complex, N-i))/(double)N, -i);
            //Positive frequencies
            xFFT.setCoef(+GSL_REAL(gsl_vector_complex_get(data_complex, i))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_complex, i))/(double)N,  i);
        }

        }
    }


    gsl_vector_complex_free(data_complex);
    gsl_vector_free(data);
    gsl_fft_real_wavetable_free(wavetable);
    gsl_fft_real_workspace_free(workspace);

}

/**
 *   \brief Get Fourier decomposition of the eccentric anomaly via Bessel functions
 **/
void eccBessel(Ofsc &xFFT, int fftN, double e)
{
    cdouble ek;
    gsl_sf_result bess;
    for(int k = 1; k < fftN; k++)
    {
        gsl_sf_bessel_Jn_e(k, k*e, &bess);
        ek = -1.0*I*bess.val/(double) k;
        xFFT.setCoef(ek, k);
        xFFT.setCoef(-1.0*ek, -k);
    }
}

//-----------------------------------------------------------------------------
// Main routine: ertbp
//-----------------------------------------------------------------------------
/**
 *  \brief Main routine to compute the Elliptic Three-Body Problem in Ofs format.
 */
void ertbp(int li_EM, int li_SEM, int coordsys)
{
    cout << "---------------------------------------------------" << endl;
    cout << "                                                   " << endl;
    cout << "              Storage of the                       " << endl;
    cout << "        Elliptic Three-Body Problem                " << endl;
    cout << "                                                   " << endl;
    cout << "---------------------------------------------------" << endl;
    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(15);

    //--------------------------------------------------------------------
    // 1. Init
    //--------------------------------------------------------------------
    //Init the ERTBP
    FBP fbp;
    init_FBP(&fbp, Csts::SUN, Csts::EARTH, Csts::MOON);

    //Init the ERTBP focused on one libration point
    FBPL fbpl;
    init_FBPL(&fbpl, &fbp, li_EM, li_SEM, Csts::ERTBP, coordsys, Csts::GRAPH, Csts::MAN_CENTER, Csts::MAN_CENTER, true, true);  //Note: PM style is NOT used

    //Parameters
    int nf       = fbpl.nf;
    double c1    = fbpl.cs_em.c1;
    double gamma = fbpl.cs_em.gamma;
    double mu_EM = fbpl.us_em.mu_EM;
    double t1    = fbpl.us_em.T;  //should be equal to 2pi here

    //--------------------------------------------------------------------
    // 2. (fac) FFT of the true anomaly
    //--------------------------------------------------------------------
    double M, E;
    double e = fbpl.us_em.lecc;

    //-------------------------------------
    //Set E[M]-M in dEv, for M = [0, ..., 2pi]
    //-------------------------------------
    int N = 2*nf+1;
    gsl_vector *dEv = gsl_vector_calloc(N);
    //Set E[i] in dEv
    for(int k = 0; k < N; k++)
    {
        M = t1*k/(double)N;  //M = t;
        E = rtnewtecc(eccDer, M, 1e-10, e, M); //E = f(M);
        gsl_vector_set(dEv, k, E-M);
    }

    //-------------------------------------
    //Actual FFT
    //-------------------------------------
    Ofsc xFFT(nf);
    eccFFT(xFFT, nf, N, dEv, 0); //E-M is odd: its a series of sinus
    //cout << "xFFT = " << endl << xFFT << endl;


    //-------------------------------------
    //With Bessel function
    //-------------------------------------
    Ofsc xFFT2(nf);
    eccBessel(xFFT2, nf, e);
    //cout << "xFFT2 = " << endl << xFFT2 << endl;

    //-------------------------------------
    //Test on a given grid
    //-------------------------------------
    //Allocation
    int fftPlot = 100;
    gsl_vector * xxL  = gsl_vector_calloc(fftPlot);
    gsl_vector * xxL2 = gsl_vector_calloc(fftPlot);

    //Test loop
    for(int i = 0; i < fftPlot; i++)
    {
        M = t1*i/fftPlot;
        E = rtnewtecc(eccDer, M, 1e-10, e, M);
        gsl_vector_set(xxL, i, fabs(creal(xFFT.evaluate(M))-(E-M)));
        gsl_vector_set(xxL2, i, fabs(creal(xFFT2.evaluate(M))-(E-M)));
    }
    cout << std::noshowpos << setiosflags(ios::scientific)  << setprecision(15);
    cout << std::showpos << "Error max on the eccentric anomaly: " << endl;
    cout << "from FFT    = " << gsl_vector_max(xxL) << endl;
    cout << "from Bessel = " << gsl_vector_max(xxL2) << endl;



    //--------------------------------------------------------------------
    // 3. FFT of Phi(nu)
    //--------------------------------------------------------------------
    double Phi, nu;
    //-------------------------------------
    //Set Phi(nu) in Phiv, for nu = [0, ..., 2pi]
    //-------------------------------------
    gsl_vector *Phiv = gsl_vector_calloc(N);
    //Set E[i] in dEv
    for(int k = 0; k < N; k++)
    {
        nu = t1*k/(double)N;  //nu = t;
        Phi = 1.0/(1.0+e*cos(nu));
        gsl_vector_set(Phiv, k, Phi);
    }

    //-------------------------------------
    //Actual FFT
    //-------------------------------------
    xFFT.zero();
    eccFFT(xFFT, nf, N, Phiv, 1);
    //cout << "xFFT = " << endl << xFFT << endl;

    //-------------------------------------
    //Test on a given grid
    //-------------------------------------
    //Test loop
    for(int i = 0; i < fftPlot; i++)
    {
        nu = t1*i/(double)fftPlot;  //nu = t;
        Phi = 1.0/(1.0+e*cos(nu));
        gsl_vector_set(xxL, i, fabs(creal(xFFT.evaluate(nu))-Phi));
        //cout << "nu = " << nu << ", eps = " << fabs(creal(xFFT.evaluate(nu))-Phi) << endl;
        //cout << "Phi = " << Phi << ", xFFT(nu) = " << creal(xFFT.evaluate(nu)) << endl;
    }
    cout << std::noshowpos << setiosflags(ios::scientific)  << setprecision(15);
    cout << std::showpos << "Error max on Phi(nu): " << endl;
    cout << "from FFT    = " << gsl_vector_max(xxL) << endl;

    //Free GSL objects
    gsl_vector_free(dEv);
    gsl_vector_free(xxL);
    gsl_vector_free(xxL2);


    //--------------------------------------------------------------------
    // 4. Set all alpha functions
    //--------------------------------------------------------------------
    Ofs<cdouble > alpha1c(nf);
    Ofs<cdouble > alpha2c(nf);
    Ofs<cdouble > alpha3c(nf);
    Ofs<cdouble > alpha4c(nf);
    Ofs<cdouble > alpha5c(nf);
    Ofs<cdouble > alpha6c(nf);
    Ofs<cdouble > alpha7c(nf);
    Ofs<cdouble > alpha8c(nf);
    Ofs<cdouble > alpha9c(nf);
    Ofs<cdouble > alpha10c(nf);
    Ofs<cdouble > alpha11c(nf);
    Ofs<cdouble > alpha12c(nf);
    Ofs<cdouble > alpha13c(nf);
    Ofs<cdouble > alpha14c(nf);
    Ofs<cdouble > alpha15c(nf);

    //Redundancy for the positions of the primaries
    Ofs<cdouble > Xe(nf);
    Ofs<cdouble > Ye(nf);
    Ofs<cdouble > Ze(nf);
    Ofs<cdouble > Xm(nf);
    Ofs<cdouble > Ym(nf);
    Ofs<cdouble > Zm(nf);
    Ofs<cdouble > Xs(nf);
    Ofs<cdouble > Ys(nf);
    Ofs<cdouble > Zs(nf);

    //Positions of the primaries in NC coordinates
    Ofs<cdouble > xe(nf);
    Ofs<cdouble > ye(nf);
    Ofs<cdouble > ze(nf);
    Ofs<cdouble > xm(nf);
    Ofs<cdouble > ym(nf);
    Ofs<cdouble > zm(nf);
    Ofs<cdouble > xs(nf);
    Ofs<cdouble > ys(nf);
    Ofs<cdouble > zs(nf);

    //-------------------------
    // If not state otherwise, alphai=0
    // i.e. alpha2,4,5 = 0
    //-------------------------
    //alpha1 = 1
    alpha1c.setCoef(1.0+0.0*I, 0);
    //alpha3 = 1
    alpha3c.setCoef(1.0+0.0*I, 0);
    //alpha6=Phi
    alpha6c.ccopy(xFFT);
    //alpha13 = -c1*Phi
    alpha13c.ofs_smult(xFFT, -c1 + 0.0*I);
    //alpha15 = Phi-1
    alpha15c.ccopy(xFFT);
    alpha15c.addCoef(-1.0+0.0*I, 0);


    //-------------------------
    // Primaries
    //-------------------------
    //Earth
    alpha9c.setCoef(+mu_EM+0.0*I, 0);
    //Moon
    alpha11c.setCoef(+mu_EM-1.0+0.0*I, 0);

    //-------------------------
    // The primaries, again
    //-------------------------
    //Earth, EM coordinates
    Xe.setCoef(+mu_EM+0.0*I, 0);
    //Earth, NC coordinates
    xe.setCoef(c1 - mu_EM/gamma + 0.0*I, 0);
    //Moon, EM coordinates
    Xm.setCoef(+mu_EM-1.0+0.0*I, 0);
    //Moon, NC coordinates
    xm.setCoef(c1 - (mu_EM-1.0)/gamma + 0.0*I, 0);

    //--------------------------
    //Put in data file
    //--------------------------
    ofs_sst(alpha1c, fbpl.cs_em.F_COEF+"alpha1", 1, "_fft");
    ofs_sst(alpha2c, fbpl.cs_em.F_COEF+"alpha2", 0, "_fft");
    ofs_sst(alpha3c, fbpl.cs_em.F_COEF+"alpha3", 1, "_fft");
    ofs_sst(alpha4c, fbpl.cs_em.F_COEF+"alpha4", 1, "_fft");
    ofs_sst(alpha5c, fbpl.cs_em.F_COEF+"alpha5", 0, "_fft");
    ofs_sst(alpha6c, fbpl.cs_em.F_COEF+"alpha6", 1, "_fft");
     //Sun
    ofs_sst(alpha7c, fbpl.cs_em.F_COEF+"alpha7", 1, "_fft");
    ofs_sst(alpha8c, fbpl.cs_em.F_COEF+"alpha8", 0, "_fft");
    //Earth
    ofs_sst(alpha9c,  fbpl.cs_em.F_COEF+"alpha9",  1, "_fft");
    ofs_sst(alpha10c, fbpl.cs_em.F_COEF+"alpha10", 0, "_fft");
    //Moon
    ofs_sst(alpha11c, fbpl.cs_em.F_COEF+"alpha11", 1, "_fft");
    ofs_sst(alpha12c, fbpl.cs_em.F_COEF+"alpha12", 0, "_fft");
    //NC additional coeffs
    ofs_sst(alpha13c, fbpl.cs_em.F_COEF+"alpha13", 1, "_fft");
    ofs_sst(alpha14c, fbpl.cs_em.F_COEF+"alpha14", 0, "_fft");
    //ERTBP additional coeff
    ofs_sst(alpha15c, fbpl.cs_em.F_COEF+"alpha15", 1, "_fft");


    //---------------
    //Primary, EM coordinates
    //Note that, at this step, the vertical motion of the primaries is undefined,
    //so we can put either Even or Odd in the ofs_sst option of gsl_Zc without much difference
    //---------------
    ofs_sst(Xs, fbpl.cs_em.F_COEF+"Ps1", 1, "_fft");
    ofs_sst(Ys, fbpl.cs_em.F_COEF+"Ps2", 0, "_fft");
    ofs_sst(Zs, fbpl.cs_em.F_COEF+"Ps3", 1, "_fft");

    ofs_sst(Xm, fbpl.cs_em.F_COEF+"Pm1", 1, "_fft");
    ofs_sst(Ym, fbpl.cs_em.F_COEF+"Pm2", 0, "_fft");
    ofs_sst(Zm, fbpl.cs_em.F_COEF+"Pm3", 1, "_fft");

    ofs_sst(Xe, fbpl.cs_em.F_COEF+"Pe1", 1, "_fft");
    ofs_sst(Ye, fbpl.cs_em.F_COEF+"Pe2", 0, "_fft");
    ofs_sst(Ze, fbpl.cs_em.F_COEF+"Pe3", 1, "_fft");


    //---------------
    //Primary, NC coordinates
    //Note that, at this step, the vertical motion of the primaries is undefined,
    //so we can put either Even or Odd in the ofs_sst option of gsl_zc without much difference
    //---------------
    ofs_sst(xs, fbpl.cs_em.F_COEF+"ps1", 1, "_fft");
    ofs_sst(ys, fbpl.cs_em.F_COEF+"ps2", 0, "_fft");
    ofs_sst(zs, fbpl.cs_em.F_COEF+"ps3", 1, "_fft");

    ofs_sst(xm, fbpl.cs_em.F_COEF+"pm1", 1, "_fft");
    ofs_sst(ym, fbpl.cs_em.F_COEF+"pm2", 0, "_fft");
    ofs_sst(zm, fbpl.cs_em.F_COEF+"pm3", 1, "_fft");

    ofs_sst(xe, fbpl.cs_em.F_COEF+"pe1", 1, "_fft");
    ofs_sst(ye, fbpl.cs_em.F_COEF+"pe2", 0, "_fft");
    ofs_sst(ze, fbpl.cs_em.F_COEF+"pe3", 1, "_fft");
}


