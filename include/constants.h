#ifndef CSTS_H
#define CSTS_H


/**
 * \file   constants.h
 * \brief  Definition of the class Csts, a simple class of numerical constants.
 * \author BLB
 * \date   2017
 */

#include <complex.h>

//----------------------------------------------------------------------------------------
//   typedef
//----------------------------------------------------------------------------------------
typedef complex double cdouble;


/**
 *  A simple class of constants.
 **/
class Csts
{
protected:
private:
public:
    //====================================================================================
    // Configuration of sempm: argument passed to the main routine (must agree with the
    // BASH constants defined in sempm/config/constants.sh
    //====================================================================================

    //------------------------------------------------------------------------------------
    // TYPE OF COMPUTATION (COMPTYPE)
    //------------------------------------------------------------------------------------
    static const int QBTBP=0;
    static const int NFO2=1;
    static const int PM=2;
    static const int PM_TEST=3;
    static const int COMPMAP=4;
    static const int COC=5;
    static const int DYNEQ=6;
    static const int TRAJ=7;
    static const int FT_TEST=8;

    //------------------------------------------------------------------------------------
    // TYPE OF MANIFOLD (MANTYPE)
    //------------------------------------------------------------------------------------
    static const int MAN_CENTER=0;
    static const int MAN_CENTER_S=1;
    static const int MAN_CENTER_U=2;
    static const int MAN_CENTER_US=3;

    //------------------------------------------------------------------------------------
    // MODEL (MODEL)
    // CRTBP = 0; QBCP = 1; BCP = 2; ERTBP = 3
    //------------------------------------------------------------------------------------
    static const int CRTBP=0;
    static const int QBCP=1;
    static const int BCP=2;
    static const int ERTBP=3;

    //------------------------------------------------------------------------------------
    // COORDINATE SYSTEM (CS) - Can also be used for UNIT SYSTEM (US) or FRAMEWORK.
    // EM = 0; SEM = 1; SE = 2
    //------------------------------------------------------------------------------------
    static const int EM=0;
    static const int SEM=1;
    static const int SE=2;


    //------------------------------------------------------------------------------------
    // PM STYLE (PMS)
    // GRAPH = 0; NORMAL FORM = 1; MIXED = 2
    //------------------------------------------------------------------------------------
    static const int GRAPH=0;
    static const int NORMFORM=1;
    static const int MIXED=2;

    //------------------------------------------------------------------------------------
    // PMAP TYPE (PMAP_TYPE)
    // PMAP = 1; TMAP = 2; EMAP = 3; IMAP = 4
    //------------------------------------------------------------------------------------
    static const int PMAP=1;
    static const int TMAP=2;
    static const int EMAP=3;
    static const int IMAP=4;
    static const int HMAP=5;
    static const int IMAPPLANAR=6;

    //------------------------------------------------------------------------------------
    // PMAP METHODS (PMAP_method)
    // DUAL_INT = 1; DUAL_INT_NO_RESET = 2;
    // DUAL_INT_STEPPED = 3; SINGLE_INT = 4.
    //------------------------------------------------------------------------------------
    static const int DUAL_INT=1;
    static const int DUAL_INT_NO_RESET=2;
    static const int DUAL_INT_STEPPED=3;
    static const int SINGLE_INT=4;

    //------------------------------------------------------------------------------------
    // LIBRATION POINT (LIBPOINT)
    //------------------------------------------------------------------------------------
    static const int EML1=1;
    static const int EML2=2;
    static const int SEL1=3;
    static const int SEL2=4;

    static const int EML3=5;
    static const int SEL3=6;

    //====================================================================================
    // Constants for NFO2 computation
    //====================================================================================
    //------------------------------------------------------------------------------------
    // Type of computation for the stable direction:
    // - STABLE_DIR_POW: using the inverse power method.
    // - STABLE_DIR_SYM: using the fact the stable direction
    //   can be obtained by symmetry from the unstable direction
    //------------------------------------------------------------------------------------
    static const int STABLE_DIR_POW = 1;
    static const int STABLE_DIR_SYM = 2;

    //------------------------------------------------------------------------------------
    // Computing the inverse of a square matrix:
    // - INVERSE_SYMP: inverse of a symplectic matrix M: inv(M) = -J * trans(M) * J
    // - INVERSE_GSL: using GSL library
    //------------------------------------------------------------------------------------
    static const int INVERSE_SYMP=1;
    static const int INVERSE_GSL=2;

    //====================================================================================
    // Environment. The rest of the constants (masses, distances...) are defined in
    // the routines of env.h
    //====================================================================================
    //------------------------------------------------------------------------------------
    // Numerotation of the Solar System planets and objects,
    // consistent with JPL's HORIZON numerotation.
    //------------------------------------------------------------------------------------
    static const int SUN=10;
    static const int MERCURY=199;
    static const int VENUS =99;
    static const int EARTH=399;
    static const int MARS=499;
    static const int JUPITER=599;
    static const int SATURN=699;
    static const int URANUS=799;
    static const int NEPTUNE=899;
    static const int PLUTO=999;
    static const int MOON=301;
    static const int EARTH_AND_MOON=700;

    //====================================================================================
    //   ORDER AND NUMBER OF VARIABLES FOR OFS AND OFTS OBJECTS
    //====================================================================================
    static const int POTENTIAL_ORDER=60; // Maximum order of the potential of the primaries
    static const int NV=6;               // Number of state variables (a priori always 6)
    static const int OFS_NV=1;           // Number of variables in the OFS object (a priori always 1)


    //====================================================================================
    //   Maps
    //====================================================================================
    static const int MAX_EVENTS=100;  // Maximum events in Poincare & Period maps

    //====================================================================================
    //   Plotting
    //====================================================================================
    static const int XY=1;  // XY plane, for plotting purposes
    static const int YZ=2;  // YZ plane, for plotting purposes
    static const int XZ=3;  // XZ plane, for plotting purposes


    //====================================================================================
    //   Aesthetics
    //====================================================================================
    static const char SSEPR[]; //= "---------------------------------";
    static const char SEPR[];  //= "---------------------------------------------------";
};

#endif // CSTS_H
