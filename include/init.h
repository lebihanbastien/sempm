#ifndef CONFIG_H_INCLUDED
#define CONFIG_H_INCLUDED

/**
 * \file init.h
 * \brief Configuration file. It allows to initialize the environment of the Sun-(Earth+Moon) Quasi-Bicircular Four-Body Problem.
 * \author BLB.
 * \date May 2015
 * \version 1.0
 */

#include <vector>
#include "ofts.h"
#include "matrix.h"
#include "env.h"

//========================================================================================
// Global structures
//========================================================================================
extern FBP SEM;   //global structure that describes the Sun-Earth-Moon system
extern FBPL SEML; //global structure that describes the Sun-Earth-Moon system around Li

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
void init_env(int t_model, int t_li, int t_man_type, int t_pms);

//========================================================================================
// Global variables for parameterization of the center manifold
//========================================================================================
extern vector<Oftsc>  CM;     //center manifold in NC coordinates
extern vector<Oftsc>  CMdot;  //time derivative of the center manifold in NC coordinates
extern vector<Oftsc> CMh;     //center manifold in TFC coordinates
extern matrix<Oftsc> JCM;     //jacobian of CM
extern vector<Oftsc>  Fh;     //reduced vector field
extern vector<Oftsc>  DWFh;   //DWFh =JCM * Fh

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
void init_inv_man(FBPL &fbpl);

/**
 *   \brief Update of the parameterization center of an invariant manifold around
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
void update_inv_man(FBPL &qbcp);


//========================================================================================
// COC: change of coordinates to achieve a diagonal form at order 2 of the Hamiltonian
//========================================================================================
extern matrix<Ofsc>  Mcoc;    //COC matrix
extern matrix<Ofsc>  Pcoc;    //COC matrix (Mcoc = Pcoc*Complex matrix)
extern matrix<Ofsc>  MIcoc;   //COC matrix = inv(Mcoc)
extern matrix<Ofsc>  Qcoc;    //COC matrix = inv(Pcoc)
extern vector<Ofsc>  Vcoc;    //COC vector

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
void init_coc(FBPL &qbcp);

/**
 *  \brief The several variables of the COC are retrieved from txt files, stored in the
 *         folder F_COC defined in fbpl.
 **/
void init_coc(matrix<Ofsc> &t_Pcoc, matrix<Ofsc> &t_Mcoc, matrix<Ofsc> &t_Qcoc,
             matrix<Ofsc> &t_MIcoc, vector<Ofsc> &t_Vcoc, FBPL& fbpl);

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
void read_coc(Ofsc& xFFT, string prefix, int k, int p);

/**
 *   \brief Number to string inner routine, using static_cast.
 *          Example: num_to_string(10) returns "10".
 **/
string num_to_string(double num);


#endif // CONFIG_H_INCLUDED
