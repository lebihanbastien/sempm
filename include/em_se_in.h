#ifndef EMINSEM_H_INCLUDED
#define EMINSEM_H_INCLUDED

/**
 * \file  eminsem.h
 * \brief Contains all the routines to perform changes of coordinates between the EM and SEM frameworks. Including
 * \author BLB.
 * \date   2016
 * \version 1.0
 */

#include "vf.h"
#include "qbtbp.h"
#include "init.h"


//----------------------------------------------------------------------------------------
// COC: Velocities <--> Momenta
//----------------------------------------------------------------------------------------
/**
 *  \brief Change the SE velocities into SE momenta
 **/
void se_v_to_se_m(double t, const double zv_se_v[], double zv_se_m[], void *params_void);

/**
 *  \brief Change the SE momenta into SE velocities
 **/
void se_m_to_se_v(double t, const double zv_se_m[], double zv_se_v[], void *params_void);

/**
 *  \brief Change the EM velocities into EM momenta
 **/
void em_v_to_em_m(double t, const double zv_em_v[], double zv_em_m[], void *params_void);

/**
 *  \brief Change the EM momenta into EM velocities
 **/
void em_m_to_em_v(double t, const double zv_em_m[], double zv_em_v[], void *params_void);

//----------------------------------------------------------------------------------------
// Change of unit system
//----------------------------------------------------------------------------------------
/**
 *   \brief From SE units to EM units for a (position, velocity) state vector and time
 *          given in IN coordinates.
 **/
void units_se_to_em(double *tc, double zv_in[], FBPL *fbpl);

/**
 *   \brief From EM units to SE units for a (position, velocity) state vector and time
 *          given in IN coordinates.
 **/
void units_em_to_se(double *tc, double zv_in[], FBPL *fbpl);

//----------------------------------------------------------------------------------------
// COC: IN <--> EM
//----------------------------------------------------------------------------------------
/**
 * \brief From EM to IN (in EM units)
 **/
void em_v_to_in(double t, const double zv_em_v[], double zv_in[], FBPL *fbpl);

/**
 * \brief From IN to EM (in EM units)
 **/
void in_to_em_v(double t, const double zv_in[], double zv_em_v[], FBPL *fbpl);

//----------------------------------------------------------------------------------------
// COC: IN <--> SE
//----------------------------------------------------------------------------------------
/**
 * \brief From SE to IN (in SE units)
 **/
void se_v_to_in(double t, const double zv_se_v[], double zv_in[], FBPL *fbpl);

/**
 * \brief From IN to SE (in SE units)
 **/
void in_to_se_v(double t, const double zv_in[], double zv_se_v[], FBPL *fbpl);

//----------------------------------------------------------------------------------------
// COC: SE <--> EM
//----------------------------------------------------------------------------------------
/**
 * \brief From SE to EM (both in position/momenta form)
 **/
void se_m_to_em_m(double t, const double zv_se_m[], double zv_em_m[], FBPL *fbpl);

/**
 * \brief From EM to SE (both in position/momenta form)
 **/
void em_m_to_se_m(double t, const double zv_em_m[], double zv_se_m[], FBPL *fbpl);

//----------------------------------------------------------------------------------------
// COC: SE <--> NCEM
//----------------------------------------------------------------------------------------
/**
 * \brief From NCEM to SE (both in position/momenta form)
 **/
void ncem_m_to_se_m(double t, const double zv_ncem_m[], double zv_se_m[], FBPL *fbpl);

//----------------------------------------------------------------------------------------
// COC: NCSEM <--> NCEM
//----------------------------------------------------------------------------------------
/**
 * \brief From NCSE to  NCEM (both in position/momenta form)
 **/
void ncse_m_to_ncem_m(double t, const double zv_ncse_m[], double zv_ncem_m[], FBPL *fbpl);

/**
 * \brief From NCEM to  NCSE (both in position/momenta form)
 **/
void ncem_m_to_ncse_m(double tEM, const double zv_ncem_m[], double zv_ncse_m[], FBPL *fbpl);

//----------------------------------------------------------------------------------------
// COC: NCSYS <--> SYS
//----------------------------------------------------------------------------------------
/**
 *  \brief COC: NC(EM or SE) coordinates to SYS (EM or SE) coordinates.
 **/
void ncsys_m_to_sys_m(double t, const double zv_ncsys_m[], double zv_sys_m[], FBPL *qbp);

/**
 *  \brief COC: from SYS (EM or SE) coordinates to NC(EM or SE) coordinates.
 **/
void sys_m_to_ncsys_m(double t, const double zv_sys_m[], double zv_ncsys_m[], FBPL *qbp);

//----------------------------------------------------------------------------------------
// COC: NCEM <--> EM
//----------------------------------------------------------------------------------------
/**
 *  \brief COC: from NCEM coordinates to EM coordinates
 **/
void ncem_m_to_em_m(double t, const double zv_ncem_m[], double zv_em_m[], FBPL *qbp);

/**
 *  \brief COC: from EM coordinates to NCEM coordinates
 **/
void em_m_to_ncem_m(double t, const double zv_em_m[], double zv_ncem_m[], FBPL *qbp);

//----------------------------------------------------------------------------------------
// COC: NCSE <--> SE
//----------------------------------------------------------------------------------------
/**
 *  \brief COC: from SE coordinates to NCSE coordinates
 **/
void se_m_to_ncse_m(double t, const double zv_se_m[], double zv_ncse_m[], FBPL *qbp);

/**
 *  \brief COC: from NCSE coordinates to SE coordinates
 **/
void ncse_m_to_se_m(double t, const double zv_ncse_m[], double zv_se_m[], FBPL *qbp);

#endif // EMINSEM_H_INCLUDED
