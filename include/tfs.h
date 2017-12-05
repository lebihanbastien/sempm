#ifndef TFS_H_INCLUDED
#define TFS_H_INCLUDED

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#include "ofs.h"

/**
 *  Fourier series structure
 **/
typedef struct TFS TFS;
struct TFS
{
    int order;
    double complex *coef;
};


/**
 *   \brief Initialize a TFS structure
 **/
void tfs_init(TFS *tfs, int order);
/**
 *   \brief From Ofs to Tfs object
 **/
void tfs_from_ofs(TFS *tfs, Ofs<double complex>& ofs, double n);
/**
 *   \brief To Ofs from Tfs object
 **/
void tfs_to_ofs(TFS *tfs, Ofs<double complex>& ofs, double n);
/**
 *   \brief Print a TFS structure
 **/
void tfs_printf(TFS *tfs);
/**
 *   \brief Kills a TFS structure
 **/
void tfs_free(TFS *tfs);

#endif // TFS_H_INCLUDED
