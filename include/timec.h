#ifndef TIMEC_H_INCLUDED
#define TIMEC_H_INCLUDED

/**
 * \file timec.h
 * \brief Time manipulation similar to MATLAB routines.
 * \author BLB, Zubin Olikara
 */

#include <stdio.h>
#include <time.h>


extern struct timespec TIC_TIME;


//Equivalent to the tic function in Matlab
void tic(void);
//Equivalent to the toc function in Matlab
double toc(void);


#endif // TIMEC_H_INCLUDED
