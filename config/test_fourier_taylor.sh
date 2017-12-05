# Configuration file for sempm - See config/constants.h for details on the constants.
#
# This file calls sempm in order to test the routines of the Fourier-Taylor algebra.
# The algebra is tested on Fourier-Taylor series of order OFTS_ORDER (for the Taylor series)
# and OFS_ORDER (for the Fourier coefficients), and of REDUCED_NV=4 variables.
# Note that REDUCED_NV is automatically set to 4 in ooftda.sh, but any other value could
# be used to test the algebra.

#-----------------------------------------------------------------------------------------
# TYPE OF COMPUTATION
#-----------------------------------------------------------------------------------------
COMPTYPE=$FT_TEST

#-----------------------------------------------------------------------------------------
# Numerical constants
#-----------------------------------------------------------------------------------------
OFTS_ORDER=8
OFS_ORDER=30
REDUCED_NV=2
STORAGE=$FALSE #(not used for FT_TEST)

#-----------------------------------------------------------------------------------------
# MODEL: RTBP = 0; QBCP = 1; BCP = 2 (not used for FT_TEST)
#-----------------------------------------------------------------------------------------
MODEL=$QBCP  

#-----------------------------------------------------------------------------------------
# LIBRATION POINT (not used for FT_TEST)
#-----------------------------------------------------------------------------------------
LIBPOINT=$EML2




