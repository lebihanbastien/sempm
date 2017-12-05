# Configuration file for sempm - See config/constants.h for details on the constants.
#
# This file calls sempm in order to solve the Quasi-Bicircular Three and Four-Body 
# Problem (QBTBP and QBCP). If STORAGE==$TRUE, the coefficients of both problems
# are stored in ./data/qbtbp and ./data/VF/QBCP/.
#

#-----------------------------------------------------------------------------------------
# TYPE OF COMPUTATION
#-----------------------------------------------------------------------------------------
COMPTYPE=$QBTBP

#-----------------------------------------------------------------------------------------
# Numerical constants
#-----------------------------------------------------------------------------------------
OFTS_ORDER=8
OFS_ORDER=30
STORAGE=$TRUE

#-----------------------------------------------------------------------------------------
# MODEL: RTBP = 0; QBCP = 1; BCP = 2
#-----------------------------------------------------------------------------------------
MODEL=$QBCP  

#-----------------------------------------------------------------------------------------
# LIBRATION POINT (not used for QBTBP)
#-----------------------------------------------------------------------------------------
LIBPOINT=$EML2




