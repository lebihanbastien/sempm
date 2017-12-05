# Configuration file for sempm - See config/constants.h for details on the constants.
#
# This file calls sempm in order to compute the change of coordinates that set the 
# Hamiltonian of the QBCP about a given libration point ($LIBPOINT) in normal diagonal
# form at order two.
#
# BLB 2017

#-----------------------------------------------------------------------------------------
# TYPE OF COMPUTATION
#-----------------------------------------------------------------------------------------
COMPTYPE=$NFO2

#-----------------------------------------------------------------------------------------
# Numerical constants
#-----------------------------------------------------------------------------------------
OFTS_ORDER=8   # (not used for NFO2)
OFS_ORDER=30
STORAGE=$TRUE

#-----------------------------------------------------------------------------------------
# MODEL: RTBP = 0; QBCP = 1; BCP = 2
#-----------------------------------------------------------------------------------------
MODEL=$QBCP  

#-----------------------------------------------------------------------------------------
# LIBRATION POINT
#-----------------------------------------------------------------------------------------
LIBPOINT=$EML2



