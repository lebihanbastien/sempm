# Configuration file for sempm

#-----------------------------------------------------------------------------------------
# TYPE OF COMPUTATION (COMPTYPE)
#-----------------------------------------------------------------------------------------
COMPTYPE=$COMPMAP

#-----------------------------------------------------------------------------------------
# Numerical constants
#-----------------------------------------------------------------------------------------
OFTS_ORDER=8
OFS_ORDER=30

#-----------------------------------------------------------------------------------------
# STORAGE
# FALSE = 0; TRUE = 1
#-----------------------------------------------------------------------------------------
STORAGE=$TRUE


#-----------------------------------------------------------------------------------------
# MODEL
# RTBP = 0; QBCP = 1; BCP = 2
#-----------------------------------------------------------------------------------------
MODEL=$QBCP  

#-----------------------------------------------------------------------------------------
# LIBRATION POINT
#-----------------------------------------------------------------------------------------
LIBPOINT=$EML2


#-----------------------------------------------------------------------------------------
# TYPE OF MANIFOLD
# MAN_CENTER    = 0
# MAN_CENTER_S  = 1
# MAN_CENTER_U  = 2
# MAN_CENTER_US = 3
#-----------------------------------------------------------------------------------------
MANTYPE=$MAN_CENTER

#-----------------------------------------------------------------------------------------
# PM STYLE
# GRAPH = 0; NORMFORM = 1; MIXED = 2
#-----------------------------------------------------------------------------------------
PMS=$GRAPH

#-----------------------------------------------------------------------------------------
# OpenMP settings
#-----------------------------------------------------------------------------------------
NUM_THREADS=2


#-----------------------------------------------------------------------------------------
# PMAP TYPE: deprecated parameter, fixed to PMAP
#-----------------------------------------------------------------------------------------
PMAP_TYPE=$PMAP

#-----------------------------------------------------------------------------------------
# PMAP METHODS: deprecated parameter, fixed to SINGLE_INT
#-----------------------------------------------------------------------------------------
PMAP_int_method=$SINGLE_INT

#-----------------------------------------------------------------------------------------
# PMAP SETTINGS
#-----------------------------------------------------------------------------------------
PMAP_ofts_order=$OFTS_ORDER #the effective order of the Taylor expansions on the map
PMAP_ofs_order=$OFS_ORDER   #the effective order of the Fourier expansions on the map
PMAP_dh_nc=0.010    #the initial relative energy at t = t0, in NC coordinates
PMAP_f_proj_T=-1  #the frequency of projection inside the map (given as a fraction of the SEM period T).  If equal to -1, no aditionnal projection is performed	
PMAP_t0_T=0  #the initial time for all trajectories on the map (x T)
PMAP_tmax_T=1000	#the maximum time for all trajectories on the map (x T) (usually set to a very large value)
PMAP_max_events=500  #the maximum number of events on the map, i.e. the number of crossings of the xy-plane, with dot(z) > 0
PMAP_n_sol=50 #number of computed trajectories
PMAP_si_min=-20.0  #minimum value for the initial conditions (s1, s3) on the map
PMAP_si_max=20.0  #maximum value for the initial conditions (s1, s3) on the map
PMAP_is_plot=1; #plot boolean
PMAP_is_par=0;  #parallel computation of not









