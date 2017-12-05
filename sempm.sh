#!/bin/sh
# script file for sempm
# RECALL: make it executable with "$ chmod +x script.sh"




#-----------------------------------------------------------------------------------------
# Source the parameters
#-----------------------------------------------------------------------------------------
source "config/constants.sh"


#-----------------------------------------------------------------------------------------
# Source the tools (bash functions)
#-----------------------------------------------------------------------------------------
source "config/tools.sh"


#-----------------------------------------------------------------------------------------
# Title
#-----------------------------------------------------------------------------------------
echo
echo "#------------------------------------------#"
echo "#          sempm configuration            #"
echo "#------------------------------------------#"

#-----------------------------------------------------------------------------------------
# Start with rm of all gnuplot_tmpdatafile_* in tmp
#-----------------------------------------------------------------------------------------
echo "The tmp folder is cleared..."
rm -f tmp/gnuplot_tmpdatafile_*


#-----------------------------------------------------------------------------------------
# Source the data folder generation
#-----------------------------------------------------------------------------------------
echo "The data and plot folders are updated if needed..."
source "config/mkdirdata.sh"
	
# DATA_DIR="data"
# if [ -d "$DATA_DIR" ]; then
#     echo "The data folder already exists"
#     
# else
# 	echo "The data folder is created"
# 	source "config/mkdirdata.sh"
# fi


#-----------------------------------------------------------------------------------------
# If no argument is passed, use the yad, if yad is installed. Else, use $1
#-----------------------------------------------------------------------------------------
# Minimum yad version accepted
yad_version_min="0.36.3"

# Check if an argument has been used 
if [ $# -eq 0 ]; then
	
	#-------------------------------------------------------------------------------------
	# Check yad version
	#-------------------------------------------------------------------------------------
	echo "No configuration file has been passed. The GUI is used instead."
	echo "Checking current Yad version...."
	
	yad_version_ok=1
	# First, check that yad exists
	command -v yad >/dev/null 2>&1 || { yad_version_ok=0; }
	
	#-------------------------------------------------------------------------------------
	# Then, the version
	#-------------------------------------------------------------------------------------
	yad_version="$(yad --version)"
	yad_version_c=$(echo $yad_version | awk 'BEGIN {FS="(" } { print $1 }')
	
	if version_gt $yad_version_c $yad_version_min; then
		 echo "... is" $yad_version
	else
		 yad_version_ok=0
		 
	fi				
	
	# Abort if necessary
	if [ $yad_version_ok -eq 0 ]; then
		 echo "Error, you need at least version " $yad_version_min " of Yad to use the GUI."
		 echo "Yad is available here: https://sourceforge.net/projects/yad-dialog/"
		 echo "Without the GUI, you must use a configuration file as argument of ooftda.sh. Exit."
		 exit 1;
	fi
	
	#-------------------------------------------------------------------------------------
	# Launch the GUI if everything is okay
	#-------------------------------------------------------------------------------------
	echo "Launching the GUI..."
    source "config/yad/ooftda_yad_notebook.sh"
    source_file="from GUI"
else
	source $1
	source_file="$1"
fi


#-----------------------------------------------------------------------------------------
# Possible Remarks (at the end):
#-----------------------------------------------------------------------------------------
REMARKS="#------------------------------------------#\n"
REMARKS+="Additional remarks:\n\n"

#-----------------------------------------------------------------------------------------
# Source file
#-----------------------------------------------------------------------------------------
echo "#------------------------------------------#"
echo "Source file:" $source_file
echo ""

#-----------------------------------------------------------------------------------------
# Additional tests
#-----------------------------------------------------------------------------------------
#Check that the variable $COMP_TYPE exist
if [ -z ${COMPTYPE+x} ]; then
	echo 'WARNING: the variable COMPTYPE is not set.'
	echo 'No computation.'
	exit	
fi

#-----------------------------------------------------------------------------------------
# Normalization (ISNORM) if FIXED to TRUE
#-----------------------------------------------------------------------------------------
ISNORM=$TRUE
REMARKS+="* ISNORM is deprecated and has been fixed to TRUE.\n"

#-----------------------------------------------------------------------------------------
# Switch to update CS, LI_EM, LI_SEM wrt LIBPOINT
#-----------------------------------------------------------------------------------------
# Check if LIBPOINT exist
if [ -z ${LIBPOINT+x} ]; then
	set_param "LIBPOINT" $EML2
fi

# Define CS, LI_EM, LI_SEM
case $LIBPOINT in
	$EML1)   
		CS=$EM
		LI_EM=1
		LI_SEM=1		
		;;
	$EML2)   
		CS=$EM
		LI_EM=2
		LI_SEM=1		
	;;
	$SEL1)   
		CS=$SEM
		LI_EM=1
		LI_SEM=1		
	;;
	$SEL2)   
		CS=$SEM
		LI_EM=1
		LI_SEM=2		
	;;
esac



#----------------------------------------
#Check that the variable $MODEL exist
#----------------------------------------
if [ -z ${MODEL+x} ]; then
	set_param "MODEL" $M_QBCP
fi

#----------------------------------------
# OFS_ORDER = Order of the Fourier series
#----------------------------------------
if [ "$MODEL" == "$CRTBP" ]; then
		OFS_ORDER=0
else  
		OFS_ORDER=30
fi

#-----------------------------------------------------------------------------------------
# Define PMS if it has not been defined
#-----------------------------------------------------------------------------------------
if [ -z ${PMS+x} ]; then
	set_param "PMS" $GRAPH
	REMARKS+="* The variable PMS, which was not defined, has been set to GRAPH.\n"

fi


#----------------------------------------
# MANTYPE    = Type of manifold
# REDUCED_NV = Number of reduced variables
#----------------------------------------
#Check that the variable $MANTYPE exist
if [ -z ${MANTYPE+x} ]; then
	#If not, default values to center manifold
	REMARKS+="* The variable MANTYPE is not set.\n"
	REMARKS+="  By default, MANTYPE=MAN_CENTER and REDUCED_NV=4.\n"
	MANTYPE=$MAN_CENTER
	REDUCED_NV=4
else
	#If yes, set the REDUCED_NV in correspondance
	case $MANTYPE in
		$MAN_CENTER)   
			        REDUCED_NV=4
			;;
		$MAN_CENTER_S)  
			        REDUCED_NV=5
			        PMS=$MIXED
			;;
		$MAN_CENTER_U)  
				REDUCED_NV=5
			        PMS=$MIXED
			;;
		$MAN_CENTER_US) REDUCED_NV=6
			        PMS=$MIXED
			;;
	esac
fi

#-----------------------------------------------------------------------------------------
# Display current global constants
#-----------------------------------------------------------------------------------------
echo "#------------------------------------------#"
echo "# GLOBAL CONSTANTS:"
echo "#------------------------------------------#"
echo "OFS_ORDER  =" $OFS_ORDER
echo "OFTS_ORDER =" $OFTS_ORDER
echo "REDUCED_NV =" $REDUCED_NV
echo


#-----------------------------------------------------------------------------------------
# Display current set of parameters
#-----------------------------------------------------------------------------------------
echo "#------------------------------------------#"
echo "# PARAMETERS: "
echo "#------------------------------------------#"
# TYPE OF COMPUTATION
case $COMPTYPE in
	$QBTBP)    echo 'COMPTYPE = QBTBP'
	;;
	$NFO2) 	   echo 'COMPTYPE = NFO2'
	;;
	$PM)       echo 'COMPTYPE = PM'
	;;
	$PM_TEST)  echo 'COMPTYPE = PM_TEST'
	;;
	$COMPMAP)  echo 'COMPTYPE = COMPMAP'
	;;
	$COC)      echo 'COMPTYPE = COC'
	;;
	$DYNEQ)    echo 'COMPTYPE = DYNEQ'
	;;
	$FT_TEST)  echo 'COMPTYPE = FT_TEST'
	;;
	*) 	  echo "COMPTYPE = " $COMPTYPE ". Unknown type."
esac
	
#LIBPOINT
case $LIBPOINT in
	$EML1)   echo 'LIBPOINT = EML1'
	;;
	$EML2)   echo 'LIBPOINT = EML2'
	;;
	$SEL1)   echo 'LIBPOINT = SEL1'
	;;
	$SEL2)   echo 'LIBPOINT = SEL2'
	;;
	*) 	  echo "LIBPOINT = " $LIBPOINT ". Unknown type."
esac

# MODEL
case $MODEL in
	$CRTBP)  echo 'MODEL   = CRTBP'
	;;
	$QBCP) 	echo 'MODEL    = QBCP'
	;;
	$BCP)   echo 'MODEL    = BCP'
	;;
	$ERTBP) echo 'MODEL    = ERTBP'
	;;
	*)      echo "MODEL    = " $MODEL ". Unknown type."
esac


# NORMALIZATION
echo "ISNORM   =" $ISNORM


# DEFAULT LIBRATION POINT
echo "LI_EM    =" $LI_EM
echo "LI_SEM   =" $LI_SEM


# PM STYLE
case $PMS in
	$GRAPH)    echo 'PMS      = GRAPH'
	;;
	$NORMFORM) echo 'PMS      = NORMFORM'
	;;
	$MIXED)    echo 'PMS      = MIXED'
	;;
	*)         echo "PMS      = " $PMS ". Unknown type."
esac

# MANIFOLD TYPE (EM)
#If yes, set the REDUCED_NV in correspondance
case $MANTYPE in
	$MAN_CENTER)    echo 'MANTYPE  = MAN_CENTER'
			;;
	$MAN_CENTER_S)  echo 'MANTYPE  = MAN_CENTER_S'
			;;
	$MAN_CENTER_U)  echo 'MANTYPE  = MAN_CENTER_U'
			;;
	$MAN_CENTER_US) echo 'MANTYPE  = MAN_CENTER_US'
			;;
	*) 	        echo "MANTYPE  = " $MANTYPE". Unknown type." 
esac


# STORAGE
echo "STORAGE  =" $STORAGE
echo


#-----------------------------------------------------------------------------------------
# Display initial conditions if necessary: PM TESTING
#-----------------------------------------------------------------------------------------
if [ $COMPTYPE == $PM_TEST ]; then

	#Check the number of initial conditions match REDUCED_NV
	if [ ${#IC[@]} != $REDUCED_NV ]; then
		echo '* The number of IC does not match REDUCED_NV. Stop.'
		exit 1
	fi
	
	# Check that the OFTS_ORDERS exist
	if [ -z ${OFTS_ORDERS+x} ]; then
		TEMP=(2 5 8)
		set_param "OFTS_ORDERS" TEMP[@]
		set_param "N_OFTS_ORDERS" 3
		REMARKS+="* The variable OFTS_ORDERS, which was not defined, has been set to [${OFTS_ORDERS[*]}] .\n"
	fi


	# Check that the OFS_ORDERS exist
	if [ -z ${OFS_ORDERS+x} ]; then
		TEMP=( 30 )
		set_param "OFS_ORDERS" $TEMP
		set_param "N_OFS_ORDERS" 1
		REMARKS+="* The variable OFS_ORDERS, which was not defined, has been set to [${OFS_ORDERS[*]}].\n"
	fi

	
	# Echo
	echo "#------------------------------------------#"
	echo "# INITIAL CONDITIONS FOR PM TESTING:"
	echo "#------------------------------------------#"
	echo "si = [" ${IC[*]} "]"
	echo

	echo "#------------------------------------------#"
	echo "# OFTS OFTS_ORDERS FOR PM TESTING:"
	echo "#------------------------------------------#"
	echo "N_OFTS_ORDERS = " $N_OFTS_ORDERS
	echo "OFTS_ORDERS   = [" ${OFTS_ORDERS[*]} "]"
	echo
	
	echo "#------------------------------------------#"
	echo "# OFS OFTS_ORDERS FOR PM TESTING:"
	echo "#------------------------------------------#"
	echo "N_OFS_ORDERS = " $N_OFS_ORDERS
	echo "OFS_ORDERS   = [" ${OFS_ORDERS[*]} "]"
	echo
	
	# Check if TMIN exist
	if [ -z ${TMIN+x} ]; then
		set_param "TMIN" 0.0
		REMARKS+="* The variable TMIN, which was not defined, has been set to $TMIN .\n"
	fi
	
	# Check if TMAX exist
	if [ -z ${TMAX+x} ]; then
		set_param "TMAX" 0.5
		REMARKS+="* The variable TMIN, which was not defined, has been set to $TMAX .\n"
	fi
	
	
	echo "#------------------------------------------#"
	echo "# TIME INVERVAL FOR PM TESTING:"
	echo "#------------------------------------------#"
	echo "[ $TMIN $TMAX ]"

	echo
fi


#-----------------------------------------------------------------------------------------
# Display initial conditions if necessary: PMAP
#-----------------------------------------------------------------------------------------
if [ $COMPTYPE == $COMPMAP ]; then


	# Check if PMAP_TYPE exist
	if [ -z ${PMAP_TYPE+x} ]; then
		set_param "PMAP_TYPE" $PMAP
		REMARKS+="* The variable PMAP_TYPE, which was not defined, has been set to $PMAP_TYPE .\n"
	fi
	
	
	echo "#------------------------------------------#"
	echo "# PMAP TYPE (deprecated):"
	echo "#------------------------------------------#"
	case $PMAP_TYPE in
		$PMAP) echo 'PMAP_TYPE = PMAP'
		;;
		$TMAP) echo 'PMAP_TYPE = TMAP'
		;;
		$EMAP) echo 'PMAP_TYPE = EMAP'
		;;
		$IMAP) echo 'PMAP_TYPE = IMAP'
		;;
		$HMAP) echo 'PMAP_TYPE = HMAP'
		;;
		$IMAPPLANAR) echo 'PMAP_TYPE = IMAPPLANAR'
		;;
		*)     echo "PMAP_TYPE = " $PMAP_TYPE ". Unknown type."
	esac
	echo

	echo "#------------------------------------------#"
	echo "# PMAP PARAMETERS:"
	echo "#------------------------------------------#"
	echo "ofts_order =" $PMAP_ofts_order
	echo "ofs_order  =" $PMAP_ofs_order
	echo "f_proj_T   =" $PMAP_f_proj_T
	echo "dh_nc      =" $PMAP_dh_nc
	echo "t0_T       =" $PMAP_t0_T
	echo "tf_T       =" $PMAP_tmax_T
	echo "max_events =" $PMAP_max_events
	echo "n_sol      =" $PMAP_n_sol
	echo "si_min     =" $PMAP_si_min
	echo "si_max     =" $PMAP_si_max
	echo "is_plot    =" $PMAP_is_plot
	echo "is_par     =" $PMAP_is_par
	echo

	echo "#------------------------------------------#"
	echo "# OpenMP SETTINGS:"
	echo "#------------------------------------------#"
	echo "NUM_THREADS =" $NUM_THREADS
	echo
	
fi

#-----------------------------------------------------------------------------------------
# Print remarks
#-----------------------------------------------------------------------------------------
echo -e $REMARKS
echo "#------------------------------------------#"

#-----------------------------------------------------------------------------------------
# Go on with the implementation?
#-----------------------------------------------------------------------------------------
echo -e "Do you want to go on with the computation (y/n)? \c "
read  ans
# bash check the answer
if [ "$ans" == "y" ]; then
	
	#-------------------------------
	#Build the array of coefficients
	#-------------------------------
	case $COMPTYPE in
	
	$QBTBP | $NFO2 | $PM | $COC | $DYNEQ | $FT_TEST) 
		COEFFS=($OFTS_ORDER $OFS_ORDER $REDUCED_NV $COMPTYPE $MODEL $LIBPOINT $PMS $MANTYPE $STORAGE)
		;;

	$PM_TEST) 
		COEFFS=($OFTS_ORDER $OFS_ORDER $REDUCED_NV $COMPTYPE $MODEL $LIBPOINT $PMS $MANTYPE $STORAGE )
		COEFFS=(${COEFFS[*]} ${IC[*]} $N_OFTS_ORDERS ${OFTS_ORDERS[*]} $N_OFS_ORDERS ${OFS_ORDERS[*]} $TMIN $TMAX)
		;;
	$COMPMAP) 
		#Building the PMAP vector
		PMAP=($PMAP_TYPE $PMAP_ofts_order $PMAP_ofs_order $PMAP_dh_nc $PMAP_f_proj_T   )
 	        PMAP=(${PMAP[*]} $PMAP_t0_T $PMAP_tmax_T $PMAP_max_events  $PMAP_n_sol $PMAP_si_min $PMAP_si_max)
		PMAP=(${PMAP[*]} $PMAP_is_plot $PMAP_is_par $PMAP_int_method)
		#Number of PMAP parameters 
		NPMAP=${#PMAP[@]} 
		#Coeffs
		COEFFS=($OFTS_ORDER  $OFS_ORDER $REDUCED_NV $COMPTYPE $MODEL)
		COEFFS=(${COEFFS[*]} $LIBPOINT $PMS $MANTYPE $STORAGE)
		COEFFS=(${COEFFS[*]} ${PMAP[*]} $NUM_THREADS ) 
		;;
	esac

	#-------------------------------
	#Call software
	#-------------------------------
	#Check the NOHUP condition
	if [ "$ISNOHUP" == "1" ]; then
			#If true, check that an output file has been set
			if [ -z ${OUT+x} ]; then
				echo 'WARNING: OUT variable has not been set.'
				echo 'OUT = default.out by default.'
				OUT='default.out'
			fi
			#Call software
			nohup bin/Release/sempm ${COEFFS[*]} > $OUT &
	else
			bin/Release/sempm ${COEFFS[*]}
	fi

	
	
else  
	echo "Stop. No computation."
fi 
