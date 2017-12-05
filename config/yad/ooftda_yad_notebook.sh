#!/bin/sh
source 'config/constants.sh'


#-----------------------------------------------------------------------------------------
# Using yad to create notebook
#-----------------------------------------------------------------------------------------

##### Random Key for the notebook
id=$(echo $[($RANDOM % ($[10000 - 32000] + 1)) + 10000] )

##### Information displayed in tab 1
tooltip="<b>Type of computation:</b>\n\n$LABEL_TYPE_OF_COMPUTATION\n\n"
tooltip+="Note that the number of reduced variables can only be user-defined for FT_TEST. "
tootlip+="It is overridden in the other cases."

####################  First tab: Computation ################### 
yad --plug="$id" --tabnum=1  --form  --columns="1" --separator="|" --item-separator='|' --focus-field="1" --dialog-sep \
--text "$tooltip"  \
--field=":CB" "$LIST_TYPE_OF_COMPUTATIONS" \
--field="<b>Numerical constants:</b>:LBL" '' \
--field="Order of the Taylor series::NUM" '8!1..40!1!' \
--field="Order of the Fourier series::NUM" '30!0..30!5!' \
--field="Number of reduced variables::NUM" '4!4..6!1!' \
--field="<b>Storage:</b>:LBL" '' \
--field="Save the results:CHK" "FALSE" \
&> ./config/yad/tmp/comp &
####################  Second tab: Model ################### 
yad --plug="$id" --tabnum=2 --form --columns="1" --separator="|" --item-separator='|' --focus-field="2" --dialog-sep  \
--field="<b>Model:</b>:RO" "QBCP" \
--field="<b>Libration Point:</b>:CB" "$LIST_LIBRATION_POINTS" \
&> ./config/yad/tmp/model &
####################  Third tab: Param Method ################### 
yad --plug="$id" --tabnum=3 --form --columns="1" --separator="|" --item-separator='|' --focus-field="2" --dialog-sep  \
--field="<b>Manifold Type:</b>:LBL" '' \
--field="$LABEL_MANIFOLD_TYPE" '' \
--field=":CB" "$LIST_MANIFOLD_TYPE" \
--field="<b>Parameterization Style:</b>:LBL" '' \
--field="$LABEL_PM_STYLE" '' \
--field=":CB" "$LIST_PM_STYLE" \
&> ./config/yad/tmp/pm &
#################### Notebook ################### 
yad --notebook --width="900" --height="700" --center --title="sempm configuration dialog" --key="$id" \
--image="logo.svg" --image-on-top  \
--tab="Computation" --tab="Model" --tab="Parameterization" \
--button="gtk-cancel:1" \
--button="gtk-ok:0"

#-----------------------------------------------------------------------------------------
# Close if desired by the user
#-----------------------------------------------------------------------------------------
closing_notebook=$?

if [ $closing_notebook -eq 1 ] 
then
	echo "Closing the GUI."
	exit
else

#-----------------------------------------------------------------------------------------
# Read the tmp data with cat
#-----------------------------------------------------------------------------------------
comp=$(cat ./config/yad/tmp/comp)
model=$(cat ./config/yad/tmp/model)
pm=$(cat ./config/yad/tmp/pm)


#-----------------------------------------------------------------------------------------
#Getting back each piece of data with awk: computation
#
# Allows to update: COMPTYPE, OFTS_ORDER, OFS_ORDER, REDUCED_NV (may be overridden), STORAGE
#
#-----------------------------------------------------------------------------------------
comp1=$(echo $comp | awk 'BEGIN {FS="|" } { print $1 }')
compt=$(echo $comp | awk 'BEGIN {FS="|" } { print $2 }')  #NULL char
comp2=$(echo $comp | awk 'BEGIN {FS="|" } { print $3 }')
comp3=$(echo $comp | awk 'BEGIN {FS="|" } { print $4 }')
comp4=$(echo $comp | awk 'BEGIN {FS="|" } { print $5 }')
compt=$(echo $comp | awk 'BEGIN {FS="|" } { print $6 }')  #NULL char
comp5=$(echo $comp | awk 'BEGIN {FS="|" } { print $7 }')

#Saving the numerical values
COMPTYPE=${!comp1}
OFTS_ORDER=$((comp2))
OFS_ORDER=$((comp3))
REDUCED_NV=$((comp4))
STORAGE=$((comp5))

#-----------------------------------------------------------------------------------------
#Getting back each piece of data with awk: model
#-----------------------------------------------------------------------------------------
model1=$(echo $model | awk 'BEGIN {FS="|" } { print $1 }')
model2=$(echo $model | awk 'BEGIN {FS="|" } { print $2 }')

#Saving the numerical values
MODEL=${!model1}
LIBPOINT=${!model2}

#-----------------------------------------------------------------------------------------
#Getting back each piece of data with awk: model
#-----------------------------------------------------------------------------------------
pmt=$(echo $pm | awk 'BEGIN {FS="|" } { print $1 }')  #NULL char
pmt=$(echo $pm | awk 'BEGIN {FS="|" } { print $2 }')  #NULL char
pm1=$(echo $pm | awk 'BEGIN {FS="|" } { print $3 }') 
pmt=$(echo $pm | awk 'BEGIN {FS="|" } { print $4 }')  #NULL char
pmt=$(echo $pm | awk 'BEGIN {FS="|" } { print $5 }')  #NULL char
pm2=$(echo $pm | awk 'BEGIN {FS="|" } { print $6 }') 

#Saving the numerical values
MANTYPE=${!pm1}
PMS=${!pm2}


#-----------------------------------------------------------------------------------------
# Switch to update CS, LI_EM, LI_SEM
#-----------------------------------------------------------------------------------------
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
	
	
#-----------------------------------------------------------------------------------------
# DEBUG
#-----------------------------------------------------------------------------------------
# Raw data from tmp
# echo "comp  = " $comp
# echo "model = " $model
# echo "pm    = " $pm


#Echo final values
# echo "COMPTYPE   =" $COMPTYPE
# echo "OFTS_ORDER =" $OFTS_ORDER
# echo "OFS_ORDER  =" $OFS_ORDER
# echo "REDUCED_NV =" $REDUCED_NV
# echo ""
# echo "MODEL     =" $MODEL
# echo "LIBPOINT  =" $LIBPOINT
# echo "STORAGE   =" $STORAGE
# echo ""
# echo "CS        =" $CS
# echo "LI_EM     =" $LI_EM
# echo "LI_SEM    =" $LI_SEM
# echo ""
# echo "MANTYPE   =" $MANTYPE
# echo "PMS       =" $PMS


fi
