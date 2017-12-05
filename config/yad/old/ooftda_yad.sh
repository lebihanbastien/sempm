#!/bin/sh
source 'config/parameters.sh'

# Example of lines
#--text="Please enter your details:" --text-align=center \
#--field="<b>Coordinate system:</b>:CB" \
#--field="$LABEL_COORDINATE_SYSTEM" \
#--field="<b>Earth-Moon libration point:</b>:CB" \
#--field="<b>Sun-Earth libration point:</b>:CB" \

#-----------------------------------------------------------------------------------------
# Using yad to create dialog
#-----------------------------------------------------------------------------------------
fields=$(yad --title="sempm configuration dialog" --columns="1" \
--geometry=400x800 --borders=10  \
--image="logo.svg" --image-on-top  \
--form --separator="," --item-separator="," \
--field="<b>Types of computation:</b>:CB" \
--field="$LABEL_TYPE_OF_COMPUTATION" \
--field="<b>Model:</b>:RO" \
--field="<b>Libration Point:</b>:CB" \
--field="<b>Storage:</b>:LBL" \
--field="Save the results in ./data:CHK" \
--button=gtk-cancel:1 --button=gtk-ok:2   \
--on-top \
"$LIST_TYPE_OF_COMPUTATIONS"  "" "QBCP" "$LIST_LIBRATION_POINTS" "" "TRUE")


#Echo the complete fields of answers
echo "$fields"

#-----------------------------------------------------------------------------------------
#Getting back each piece of data with awk.
#-----------------------------------------------------------------------------------------
field1=$(echo $fields | awk 'BEGIN {FS="," } { print $1 }')
fieldt=$(echo $fields | awk 'BEGIN {FS="," } { print $2 }')
field2=$(echo $fields | awk 'BEGIN {FS="," } { print $3 }')
field3=$(echo $fields | awk 'BEGIN {FS="," } { print $4 }')
fieldt=$(echo $fields | awk 'BEGIN {FS="," } { print $5 }')
field4=$(echo $fields | awk 'BEGIN {FS="," } { print $6 }')


#-----------------------------------------------------------------------------------------
#Saving the numerical values
#-----------------------------------------------------------------------------------------
COMP_TYPE=${!field1}
MODEL=${!field2}
LIBPOINT=${!field3}
STORAGE=$((field4))

#CS=${!field3}
#LI_EM=$((field4))
#LI_SEM=$((field5))

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
	
	
#Echo
echo "COMP_TYPE =" $COMP_TYPE
echo "MODEL     =" $MODEL
echo "LIBPOINT  =" $LIBPOINT
echo "STORAGE   =" $STORAGE
echo ""
echo "CS        =" $CS
echo "LI_EM     =" $LI_EM
echo "LI_SEM    =" $LI_SEM


