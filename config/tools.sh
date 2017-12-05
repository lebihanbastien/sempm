#--------------------------------------------------------------
# Function to update unset parameters
#--------------------------------------------------------------
function set_param
{
	# argAry2 contains the arguments after $1, 
	# which may correspond to either an array or a scalar
	declare -a argAry2=("${!2}")
	local size=${#argAry2[@]}
    
	
	# If the size of argAry2 is equal to 1, we the second argument is a scalar.
	# if not, it is an array and we use a specific call
	if [ "$size" == "1" ]; then
	
		#Evaluate
		eval "$1=$2"
		
		# Warn the user
		#echo "#------------------------------------------#"
		#echo 'WARNING: the variable' $1 'is not set.'
		#echo $2 'is used by default.'
		#echo "#------------------------------------------#"
	
	
	else
	
		#Evaluate
		eval "$1=( ${argAry2[@]} )"  
		
		# Warn the user
		#echo "#------------------------------------------#"
		#echo 'WARNING: the variable' $1 'is not set.'
		#echo '(' ${argAry2[@]} ') is used by default.'
		#echo "#------------------------------------------#"
	fi

}

#--------------------------------------------------------------
# Function to compare software versions
#--------------------------------------------------------------
function version_gt() { test "$(printf '%s\n' "$@" | sort -V | head -n 1)" != "$1"; }
