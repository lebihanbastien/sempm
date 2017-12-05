#-----------------------------------------------------------------------------------------
# Creating data folder tree architecture
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Functions
#-----------------------------------------------------------------------------------------
function build_generic_folder
{
	#QBCP
	mkdir -p data/$1/QBCP/EM/L1
	mkdir -p data/$1/QBCP/EM/L2
	
	mkdir -p data/$1/QBCP/SEM/L1
	mkdir -p data/$1/QBCP/SEM/L2
	
	# CRTBP
	mkdir -p data/$1/CRTBP/EM/L1
	mkdir -p data/$1/CRTBP/EM/L2
	
	mkdir -p data/$1/CRTBP/SEM/L1
	mkdir -p data/$1/CRTBP/SEM/L2
}

function build_pm_subfolder
{
	mkdir -p data/$1/$2/$3/$4/DWf
	mkdir -p data/$1/$2/$3/$4/FW
	mkdir -p data/$1/$2/$3/$4/Potential
	mkdir -p data/$1/$2/$3/$4/rvf
	mkdir -p data/$1/$2/$3/$4/Test
	mkdir -p data/$1/$2/$3/$4/W
	mkdir -p data/$1/$2/$3/$4/Wdot
}


function build_pm_folder
{
	build_pm_subfolder $1 QBCP EM L1
	build_pm_subfolder $1 QBCP EM L2
	build_pm_subfolder $1 QBCP SEM L1
	build_pm_subfolder $1 QBCP SEM L2
	
	build_pm_subfolder $1 CRTBP EM L1
	build_pm_subfolder $1 CRTBP EM L2
	build_pm_subfolder $1 CRTBP SEM L1
	build_pm_subfolder $1 CRTBP SEM L2
}


#-----------------------------------------------------------------------------------------
# data folder
#-----------------------------------------------------------------------------------------
# QBTBP
mkdir -p data/qbtbp

# COC
build_generic_folder COC

#VF
build_generic_folder VF

#Param method
build_pm_folder CS
build_pm_folder CU
build_pm_folder CUS
build_pm_folder MS
build_pm_folder NF
build_pm_folder CM

#-----------------------------------------------------------------------------------------
# print folder
#-----------------------------------------------------------------------------------------
#QBCP
mkdir -p print/QBCP/EM/L1/orbits
mkdir -p print/QBCP/EM/L2/orbits

mkdir -p print/QBCP/SEM/L1/orbits
mkdir -p print/QBCP/SEM/L2/orbits

# CRTBP
mkdir -p print/CRTBP/EM/L1/orbits
mkdir -p print/CRTBP/EM/L2/orbits

mkdir -p print/CRTBP/SEM/L1/orbits
mkdir -p print/CRTBP/SEM/L2/orbits
