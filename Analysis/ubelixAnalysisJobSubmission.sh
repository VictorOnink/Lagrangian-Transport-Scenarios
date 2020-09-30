#This will be a general code that will allow easy submission of large beaching simulation jobs. It is supposed to work in unison with generalBeachingScenarios.py, which is the python code that contains the actually parcels to run the different simulations

#####################################################################################
# First we define the general parameters of the run                                 #
#####################################################################################
#0=first order, 1=coastal, 2=stochastic beaching/resuspension, 3=coast type dependent
SCENARIO=3
export SCENARIO
#for scenario 1, the time a particle must be near the coast to beach (in days)
VICINITY=2
export VICINITY
#for scenario 2, the beaching and resuspension timescales (in days)
SHORETIME=26
export SHORETIME
RESUSTIME=69
export RESUSTIME
#For scenario 3, we need to indicate if beaching is more likely with sand or not sand. 0 = more sand is less likely beaching and resuspension, 1 = more sand is more likely beaching and resuspension
SHOREDEPEN=0
export SHOREDEPEN
#the starting year of the simulation, and how many years the simulation will take
STARTYEAR=2010
export STARTYEAR
SIMLEN=5 #Number of years the simulation runs
export SIMLEN
STOKES=0 #0 = include stokes, 1 = do not include stokes
export STOKES
INPUT=0 #0=Jambeck, higher values aren't relevant yet but can be added
export INPUT

ENSEMBLE=1 #For when we want to run multiple iterations of the same parameter setup
export ENSEMBLE

#Now, we can set the job name prefix

if [ "$SCENARIO" -eq "0" ]; then
	RUNNAMEPREFIX="Analysis_AdvDifOnly_y="${STARTYEAR}"_"
	if [ "$STOKES" -eq "1" ]; then
	    RUNNAMEPREFIX=${RUNNAMEPREFIX}"NS_"
        fi
elif [ "$SCENARIO" -eq "1" ]; then
	RUNNAMEPREFIX="Analysis_Prox_vic="${VICINITY}"_y="${STARTYEAR}"_"
	if [ "$STOKES" -eq "1" ]; then
	    RUNNAMEPREFIX=${RUNNAMEPREFIX}"NS_"
        fi
elif [ "$SCENARIO" -eq "2" ]; then
	RUNNAMEPREFIX="Analysis_Stochastic_ST="${SHORETIME}"_RT="${RESUSTIME}"_y"${STARTYEAR}"_"
	if [ "$STOKES" -eq "1" ]; then
	    RUNNAMEPREFIX=${RUNNAMEPREFIX}"NS_"
        fi
elif [ "$SCENARIO" -eq "3" ]; then
        RUNNAMEPREFIX="Analysis_SDResus_SD="${shoreDepen}"_ST="${SHORETIME}"_RT="${RESUSTIME}"_y"${STARTYEAR}"_"
        if [ "$STOKES" -eq "1" ]; then
	    RUNNAMEPREFIX=${RUNNAMEPREFIX}"NS_"
        fi
elif [ "$SCENARIO" -eq "4" ]; then
        RUNNAMEPREFIX="Analysis_Turrell_Wmin="${WMIN}"_ST="${SHORETIME}"_RT="${RESUSTIME}"_y"${STARTYEAR}"_"
        if [ "$STOKES" -eq "1" ]; then
	    RUNNAMEPREFIX=${RUNNAMEPREFIX}"NS_"
        fi
fi


RUNNAMEPREFIX=${RUNNAMEPREFIX}"ENSEMBLE="${ENSEMBLE}
echo $RUNNAMEPREFIX
#####################################################################################
# Now the part where we create the submission file                                  #
#####################################################################################
runname=$RUNNAMEPREFIX
part1="#!/bin/sh"
part2="#SBATCH --mail-type=begin,end,fail"
part3="#SBATCH --mail-user=victor.onink@climate.unibe.ch"
part4="#SBATCH --job-name="$runname
part5="#SBATCH --output="runOutput/$runname".o%j"m
part6="#SBATCH --mem-per-cpu=6G"
part7="#SBATCH --time=00:10:00"
part8='#SBATCH --partition=debug'
#loading the bash and setting the environment
part9="source /home/ubelix/climate/vo18e689/.bash_profile"
part10="source /home/ubelix/climate/vo18e689/anaconda3/bin/activate py3_parcels_v2_2"
part11='cd "/home/ubelix/climate/vo18e689/codes/Next-Stage-Plastic-Beaching/Analysis"'
#And now the actual running of the code
part12="python masterProcessing.py -p 10 -v"
#and now the creation of the submission file
for i in {1..12} 
do
partGrab="part"$i
echo ${!partGrab} >> jobsubmissionFile_${run}_${restartnum}.sh
done

sbatch jobsubmissionFile_${run}_${restartnum}.sh

rm jobsubmissionFile_${run}_${restartnum}.sh
