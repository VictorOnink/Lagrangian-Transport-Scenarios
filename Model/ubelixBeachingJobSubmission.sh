#This will be a general code that will allow easy submission of large beaching simulation jobs. It is supposed to work in unison with generalBeachingScenarios.py, which is the python code that contains the actually parcels to run the different simulations

#####################################################################################
# General run parameters                                                            #
#####################################################################################
#0=first order, 1=coastal, 2=stochastic beaching/resuspension, 3=coast type dependent, 4=Turrell (2020)
SCENARIO=4
#for scenario 1, the time a particle must be near the coast to beach (in days)
VICINITY=2
#for scenario 2, the beaching and resuspension timescales (in days)
SHORETIME=10
RESUSTIME=69
#for scenario 3, the shore dependence scenario.
SHOREDEPEN=0
#for scenario 4, the minimum wind speed for resusplension. Divide by 10 for actual value
WMIN=3
#the starting year of the simulation, and how many years the simulation will take
STARTTIME=2010
#Which input distribution do we want to use? 0=Jambeck, 1=lebreton
INPUT=0
#Start year of the simulation. 0 = new simulation, otherwise it picks up from a previous simulation
START=0 
#Number of years the simulation runs
SIMLEN=1
#Inclusion of Stokes drift. 0 = include stokes, 1 = do not include stokes
STOKES=0 
#Ensemble member
ENSEMBLE=1

export SCENARIO
export VICINITY
export ENSEMBLE
export SHORETIME
export RESUSTIME
export SHOREDEPEN
export WMIN
export STARTTIME
export INPUT
export START
export SIMLEN
export STOKES
export ENSEMBLE

#Setting the name of the job
if [ "$SCENARIO" -eq "0" ]; then
	RUNNAMEPREFIX="AdvDifOnly_y="${STARTTIME}"_"
	if [ "$STOKES" -eq "1" ]; then
	    RUNNAMEPREFIX=${RUNNAMEPREFIX}"NS_"
        fi
elif [ "$SCENARIO" -eq "1" ]; then
	RUNNAMEPREFIX="Prox_vic="${VICINITY}"_y="${STARTTIME}"_"
	if [ "$STOKES" -eq "1" ]; then
	    RUNNAMEPREFIX=${RUNNAMEPREFIX}"NS_"
        fi
elif [ "$SCENARIO" -eq "2" ]; then
	RUNNAMEPREFIX="Stochastic_ST="${SHORETIME}"_RT="${RESUSTIME}"_y"${STARTTIME}"_"
	if [ "$STOKES" -eq "1" ]; then
	    RUNNAMEPREFIX=${RUNNAMEPREFIX}"NS_"
        fi
elif [ "$SCENARIO" -eq "3" ]; then
        RUNNAMEPREFIX="SDResus_SD="${shoreDepen}"_ST="${SHORETIME}"_RT="${RESUSTIME}"_y"${STARTTIME}"_"
        if [ "$STOKES" -eq "1" ]; then
	    RUNNAMEPREFIX=${RUNNAMEPREFIX}"NS_"
        fi
elif [ "$SCENARIO" -eq "4" ]; then
        RUNNAMEPREFIX="Turrell_Wmin="${WMIN}"_ST="${SHORETIME}"_RT="${RESUSTIME}"_y"${STARTTIME}"_"
        if [ "$STOKES" -eq "1" ]; then
	    RUNNAMEPREFIX=${RUNNAMEPREFIX}"NS_"
        fi
fi


RUNNAMEPREFIX=${RUNNAMEPREFIX}"ENSEMBLE="${ENSEMBLE}
echo $RUNNAMEPREFIX
#####################################################################################
# Now the part where we create the submission file                                  #
#####################################################################################

#The number of runs we do, dependent on the input scenario.
if [ "$INPUT" -eq "0" ]; then
    runlength=8
elif [ "$INPUT" -eq "1" ]; then
    runlength=3
fi

#Create a loop to run over all 
for ((RUN=0; RUN<=$runlength;RUN++))
do
	export RUN           
	for ((RESTARTNUM=$START; RESTARTNUM<$SIMLEN; RESTARTNUM++))
	do
	   export RESTARTNUM
	   runname=$RUNNAMEPREFIX"_run="$RUN"_restart="$RESTARTNUM
	   part1="#!/bin/sh"
	   part2="#SBATCH --mail-type=begin,end,fail"
	   part3="#SBATCH --mail-user=victor.onink@climate.unibe.ch"
	   part4="#SBATCH --job-name="$runname
	   part5="#SBATCH --output="runOutput/$runname".o%j"
	   part6="#SBATCH --mem-per-cpu=6G"
           part7="#SBATCH --time=00:10:00"
           part8="#SBATCH --partition=debug"
	   #loading the bash and setting the environment
	   part9="source /home/ubelix/climate/vo18e689/.bash_profile"
	   part10="source /home/ubelix/climate/vo18e689/anaconda3/bin/activate py3_parcels_v2_2"
	   part11='cd "/home/ubelix/climate/vo18e689/codes/Next-Stage-Plastic-Beaching/Model/"'
	   #And now the actual running of the code
	   part12="python generalBeachingScenarios.py -p 10 -v"
	   #and now the creation of the submission file
	   for i in {1..12} 
	   do
		partGrab="part"$i
		echo ${!partGrab} >> jobsubmissionFile_${RUN}_${RESTARTNUM}.sh
	   done
	   if [ "$RESTARTNUM" -eq "$START" ]; then
   	      jobid=$(sbatch --parsable jobsubmissionFile_${RUN}_${RESTARTNUM}.sh)
           else
	      jobid=$(sbatch --parsable --dependency=afterok:${jobid} jobsubmissionFile_${RUN}_${RESTARTNUM}.sh)
	   fi
	   rm jobsubmissionFile_${RUN}_${RESTARTNUM}.sh
	done
done
