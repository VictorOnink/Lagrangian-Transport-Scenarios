#####################################################################################
# General run parameters                                                            #
#####################################################################################
SUBMISSION='simulation'
DEBUG=0 # 0 = Not a debug run, 1 = a debug run
#0=first order, 1=coastal, 2=stochastic beaching/resuspension, 3=coast type dependent, 4 = Turrell (2020)
#5 = Cozar based fragmentation, 6 = Size dependent transport
SCENARIO=6
#for scenario 1, the time a particle must be near the coast to beach (in days)
VICINITY=2
#for scenario 2, the beaching and resuspension timescales (in days)
SHORETIME=20
RESUSTIME=69
#for scenario 3, the shore dependence scenario.
SHOREDEPEN=0
#for scenario 4, the minimum wind speed for resusplension. Divide by 10 for actual value
WMIN=3
#for scenario 6, the initial size of the particle in 1e-6 m
PARTICLE_SIZE=5
#for scenario 6, the critical bottom shear stress for particle resuspension (x1e-3)
SEABED_CRIT=140
#the starting year of the simulation, and how many years the simulation will take
STARTYEAR=2010
#Which input distribution do we want to use? 0=Jambeck, 1=lebreton, 2=point release, 3=uniform release
INPUT=1
#Which advection data do we want to use?
# 0 = Global HYCOM, 1 = Caribbean HYCOM, 2 = Mediterranean CMEMS
ADVECTION_DATA=2
#Start year of the simulation. 0 = new simulation, otherwise it picks up from a previous simulation
START=0
#Number of years the simulation runs
SIMLEN=3
#Inclusion of Stokes drift. 0 = include stokes, 1 = do not include stokes
STOKES=0 
#Ensemble member
ENSEMBLE=1
#Ubelix server, so server==1
SERVER=1

export SUBMISSION
export SCENARIO
export VICINITY
export SHORETIME
export RESUSTIME
export SHOREDEPEN
export WMIN
export PARTICLE_SIZE
export SEABED_CRIT
export STARTYEAR
export INPUT
export ADVECTION_DATA
export START
export SIMLEN
export STOKES
export ENSEMBLE
export SERVER

#Setting the name of the job
if [ "$SCENARIO" -eq "0" ]; then
	RUNNAMEPREFIX="AdvDifOnly_y="${STARTYEAR}"_"
	if [ "$STOKES" -eq "1" ]; then
	    RUNNAMEPREFIX=${RUNNAMEPREFIX}"NS_"
        fi
elif [ "$SCENARIO" -eq "1" ]; then
	RUNNAMEPREFIX="Prox_vic="${VICINITY}"_y="${STARTYEAR}"_"
	if [ "$STOKES" -eq "1" ]; then
	    RUNNAMEPREFIX=${RUNNAMEPREFIX}"NS_"
        fi
elif [ "$SCENARIO" -eq "2" ]; then
	RUNNAMEPREFIX="Stochastic_ST="${SHORETIME}"_RT="${RESUSTIME}"_y"${STARTYEAR}"_"
	if [ "$STOKES" -eq "1" ]; then
	    RUNNAMEPREFIX=${RUNNAMEPREFIX}"NS_"
        fi
elif [ "$SCENARIO" -eq "3" ]; then
        RUNNAMEPREFIX="SDResus_SD="${SHOREDEPEN}"_ST="${SHORETIME}"_RT="${RESUSTIME}"_y"${STARTYEAR}"_"
        if [ "$STOKES" -eq "1" ]; then
	    RUNNAMEPREFIX=${RUNNAMEPREFIX}"NS_"
        fi
elif [ "$SCENARIO" -eq "4" ]; then
        RUNNAMEPREFIX="Turrell_Wmin="${WMIN}"_ST="${SHORETIME}"_y"${STARTYEAR}"_"
        if [ "$STOKES" -eq "1" ]; then
	    RUNNAMEPREFIX=${RUNNAMEPREFIX}"NS_"
        fi
elif [ "$SCENARIO" -eq "5" ]; then
        RUNNAMEPREFIX="KaandorpFrag_ST="${SHORETIME}"_RT="${RESUSTIME}"_y"${STARTYEAR}"_"
        if [ "$STOKES" -eq "1" ]; then
	    RUNNAMEPREFIX=${RUNNAMEPREFIX}"NS_"
        fi
elif [ "$SCENARIO" -eq "6" ]; then
        RUNNAMEPREFIX="SizeTransport_SIZE="${PARTICLE_SIZE}"_ST="${SHORETIME}"_RT="${RESUSTIME}"_y"${STARTYEAR}"_"
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
    runlength=0 #8
elif [ "$INPUT" -eq "1" ]; then
    runlength=0 #3
elif [ "$INPUT" -eq "2" ]; then
    runlength=0
elif [ "$INPUT" -eq "3" ]; then
    runlength=0
fi

#Create a loop to run over all 
for ((RUN=0; RUN<=$runlength; RUN++))
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
	   part6="#SBATCH --mem-per-cpu=10G"
	   if [ "$DEBUG" -eq "0" ]; then
	          part7="#SBATCH --time=95:59:00"
            part8="#SBATCH --partition=epyc2"
            part9='#SBATCH --qos=job_epyc2'
     else
            part7="#SBATCH --time=00:10:00"
            part8="#SBATCH --partition=epyc2"
            part9='#SBATCH --qos=job_epyc2_debug'
     fi
	   #loading the bash and setting the environment
	   part10="source /storage/homefs/vo18e689/.bash_profile"
	   part11="source /storage/homefs/vo18e689/anaconda3/bin/activate py3_parcels"
	   part12='cd "/storage/homefs/vo18e689/codes/Next-Stage-Plastic-Beaching/"'
	   #And now the actual running of the code
	   part13="python src/main.py -p 10 -v"
	   #and now the creation of the submission file
	   for i in {1..13}
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
