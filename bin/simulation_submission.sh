#####################################################################################
# General run parameters                                                            #
#####################################################################################
SUBMISSION='simulation'
DEBUG=0 # 0 = Not a debug run, 1 = a debug run
#0=first order, 1=coastal, 2=stochastic beaching/resuspension, 3=coast type dependent, 4 = Turrell (2020)
#5 = Kaandorp based fragmentation, 6 = Size dependent transport
SCENARIO=6
#for scenario 1, the time a particle must be near the coast to beach (in days)
VICINITY=2
#for scenario 2, the beaching and resuspension timescales (in days)
SHORETIME_list=(20)
RESUSTIME_list=(69)
#for scenario 3, the shore dependence scenario.
SHOREDEPEN=0
#for scenario 4, the minimum wind speed for resuspension. Divide by 10 for actual value
WMIN=3
#for scenario 5 and 6, the initial size of the particle in 1e-6 m
PARTICLE_SIZE_list=(5000 1000 500 100 90 80 70 60 50 40 30 20 10 5 1)
#for scenario 5 and 6, the critical bottom shear stress for particle resuspension (x1e-3)
SEABED_CRIT=140
# For scenario 5, the fragmentation parameters p (x1e-1), DN (x1e-1), the number of size classes and the fragmentation
# timescale (DAYS)
P=4
DN=25
SIZE_CLASS_NUMBER=6
LAMBDA_FRAG=385
#the starting year of the simulation, and how many years the simulation will take
STARTYEAR=2010
#Which input distribution do we want to use? 0=Jambeck, 1=lebreton, 2=point release, 3=uniform release
INPUT=1
#Which advection data do we want to use?
# 0 = Global HYCOM, 1 = Caribbean HYCOM, 2 = Mediterranean CMEMS
ADVECTION_DATA=2
#Start year of the simulation. 0 = new simulation, otherwise it picks up from a previous simulation
START=1
#Number of years the simulation runs
SIMLEN=3
#Inclusion of Stokes drift. 0 = include stokes, 1 = do not include stokes
STOKES=0 
#Ensemble member
ENSEMBLE=1
#Ubelix server, so server==1
SERVER=1

#Exporting the variables
export SUBMISSION
export SCENARIO
export VICINITY
export SHOREDEPEN
export WMIN
export SEABED_CRIT
export P
export DN
export SIZE_CLASS_NUMBER
export LAMBDA_FRAG
export STARTYEAR
export INPUT
export ADVECTION_DATA
export START
export SIMLEN
export STOKES
export ENSEMBLE
export SERVER

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

#####################################################################################
# Creating the submission file and submitting the job                               #
#####################################################################################
for SHORETIME in "${SHORETIME_list[@]}"; do
  export SHORETIME
  for RESUSTIME in "${RESUSTIME_list[@]}"; do
    export RESUSTIME
    for PARTICLE_SIZE in "${PARTICLE_SIZE_list[@]}"; do
      export PARTICLE_SIZE

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
        RUNNAMEPREFIX="Stochastic_ST="${SHORETIME}"_RT="${RESUSTIME}"_y="${STARTYEAR}"_"
        if [ "$STOKES" -eq "1" ]; then
            RUNNAMEPREFIX=${RUNNAMEPREFIX}"NS_"
        fi
      elif [ "$SCENARIO" -eq "3" ]; then
        RUNNAMEPREFIX="SDResus_SD="${SHOREDEPEN}"_ST="${SHORETIME}"_RT="${RESUSTIME}"_y="${STARTYEAR}"_"
        if [ "$STOKES" -eq "1" ]; then
          RUNNAMEPREFIX=${RUNNAMEPREFIX}"NS_"
        fi
      elif [ "$SCENARIO" -eq "4" ]; then
        RUNNAMEPREFIX="Turrell_Wmin="${WMIN}"_ST="${SHORETIME}"_y="${STARTYEAR}"_"
        if [ "$STOKES" -eq "1" ]; then
          RUNNAMEPREFIX=${RUNNAMEPREFIX}"NS_"
        fi
      elif [ "$SCENARIO" -eq "5" ]; then
        RUNNAMEPREFIX="KaandorpFrag_ST="${SHORETIME}"_RT="${RESUSTIME}"_y="${STARTYEAR}"_"
        if [ "$STOKES" -eq "1" ]; then
          RUNNAMEPREFIX=${RUNNAMEPREFIX}"NS_"
        fi
      elif [ "$SCENARIO" -eq "6" ]; then
        RUNNAMEPREFIX="SizeTransport_SIZE="${PARTICLE_SIZE}"_ST="${SHORETIME}"_y="${STARTYEAR}"_tau="${SEABED_CRIT}"_"
        if [ "$STOKES" -eq "1" ]; then
          RUNNAMEPREFIX=${RUNNAMEPREFIX}"NS_"
        fi
      fi

      RUNNAMEPREFIX=${RUNNAMEPREFIX}"ENSEMBLE="${ENSEMBLE}"_"
      echo $RUNNAMEPREFIX

      #Looping over all the runs based on the input scenario
      for ((RUN=0; RUN<=$runlength; RUN++)); do
        export RUN
        # looping over all the simulation years
        for ((RESTARTNUM=$START; RESTARTNUM<$SIMLEN; RESTARTNUM++)); do
           export RESTARTNUM
           runname=$RUNNAMEPREFIX"_run="$RUN"_restart="$RESTARTNUM
           part1="#!/bin/sh"
           part2="#SBATCH --mail-type=begin,end,fail"
           part3="#SBATCH --mail-user=victor.onink@climate.unibe.ch"
           part4="#SBATCH --job-name="$runname
           part5="#SBATCH --output="runOutput/$runname".o%j"
           part6="#SBATCH --mem-per-cpu=40G"
           if [ "$DEBUG" -eq "0" ]; then
            part7="#SBATCH --time=95:59:00"
            part8="#SBATCH --partition=epyc2"
            part9='#SBATCH --qos=job_epyc2'
           else
            part7="#SBATCH --time=00:10:00"
            part8="#SBATCH --partition=epyc2"
            part9='#SBATCH --qos=job_epyc2_debug'
           fi
           part10="source /storage/homefs/vo18e689/.bash_profile"
           part11="source /storage/homefs/vo18e689/anaconda3/bin/activate py3_parcels"
           part12='cd "/storage/homefs/vo18e689/codes/Next-Stage-Plastic-Beaching/"'
           part13="python src/main.py -p 10 -v"

           #and now the creation of the submission file
           for i in {1..13}; do
              partGrab="part"$i
              echo ${!partGrab} >> jobsubmissionFile_${RUN}_${RESTARTNUM}.sh
           done

           #Submitting the job
           if [ "$RESTARTNUM" -eq "$START" ]; then
              jobid=$(sbatch --parsable jobsubmissionFile_${RUN}_${RESTARTNUM}.sh)
           else
              jobid=$(sbatch --parsable --dependency=afterok:${jobid} jobsubmissionFile_${RUN}_${RESTARTNUM}.sh)
           fi

           #Removing the used submission file
           rm jobsubmissionFile_${RUN}_${RESTARTNUM}.sh
        done
      done
    done
  done
done