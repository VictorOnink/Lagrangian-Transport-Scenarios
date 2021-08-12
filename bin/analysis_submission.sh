#####################################################################################
# First we define the general parameters of the run                                 #
#####################################################################################
SUBMISSION='analysis'
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
#for scenario 4, the minimum wind speed for resusplension. Divide by 10 for actual value
WMIN=3
#for scenario 5 and 6, the initial size of the particle in 1e-6 m
PARTICLE_SIZE_list=(5000 1000 500 100 90 80 70 60 50 40 30 20 10 5 1)
#for scenario 5 and 6, the critical bottom shear stress for particle resuspension (x1e-3)
SEABED_CRIT=140
# For scenario 5, the fragmentation parameters p (x1e-1), DN (x1e-1), the number of size classes and fragmentation
# timescale (days)
P=4
DN=25
SIZE_CLASS_NUMBER=6
LAMBDA_FRAG_list=(388)
#the starting year of the simulation, and how many years the simulation will take
STARTYEAR=2010
#Which input distribution do we want to use? 0=Jambeck, 1=lebreton, 2=point release, 3=uniform release
INPUT=1
#Which advection data do we want to use?
# 0 = Global HYCOM, 1 = Caribbean HYCOM, 2 = Mediterranean CMEMS
ADVECTION_DATA=2
#Number of years the simulation runs
SIMLEN=1
#Inclusion of Stokes drift. 0 = include stokes, 1 = do not include stokes
STOKES=0
#Ensemble member
ENSEMBLE=1
#Ubelix server, so server==1
SERVER=1

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
export SIMLEN
export STOKES
export INPUT
export ENSEMBLE
export STARTYEAR
export ADVECTION_DATA
export SERVER

#A number of switches to indicate which analysis steps we want to run.
#0 = off, 1 = on
CONCENTRATION=0
VERTICAL_CONCENTRATION=0
TIMESERIES=0
MAX_DISTANCE=0
TIMESLICING=0
STATISTICS=0
SEPARATION=1
SIZE_SPECTRUM=0

export CONCENTRATION
export VERTICAL_CONCENTRATION
export TIMESERIES
export MAX_DISTANCE
export TIMESLICING
export STATISTICS
export SEPARATION
export SIZE_SPECTRUM

#####################################################################################
# Now the part where we create the submission file and submit the job               #
#####################################################################################

for SHORETIME in "${SHORETIME_list[@]}"; do
  export SHORETIME
  for RESUSTIME in "${RESUSTIME_list[@]}"; do
    export RESUSTIME
    for PARTICLE_SIZE in "${PARTICLE_SIZE_list[@]}"; do
      export PARTICLE_SIZE
        for LAMBDA_FRAG in "${LAMBDA_FRAG_list[@]}"; do
          export  LAMBDA_FRAG

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
            RUNNAMEPREFIX="Analysis_Stochastic_ST="${SHORETIME}"_RT="${RESUSTIME}"_y="${STARTYEAR}"_"
            if [ "$STOKES" -eq "1" ]; then
              RUNNAMEPREFIX=${RUNNAMEPREFIX}"NS_"
            fi
          elif [ "$SCENARIO" -eq "3" ]; then
            RUNNAMEPREFIX="Analysis_SDResus_SD="${SHOREDEPEN}"_ST="${SHORETIME}"_RT="${RESUSTIME}"_y="${STARTYEAR}"_"
            if [ "$STOKES" -eq "1" ]; then
              RUNNAMEPREFIX=${RUNNAMEPREFIX}"NS_"
            fi
          elif [ "$SCENARIO" -eq "4" ]; then
            RUNNAMEPREFIX="Analysis_Turrell_Wmin="${WMIN}"_ST="${SHORETIME}"_y="${STARTYEAR}"_"
            if [ "$STOKES" -eq "1" ]; then
              RUNNAMEPREFIX=${RUNNAMEPREFIX}"NS_"
            fi
          elif [ "$SCENARIO" -eq "5" ]; then
            RUNNAMEPREFIX="Analysis_KaandorpFrag_ST="${SHORETIME}"_RT="${RESUSTIME}"_y="${STARTYEAR}"_"
            if [ "$STOKES" -eq "1" ]; then
              RUNNAMEPREFIX=${RUNNAMEPREFIX}"NS_"
            fi
          elif [ "$SCENARIO" -eq "6" ]; then
            RUNNAMEPREFIX="Analysis_SizeTransport_SIZE="${PARTICLE_SIZE}"_ST="${SHORETIME}"_y="${STARTYEAR}"_tau="${SEABED_CRIT}"_"
            if [ "$STOKES" -eq "1" ]; then
              RUNNAMEPREFIX=${RUNNAMEPREFIX}"NS_"
            fi
          fi
          RUNNAMEPREFIX=${RUNNAMEPREFIX}"ENSEMBLE="${ENSEMBLE}"_"
          echo $RUNNAMEPREFIX

          # specifying the parts of the submission file
          part1="#!/bin/sh"
          part2="#SBATCH --mail-type=begin,end,fail"
          part3="#SBATCH --mail-user=victor.onink@climate.unibe.ch"
          part4="#SBATCH --job-name="${RUNNAMEPREFIX}
          part5="#SBATCH --output="runOutput/${RUNNAMEPREFIX}".o%j"
          part6="#SBATCH --mem-per-cpu=20G"
          if [ "$DEBUG" -eq "0" ]; then
            part7="#SBATCH --time=04:00:00"
            part8="#SBATCH --partition=epyc2"
            part9='#SBATCH --qos=job_epyc2'
          else
            part7="#SBATCH --time=00:19:00"
            part8="#SBATCH --partition=epyc2"
            part9='#SBATCH --qos=job_epyc2_debug'
          fi
          part10="source /storage/homefs/vo18e689/.bash_profile"
          part11="source /storage/homefs/vo18e689/anaconda3/bin/activate py3_parcels"
          part12='cd "/storage/homefs/vo18e689/codes/Next-Stage-Plastic-Beaching/"'
          part13="python src/main.py -p 10 -v"

          # Putting all the parts into the submission file
          for i in {1..13}; do
            partGrab="part"$i
            echo ${!partGrab} >>jobsubmissionFile.sh
          done

          # submitting the job
          sbatch jobsubmissionFile.sh

          # deleting the submission file
          rm jobsubmissionFile.sh
      done
    done
  done
done
