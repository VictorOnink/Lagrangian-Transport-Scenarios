#!/bin/sh
module load Workspace
#####################################################################################
# General run parameters                                                            #
#####################################################################################
SUBMISSION='analysis'
export SUBMISSION

# 0 = Not a debug run, 1 = a debug run
DEBUG=0

# Which scenario would you like to run? The current options are listed below:
# 'AdvectionDiffusionOnly', 'CoastalProximity', 'Stochastic', 'ShoreDependentResuspension', 'TurrellResuspension',
# 'SizeTransport', 'FragmentationKaandorp', 'FragmentationKaandorpPartial', 'BlueCloud'
SCENARIO='BlueCloud'
export SCENARIO

#for scenario 1, the time a particle must be near the coast to beach (in days)
VICINITY=2
export VICINITY

#for scenario 2, the beaching and resuspension timescales (in days)
SHORETIME_list=(26)
RESUSTIME_list=(69)

#for scenario 3, the shore dependence scenario.
SHOREDEPEN=0
export SHOREDEPEN

#for scenario 4, the minimum wind speed for resusplension. Divide by 10 for actual value
WMIN=3
export WMIN

#for scenario 5 and 6, the initial size of the particle in 1e-6 m and the particle rho
PARTICLE_SIZE_list=(5000)  #(5000 2500 1250 625 313 156 78 39 20 10 5 2)
INIT_DENSITY_list=(920) # (30 920 980 1020)
export INIT_DENSITY

#for scenario 5 and 6, the critical bottom shear stress for particle resuspension
SEABED_CRIT='0'
export SEABED_CRIT

#for scenarios 5 - 7, if fixed == 1 we have size-independent resuspension, otherwise the resuspension timescale is a
#function of the particle rise velocity
FIXED_RESUS=0
export FIXED_RESUS

# For scenario 7, the fragmentation parameters p (x1e-1), DN (x1e-1), the number of size classes and fragmentation
# timescale (days). Also, the option of including ocean fragmentation or not
P='0.4'
DN='2.5'
SIZE_CLASS_NUMBER=6
LAMBDA_FRAG_list=(388)
OCEAN_FRAG=0
LAMBDA_OCEAN_FRAG_LIST=(388)
export P
export DN
export SIZE_CLASS_NUMBER
export OCEAN_FRAG

# For scenario 8 & 9, which release site do we want to do the backwards simulation for
RELEASE_SITE=0
export RELEASE_SITE

#the starting year of the simulation, how many years the simulation will take, and whether the simulation is forwards
#or backwards in time
YEAR=2012
STARTMONTH_list=(1)
STARTDAY=1
BACKWARD=1
export STARTDAY
export BACKWARD

# Which input distribution do we want to use? 0=Jambeck, 1=lebreton, 2=lebretondivision, 3=lebretonKaandorpInit,
# 4=point, 5=uniform
INPUT=3
export INPUT

#Which advection data do we want to use?
# 0 = Global HYCOM, 1 = Caribbean HYCOM, 2 = Mediterranean CMEMS
ADVECTION_DATA=2
export ADVECTION_DATA

#Number of years the simulation runs
SIMLEN=3
export SIMLEN

# For the analysis, if we have multiple years that we want to combine into one analysis set (so if we have continuous
# particle release), then this how many years we want to include
if [ "$SCENARIO" -eq "7" ]; then
  COMBINE_YEARS=${SIMLEN}
else
  COMBINE_YEARS=1
fi

#Inclusion of Stokes drift. 0 = include stokes, 1 = do not include stokes
STOKES=0
export STOKES

#Ensemble member
ENSEMBLE=1
export ENSEMBLE

#What server everything is running on
SERVER='UBELIX'
export SERVER

#A number of switches to indicate which analysis steps we want to run.
#0 = off, 1 = on
CONCENTRATION=0
MONTHLY_CONCENTRATION=0
VERTICAL_CONCENTRATION=0
SPATIAL_VERTICAL_PROFILES=0
LONLAT_CONCENTRATION=0
TIMESERIES=0
MAX_DISTANCE=0
TIMESLICING=1
STATISTICS=0
SEPARATION=0
SIZE_SPECTRUM=0
BAYESIAN=0

export CONCENTRATION
export MONTHLY_CONCENTRATION
export VERTICAL_CONCENTRATION
export LONLAT_CONCENTRATION
export TIMESERIES
export MAX_DISTANCE
export TIMESLICING
export STATISTICS
export SEPARATION
export SIZE_SPECTRUM
export BAYESIAN
export SPATIAL_VERTICAL_PROFILES

#####################################################################################
# A number of tests to make sure that I don't submit too many jobs by accident      #
#####################################################################################
if [ $SCENARIO == "SizeTransport" -a ${#LAMBDA_FRAG_list[@]} -gt 1 ]; then
  echo 'For SizeTransport, do not submit multiple LAMBDA_FRAG values.'
  exit
fi
if [ $SCENARIO == "SizeTransport" -a ${#LAMBDA_OCEAN_FRAG_LIST[@]} -gt 1 ]; then
  echo 'For SizeTransport, do not submit multiple LAMBDA_OCEAN_FRAG values.'
  exit
fi
if [ $SCENARIO == "SizeTransport" -a  "$INPUT" -ne "1" ]; then
  echo 'For SizeTransport, make sure to use input scenario 1!!!!.'
  exit
fi
if [ $SCENARIO == "FragmentationKaandorpPartial" -a ${#PARTICLE_SIZE_list[@]} -gt 1 ]; then
  echo 'For KaandorpFragmentationPartial, please only submit one PARTICLE_SIZE.'
  exit
fi
if [ $SCENARIO == "FragmentationKaandorpPartial" -a ${#LAMBDA_OCEAN_FRAG_LIST[@]} -gt 1 -a $OCEAN_FRAG -eq 0 ]; then
  echo 'Without OCEAN_FRAG, do not submit more than one LAMBDA_OCEAN_FRAG value'
  exit
fi
if [ $SCENARIO != "FragmentationKaandorpPartial" -a  "$POST_PROCESS" -ne "0" ]; then
  echo 'Postprocessing only applies for KaandorpFragmentationPartial, not any other scenario'
  exit
fi

#The number of runs we do, dependent on the input scenario.
if [ "$INPUT" -eq "0" ]; then
  RUNLENGTH=0
elif [ "$INPUT" -eq "1" ]; then
  RUNLENGTH=0
elif [ "$INPUT" -eq "2" ]; then
  RUNLENGTH=9
elif [ "$INPUT" -eq "3" ]; then
  RUNLENGTH=9
elif [ "$INPUT" -eq "4" ]; then
  RUNLENGTH=0
elif [ "$INPUT" -eq "5" ]; then
  RUNLENGTH=0
fi
#####################################################################################
# Now the part where we create the submission file and submit the job               #
#####################################################################################

for SHORETIME in "${SHORETIME_list[@]}"; do
  export SHORETIME
  for RESUSTIME in "${RESUSTIME_list[@]}"; do
    export RESUSTIME
    for PARTICLE_SIZE in "${PARTICLE_SIZE_list[@]}"; do
      export PARTICLE_SIZE
      for INIT_DENSITY in "${INIT_DENSITY_list[@]}"; do
        export INIT_DENSITY
        for LAMBDA_FRAG in "${LAMBDA_FRAG_list[@]}"; do
          export LAMBDA_FRAG
          for LAMBDA_OCEAN_FRAG in "${LAMBDA_OCEAN_FRAG_LIST[@]}"; do
            export LAMBDA_OCEAN_FRAG
            for STARTMONTH in "${STARTMONTH_list[@]}"; do
              export STARTMONTH

              #Now, we can set the job name prefix
              if [ $SCENARIO == "AdvectionDiffusionOnly" ]; then
                RUNNAMEPREFIX="AdvDifOnly_y="${STARTYEAR}"_"
              elif [ $SCENARIO == "CoastalProximity" ]; then
                RUNNAMEPREFIX="Prox_vic="${VICINITY}"_y="${STARTYEAR}"_"
              elif [ $SCENARIO == "Stochastic" ]; then
                RUNNAMEPREFIX="Stochastic_ST="${SHORETIME}"_RT="${RESUSTIME}"_y="${STARTYEAR}"_"
              elif [ $SCENARIO == "ShoreDependentResuspension" ]; then
                RUNNAMEPREFIX="SDResus_SD="${SHOREDEPEN}"_ST="${SHORETIME}"_RT="${RESUSTIME}"_y="${STARTYEAR}"_"
              elif [ $SCENARIO == "TurrellResuspension" ]; then
                RUNNAMEPREFIX="Turrell_Wmin="${WMIN}"_ST="${SHORETIME}"_y="${STARTYEAR}"_"
              elif [ $SCENARIO == "FragmentationKaandorp" ]; then
                RUNNAMEPREFIX="KaandorpFrag_ST="${SHORETIME}"_RT="${RESUSTIME}"_y="${STARTYEAR}"_"
              elif [ $SCENARIO == "FragmentationKaandorpPartial" ]; then
                RUNNAMEPREFIX="PartialKaandorpFrag_ST="${SHORETIME}"_RT="${RESUSTIME}"_y="${STARTYEAR}"_"
                if [ "$POST_PROCESS" -eq "1" ]; then
                  RUNNAMEPREFIX='PP_'${RUNNAMEPREFIX}
                fi
              elif [ $SCENARIO == "SizeTransport" ]; then
                RUNNAMEPREFIX="SizeTransport_SIZE="${PARTICLE_SIZE}"_ST="${SHORETIME}"_y="${STARTYEAR}"_tau="${SEABED_CRIT}"_"
              elif [ $SCENARIO == "BlueCloud" ]; then
                RUNNAMEPREFIX="BlueCloud_ST="${SHORETIME}"_RT="${RESUSTIME}"_y="${STARTYEAR}"_"
              fi
              echo $RUNNAMEPREFIX

              #Initializing a string used for keeping track of the job dependencies
              PARALLEL_STEP=1
              export PARALLEL_STEP
              #First we are going to submit all the jobs for the individual run/restart files, so each runs the analysis
              #code just for that specific parcels output file. We also need to consider the various starting years
              RESTART_REMOVE=0
              for ((STARTYEAR=${YEAR}; STARTYEAR<$((YEAR+COMBINE_YEARS)); STARTYEAR++)); do
                export STARTYEAR
                for ((RUN=0; RUN<=${RUNLENGTH}; RUN++)); do
                  export RUN
                  # looping over all the simulation years
                  for ((RESTARTNUM=0; RESTARTNUM<$((SIMLEN-RESTART_REMOVE)); RESTARTNUM++)); do
                    export RESTARTNUM
                    # specifying the parts of the submission file
                    part1="#!/bin/sh"
                    part2="#SBATCH --mail-type=fail"
                    part3="#SBATCH --mail-user=victor.onink@climate.unibe.ch"
                    part4="#SBATCH --job-name="${RUNNAMEPREFIX}
                    part5="#SBATCH --output="runOutput/${RUNNAMEPREFIX}".o%j"
                    part6="#SBATCH --mem-per-cpu=20G"
                    if [ "$DEBUG" -eq "0" ]; then
                      part7="#SBATCH --time=02:00:00"
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
                      echo ${!partGrab} >> jobsubmissionFile_${RUN}_${RESTARTNUM}.sh
                    done

                    # submitting the job
                    jobid=$(sbatch --parsable jobsubmissionFile_${RUN}_${RESTARTNUM}.sh)
                    echo -n ':'${jobid} >> job_id.txt
                    # deleting the submission file
                    rm jobsubmissionFile_${RUN}_${RESTARTNUM}.sh
                  done
                done
                RESTART_REMOVE=$((RESTART_REMOVE+1))
              done
            done

            # Remove character of the JOB_TRACKER so that we don't end with :
            PARALLEL_STEP=2
            export PARALLEL_STEP
            STARTYEAR=${YEAR}
            export STARTYEAR
            RUN=${RUNLENGTH}
            RESTARTNUM=0
            export RUN
            export RESTARTNUM
            # specifying the parts of the submission file
            part1="#!/bin/sh"
            part2="#SBATCH --mail-type=fail"
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
              echo ${!partGrab} >> jobsubmissionFile.sh
            done

            # Submitting the job that will join all the various analysis files together
            dependence=$( cat job_id.txt )
            sbatch --dependency=afterok${dependence} jobsubmissionFile.sh
  #          sbatch jobsubmissionFile.sh
            rm jobsubmissionFile.sh
            rm job_id.txt
          done
        done
      done
    done
  done
done
