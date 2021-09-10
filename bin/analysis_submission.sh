#!/bin/sh
while IFS= read -r line; do
    echo "Text read from file: $line"
done < "job_id.txt"

######################################################################################
## First we define the general parameters of the run                                 #
######################################################################################
#SUBMISSION='analysis'
#export SUBMISSION
#DEBUG=0 # 0 = Not a debug run, 1 = a debug run
##0=first order, 1=coastal, 2=stochastic beaching/resuspension, 3=coast type dependent, 4 = Turrell (2020)
##5 = Size dependent transport, 6 = Kaandorp based fragmentation, 7 = alternate Kaandorp fragmentation
#SCENARIO=7
#export SCENARIO
##for scenario 1, the time a particle must be near the coast to beach (in days)
#VICINITY=2
#export VICINITY
##for scenario 2, the beaching and resuspension timescales (in days)
#SHORETIME_list=(20)
#RESUSTIME_list=(69)
##for scenario 3, the shore dependence scenario.
#SHOREDEPEN=0
#export SHOREDEPEN
##for scenario 4, the minimum wind speed for resusplension. Divide by 10 for actual value
#WMIN=3
#export WMIN
##for scenario 5 and 6, the initial size of the particle in 1e-6 m and the particle rho
#PARTICLE_SIZE_list=(5000)
#INIT_DENSITY=920
#export INIT_DENSITY
##for scenario 5 and 6, the critical bottom shear stress for particle resuspension (x1e-3)
#SEABED_CRIT=0
#export SEABED_CRIT
## For scenario 5, the fragmentation parameters p (x1e-1), DN (x1e-1), the number of size classes and fragmentation
## timescale (days). Also, the option of including ocean fragmentation or not
#P=4
#DN=25
#SIZE_CLASS_NUMBER=6
#LAMBDA_FRAG_list=(388)
#OCEAN_FRAG=0
#LAMBDA_OCEAN_FRAG_LIST=(388)
#export P
#export DN
#export SIZE_CLASS_NUMBER
#export OCEAN_FRAG
##(1 10 100 200 300 388)
##the starting year of the simulation, and how many years the simulation will take
#YEAR=2010
#STARTMONTH_list=(1 2 3 4 5 6 7 8 9 10 11 12)
#STARTDAY=1
#export STARTDAY
##Which input distribution do we want to use? 0=Jambeck, 1=lebreton, 2=lebretondivision, 3=point, 4=uniform
#INPUT=2
#export INPUT
##Which advection data do we want to use?
## 0 = Global HYCOM, 1 = Caribbean HYCOM, 2 = Mediterranean CMEMS
#ADVECTION_DATA=2
#export ADVECTION_DATA
##Number of years the simulation runs
#SIMLEN=2
#export SIMLEN
## For the analysis, if we have multiple years that we want to combine into one analysis set (so if we have continuous
## particle release), then this how many years we want to include
#if [ "$SCENARIO" -eq "7" ]; then
#  COMBINE_YEARS=${SIMLEN}
#else
#  COMBINE_YEARS=1
#fi
##Inclusion of Stokes drift. 0 = include stokes, 1 = do not include stokes
#STOKES=0
#export STOKES
##Ensemble member
#ENSEMBLE=1
#export ENSEMBLE
##Ubelix server, so server==1
#SERVER=1
#export SERVER
#
##A number of switches to indicate which analysis steps we want to run.
##0 = off, 1 = on
#CONCENTRATION=1
#VERTICAL_CONCENTRATION=0
#TIMESERIES=0
#MAX_DISTANCE=0
#TIMESLICING=0
#STATISTICS=0
#SEPARATION=0
#SIZE_SPECTRUM=0
#
#export CONCENTRATION
#export VERTICAL_CONCENTRATION
#export TIMESERIES
#export MAX_DISTANCE
#export TIMESLICING
#export STATISTICS
#export SEPARATION
#export SIZE_SPECTRUM
#
######################################################################################
## A number of tests to make sure that I don't submit too many jobs by accident      #
######################################################################################
#if [ "$SCENARIO" -eq "5" -a ${#LAMBDA_FRAG_list[@]} -gt 1 ]; then
#  echo 'For SizeTransport, do not submit multiple LAMBDA_FRAG values.'
#  exit
#fi
#if [ "$SCENARIO" -eq "5" -a ${#LAMBDA_OCEAN_FRAG_LIST[@]} -gt 1 ]; then
#  echo 'For SizeTransport, do not submit multiple LAMBDA_OCEAN_FRAG values.'
#  exit
#fi
#if [ "$SCENARIO" -eq "7" -a ${#PARTICLE_SIZE_list[@]} -gt 1 ]; then
#  echo 'For KaandorpFragmentationPartial, please only submit one PARTICLE_SIZE.'
#  exit
#fi
#if [ "$SCENARIO" -eq "7" -a ${#LAMBDA_OCEAN_FRAG_LIST[@]} -gt 1 -a $OCEAN_FRAG -eq 0 ]; then
#  echo 'Without OCEAN_FRAG, do not submit more than one LAMBDA_OCEAN_FRAG value'
#  exit
#fi
#
##The number of runs we do, dependent on the input scenario.
#if [ "$INPUT" -eq "0" ]; then
#  runlength=0 #8
#elif [ "$INPUT" -eq "1" ]; then
#  runlength=0 #3
#elif [ "$INPUT" -eq "2" ]; then
#  runlength=9
#elif [ "$INPUT" -eq "3" ]; then
#  runlength=0
#fi
#
######################################################################################
## Now the part where we create the submission file and submit the job               #
######################################################################################
#
#for SHORETIME in "${SHORETIME_list[@]}"; do
#  export SHORETIME
#  for RESUSTIME in "${RESUSTIME_list[@]}"; do
#    export RESUSTIME
#    for PARTICLE_SIZE in "${PARTICLE_SIZE_list[@]}"; do
#      export PARTICLE_SIZE
#      for LAMBDA_FRAG in "${LAMBDA_FRAG_list[@]}"; do
#        export LAMBDA_FRAG
#        for LAMBDA_OCEAN_FRAG in "${LAMBDA_OCEAN_FRAG_LIST[@]}"; do
#          export LAMBDA_OCEAN_FRAG
#          for STARTMONTH in "${STARTMONTH_list[@]}"; do
#            export STARTMONTH
#
#            #Now, we can set the job name prefix
#            if [ "$SCENARIO" -eq "0" ]; then
#              RUNNAMEPREFIX="Analysis_AdvDifOnly_y="${YEAR}"_"
#            elif [ "$SCENARIO" -eq "1" ]; then
#              RUNNAMEPREFIX="Analysis_Prox_vic="${VICINITY}"_y="${YEAR}"_"
#            elif [ "$SCENARIO" -eq "2" ]; then
#              RUNNAMEPREFIX="Analysis_Stochastic_ST="${SHORETIME}"_RT="${RESUSTIME}"_y="${YEAR}"_"
#            elif [ "$SCENARIO" -eq "3" ]; then
#              RUNNAMEPREFIX="Analysis_SDResus_SD="${SHOREDEPEN}"_ST="${SHORETIME}"_RT="${RESUSTIME}"_y="${YEAR}"_"
#            elif [ "$SCENARIO" -eq "4" ]; then
#              RUNNAMEPREFIX="Analysis_Turrell_Wmin="${WMIN}"_ST="${SHORETIME}"_y="${YEAR}"_"
#            elif [ "$SCENARIO" -eq "6" ]; then
#              RUNNAMEPREFIX="Analysis_KaandorpFrag_ST="${SHORETIME}"_RT="${RESUSTIME}"_y="${YEAR}"_"
#            elif [ "$SCENARIO" -eq "7" ]; then
#              RUNNAMEPREFIX="Analysis_PartialKaandorpFrag_ST="${SHORETIME}"_RT="${RESUSTIME}"_y="${YEAR}"_lamf="${LAMBDA_FRAG}"_lamfO="${LAMBDA_FRAG_LAMBDA_OCEAN_FRAG}"_"
#            elif [ "$SCENARIO" -eq "5" ]; then
#              RUNNAMEPREFIX="Analysis_SizeTransport_SIZE="${PARTICLE_SIZE}"_ST="${SHORETIME}"_y="${YEAR}"_"
#            fi
#            if [ "$STOKES" -eq "1" ]; then
#              RUNNAMEPREFIX=${RUNNAMEPREFIX}"NS_"
#            fi
#            RUNNAMEPREFIX=${RUNNAMEPREFIX}"ENSEMBLE="${ENSEMBLE}"_"
#
#            #Initializing a string used for keeping track of the job dependencies
##            JOB_TRACKER=''
#            PARALLEL_STEP=1
#            export PARALLEL_STEP
#            #First we are going to submit all the jobs for the individual run/restart files, so each runs the analysis
#            #code just for that specific parcels output file. We also need to consider the various starting years
#            RESTART_REMOVE=0
#            for ((STARTYEAR=${YEAR}; STARTYEAR<$((YEAR+COMBINE_YEARS)); STARTYEAR++)); do
#              export STARTYEAR
#              for ((RUN=0; RUN<=$runlength; RUN++)); do
#                export RUN
#                # looping over all the simulation years
#                for ((RESTARTNUM=0; RESTARTNUM<$((SIMLEN-RESTART_REMOVE)); RESTARTNUM++)); do
#                  export RESTARTNUM
#                  # specifying the parts of the submission file
#                  part1="#!/bin/sh"
#                  part2="#SBATCH --mail-type=begin,end,fail"
#                  part3="#SBATCH --mail-user=victor.onink@climate.unibe.ch"
#                  part4="#SBATCH --job-name="${RUNNAMEPREFIX}
#                  part5="#SBATCH --output="runOutput/${RUNNAMEPREFIX}".o%j"
#                  part6="#SBATCH --mem-per-cpu=20G"
#                  if [ "$DEBUG" -eq "0" ]; then
#                    part7="#SBATCH --time=01:30:00"
#                    part8="#SBATCH --partition=epyc2"
#                    part9='#SBATCH --qos=job_epyc2'
#                  else
#                    part7="#SBATCH --time=00:19:00"
#                    part8="#SBATCH --partition=epyc2"
#                    part9='#SBATCH --qos=job_epyc2_debug'
#                  fi
#                  part10="source /storage/homefs/vo18e689/.bash_profile"
#                  part11="source /storage/homefs/vo18e689/anaconda3/bin/activate py3_parcels"
#                  part12='cd "/storage/homefs/vo18e689/codes/Next-Stage-Plastic-Beaching/"'
#                  part13="python src/main.py -p 10 -v"
#
#                  # Putting all the parts into the submission file
#                  for i in {1..13}; do
#                    partGrab="part"$i
#                    echo ${!partGrab} >> jobsubmissionFile_${RUN}_${RESTARTNUM}.sh
#                  done
#
#                  # submitting the job
#                  jobid=$(sbatch --parsable jobsubmissionFile_${RUN}_${RESTARTNUM}.sh)
#                  echo -n ':'${jobid} >> job_id.txt
#                  #JOB_TRACKER=${JOB_TRACKER}':'${jobid}
#                  #JOB_TRACKER[${#JOB_TRACKER[@]}]=${jobid}
#                  scancel ${jobid}
#                  # deleting the submission file
#                  rm jobsubmissionFile_${RUN}_${RESTARTNUM}.sh
#                done
#              done
#              RESTART_REMOVE=$((RESTART_REMOVE+1))
#            done
#          done
#
#          # Remove character of the JOB_TRACKER so that we don't end with :
#          PARALLEL_STEP=2
#          export PARALLEL_STEP
#          STARTYEAR=${YEAR}
#          export STARTYEAR
#          # specifying the parts of the submission file
#          part1="#!/bin/sh"
#          part2="#SBATCH --mail-type=begin,end,fail"
#          part3="#SBATCH --mail-user=victor.onink@climate.unibe.ch"
#          part4="#SBATCH --job-name="${RUNNAMEPREFIX}
#          part5="#SBATCH --output="runOutput/${RUNNAMEPREFIX}".o%j"
#          part6="#SBATCH --mem-per-cpu=20G"
#          if [ "$DEBUG" -eq "0" ]; then
#            part7="#SBATCH --time=01:30:00"
#            part8="#SBATCH --partition=epyc2"
#            part9='#SBATCH --qos=job_epyc2'
#          else
#            part7="#SBATCH --time=00:19:00"
#            part8="#SBATCH --partition=epyc2"
#            part9='#SBATCH --qos=job_epyc2_debug'
#          fi
#          part10="source /storage/homefs/vo18e689/.bash_profile"
#          part11="source /storage/homefs/vo18e689/anaconda3/bin/activate py3_parcels"
#          part12='cd "/storage/homefs/vo18e689/codes/Next-Stage-Plastic-Beaching/"'
#          part13="python src/main.py -p 10 -v"
#
#          # Putting all the parts into the submission file
#          for i in {1..13}; do
#            partGrab="part"$i
#            echo ${!partGrab} >> jobsubmissionFile.sh
#          done
#
#          # Submitting the job that will join all the various analysis files together
##          sbatch --dependency=afterok${JOB_TRACKER} jobsubmissionFile.sh
#          rm jobsubmissionFile.sh
#        done
#      done
#    done
#  done
#done
