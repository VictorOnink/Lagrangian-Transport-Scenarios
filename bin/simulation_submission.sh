#!/bin/sh
module load Workspace
#####################################################################################
# General run parameters                                                            #
#####################################################################################
SUBMISSION='simulation'
export SUBMISSION

# 0 = Not a debug run, 1 = a debug run
DEBUG=1

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

#for scenario 4, the minimum wind speed for resuspension.
WMIN='0.3'
export WMIN

#for scenario 5 and 6, the initial size of the particle in 1e-6 m and the rho of the particle
PARTICLE_SIZE_list=(5000) # (5000 2500 1250 625 313 156 78 39 20 10 5 2)
INIT_DENSITY_list=(920) # (30 920 980 1020)

#for scenarios 5 - 7, the critical bottom shear stress for particle resuspension
SEABED_CRIT='0'
export SEABED_CRIT

#for scenarios 5 - 7, if fixed == 1 we have size-independent resuspension, otherwise the resuspension timescale is a
#function of the particle rise velocity
FIXED_RESUS=0
export FIXED_RESUS

# For scenario 7, the fragmentation parameters p (x1e-1), DN (x1e-1), the number of size classes and the fragmentation
# timescale (DAYS). Also, the option of including ocean fragmentation or not. Finally, the sink removal rate
P='0.4'
DN='2.5'
SIZE_CLASS_NUMBER=6
LAMBDA_FRAG_list=(388)  # (388 1000 10000 35000 50000)
OCEAN_FRAG=0
LAMBDA_OCEAN_FRAG_LIST=(388)
export P
export DN
export SIZE_CLASS_NUMBER
export OCEAN_FRAG

# For scenario 7, we are either running a lagrangian simulation OR we are running postprocessing on the parcels output
# to calculate the particle numbers.
# POST_PROCESS == 0 -> run lagrangian simulation
# POST_PROCESS == 1 -> run post processing
POST_PROCESS=0
export POST_PROCESS

# For scenario 8, which release site do we want to do the backwards simulation for
RELEASE_SITE=0
export RELEASE_SITE

#the starting year of the simulation, how many years the simulation will take, and whether the simulation is forwards
#or backwards in time
STARTYEAR=2012
STARTMONTH_list=(1)
STARTDAY=1
BACKWARD=1
export STARTYEAR
export STARTDAY
export BACKWARD

# Which input distribution do we want to use? 0=Jambeck, 1=lebreton, 2=lebretondivision, 3=lebretonKaandorpInit,
# 4=point, 5=uniform
INPUT=4
export INPUT

#Which advection data do we want to use?
# 0 = Global HYCOM, 1 = Caribbean HYCOM, 2 = Mediterranean CMEMS
ADVECTION_DATA=2
export ADVECTION_DATA

#Start year of the simulation. 0 = new simulation, otherwise it picks up from a previous simulation
START=0

#Number of years the simulation runs
SIMLEN=1
export SIMLEN

#Inclusion of Stokes drift. 0 = include stokes, 1 = do not include stokes
STOKES=0
export STOKES

#Ensemble member
ENSEMBLE=1
export ENSEMBLE

#What server everything is running on
SERVER='UBELIX'
export SERVER

#The number of runs we do, dependent on the input scenario.
if [ "$INPUT" -eq "0" ]; then
  runlength=0
elif [ "$INPUT" -eq "1" ]; then
  runlength=0
elif [ "$INPUT" -eq "2" ]; then
  runlength=9
elif [ "$INPUT" -eq "3" ]; then
  runlength=9
elif [ "$INPUT" -eq "4" ]; then
  runlength=0
elif [ "$INPUT" -eq "5" ]; then
  runlength=0
fi

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


#####################################################################################
# Creating the submission file and submitting the job                               #
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

              #Setting the name of the job
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
                   if [ "$DEBUG" -eq "0" ]; then
                    if [ $POST_PROCESS -eq 0 ]; then
                     part6="#SBATCH --mem-per-cpu=40G"
                     part7="#SBATCH --time=95:59:00"
                    else
                     part6="#SBATCH --mem-per-cpu=40G"
                     part7="#SBATCH --time=95:59:00"
                    fi
                    part8="#SBATCH --partition=epyc2"
                    part9='#SBATCH --qos=job_epyc2'
                   else
                    part6="#SBATCH --mem-per-cpu=40G"
                    part7="#SBATCH --time=00:20:00"
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
      done
    done
  done
done
