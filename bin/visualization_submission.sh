#!/bin/sh
module load Workspace
#####################################################################################
# First we define the general parameters of the run                                 #
#####################################################################################
SUBMISSION='visualization'

# 0 = Not a debug run, 1 = a debug run
DEBUG=1

# Which scenario would you like to run? The current options are listed below:
# 'AdvectionDiffusionOnly', 'CoastalProximity', 'Stochastic', 'ShoreDependentResuspension', 'TurrellResuspension',
# 'SizeTransport', 'FragmentationKaandorp', 'FragmentationKaandorpPartial', 'BlueCloud'
SCENARIO='BlueCloud'

#What server everything is running on
SERVER='UBELIX'
export SERVER

export SUBMISSION
export SCENARIO
export STOKES
export BACKWARD
export SERVER
export INPUT
export ADVECTION_DATA

#Now, we can set the job name prefix
if [ $SCENARIO == "AdvectionDiffusionOnly" ]; then
  RUNNAMEPREFIX="Visualization_AdvDifOnly"
elif [ $SCENARIO == "CoastalProximity" ]; then
  RUNNAMEPREFIX="Visualization_Prox"
elif [ $SCENARIO == "Stochastic" ]; then
  RUNNAMEPREFIX="Visualization_Stochastic"
elif [ $SCENARIO == "ShoreDependentResuspension" ]; then
  RUNNAMEPREFIX="Visualization_SDResus"
elif [ $SCENARIO == "TurrellResuspension" ]; then
  RUNNAMEPREFIX="Visualization_Turrell"
elif [ $SCENARIO == "FragmentationKaandorp" ]; then
  RUNNAMEPREFIX="Visualization_KaandorpFrag"
elif [ $SCENARIO == "FragmentationKaandorpPartial" ]; then
  RUNNAMEPREFIX="Visualization_PartialKaandorpFrag"
elif [ $SCENARIO == "SizeTransport" ]; then
  RUNNAMEPREFIX="Visualization_SizeTransport"
elif [ $SCENARIO == "BlueCloud" ]; then
  RUNNAMEPREFIX="Visualization_BlueCloud"
fi

echo $RUNNAMEPREFIX

#####################################################################################
# Now the part where we create the submission file                                  #
#####################################################################################
runname=$RUNNAMEPREFIX
if [ $SERVER == "UBELIX" ]; then
  #Submission for ubelix
  part1="#!/bin/sh"
  part2="#SBATCH --mail-type=fail"
  part3="#SBATCH --mail-user=victor.onink@climate.unibe.ch"
  part4="#SBATCH --job-name="$runname
  part5="#SBATCH --output="runOutput/$runname".o%j"
  part6="#SBATCH --mem-per-cpu=20G"
  if [ "$DEBUG" -eq "0" ]; then
        part7="#SBATCH --time=48:00:00"
        part8="#SBATCH --partition=epyc2"
        part9='#SBATCH --qos=job_epyc2'
  else
        part7="#SBATCH --time=00:19:00"
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
  echo ${!partGrab} >> jobsubmissionFile.sh
  done
elif [ $SERVER == "KUPHAVEN" ]; then
  #Submission for kuphaven
  part1="#!/bin/sh"
  part2="#SBATCH --mail-type=begin,end,fail"
  part3="#SBATCH --mail-user=victor.onink@climate.unibe.ch"
  part4="#SBATCH --job-name="$runname
  part5="#SBATCH --output="runOutput/$runname".o%j"
  part6="#SBATCH --mem-per-cpu=20G"
  if [ "$DEBUG" -eq "0" ]; then
    part7="#SBATCH --time=48:00:00"
    part8="#SBATCH --partition=long"
  else
    part7="#SBATCH --time=00:29:59"
    part8="#SBATCH --partition=debug"
  fi
  part9="source /storage/climatestor/Bern3dLPX/onink/alphadata04/.bash_profile"
  part10="source /storage/climatestor/Bern3dLPX/onink/alphadata04/anaconda3/bin/activate py3_parcels"
  part11='cd "/storage/climatestor/Bern3dLPX/onink/alphadata04/lagrangian_sim/BeachingSim/Next-Stage-Plastic-Beaching/"'
  part12="python src/main.py -p 10 -v"
  #and now the creation of the submission file
  for i in {1..12}
  do
  partGrab="part"$i
  echo ${!partGrab} >> jobsubmissionFile.sh
  done
fi

sbatch jobsubmissionFile.sh

rm jobsubmissionFile.sh
