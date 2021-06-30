# bin 
This contains bash scripts used to submit simulation, analysis and visualization jobs 
on ubelix. Within the bash script it is possible to set simulation parameters such as
beaching timescales, input scenarios and simulation scenarios, and these parameters are
then loaded into the ```src/settings.py``` file.

The files are:
- `simulation_submission.sh`: This submits the jobs for the Lagrangian parcels simulations
- `analysis_submission.sh`: This submits jobs for post-processing of the parcels output files
- `visualization_submission.sh`: This submits jobs for making figures of the post-processed
parcels output.