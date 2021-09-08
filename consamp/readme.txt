it is assumed that you are using anaconda for the environment management.

from root folder run

1. create a new environment
conda create -n constrained_sampling

2. active new environment
conda activate constrained_sampling

3. install all dependencies
conda install --file requirements.txt

4. run app (from constrained-sampling/consamp/plotting)
python app_3d.py

5. click link in the console output to open a browser and view the app



folder/experiment naming:
-globalsampling_globaloptimization_localsampling_localoptimization_additionalinfo
-e.g uniform_biased_RRT_biased_no_filter

files:
-summary.xlsx (metrics for 3d examples on sphere) 
