# ucn2png
Convert ucn to png and shapefile by hpham@intera.com
Last updated: 03/15/2022

INSTALLATION
1. Install Anaconda Individual Edition https://repo.anaconda.com/archive/Anaconda3-2021.11-Windows-x86_64.exe
2. Install Anaconda to your Windows PC
3. Setup a new conda environment
- Open Anaconda Prompt (Anaconda3)
- Check list of current conda environments installed in your PC:
   conda env list
- Install a new conda environment
   conda create -n flopy python=3
- Activate the just created conda environment
   conda activate flopy
- Install some libraries needed for the script (all libs at once to avoid conflict)
   conda install -c conda-forge flopy geopandas numpy matplotlib pandas
4. Or you can install a new conda env using flopy.yml included in this repo
   conda env create -f flopy.yml
   
   
RUN SCRIPT
python main.py input/input.csv

