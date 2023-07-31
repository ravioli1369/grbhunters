# Statistical Siginificance of Gamma Ray Bursts

## Getting Started
This repository contains the code used to perform the statistical analysis of the data from the Cadmium-Zinc-Telluride (CZT) detector on Astrosat to determine the significance of the Gamma Ray Bursts (GRBs) detected by the detector. The pipeline performs the tasks of finding the GRB as well as determining its Signal to Noise Ratio (SNR). 

### Prerequisites

The code is written in Python 3.9.1 and uses the following packages:
* [Astropy](https://www.astropy.org/)
* [Numpy](https://numpy.org/)
* [Scipy](https://www.scipy.org/)
* [Matplotlib](https://matplotlib.org/)

Note that you will require a working version of cztpipeline version 3 to run the code. The code was tested on czt_pipeline_20221209_v3.0. The pipeline can be downloaded from [here](http://astrosat-ssc.iucaa.in/cztiData).

### Running the code
The relevant codes are two python scripts: `final_script.py` and `pipelinev3.py`. The former is the main script which calls the latter. To run the code, simply run the following command in the terminal:
```
python3 final_script.py -d <path to data> -t <trigger time> -n <name of GRB> --timebin <Optional: Binsize of time series>
```
Note that the `.evt`, `.mkf`, `livetime.fits` and `mkfThresholds.txt` files must be in the same directory as given by the path to data. The trigger time is the time of the GRB in Astrosat seconds. The name of the GRB is used to name the output files (this can be anything but it's recommended to choose something sensible). The timebin is the binsize of the time series used to find the GRB and calculate the SNR. If this is not specified, the script will loop through a range of binsizes and find the one which gives the highest SNR. The final output is a `.pdf` file containing various plots and a `.txt` file containing the SNR and the significance of the GRB.


Happy Hunting :)