# VisualSearch
Repository to migrate Matlab code to Python.

# Run Ideal Bayesian Searcher in Matlab
To run the IBS or any of the variants, run the *run\_main.m* script inside the dynamic_model folder (in MatlabCode folder). This script allows to set different search strategies and priors, the options (check BayesianSaliencyModel for more info) are:

* Searcher
	1. Correlation IBS 
	2. Geisler's IBS
	3. Greedy
* Prior
	1. DeepGaze2
	2. MLNet
	3. SAM-ResNet
	4. SAM-VGG
	5. Center
	6. Flat
	7. Noisy

# Matlab Version
The model was tested on Matlab 2016 and higher.
