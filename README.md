# Visual Search on Natural Scenes
This repository contains the code for the paper **Modeling Human Visual Search in Natural Scenes: A Combined Bayesian Searcher and Saliency Map Approach**

# Run Visual Search
To run the IBS model or any of the variants, run the *run\_main.m* script inside the dynamic_model folder (in MatlabCode folder). This script allows to set different search strategies and priors, the possible parameters (check BayesianSaliencyModel.m for more info) are:

* **Searcher**
	1. Ideal Bayesian Searcher (IBS)
	2. Structural Similarity IBS (sIBS)
	3. Correlation IBS (cIBS)
	4. Greedy 
* **Saliency Prior**
	1. DeepGaze2
	2. SAM-ResNet
	3. SAM-VGG
	4. Center
	5. Flat
	6. Noisy

For more info about the parameters check `BayesianSaliencyModel.m` and `run_main.m`.

## Requierements
This models was developed and tested on Matlab 2018b. To use any searcher the saliency prior distribution must be provided. We use the code provided by each respective paper to calculate in the case of the data presented here.

## sIBS
To use the sIBS model or repliclate the results of the paper, you can precalculate the in the working directory. 

## Paper's results
To reproduce the Frontiers' paper results checkout to the frontiers branch, donwload the data at this link and first execute `run_paper.m` to calculate all models' scanpaths and then run `run_analysis.m` to replicate figures from the results\scripts folder. Depending on your PC this may take a while (more than a day).

# Cite

Please cite using the following BibTex:

@article{10.3389/fnsys.2022.882315,
   author={Bujia, Gaston and Sclar, Melanie and Vita, Sebastian and Solovey, Guillermo and Kamienkowski, Juan Esteban},   
   title={Modeling Human Visual Search in Natural Scenes: A Combined Bayesian Searcher and Saliency Map Approach},      
   journal={Frontiers in Systems Neuroscience},
   volume={16},
   year={2022},
   url={https://www.frontiersin.org/articles/10.3389/fnsys.2022.882315},       
   doi={10.3389/fnsys.2022.882315},
   issn={1662-5137}   
}

This work was presented as [*Modeling human visual search: A combined Bayesian searcher and saliency map approach for eye movement guidance in natural scenes*](https://arxiv.org/pdf/2009.08373) in [Shared Visual Representations in Human & Machine Intelligence,
2020 NeurIPS Workshop (Dec 12, 2020)](https://www.svrhm.com/). You can find the presentation for that version at the [CBMM MIT's web](https://cbmm.mit.edu/video/modeling-human-visual-search-combined-bayesian-searcher-and-saliency-map-approach-eye-movement).
