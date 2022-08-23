# VisualSearch
This repository contains the code for the paper [*Modeling human visual search: A combined Bayesian searcher and saliency map approach for eye movement guidance in natural scenes*](https://arxiv.org/pdf/2009.08373) presented in [Shared Visual Representations in Human & Machine Intelligence,
2020 NeurIPS Workshop (Dec 12, 2020)](https://www.svrhm.com/)

Please cite with the following BibTeX: 

@article{sclar2020modeling,
  title={Modeling human visual search: A combined Bayesian searcher and saliency map approach for eye movement guidance in natural scenes},
  author={Sclar, Melanie and Bujia, Gast{\'o}n and Vita, Sebasti{\'a}n and Solovey, Guillermo and Kamienkowski, Juan Esteban},
  journal={arXiv preprint arXiv:2009.08373},
  year={2020}
}

# Environment
The model was tested on Matlab 2016 and higher.

# Run Visual Search
To run the IBS model or any of the variants, run the *run\_main.m* script inside the dynamic_model folder (in MatlabCode folder). This script allows to set different search strategies and priors, the possible parameters (check BayesianSaliencyModel.m for more info) are:

* **Searcher**
	1. Ideal Bayesian Searcher (IBS)
	2. Structural Similarity IBS (sIBS)
	3. Correlation IBS (cIBS)
	4. Greedy 
* **Saliency Prior**
	1. DeepGaze2
	2. MLNet
	3. SAM-ResNet
	4. SAM-VGG
	5. Center
	6. Flat
	7. Noisy

For more info about the parameters check *BayesianSaliencyModel.m*.

## Requierements
To use any searcher the saliency prior distribution must be provided. We use the code provided by each respective paper to calculate in the case of the data presented here. Results can be found at TBD. 

### sIBS
To use the sIBS model or repliclate the results of the paper, in the working directory. 

# Paper's results
To reproduce the Frontiers' paper results checkout to the frontiers branch, and first execute `run_paper.m` to calculate all models' scanpaths and then run `run_analysis.m` to replicate figures from the results\scripts folder. Depending on your PC this may take a while (more than a day).

# TODO MATLAB CODE
- Add analysis code
- Update readme
- Add images to gitignore
