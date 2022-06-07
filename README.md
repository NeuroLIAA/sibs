# Modeling Human Visual Search in Natural Scenes: A Combined Bayesian Searcher and Saliency Map Approach

This repository contains the code for the paper [Modeling Human Visual Search in Natural Scenes: A Combined Bayesian Searcher and Saliency Map Approach](https://www.frontiersin.org/articles/10.3389/fnsys.2022.882315/full). Previously the [preprint](https://arxiv.org/pdf/2009.08373) was presented in [Shared Visual Representations in Human & Machine Intelligence,
2020 NeurIPS Workshop (Dec 12, 2020)](https://www.svrhm.com/)

# Environment
The model was tested on Matlab 2016 and higher.

# Run Searcher in Matlab
To run the IBS model or any of the variants, run the *run\_main.m* script inside the dynamic_model folder (in MatlabCode folder). This script allows to set different search strategies and priors, the possible parameters (check BayesianSaliencyModel.m for more info) are:

* **Searcher**
	1. Structural Similarity IBS 	
	1. Correlation IBS 
	2. Geisler's IBS
	3. Greedy
* **Prior**
	1. DeepGazII
	2. MLNet
	3. SAM-ResNet
	4. SAM-VGG
	5. Center
	6. Flat
	7. Noisy

For more info about the parameters check *BayesianSaliencyModel.m*.

## Replicate paper results

## Example

# Cite us

Please cite with the following BibTeX: 

```
@article{Bujia2022vsearch,
author = {Bujia, Gaston and Sclar, Melanie and Vita, Sebastian and Solovey, Guillermo and Kamienkowski, Juan Esteban},
doi = {10.3389/fnsys.2022.882315},
issn = {1662-5137},
journal = {Frontiers in Systems Neuroscience},
title = {{Modeling Human Visual Search in Natural Scenes: A Combined Bayesian Searcher and Saliency Map Approach}},
url = {https://www.frontiersin.org/article/10.3389/fnsys.2022.882315},
volume = {16},
year = {2022}
}
```
