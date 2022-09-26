# Visual Search on Natural Scenes
This repository contains the code for the paper [**Modeling Human Visual Search in Natural Scenes: A Combined Bayesian Searcher and Saliency Map Approach**](https://www.frontiersin.org/articles/10.3389/fnsys.2022.882315/full). 

## Abstract
Finding objects is essential for almost any daily-life visual task. Saliency models have been useful to predict fixation locations in natural images during a free-exploring task. However, it is still challenging to predict the sequence of fixations during visual search. Bayesian observer models are particularly suited for this task because they represent visual search as an active sampling process. Nevertheless, how they adapt to natural images remains largely unexplored. Here, we propose a unified Bayesian model for visual search guided by saliency maps as prior information. We validated our model with a visual search experiment in natural scenes. We showed that, although state-of-the-art saliency models performed well in predicting the first two fixations in a visual search task ( 90% of the performance achieved by humans), their performance degraded to chance afterward. Therefore, saliency maps alone could model bottom-up first impressions but they were not enough to explain scanpaths when top-down task information was critical. In contrast, our model led to human-like performance and scanpaths as revealed by: first, the agreement between targets found by the model and the humans on a trial-by-trial basis; and second, the scanpath similarity between the model and the humans, that makes the behavior of the model indistinguishable from that of humans. Altogether, the combination of deep neural networks based saliency models for image processing and a Bayesian framework for scanpath integration probes to be a powerful and flexible approach to model human behavior in natural scenarios.

## Run Visual Search
To run the IBS model or any of the variants, run the `run_main.m` script inside the dynamic_model folder (in MatlabCode folder). This script allows to set different search strategies and priors, the possible parameters are:

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

To run the model individually check `BayesianSaliencyModel.m` and `main.m`.

## Requierements
These models were developed and tested on Matlab 2018b. To use any searcher the saliency prior distribution must be provided in the `MatlabCode/saliency/SALIENCY_MAP`, images and templates sould be placed on `MatlabCode/data_images/images`, `MatlabCode/data_images/templates` 
respectively. We use the code provided by each respective paper to calculate saliency models for SAM and DeepGazeII.

## sIBS
To use the sIBS model the Structural Similarity map for each image should be provided at `MatlabCode/ssim/`.

## Paper's results
To reproduce the Frontiers' paper results checkout to the frontiers branch, first execute `MatlabCode/dynamic_models/run_paper.m` to calculate all models' scanpaths and then run `run_analysis.m` to replicate figures from the `MatlabCode/results/scripts` folder. Depending on your PC this may take a while (more than a day).

## Cite

Please cite using the following BibTex:

```
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
```

This work was presented as [*Modeling human visual search: A combined Bayesian searcher and saliency map approach for eye movement guidance in natural scenes*](https://arxiv.org/pdf/2009.08373) in [Shared Visual Representations in Human & Machine Intelligence,
2020 NeurIPS Workshop (Dec 12, 2020)](https://www.svrhm.com/). You can find the presentation for that version at the [CBMM MIT's web](https://cbmm.mit.edu/video/modeling-human-visual-search-combined-bayesian-searcher-and-saliency-map-approach-eye-movement).
