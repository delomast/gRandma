# gRandma
grandparentage inference </br>
</br> </br> 
Inference of grandparentage based on genetic markers </br>

Citation for grandparent inference/error estimation algorithms: Delomas, T.A. and Campbell, M.R. 2022. 
Grandparent inference from genetic data: The potential for parentage-based tagging programs to 
identify offspring of hatchery strays. North American Journal of Fisheries Management. https://doi.org/10.1002/nafm.10714

Citation for single parent inference/error estimation algorithms: Steele, C.A., et al. 2022. Single-parentage 
assignments reveal negative-assortative mating in an endangered salmonid. Ecology and Evolution. https://doi.org/10.1002/ece3.8846

Potentially some other relationships as well, but main goal is grandparentage </br>
Written with mixed stock analyses in mind (i.e. descendants of multiple populations are mixed together at sampling) </br>
Performs inference and error estimation </br>
Current options: </br>
* "single-sided" grandparentage: a trio of grandchild + both maternal grandparents (or both paternal grandparents) </br>
* single parentage: a pair of parent + offspring </br>

Install with:
```
devtools::install_github("delomast/gRandma")
```

To install and view the vignette:
```
devtools::install_github("delomast/gRandma", build_vignettes = TRUE)
browseVignettes("gRandma")
```

