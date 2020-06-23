# gRandma
grandparentage inference </br>
</br> </br> 
Inference of grandparentage based on genetic markers </br>
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

