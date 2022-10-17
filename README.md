# RainbowGraph
Estimate underlying regeneration characteristics in rainbow Cre loxP lineage tracing using.
Given a dataset representing a tissue as a network of cells with given colours, predict what fraction of observed cells were present during Cre activation and how many among them were proliferating.
Package contains tools for:  
-efficiently simulating proliferation/regeneration in a given cell network  
-training a machine learning model on the simulated data  
-using the trained model to predict characteristics of the original data  


``` r

# package installation
# install.packages("devtools")
devtools::install_github("jmorys/RainbowGraph")

# This package relies on an R interface to Keras and Tensorflow, which have to be installed separately.
# To install them directly from R, run:
# keras::install_keras()
# Contrary to the documentation, the above may attempt to install the gpu version. To avoid that use:
# keras::install_keras(version = "cpu")

# To check package functionality, run (it will take several minutes):
test_result <- RainbowGraph::test_functionality()

# This runs the "get_complete_results_bayes" function, which contains the entire pipeline, on a test_graph supplied with the package.

```
