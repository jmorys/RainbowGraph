# RainbowGraph
Package for estimating simple regeneration parameters of graphs representing tissue tracked with random colour Cre loxP lineage tracing.
Provides tools for simulating regeneration on a given structure, creation of machine learning model for inferring regeneration parameters using simulated data, and using created model to predict parameters for original labeling.


``` r
# package installation
# install.packages("devtools")
devtools::install_github("jmorys/RainbowGraph")

# This package relies on an R interface to Keras and Tensorflow, and as such those tools have to be installed.
# this can be achieved by running
# Keras::install_keras()
# contrary to documentation keras may attempt to install gpu version. If it causes problems during installation try
# Keras::install_keras(version = "cpu")

# to check package functionality run
test_result - RainbowGraph::test_functionality()

# this runs the get_complete_results_bayes function, which performs all steps of regeneration characteristics prediction, on a test_graph supplied in package.
# running this function is times consuming and will take several minutes

```
