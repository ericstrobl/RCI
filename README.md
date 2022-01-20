# Root Causal Inference (RCI)

RCI is an algorithm that discovers the root causes of disease in a patient-specific manner. The algorithm extracts mutually independent errors from the data assuming that the causal relations can be summarized using a linear structural equation model. The algorithm then outputs patient-specific Shapley values of the log-odds using logistic regression.

# Installation

> library(devtools)

> install_github("ericstrobl/RCI")

> library(RCI)

# Toy Example

Generate a random DAG with 10 variables and an expected neighborhood size of 2:

> G = generate_DAG(ps[p],en=2)

Create a dataset of 1000 samples from the DAG:

> data = sample_DAG_Y(1000,G)

> X = data$data[,-X$Y]; Y = data$data[,X$Y]

Run the RCI algorithm and print the Shapley values. Note that the column names correspond to the variable number in X:

> out = RCI(X,Y)

> print(out$scores)


