# Astronomical Diagrams project
The project consisted of the analysis of astronomical diagrams from manuscripts dating from the 12th century. 
The database is not available on GitHub due to its weight, however I can email it to you if needed.
## Content of the project
The objective was to implement a method for extracting and clustering the diagrams to facilitate their analysis by astronomers. 
The extraction of diagrams from manuscripts was previously done thanks to a segmentation algorithm using DocExtrator. 
Then a convolutional neural network architecture has been implemented with PyTorch to analyse the deep features of the diagrams. 
To generate the diagram clusters, several dimension reduction and clustering methods were tested (GaussianMixtures, etc..).
Those retained here are PCA and Kmeans. 
The data visualization was done via Matplotlib and Gephi software allowing a better visualization of the network.
