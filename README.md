# CodFleEcoEvo
Code and data used for the manuscript **Competitive traps from eco-evolutionary feedbacks hinder recovery in harvested Fish populations**.

Authors: Mikael Ohlsson\*, György Barabás, Max Lindmark, Michele Casini, Viktor Thunell, Valerio Bartolino, Mattias Sköld, Anna Eklöf.

\*Responsible for repository: mikael.ohlsson@liu.se

Before cloning the repository, make sure to have Git LFS (Large File Storage) installed (`git lfs install`). If you have already cloned the repository without Git LFS installed, install it as previously and run `git lfs pull`.
 
The code is included in two separate R-markdown files. Simply open the respective files in RStudio and knit them to render the HTML output and generate the figures used in the manuscript.

`R/ecoevo.Rmd` contains the complete code for running the eco-evolutionary model. Code for generating related figures for different scenarios is also included.

`R/empirical_data.Rmd` contains the code to read in used data from `data/`, run analyses and generate the figures included in the manuscript.
