# CodFleEcoEvo

Code and data used for the manuscript **Competitive traps from eco-evolutionary feedbacks hinder recovery in harvested Fish populations**.

Authors: Mikael Ohlsson\*, György Barabás, Max Lindmark, Michele Casini, Viktor Thunell, Valerio Bartolino, Mattias Sköld, and Anna Eklöf.

\*Responsible for repository: mikael.ohlsson@liu.se

The code is included in the `markdown/` folder, in separate Quarto Markdown files. Simply open and render them in RStudio to create the HTML output and generate the figures used in the manuscript. Rendering has been tested under R 4.3.3 with `tidyverse` 2.0.0.

- `markdown/cod-flounder-basic.qmd` contains the complete code for running the eco-evolutionary model. Code for generating related figures for different scenarios is also included.
- `markdown/cod-flounder-feedback.qmd` contains the complete code for running the eco-evolutionary model with cod-to-flounder feedbacks included. It reproduces Figures S1 and S2 from the Supporting Information. *Note:* rendering this file may take several minutes.
- `markdown/empirical-data.qmd` contains the code to read in data from `data/`, run analyses, and generate the figures included in the manuscript.

The empirical data are in the `data/` folder, in compressed `.rds` format. To load them, use either the base R function `readRDS`, or the `tidyverse` version `read_rds` (the markdown file `markdown/empirical-data.qmd` uses the latter).

The `manuscript/` folder contains the source files for the main text and the Supporting Information.
