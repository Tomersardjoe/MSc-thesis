<!-- README based on template: https://github.com/othneildrew/Best-README-Template -->

<a id="readme-top"></a>

<!-- PROJECT SHIELDS -->
[![Unlicense License][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url]



<!-- PROJECT HEADER -->
<br />
<div align="center">

  <h1 align="center">MSc thesis repository - Tomer Sardjoe</h1>

  <p align="center">
    <i>Assessing pangenome structure influence of gene co-occurrence estimation approaches in bacterial pangenomes</i>
    <br />
    <a href="https://github.com/Tomersardjoe/MSc-thesis/blob/main/MSc_thesis_Tomer_Sardjoe.pdf"><strong>Read the manuscript »</strong></a>
    <br />
    <br />
  </p>
</div>



<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#project-description">Project description</a></li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li>
      <a href="#usage">Usage</a>
      <ul>
        <li><a href="#preparation">Preparation</a></li>
        <li><a href="#gene-co-occurrence-estimation">Gene co-occurrence estimation</a></li>
        <li><a href="#pseudo-simulation">Pseudo-simulation</a></li>
        <li><a href="#aggregated-results">Aggregated results</a></li>
      </ul>
    </li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>



<!-- PROJECT DESCRIPTION -->
## Project description

This repository contains the raw data, appendices in my manuscript, and all scripts that were written or modified during my MSc thesis project at the <a href="https://www.wur.nl/en/chair-groups/bioinformatics"> Bioinformatics group</a> at <a href="https://www.wur.nl/en"> WUR</a>. 

The main objective of my thesis was to quantify the effect of pangenome structure on gene co-occurrence esetimation approaches. The main objectives of the project were:
* Create a real and two pseudo-simulated pangenomic datasets
* Compare <a href="https://github.com/fwhelan/coinfinder"> Coinfinder</a>, <a href="https://github.com/fbaumdicker/goldfinder"> Goldfinder</a> and <a href="https://github.com/alanbeavan/PanForest"> PanForest</a> gene co-occurrence estimation approaches
* Evaluate the influence of pangenome structure and phylogenetic independence (D-value)

For more information and a closer look at the findings, please read the <a href="https://github.com/Tomersardjoe/MSc-thesis/blob/main/MSc_thesis_Tomer_Sardjoe.pdf"> manuscript</a>.

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- GETTING STARTED -->
## Getting Started

To be able to run the scripts and recreate the analyses, follow the steps below.

### Prerequisites

First, clone the repository onto your machine (note that the scripts were written and tested in a Linux environment, I will assume that you will clone the repository to a Linux system as well).
* Cloning the repository
```sh
git clone https://github.com/Tomersardjoe/MSc-thesis.git
```

### Installation

The first thing you will want to do is to retrieve all the dependencies. Luckily, the project uses <a href="https://docs.conda.io/projects/conda/en/latest/index.html"> Conda</a> environments for this.

1. Create the Conda environments from the environments .yaml files. These files can be found in the conda_envs directory and the environments can be created like so:
```sh
conda env create -f coinfinder.yaml
conda env create -f goldfinder.yaml
conde env create -f panforest.yaml
conda env create -f analysis.yaml
```
**Please DO NOT change the names of the Conda environments, the scripts rely on their exact naming for environment detection**.

Next, download and install the Coinfinder, Goldfinder and PanForest software in the root of the directory (**Not for thesis assessors as the software is already installed in my Lustre environment**).
Make sure that the installed tools produced their own directory in the root of our directory and they are named "goldfinder" and "panforest"
(Coinfinder does not produce its own directory, as it is installed through miniforge3).
   
2. Download and install the three tools according to their instructions:

* <a href="https://github.com/fwhelan/coinfinder"> Coinfinder</a>
* <a href="https://github.com/fbaumdicker/goldfinder"> Goldfinder</a>
* <a href="https://github.com/alanbeavan/PanForest"> PanForest</a>

A few of Coinfinder and PanForest's scripts were modified during my thesis to make them compatible with my data.
This means that the modified scripts have to be overwritten in the tool directories (**Not for thesis assessors, my Lustre environment already has the correct files**).
You can locate these files in the "Modified scripts" directory. Below you will find the destination of each file.
  
3. Move the modified scripts to the correct directories:
```sh
mv Modified\ scripts/lineage.R miniforge3/envs/coinfinder/bin/coinfind-code/
mv Modified\ scripts/calculate_d.R panforest/
mv Modified\ scripts/calculate_d_100.R panforest/
mv Modified\ scripts/process_matrix.py panforest/
mv Modified\ scripts/rf_module.py panforest/
```
   After following the above steps, the project environment is ready for the analyses. 

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- USAGE -->
## Usage

Following all the steps described below allow the complete recreation of my work.

### Preparation

Some preparatory scripts need to be run to pre-process the data and create a selection of bacterial pangenomes for subsequent analyses.

1. Rename the "Data" directory to "real_data" and "Own scripts" to "scripts":

 ```sh
 mv Data/ real_pangenomes
 mv Own\ scripts/ scripts
 ```

The intersection in bacterial species is detected and the required files for the gene co-occurrence methods are created by running the match_all.py script.
Make sure to first activate the "analysis" Conda environment:

```sh
conda activate analysis
```

2.  Run match_all.py:
```sh
python3 scripts/prep/match_all.py --dewar real_pangenomes/Dewar_data.csv --barth real_pangenomes/Barth_Baumdicker_Weigel_data.csv --gpa real_pangenomes/GPA_matrices_red/ --tree real_pangenomes/Trees_red/ --outdir real_pangenomes
```

The script creates gpa_matches_all and tree_matches_all directories in the real_pangenomes directory. The gpa_matches_all directory contains the gene presence/absence matrices
in .csv and .tab format of the intersection of bacteria between the two input datasets. The tree_matches_all directory contains the phylogenetic trees for the pangenomes in .nwk format.
A matched_all.csv file is created that contains the intersected pangenomes between the two datasets, merged into one .csv file.
From here, the genomic fluidity and openness of the intersected pangenomes can be calculated by running the calc_alpha_fluidity.sh script.

3. Run calc_alpha_fluidity.sh:
```sh
./scripts/prep/calc_alpha_fluidity.sh --gpa real_pangenomes/gpa_matches_all/ --outdir real_pangenomes
```
This driver script calls the fluidity_calc.R script for each pangenome in the gpa_matches_all directory and calculates their genomic fluidity and openness. The output file
is called alpha_fluidity_all.csv.

Now, a selection of pangenomes can be made to continue the analyses only for the selected species (this was done during the thesis project), however, the pipeline can also be executed
on all the remaining species (114 bacterial species). For the sake of replicability I will show you how to make a selection.

4. Selecting bacterial species

Create a selection.csv file:

```sh
touch real_pangenomes/selection.csv
```

Use your favourite text editor to insert a "species_taxid" column header, followed by a species taxonomic ID on each line (these can be found in matched_all.csv). The resulting file should look like this:

```
species_taxid
632
29459
1719
485
35814
75985
287
197
520
210
1313
470
562
1396
1352
40324
```

Then run the select_species.py script to create a subset of matched_all.csv, only for the selected species, called species_categories.csv, 
a gpa_matches and tree_matches directory that contain the GPA and tree subsets for the selected species.

5. Run select_species.py
```sh
python3 scripts/prep/select_species.py --taxids real_pangenomes/selection.csv --matched real_pangenomes/matched_all.csv --gpa real_pangenomes/gpa_matches_all/ --tree real_pangenomes/tree_matches_all/ --outdir real_pangenomes
```

The pangenome structure characterisics are now know for each pangenome (fluidity and openness). The Python script categorise.py inserts the fluidity and openness column in match_all.csv
and performs k-means clustering on the fluidity and openness values to define panagenome structure classes. Optionally, it will mark the selected species in the output clustering visualisation.

6. Run categorise.py
```sh
python3 scripts/prep/categorise.py --matched real_pangenomes/matched_all.csv --alpha real_pangenomes/alpha_fluidity_all.csv [--species real_pangenomes] --outdir real_pangenomes
```

<p align="right">(<a href="#readme-top">back to top</a>)</p>

### Gene co-occurrence estimation

The data is now fully prepared for the gene co-occurrence analyses. Note that the following scripts will take a long time (over 24 hours) to run, depending on your system and number of selected pangenomes.
The steps below will again show instructions for a selection of pangenomes, however, all pangenomes can be run as well. Note that when all pangenome are run parameters in the driver
scripts are changed to speed up computation at the cost of accurracy or significance of the analyses, as well as built-in concurrency.

1. Coinfinder is run first. Activate the coinfinder Conda environment like this:
```sh
conda activate coinfinder
```

2. Then, run the run_coinfinder.sh driver script that executes Coifinder on the selected (or all) pangenomes:
```sh
./scripts/coinfinder/run_coinfinder.sh --dataset real --scope selected (or "all" if using all pangenomes)
```
This script preprocesses the phylogenetic trees by substituting zero-length branches with very small (1e-6) numbers. The driver script runs Coinfinder so that it performs no 
multiple testing and associating genes are reported regardless of significance, which is required to apply manual Benjamini-Hochberg False Discover Rate correction. Because
significance is not taken into account, the resulting network files become huge. Coinfinder currently does not have a parameter that disables the creation of these network files,
so they are simply removed.

3. Switch to the "analysis" Conda environment:
```sh
conda activate analysis
```
The output files are created in real_pangenomes/coinfinder_runs_selected. This directory contains subdirectories that are named after the species taxonomic ID of each pangenome
and at this point contain the raw Coinfinder output. Next up, BH FDR is applied by running the bh_coinfinder.sh driver script that calls the p_adj_coinfinder.R script on 
every pairs.tsv file. This script performs the BH FDR correction and applies removes any gene pairs below the 0.05 significance threshold. Then, the Python script filter_sig.py
is called, which first makes a copy of the nodes.tsv file that contains the D-value calculations for all genes and calls that file nodes_all.tsv. Then it filters the insignificant
entries from all other files that contain them (edges.tsv, nodes.tsv, and components.tsv).

4. Apply manual BH FDR correction and significance filtering with the bh_coinfinder.sh script:
```sh
./scripts/coinfinder/bh_coinfinder.sh real_pangenomes/coinfinder_runs_selected/
```

Finally, the calculation of D-value cutoffs and pagenome specific visualisations are created by running the get_d_coinfinder.sh driver script. This script calls d_distribution_coifinder.R
on all the, now BH FDR corrected significantly associated genes and gene pairs, in the coinfinder_runs_selected directory.

5. Calculate D-value cutoffs and generate visualisations with the get_d_coinfinder.sh script:
```sh
./scripts/coinfinder/get_d_coinfinder.sh --dataset real --scope selected
```
All outputs are generated for each coinfinder_runs_selected subdirectory in a "d_cutoff" directory. The output consists of the calculated D-value cutoff in a .txt file,
the gene pairs and genes that are significantly co-occurring above the D-value cutoff in a .csv file and visualisations.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

6. Goldfinder is run next, make sure to activate its Conda environment like this:
```sh
conda activate goldfinder
```

7. Then, run the run_goldfinder.sh driver script, which will execute Goldfinder on the selected pangenomes:
```sh
./scripts/goldfinder/run_goldfinder.sh --dataset real --scope selected
```
This script creates a goldfinder_runs_selected output directory, which after the script has processed the pangenomes, will contain all the output files of Goldfinder per species
taxnonomic ID subdirectory.

8. Switch to the "analysis" Conda environment:
```sh
conda activate analysis
```
Goldfinder has built-in BH FDR correction and no further filtering steps are required, so all that is left to do is to run the get_d_goldfinder.sh driver script, which will
call d_distribution_goldfinder.R that saves the D-values for each significantly associated gene pair and the genes, as well as generates visualisations per pangenome.

9. Run get_d_goldfinder.sh
```sh
./scripts/goldfinder/get_d_goldfinder.sh --dataset real --scope selected
```
All outputs are generated for each goldfinder_runs_selected subdirectory in a "d_distribution" directory. The output consists of the the gene pairs and genes that are significantly
co-occurring in a .csv file and visualisations.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

10. PanForest is run last. Don't forget to activate its Conda environment:
```sh
conda activate panforest
```
11. Run_panforest.sh driver script is run like this:
```sh
./scripts/panforest/run_panforest.sh --dataset real --mode filtered --scope selected
```

The `--mode` flag is used to determine whether D-values below zero should be filtered out or not (with "unfiltered" not filtering them out). The driver script first calls
pre_forest_d_filter.py, which filters out all genes in the input GPA matrix below a D-value of zero. Then it calls process_matrix.py, which is a script from PanForest that
removes all genes present in less than 1% amd genes present in more than 99% of genomes, as well as collapsing identical gene presence absence patterns between genes into family
groups. After processing all the pangenomes, the script calls the PanForest.py script that runs the PanForest random forest analyses on the pangenomes. The output of the analyses
are written to "panforest_runs_selected", where each subdirectory is representative of a species taxonomic ID.

12. Switch to the "analysis" Conda environment:
```sh
conda activate analysis
```

Because PanForest collapses identical genes into family groups, D-values will have to be recalculated. The calculate_d_panforest.sh script provides this functionality.

13. Run the calculate_d_panforest.sh script:
```sh
./scripts/panforest/calculate_d_panforest.sh --dataset real --mode filtered --scope selected
```

This driver script first calls prep_d_calc.R, which removes any empty columns from the collapsed matrix and creates a one-row .csv file of all the gene IDs that PanForest needs.
It then finds the associated zero-length branches replaced phylogenetic trees that were generated for Coinfinder for each pangenome. Lastly, it runs calculate_d.R, which is a
script provided by PanForest. This script calculates the D-values of all the genes and family groups. The output files are generated in an "imp_cutoff" subdirectory for every
pangenome, which at this point consists of the cleaned collapsed matrix, a nodes.tsv file, and the nodes_in.csv file with the one-row of gene IDs.

Finally, the get_d_panforest.sh driver script is run that calls imp_distribution.R for every pangenome in the panforest_runs_selected directory:

14. Run the get_d_panforest.sh script:
```sh
./scripts/panforest/get_d_panforest.sh --dataset real --mode filtered --scope selected
```
The imp_distribution.R script symmetrizes the output Gini importance score matrix, visualises what the difference in D-value cutoff would look like between an elbow detection
method and a simple Q3 cutoff (**These cutoffs are not applied to the data but were used to study their effects earlier in the project**). The only cutoffs that are applied are the fixed threshold of
0.01 in the Gini importance score matrix, as well as a 0.9 accuracy and 0.9 F1-score cutoff of the metrics that were calculated per gene in the random forest analysis.
The remaining gene pairs and genes are written to the "imp_cutoff" directory as .csv file. Other outputs in this directory are visualisations.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

### Pseudo-simulation

The steps above are all run on the real pangenomic data. Here, the steps to recreate the pseudo-simulated datasets is explained.

1. Make sure the "analysis" Conda environment is activated:
```sh
conda activate analysis
```

2. Run make_sim.sh:
```sh
./scripts/pseudo-sim/make_sim.sh --dataset perfect --scope selected
./scripts/pseudo-sim/make_sim.sh --dataset flip --scope selected
```
This driver script calls sim_pan.py that generates the pseudo-simulated datasaet. Depending on the chosen dataset (perfect = perfect co-occurrence or flip = one-flip co-occurrence),
the script takes genes between the 5th and 95th percentile of the pangenome's D-value distribution and picks 10% of the total genes, evenly spaced, for either perfect duplication
or one-flip duplication. Under perfect duplication, 10% of the total gene pool gets a "_dup" gene inserted in random positions in the GPA matrix. Under one-flip duplication,
the same process happens, but now one presence or absence of a gene in a single genome is randomly flipped to the other state. This creates a pattern of near co-occurrence. The
resulting output files are written to a simulated_pangenome_perfect or simulated_pangenome_flip in the root directory. These directories contain a gpa_matches subdirectory that
contains the GPA matrices with the inserted duplicated genes.

3. Repeat the co-occurrence analyses steps <a href=#gene-co-occurrence-estimation>above</a> (substituting --dataset real with --dataset perfect or --dataset flip).

<p align="right">(<a href="#readme-top">back to top</a>)</p>

### Aggregated results

So far, only individual pangenome results have been generated but my manuscript reports findings and interpretations of the aggregated results. This section explains how to
aggregate the results for real pangenomic data and the pseudo-simulated datasets and how to generate the visualisations used in my manuscript.

1. Make sure the "analysis" Conda environment is activated:
```sh
conda activate analysis
```

2. The real pangenomic data is summarised first. The run_real_match.sh scripts calls real_match.py:
```sh
./scripts/combined/run_real_match.sh --dataset real --mode filtered --scope selected
```
This script counts all the co-occurring gene pairs in the output file generate by each tool.
It then looks up the fluidity and openness values associated with the pangenome and appends these values to a .tsv file in the root of the real_pangenomes directory.

3. The pseudo-simulated data is summarised in a similar fashion but with more functionality. The run_dup_match.sh script is used here to call dup_match.py:
```sh
./scripts/pseudo-sim/run_dup_match.sh --dataset perfect --mode filtered --scope selected
./scripts/pseudo-sim/run_dup_match.sh --dataset flip --mode filtered --scope selected
```
This script counts the total duplicated genes that were inserted and how many of them were found. It then calculates performance metrics (precision, recall and F1-score).
It looks up the fluidity, openness and structure category of each pangenomes and writes these values to a .tsv file in the root of each of the pseudo-simulated datasets.

4. Finally, the upset_all.sh script takes care of all the statistics and visualisations used in my manuscript:
```sh
./scripts/combined/upset_all.sh --dataset real --mode filtered --scope selected
./scripts/combined/upset_all.sh --dataset perfect --mode filtered --scope selected
./scripts/combined/upset_all.sh --dataset flip --mode filtered --scope selected
```
The driver script calls aggregated_upset.R for every pangenome in the real/perfect/one-flip dataset. It loads in each tool's .csv output files, classifies overlaps between methods,
and prepares gene- and pair-level data for UpSet plots, aggregated by pangenome structure category.
For Goldfinder, the gene pairs and genes below the D-value cutoff of Coinfinder where identified and overlayed on the UpSetplots.
If a summary file is provided, which only happens if a pseudo-simulated dataset is passed to the script, the script inlcudes precision, recall and F1-score metrics in as panels
underneath the UpSet plots. It also generates a plot that shows at which D-value duplicated genes were recovered by each tool.
Statistical tests (chi-squared, Fisher's exact, proportion tests, ANOVA, and EMMeans posthoc tests) are performed to assess whether proportions of agreement differs significantly
between proportions. All output statistics and visualisations are written to the combined_results directory in the root of the project directory.



<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- LICENSE -->
## License

Distributed under the GPL-3.0 license. See <a href="https://github.com/Tomersardjoe/MSc-thesis/blob/main/LICENSE"> LICENSE</a> for more information.

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- CONTACT -->
## Contact

Tomer Sardjoe - tomer.sardjoe@wur.nl

Project Link: [https://github.com/Tomersardjoe/MSc-thesis](https://github.com/Tomersardjoe/MSc-thesis)

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- ACKNOWLEDGMENTS -->
## Acknowledgments
Franz Baumdicker for providing the data that was foundational to this work.

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- MARKDOWN LINKS & IMAGES -->
[license-shield]: https://img.shields.io/badge/License-GPL_3.0-green
[license-url]: https://github.com/Tomersardjoe/MSc-thesis/blob/main/LICENSE
[linkedin-shield]: https://img.shields.io/badge/LinkedIn-blue
[linkedin-url]: https://www.linkedin.com/in/tomersardjoe/
