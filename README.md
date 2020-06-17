# Personality and shell morph in *Cepaea nemoralis* snails

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3899042.svg)](https://doi.org/10.5281/zenodo.3899042)

This repo contains all data and code needed to re-do the analyses and figures in our manuscript

"Boldness and exploration are linked to shell morph but not environmental context in the snail *Cepaea nemoralis*"  
(by Maxime Dahirel, Valentin Gaudu, Armelle Ansart)

(link to bioRxiv preprint [here](https://doi.org/10.1101/866947))

data in `.csv` format are in the `data` folder, R script (including detailed information about the analysis) in the `R` folder.

This folder is a RStudio project folder, and the script uses the [`here` package](https://here.r-lib.org/) (see also [here](https://github.com/jennybc/here_here)). 
This means all file paths are relative, and the analysis *should* work on your computer no questions asked, whether you use the R project or not, no code line to change as long as you download the entire repo (you just need to install all the needed packages first, of course).

If you run the script for the first time, models and some other time-consuming outputs will be saved in the `R_output` folder so you don't have to re-run them everytime. 

If you don't want/don't have the time to run them, copy and unzip into the `R_output` folder the content of the zipped folder named `cepaea-personality-2017-rdata-output.zip` that is bundled with each release [here](https://github.com/mdahirel/cepaea-personality-2017/releases) or on [Zenodo](https://doi.org/10.5281/zenodo.3899042). It is an exact copy of the model files I saved in my own `R_output` folder (which were added to `gitignore` because of their size).
