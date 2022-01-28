# Association of diet and antimicrobial resistance
 Scripts and docker images used in paper: Association of diet and antimicrobial resistance in healthy US adults. Paper has been published in XXXXXX and can be found here: XXXXXX

### **Data Availibility**
- Metagenomes are deposited in NCBI Sequence Read Archive (SRA) under the [study accession SRP354271](https://dataview.ncbi.nlm.nih.gov/object/PRJNA795985). Requests for non-metagenomic data from the USDA ARS WHNRC Nutritional Phenotyping Study used in this analysis should be made via an email to the senior WHNRC author on the publication of interest. Requests will be reviewed quarterly by a committee consisting of the study investigators.

### **Containers for reproducibility**
- A docker container for analysis in R and ML analysis are provided. Docker must be installed to run. These images were built and run on Docker v4.1.1. To run, inside the directory of the appropriate folder (R or python), run: ```docker build -t [tag_name] .```
- In order to run the R container, run: ```docker run --rm -it -p 8787:8787 -e PASSWORD=yourpasswordhere -v path/to/ARG_scripts_and_data:/home amr_r_env``` (or whatever you named the container)
- In order to run the ML container: ```docker run --rm -it -p 8888:8888 -v path/to/ARG_scripts_and_data/:/home/ amr_py_env bash``` (or whatever you named the ML container)
    - then run: ```jupyter notebook --ip 0.0.0.0 --p 8888 --no-browser --allow-root```
- both of these ```docker run``` commands will spin up a local host server. If  you are using chrome, navigate to the appropriate site (i.e., localhost:8787)

### **Code**
- all code used in this analysis is provided
- R analyses were run locally on a intel-based macbook pro
- ML analyses were run on Ceres, a supercomputer for use by USDA researchers
- Sequence preprocessesing was run on Spitfire, a slurm-based HPC cluster managed by the UC Davis Genome Center

### **Figures**
- to generate most of the main/supp figures and tables in the manuscript, run generate_figures.R script with the necessary scripts for sourcing and necessary input data