# Association of diet and antimicrobial resistance
 Scripts and docker images used in paper: Association of diet and antimicrobial resistance in healthy US adults. Paper has been published in XXXXXX and can be found here: XXXXXX

### **Data Availibility**
- Metagenomes are deposited in NCBI Sequence Read Archive (SRA) under the [study accession SRP354271](https://dataview.ncbi.nlm.nih.gov/object/PRJNA795985). Requests for non-metagenomic data from the USDA ARS WHNRC Nutritional Phenotyping Study used in this analysis should be made via an email to the senior WHNRC author on the publication of interest. Requests will be reviewed quarterly by a committee consisting of the study investigators.

________________________________________
### **Containers for reproducibility**
_________________________________________________
- A docker container for analysis in R and ML analysis are provided. Docker must be installed to run. These images were built and run on Docker v4.1.1. To run, inside the directory of the appropriate folder (R or python), run: 
    - ```docker build -t [tag_name] .```
    - Example: ```docker build -t amr_r_env:1.0 .```
- In order to run the R container: 
    1. first clone the repo to your local machine
    2. ```docker run --rm -it -p 8787:8787 -e PASSWORD=yourpasswordhere -v path/to/ARGs_and_Diet/Rcode:/home/Rcode/ -v path/to/ARGs_and_Diet/data:/home/data/ amr_r_env``` (or whatever you named the container)
    3. navigate to http://localhost:8787/ in a browser window
    4. log into the Rstudio local server
        - username: rstudio
        - password: yourpasswordhere (if you didnt set one in docker run command)
    5. change to the scripts working directory inside R.
        - setwd("/home/Rcode")
    6. navigate the filesystem to the working directory.
        ![plot showing changing working directory in the file pane](https://github.com/aoliver44/ARGs_and_Diet/blob/main/utilities/readme_picture.png)


- In order to run the ML container: 
    1. ```docker run --rm -it -p 8888:8888 -v path/to/ARG_scripts_and_data/:/home/ amr_py_env /bin/bash -c "cd /home && jupyter notebook --ip='0.0.0.0' --port=8888 --no-browser --allow-root"``` (or whatever you named the ML container)
    2. navigate to http://localhost:8888/ in a browser window
    3. copy a token from the terminal and paste it in the browser window where it asks for one.

### **Code**
- all code used in this analysis is provided
- R analyses were run locally on a intel-based macbook pro
- ML analyses were run on Ceres, a supercomputer for use by USDA researchers
- Sequence preprocessing was run on Spitfire, a slurm-based HPC cluster managed by the UC Davis Genome Center
    - code for sequencing preprocessing can be found [here]()

### **Figures**
- to generate most of the main/supp figures and tables in the manuscript, run generate_figures.R script with the necessary scripts for sourcing and necessary input data