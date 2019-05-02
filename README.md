This snakemake pipeline contains all the necessary scripts to recreate the analysis for (eventual paper).

Writen by Kristen Wells April 11, 2019. Last modified April 23, 2019

To use:

1. Download and install miniconda3: For Linux

'''
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh bash Miniconda3-latest-Linux-x86_64.sh
'''

2. Copy this entire directory to where you wish to run this pipeline. Inside this directory:
  - make a directory for each file

'''
mkdir {sample}
cd {sample}
mkdir {sample}_raw_data
'''

  - move sample tarball and meta data tarball into the raw_data folder and untar.
  - Place the sample sheet named in the location {sample}/sample_sheet_{sample}.csv This must be repeated for all samples in the experiment.
3. Install Snakemake:

'''
conda install snakemake -c bioconda -c conda-forge
'''
4. Download and extract the mouse reference genome

'''
wget http://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-mm10-1.2.0.tar.gz
tar -xzvf refdata-cellranger-mm10-1.2.0.tar.gz
'''

5. Update the config file (config.yaml)
  - mapping_sample: The sample used to map clusters for the combined analysis
  - samples: a list of all samples in the experiment. Include path to the reference and gtf file for each sample. Also include if the sample should be included in the allSamples run and the Controls run
  - project: The name of the project
  - data_dir: The path to the directory containing any extra data files
  - mem: The amount of memory to give the mapping rules
  - fig_files: A list of all files needed for the final analysis. These all must be made at some point during the snakemake pipeline
6. If you are running on a cluster, update the cluster.json for with your specific specs

To make a seurat object for all samples from the 10x run

'''
snakemake --cores 30
'''

or submit with

'''
sbatch analyze_individ.sh
'''

To make figures, it's best to submit this to the cluster using

'''
sbatch submit_figures.sh
'''

After running the submit_figures.sh script, there will be new directories for your combined analysis, and analysis_out directories in each sample directory containing images from the analysis. There will also be a "figure_output" director containing all figures.

This pipeline has been developed and tested with Snakemake 5.2.2 or higher. Check your version by typing: snakemake --version

If you are running a version of Snakemake prior to 5.2.2, update to the latest version: conda update snakemake -c bioconda -c conda-forge

The envrionments for each rule are included under the rule header

Known problems: All samples will run in parallel. Sometimes one sample will crash during the count and mkfastq rules. Just rerunning the snakefile once the other samples have completed will likely fix this.

Additionally, sometimes loading the conda environment fails. Again, this is seemingly random, and running the submit_figures.sh again once the original run has ended will likely fix the problem.
