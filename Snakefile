from snakemake.utils import min_version

##### set minimum snakemake version #####
min_version("5.1.4")

##### load config file #####
configfile: "config.yaml"

##### determine cellranger output directory name #####
def get_reference(wildcards):
    output_dir=config['samples'][wildcards.sample]['reference'].split("/")[-1]
    if "cellranger" in output_dir:
        output_dir=output_dir.split("-")[2]
    full_output_dir = "/outs/filtered_gene_bc_matrices/" + output_dir
    return full_output_dir

SAMPLES = [key for key in config["samples"]]

ALL_SAMPLES = []
CONTROLS = []
for i in SAMPLES:
    if config["samples"][i]["allSamples"]:
        ALL_SAMPLES.append(i)
    if config["samples"][i]["Controls"]:
        CONTROLS.append(i)

rule all:
    input:
        input_list = expand("{data_set}/analysis_outs/seurat_{data_set}.rda", data_set = config["samples"])

rule mkfastq:
    input:
        runDir="{sample}/{sample}_raw_data",
        sampleSheet="{sample}/sample_sheet_{sample}.csv"
    output:
        directory("{sample}/{sample}_mkfastq/outs/fastq_path")
    params:
        sampleID="{sample}_mkfastq",
        sampleName="{sample}",
        runDirP="{sample}_raw_data",
        sampleSheetP="sample_sheet_{sample}.csv"
    threads: 30
    shell:
        """
        module load bcl2fastq
        cd {params.sampleName}
        rm -r {params.sampleID}
        cellranger mkfastq --run={params.runDirP} --samplesheet={params.sampleSheetP} --id={params.sampleID} --ignore-dual-index --localcores={threads}
        cd ..
        """

rule count:
    input:
        fastqDir="{sample}/{sample}_mkfastq/outs/fastq_path",
        ref=lambda wildcards: config['samples'][wildcards.sample]['reference']
    output:
        directory("{sample}/{sample}_count")
    params:
        sampleID="{sample}_count",
        sampleName="{sample}",
        mem=config['mem'],
        fastqDirP="{sample}_mkfastq/outs/fastq_path"
    threads: 30
    shell:
        """
        cd {params.sampleName}
        cellranger count --id={params.sampleID} --fastqs={params.fastqDirP} --sample={params.sampleName} --transcriptome={input.ref} --localcores={threads} --localmem={params.mem}
        cd ..
        """

rule seurat_object:
    input:
        "{sample}/{sample}_count"
    output:
        "{sample}/analysis_outs/seurat_{sample}_empty.rda"
    params:
        proj=config['project'],
        ref_dir=get_reference
    conda:
        "envs/create_seurat.yaml"
    script:
        "scripts/create_seurat.R"

rule final_analysis:
    input:
        "{sample}/analysis_outs/seurat_{sample}_empty.rda"
    output:
        pdf="{sample}/analysis_outs/{sample}_qc_images.pdf",
        seurat="{sample}/analysis_outs/seurat_{sample}.rda"
    params:
        image_pdf="{sample}/analysis_outs/{sample}_images.pdf",
        data_dir=config['data_dir'],
        pseudotime="{sample}/analysis_outs/{sample}_pseudotime.rda",
        script = "{sample}/scripts/analysis_driver.R"
    conda:
        "envs/analysis_driver.yaml"
    script:
        "{params.script}"


rule combine_samples:
    input:
        data_list = expand("{data_set}/analysis_outs/seurat_{data_set}_empty.rda", data_set = ALL_SAMPLES)
    output:
        "{sample}/analysis_outs/{sample}_combined_empty.rda"
    params:
        qual_pdf = "{sample}/analysis_outs/{sample}_combined_quality.pdf",
        data_names = SAMPLES,
        script = "{sample}/scripts/create_combined_seurat.R"
    conda:
        "envs/create_seurat.yaml"
    script:
        "{params.script}"

rule combine_controls:
    input:
        data_list = expand("{data_set}/analysis_outs/seurat_{data_set}_empty.rda", data_set = CONTROLS)
    output:
        "{sample}/analysis_outs/{sample}_merged_empty.rda"
    params:
        qual_pdf = "{sample}/analysis_outs/{sample}_merged_quality.pdf",
        data_names = CONTROLS,
        script = "{sample}/scripts/create_combined_seurat.R"
    conda:
        "envs/create_seurat.yaml"
    script:
        "{params.script}"

rule combined_analysis:
    input:
        seurat_object = "{sample}/analysis_outs/{sample}_{merge_type}_empty.rda",
        mapping_object = "{mapping}/analysis_outs/seurat_{mapping}.rda".format(mapping = config['mapping_sample'])
    params:
        image_pdf="{sample}/analysis_outs/{sample}_images.pdf",
        data_dir = config['data_dir'],
        script = "{sample}/scripts/analysis_driver.R"
    output:
        seurat_object = "{sample}/analysis_outs/seurat_{sample}_{merge_type}.rda",
        slingshot_object = "{sample}/analysis_outs/{sample}_{merge_type}_slingshot.rda"
    conda:
        "envs/analysis_driver.yaml"
    script:
        "{params.script}"

rule subset_combined:
    input:
        "{sample}/analysis_outs/seurat_{sample}_combined.rda"
    params:
        subset_name = "{subset_name}",
        pdf_file = "{sample}/analysis_outs/{sample}_{subset_name}_combined_images.pdf"
    output:
        "{sample}/analysis_outs/seurat_{sample}_{subset_name}_combined.rda"
    conda:
        "envs/create_seurat.yaml"
    script:
        "scripts/subset_analysis.R"

rule make_figures:
    input:
        config['fig_files']
    params:
        data_dir = config['data_dir'],
        save_dir = "figure_output"
    output:
        "figure_output/complete_figs.txt"
    conda:
        "envs/figures.yaml"
    script:
        "scripts/figures.R"
