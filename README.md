

# CRISPR_hands_on

<h2>Introduction</h2>

**Why Nextflow?:** I think its easier to use and simplifies the design
and execution of bioinformatics complex workflows, enabling
reproducibility and scalability across diverse computing environments.

## Summary of some of the Beneifits:

-   **Reproducibility:** Ensures consistent results across different
    environments.

-   **Scalability:** Easily scales workflows on local machines,
    clusters, or cloud platforms.

-   **Portability:** Supports multiple containerization technologies
    like Docker and Singularity.

-   **Parallelization:** Enables concurrent execution of tasks for
    faster analysis.

-   **Modularity:** Allows for flexible workflow design and reuse of
    components.

-   **Robust error handling:** Automatically manages task retries and
    recovery.

    ## Prerequisites

    Before we start ensure that you have the following:

-   conda / miniforge

-   github account

-   Linux terminal access \| wsl for windows \| virtual machine

-   Text editor e.g Vscode, nano, vim etc

<h3>Step 01: setting up our environments</h3>

-   Open your terminal and clone this repo in your home directory <br>

        cd ~

Then \`\`\` git clone '<https://github.com/jgoke1/CRISPR_hands_on.git>'

cd CRISPR_hands_on/scripts

    ***create a conda environment***

Let's name this new environment crispr_env *conda create --name
crispr_env*

activate the new environment *conda activate crispr_env*

install all the softwares to be used

*conda install -c bioconda fastplong - for quality check and trimming*

*conda install -c bioconda flye - for assembly*

*conda install -c bioconda minimap2 - for Alignment*

*conda install -c bioconda samtools - for Alignment visualization and
manipulations of alignment file format*

*conda install -c bioconda bedtools - for manipulations of file format*

*conda install -c bioconda card - for AMR genes*

      #### Download your data into our data_dir 

mkdir -p data_dir cd data_dir

    ### create a text file for our samples id

*touch sample_id.txt* *nano sample_id.txt*

ensure your sample ids are in the text file

\### create a text file for our samples id

    ### step 02: create the nexflow script
    - Lets start by navigating to the scripts dir in our cloned repo <br>

// activate the nextflow env conda activate nextflow

cd ../scripts touch my_nextflow_script.nf

    - open the script in your preferred text editor In my case i will use Vs code

code my_nextflow_script.nf \# Vscode users vim my_nextflow_script.nf \#
vim users nano my_nextflow_script.nf \# nano users

    - Basic concepts <br>
    -- Nextflow workflows are built by connecting different processes.<br>
    -- Each process can use any scripting language executable on Linux (e.g., Bash, Python).<br>
    -- Processes run independently and are isolated from one another.<br>
    -- Communication between processes occurs through asynchronous channels (FIFO queues).<br>
    -- Inputs and outputs of processes are defined as channels.<br>
    -- Workflow execution flow is determined by the input/output declarations.<br>

// 3 main parts that will make up our script - Parameters - Processes -
Workflow execution block

    ### step 03: Define parameters for our workflow in a nextflow.config 

    Create a nextflow.config file in the scripts dir

touch nextflow.config

    We shall then open the file and in put parameters of our nexflow pipeline

#!/usr/bin/env nextflow

# // ================================ // Define the workflow parameters //

// 'params' block defines user-provided inputs and options for the
workflow. // These can be specified by the user when running the
workflow or set to default values. params { // Path to a file containing
sample IDs. Each ID will be processed in the workflow. sampleID =
'/home/jovyan/Workflow-Mg-t-Nextflow-Tutorial-for-Bioinformatics/data/sample.txt'

    // Path to a file containing sample IDs. Each ID will be processed in the workflow.
        sample_id = '/data/home/kemi_adepeju/nex_crispr/sample.txt'
        
    // Path to the directory containing MiSeq sequencing reads. This will be used for analysis.
        reads = '/data/home/kemi_adepeju/nex_crispr/test_data'
        
        outdir = '/data/home/kemi_adepeju/nex_crispr/test_data/output'
        entero_ref = '/data/home/kemi_adepeju/crispr/enterobacterales_reference.mmi'
        crispr_target = '/data/home/kemi_adepeju/crispr/crispr_target.mmi'

}

# // ================================ // Executor Settings //

// The 'executor' block defines the settings for how tasks will be
executed. // This includes resource allocation and queue management for
jobs in the workflow. executor { // The 'queueSize' defines the maximum
number of tasks that can be queued at once. queueSize = 100 // Set to
allow a maximum of 100 tasks in the queue.

    // 'cpus' specifies the number of CPU cores allocated for each task.
    cpus = 4  // Each task will get 4 CPU cores.

    // 'memory' defines the amount of memory (RAM) allocated per task.
    memory = '24 GB'  // Each task is allocated 24 GB of RAM.

}

// ================================ // Additional workflow steps will
follow here // ================================

// Typically, the workflow's processes would be defined after these
parameters. // For example, processes like FASTQC, Fastp, and MultiQC
would be described, // with conditions for skipping the steps based on
the parameters defined above.

    ### step 04: Create workflow processes

    We can now go back to the nextflow script we created to define the processes for our pipeline

    #### FASTPLONG process

process FASTPLONG { // Specify the Conda environment to use for this
process conda 'fastqc'

    input:
        // Input file path for the MiSeq reads 
        path reads 
        // Input file containing a list of sample IDs
        path sample_id

        output:
        // Directory where FastQC results will be saved
        path outdir

        script:

        """
        # Create the output directory if it doesn't exist
        mkdir -p outdir/fastp_results
       
        # Loop through each sample ID from the sampleID file
        for sample in \$(cat "${sample_id}"); do
            # Run Fastplong for each sample's fastq files and save the output in the fastp_results directory
            fastplong \
                -i "${reads}"/"\${sample}".fastq \
                -o outdir/fastp_results/${sample_id}_filt_.fastq \
                -h outdir/fastp_results/${sample_id}_filt.html \
                -q 15 -u 40 -m 20 -l 500 --verbose
        done
        
           """

}

    #### MULTIQC process

process MULTIQC_01 { // Specify the Conda environment to use for this
process conda 'multiqc'

    i  input:
        // Input file path for the fastqc results folder
        path fastp_results

        output:
        // Directory where MultiQC results will be saved
        path outdir



        script:
        """
        # Create the output directory if it doesn't exist 
        mkdir -p outdir/multiq_output
        multiqc outdir/fastp_results/* --outdir outdir/multiqc_output
        """

}

    #### ASSEMBLY process

process ASSEMBLY { /

    input:
       input:
        path fastp_results
        path sample_id

        output:
        path "assembly"

        script:
        """
       mkdir -p assembly
        for sample in \$(cat "${sample_id}"); do
            flye --nano-raw fastp_results/\${sample}_filt.fastq\ \
            --out-dir ${assembly}/\${sample}.fasta\
            --genome-size 5m\
            --threads 4
        """

}

    #### ALIGNMENT process

process ALIGNMENT { // Specify the Conda environment to use for this
process conda 'fastqc'

    iinput:
        path fastp_results
        path sample_id

        output:
        path "assembly"

    script:
    """
    mkdir
    for sample in \$(cat "${sample_id}"); do
        minimap2 -ax map-ont ${entero_ref} fastp_results/\${sample}_filt_.fastq > alignment/\${sample}.bam
           samtools view -bS -F 4 alignment/\${sample}.bam > alignment/\${sample}.sam
            samtools fastq alignment/\${sample}.bam > alignment/\${sample}.fastq
        done
    """

}

    #### CRISPR_TARGET process

process CRISPR_TARGET {

    iinput:
        path assembly
        path sample_id

        output:
        path "crispr"

    script:
    """
    mkdir fast
    for sample in \$(cat "${sample_id}"); do
        minimap2 -ax map-ont ${entero_ref} fastp_results/\${sample}_filt_.fastq > alignment/\${sample}.bam
           samtools view -bS -F 4 alignment/\${sample}.bam > alignment/\${sample}.sam
            samtools fastq alignment/\${sample}.bam > alignment/\${sample}.fastq
        done
    """

}

    
    ### step 05:Create  workflow execution block
    workflow {  
        fastp_results = FASTPLONG(params.read, params.sample_id)
        multiqc_report = MULTIQC_01(fastqc_results)
        assembly = ASSEMBLY(fastp_results)
        alignment = ALIGNMENT(fastp_results, params.sample_id, params.entero_ref)
        crispr = CRISPR_TARGET(fastp_results, params.crispr_target, params.sample_id)
    }

nextflow run <your-script> -c nextflow.config \`\`\` You can also add
***-resume*** incase any of the processes fail to restart the pipeline
without re running the processes that were successful

## Resources

-   [nextflow
    training](https://training.nextflow.io/basic_training/intro/) <br>
-   [awesome nextflow
    repo](https://github.com/nextflow-io/awesome-nextflow)
