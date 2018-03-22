configfile: "Snake.config.json"
import os.path


wildcard_constraints:
    binsize  = "\d+",
    seed     = "\d+",
    windows  = "\d+_[a-zA-Z]+",
    chrom    = "[chr0-9XY]+",
    segments = "[a-zA-Z]+",
    method   = "[a-zA-Z0-9_]+"



SIMUL_WINDOW = [50000]
SIMUL_SEEDS  = [5]
METHODS      = ["maryam", "simple"]
CHROMOSOMES  = config['chromosomes']


rule simul:
    input:
        expand("plots/simulation{seed}-{binsize}/{binsize}_fixed.pdf",
                seed   = SIMUL_SEEDS,
                binsize = SIMUL_WINDOW),
        expand("evaluation/simulation{seed}_{binsize}.{segments}.{method}.pdf",
                seed  = SIMUL_SEEDS,
                binsize = SIMUL_WINDOW,
                segments = ["few","medium"],
                method = METHODS),
        expand("sv_calls/simulation{seed}-{binsize}/{binsize}_fixed.{segments}/{method}.{chrom}.pdf",
                seed   = SIMUL_SEEDS,
                binsize = SIMUL_WINDOW,
                segments = ["few","medium"],
                chrom = CHROMOSOMES,
                method = METHODS)


################################################################################
# Simulation of count data                                                     #
################################################################################

rule evaluate_simulation:
    input:
        prob = "sv_calls/simulation{seed}-{binsize}/{binsize}_fixed.{segments}/{method}.txt",
        simul = "simulation/variants/genome{seed}-{binsize}.txt"
    output:
        "evaluation/simulation{seed}_{binsize}.{segments}.{method}.pdf"
    script:
        "utils/evaluate_simulation.R"

rule simulate_genome:
    output:
        tsv="simulation/genome/genome{seed}.tsv"
    log:
        "log/simulate_genome/genome{seed}.tsv"
    params:
        svcount     =     200,
        minsize     =  100000,
        maxsize     = 5000000,
        mindistance = 1000000,
    shell:
        "utils/simulate_SVs.R {wildcards.seed} {params.svcount} {params.minsize} {params.maxsize} {params.mindistance} {output.tsv} > {log} 2>&1"

rule add_vafs_to_simulated_genome:
    input:
        tsv="simulation/genome/genome{seed}.tsv"
    output:
        tsv="simulation/genome-with-vafs/genome{seed}.tsv"
    params:
        min_vaf = config["simulation_min_vaf"],
        max_vaf = config["simulation_max_vaf"],
    shell:
        """
        awk -v min_vaf={params.min_vaf} -v max_vaf={params.max_vaf} -v seed={wildcards.seed} \
        'BEGIN {{srand(seed); OFS="\\t"}} {{vaf=min_vaf+rand()*(max_vaf-min_vaf); print $0, vaf}}' {input.tsv} > {output.tsv}
        """

def min_coverage(wildcards):
    return round(float(config["simulation_min_reads_per_library"]) * int(wildcards.binsize) / float(config["genome_size"]))

def max_coverage(wildcards):
    return round(float(config["simulation_max_reads_per_library"]) * int(wildcards.binsize) / float(config["genome_size"]))

def neg_binom_p(wildcards):
    return float(config["simulation_neg_binom_p"][wildcards.binsize])

rule simulate_counts:
    input:
        config="simulation/genome-with-vafs/genome{seed}.tsv",
    output:
        counts="simulation/counts/genome{seed}-{binsize}.txt.gz",
        segments="simulation/segments/genome{seed}-{binsize}.txt",
        phases="simulation/phases/genome{seed}-{binsize}.txt",
        info="simulation/info/genome{seed}-{binsize}.txt",
        sce="simulation/sce/genome{seed}-{binsize}.txt",
        variants="simulation/variants/genome{seed}-{binsize}.txt",
    params:
        mc_command   = config["mosaicatcher"],
        neg_binom_p  = neg_binom_p,
        min_coverage = min_coverage,
        max_coverage = max_coverage,
        cell_count   = config["simulation_cell_count"],
        alpha        = config["simulation_alpha"],
    log:
        "log/simulate_counts/genome{seed}-{binsize}.txt"
    shell:
        """
            {params.mc_command} simulate \
            -w {wildcards.binsize} \
            --seed {wildcards.seed} \
            -n {params.cell_count} \
            -p {params.neg_binom_p} \
            -c {params.min_coverage} \
            -C {params.max_coverage} \
            -a {params.alpha} \
            -V {output.variants} \
            -i {output.info} \
            -o {output.counts} \
            -U {output.segments} \
            -P {output.phases} \
            -S {output.sce} \
            --sample-name simulation{wildcards.seed}-{wildcards.binsize} \
            {input.config} > {log} 2>&1
        """

rule link_to_simulated_counts:
    input:
        counts="simulation/counts/genome{seed}-{binsize}.txt.gz",
        info="simulation/info/genome{seed}-{binsize}.txt",
    output:
        counts = "counts/simulation{seed}-{binsize}/{binsize}_fixed.txt.gz",
        info   = "counts/simulation{seed}-{binsize}/{binsize}_fixed.info"
    run:
        d = os.path.dirname(output.counts)
        count_file = os.path.basename(output.counts)
        info_file = os.path.basename(output.info)
        shell("cd {d} && ln -s ../../{input.counts} {count_file} && ln -s ../../{input.info} {info_file} && cd ../..")


rule link_to_simulated_strand_states:
    input:
        sce="simulation/sce/genome{seed}-{binsize}.txt",
    output:
        states="strand_states/simulation{seed}-{binsize}/final.txt",
    run:
        d = os.path.dirname(output.states)
        f = os.path.basename(output.states)
        shell("cd {d} && ln -s ../../{input.sce} {f} && cd ../..")




################################################################################
# Plots                                                                        #
################################################################################

rule plot_mosaic_counts:
    input:
        counts = "counts/{sample}/{file_name}.txt.gz",
        info   = "counts/{sample}/{file_name}.info"
    output:
        "plots/{sample}/{file_name}.pdf"
    log:
        "log/plot_mosaic_counts/{sample}.{file_name}.txt"
    params:
        plot_command = "Rscript " + config["plot_script"]
    shell:
        """
        {params.plot_command} {input.counts} {input.info} {output} > {log} 2>&1
        """

rule plot_SV_calls:
    input:
        counts = "counts/{sample}/{windows}.txt.gz",
        probs  = "sv_calls/{sample}/{windows}.{segments}/{method}.txt"
    output:
        expand("sv_calls/{{sample}}/{{windows}}.{{segments}}/{{method}}.{chrom}.pdf", chrom = config['chromosomes'])
    log:
        "log/{sample}/plot_SV_call.{windows}.{segments}.txt"
    params:
        plot_command = "Rscript " + config["sv_plot_script"]
    shell:
        """
        {params.plot_command} \
            {input.counts} {input.probs} \
            sv_calls/{wildcards.sample}/{wildcards.windows}.{wildcards.segments} \
            > {log} 2>&1
        """


################################################################################
# Segmentation                                                                 #
################################################################################

rule segmentation:
    input:
        "counts/{sample}/{file_name}.txt.gz"
    output:
        "segmentation/{sample}/{file_name}.txt"
    log:
        "log/{sample}/segmentation.{file_name}.txt"
    params:
        mc_command = config["mosaicatcher"]
    shell:
        """
        {params.mc_command} segment \
        -o {output} \
        {input} > {log} 2>&1
        """

rule plot_segmentation:
    input:
        counts = "counts/{sample}/{file_name}.txt.gz",
        segments = "segmentation/{sample}/{file_name}.txt"
    output:
        "segmentation/{sample}/{file_name}/{chrom}.pdf"
    log:
        "log/{sample}/plot_segmentation.{file_name}.{chrom}.txt"
    params:
        command = config["plot_segments"]
    shell:
        """
        Rscript {params.command} {input.counts} {input.segments} {wildcards.chrom} {output} > {log} 2>&1
        """



# Pick a few segmentations and prepare the input files for SV classification
rule prepare_segments:
    input:
        "segmentation/{sample}/{windows}.txt"
    output:
        "segmentation2/{sample}/{windows}.{bpdens}.txt"
    log:
        "log/{sample}/prepare_segments.{windows}.{bpdens}.txt"
    params:
        quantile = lambda wc: config["bp_density"][wc.bpdens]
    script:
        "utils/helper.prepare_segments.R"


################################################################################
# SV classification                                                            #
################################################################################

rule install_MaRyam:
    output:
        "utils/R-packages2/MaRyam/R/MaRyam"
    log:
        "log/install_MaRyam.log"
    shell:
        """
        TAR=$(which tar) Rscript utils/install_maryam.R > {log} 2>&1
        """

rule run_sv_classification:
    input:
        maryam = "utils/R-packages2/MaRyam/R/MaRyam",
        counts = "counts/{sample}/{windows}.txt.gz",
        info   = "counts/{sample}/{windows}.info",
        states = "strand_states/{sample}/final.txt",
        bp     = "segmentation2/{sample}/{windows}.{segments}.txt"
    output:
        outdir = "sv_probabilities/{sample}/{windows}.{segments}/",
        out1   = "sv_probabilities/{sample}/{windows}.{segments}/allSegCellProbs.table",
        out2   = "sv_probabilities/{sample}/{windows}.{segments}/allSegCellGTprobs.table",
        bamNames = "sv_probabilities/{sample}/{windows}.{segments}/bamNames.txt"
    log:
        "log/{sample}/run_sv_classification.{windows}.{segments}.txt"
    params:
        windowsize    = lambda wc: wc.windows.split("_")[0]
    shell:
        """
        set -x
        # set haplotypeInfo if phasing info is available
        Rscript utils/MaRyam_pipeline.R \
                binRCfile={input.counts} \
                BRfile={input.bp} \
                infoFile={input.info} \
                stateFile={input.states} \
                outputDir={output.outdir} \
                bin.size={params.windowsize} \
                K=22 \
                maximumCN=4 \
                utils/R-packages2/ > {log} 2>&1
        """

rule convert_SVprob_output:
    input:
        probs    = "sv_probabilities/{sample}/{windows}.{segments}/allSegCellProbs.table",
        info     = "counts/{sample}/{windows}.info",
        bamNames = "sv_probabilities/{sample}/{windows}.{segments}/bamNames.txt"
    output:
        "sv_calls/{sample}/{windows}.{segments}/maryam.txt"
    params:
        sample_name = lambda wc: wc.sample
    log:
        "log/{sample}/convert_SVprob_output.{windows}.{segments}.txt"
    script:
        "utils/helper.convert_svprob_output.R"

### New SV classification

rule run_sv_classification_preparation:
    input:
        counts = "counts/{sample}/{windows}.txt.gz",
        info   = "counts/{sample}/{windows}.info",
        states = "strand_states/{sample}/final.txt",
        bp     = "segmentation2/{sample}/{windows}.{bpdens}.txt"
    output:
        "sv_probabilities/{sample}/{windows}.{bpdens}/raw_probabilities.Rdata"
    script:
        "utils/sv_classifier.R"


rule sv_classifier_simple:
    input:
        "sv_probabilities/{sample}/{windows}.{bpdens}/raw_probabilities.Rdata"
    output:
        prob = "sv_calls/{sample}/{windows}.{bpdens}/simple.txt"
    shell:
        "echo 'work in progress'"

