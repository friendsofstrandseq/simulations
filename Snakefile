configfile: "Snake.config.json"
import os.path


wildcard_constraints:
    binsize   = "\d+",
    seed      = "\d+",
    windows   = "\d+_[a-zA-Z]+",
    chrom     = "[chr0-9XY]+",
    segments  = "[a-zA-Z0-9]+",
    method    = "[a-zA-Z0-9_-]+",
    minsvsize = "\d+",
    maxsvsize = "\d+",
    minvaf    = "\d+",
    maxvaf    = "\d+",
    mincells  = "\d+",
    llr       = "\d"


localrules:
    simul,
    simulate_genome,
    add_vafs_to_simulated_genome,
    link_to_simulated_counts,
    link_to_simulated_strand_states,
    new_simulate_genome,
    new_link_to_simulated_counts,
    new_link_to_simulated_strand_states,
    prepare_segments,
    install_mosaiClassifier,
    convert_SVprob_output,
    sv_classifier_filter


### Global settings
NUM_SVS         = config["num_SVs_per_genome"]
SIMUL_SEEDS     = [1]
SIMUL_WINDOW    = [50000]
NUM_CELLS       = config["num_cells"]
METHODS         = ["simple_llr2",
                   "simple_llr4",
                   "mc_simple",
                   "mc_biallelic"]
CHROMOSOMES     = config['chromosomes']
SEGMENTS        = ["fraction10",
                   "fraction20",
                   "fraction30",
                   "fraction50",
                   "fraction70",
                   "fraction100"]
SIZE_RANGES     = ["200000-500000",
                   "500000-1000000",
                   "1000000-2000000",
                   "2000000-5000000"]
VAF_RANGES      = ["2-5", "5-10", "10-20", "20-40", "40-80", "95-100"]
SV_CLASSES      = ["het_del", "het_dup", "het_inv"]


### Main rule
rule simul:
    input:
        # New SV evaluation curves
        # (based on separate simulations for SV sizes and VAFs)
        #
        # Do not Plot SV calls of some of the new simulations
        # expand("sv_plots/seed{seed}_size{sizerange}_vaf{vafrange}_{svclass}-{binsize}/{binsize}_fixed.{segments}/{method}.{chrom}.pdf",
        #                seed      = SIMUL_SEEDS[1:2],
        #                sizerange = SIZE_RANGES,
        #                vafrange  = VAF_RANGES,
        #                segments  = SEGMENTS,
        #                binsize   = SIMUL_WINDOW,
        #                method    = METHODS,
        #                chrom     = CHROMOSOMES,
        #                svclass   = SV_CLASSES),
        expand("results/meta-{binsize}.pdf",
                binsize = SIMUL_WINDOW)


def min_coverage(wildcards):
    return round(float(config["simulation_min_reads_per_library"]) * int(wildcards.binsize) / float(config["genome_size"]))

def max_coverage(wildcards):
    return round(float(config["simulation_max_reads_per_library"]) * int(wildcards.binsize) / float(config["genome_size"]))

def neg_binom_p(wildcards):
    return float(config["simulation_neg_binom_p"][wildcards.binsize])




################################################################################
# New Simulation of count data                                                 #
################################################################################

rule new_simulate_genome:
    output:
        "simulation_new/seed{seed}_size{minsvsize}-{maxsvsize}_vaf{minvaf}-{maxvaf}_{svclass}/genome.tsv"
    log:
        "log/simulation_new/seed{seed}_size{minsvsize}-{maxsvsize}_vaf{minvaf}-{maxvaf}_{svclass}.txt"
    params:
        svcount     = NUM_SVS,
        mindistance = 1000000
    shell:
        """
        tmp=$(mktemp)
        Rscript utils/simulate_SVs.R {wildcards.seed} {params.svcount} \
            {wildcards.minsvsize} {wildcards.maxsvsize} {params.mindistance} \
            {wildcards.svclass} $tmp \
            > {log} 2>&1
        awk -v minvaf={wildcards.minvaf} -v maxvaf={wildcards.maxvaf} \
            -v seed={wildcards.seed} \
            'BEGIN {{srand(seed); OFS="\\t"}}
            {{ vaf = int(minvaf + 1.0 * rand() * (maxvaf - minvaf))/100.0; print $0, vaf}}' $tmp \
        > {output} 2> {log}
        """

rule new_simulate_counts:
    input:
        config = "simulation_new/seed{seed}_size{minsvsize}-{maxsvsize}_vaf{minvaf}-{maxvaf}_{svclass}/genome.tsv"
    output:
        counts   = "simulation_new/seed{seed}_size{minsvsize}-{maxsvsize}_vaf{minvaf}-{maxvaf}_{svclass}/counts-{binsize}.txt.gz",
        segments = "simulation_new/seed{seed}_size{minsvsize}-{maxsvsize}_vaf{minvaf}-{maxvaf}_{svclass}/segments-{binsize}.txt",
        phases   = "simulation_new/seed{seed}_size{minsvsize}-{maxsvsize}_vaf{minvaf}-{maxvaf}_{svclass}/phases-{binsize}.txt",
        info     = "simulation_new/seed{seed}_size{minsvsize}-{maxsvsize}_vaf{minvaf}-{maxvaf}_{svclass}/info-{binsize}.txt",
        sce      = "simulation_new/seed{seed}_size{minsvsize}-{maxsvsize}_vaf{minvaf}-{maxvaf}_{svclass}/sce-{binsize}.txt",
        variants = "simulation_new/seed{seed}_size{minsvsize}-{maxsvsize}_vaf{minvaf}-{maxvaf}_{svclass}/variants-{binsize}.txt",
    params:
        mc_command   = config["mosaicatcher"],
        neg_binom_p  = neg_binom_p,
        min_coverage = min_coverage,
        max_coverage = max_coverage,
        cell_count   = NUM_CELLS,
        alpha        = config["simulation_alpha"]
    log:
        "log/simulate_counts/seed{seed}_size{minsvsize}-{maxsvsize}_vaf{minvaf}-{maxvaf}_{svclass}.txt"
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
            --sample-name seed{wildcards.seed}_size{wildcards.minsvsize}-{wildcards.maxsvsize}_vaf{wildcards.minvaf}-{wildcards.maxvaf}-{wildcards.binsize}_{wildcards.svclass} \
            {input.config} > {log} 2>&1
        """

rule new_link_to_simulated_counts:
    input:
        counts = "simulation_new/seed{seed}_size{minsvsize}-{maxsvsize}_vaf{minvaf}-{maxvaf}_{svclass}/counts-{binsize}.txt.gz",
        info   = "simulation_new/seed{seed}_size{minsvsize}-{maxsvsize}_vaf{minvaf}-{maxvaf}_{svclass}/info-{binsize}.txt",
    output:
        counts = "counts/seed{seed}_size{minsvsize}-{maxsvsize}_vaf{minvaf}-{maxvaf}_{svclass}-{binsize}/{binsize}_fixed.txt.gz",
        info   = "counts/seed{seed}_size{minsvsize}-{maxsvsize}_vaf{minvaf}-{maxvaf}_{svclass}-{binsize}/{binsize}_fixed.info"
    run:
        d = os.path.dirname(output.counts)
        count_file = os.path.basename(output.counts)
        info_file  = os.path.basename(output.info)
        shell("cd {d} && ln -s ../../{input.counts} {count_file} && ln -s ../../{input.info} {info_file} && cd ../..")


rule new_link_to_simulated_strand_states:
    input:
        sce    = "simulation_new/seed{seed}_size{minsvsize}-{maxsvsize}_vaf{minvaf}-{maxvaf}_{svclass}/sce-{binsize}.txt",
    output:
        states = "strand_states/seed{seed}_size{minsvsize}-{maxsvsize}_vaf{minvaf}-{maxvaf}_{svclass}-{binsize}/final.txt"
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

rule plot_SV_calls_new:
    input:
        counts = "counts/seed{seed}_size{minsvsize}-{maxsvsize}_vaf{minvaf}-{maxvaf}-{binsize}_{svclass}/{binsize}_fixed.txt.gz",
        segs   = "segmentation2/seed{seed}_size{minsvsize}-{maxsvsize}_vaf{minvaf}-{maxvaf}-{binsize}_{svclass}/{binsize}_fixed.{segments}.txt",
        svs    = "sv_calls/seed{seed}_size{minsvsize}-{maxsvsize}_vaf{minvaf}-{maxvaf}-{binsize}_{svclass}/{binsize}_fixed.{segments}/{method}.txt",
        simul  = "simulation_new/seed{seed}_size{minsvsize}-{maxsvsize}_vaf{minvaf}-{maxvaf}_{svclass}/variants-{binsize}.txt"
    output:
        expand("sv_plots/seed{{seed}}_size{{minsvsize}}-{{maxsvsize}}_vaf{{minvaf}}-{{maxvaf}}_{{svclass}}-{{binsize}}/{{binsize}}_fixed.{{segments}}/{{method}}.{chrom}.pdf", chrom = config['chromosomes'])
    script:
        "utils/plot_sv_calls.R"



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
        -m 0.2 -M 50000000 -o {output} \
        {input} > {log} 2>&1
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

### New SV classification (Sascha)

rule sv_classifier_preparation:
    input:
        counts = "counts/{sample}/{windows}.txt.gz",
        info   = "counts/{sample}/{windows}.info",
        states = "strand_states/{sample}/final.txt",
        bp     = "segmentation2/{sample}/{windows}.{bpdens}.txt"
    output:
        "sv_probabilities/{sample}/{windows}.{bpdens}/raw_probabilities.Rdata"
    script:
        "utils/sv_classifier.R"


rule sv_classifier_mostsimple:
    input:
        prob = "sv_probabilities/{sample}/{windows}.{bpdens}/raw_probabilities.Rdata"
    output:
        "sv_calls/{sample}/{windows}.{bpdens}/simple_llr{llr}.txt"
    params:
        llr_cutoff = lambda wc: wc.llr
    script:
        "utils/sv_classifier_mostsimple.R"


rule sv_classifier_biallelic:
    input:
        "sv_probabilities/{sample}/{windows}.{bpdens}/raw_probabilities.Rdata"
    output:
        prob = "sv_calls/{sample}/{windows}.{bpdens}/biallelic_min{mincells}cells.txt"
    params:
        mincells = lambda wc: wc.mincells
    script:
        "utils/sv_classifier_biallelic.R"



### Newest classifier: "MosaiClassifier"

rule install_mosaiClassifier:
    output:
        "utils/mosaiClassifier.snakemake.R"
    shell:
        """
        git clone https://github.com/friendsofstrandseq/pipeline.git mosaiClassifier
        cd mosaiClassifier
        cd ../utils
        ln -s ../mosaiClassifier/utils/mosaiClassifier
        ln -s ../mosaiClassifier/utils/mosaiClassifier.snakemake.R
        ln -s ../mosaiClassifier/utils/mosaiClassifier_call.snakemake.R
        ln -s ../mosaiClassifier/utils/mosaiClassifier_call_biallelic.snakemake.R
        """

rule mosaiClassifier_calc_probs:
    input:
        installer = "utils/mosaiClassifier.snakemake.R",
        counts = "counts/{sample}/{windows}.txt.gz",
        info   = "counts/{sample}/{windows}.info",
        states = "strand_states/{sample}/final.txt",
        bp     = "segmentation2/{sample}/{windows}.{bpdens}.txt"
    output:
        output = "sv_probabilities/{sample}/{windows}.{bpdens}/probabilities.Rdata"
    log:
        "log/{sample}/mosaiClassifier_calc_probs.{windows}.{bpdens}.txt"
    script:
        "utils/mosaiClassifier.snakemake.R"

rule mosaiClassifier_make_call:
    input:
        probs = "sv_probabilities/{sample}/{windows}.{bpdens}/probabilities.Rdata"
    output:
        "sv_calls/{sample}/{windows}.{bpdens}/mc_simple.txt"
    script:
        "utils/mosaiClassifier_call.snakemake.R"

rule mosaiClassifier_make_call_biallelic:
    input:
        probs = "sv_probabilities/{sample}/{windows}.{bpdens}/probabilities.Rdata"
    output:
        "sv_calls/{sample}/{windows}.{bpdens}/mc_biallelic.txt"
    script:
        "utils/mosaiClassifier_call_biallelic.snakemake.R"





################################################################################
# Evaluation of MosaiCatcher                                                   #
################################################################################


rule evaluation_newer_version:
    input:
        truth = expand("simulation_new/seed{seed}_size{sizerange}_vaf{vafrange}_{svclass}/variants-{{binsize}}.txt",
                       seed = SIMUL_SEEDS,
                       sizerange = SIZE_RANGES,
                       vafrange  = VAF_RANGES,
                       svclass   = SV_CLASSES),
        calls = expand("sv_calls/seed{seed}_size{sizerange}_vaf{vafrange}_{svclass}-{{binsize}}/{{binsize}}_fixed.{segments}/{{method}}.txt",
                       seed = SIMUL_SEEDS,
                       sizerange = SIZE_RANGES,
                       vafrange  = VAF_RANGES,
                       segments  = SEGMENTS,
                       svclass   = SV_CLASSES)
    output:
        "results/{binsize}_{method}.pdf",
        "results/{binsize}_{method}.pdf.txt"
    script:
        "utils/evaluation.stratified.R"


rule evaluation_meta:
    input:
        expand("results/{{binsize}}_{method}.pdf.txt",
               method = METHODS)
    output:
        "results/meta-{binsize}.pdf"
    script:
        "utils/evaluation.meta.R"

