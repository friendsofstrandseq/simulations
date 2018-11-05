configfile: "Snake.config.json"
import os.path

wildcard_constraints:
    binsize   = "\d+",
    seed      = "\d+",
    windows   = "\d+_[a-zA-Z]+",
    chrom     = "[chr0-9XY]+",
    segments  = "[a-zA-Z0-9]+",
    method    = "[a-zA-Z0-9_\.]+",
    minsvsize = "\d+",
    maxsvsize = "\d+",
    minvaf    = "\d+",
    maxvaf    = "\d+",
    mincells  = "\d+",
    llr       = "\d"


localrules:
    simul,
    new_simulate_genome,
    new_link_to_simulated_counts,
    new_link_to_simulated_strand_states,
    prepare_segments,
    install_mosaiClassifier,
	evaluation_newer_version,
    evaluation_meta



### Global settings
NUM_SVS         = config["num_SVs_per_genome"]
CHROM           = config["chromosomes"]
SIMUL_SEEDS     = config["seeds"]
SIMUL_WINDOW    = config["window_sizes"]
NUM_CELLS       = config["num_cells"]
SIZE_RANGES     = config["size_ranges"]
VAF_RANGES      = config["vaf_ranges"]
SV_CLASSES      = config["sv_classes"]
METHODS         = ["simpleCalls_llr4_poppriorsTRUE_haplotagsFALSE_gtcutoff0.05_regfactor6_filterTRUE"]
SEGMENTS        = ["fraction100"]


### Main rule
rule simul:
    input:
        expand("results/meta-{binsize}.pdf", binsize = SIMUL_WINDOW),

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
        "log/simulate_counts/seed{seed}_size{minsvsize}-{maxsvsize}_vaf{minvaf}-{maxvaf}_{svclass}-{binsize}.txt"
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
# Segmentation                                                                 #
################################################################################

rule segmentation:
    input:
        "counts/{sample}/{file_name}.txt.gz"
    output:
        "segmentation/{sample}/{file_name}.txt.fixme"
    log:
        "log/segmentation/{sample}/{file_name}.log"
    params:
        mc_command = config["mosaicatcher"]
    shell:
        """
        {params.mc_command} segment \
        --remove-none \
        --forbid-small-segments 4 \
        -M 50000000 \
        -o {output} \
        {input} > {log} 2>&1
        """

# TODO: This is a workaround because latest versions of "mosaic segment" don't compute the "bps"
# TODO: column properly. Remove once fixed in the C++ code.
rule fix_segmentation:
    input:
        "segmentation/{sample}/{window}_{file_name}.txt.fixme"
    output:
        "segmentation/{sample}/{window,\d+}_{file_name}.txt"
    shell:
        'awk \'BEGIN {{OFS="\\t"}} {{if ($1=="{wildcards.sample}") $12=int(($14-1)/100000); print}}\' {input} > {output}'

# Pick a few segmentations and prepare the input files for SV classification
rule prepare_segments:
    input:
        "segmentation/{sample}/{windows}.txt",
    output:
        "segmentation2/{sample}/{windows}.fraction100.txt"
    log:
        "log/prepare_segments/{sample}/{windows}.fraction100.log"
    shell:
        './utils/select_segmentation.py --output_jointseg {output} {input}'


################################################################################
# SV classification                                                            #
################################################################################

rule install_mosaiClassifier:
    output:
        "utils/mosaiClassifier.snakemake.R",
        "utils/helper.prepare_segments.R"
    shell:
        """
        git clone https://github.com/friendsofstrandseq/pipeline.git mosaiClassifier
        cd mosaiClassifier
        cd ../utils
        ln -s ../mosaiClassifier/utils/mosaiClassifier
        ln -s ../mosaiClassifier/utils/mosaiClassifier.snakemake.R
        ln -s ../mosaiClassifier/utils/mosaiClassifier_call.snakemake.R
        ln -s ../mosaiClassifier/utils/mosaiClassifier_call_biallelic.snakemake.R
        ln -s ../mosaiClassifier/utils/helper.prepare_segments.R
	ln -s ../mosaiClassifier/utils/sv_consistency_barplot.R
	ln -s ../mosaiClassifier/utils/sv_consistency_barplot.snakemake.R
	ln -s ../mosaiClassifier/utils/select_segmentation.py
	ln -s ../mosaiClassifier/utils/plot-sv-calls.R
        """

rule mosaiClassifier_make_call:
    input:
        probs = "sv_probabilities/{sample}/{windows}.{bpdens}/probabilities.Rdata"
    output:
        "sv_calls/{sample}/{windows}.{bpdens}/simpleCalls_llr{llr}_poppriors{pop_priors,(TRUE|FALSE)}_haplotags{use_haplotags, (TRUE|FALSE)}_gtcutoff{gtcutoff,[0-9\\.]+}_regfactor{regfactor,[0-9]+}_filter{use_filter, (TRUE|FALSE)}.txt"
    log:
        "log/mosaiClassifier_make_call/{sample}/{windows}.{bpdens}.llr{llr}.poppriors{pop_priors}.haplotags{use_haplotags}.gtcutoff{gtcutoff}.regfactor{regfactor}.filter{use_filter}.log"
    script:
        "utils/mosaiClassifier_call.snakemake.R"

rule mosaiClassifier_calc_probs:
    input:
        counts = "counts/{sample}/{windows}.txt.gz",
        info   = "counts/{sample}/{windows}.info",
        states = "strand_states/{sample}/final.txt",
        bp     = "segmentation2/{sample}/{windows}.{bpdens}.txt"
    output:
        output = "sv_probabilities/{sample}/{windows}.{bpdens}/probabilities.Rdata"
    log:
        "log/mosaiClassifier_calc_probs/{sample}/{windows}.{bpdens}.log"
    script:
        "utils/mosaiClassifier.snakemake.R"

rule mosaiClassifier_make_call_biallelic:
    input:
        probs = "sv_probabilities/{sample}/{windows}.{bpdens}/probabilities.Rdata"
    output:
	"sv_calls/{sample}/{windows}.{bpdens}/biAllelic_llr{llr}_poppriors{poppriors}_haplotags{use_haplotags, (TRUE|FALSE)}_gtcutoff{gtcutoff,[0-9\\.]+}_regfactor{regfactor}_filter{use_filter, (TRUE|FALSE)}.txt"
    log:
	"log/mosaiClassifier_make_call_biallelic/{sample}/{windows}.{bpdens}.llr{llr}.poppriors{poppriors}.haplotags{use_haplotags}.gtcutoff{gtcutoff}.regfactor{regfactor}.filter{use_filter}.log"
    script:
        "utils/mosaiClassifier_call_biallelic.snakemake.R"





################################################################################
# Evaluation of MosaiCatcher                                                   #
################################################################################


rule evaluation_newer_version:
    input:
        truth = expand("simulation_new/seed{seed}_size{sizerange}_vaf{vafrange}_{svclass}/variants-{binsize}.txt",
                       seed = SIMUL_SEEDS,
                       sizerange = SIZE_RANGES,
                       vafrange  = VAF_RANGES,
                       svclass   = SV_CLASSES,
                       binsize   = SIMUL_WINDOW),
        calls = expand("sv_calls/seed{seed}_size{sizerange}_vaf{vafrange}_{svclass}-{binsize}/{binsize}_fixed.{segments}/{method}.txt",
                       seed = SIMUL_SEEDS,
                       sizerange = SIZE_RANGES,
                       vafrange  = VAF_RANGES,
                       segments  = SEGMENTS,
                       svclass   = SV_CLASSES,
                       binsize = SIMUL_WINDOW,
                       method = METHODS)
    output:
        table      = "results/{binsize, [0-9]+}_{method}.pdf.txt",
    log:
        "log/results/{binsize, [0-9]+}_{method}.log"
    script:
        "utils/evaluation.stratified.R"


rule evaluation_meta:
    input:
        expand("results/{binsize}_{method}.pdf.txt",
               binsize = SIMUL_WINDOW,
               method = METHODS)
    output:
        "results/meta-{binsize, [0-9]+}.pdf"
    script:
        "utils/evaluation.meta.R"


