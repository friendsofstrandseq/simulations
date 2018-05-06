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
    install_MaRyam,
    convert_SVprob_output,
    sv_classifier_filter


NUM_SVS         = 150
NUM_CELLS       = 200
SIMUL_WINDOW    = [50000]
SIMUL_SEEDS_OLD = [1,2,3,4,5,6,7,8,9,10]
SIMUL_SEEDS_NEW = [5,6,7]
METHODS         = ["maryam",
                   "simple_llr0", "simple_llr2", "simple_llr4",
				   "merge_llr0",  "merge_llr2",  "merge_llr4",
                   "simple_llr2___size300000-vaf2", "simple_llr2___size300000-vaf4",
                   "simple_llr2___size500000-vaf2", "simple_llr2___size500000-vaf4",
				   "merge_llr2___size500000-vaf2", "merge_llr2___size500000-vaf4"]
CHROMOSOMES     = config['chromosomes']
SEGMENTS        = ["fraction15","fraction30","fraction55"]
ALL_SEGMENTS    = ["fraction08", "fraction15", \
                   "fraction20", "fraction25", "fraction30", "fraction35", \
                   "fraction45","fraction55","fraction65"]
SIZE_RANGES     = ["100000-200000", "200000-400000", "400000-800000", "800000-1600000", "1600000-3200000","3200000-6400000"]
VAF_RANGES      = ["1-5", "5-10", "10-20", "20-50", "50-100"]

rule simul:
    input:
        # Overview plots of cells
        # expand("plots/simulation{seed}-{binsize}/{binsize}_fixed.pdf",
        #         seed   = SIMUL_SEEDS,
        #         binsize = SIMUL_WINDOW),
        #
        # Plot SV calls together with simulated variants per chromosome
        # expand("sv_plots/simulation{seed}-{binsize}/{binsize}_fixed.{segments}/{method}.{chrom}.pdf",
        #         seed      = SIMUL_SEEDS,
        #         binsize   = SIMUL_WINDOW,
        #         segments  = SEGMENTS,
        #         chrom     = CHROMOSOMES,
        #         method    = METHODS),
        #
        # Evaluation curve that I designed during the hackathon
        # (based on a simulation of mixed SV sizes and VAFs)
        expand("results/evaluation_hackathon/{binsize}_{method}.pdf",
                binsize = [50000],
                method  = METHODS),
        #
        # New SV evaluation curves
        # (based on separate simulations for SV sizes and VAFs)
        expand("results/evaluation_stratified/{binsize}_{method}.pdf",
                binsize = [50000],
                method  = METHODS),
        #
        # Plot SV calls of some of the new simulations
        expand("sv_plots/seed{seed}_size{sizerange}_vaf{vafrange}-{binsize}/{binsize}_fixed.{segments}/{method}.{chrom}.pdf",
                       seed = [5],
                       sizerange = SIZE_RANGES,
                       vafrange  = VAF_RANGES,
                       segments  = SEGMENTS,
                       binsize   = SIMUL_WINDOW,
                       method    = METHODS,
                       chrom     = CHROMOSOMES),
        expand("results/evaluation_stratified/meta-{binsize}.pdf",
                binsize = SIMUL_WINDOW)



################################################################################
# Simulation of count data                                                     #
################################################################################

rule simulate_genome:
    output:
        tsv="simulation/genome/genome{seed}.tsv"
    log:
        "log/simulate_genome/genome{seed}.tsv"
    params:
        svcount     = NUM_SVS,
        minsize     =  100000,
        maxsize     = 5000000,
        mindistance = 1000000,
    shell:
        "Rscript utils/simulate_SVs.R {wildcards.seed} {params.svcount} {params.minsize} {params.maxsize} {params.mindistance} {output.tsv} > {log} 2>&1"

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
        cell_count   = NUM_CELLS,
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
# New Simulation of count data                                                 #
################################################################################

rule new_simulate_genome:
    output:
        "simulation_new/seed{seed}_size{minsvsize}-{maxsvsize}_vaf{minvaf}-{maxvaf}/genome.tsv"
    log:
        "log/simulation_new/seed{seed}_size{minsvsize}-{maxsvsize}_vaf{minvaf}-{maxvaf}.txt"
    params:
        svcount     = NUM_SVS,
        mindistance = 1000000
    shell:
        """
        tmp=$(mktemp)
        Rscript utils/simulate_SVs.R {wildcards.seed} {params.svcount} \
            {wildcards.minsvsize} {wildcards.maxsvsize} {params.mindistance} $tmp \
            > {log} 2>&1
        awk -v minvaf={wildcards.minvaf} -v maxvaf={wildcards.maxvaf} \
            -v seed={wildcards.seed} \
            'BEGIN {{srand(seed); OFS="\\t"}}
            {{ vaf = int(minvaf + 1.0 * rand() * (maxvaf - minvaf))/100.0; print $0, vaf}}' $tmp \
        > {output} 2> {log}
        """

rule new_simulate_counts:
    input:
        config = "simulation_new/seed{seed}_size{minsvsize}-{maxsvsize}_vaf{minvaf}-{maxvaf}/genome.tsv"
    output:
        counts   = "simulation_new/seed{seed}_size{minsvsize}-{maxsvsize}_vaf{minvaf}-{maxvaf}/counts-{binsize}.txt.gz",
        segments = "simulation_new/seed{seed}_size{minsvsize}-{maxsvsize}_vaf{minvaf}-{maxvaf}/segments-{binsize}.txt",
        phases   = "simulation_new/seed{seed}_size{minsvsize}-{maxsvsize}_vaf{minvaf}-{maxvaf}/phases-{binsize}.txt",
        info     = "simulation_new/seed{seed}_size{minsvsize}-{maxsvsize}_vaf{minvaf}-{maxvaf}/info-{binsize}.txt",
        sce      = "simulation_new/seed{seed}_size{minsvsize}-{maxsvsize}_vaf{minvaf}-{maxvaf}/sce-{binsize}.txt",
        variants = "simulation_new/seed{seed}_size{minsvsize}-{maxsvsize}_vaf{minvaf}-{maxvaf}/variants-{binsize}.txt",
    params:
        mc_command   = config["mosaicatcher"],
        neg_binom_p  = neg_binom_p,
        min_coverage = min_coverage,
        max_coverage = max_coverage,
        cell_count   = NUM_CELLS,
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
            --sample-name seed{wildcards.seed}_size{wildcards.minsvsize}-{wildcards.maxsvsize}_vaf{wildcards.minvaf}-{wildcards.maxvaf}-{wildcards.binsize} \
            {input.config} > {log} 2>&1
        """

rule new_link_to_simulated_counts:
    input:
        counts = "simulation_new/seed{seed}_size{minsvsize}-{maxsvsize}_vaf{minvaf}-{maxvaf}/counts-{binsize}.txt.gz",
        info   = "simulation_new/seed{seed}_size{minsvsize}-{maxsvsize}_vaf{minvaf}-{maxvaf}/info-{binsize}.txt",
    output:
        counts = "counts/seed{seed}_size{minsvsize}-{maxsvsize}_vaf{minvaf}-{maxvaf}-{binsize}/{binsize}_fixed.txt.gz",
        info   = "counts/seed{seed}_size{minsvsize}-{maxsvsize}_vaf{minvaf}-{maxvaf}-{binsize}/{binsize}_fixed.info"
    run:
        d = os.path.dirname(output.counts)
        count_file = os.path.basename(output.counts)
        info_file  = os.path.basename(output.info)
        shell("cd {d} && ln -s ../../{input.counts} {count_file} && ln -s ../../{input.info} {info_file} && cd ../..")


rule new_link_to_simulated_strand_states:
    input:
        sce    = "simulation_new/seed{seed}_size{minsvsize}-{maxsvsize}_vaf{minvaf}-{maxvaf}/sce-{binsize}.txt",
    output:
        states = "strand_states/seed{seed}_size{minsvsize}-{maxsvsize}_vaf{minvaf}-{maxvaf}-{binsize}/final.txt"
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

rule plot_SV_calls_old:
    input:
        counts = "counts/simulation{seed}-{binsize}/{binsize}_fixed.txt.gz",
        segs   = "segmentation2/simulation{seed}-{binsize}/{binsize}_fixed.{segments}.txt",
        svs    = "sv_calls/simulation{seed}-{binsize}/{binsize}_fixed.{segments}/{method}.txt",
        simul  = "simulation/variants/genome{seed}-{binsize}.txt"
    output:
        expand("sv_plots/simulation{{seed}}-{{binsize}}/{{binsize}}_fixed.{{segments}}/{{method}}.{chrom}.pdf", chrom = config['chromosomes'])
    script:
        "utils/plot_sv_calls.R"

rule plot_SV_calls_new:
    input:
        counts = "counts/seed{seed}_size{minsvsize}-{maxsvsize}_vaf{minvaf}-{maxvaf}-{binsize}/{binsize}_fixed.txt.gz",
        segs   = "segmentation2/seed{seed}_size{minsvsize}-{maxsvsize}_vaf{minvaf}-{maxvaf}-{binsize}/{binsize}_fixed.{segments}.txt",
        svs    = "sv_calls/seed{seed}_size{minsvsize}-{maxsvsize}_vaf{minvaf}-{maxvaf}-{binsize}/{binsize}_fixed.{segments}/{method}.txt",
        simul  = "simulation_new/seed{seed}_size{minsvsize}-{maxsvsize}_vaf{minvaf}-{maxvaf}/variants-{binsize}.txt"
    output:
        expand("sv_plots/seed{{seed}}_size{{minsvsize}}-{{maxsvsize}}_vaf{{minvaf}}-{{maxvaf}}-{{binsize}}/{{binsize}}_fixed.{{segments}}/{{method}}.{chrom}.pdf", chrom = config['chromosomes'])
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
        -m 0.5 -M 50000000 -o {output} \
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

rule install_MaRyam:
    output:
        "utils/R-packages2/MaRyam/R/MaRyam"
    log:
        "log/install_MaRyam.log"
    shell:
        """
        TAR=$(which tar) Rscript utils/install_maryam.R > {log} 2>&1
        """

rule sv_classifier_maryam:
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

rule sv_classifier_merger:
    input:
        prob = "sv_probabilities/{sample}/{windows}.{bpdens}/raw_probabilities.Rdata"
    output:
        "sv_calls/{sample}/{windows}.{bpdens}/merge_llr{llr}.txt"
    params:
        llr_cutoff = lambda wc: wc.llr
    script:
        "utils/sv_classifier_merger.R"


rule sv_classifier_biallelic:
    input:
        "sv_probabilities/{sample}/{windows}.{bpdens}/raw_probabilities.Rdata"
    output:
        prob = "sv_calls/{sample}/{windows}.{bpdens}/biallelic_min{mincells}cells.txt"
    params:
        mincells = lambda wc: wc.mincells
    script:
        "utils/sv_classifier_biallelic.R"


rule sv_classifier_filter:
    input:
        "sv_calls/{sample}/{windows}.{bpdens}/{method}.txt"
    output:
        "sv_calls/{sample}/{windows}.{bpdens}/{method}___size{minsvsize}-vaf{minvaf}.txt"
    params:
        minsvsize = lambda wc: wc.minsvsize,
        minvaf    = lambda wc: wc.minvaf,
        numcells  = NUM_CELLS
    script:
        "utils/filter_sv_calls.R"




################################################################################
# Evaluation of MosaiCatcher                                                   #
################################################################################


rule evaluation_hackathon:
    input:
        truth = expand("simulation/variants/genome{seed}-{{binsize}}.txt", \
                       seed = SIMUL_SEEDS_OLD),
        calls = expand("sv_calls/simulation{seed}-{{binsize}}/{{binsize}}_fixed.{segments}/{{method}}.txt", \
                       seed = SIMUL_SEEDS_OLD, segments = ALL_SEGMENTS)
    params:
        N_cells = NUM_CELLS
    output:
        "results/evaluation_hackathon/{binsize}_{method}.pdf"
    script:
        "utils/evaluation.R"

rule evaluation_newer_version:
    input:
        truth = expand("simulation_new/seed{seed}_size{sizerange}_vaf{vafrange}/variants-{{binsize}}.txt",
                       seed = SIMUL_SEEDS_NEW,
                       sizerange = SIZE_RANGES,
                       vafrange  = VAF_RANGES),
        calls = expand("sv_calls/seed{seed}_size{sizerange}_vaf{vafrange}-{{binsize}}/{{binsize}}_fixed.{segments}/{{method}}.txt",
                       seed = SIMUL_SEEDS_NEW,
                       sizerange = SIZE_RANGES,
                       vafrange  = VAF_RANGES,
                       segments = ALL_SEGMENTS),
    output:
        "results/evaluation_stratified/{binsize}_{method}.pdf",
        "results/evaluation_stratified/{binsize}_{method}.pdf.txt"
    script:
        "utils/evaluation.stratified.R"


rule evaluation_meta:
    input:
        expand("results/evaluation_stratified/{{binsize}}_{method}.pdf.txt",
               method = METHODS)
    output:
        "results/evaluation_stratified/meta-{binsize}.pdf"
    script:
        "utils/evaluation.meta.R"
