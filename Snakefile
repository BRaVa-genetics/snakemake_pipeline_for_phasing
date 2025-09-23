"""
Snakefile pipeline for the phasing common and rare variants with SHAPEIT5

The working directory will be the one from where this script will be executed
The pipeline consists of 4 small steps (small conceptually, might be tough computationally)
    1) prepare BCF files (common & rare)
    2) phase common for the whole cohort (the scaffold)
    3) phase rare for the whole cohort, in chunks
    4) prepare and phase trios (just common variants; phase using the script for step-2)
    5) assess phasing and get files with the PP distribution

#### Required input ####
1. {TAG}.100trios.pedigree: a list of {child,father,mother} triplets for trio phasing
2. Several lists of sample IDs, such as mapping between array and WES (optional)
3. Genetic maps and lists of chunks, per chromosome (see the corresponding scripts)
4. A few parameters for Shapeit5 which are already set

Check the `README.md` for more details.

GK - Sep 26th, 2023, udpated Sept 2025
"""

## General input and parameters ##
TAG=config['tag']
WORK_DIR=config["work_dir"]
input_wes_prefix=config["input_wes_prefix"]
input_bed_prefix=config["input_bed_prefix"]
genet_map_prefix=config["genet_map_prefix"]
chunk_list_prefix=config["chunk_list_prefix"]
# PEDIGREE="{WORK_DIR}/sample_lists/{TAG}.100trios.pedigree"

rule prepare_bcfs:
    threads:5
    resources:
        mem_mb=16000,
    params:
        workdir = config["work_dir"],
        bed_prefix = config["input_bed_prefix"],
        pedigree = config["pedigree"]
    input:
        raw_wes = input_wes_prefix+"/chr{chrom}_hard_filters.tidy.vep.vcf.gz",
    output:
        bcf_common = "{wd}/sandbox/{tag}.notrios.common_merged.chr{chrom}.bcf",
        bcf_rare   = "{wd}/sandbox/{tag}.notrios.rare_prepared.chr{chrom}.bcf"
    run:
        shell("bash smk_01_prep_all_bcf.sh {wildcards.chrom} {wildcards.tag} {params.workdir} {input.raw_wes} {params.bed_prefix}")

rule phase_common:
    threads:10
    resources:
        mem_mb=16000,
    params:
        out_pref = config["work_dir"] + "/phased_common/{tag}.notrios"
    input:
        input_bcf = rules.prepare_bcfs.output.bcf_common,
        genet_map = genet_map_prefix + "/chr{chrom}.b38.gmap.gz"
    output:
        phased_bcf = "{wd}/phased_common/{tag}.notrios.chr{chrom}.bcf"
    shell:
        """
        bash smk_02_phase_common.sh {wildcards.chrom} {input.input_bcf} {input.genet_map} {params.out_pref}
        """

rule phase_rare:
    resources:
        mem_mb=16000
    input:
        input_bcf = rules.prepare_bcfs.output.bcf_rare,
        scaffold = rules.phase_common.output.phased_bcf
    params:
        workdir = config["work_dir"],
        genet_map = genet_map_prefix + "/chr{chrom}.b38.gmap.gz",
        chunk_list = chunk_list_prefix + "chr{chrom}.txt",
        tag = config["tag"] + ".notrios"
    output:
        phased_bcf = "{wd}/phased_rare/{tag}.notrios.phased_all.chr{chrom}.bcf"
    shell:
        """
        bash smk_03_phase_rare.sh {wildcards.chrom} {input.scaffold} {input.input_bcf} {params.tag} {params.genet_map} {params.chunk_list} {params.workdir}
        """

rule phase_trios:
    resources:
        mem_mb=16000
    threads:5
    params:
        workdir = config["work_dir"],
        bed_prefix = config["input_bed_prefix"],
        outpref = config["work_dir"] + "/phased_common/{tag}.trios",
        pedigree = config["pedigree"]
    input:
        raw_wes = input_wes_prefix + "/chr{chrom}_hard_filters.tidy.vep.vcf.gz",
        genet_map = genet_map_prefix + "/chr{chrom}.b38.gmap.gz"
    output:
        trios_prep = "{wd}/sandbox/{tag}.trios.prepared.chr{chrom}.bcf",
        trios_csi  = "{wd}/sandbox/{tag}.trios.prepared.chr{chrom}.bcf.csi",
        phased_bcf = "{wd}/phased_common/{tag}.trios.chr{chrom}.bcf",
        phased_csi = "{wd}/phased_common/{tag}.trios.chr{chrom}.bcf.csi"
    shell:
        """
        bash smk_04_phase_trios.sh {wildcards.chrom} {wildcards.tag} {params.workdir} {input.raw_wes} {params.bed_prefix}
        bash smk_02_phase_common.sh {wildcards.chrom} {output.trios_prep} {input.genet_map} {params.outpref} {params.pedigree}
        """

rule assess_phasing:
    input:
        phased_common= rules.phase_common.output.phased_bcf,
        phased_rare  = rules.phase_rare.output.phased_bcf,
        phased_trios = rules.phase_trios.output.phased_bcf
    output:
        ser_comn = "{wd}/phasing_assessment/{tag}.assess_common.chr{chrom}.variant.switch.txt.gz",
        ser_rare = "{wd}/phasing_assessment/{tag}.assess_rare.pp0.50.chr{chrom}.variant.switch.txt.gz"
    run:
        prefix = WORK_DIR + "/phasing_assessment/" + TAG
        pedigree = config["pedigree"]
        shell("bash smk_05_assess_phasing.sh common {wildcards.chrom} {input.phased_trios} {input.phased_common} {pedigree} {prefix}")
        shell("bash smk_05_assess_phasing.sh rare {wildcards.chrom} {input.phased_trios} {input.phased_rare} {pedigree} {prefix}")

rule get_pp_distribution:
    input:
        phased_rare  = rules.phase_rare.output.phased_bcf
    output:
        pp_mac0_pp90 = "{wd}/phasing_assessment/{tag}.pp.chr{chrom}.maf00015.pp0.90.gz",
        pp_mac0_pp50 = "{wd}/phasing_assessment/{tag}.pp.chr{chrom}.maf00015.pp0.50.gz"
    run:
        prefix = WORK_DIR + "/phasing_assessment/" + TAG
        shell("bash smk_05_assess_phasing.sh getpp {wildcards.chrom} {input.phased_rare} {prefix}.pp")


rule run_all:
    input:
        expand(WORK_DIR+"/phased_common/{tag}.trios.chr{chrom}.bcf", chrom=range(1,23), tag={TAG}),
        expand(WORK_DIR+"/phased_rare/{tag}.notrios.phased_all.chr{chrom}.bcf", chrom=range(1,23), tag={TAG}),   
        expand(WORK_DIR+"/phasing_assessment/{tag}.assess_rare.pp0.50.chr{chrom}.variant.switch.txt.gz", chrom=range(1,23), wd={WORK_DIR}, tag={TAG}),
        # expand(WORK_DIR+"/phasing_assessment/{tag}.assess_common.chr{chrom}.variant.switch.txt.gz", chrom=range(21,22), tag={TAG}),
        # expand(WORK_DIR+"/phasing_assessment/{tag}.pp.chr{chrom}.maf00015.pp0.90.gz", chrom=range(1,23), tag={TAG}),
