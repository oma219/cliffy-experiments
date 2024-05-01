##################################################
# Name: exp_2.smk
# Description: Contains the workflow and methods
#              needed for experiment 2.
#
#              This experiment builds the 16S
#              rRNA read dataset using 
#              MicrobeMixer code.
#
# Date: Mar 18th, 2024
##################################################

####################################################
# Section 1: Helper functions needed for these
#            experiment rules.
####################################################

def get_input_fasta_file_exp2(wildcards):
    """ Returns path to SILVA database with taxonomic path in each header """
    file_list = []
    for data_file in os.listdir(f"exp2_data/"):
        if data_file.endswith(".fasta") and data_file.startswith("SILVA"):
            file_list.append(f"exp2_data/" + data_file)
    assert len(file_list) == 1
    return file_list[0]

def get_input_taxonomy_exp2(wildcards):
    """ Returns path to SILVA database with taxonomic path in each header """
    file_list = []
    for data_file in os.listdir(f"exp2_data/"):
        if data_file.endswith(".txt") and data_file.startswith("tax_slv_ssu"):
            file_list.append(f"exp2_data/" + data_file)
    assert len(file_list) == 1
    return file_list[0]

def get_input_primers_exp2(wildcards):
    """ Returns path to primers file """
    file_list = []
    for data_file in os.listdir(f"exp2_data/"):
        if data_file.endswith("_primers.txt") and data_file.startswith(wildcards.primer):
            file_list.append(f"exp2_data/" + data_file)
    assert len(file_list) == 1
    return file_list[0]

####################################################
# Section 2: Rules needed for this experiment type
####################################################

# Section 2.1: Generate biome abundance files for each biome

rule generate_biome_abundance_files_exp2:
    input:
        get_input_taxonomy_exp2
    output:
        "exp2_abundance_files/{biome}/abundance.csv",
        "exp2_abundance_files/{biome}/stdout.log"
    shell:
        """
        # get the EBI lineage string
        lineage=""
        if [ "{wildcards.biome}" == "waste_water" ]; then
            lineage="root:Engineered:Wastewater"
        elif [ "{wildcards.biome}" == "human_gut" ]; then
            lineage="root:Host-associated:Human:Digestive system"
        elif [ "{wildcards.biome}" == "aquatic" ]; then
            lineage="root:Environmental:Aquatic"
        elif [ "{wildcards.biome}" == "soil" ]; then
            lineage="root:Environmental:Terrestrial:Soil"
        else
            exit 1
        fi

        # run command
        python3 {microbe_mixer_prog} stats --biome "$lineage" --taxonomy {input} --output {output[0]} > {output[1]}
        """

# Section 2.2: Simulate reads for each biome for all primers

rule generate_reads_for_biome_exp2:
    input:
        get_input_fasta_file_exp2,
        get_input_taxonomy_exp2,
        get_input_primers_exp2,
        "exp2_abundance_files/{biome}/abundance.csv"
    output:
        "exp2_reads/{biome}/{primer}/{primer}_mate_1.fq",
        "exp2_reads/{biome}/{primer}/{primer}_mate_2.fq",
        "exp2_reads/{biome}/{primer}/stdout.log"
    shell:
        """
        # run command
        python3 {microbe_mixer_prog} simulate \
                --biome-abundance exp2_abundance_files/{wildcards.biome}/abundance.csv \
                --primers {input[2]} \
                --silva-taxonomy {input[1]} \
                --silva-ref {input[0]} \
                --num-reads {num_reads_to_simulate_exp2} \
                --output-name {wildcards.primer} \
                --temp-dir exp2_reads/{wildcards.biome}/{wildcards.primer}/ > {output[2]}
        
        # remove extra files
        rm exp2_reads/{wildcards.biome}/{wildcards.primer}/genus_*_reads_1.fq 
        rm exp2_reads/{wildcards.biome}/{wildcards.primer}/genus_*reads_2.fq 
        rm exp2_reads/{wildcards.biome}/{wildcards.primer}/genus_*_seqs.fna
        """

# Section 2.3: Overall rule to run experiment 2

rule run_exp2:
    input:
        expand("exp2_reads/{biome}/{primer}/{primer}_mate_1.fq", primer=["V1_V2", "V3_V4", "V4_V4", "V4_V5"], biome=["waste_water", "human_gut", "aquatic", "soil"])
