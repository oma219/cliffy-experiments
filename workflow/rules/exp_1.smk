##################################################
# Name: exp_1.smk
# Description: Contains the workflow and methods
#              needed for experiment 1.
#
#              This experiment builds Cliffy 
#              indexes over SILVA. For 
#              Cliffy, it builds the index
#              with and without taxonomic compression.
#
#              Output are statistics on index 
#              size for Cliffy, as well as analysis 
#              over the number of monotonic increases 
#              in each profile.
#
# Date: Mar 18th, 2024
##################################################

####################################################
# Section 1: Helper functions needed for these
#            experiment rules.
####################################################

def get_input_fasta_file_exp1(wildcards):
    """ Returns path to SILVA database with taxonomic path in each header """
    file_list = []
    for data_file in os.listdir(f"exp1_data/"):
        if data_file.endswith(".fasta") and data_file.startswith("SILVA"):
            file_list.append(f"exp1_data/" + data_file)
    assert len(file_list) == 1
    return file_list[0]

####################################################
# Section 2: Rules needed for this experiment type
####################################################

##############################################################
# Section 2.1: Generate FASTA files for every genera in SILVA
##############################################################
rule convert_to_uppercase_and_samelines_exp1:
    input:
        get_input_fasta_file_exp1
    output:
        "exp1_silva_database/silva_database.fa"
    shell:
        "seqtk seq -U {input} > {output}"

rule seperate_SILVA_ref_into_genera_exp1:
    input:
        "exp1_silva_database/silva_database.fa"
    output:
        "exp1_rna_input_files/filelist.txt"
    shell:
        """
        python3 {repo_dir}/src/write_silva_genera.py \
        -i {input} \
        -o {work_dir}/exp1_rna_input_files/ \
        --tree exp1_data/tax_slv_ssu_138.1.tre \
        --tree-rank exp1_data/tax_slv_ssu_138.1.txt \
        -n {num_genera_exp1}
        """

rule generate_dna_version_of_fasta_files_exp1:
    input:
        "exp1_rna_input_files/filelist.txt"
    output:
        "exp1_dna_input_files/filelist.txt"
    shell:
        """
        work_dir=$(pwd)
        files=(exp1_rna_input_files/*.fa)
        num_docs=${{#files[@]}}
        
        touch {output}
        for i in $(seq 1 $num_docs); do
            input_file="${{work_dir}}/exp1_rna_input_files/doc_${{i}}_seq.fa"
            output_file="${{work_dir}}/exp1_dna_input_files/doc_${{i}}_seq.fa"

            seqtk seq -r $input_file > $output_file
            printf "%s %d\n" $output_file $i >> {output}
        done
        """

########################################################
# Section 2.2: Build the indexes for Cliffy on full text
########################################################

rule build_index_with_two_pass_full_dap_exp1:
    input:
        "exp1_dna_input_files/filelist.txt"
    output:
        "exp1_index/full_dap_full_text/filelist.txt",
        "exp1_index/full_dap_full_text/build.log",
        "exp1_index/full_dap_full_text/output.fna.sdap",
        "exp1_index/full_dap_full_text/output.fna.edap",
        "exp1_index/full_dap_full_text/time_and_mem.log"
    shell:
        """
        cp {input} {output[0]}
        {time_prog} {time_format} --output={output[4]} \
        pfp_doc64 build --filelist {output[0]} \
                        --output exp1_index/full_dap_full_text/output \
                        --revcomp \
                        --no-ftab \
                        --two-pass exp1_index/full_dap_full_text/temp  \
                        --tmp-size {tmp_mem_used_exp1} 2> {output[1]}
        """

rule build_index_with_two_pass_taxonomic_compressed_exp1:
    input:
        "exp1_dna_input_files/filelist.txt"
    output:
        "exp1_index/tax_compressed_full_text/filelist.txt",
        "exp1_index/tax_compressed_full_text/build.log",
        "exp1_index/tax_compressed_full_text/output.fna.taxcomp.sdap",
        "exp1_index/tax_compressed_full_text/output.fna.taxcomp.edap",
        "exp1_index/tax_compressed_full_text/output.fna.taxcomp.of.sdap",
        "exp1_index/tax_compressed_full_text/output.fna.taxcomp.of.edap",
        "exp1_index/tax_compressed_full_text/time_and_mem.log"
    shell:
        """
        cp {input} {output[0]}
        {time_prog} {time_format} --output={output[6]} \
        pfp_doc64 build --filelist {output[0]} \
                        --output exp1_index/tax_compressed_full_text/output \
                        --revcomp \
                        --taxcomp \
                        --num-col 7 \
                        --two-pass exp1_index/tax_compressed_full_text/temp  \
                        --tmp-size {tmp_mem_used_exp1} 2> {output[1]}
        """

#############################################################################
# Section 2.3: Build the indexes for Cliffy on DNA minimizer version of text
#############################################################################

rule build_index_with_two_pass_full_dap_with_dna_minimizers_exp1:
    input:
        "exp1_dna_input_files/filelist.txt"
    output:
        "exp1_index/full_dap_dna_minimizers/filelist.txt",
        "exp1_index/full_dap_dna_minimizers/build.log",
        "exp1_index/full_dap_dna_minimizers/output.fna.sdap",
        "exp1_index/full_dap_dna_minimizers/output.fna.edap",
        "exp1_index/full_dap_dna_minimizers/time_and_mem.log"
    shell:
        """
        cp {input} {output[0]}
        {time_prog} {time_format} --output={output[4]} \
        pfp_doc64 build --filelist {output[0]} \
                        --output exp1_index/full_dap_dna_minimizers/output \
                        --revcomp \
                        --no-ftab \
                        --dna-minimizers  \
                        --small-window 4 \
                        --large-window 11 \
                        --two-pass exp1_index/full_dap_dna_minimizers/temp  \
                        --tmp-size {tmp_mem_used_exp1} 2> {output[1]}
        """

rule build_index_with_two_pass_taxonomic_compressed_with_dna_minimizers_exp1:
    input:
        "exp1_dna_input_files/filelist.txt"
    output:
        "exp1_index/tax_compressed_dna_minimizers/filelist.txt",
        "exp1_index/tax_compressed_dna_minimizers/build.log",
        "exp1_index/tax_compressed_dna_minimizers/output.fna.taxcomp.sdap",
        "exp1_index/tax_compressed_dna_minimizers/output.fna.taxcomp.edap",
        "exp1_index/tax_compressed_dna_minimizers/output.fna.taxcomp.of.sdap",
        "exp1_index/tax_compressed_dna_minimizers/output.fna.taxcomp.of.edap",
        "exp1_index/tax_compressed_dna_minimizers/time_and_mem.log"
    shell:
        """
        cp {input} {output[0]}
        {time_prog} {time_format} --output={output[6]} \
        pfp_doc64 build --filelist {output[0]} \
                        --output exp1_index/tax_compressed_dna_minimizers/output \
                        --revcomp \
                        --taxcomp \
                        --num-col 7 \
                        --dna-minimizers  \
                        --small-window 4 \
                        --large-window 11 \
                        --two-pass exp1_index/tax_compressed_dna_minimizers/temp  \
                        --tmp-size {tmp_mem_used_exp1} 2> {output[1]}
        """

###################################################################
# Section 2.2: Build the indexes for Cliffy on minimizer digested
###################################################################

rule build_index_with_two_pass_full_dap_with_minimizers_exp1:
    input:
        "exp1_dna_input_files/filelist.txt"
    output:
        "exp1_index/full_dap_minimizers/filelist.txt",
        "exp1_index/full_dap_minimizers/build.log",
        "exp1_index/full_dap_minimizers/output.fna.sdap",
        "exp1_index/full_dap_minimizers/output.fna.edap",
        "exp1_index/full_dap_minimizers/time_and_mem.log"
    shell:
        """
        cp {input} {output[0]}
        {time_prog} {time_format} --output={output[4]} \
        pfp_doc64 build --filelist {output[0]} \
                        --output exp1_index/full_dap_minimizers/output \
                        --revcomp \
                        --no-ftab \
                        --minimizers  \
                        --small-window 4 \
                        --large-window 11 \
                        --two-pass exp1_index/full_dap_minimizers/temp  \
                        --tmp-size {tmp_mem_used_exp1} 2> {output[1]}
        """

rule build_index_with_two_pass_taxonomic_compressed_with_minimizers_exp1:
    input:
        "exp1_dna_input_files/filelist.txt"
    output:
        "exp1_index/tax_compressed_minimizers/filelist.txt",
        "exp1_index/tax_compressed_minimizers/build.log",
        "exp1_index/tax_compressed_minimizers/output.fna.taxcomp.sdap",
        "exp1_index/tax_compressed_minimizers/output.fna.taxcomp.edap",
        "exp1_index/tax_compressed_minimizers/output.fna.taxcomp.of.sdap",
        "exp1_index/tax_compressed_minimizers/output.fna.taxcomp.of.edap",
        "exp1_index/tax_compressed_minimizers/time_and_mem.log"
    shell:
        """
        cp {input} {output[0]}
        {time_prog} {time_format} --output={output[6]} \
        pfp_doc64 build --filelist {output[0]} \
                        --output exp1_index/tax_compressed_minimizers/output \
                        --revcomp \
                        --taxcomp \
                        --num-col 5 \
                        --minimizers  \
                        --small-window 4 \
                        --large-window 11 \
                        --two-pass exp1_index/tax_compressed_minimizers/temp  \
                        --tmp-size {tmp_mem_used_exp1} 2> {output[1]}
        """

##############################################################################
# Section 2.3: Build the full DAP indexes for Cliffy using no reverse complement 
#              in order to load the full DAPs and print them in csv format
##############################################################################

rule build_index_with_two_pass_full_dap_no_ftab_no_revcomp_exp1:
    input:
        "exp1_dna_input_files/filelist.txt"
    output:
        "exp1_index_no_revcomp/full_dap_full_text/filelist.txt",
        "exp1_index_no_revcomp/full_dap_full_text/build.log",
        "exp1_index_no_revcomp/full_dap_full_text/output.fna.sdap",
        "exp1_index_no_revcomp/full_dap_full_text/output.fna.edap",
        "exp1_index_no_revcomp/full_dap_full_text/time_and_mem.log"
    shell:
        """
        cp {input} {output[0]}
        {time_prog} {time_format} --output={output[4]} \
        pfp_doc64 build --filelist {output[0]} \
                        --output exp1_index_no_revcomp/full_dap_full_text/output \
                        --no-ftab \
                        --two-pass exp1_index_no_revcomp/full_dap_full_text/temp  \
                        --tmp-size {tmp_mem_used_exp1} 2> {output[1]}
        """
rule build_index_with_two_pass_full_dap_with_dna_minimizers_no_ftab_no_revcomp_exp1:
    input:
        "exp1_dna_input_files/filelist.txt"
    output:
        "exp1_index_no_revcomp/full_dap_dna_minimizers/filelist.txt",
        "exp1_index_no_revcomp/full_dap_dna_minimizers/build.log",
        "exp1_index_no_revcomp/full_dap_dna_minimizers/output.fna.sdap",
        "exp1_index_no_revcomp/full_dap_dna_minimizers/output.fna.edap",
        "exp1_index_no_revcomp/full_dap_dna_minimizers/time_and_mem.log"
    shell:
        """
        cp {input} {output[0]}
        {time_prog} {time_format} --output={output[4]} \
        pfp_doc64 build --filelist {output[0]} \
                        --output exp1_index_no_revcomp/full_dap_dna_minimizers/output \
                        --no-ftab \
                        --dna-minimizers  \
                        --small-window 4 \
                        --large-window 11 \
                        --two-pass exp1_index_no_revcomp/full_dap_dna_minimizers/temp  \
                        --tmp-size {tmp_mem_used_exp1} 2> {output[1]}
        """

rule build_index_with_two_pass_full_dap_with_minimizers_no_ftab_no_revcomp_exp1:
    input:
        "exp1_dna_input_files/filelist.txt"
    output:
        "exp1_index_no_revcomp/full_dap_minimizers/filelist.txt",
        "exp1_index_no_revcomp/full_dap_minimizers/build.log",
        "exp1_index_no_revcomp/full_dap_minimizers/output.fna.sdap",
        "exp1_index_no_revcomp/full_dap_minimizers/output.fna.edap",
        "exp1_index_no_revcomp/full_dap_minimizers/time_and_mem.log"
    shell:
        """
        cp {input} {output[0]}
        {time_prog} {time_format} --output={output[4]} \
        pfp_doc64 build --filelist {output[0]} \
                        --output exp1_index_no_revcomp/full_dap_minimizers/output \
                        --no-ftab \
                        --minimizers  \
                        --small-window 4 \
                        --large-window 11 \
                        --two-pass exp1_index_no_revcomp/full_dap_minimizers/temp  \
                        --tmp-size {tmp_mem_used_exp1} 2> {output[1]}
        """

################################################################
# Section 2.4: Analyze the number of increases in each profile
################################################################

rule print_out_all_document_array_profiles_exp1:
    input:
        "exp1_index_no_revcomp/full_dap_{digestion_type}/filelist.txt",
        "exp1_index_no_revcomp/full_dap_{digestion_type}/output.fna.sdap",
        "exp1_index_no_revcomp/full_dap_{digestion_type}/output.fna.edap"
    output:
        "exp1_results/full_dap_{digestion_type}_csv/output.sdap.csv",
        "exp1_results/full_dap_{digestion_type}_csv/output.edap.csv"
    shell:
        """
        pfp_doc64 info --ref exp1_index_no_revcomp/full_dap_{wildcards.digestion_type}/output \
                       --output exp1_results/full_dap_{wildcards.digestion_type}_csv/output
        """

rule analyze_document_array_profiles_exp1:
    input:
        "exp1_results/full_dap_{digestion_type}_csv/output.sdap.csv",
        "exp1_results/full_dap_{digestion_type}_csv/output.edap.csv"
    output:
        "exp1_results/full_dap_{digestion_type}_analysis/increases.csv",
        "exp1_results/full_dap_{digestion_type}_analysis/sig_vals.csv"
    shell:
        """
        {dap_read_prog} -i {input[0]} -o {output[0]} -s {output[1]}
        """

#######################################################
# Section 2.6: Combine results for monotonic analysis
#######################################################

rule combine_all_monotonic_increases_exp1:
    input:
        expand("exp1_results/full_dap_{digestion_type}_analysis/increases.csv", digestion_type=["full_text", "minimizers", "dna_minimizers"])
    output:
        "exp1_final_results/combined_increases.csv"
    run:
        with open(output[0], "w") as out_fd:
            for digestion_type in ["full_text", "minimizers", "dna_minimizers"]:
                with open(f"exp1_results/full_dap_{digestion_type}_analysis/increases.csv") as in_fd:
                    lines = [x.strip() for x in in_fd.readlines()]
                    for line in lines:
                        out_fd.write(f"{digestion_type},{line}\n")


#######################################################
# Section 2.6: Combine index sizees into csv file
#######################################################

rule combine_all_index_sizes_into_csv_file_exp1:
    input:
        expand("exp1_index/{index_type}_{digestion}/output.fna.sdap",
                           index_type=["full_dap"],
                           digestion=["full_text", "minimizers", "dna_minimizers"]),
        expand("exp1_index/{index_type}_{digestion}/output.fna.taxcomp.sdap",
                           index_type=["tax_compressed"],
                           digestion=["full_text", "minimizers", "dna_minimizers"])
        
    output:
        "exp1_final_results/results.csv"
    run:
        file_extensions_full = [".bwt.cliffy", ".F.cliffy", ".sdap", ".edap", ".runcnt"]
        file_extensions_tax = [".bwt.cliffy", ".F.cliffy", ".taxcomp.sdap", ".taxcomp.edap", 
                               ".taxcomp.of.sdap", ".taxcomp.of.edap", ".taxcomp.ofptr.sdap", ".taxcomp.ofptr.edap",
                               ".runcnt", ".ftab"]
        
        with open(output[0], "w") as out_fd:
            for index_type in ["full_dap", "tax_compressed"]:
                for digestion in ["full_text", "minimizers", "dna_minimizers"]:
                    curr_prefix = f"exp1_index/{index_type}_{digestion}/output.fna"
                    if index_type == "full_dap":
                        file_extensions = file_extensions_full
                    else:
                        file_extensions = file_extensions_tax
                    
                    # iterate through files and sum up size
                    total_size = 0
                    for file_extension in file_extensions:
                        curr_file = f"{curr_prefix}{file_extension}"
                        if os.path.exists(curr_file):
                            total_size += os.path.getsize(curr_file)
                        else:
                            raise FileNotFoundError(f"File {curr_file} does not exist")
                    
                    # write to output file
                    out_fd.write(f"{index_type},{digestion},{total_size},{round(total_size/1073741824, 4)}\n")

        



##############################################
# Section 2.6: Overall rule of experiment 1
##############################################

rule run_exp1:
    input:
        expand("exp1_results/full_dap_{digestion_type}_analysis/increases.csv", digestion_type=["full_text", "minimizers", "dna_minimizers"]),
        expand("exp1_results/full_dap_{digestion_type}_analysis/sig_vals.csv", digestion_type=["full_text", "minimizers", "dna_minimizers"]),
        expand("exp1_index/tax_compressed_{digestion_type}/output.fna.taxcomp.sdap", digestion_type=["full_text", "minimizers", "dna_minimizers"]),
        expand("exp1_index/tax_compressed_{digestion_type}/output.fna.taxcomp.edap", digestion_type=["full_text", "minimizers", "dna_minimizers"]),
        expand("exp1_index/full_dap_{digestion_type}/output.fna.sdap", digestion_type=["full_text", "minimizers", "dna_minimizers"]),
        expand("exp1_index/full_dap_{digestion_type}/output.fna.edap", digestion_type=["full_text", "minimizers", "dna_minimizers"]),
