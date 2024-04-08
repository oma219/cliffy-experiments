##################################################
# Name: exp_3.smk
# Description: Contains the workflow and methods
#              needed for experiment 2.
#
#              This experiment queries the
#              rRNA read dataset using both Kraken
#              and Cliffy.
#
#              Measure the run time and accuracy
#              for both approaches.
#
# Date: Mar 18th, 2024
##################################################

####################################################
# Section 1: Helper functions needed for these
#            experiment rules.
####################################################

####################################################
# Section 2: Rules needed for this experiment type
####################################################

######################################################
# Section 2.1: Run cliffy using index over full text
######################################################

rule run_cliffy_full_text_query:
    output:
        "exp3_cliffy_results/no_digestion/{dataset}/{region}/trial_1/mate_1.listings",
        "exp3_cliffy_results/no_digestion/{dataset}/{region}/trial_1/mate_1.time",
        "exp3_cliffy_results/no_digestion/{dataset}/{region}/trial_1/mate_1.log",
        "exp3_cliffy_results/no_digestion/{dataset}/{region}/trial_1/mate_2.listings",
        "exp3_cliffy_results/no_digestion/{dataset}/{region}/trial_1/mate_2.time",
        "exp3_cliffy_results/no_digestion/{dataset}/{region}/trial_1/mate_2.log",
        "exp3_cliffy_results/no_digestion/{dataset}/{region}/trial_2/mate_1.listings",
        "exp3_cliffy_results/no_digestion/{dataset}/{region}/trial_2/mate_1.time",
        "exp3_cliffy_results/no_digestion/{dataset}/{region}/trial_2/mate_1.log",
        "exp3_cliffy_results/no_digestion/{dataset}/{region}/trial_2/mate_2.listings",
        "exp3_cliffy_results/no_digestion/{dataset}/{region}/trial_2/mate_2.time",
        "exp3_cliffy_results/no_digestion/{dataset}/{region}/trial_2/mate_2.log"
    shell:
        """
        # Run mate 1 and mate 2 in the first trial
        pfp_doc64 run --ref exp3_cliffy_indexes/tax_compressed_full_text/output \
				      --pattern exp3_read_files/{wildcards.dataset}/{wildcards.region}_mate_1.fq \
                      --output exp3_cliffy_results/no_digestion/{wildcards.dataset}/{wildcards.region}/trial_1/mate_1 \
				      --taxcomp \
                      --ftab \
				      --num-col 7 2> {output[2]}
        grep 'querying the patterns' {output[2]} | awk '{{ print substr($7, 2) }}' > {output[1]}

        pfp_doc64 run --ref exp3_cliffy_indexes/tax_compressed_full_text/output \
				      --pattern exp3_read_files/{wildcards.dataset}/{wildcards.region}_mate_2.fq \
                      --output exp3_cliffy_results/no_digestion/{wildcards.dataset}/{wildcards.region}/trial_1/mate_2 \
				      --taxcomp \
                      --ftab \
				      --num-col 7 2> {output[5]}
        grep 'querying the patterns' {output[5]} | awk '{{ print substr($7, 2) }}' > {output[4]}

        # Run mate 1 and mate 2 in the second trial
        pfp_doc64 run --ref exp3_cliffy_indexes/tax_compressed_full_text/output \
				      --pattern exp3_read_files/{wildcards.dataset}/{wildcards.region}_mate_1.fq \
                      --output exp3_cliffy_results/no_digestion/{wildcards.dataset}/{wildcards.region}/trial_2/mate_1 \
				      --taxcomp \
                      --ftab \
				      --num-col 7 2> {output[8]}
        grep 'querying the patterns' {output[8]} | awk '{{ print substr($7, 2) }}' > {output[7]}

        pfp_doc64 run --ref exp3_cliffy_indexes/tax_compressed_full_text/output \
				      --pattern exp3_read_files/{wildcards.dataset}/{wildcards.region}_mate_2.fq \
                      --output exp3_cliffy_results/no_digestion/{wildcards.dataset}/{wildcards.region}/trial_2/mate_2 \
				      --taxcomp \
                      --ftab \
				      --num-col 7 2> {output[11]}
        grep 'querying the patterns' {output[11]} | awk '{{ print substr($7, 2) }}' > {output[10]}
        """

########################################################################
# Section 2.2: Run cliffy using index over DNA minimizer digested text
########################################################################

rule run_cliffy_dna_minimizers_query:
    output:
        "exp3_cliffy_results/dna_minimizers/{dataset}/{region}/trial_1/mate_1.listings",
        "exp3_cliffy_results/dna_minimizers/{dataset}/{region}/trial_1/mate_1.time",
        "exp3_cliffy_results/dna_minimizers/{dataset}/{region}/trial_1/mate_1.log",
        "exp3_cliffy_results/dna_minimizers/{dataset}/{region}/trial_1/mate_2.listings",
        "exp3_cliffy_results/dna_minimizers/{dataset}/{region}/trial_1/mate_2.time",
        "exp3_cliffy_results/dna_minimizers/{dataset}/{region}/trial_1/mate_2.log",
        "exp3_cliffy_results/dna_minimizers/{dataset}/{region}/trial_2/mate_1.listings",
        "exp3_cliffy_results/dna_minimizers/{dataset}/{region}/trial_2/mate_1.time",
        "exp3_cliffy_results/dna_minimizers/{dataset}/{region}/trial_2/mate_1.log",
        "exp3_cliffy_results/dna_minimizers/{dataset}/{region}/trial_2/mate_2.listings",
        "exp3_cliffy_results/dna_minimizers/{dataset}/{region}/trial_2/mate_2.time",
        "exp3_cliffy_results/dna_minimizers/{dataset}/{region}/trial_2/mate_2.log"
    shell:
        """
        # Run mate 1 and mate 2 in the first trial
        pfp_doc64 run --ref exp3_cliffy_indexes/tax_compressed_dna_minimizers/output \
				      --pattern exp3_read_files/{wildcards.dataset}/{wildcards.region}_mate_1.fq \
                      --output exp3_cliffy_results/dna_minimizers/{wildcards.dataset}/{wildcards.region}/trial_1/mate_1 \
				      --taxcomp \
                      --ftab \
                      --dna-minimizers \
                      --small-window 4 \
                      --large-window 11 \
				      --num-col 7 2> {output[2]}
        grep 'querying the patterns' {output[2]} | awk '{{ print substr($7, 2) }}' > {output[1]}

        pfp_doc64 run --ref exp3_cliffy_indexes/tax_compressed_dna_minimizers/output \
				      --pattern exp3_read_files/{wildcards.dataset}/{wildcards.region}_mate_2.fq \
                      --output exp3_cliffy_results/dna_minimizers/{wildcards.dataset}/{wildcards.region}/trial_1/mate_2 \
				      --taxcomp \
                      --ftab \
                      --dna-minimizers \
                      --small-window 4 \
                      --large-window 11 \
				      --num-col 7 2> {output[5]}
        grep 'querying the patterns' {output[5]} | awk '{{ print substr($7, 2) }}' > {output[4]}

        # Run mate 1 and mate 2 in the second trial
        pfp_doc64 run --ref exp3_cliffy_indexes/tax_compressed_dna_minimizers/output \
				      --pattern exp3_read_files/{wildcards.dataset}/{wildcards.region}_mate_1.fq \
                      --output exp3_cliffy_results/dna_minimizers/{wildcards.dataset}/{wildcards.region}/trial_2/mate_1 \
				      --taxcomp \
                      --ftab \
                      --dna-minimizers \
                      --small-window 4 \
                      --large-window 11 \
				      --num-col 7 2> {output[8]}
        grep 'querying the patterns' {output[8]} | awk '{{ print substr($7, 2) }}' > {output[7]}

        pfp_doc64 run --ref exp3_cliffy_indexes/tax_compressed_dna_minimizers/output \
				      --pattern exp3_read_files/{wildcards.dataset}/{wildcards.region}_mate_2.fq \
                      --output exp3_cliffy_results/dna_minimizers/{wildcards.dataset}/{wildcards.region}/trial_2/mate_2 \
				      --taxcomp \
                      --ftab \
                      --dna-minimizers \
                      --small-window 4 \
                      --large-window 11 \
				      --num-col 7 2> {output[11]}
        grep 'querying the patterns' {output[11]} | awk '{{ print substr($7, 2) }}' > {output[10]}
        """

###################################################################
# Section 2.3: Run cliffy using index over minimizer digested text
###################################################################

rule run_cliffy_minimizers_query:
    output:
        "exp3_cliffy_results/minimizers/{dataset}/{region}/trial_1/mate_1.listings",
        "exp3_cliffy_results/minimizers/{dataset}/{region}/trial_1/mate_1.time",
        "exp3_cliffy_results/minimizers/{dataset}/{region}/trial_1/mate_1.log",
        "exp3_cliffy_results/minimizers/{dataset}/{region}/trial_1/mate_2.listings",
        "exp3_cliffy_results/minimizers/{dataset}/{region}/trial_1/mate_2.time",
        "exp3_cliffy_results/minimizers/{dataset}/{region}/trial_1/mate_2.log",
        "exp3_cliffy_results/minimizers/{dataset}/{region}/trial_2/mate_1.listings",
        "exp3_cliffy_results/minimizers/{dataset}/{region}/trial_2/mate_1.time",
        "exp3_cliffy_results/minimizers/{dataset}/{region}/trial_2/mate_1.log",
        "exp3_cliffy_results/minimizers/{dataset}/{region}/trial_2/mate_2.listings",
        "exp3_cliffy_results/minimizers/{dataset}/{region}/trial_2/mate_2.time",
        "exp3_cliffy_results/minimizers/{dataset}/{region}/trial_2/mate_2.log"
    shell:
        """
        # Run mate 1 and mate 2 in the first trial
        pfp_doc64 run --ref exp3_cliffy_indexes/tax_compressed_minimizers/output \
				      --pattern exp3_read_files/{wildcards.dataset}/{wildcards.region}_mate_1.fq \
                      --output exp3_cliffy_results/minimizers/{wildcards.dataset}/{wildcards.region}/trial_1/mate_1 \
				      --taxcomp \
                      --ftab \
                      --minimizers \
                      --small-window 4 \
                      --large-window 11 \
				      --num-col 7 2> {output[2]}
        grep 'querying the patterns' {output[2]} | awk '{{ print substr($7, 2) }}' > {output[1]}

        pfp_doc64 run --ref exp3_cliffy_indexes/tax_compressed_minimizers/output \
				      --pattern exp3_read_files/{wildcards.dataset}/{wildcards.region}_mate_2.fq \
                      --output exp3_cliffy_results/minimizers/{wildcards.dataset}/{wildcards.region}/trial_1/mate_2 \
				      --taxcomp \
                      --ftab \
                      --minimizers \
                      --small-window 4 \
                      --large-window 11 \
				      --num-col 7 2> {output[5]}
        grep 'querying the patterns' {output[5]} | awk '{{ print substr($7, 2) }}' > {output[4]}

        # Run mate 1 and mate 2 in the second trial
        pfp_doc64 run --ref exp3_cliffy_indexes/tax_compressed_minimizers/output \
				      --pattern exp3_read_files/{wildcards.dataset}/{wildcards.region}_mate_1.fq \
                      --output exp3_cliffy_results/minimizers/{wildcards.dataset}/{wildcards.region}/trial_2/mate_1 \
				      --taxcomp \
                      --ftab \
                      --minimizers \
                      --small-window 4 \
                      --large-window 11 \
				      --num-col 7 2> {output[8]}
        grep 'querying the patterns' {output[8]} | awk '{{ print substr($7, 2) }}' > {output[7]}

        pfp_doc64 run --ref exp3_cliffy_indexes/tax_compressed_minimizers/output \
				      --pattern exp3_read_files/{wildcards.dataset}/{wildcards.region}_mate_2.fq \
                      --output exp3_cliffy_results/minimizers/{wildcards.dataset}/{wildcards.region}/trial_2/mate_2 \
				      --taxcomp \
                      --ftab \
                      --minimizers \
                      --small-window 4 \
                      --large-window 11 \
				      --num-col 7 2> {output[11]}
        grep 'querying the patterns' {output[11]} | awk '{{ print substr($7, 2) }}' > {output[10]}
        """

rule run_exp3_queries:
    input:
        expand("exp3_cliffy_results/{index_type}/{dataset}/{region}/trial_{trial_num}/mate_{mate}.listings", 
                index_type=["no_digestion", "dna_minimizers", "minimizers"],
                dataset=["human_gut", "aquatic", "soil"],
                region=["V1_V2", "V3_V4", "V4_V4", "V4_V5"],
                trial_num=[1, 2],
                mate=[1, 2])

###################################################################
# Section 2.4: Run kraken2 on the all read files
###################################################################

rule run_default_kraken2_query:
    output:
        "exp3_kraken_results/{database}/{dataset}/{region}/trial_1/output.kreport2", 
        "exp3_kraken_results/{database}/{dataset}/{region}/trial_1/output.kraken2",
        "exp3_kraken_results/{database}/{dataset}/{region}/trial_1/output.kraken2.log", 
        "exp3_kraken_results/{database}/{dataset}/{region}/trial_1/output.bracken", 
        "exp3_kraken_results/{database}/{dataset}/{region}/trial_1/output.bracken.time", 
        "exp3_kraken_results/{database}/{dataset}/{region}/trial_1/output.total.time", 
        "exp3_kraken_results/{database}/{dataset}/{region}/trial_2/output.kreport2", 
        "exp3_kraken_results/{database}/{dataset}/{region}/trial_2/output.kraken2", 
        "exp3_kraken_results/{database}/{dataset}/{region}/trial_2/output.kraken2.log", 
        "exp3_kraken_results/{database}/{dataset}/{region}/trial_2/output.bracken", 
        "exp3_kraken_results/{database}/{dataset}/{region}/trial_2/output.bracken.time",
        "exp3_kraken_results/{database}/{dataset}/{region}/trial_2/output.total.time", 
    shell:
        """
        # Run trial 1
        kraken2 --db {wildcards.database} \
                --threads 1 \
                --report {output[0]} \
                --paired exp3_read_files/{wildcards.dataset}/{wildcards.region}_mate_1.fq  exp3_read_files/{wildcards.dataset}/{wildcards.region}_mate_2.fq > {output[1]} 2> {output[2]}
        
        {time_prog} {time_format} --output={output[4]} \
        bracken -d /vast/blangme2/oahmed/kraken2_dbs/{wildcards.database} \
				-r 250 \
				-l G  \
				-i {output[0]}\
				-o {output[3]}

        kraken_time=$(awk '/processed in/ {{sub(/s/, "", $7); print $7}}' {output[2]})
        bracken_time=$(awk '{{print $6}}' {output[4]})
        
        echo "$kraken_time + $bracken_time" > {output[5]}
        
        # Run trial 2
        kraken2 --db {wildcards.database}\
                --threads 1 \
                --report {output[6]} \
                --paired exp3_read_files/{wildcards.dataset}/{wildcards.region}_mate_1.fq  exp3_read_files/{wildcards.dataset}/{wildcards.region}_mate_2.fq > {output[7]} 2> {output[8]}
        
        {time_prog} {time_format} --output={output[10]} \
        bracken -d /vast/blangme2/oahmed/kraken2_dbs/{wildcards.database} \
				-r 250 \
				-l G  \
				-i {output[6]} \
				-o {output[9]}
        
        kraken_time=$(awk '/processed in/ {{sub(/s/, "", $7); print $7}}' {output[8]})
        bracken_time=$(awk '{{print $6}}' {output[10]})
        
        echo "$kraken_time + $bracken_time" > {output[11]}
        """

rule run_all_kraken2_queries_exp3:
    input:
        expand("exp3_kraken_results/{database}/{dataset}/{region}/trial_{trial_num}/output.kreport2",
                database=["silva_m31", "silva_m27", "silva_m23", "silva_m15", 
                          "silva_m31_nospaces", "silva_m27_nospaces", "silva_m23_nospaces", "silva_m15_nospaces"],
                dataset=["human_gut", "aquatic", "soil"],
                region=["V1_V2", "V3_V4", "V4_V4", "V4_V5"],
                trial_num=[1, 2])


###########################################################
# Section 2.5: Compute accuracy of cliffy output
###########################################################

rule analyze_results_from_docprof_approaches_exp3:
    input:
        "exp3_cliffy_results/{digestion}/{dataset}/{region}/trial_2/mate_1.listings",
        "exp3_cliffy_results/{digestion}/{dataset}/{region}/trial_2/mate_2.listings",
    output:
        "exp3_cliffy_analysis/{digestion}/{dataset}/{region}/output.classification_results.csv", 
        "exp3_cliffy_analysis/{digestion}/{dataset}/{region}/output.abundance_results.csv",
        "exp3_cliffy_analysis/{digestion}/{dataset}/{region}/TP_results.csv"
    shell:
        """
        python3 {repo_dir}/src/classify_silva_readset.py \
                        --mate1-listings {input[0]} \
                        --mate2-listings {input[1]} \
                        --doc-id-to-traversal {doc_to_trav_exp3} \
                        --silva-tax-ranks {silva_tax_rank_exp3} \
                        --readset-truthset exp3_read_files/{wildcards.dataset}/{wildcards.region}_seqtax.txt \
                        --output-dir exp3_cliffy_analysis/{wildcards.digestion}/{wildcards.dataset}/{wildcards.region}/ \
                        --trav-to-length {trav_to_seq_exp3} \
                        --sample-name output \
                        --classify-docprof 
        """

rule analyze_results_from_kraken2_approaches_exp3:
    input:
        "exp3_kraken_results/{database}/{dataset}/{region}/trial_2/output.kreport2", 
        "exp3_kraken_results/{database}/{dataset}/{region}/trial_2/output.kraken2", 
        "exp3_kraken_results/{database}/{dataset}/{region}/trial_2/output.bracken", 
    output:
        "exp3_kraken_analysis/{database}/{dataset}/{region}/output.classification_results.csv",
        "exp3_kraken_analysis/{database}/{dataset}/{region}/output.abundance_results.csv",
        "exp3_kraken_analysis/{database}/{dataset}/{region}/TP_FP_VP_results.txt",
    shell:
        """
        python3 {repo_dir}/src/classify_silva_readset.py \
                        --doc-id-to-traversal {doc_to_trav_exp3} \
                        --silva-tax-ranks {silva_tax_rank_exp3} \
                        --readset-truthset exp3_read_files/{wildcards.dataset}/{wildcards.region}_seqtax.txt \
                        --output-dir exp3_kraken_analysis/{wildcards.database}/{wildcards.dataset}/{wildcards.region}/ \
 				        --kraken2-read-class {input[1]} \
 				        --bracken-output {input[2]} \
                        --trav-to-length {trav_to_seq_exp3} \
                        --sample-name output \
                        --classify-kraken
        """

rule run_cliffy_analyses_on_all_datasets_exp3:
    input:
        expand("exp3_cliffy_analysis/{digestion}/{dataset}/{region}/output.classification_results.csv", 
                digestion=["no_digestion", "dna_minimizers", "minimizers"],
                dataset=["human_gut", "aquatic", "soil"],
                region=["V1_V2", "V3_V4", "V4_V4", "V4_V5"])

rule run_kraken_analyses_on_all_datasets_exp3:
    input:
        expand("exp3_kraken_analysis/{database}/{dataset}/{region}/output.classification_results.csv",
                database=["silva_m31", "silva_m27", "silva_m23", "silva_m15", 
                          "silva_m31_nospaces", "silva_m27_nospaces", "silva_m23_nospaces", "silva_m15_nospaces"],
                dataset=["human_gut", "aquatic", "soil"],
                region=["V1_V2", "V3_V4", "V4_V4", "V4_V5"])


###########################################################
# Section 2.6: Gather the time results for each method
###########################################################

rule generate_time_results_file_exp3:
    input:
        expand("exp3_cliffy_results/minimizers/{dataset}/{region}/trial_2/mate_{mate_num}.time",
                dataset=["human_gut", "aquatic", "soil"],
                region=["V1_V2", "V3_V4", "V4_V4", "V4_V5"],
                mate_num=[1, 2]),
        expand("exp3_kraken_results/silva_m31/{dataset}/{region}/trial_2/output.total.time",
                dataset=["human_gut", "aquatic", "soil"],
                region=["V1_V2", "V3_V4", "V4_V4", "V4_V5"])
    output:
        "exp3_final_results/time.csv"
    run:
        def get_time_from_cliffy_file(path):
            with open(path, "r") as in_fd:
                lines = [x.strip() for x in in_fd.readlines()]
                assert len(lines) == 1, f"Error: {path} has more than one line = {lines}"
                return float(lines[0])

        def get_time_from_kraken_file(path):
            with open(path, "r") as in_fd:
                lines = [x.strip() for x in in_fd.readlines()]
                assert len(lines) == 1, f"Error: {path} has more than one line = {lines}"
                time1, time2 = lines[0].split(" +" )
                return round(float(time1) + float(time2), 4)

        with open(output[0], "w") as out_fd:
            out_fd.write("dataset,region,method,digestion,time\n")
            for dataset in ["human_gut", "aquatic", "soil"]:
                for region in ["V1_V2", "V3_V4", "V4_V4", "V4_V5"]:

                    # summarize results for all types of cliffy 
                    for digestion_type in ["no_digestion", "dna_minimizers", "minimizers"]:
                        cliffy_files= [f"exp3_cliffy_results/{digestion_type}/{dataset}/{region}/trial_2/mate_1.time",
                                       f"exp3_cliffy_results/{digestion_type}/{dataset}/{region}/trial_2/mate_1.time"]

                        cliffy_time = get_time_from_cliffy_file(cliffy_files[0]) + get_time_from_cliffy_file(cliffy_files[1])
                        out_fd.write(f"{dataset},{region},cliffy,{digestion_type},{cliffy_time}\n")
                    
                    # extract the corresponding result for kraken
                    kraken_time = get_time_from_kraken_file(f"exp3_kraken_results/silva_m31/{dataset}/{region}/trial_2/output.total.time")
                    out_fd.write(f"{dataset},{region},kraken,na,{kraken_time}\n")

###########################################################
# Section 2.7: Gather the read classification results
###########################################################

rule generate_read_classification_results_exp3:
    input:
        expand("exp3_cliffy_analysis/{digestion}/{dataset}/{region}/output.classification_results.csv", 
                digestion=["no_digestion", "dna_minimizers", "minimizers"],
                dataset=["human_gut", "aquatic", "soil"],
                region=["V1_V2", "V3_V4", "V4_V4", "V4_V5"]),
        expand("exp3_kraken_analysis/{database}/{dataset}/{region}/output.classification_results.csv",
                database=["silva_m31", "silva_m27", "silva_m23", "silva_m15", 
                          "silva_m31_nospaces", "silva_m27_nospaces", "silva_m23_nospaces", "silva_m15_nospaces"],
                dataset=["human_gut", "aquatic", "soil"],
                region=["V1_V2", "V3_V4", "V4_V4", "V4_V5"])
    output:
        "exp3_final_results/accuracy.csv"
    run:
        def read_cliffy_classification_results(path):
            with open(path, "r") as in_fd:
                lines = [x.strip() for x in in_fd.readlines()]
                assert len(lines) == 12, f"Error: {path} has {len(lines)} lines."
                return lines

        def read_kraken_classification_results(path):
            with open(path, "r") as in_fd:
                lines = [x.strip() for x in in_fd.readlines()]
                assert len(lines) == 12, f"Error: {path} has {len(lines)} lines."
                return lines

        with open(output[0], "w") as out_fd:
            out_fd.write("dataset,region,digestion,method,query_approach,clade,sensitivity\n")
            for dataset in ["human_gut", "aquatic", "soil"]:
                for region in ["V1_V2", "V3_V4", "V4_V4", "V4_V5"]:

                    # copy all results from the individual files
                    for digestion_type in ["no_digestion", "dna_minimizers", "minimizers"]:
                        cliffy_files= [f"exp3_cliffy_analysis/{digestion_type}/{dataset}/{region}/output.classification_results.csv"]
                        cliffy_results = read_cliffy_classification_results(cliffy_files[0])

                        out_fd.writelines(f"{dataset},{region},{digestion_type},{s}\n" for s in cliffy_results)

                    # extract the corresponding result for kraken
                    for database_type in ["silva_m31", "silva_m27", "silva_m23", "silva_m15", "silva_m31_nospaces", "silva_m27_nospaces", "silva_m23_nospaces", "silva_m15_nospaces"]:
                        database_desc = database_type.split("silva_")[1]
                        kraken_results = read_kraken_classification_results(f"exp3_kraken_analysis/{database_type}/{dataset}/{region}/output.classification_results.csv")
                        
                        new_kraken_results = new_strings = [s.replace('kraken2_with_vp', 'kraken2_with_vp_' + database_desc).replace('kraken2_without_vp', 'kraken2_without_vp_' + database_desc) for s in kraken_results]
                        out_fd.writelines(f"{dataset},{region},na,{s}\n" for s in new_kraken_results)

###########################################################
# Section 2.7: Process classifications in order to compare
#              classifications
###########################################################

rule process_per_read_classifications_exp3:
    input:
        "exp3_cliffy_analysis/no_digestion/{dataset}/{region}/TP_results.csv",
        "exp3_kraken_analysis/{database}/{dataset}/{region}/TP_FP_VP_results.txt"
    output:
        "exp3_compare_acc/{database}/{dataset}/{region}/FP_Kraken.txt",
        "exp3_compare_acc/{database}/{dataset}/{region}/TP_Kraken.txt",
        "exp3_compare_acc/{database}/{dataset}/{region}/FP_in_kraken_TP_in_cliffy.txt",
        "exp3_compare_acc/{database}/{dataset}/{region}/TP_in_kraken_TP_in_cliffy.txt"
    shell:
        """
        awk -F, '{{ if ($2 == "FP") {{print $1}} }}' {input[1]} > {output[0]}
        awk -F, '{{ if ($2 == "TP") {{print $1}} }}' {input[1]} > {output[1]}

        awk -F, '{{ print $1 }}' {input[0]} {output[0]} | sort | uniq -d > {output[2]}
        awk -F, '{{ print $1 }}' {input[0]} {output[1]} | sort | uniq -d > {output[3]}
        """

rule analyze_cliffy_and_kraken_classifications_comparisons_exp3:
    input:
        "exp3_compare_acc/{database}/{dataset}/{region}/FP_in_kraken_TP_in_cliffy.txt",
        "exp3_compare_acc/{database}/{dataset}/{region}/TP_in_kraken_TP_in_cliffy.txt",
        "exp3_cliffy_results/no_digestion/{dataset}/{region}/trial_2/mate_1.listings",
        "exp3_cliffy_results/no_digestion/{dataset}/{region}/trial_2/mate_2.listings",
        "exp3_kraken_results/{database}/{dataset}/{region}/trial_2/output.kraken2"
    output:
        "exp3_compare_acc/{database}/{dataset}/{region}/FP.csv",
        "exp3_compare_acc/{database}/{dataset}/{region}/TP.csv",
        "exp3_compare_acc/{database}/{dataset}/{region}/FP_TP.csv"
    shell:
        """
        python3 {repo_dir}/src/analyze_classifications.py \
                --fp-names {input[0]} \
                --mate1-listings {input[2]} \
                --mate2-listings {input[3]} \
                --kraken-results {input[4]} \
                --truthset  exp3_read_files/{wildcards.dataset}/{wildcards.region}_seqtax.txt\
                --doc-to-trav {doc_to_trav_exp3} \
                --silva-taxonomy {silva_tax_rank_exp3} \
                --output-dir exp3_compare_acc/{wildcards.database}/{wildcards.dataset}/{wildcards.region}/  \
                --output-file-prefix FP \
                --trav-to-seqlength {trav_to_seq_exp3}

        python3 {repo_dir}/src/analyze_classifications.py \
                --fp-names {input[1]} \
                --mate1-listings {input[2]} \
                --mate2-listings {input[3]} \
                --kraken-results {input[4]} \
                --truthset  exp3_read_files/{wildcards.dataset}/{wildcards.region}_seqtax.txt \
                --doc-to-trav {doc_to_trav_exp3} \
                --silva-taxonomy {silva_tax_rank_exp3} \
                --output-dir exp3_compare_acc/{wildcards.database}/{wildcards.dataset}/{wildcards.region}/ \
                --output-file-prefix TP \
                --trav-to-seqlength {trav_to_seq_exp3}
        
        cat {output[0]} {output[1]} > {output[2]}
        """

###########################################################
# Section 2.8: Combine the results from Section 2.7
###########################################################

rule generate_read_classification_comparison_results_exp3:
    input:
        expand("exp3_compare_acc/{database}/{dataset}/{region}/FP_TP.csv", 
                database=["silva_m31", "silva_m27", "silva_m23", "silva_m15", 
                          "silva_m31_nospaces", "silva_m27_nospaces", "silva_m23_nospaces", "silva_m15_nospaces"],
                dataset=["aquatic", "soil"],
                region=["V1_V2", "V4_V4"]),
    output:
        "exp3_final_results/accuracy_debugging.csv"
    run:
        with open(output[0], "w") as out_fd:
            out_fd.write("dataset,region,type,x,no_hits,hits_to_correct_genus,lengths_contain_doc,lengths_exclusive_doc,percent_to_clade,avg_clade_size\n")
            for database in ["silva_m31", "silva_m27", "silva_m23", "silva_m15", 
                             "silva_m31_nospaces", "silva_m27_nospaces", "silva_m23_nospaces", "silva_m15_nospaces"]:
                for dataset in ["aquatic", "soil"]:
                    for region in ["V1_V2", "V4_V4"]:
                        with open(f"exp3_compare_acc/{database}/{dataset}/{region}/FP_TP.csv", "r") as in_fd:
                            lines = [x.strip() for x in in_fd.readlines()]
                            for line in lines:
                                out_fd.write(f"{database},{dataset},{region},{line}\n")


###########################################################
# Section 2.9: Combine the abundance results into one file
###########################################################

rule generate_final_abundance_file_exp3:
    input:
        expand("exp3_cliffy_analysis/minimizers/{dataset}/{region}/output.abundance_results.csv", 
                dataset=["human_gut", "aquatic", "soil"],
                region=["V1_V2", "V3_V4", "V4_V4", "V4_V5"]),
        expand("exp3_kraken_analysis/silva_m31/{dataset}/{region}/output.abundance_results.csv",
                dataset=["human_gut", "aquatic", "soil"],
                region=["V1_V2", "V3_V4", "V4_V4", "V4_V5"])
    output:
        "exp3_final_results/abundance.csv",
        "exp3_final_results/abundance_bc_distances.csv"
    run:
        with open(output[0], "w") as raw_data_fd, open(output[1], "w") as bc_data_fd:
            raw_data_fd.write("dataset,region,method,genus,abundance\n")
            bc_data_fd.write("dataset,region,method,bcdist\n")

            for dataset in ["human_gut", "aquatic", "soil"]:
                for region in ["V1_V2", "V3_V4", "V4_V4", "V4_V5"]:

                    # copy results from cliffy file
                    cliffy_file = f"exp3_cliffy_analysis/minimizers/{dataset}/{region}/output.abundance_results.csv"
                    with open(cliffy_file, "r") as in_fd:
                        lines = [x.strip() for x in in_fd.readlines()]
                        assert "Bray-Curtis" in lines[-1]
                        for line in lines[:-1]:
                            raw_data_fd.write(f"{dataset},{region},{line}\n")
                        bc_data_fd.write(f"{dataset},{region},cliffy,{lines[-1].split()[-1]}\n")
                    
                    # copy results from kraken file
                    kraken_file = f"exp3_kraken_analysis/silva_m31/{dataset}/{region}/output.abundance_results.csv"
                    with open(kraken_file, "r") as in_fd:
                        lines = [x.strip() for x in in_fd.readlines()]
                        assert "Bray-Curtis" in lines[-1]
                        for line in lines[:-1]:
                            if line.split(",")[0] != "truth":
                                raw_data_fd.write(f"{dataset},{region},{line}\n")
                        bc_data_fd.write(f"{dataset},{region},kraken2,{lines[-1].split()[-1]}\n")
