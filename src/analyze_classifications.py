
# Name: analyze_classifications.py
# Description: Compare the classifications produced by 
#              by cliffy and kraken to see how they differ
#              and explain when they agree and when they don't
# Date: March 27th, 2024

import argparse
import os



#################################################
# Helper methods
#################################################

def create_truthset_dict(truthset_file):
    """ parse the two-column file into dictionary mapping read name to traversal """
    genus_id_to_taxa = {}
    with open(truthset_file, "r") as in_fd:
        for line in in_fd:
            line_split = [x.strip() for x in line.split()]

            # parse out the traversal and genus num
            last_pos_of_traversal = max([i for i in range(len(line_split)) if ";" in line_split[i]])  
            traversal = ' '.join(line_split[0:last_pos_of_traversal+1])
            curr_genus_num = int(line_split[last_pos_of_traversal+1])

            genus_id_to_taxa[curr_genus_num] = traversal
    return genus_id_to_taxa

def parse_silva_taxonomy(silva_taxonomy_file):
    taxa_to_silva_id = {}; silva_id_to_taxa = {}
    with open(silva_taxonomy_file, "r") as in_fd:
        for line in in_fd:
            last_pos_of_traversal, traversal = get_traversal_from_full_line(line)
            
            line_split = line.split()
            node_id = int(line_split[last_pos_of_traversal+1])
            rank = line_split[last_pos_of_traversal+2]

            taxa_to_silva_id[traversal] = (node_id, rank)
            silva_id_to_taxa[node_id] = (traversal, rank)
    return taxa_to_silva_id, silva_id_to_taxa

def get_traversal_from_full_line(line):
    """
    Get the traversal part of line, the presence of spaces in
    in traversal string makes it slightly less than trivial

    line: string with entire lines contents from silva taxonomy file
    """
    line_split = line.split()
    last_pos_of_traversal = max([i for i in range(len(line_split)) if ";" in line_split[i]])
    traversal = ' '.join(line_split[0:last_pos_of_traversal+1])
    
    assert (len(line_split) - last_pos_of_traversal) <= 4   
    return last_pos_of_traversal, traversal

def parse_doc_to_trav(doc_to_trav_file):
    full_traversal_to_id = {}
    with open(doc_to_trav_file, "r") as in_fd:
        for line in in_fd:
            line_split = [x.strip() for x in line.split()]
            curr_traversal = " ".join(line_split[1:])
            full_traversal_to_id[curr_traversal] = int(line_split[0])-1
    return full_traversal_to_id

def parse_taxa_to_seqlength(trav_to_seqlength_file):
    with open(trav_to_seqlength_file, "r") as in_fd:
        taxa_to_seqlength = {}
        for line in in_fd:
            last_pos_of_traversal, traversal = get_traversal_from_full_line(line)
            seqlength = int(line.split()[last_pos_of_traversal+1])

            assert (last_pos_of_traversal+1+1) == len(line.split())
            taxa_to_seqlength[traversal] = seqlength
        return taxa_to_seqlength

def generate_total_seqlength_for_any_taxa(genus_to_seqlength, taxa_to_silva_id):
    taxa_to_seqlength = {taxa: 0 for taxa in taxa_to_silva_id.keys()}
    for genus_taxa, genus_length in genus_to_seqlength.items():
        for taxa in taxa_to_seqlength.keys():
            if genus_taxa.startswith(taxa):
                taxa_to_seqlength[taxa] += genus_length
    return taxa_to_seqlength

def create_fp_reads_set(fp_names_file):
    fp_reads = set()
    with open(fp_names_file, "r") as fp_names:
        for line in fp_names:
            read_name = line.strip()
            fp_reads.add(read_name)
    return fp_reads

def load_kraken_results(kraken_results_file):
    read_name_to_kraken_results = {}
    with open(kraken_results_file, "r") as kraken_results:
        for line in kraken_results:
            line_split = line.split("\t")
            assert len(line_split) == 5, f"{line} had an unexpected number of columns"

            read_name = line_split[1]
            classification = int(line_split[2])
            results = line_split[4]

            read_name_to_kraken_results[read_name] = (classification, results)
    return read_name_to_kraken_results

def load_cliffy_results(mate1_path, mate2_path):
    read_name_to_cliffy_results = {}
    with open(mate1_path, "r") as mate1_fd, open(mate2_path, "r") as mate2_fd:
        mate1_lines = [x.strip() for x in  mate1_fd.readlines()]
        mate2_lines = [x.strip() for x in  mate2_fd.readlines()]
    
    header_m1 = ""; listing_m1 = ""; 
    header_m2 = ""; listing_m2 = ""; pos = 0

    for mate1_line, mate2_line in zip(mate1_lines, mate2_lines):
        if mate1_line.startswith(">") and mate2_line.startswith(">"):
            header_m1 = mate1_line[1:]; header_m2 = mate2_line[1:]
            pos += 1
        elif pos == 1:
            listing_m1 = mate1_line; listing_m2 = mate2_line 
            pos = 0

            # make sure the names are the same except the mate # at end
            assert header_m1[:-1] == header_m2[:-1]

            # map the read name for mate 1 to both listings
            read_name_to_cliffy_results[header_m1] = (listing_m1, listing_m2)
    return read_name_to_cliffy_results

def combine_results(fp_reads, read_name_to_kraken_results, read_name_to_cliffy_results, genus_id_to_taxa, taxa_to_silva_id, taxa_to_doc_id):
    full_read_set = []
    for read_name in fp_reads:
        # get the kraken/cliffy results
        kraken_classification, kraken_results = read_name_to_kraken_results[read_name]
        cliffy_mate1_listing, cliffy_mate2_listing = read_name_to_cliffy_results[read_name]

        # get the genus number and correct taxa
        genus_num = int(read_name.split("_")[1])
        assert genus_num in genus_id_to_taxa, f"Genus number {genus_num} not found in truthset"
        true_taxa = genus_id_to_taxa[genus_num]

        # look up true taxa in SILVA taxonomy
        assert true_taxa in taxa_to_silva_id, f"True taxa {true_taxa} not found in SILVA taxonomy"
        silva_correct_id, true_rank = taxa_to_silva_id[true_taxa]
        assert true_rank == "genus", f"True rank={true_rank} is not genus for silva node id #{silva_correct_id}"

        # look up the correct document in cliffy index
        assert true_taxa in taxa_to_doc_id, f"True taxa {true_taxa} not found in doc to trav mapping"
        cliffy_correct_id = taxa_to_doc_id[true_taxa]
    
        # create the object
        read_obj = {
            "read_name": read_name,
            "kraken_classification": kraken_classification,
            "kraken_results": kraken_results,
            "cliffy_mate1_listing": cliffy_mate1_listing,
            "cliffy_mate2_listing": cliffy_mate2_listing,
            "genus_num": genus_num,
            "true_taxa": true_taxa,
            "silva_correct_id": silva_correct_id,
            "cliffy_correct_id": cliffy_correct_id
        }

        # add the object to the full set
        full_read_set.append(read_obj)
    return full_read_set

def analyze_combined_results(full_read_set, silva_id_to_taxa, taxa_to_seqlength):
    # define overall lists to store results
    list_num_no_hits = []; list_num_hits_to_correct_genus = [];
    list_lengths_containing_doc = []; list_lengths_exclusive_to_doc = []

    list_clade_level_hits = [0, 0, 0, 0, 0, 0] # genus, family, order, class, phylum, kingdom
    list_clade_hit_seqlengths = [{"length_sum": 0, "num_seqs": 0} for _ in range(6)]

    # methods to analyze each method's results
    def process_kraken_result(input_str):
        output_pairs = input_str.split()
        counts = {}; not_found_count = 0

        for pair in output_pairs:
            if pair == "|:|":
                continue
            x, count = pair.split(":")

            count = int(count)
            node_id = 0 if x == "A" else int(x)   

            counts[node_id] = counts.get(node_id, 0) + count
        return counts

    def process_doc_listings_results(list_1, list_2, correct_doc):
        list_1_pairs = list_1.split(); list_2_pairs = list_2.split()
        assert len(list_1_pairs) % 2 == 0 and len(list_2_pairs) % 2 == 0
        
        complete_listing = list_1_pairs + list_2_pairs
        lengths_containing_doc = []
        lengths_exclusive_to_doc = []

        for i in range(0, len(complete_listing), 2):
            length = int(complete_listing[i].split(',')[1].strip(']')) - int(complete_listing[i].split(',')[0].strip('[')) + 1
            docs = list(map(int, complete_listing[i+1].strip('{}').split(',')))
            left_doc = docs[0]; right_doc = docs[-1]

            if correct_doc >= left_doc and correct_doc <= right_doc:
                lengths_containing_doc.append(length)
            if len(docs) == 1 and correct_doc == docs[0]:
                lengths_exclusive_to_doc.append(length)

        return lengths_containing_doc, lengths_exclusive_to_doc

    # go through each read and process it
    for read_obj in full_read_set:

        # kraken-related analysis
        kraken_analysis = process_kraken_result(read_obj["kraken_results"])
        num_no_hits = kraken_analysis.get(0, 0)
        num_hits_to_correct_genus = kraken_analysis.get(read_obj["silva_correct_id"], 0)

        ranks_hits = [silva_id_to_taxa[taxa_id][1] for taxa_id, count in kraken_analysis.items() if taxa_id not in [0, 1] for _ in range(count)]
        list_clade_level_hits[0] += ranks_hits.count("genus")
        list_clade_level_hits[1] += ranks_hits.count("family")
        list_clade_level_hits[2] += ranks_hits.count("order")
        list_clade_level_hits[3] += ranks_hits.count("class")
        list_clade_level_hits[4] += ranks_hits.count("phylum")
        list_clade_level_hits[5] += ranks_hits.count("domain") 

        for taxa_id, _ in kraken_analysis.items():
            if taxa_id not in [0, 1]:
                curr_taxa, curr_rank = silva_id_to_taxa[taxa_id]
                curr_taxa_length = taxa_to_seqlength[curr_taxa]

                if curr_rank == "genus":
                    list_clade_hit_seqlengths[0]["length_sum"] += curr_taxa_length
                    list_clade_hit_seqlengths[0]["num_seqs"] += 1
                elif curr_rank == "family":
                    list_clade_hit_seqlengths[1]["length_sum"] += curr_taxa_length
                    list_clade_hit_seqlengths[1]["num_seqs"] += 1
                elif curr_rank == "order":
                    list_clade_hit_seqlengths[2]["length_sum"] += curr_taxa_length
                    list_clade_hit_seqlengths[2]["num_seqs"] += 1
                elif curr_rank == "class":
                    list_clade_hit_seqlengths[3]["length_sum"] += curr_taxa_length
                    list_clade_hit_seqlengths[3]["num_seqs"] += 1
                elif curr_rank == "phylum":
                    list_clade_hit_seqlengths[4]["length_sum"] += curr_taxa_length
                    list_clade_hit_seqlengths[4]["num_seqs"] += 1
                elif curr_rank == "domain":
                    list_clade_hit_seqlengths[5]["length_sum"] += curr_taxa_length
                    list_clade_hit_seqlengths[5]["num_seqs"] += 1

        list_num_no_hits.append(num_no_hits)
        list_num_hits_to_correct_genus.append(num_hits_to_correct_genus)

        # debugging print statement ...
        # print(f"\nRead {read_obj['read_name']} had...\n\t{num_no_hits} no hits\n\t{num_hits_to_correct_genus}"
        #       f" hits to correct genus of {read_obj['silva_correct_id']}\n\t"
        #       f"it was classified to {read_obj['kraken_classification']}\n")

        lengths_containing_doc, lengths_exclusive_to_doc = process_doc_listings_results(read_obj["cliffy_mate1_listing"], read_obj["cliffy_mate2_listing"], read_obj["cliffy_correct_id"])
        list_lengths_containing_doc.extend(lengths_containing_doc)
        list_lengths_exclusive_to_doc.extend(lengths_exclusive_to_doc)

        # debugging print statement ...
        # print(read_obj["cliffy_mate1_listing"], read_obj["cliffy_mate2_listing"], read_obj["cliffy_correct_id"])
        # print(f"Lengths containing doc: {lengths_containing_doc}\nLengths exclusive to doc: {lengths_exclusive_to_doc}\n")

    return list_num_no_hits, list_num_hits_to_correct_genus, list_lengths_containing_doc, list_lengths_exclusive_to_doc, list_clade_level_hits, list_clade_hit_seqlengths

def write_results_to_disk(output_dir, output_file_prefix, list_num_no_hits, list_num_hits_to_correct_genus, list_lengths_containing_doc, list_lengths_exclusive_to_doc, list_clade_level_hits, list_clade_hit_seqlengths):
    print(f"Here are the max-values: {max(list_num_no_hits)}, {max(list_num_hits_to_correct_genus)}, {max(list_lengths_containing_doc)}, {max(list_lengths_exclusive_to_doc)}, {len(list_clade_level_hits)},{len(list_clade_hit_seqlengths)}")
    with open(output_dir + output_file_prefix + ".csv", "w") as out_fd:
        for i in range(400):
            num_hits_to_clade = round(list_clade_level_hits[i]/sum(list_clade_level_hits),3) if i < len(list_clade_level_hits) else 0
            avg_length = round(list_clade_hit_seqlengths[i]["length_sum"]/list_clade_hit_seqlengths[i]["num_seqs"], 3) if i < len(list_clade_hit_seqlengths) else 0
            out_fd.write(f"{output_file_prefix},{i},{list_num_no_hits.count(i)},{list_num_hits_to_correct_genus.count(i)},{list_lengths_containing_doc.count(i)},{list_lengths_exclusive_to_doc.count(i)},{num_hits_to_clade},{avg_length}\n")

#################################################
# Main method/Argument parsing
#################################################

def main(args):
    # create a dictionaries to help with translating between ids and taxa
    genus_id_to_taxa = create_truthset_dict(args.truthset_file)
    taxa_to_silva_id, silva_id_to_taxa = parse_silva_taxonomy(args.silva_taxonomy_file)
    taxa_to_doc_id = parse_doc_to_trav(args.doc_to_trav_file)
    genus_to_seqlength = parse_taxa_to_seqlength(args.trav_to_seqlength_file)
    taxa_to_seqlength = generate_total_seqlength_for_any_taxa(genus_to_seqlength, taxa_to_silva_id)

    # create a set for the reads of interest
    fp_reads = create_fp_reads_set(args.fp_names_file)

    # load the results from the kraken classification file
    read_name_to_kraken_results = load_kraken_results(args.kraken_results_file)

    # load the cliffy listings from the two listings file
    read_name_to_cliffy_results = load_cliffy_results(args.mate1_listings_file, args.mate2_listings_file)

    # assert that the read names are the same in all three files
    for read_name in fp_reads:
        assert read_name in read_name_to_kraken_results, f"Read name {read_name} not found in kraken results"
        assert read_name in read_name_to_cliffy_results, f"Read name {read_name} not found in cliffy results"
    
    # combine the results and create an object for each read
    full_read_set = combine_results(fp_reads, read_name_to_kraken_results, read_name_to_cliffy_results, genus_id_to_taxa, taxa_to_silva_id, taxa_to_doc_id)
    print(f"[log] loaded {len(full_read_set)} reads ...")

    # get the final analysis results and write them to disk
    list_num_no_hits, list_num_hits_to_correct_genus, list_lengths_containing_doc, list_lengths_exclusive_to_doc, list_clade_level_hits, list_clade_hit_seqlengths = analyze_combined_results(full_read_set, silva_id_to_taxa, taxa_to_seqlength)
    write_results_to_disk(args.output_dir, args.output_file_prefix, list_num_no_hits, list_num_hits_to_correct_genus, list_lengths_containing_doc, list_lengths_exclusive_to_doc, list_clade_level_hits, list_clade_hit_seqlengths)
    print(f"[log] wrote results to {args.output_dir + args.output_file_prefix + '.csv'}")

def parse_arguments():
    parser = argparse.ArgumentParser(description="Compare the classifications produced by cliffy and kraken")

    parser.add_argument("--fp-names", dest="fp_names_file", type=str, required=True, help="*.txt file with a read name on each line that represent FPs by Kraken but TPs for Cliffy")
    parser.add_argument("--mate1-listings", dest="mate1_listings_file", type=str, required=True, help="listings for mate 1 produced by cliffy")
    parser.add_argument("--mate2-listings", dest="mate2_listings_file", type=str, required=True, help="listings for mate 2 produced by cliffy")
    parser.add_argument("--kraken-results", dest="kraken_results_file", type=str, required=True, help="classification file from kraken")
    parser.add_argument("--truthset", dest="truthset_file", type=str, required=True, help="truthset file: maps genus ids to their taxa")
    parser.add_argument("--doc-to-trav", dest="doc_to_trav_file", type=str, required=True, help="mapping file from document ids to taxa traversal")
    parser.add_argument("--silva-taxonomy", dest="silva_taxonomy_file", type=str, required=True, help="taxonomy file from SILVA")
    parser.add_argument("--output-dir", dest="output_dir", type=str, required=True, help="directory to write output files to")
    parser.add_argument("--output-file-prefix", dest="output_file_prefix", type=str, required=True, help="prefix for output files")
    parser.add_argument("--trav-to-seqlength", dest="trav_to_seqlength_file", type=str, required=True, help="taxa to sequence length")

    args = parser.parse_args()
    return args

def check_arguments(args):
    # check the file arguments
    for file_path in [args.fp_names_file, args.mate1_listings_file, args.mate2_listings_file, args.kraken_results_file, args.truthset_file, args.doc_to_trav_file, args.silva_taxonomy_file, args.trav_to_seqlength_file]:
        if not os.path.exists(file_path):
            raise ValueError(f"File {file_path} does not exist")
    
    # check the output directory
    if not os.path.exists(args.output_dir):
        raise ValueError(f"Directory {args.output_dir} does not exist")
    
    # check if output dir has a trailing slash
    if args.output_dir[-1] != "/":
        args.output_dir += "/"

if __name__ == "__main__":
    args = parse_arguments()
    check_arguments(args)

    main(args)
    