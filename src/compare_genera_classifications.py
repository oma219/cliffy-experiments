
# Name: compare_genera_classifications.py
# Description: Compare the classifications across different genera
# Date: March 27th, 2024

import argparse
import os
import numpy as np

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

def process_input_file(input_file, genus_id_to_tp_fp):
    with open(input_file, "r") as in_fd:
        lines = [x.strip() for x in in_fd.readlines()]
        for line in lines:
            curr_genus_num = int(line.split(",")[0].split("_")[1])
            curr_read_result = line.split(",")[1]
            assert curr_genus_num in genus_id_to_tp_fp

            if curr_read_result == "TP":
                genus_id_to_tp_fp[curr_genus_num][0] += 1
            else:
                genus_id_to_tp_fp[curr_genus_num][1] += 1
    return genus_id_to_tp_fp


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

    # compute the number of TPs vs (FPs + VPs) for each genera id in readset
    genus_id_to_tp_fp = {x: [0, 0] for x in range(1, 101)}
    genus_id_to_tp_fp = process_input_file(args.input_file, genus_id_to_tp_fp)

    # sort this dictionary by genera that we get the highest % correct
    genus_id_to_tp_fp_sorted = dict(sorted(genus_id_to_tp_fp.items(), key=lambda item: item[1][0] / (item[1][0] + item[1][1]) if (item[1][0] + item[1][1]) != 0 else 0, reverse=True))
    x = []; y = []
    for key, value in genus_id_to_tp_fp_sorted.items():
        if sum(value) == 0:
            continue

        curr_genus_num = key
        curr_genus_traversal = genus_id_to_taxa[curr_genus_num]
        curr_phylum_taxa = ";".join(curr_genus_traversal.split(";")[:2]) + ";"
        curr_class_taxa = ";".join(curr_genus_traversal.split(";")[:3]) + ";"
        curr_order_taxa = ";".join(curr_genus_traversal.split(";")[:4]) + ";"
        curr_family_taxa = ";".join(curr_genus_traversal.split(";")[:5]) + ";"

        total_size = taxa_to_seqlength[curr_phylum_taxa] + taxa_to_seqlength[curr_class_taxa] + taxa_to_seqlength[curr_order_taxa] + taxa_to_seqlength[curr_family_taxa] + taxa_to_seqlength[curr_genus_traversal]

        print(key, value, round(value[0]/sum(value), 2), round(total_size/1000000, 2))
        x.append(value[0]/sum(value)); y.append(total_size/1000000)
        # print(key, value, curr_phylum_taxa, round(taxa_to_seqlength[curr_phylum_taxa]/1000000))
        # print(key, value, curr_genus_traversal)
    
    print(np.corrcoef(x, y))


def parse_arguments():
    parser = argparse.ArgumentParser(description="Compare the classifications produced by cliffy and kraken")

    parser.add_argument("--input", dest="input_file", type=str, required=True, help="*.txt file with a read name on each line that represent FPs by Kraken but TPs for Cliffy")
    parser.add_argument("--truthset", dest="truthset_file", type=str, required=True, help="truthset file: maps genus ids to their taxa")
    parser.add_argument("--doc-to-trav", dest="doc_to_trav_file", type=str, required=True, help="mapping file from document ids to taxa traversal")
    parser.add_argument("--silva-taxonomy", dest="silva_taxonomy_file", type=str, required=True, help="taxonomy file from SILVA")
    parser.add_argument("--trav-to-seqlength", dest="trav_to_seqlength_file", type=str, required=True, help="taxa to sequence length")

    args = parser.parse_args()
    return args

def check_arguments(args):
    # check the file arguments
    for file_path in [args.input_file, args.truthset_file, args.doc_to_trav_file, args.silva_taxonomy_file, args.trav_to_seqlength_file]:
        if not os.path.exists(file_path):
            raise ValueError(f"File {file_path} does not exist")
    
if __name__ == "__main__":
    args = parse_arguments()
    check_arguments(args)

    main(args)