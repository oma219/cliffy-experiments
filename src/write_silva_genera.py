#!/usr/bin/env python3

# Name: write_silva_genera.py
# Description: This script writes out each genera of SILVA separately
#              to a file. 
# Date: Mar 18th, 2024

import argparse
import os
import multiprocess as mp

#######################################################
# Section 1: Parse SILVA Reference and Taxonomy Files
#######################################################

def get_tax_lineage(line_split):
    """
    return: full taxonomic lineage from the SILVA taxonomy file

    The main special case we need to handle is when there are spaces
    in the taxonomic lineage, so we need to make sure we take care
    of that.

    e.g. Archaea;Crenarchaeota;Nitrososphaeria;Nitrosopumilales;Nitrosopumilaceae;Candidatus Nitrosopelagicus; 42941 genus 138
         should return ...
         Archaea;Crenarchaeota;Nitrososphaeria;Nitrosopumilales;Nitrosopumilaceae;Candidatus Nitrosopelagicus;
    """
    last_pos_of_traversal = max([i for i in range(len(line_split)) if ";" in line_split[i]])
    traversal = ' '.join(line_split[0:last_pos_of_traversal+1])
    return traversal

def get_tax_rank(line_split):
    """
    return: the taxonomic rank for a given line from SILVA file

    This is only slighly complicated by the fact that if you split
    every line by spaces, you could have a variable number. So you
    have go from the back BUT for some lines like Archaea; it looks like 
    this ...

    Archaea; 2 domain 

    Usually the is an additional field at the end that specifies the version.

    So the rule I will use is to scan for first non-number and return that 
    as the rank.
    """
    last_pos_of_traversal = max([i for i in range(len(line_split)) if ";" in line_split[i]])
    tax_rank = line_split[last_pos_of_traversal+2]
    return tax_rank 

def get_tax_id(line_split):
    """
    return: the taxonomic id for a given line from SILVA file
    """
    last_pos_of_traversal = max([i for i in range(len(line_split)) if ";" in line_split[i]])
    tax_id = line_split[last_pos_of_traversal+1]
    return int(tax_id)

def parse_taxonomy_file(taxonomy_file):
    """ parse the SILVA taxonomy (.txt) file """

    # read in the SILVA taxonomy file and build a dictionary mapping lineage to rank
    with open(taxonomy_file, 'r') as in_fd:
        all_lines = [line.strip().split() for line in in_fd]
        full_dict = {get_tax_lineage(x): get_tax_rank(x) for x in all_lines}
    assert len(all_lines) == len(full_dict), "issue with parsing the SILVA taxonomy file"

    # get a list of the genera from the list of lineages
    genera_list = [x for x in full_dict.keys() if full_dict[x] == "genus"]

    log_message("info", f"found {len(genera_list)} genera in the SILVA taxonomy file.")
    log_message("info", f"found {len(full_dict)} different lineages in the SILVA taxonomy file.")

    return genera_list, full_dict

def process_reference_file(ref_path, genera_to_seqs, full_lineage_dict):
    """ take reference sequences in SILVA file and parse into dictionary indexed by genera """

    def parse_taxonomy_lineage_to_genus(lineage):
        """
        Note: this method is exhaustive. It would be easy to just remove that last
              part of the lineage, but I want to be 100% that is the genus.

        e.g. Bacteria;Firmicutes;Bacilli;Bacillales;Planococcaceae;Planococcus;Planococcus citreus
             would return ...
             Bacteria;Firmicutes;Bacilli;Bacillales;Planococcaceae;Planococcus;
        """
        lineage_split = lineage.split(";")
        for i in range(len(lineage_split)-1, 0, -1):
            curr_lineage = ";".join(lineage_split[:i]) + ";"
            if full_lineage_dict[curr_lineage] == "genus":
                return curr_lineage

        # sequence is not specified to genus level
        return ""

    # go through each sequences in SILVA file
    with open(ref_path, "r") as ref_file:
        curr_seq = ""; curr_lineage = ""; curr_header = ""

        included_count = 0; not_included_count = 0
        for line in ref_file:
            if line[0] == ">" and curr_seq != "":
                # handle previous seq
                if curr_lineage in genera_to_seqs:
                    genera_to_seqs[curr_lineage].append((curr_header, curr_seq))
                    included_count += 1
                else:
                    not_included_count += 1

                # reset for next seq
                curr_seq = ""
                curr_lineage = parse_taxonomy_lineage_to_genus(" ".join(line.split()[1:])) # removes accession num at beginning
                curr_header = line.strip()
            elif line[0] == ">" and len(curr_seq) == 0:
                curr_lineage = parse_taxonomy_lineage_to_genus(" ".join(line.split()[1:])) # removes accession num at beginning
                curr_header = line.strip()
            else:
                curr_seq += line.strip()
        
        # process last seq
        if curr_lineage in genera_to_seqs:
            genera_to_seqs[curr_lineage].append((curr_header, curr_seq))
            included_count += 1
        else:
            not_included_count += 1

    log_message("results", f"{included_count} sequences were specified" 
                           f" to genus level, {not_included_count} were not.\n")
    return genera_to_seqs

#######################################################
# Section 2: Parse the SILVA taxonomy
#######################################################

class TreeNode(object):   
    def __init__(self, data, traversal="", rank=""):
        self.data = data
        self.children = []
        self.traversal = traversal
        self.rank = rank

    def add_child(self, obj):
        self.children.append(obj)
    
    def add_traversal(self, traversal):
        self.traversal = traversal
    
    def add_rank(self, rank):
        self.rank = rank

def get_silva_newick_tree(input_filepath):
    """ Loads the input in Newick Tree format """
    in_fd = open(input_filepath, "r")
    lines = [x.strip() for x in in_fd.readlines()]
    assert len(lines) == 1
    
    # Remove two extra characters at end
    silva_tree_str = lines[0][:-2]
    return silva_tree_str

def generate_node_id_to_traversal_dict(traversal_to_level_path):
    """ create a dictionary mapping node id to taxonomy traversal """
    node_map = {}
    with open(traversal_to_level_path, "r") as in_fd:
        for line in in_fd:
            line_split = line.split()

            traversal = get_tax_lineage(line_split)
            node_id = get_tax_id(line_split)
            tax_rank = get_tax_rank(line_split)
            
            node_map[node_id] = (traversal, tax_rank)
    return node_map

def build_newick_tree(curr_node, tree_str, node_map):
    """ builds the tree based on newick format """
    validate_subtree(tree_str)
    if "(" not in tree_str and ")" not in tree_str:
        children = tree_str.split(",")
        for child in children:
            assert int(child) in node_map
            curr_traversal, curr_rank = node_map[int(child)]
            curr_node.add_child(TreeNode(child, curr_traversal, curr_rank))
    else:
        interior_children = node_split(tree_str)
        for child in interior_children:
            if "(" not in child:
                assert int(child) in node_map
                curr_traversal, curr_rank = node_map[int(child)]
                curr_node.add_child(TreeNode(child, curr_traversal, curr_rank))
            else:
                curr_node.add_child(build_newick_tree(TreeNode("Interior"), strip_parentheses(child), node_map))
    return curr_node

def strip_parentheses(tree_str):
    return tree_str[1:-1]

def validate_subtree(tree_str):
    assert tree_str.count("(") == tree_str.count(")")

def node_split(input_str):
    """ 
    Split str based on comma but go around interior nodes so 
    ((One,Two),(Three,Four)),Five would return a list of 
    [((One,Two),(Three,Four)), Five]
    """
    out_list = []; pos = 0
    while pos < len(input_str):
        curr_str = ""
        # Case 1: interior node 
        if input_str[pos] == "(":
            curr_str += input_str[pos]
            status = 1
            pos += 1

            while status != 0:
                if input_str[pos] == "(":
                    status += 1
                elif input_str[pos] == ")":
                    status -= 1
                curr_str += input_str[pos]
                pos += 1
            out_list.append(curr_str)
        # Case 2: comma
        elif input_str[pos] == ",":
              pos += 1
        # Case 3: tip
        else:
            while pos < len(input_str) and input_str[pos] != ",":
                curr_str += input_str[pos]
                pos += 1
            out_list.append(curr_str)
    return out_list

def traversal_newick_tree(curr_node, tree_traversal):
    """ performs a pre-order traversal of tree """
    tree_traversal.append((curr_node.traversal, curr_node.rank))
    if len(curr_node.children) == 0:
        return tree_traversal
    else:
        for child in curr_node.children:
            tree_traversal = traversal_newick_tree(child, tree_traversal)
    return tree_traversal

#######################################################
# Section 3: Output data to files
#######################################################

def write_individual_genus_to_file(input_data):
    for genera_tuple in input_data:
        genus_num, traversal, genus_to_seq_dict, output_dir = genera_tuple
        with open(output_dir+f"doc_{genus_num+1}_seq.fa", "w") as out_fd:
            assert traversal in genus_to_seq_dict
            for header, seq in genus_to_seq_dict[traversal]:
                out_fd.write(f"{header}\n{seq}\n")
    return len(input_data)

def create_filelist(output_dir, num_genera_written):
    with open(output_dir+"filelist.txt", "w") as out_fd:
        for i in range(1, num_genera_written+1):
            out_fd.write(f"{output_dir}doc_{i}_seq.fa {i}\n")

def write_out_doc_to_traversal_file(output_dir, genera_written):
    with open(output_dir+"doc_to_traversal.txt", "w") as out_fd:
        for i,traversal,_,_ in genera_written:
            out_fd.write(f"{i+1} {traversal}\n")

#######################################################
# Main method & Argument parsing
#######################################################

def main(args):

    # Step 1: get a list of the genera and all sequences associated with them
    genera_list, full_lineage_dict = parse_taxonomy_file(args.tree_rank_path)
    genera_to_seqs = {genus: [] for genus in genera_list}
    genera_to_seqs = process_reference_file(args.input_file, genera_to_seqs, full_lineage_dict)

    # Step 2: build the SILVA taxonomy
    newick_tree = get_silva_newick_tree(args.tree_path)
    node_map = generate_node_id_to_traversal_dict(args.tree_rank_path)
    root = build_newick_tree(TreeNode("Root"), strip_parentheses(newick_tree), node_map)
    log_message("info", f"successfully built the SILVA taxonomy tree.\n")

    # Step 3: get the pre-order traversal of the tree
    traversal = traversal_newick_tree(root, [])
    genera_in_order = [tup[0] for tup in traversal if tup[1] == "genus"]
    log_message("info", f"computed the pre-order traversal of tree with {len(genera_in_order)} genera.\n")
    assert all([key in genera_in_order for key in genera_to_seqs.keys()]), "issue with parsing the SILVA taxonomy file"

    # Step 4: focus on genera with sequences (will be most of them, maybe only 1 or 2 without)
    genera_in_order_with_seq = [genus for genus in genera_in_order if len(genera_to_seqs[genus]) > 0]
    log_message("info", f"found {len(genera_in_order_with_seq)} genera with sequences.\n")

    # Step 5: determine how many genera to write and prepare data for parallel processing
    genera_to_write = args.num_genera if args.num_genera < len(genera_in_order_with_seq) else len(genera_in_order_with_seq)
    genera_write_list = [(i, genera_in_order_with_seq[i], genera_to_seqs, args.output_dir) for i in range(genera_to_write)]

    # Step 6: write out genera references in paralell
    batch_size = 250
    batches = [genera_write_list[i:i + batch_size] for i in range(0, len(genera_write_list), batch_size)]

    with mp.Pool() as pool:
        res = pool.map(write_individual_genus_to_file, batches)
        log_message("info", f"used {pool._processes} out of {mp.cpu_count()} processors to write out genera.")
    
    log_message("info", f"successfully wrote out {sum(res)} genera to separate files.\n")
    assert sum(res) == len(genera_write_list), "issue with writing out genera to separate files"

    # Step 7: Write out the metadata files
    create_filelist(args.output_dir, sum(res))
    log_message("info", "successfully wrote out filelist.txt")

    write_out_doc_to_traversal_file(args.output_dir, genera_write_list)
    log_message("info", "successfully wrote out doc_to_traversal.txt\n")


def parse_arguments():
    """ Defines the command-line argument parser, and return arguments """
    parser = argparse.ArgumentParser(description="This script allows users to separate a SILVA database at genus-level.")
    parser.add_argument("-i", "--input", dest="input_file", required=True, help="input SILVA database")
    parser.add_argument("-o", "--output-dir", dest="output_dir", required=True, help="output directory for files")
    parser.add_argument("--tree", dest="tree_path", default="", required=True, help="path to *.tre file for the SILVA database")
    parser.add_argument("--tree-rank", dest="tree_rank_path", default="", required=True, help="path to *.txt file for the SILVA database")
    parser.add_argument("-n", "--num", dest="num_genera", default=10, type=int, help="number of genera to write out to separate files")
    args = parser.parse_args()
    return args

def check_arguments(args):
    """ Validate the command-line arguments """

    # make sure input file and output directory are valid
    if not os.path.isfile(args.input_file):
        print("Error: the path for the provided input file is not valid.")
        exit(1)
    if not os.path.isdir(args.output_dir):
        print(f"Error: the output directory path is not valid = {args.output_dir}")
        exit(1)

    # modify output directory argument
    elif args.output_dir[-1] != '/':
        args.output_dir += "/"
    
    # make sure SILVA files look valid
    if not os.path.isfile(args.tree_path) or ".tre" not in args.tree_path:
        print("Error: *.tre file is not valid")
        exit(1)
    if not os.path.isfile(args.tree_rank_path) or ".txt" not in args.tree_rank_path:
        print("Error: *.txt file is not valid")
        exit(1)

def log_message(level="info", message=""):
    if level == "error":
        print(f"\033[31m[log::{level}]\033[0m {message}")
    elif level == "warning":
        print(f"\033[33m[log::{level}]\033[0m {message}")
    else:
        print(f"\033[32m[log::{level}]\033[0m {message}")

if __name__ == "__main__":
    args = parse_arguments()
    check_arguments(args)
    main(args)