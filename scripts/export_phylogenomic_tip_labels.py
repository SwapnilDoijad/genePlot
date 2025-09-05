import argparse
from Bio import Phylo

def export_tip_labels(input_tree, output_file):
    # Load the tree from the input file
    tree = Phylo.read(input_tree, "newick")

    # Get all the terminal nodes (tips)
    tips = tree.get_terminals()

    # Extract labels sequentially
    tip_labels = [tip.name for tip in tips]

    # Write the labels to the output file
    with open(output_file, "w") as f:
        f.write("\n".join(tip_labels))

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Export phylogenomic tip labels sequentially")
    parser.add_argument("-i", "--input", help="Input tree file (Newick format)", required=True)
    parser.add_argument("-o", "--output", help="Output file for tip labels", required=True)
    args = parser.parse_args()

    # Call the function to export tip labels
    export_tip_labels(args.input, args.output)

#  source /home/groups/VEO/tools/biopython/myenv/bin/activate