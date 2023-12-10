from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo
import io

from matplotlib import pyplot as plt

# Set this to the absolute path of your clustalo.exe
CLUSTAL_PATH = "C:/Users/price/PycharmProjects/CS_123A_Project/Clustal_Windows/clustalo.exe"


# Perform MSA using Clustal Omega
def perform_msa(input_file, output_file):
    clustalomega_cline = ClustalOmegaCommandline(infile=input_file, outfile=output_file, verbose=True, auto=True,
                                                 cmd=CLUSTAL_PATH)
    stdout, stderr = clustalomega_cline()


# Calculate Distance Matrix from Aligned Sequences
def calculate_distance_matrix(aligned_sequences_file):
    alignment = AlignIO.read(aligned_sequences_file, "fasta")
    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(alignment)
    return distance_matrix


# Construct Phylogenetic Tree
def construct_tree(distance_matrix, algorithm='upgma'):
    constructor = DistanceTreeConstructor()
    if algorithm.lower() == 'upgma':
        tree = constructor.upgma(distance_matrix)
    elif algorithm.lower() == 'nj':
        tree = constructor.nj(distance_matrix)
    else:
        raise ValueError("Unsupported algorithm. Please choose either 'upgma' or 'nj'.")
    return tree


# Output tree in Newick Format
def tree_to_newick(tree):
    newick_io = io.StringIO()
    Phylo.write(tree, newick_io, 'newick')
    return newick_io.getvalue().strip()


# Replace 'input_sequences.fasta' with your actual input file containing the sequences
input_sequences_file = 'sampleseq.txt'
aligned_sequences_file = 'outseq7.txt'

# Step 1: Perform MSA
perform_msa(input_sequences_file, aligned_sequences_file)

# Step 2: Calculate the distance matrix
distance_matrix = calculate_distance_matrix(aligned_sequences_file)
print(distance_matrix)

# Step 3: Construct the UPGMA tree
upgma_tree = construct_tree(distance_matrix, 'upgma')

# Step 4: Construct the NJ tree
nj_tree = construct_tree(distance_matrix, 'nj')

# Convert trees to Newick Format
upgma_newick = tree_to_newick(upgma_tree)
nj_newick = tree_to_newick(nj_tree)

# Output the results
print("UPGMA Tree in Newick Format:")
print(upgma_newick)
print("\nNJ Tree in Newick Format:")
print(nj_newick)

print("UPGMA Tree:")
Phylo.draw_ascii(upgma_tree)

print("\nNJ Tree:")
Phylo.draw_ascii(nj_tree)

fig1 = plt.figure(figsize=(13, 5), dpi=100)
axes1 = fig1.add_subplot(1, 1, 1)
Phylo.draw(upgma_tree, axes=axes1)

fig2 = plt.figure(figsize=(13, 5), dpi=100)
axes2 = fig2.add_subplot(1, 1, 1)
Phylo.draw(nj_tree, axes=axes2)
