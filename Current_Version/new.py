from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import AlignIO, Phylo
from Bio.Phylo import Consensus, TreeConstruction
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from matplotlib import pyplot as plt
import io

# Set this to the absolute path of your clustalo.exe
CLUSTAL_PATH = "C:/Users/price/PycharmProjects/CS_123A_Project/Clustal_Windows/clustalo.exe"
TREES_MADE_PER_CONSENSUS = 10


# Perform MSA using Clustal Omega
def perform_msa(input_file, output_file):
    clustalomega_cline = ClustalOmegaCommandline(infile=input_file, outfile=output_file, verbose=True, auto=True,
                                                 cmd=CLUSTAL_PATH)
    stdout, stderr = clustalomega_cline()


# Retrieve alignment from given file
def create_alignment(aligned_sequences_file):
    alignment = AlignIO.read(aligned_sequences_file, "fasta")
    return alignment


# Calculate Distance Matrix from Aligned Sequences
def calculate_distance_matrix(alignment):
    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(alignment)
    return distance_matrix


# Construct Phylogenetic Tree
def construct_tree(alignment, distance_matrix, algorithm='upgma'):
    calculator = DistanceCalculator('identity')
    if algorithm.lower() == 'upgma':
        x = TreeConstruction.DistanceTreeConstructor(distance_calculator=calculator, method='upgma')
        tree_collection = Consensus.bootstrap_trees(alignment, TREES_MADE_PER_CONSENSUS, x)
        tree = Consensus.majority_consensus(tree_collection)
        clade_nums = []
        for clade in tree.find_clades():
            if clade.confidence is not None:
                clade_nums.append(clade.confidence)

    elif algorithm.lower() == 'nj':
        x = TreeConstruction.DistanceTreeConstructor(distance_calculator=calculator, method='nj')
        tree_collection = Consensus.bootstrap_trees(alignment, TREES_MADE_PER_CONSENSUS, x)
        tree = Consensus.majority_consensus(tree_collection)
        clade_nums = []
        for clade in tree.find_clades():
            if clade.confidence is not None:
                clade_nums.append(clade.confidence)
    else:
        raise ValueError("Unsupported algorithm. Please choose either 'upgma' or 'nj'.")
    return tree


# Replace 'input_sequences.fasta' with your actual input file containing the sequences
input_sequences_file = 'sampleseq.txt'
aligned_sequences_file = 'outseq7.txt'

# Step 1: Perform MSA
# perform_msa(input_sequences_file, aligned_sequences_file)

# Step 2: Calculate the distance matrix
alignment = create_alignment('outseq7.txt')
distance_matrix = calculate_distance_matrix(alignment)
print(distance_matrix, "\n")

# Step 3: Construct the UPGMA tree
upgma_tree = construct_tree(alignment, distance_matrix, 'upgma')

# Step 4: Construct the NJ tree
nj_tree = construct_tree(alignment, distance_matrix, 'nj')

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