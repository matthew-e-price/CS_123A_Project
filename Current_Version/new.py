import time

from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import AlignIO, Phylo
from Bio.Phylo import Consensus, TreeConstruction
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from matplotlib import pyplot as plt
import io

# Set this to the absolute path of your clustalo.exe
CLUSTAL_PATH = "C:/Users/price/PycharmProjects/CS_123A_Project/Clustal_Windows/clustalo.exe"
TREES_MADE_PER_CONSENSUS = 5


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
def construct_tree(distance_matrix, algorithm='upgma'):
    constructor = DistanceTreeConstructor()
    if algorithm.lower() == 'upgma':
        tree = constructor.upgma(distance_matrix)
    elif algorithm.lower() == 'nj':
        tree = constructor.nj(distance_matrix)
    else:
        raise ValueError("Unsupported algorithm. Please choose either 'upgma' or 'nj'.")
    return tree


# Compare generated tree to bootstrap trees
def compare_to_bootstrap(tree, alignment, algorithm="upgma"):
    calculator = DistanceCalculator('identity')
    if algorithm.lower() == 'upgma':
        tree_constructor = TreeConstruction.DistanceTreeConstructor(distance_calculator=calculator, method='upgma')
        tree_collection = Consensus.bootstrap_trees(alignment, TREES_MADE_PER_CONSENSUS, tree_constructor)
        new_tree = Consensus.get_support(tree, tree_collection, len_trees=TREES_MADE_PER_CONSENSUS)
        clade_nums = []
        for clade in new_tree.find_clades():
            if clade.confidence is not None:
                clade_nums.append(clade.confidence)
    elif algorithm.lower() == 'nj':
        tree_constructor = TreeConstruction.DistanceTreeConstructor(distance_calculator=calculator, method='nj')
        tree_collection = Consensus.bootstrap_trees(alignment, TREES_MADE_PER_CONSENSUS, tree_constructor)
        new_tree = Consensus.get_support(tree, tree_collection, len_trees=TREES_MADE_PER_CONSENSUS)
        clade_nums = []
        for clade in new_tree.find_clades():
            if clade.confidence is not None:
                clade_nums.append(clade.confidence)
    else:
        raise ValueError("Unsupported algorithm. Please choose either 'upgma' or 'nj'.")
    return tree


def total_confidence(tree):
    clade_nums = []
    confidence = 0
    for clade in tree.find_clades():
        if clade.confidence is not None:
            clade_nums.append(clade.confidence)
            confidence += clade.confidence
    confidence = confidence / len(clade_nums)


# Replace 'input_sequences.fasta' with your actual input file containing the sequences
input_sequences_file = 'sampleseq.txt'
aligned_sequences_file = 'outseq7.txt'

# Step 1: Perform MSA
# perform_msa(input_sequences_file, aligned_sequences_file)

# Step 2: Calculate the distance matrix
alignment = create_alignment('outseq7.txt')
distance_matrix = calculate_distance_matrix(alignment)
# print(distance_matrix, "\n")

# Step 3: Construct the UPGMA tree and compare it to bootstrap trees
upgma_start = time.perf_counter()
upgma_tree = construct_tree(distance_matrix, 'upgma')
upgma_end = time.perf_counter()

upgma_tree = compare_to_bootstrap(upgma_tree, alignment, "upgma")

# Step 4: Construct the NJ tree
nj_start = time.perf_counter()
nj_tree = construct_tree(distance_matrix, 'nj')
nj_end = time.perf_counter()

nj_tree = compare_to_bootstrap(nj_tree, alignment, "nj")

print("UPGMA Tree:")
print((upgma_end - upgma_start) * 1000, "ms")
Phylo.draw_ascii(upgma_tree)

print("\nNJ Tree:")
print((nj_end - nj_start) * 1000, "ms")
Phylo.draw_ascii(nj_tree)

fig1 = plt.figure(figsize=(13, 5), dpi=100)
axes1 = fig1.add_subplot(1, 1, 1)
Phylo.draw(upgma_tree, axes=axes1)

fig2 = plt.figure(figsize=(13, 5), dpi=100)
axes2 = fig2.add_subplot(1, 1, 1)
Phylo.draw(nj_tree, axes=axes2)
