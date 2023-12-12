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
                                                 cmd=CLUSTAL_PATH, force=True)
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
        for clade in new_tree.get_nonterminals():
            clade.name = ""
    elif algorithm.lower() == 'nj':
        tree_constructor = TreeConstruction.DistanceTreeConstructor(distance_calculator=calculator, method='nj')
        tree_collection = Consensus.bootstrap_trees(alignment, TREES_MADE_PER_CONSENSUS, tree_constructor)
        new_tree = Consensus.get_support(tree, tree_collection, len_trees=TREES_MADE_PER_CONSENSUS)
        for clade in new_tree.get_nonterminals():
            clade.name = ""
    else:
        raise ValueError("Unsupported algorithm. Please choose either 'upgma' or 'nj'.")
    return new_tree


def average_confidence(tree):
    clade_nums = []
    confidence = 0
    for clade in tree.get_nonterminals():
        if clade.confidence is not None:
            clade_nums.append(clade.confidence)
            confidence += clade.confidence
        else:
            clade_nums.append(0)
    return confidence / len(clade_nums)


def draw_tree(tree, file, algorithm="upgma"):
    if algorithm.lower() == 'upgma':
        fig = plt.figure(figsize=(15, 8), dpi=100)
        plt.rc("font", size=12)
        plt.rc("xtick", labelsize=10)
        plt.rc("ytick", labelsize=10)
        axes = fig.add_subplot(1, 1, 1)
        Phylo.draw(tree, axes=axes, show_confidence=False, do_show=False)
        filename = "Trees/" + file[:file.index(".txt")] + "_upgma"
        fig.savefig(filename)
        plt.close(fig)
    elif algorithm.lower() == 'nj':
        fig = plt.figure(figsize=(15, 8), dpi=100)
        plt.rc("font", size=12)
        plt.rc("xtick", labelsize=10)
        plt.rc("ytick", labelsize=10)
        axes = fig.add_subplot(1, 1, 1)
        Phylo.draw(tree, axes=axes, show_confidence=False, do_show=False)
        filename = "Trees/" + file[:file.index(".txt")] + "_nj"
        fig.savefig(filename)
        plt.close(fig)
    else:
        raise ValueError("Unsupported algorithm. Please choose either 'upgma' or 'nj'.")


def upgma_test_correct(distance_matrix, tree):
    for a in tree.get_terminals():
        for b in tree.get_terminals():
            for c in tree.get_terminals():
                if a != b and b != c and c != a:
                    distances = [round(tree.distance(a.name, b.name), 6),
                                 round(tree.distance(b.name, c.name), 6),
                                 round(tree.distance(c.name, a.name), 6)]
                    for i in range(3):
                        if distances[i] > max(distances[i-1], distances[i-2]):
                            return False
    return True


def make_trees(input, output):
    # Replace 'input_sequences.fasta' with your actual input file containing the sequences
    input_sequences_file = "Sequences/" + input
    aligned_sequences_file = "Alignments/" + output

    # Step 1: Perform MSA
    #perform_msa(input_sequences_file, aligned_sequences_file)

    # Step 2: Calculate the distance matrix
    alignment = create_alignment(aligned_sequences_file)
    distance_matrix = calculate_distance_matrix(alignment)

    # Step 3: Construct the UPGMA tree and compare it to bootstrap trees
    upgma_start = time.perf_counter()
    upgma_tree = construct_tree(distance_matrix, 'upgma')
    upgma_end = time.perf_counter()

    upgma_tree = compare_to_bootstrap(upgma_tree, alignment, "upgma")
    print(distance_matrix)
    print(upgma_test_correct(distance_matrix, upgma_tree))

    # Step 4: Construct the NJ tree
    nj_start = time.perf_counter()
    nj_tree = construct_tree(distance_matrix, 'nj')
    nj_end = time.perf_counter()

    nj_tree = compare_to_bootstrap(nj_tree, alignment, "nj")

    draw_tree(upgma_tree, input, "upgma")
    draw_tree(nj_tree, input, "nj")

    stats_dic = {
        "upgma_time": (upgma_end - upgma_start) * 1000,
        "upgma_confidence": average_confidence(upgma_tree),
        "nj_time": (nj_end - nj_start) * 1000,
        "nj_confidence": average_confidence(nj_tree)
    }
    return stats_dic

