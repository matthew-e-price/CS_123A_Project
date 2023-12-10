import time

import Bio.Phylo.Consensus
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from matplotlib import pyplot as plt


# Creates an alignment object using a clw file
def create_alignment(file_path):
    alignment = AlignIO.read(file_path, "clustal")
    return alignment


# Creates a distance matrix using the blosum62 table
def create_distance_matrix(alignment):
    calculator = DistanceCalculator("blosum62")
    distance_matrix = calculator.get_distance(alignment)
    return distance_matrix


# Creates a UPGMA tree given a distance matrix
def create_upgma_tree(distance_matrix, alignment):
    tree_constructor = DistanceTreeConstructor(method="upgma")
    upgma_tree = tree_constructor.upgma(distance_matrix)
    Bio.Phylo.Consensus.bootstrap_trees(alignment, 1000, tree_constructor)
    return upgma_tree


# Creates a NJ tree given a distance matrix
def create_nj_tree(distance_matrix):
    tree_constructor = DistanceTreeConstructor()
    nj_tree = tree_constructor.nj(distance_matrix)
    return nj_tree


# Creates trees through UPGMA and NJ and compares
# both the time taken to make them and the trees
# that are generated.
def main():
    alignment = create_alignment("sequence_example_1.clw")

    distance_matrix = create_distance_matrix(alignment)

    upgma_start = time.perf_counter()
    upgma_tree = create_upgma_tree(distance_matrix, alignment)
    upgma_end = time.perf_counter()

    nj_start = time.perf_counter()
    nj_tree = create_nj_tree(distance_matrix)
    nj_end = time.perf_counter()

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


if __name__ == "__main__":
    main()
