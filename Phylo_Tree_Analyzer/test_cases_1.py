from matplotlib import pyplot as plt

import phylo_wrk


def main():
    # List of dictionaries
    stats = []

    # Test 1
    result = phylo_wrk.make_trees("FOX03.txt")
    stats.append([result["upgma_time"], result["upgma_confidence"], result["nj_time"], result["nj_confidence"]])

    # Test 2
    result = phylo_wrk.make_trees("MLH1_gene_animals.txt")
    stats.append([result["upgma_time"], result["upgma_confidence"], result["nj_time"], result["nj_confidence"]])

    # Test 3
    result = phylo_wrk.make_trees("MLH1_gene_plants.txt")
    stats.append([result["upgma_time"], result["upgma_confidence"], result["nj_time"], result["nj_confidence"]])

    # Test 4
    result = phylo_wrk.make_trees("MLH1_gene_both.txt")
    stats.append([result["upgma_time"], result["upgma_confidence"], result["nj_time"], result["nj_confidence"]])

    # Test 5
    result = phylo_wrk.make_trees("MLH1_protein_animals.txt")
    stats.append([result["upgma_time"], result["upgma_confidence"], result["nj_time"], result["nj_confidence"]])

    # Test 6
    result = phylo_wrk.make_trees("MLH1_protein_plants.txt")
    stats.append([result["upgma_time"], result["upgma_confidence"], result["nj_time"], result["nj_confidence"]])

    # Test 7
    result = phylo_wrk.make_trees("MLH1_protein_both.txt")
    stats.append([result["upgma_time"], result["upgma_confidence"], result["nj_time"], result["nj_confidence"]])

    # Plot table of stats
    col_labels = ["UPGMA Time", "UPGMA Confidence", "NJ Time", "NJ Confidence"]
    row_labels = ["Test 1", "Test 2", "Test 3", "Test 4", "Test 5", "Test 6", "Test 7"]
    table, ax = plt.subplots(1, 1)
    ax.axis("tight")
    ax.axis("off")
    table = ax.table(cellText=stats,
                     rowLabels=row_labels,
                     colLabels=col_labels,
                     loc="center")
    plt.show()


if __name__ == "__main__":
    main()
