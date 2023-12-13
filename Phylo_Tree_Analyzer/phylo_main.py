# Class is written by Matthew Price

import phylo_wrk
import test_cases_1
from matplotlib import pyplot as plt


# Set this to the absolute path of your clustalo.exe
# MUST be set as absolute path, relative is not accepted for some reason
CLUSTAL_PATH = "ENTER ABS PATH HERE"

# Name of your txt file (MUST BE IN "Phylo_Tree_Analyzer/Sequences")
# EXAMPLE: "demo_sequence.txt"
INPUT_FILE = "FOX03.txt"

# Number of bootstrap trees to be generated, recommended 10
BOOT_TREES = 10

# Toggles the use of ClustalOmegaCommandLine
# Disable if you already have the proper alignment file in "/Alignments"
# Disabling will save a lot of time
USE_CLINE = True


# Plot table of stats and generate phylogenetic trees
def main():
    # List of dictionaries
    stats = []

    # Test 1
    result = phylo_wrk.make_trees(INPUT_FILE, BOOT_TREES, USE_CLINE, CLUSTAL_PATH)
    stats.append([result["upgma_time"], result["upgma_confidence"], result["nj_time"], result["nj_confidence"]])

    col_labels = ["UPGMA Time", "UPGMA Confidence", "NJ Time", "NJ Confidence"]
    row_labels = ["Result"]

    fig = plt.figure(figsize=(16, 12), dpi=100)
    ax = fig.add_subplot(1, 1, 1)

    ax.axis("tight")
    ax.axis("off")
    table = ax.table(cellText=stats,
                     rowLabels=row_labels,
                     colLabels=col_labels,
                     colWidths=[0.2 for x in col_labels],
                     loc="center")
    table.scale(1.0, 2.0)

    filename = "Results/" + INPUT_FILE[:INPUT_FILE.index(".txt")] + "_results"
    fig.savefig(filename)

    plt.show()


if __name__ == "__main__":
    main()
