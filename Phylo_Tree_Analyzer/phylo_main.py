import phylo_wrk
import test_cases_1
from matplotlib import pyplot as plt


# Name of your txt file (MUST BE IN "Phylo_Tree_Analyzer/Sequences")
# EXAMPLE: "demo_sequence.txt"
INPUT_FILE = "FOX03.txt"


# 0 - Run with provided file
# 1 - Run using "test_cases_1"
# 2 - Run using "test_cases_2"
METHOD = 0


# Plot table of stats and generate phylogenetic trees
def main():
    if METHOD == 0:
        # List of dictionaries
        stats = []

        # Test 1
        result = phylo_wrk.make_trees(INPUT_FILE)
        stats.append([result["upgma_time"], result["upgma_confidence"], result["nj_time"], result["nj_confidence"]])

        col_labels = ["UPGMA Time", "UPGMA Confidence", "NJ Time", "NJ Confidence"]
        row_labels = ["Result"]
        table, ax = plt.subplots(1, 1)
        ax.axis("tight")
        ax.axis("off")
        table = ax.table(cellText=stats,
                         rowLabels=row_labels,
                         colLabels=col_labels,
                         loc="center")
        plt.show()
    elif METHOD == 1:
        test_cases_1.main()
    else:
        test_cases_1.main()


if __name__ == "__main__":
    main()
