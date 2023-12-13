from matplotlib import pyplot as plt
import phylo_wrk


# Set this to the absolute path of your clustalo.exe
# MUST be set as absolute path, relative is not accepted for some reason
CLUSTAL_PATH = "C:/Users/price/PycharmProjects/CS_123A_Project/ClustalO/clustalo.exe"

# Number of bootstrap trees to be generated, recommended 10
BOOT_TREES = 100

# Toggles the use of ClustalOmegaCommandLine
# Disable if you already have the proper alignment file in "/Alignments"
# Disabling will save a lot of time
USE_CLINE = False


def main():
    # List of dictionaries
    stats = []

    # Test 1
    result = phylo_wrk.make_trees("Convergent_Animal.txt", BOOT_TREES, USE_CLINE, CLUSTAL_PATH)
    stats.append([result["upgma_time"], result["upgma_confidence"], result["nj_time"], result["nj_confidence"]])

    # Test 2
    result = phylo_wrk.make_trees("Convergent_Plants.txt", BOOT_TREES, USE_CLINE, CLUSTAL_PATH)
    stats.append([result["upgma_time"], result["upgma_confidence"], result["nj_time"], result["nj_confidence"]])

    # Test 3
    result = phylo_wrk.make_trees("Divergent_Animal.txt", BOOT_TREES, USE_CLINE, CLUSTAL_PATH)
    stats.append([result["upgma_time"], result["upgma_confidence"], result["nj_time"], result["nj_confidence"]])

    # Test 4
    result = phylo_wrk.make_trees("Divergent_Plants.txt", BOOT_TREES, USE_CLINE, CLUSTAL_PATH)
    stats.append([result["upgma_time"], result["upgma_confidence"], result["nj_time"], result["nj_confidence"]])

    # Test 5
    result = phylo_wrk.make_trees("Homologous_Animals.txt", BOOT_TREES, USE_CLINE, CLUSTAL_PATH)
    stats.append([result["upgma_time"], result["upgma_confidence"], result["nj_time"], result["nj_confidence"]])

    # Test 6
    result = phylo_wrk.make_trees("Homologous_Plants.txt", BOOT_TREES, USE_CLINE, CLUSTAL_PATH)
    stats.append([result["upgma_time"], result["upgma_confidence"], result["nj_time"], result["nj_confidence"]])

    # Test 7
    result = phylo_wrk.make_trees("Non_Homologous_Animals.txt", BOOT_TREES, USE_CLINE, CLUSTAL_PATH)
    stats.append([result["upgma_time"], result["upgma_confidence"], result["nj_time"], result["nj_confidence"]])

    # Test 8
    result = phylo_wrk.make_trees("Non_Homologous_Plants.txt", BOOT_TREES, USE_CLINE, CLUSTAL_PATH)
    stats.append([result["upgma_time"], result["upgma_confidence"], result["nj_time"], result["nj_confidence"]])

    # Test 9
    result = phylo_wrk.make_trees("MLH1_Animal_Gene_Dataset.txt", BOOT_TREES, USE_CLINE, CLUSTAL_PATH)
    stats.append([result["upgma_time"], result["upgma_confidence"], result["nj_time"], result["nj_confidence"]])

    # Test 10
    result = phylo_wrk.make_trees("MLH1_Plant_Gene_Dataset.txt", BOOT_TREES, USE_CLINE, CLUSTAL_PATH)
    stats.append([result["upgma_time"], result["upgma_confidence"], result["nj_time"], result["nj_confidence"]])

    # Plot table of stats
    col_labels = ["UPGMA Time", "UPGMA Confidence", "NJ Time", "NJ Confidence"]
    row_labels = ["Convergent Animals", "Convergent Plants", "Divergent Animals", "Divergent Plants",
                  "Homologous Animals", "Homologous Plants", "Non-Homologous Animals", "Non-Homologous Plants",
                  "MLH1 Animals (All)", "MLH1 Plants (All)"]

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

    filename = "Results/" + "test_cases_2" + "_results"
    fig.savefig(filename)

    plt.show()


if __name__ == "__main__":
    main()
