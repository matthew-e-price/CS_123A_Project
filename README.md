# CS_123A_Project Instructions

1) Load the code into your Python IDE of choice (I used PyCharm)
    - Unfortunately, we had issues making an executable file due to issues with Biopython
2) Install all libraries and dependencies
    - matplotlib
    - Bio
3) Fix the UPGMA bug by copying the code in TreeConstruction(FIXED).py into "Bio/Phylo/TreeConstruction"
    - In PyCharm, it's located at "venv/Lib/site-packages/Bio/Phylo/TreeConstruction"
    - You just need to copy and paste all the code
    - Everything will work if you don't do this, but UPGMA trees won't be ultrametric and correct
4) In phylo_main.py, provide the absolute path to ClustalO/clustalo.exe in the variable CLUSTAL_PATH
    - Make sure all instances of "\" are replaced with "/"
    - For example, "C:/Users/price/PycharmProjects/CS_123A_Project/ClustalO/clustalo.exe"
    - Sorry for using an absolute path, we couldn't find any other way to make the command line work properly
5) Put your text file in the Sequences directory
    - Must be ".txt"
    - Must be Fasta format
    - If any of the sequence names have spaces, then there could be issues
      - For example, "Black Bear" and "Black Cat" would result in an error, since the name for both would be "Black"
      - A solution would be to replace each space with "_", so "Black_Bear"
6) Put the name of your input file in the variable INPUT_FILE
    - For example, "FOX03.txt"
8) Set the number of bootstrap trees you want to be generated for each actual tree created in the variable BOOT_TREES
    - The more bootstrap trees generated, the longer the program will take
    - The default is 10, though you can use 100 if you want
9) Set if you wish to generate a new alignment in the variable USE_CLINE
    - True: Will generate a new alignment file
    - False: Will not generate a new alignment file
      - This will only work if there is already a generated alignment file for the input file, stored in the directory "Alignments"
      - For example, if you run the program with an input file, then you can set this to False and run it again since the alignment file has been created
10) Now that everything's set up, run phylo_main.py
    - If you are generating an alignment file, then this step could take several minutes
    - After the program finishes, you will be shown a table of statistics
    - The produced trees will be stored in the directory "Trees"
    - The table will be stored in the directory "Results"

## Potential matplotlib Issue with PyCharm
If you get an error when running phylo_main.py where there "is no attribute 'FigureCanvas'", then use the following fix
- Go to File | Settings | Tools | Python Scientific
- Disable the option "Show plots in tools window"
