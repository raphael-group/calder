CALDER
=======================

CALDER (Cancer Analysis of Longitudinal Data through Evolutionary Reconstruction) is an algorithm for inferring evolutionary phylogenies using multiple longitudinal bulk DNA sequencing samples from the same patient. CALDER improves upon previous methods by enforcing the evolutionary relationships that are expected between temporally ordered samples. 

Setup
------------------------
The setup process for CALDER requires the following steps:

### Download
The following command clones the current CALDER repository from GitHub:

    git clone https://github.com/raphael-group/calder.git

### Requirements
The following software is required for CALDER:

* Linux/Unix or Windows
* [Java 8 Runtime Environment (JRE)](https://www.oracle.com/technetwork/java/javase/downloads/jre8-downloads-2133155.html)
* [ILP solver](#ilp-solvers) (see below)
* Optional: [Absence-Aware Clustering](https://github.com/raphael-group/Absence-Aware-Clustering) to cluster mutations.

### Testing
With the dependencies set up correctly, the following command will run CALDER on the provided test input and write the results to a subdirectory called "output":

    java -jar calder.jar -i CLL003_clustered.txt -o output

This should take no more than a few seconds to run and the output should match the contents of the sample_output folder.

Use
----------------
CALDER has the following steps.

### Input
The input file is a tab-separated text file representing a matrix of read counts. Each row corresponds to a longitudinal sample, and alternating columns designate the reference reads and variant reads covering each mutation, respectively. For example, an instance with 3 samples and 4 mutations could be like so:

        a     a   b   b   c   c   d   d
    t1  700   300 0   0   0   0   0   0
    t2  700   300 800 200 900 100 900 100
    t3  600   400 800 200 900 100 900 100
    
For real datasets with a considerable number of mutations (more than 40), we recommend using [Absence-Aware Clustering](https://github.com/raphael-group/Absence-Aware-Clustering) to cluster mutations.

CALDER assumes that input mutations are in copy-neutral regions, i.e., that the number of reads with a mutation is proportional to the number of cells with that mutation. If you suspect this assumption does not hold for your data, consider excluding mutations that may be affected by CNA; alternatively, if you have copy number calls (e.g., from [HATCHet](https://github.com/raphael-group/hatchet)), you could correct the read counts to represent the true CCF.

### Running
The command to run CALDER is simply "java -jar calder.jar" followed by command line arguments. The option -i to designate the input file is required.

### Output
For each solution, CALDER produces 2 text files: a DOT file containing the inferred phylogenetic tree T, and a CSV file containing the inferred frequency matrix Fhat and the clone proportion matrix U. DOT files can be visualized using standard tools such as [`graphviz`](https://www.graphviz.org/) (see below for an example), and the matrices in CSV format can also be manipulated using standard tools -- (see soln_to_timescape.py for an example that does so using the [`pandas`](https://pandas.pydata.org/) library in Python).

### Visualizing output
To visualize a tree using [Graphviz](https://www.graphviz.org/) (after installing it), you can navigate to the output directory and run the following command:
```    
dot -Tpng CLL003_tree1.dot > CLL003_soln1.png
```
See the Graphviz documentation for more options.

We provide a script to support visualizing clone mixture proportions using the [Timescape R package](https://bioconductor.org/packages/release/bioc/html/timescape.html). This requires the following dependencies:
* Python 3
* Python packages: `networkx`, `pandas`, and `pydot`
* R >= 3.3
* R package: [Timescape](https://bioconductor.org/packages/release/bioc/html/timescape.html) (and its dependencies)

First, run the Python script to convert the solution DOT and CSV files to Timescape-formatted files (assuming that `python` refers to Python 3):
```
python soln_to_timescape.py outdir/CLL003_soln1.csv outdir/CLL003_tree1.dot CLL003
```
Then, run the following commands in R to generate the visualization:
```
library(timescape)
prev <- read.table("CLL003_prev.txt", header=TRUE)
edges <- read.table("CLL003_edges.txt", header=TRUE)
timescape(clonal_prev = prev, tree_edges = edges)
```
See the Timescape documentation for more options.

### Clustering mutations
We recommend clustering mutations by frequency before running CALDER - primarily because we generally expect to have multiple mutations distinguishing between any two clonal expansion events, and therefore between any two clones. We recommend using [Absence-Aware Clustering](https://github.com/raphael-group/Absence-Aware-Clustering), a clustering algorithm that pays particular attention to the distinction between mutation presence and absence. Python scripts are included to convert a CALDER input file to the format required by the clustering software, and to convert the clustering output back to CALDER input format.

Requirements:
* [Absence-Aware Clustering](https://github.com/raphael-group/Absence-Aware-Clustering) (see the link for dependencies as well as instructions for installation and usage)
* [Python 3](https://www.python.org/downloads/)

The following command converts CALDER-formatted input to clustering input (assuming that `python` refers to Python 3, otherwise use `python3` explicitly):
```
python calder_to_clustering.py calder_input.txt clustering_input.txt
```

Then, after running Absence-Aware Clustering, use the following command to apply the cluster assignments to the original data (where `clustering_assignments.txt` is the output file from the top level of the clustering output directory):
```
python apply_clustering.py clustering_input.txt cluster_assignments.txt calder_input_clustered.txt
```

### CALDER Command line options
    Required
    -i,--input <arg>       input file path
    -o,--output <arg>      output directory
    
    Additional options
    -a,--alpha <arg>       confidence level alpha (default 0.9)
    -h,--threshold <arg>   detection threshold h (default 0.01) 
    -O,--objective <arg>   objective function (l0 or l1)
    -N,--nonlongitudinal   do not enforce longitudinal constraints (for non-longitudinal data)
    -c,--printconf         print effective confidence level
    -g,--print-graph       print ancestry graph
    -n,--intervals         print confidence intervals (default false)
    -s,--solutions <arg>   maximum number of optimal solutions to return (default 1)
    -t,--time              track and output timing information
    -v,--solver <arg>      MILP solver back-end (default gurobi)
    
### ILP solvers
CALDER requires a specialized ILP solver. We recommend the [Gurobi optimizer](http://www.gurobi.com/index) (version 8.0 required), as it is fast, easy to install, and supported on all platforms (website includes instructions for obtaining a license, downloading, and installing). 

If for some reason you would prefer not to use Gurobi (e.g., if you are part of a non-academic entity and not interested in purchasing a license), we also support the [GLPK solver](https://www.gnu.org/software/glpk/) with the [GLPK for Java interface](http://glpk-java.sourceforge.net/), or the [lp-solve solver](https://sourceforge.net/projects/lpsolve/files/lpsolve/5.5.2.5/). We generally found GLPK to be faster and easier to use on all platforms. For more details, see the specific installation instructions for each solver. Installation tends to be easier on Linux systems than on Mac or Windows systems. Note that you will need to **specify the alternate solver** using the -v option.

Additional information
----------------
###
For assistance with running CALDER, interpreting the results, or other related questions, please email me (Matt Myers) at this address: [mm63@cs.princeton.edu](mailto:mm63@cs.princeton.edu)

### License
See `LICENSE` for license information.

### Citation
If you use CALDER in your work, please cite the following paper ([available here](https://www.cell.com/cell-systems/fulltext/S2405-4712(19)30191-7)):


Myers, M.A., Satas, G. and Raphael, B.J., 2019. CALDER: Inferring Phylogenetic Trees from Longitudinal Tumor Samples. *Cell Systems.*
