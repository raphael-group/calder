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
* [Gurobi optimizer](http://www.gurobi.com/index) version 8.0, or any of the first 5 [ILP solvers listed here](http://javailp.sourceforge.net/)

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

### Running
The command to run CALDER is simply "java -jar calder.jar" followed by command line arguments. The option -i to designate the input file is required.

### Output
The main output of CALDER is a text file for each solution, <input filename>_tree<solution index>.txt. By default, this file lists the inferred clone proportion matrix U and tree (represented as a list of edges with 0-indexed ids corresponding to columns of U). See the sample_output folder for an example.

Stay tuned for additional utilities to visualize output.

### Command line options
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
    
Additional information
----------------

### License
See `LICENSE` for license information.

### Citation
If you use CALDER in your work, please cite the following preprint: https://www.biorxiv.org/content/10.1101/526814v1.
