CALDER
=======================

CALDER (Cancer Analysis of Longitudinal Data through Evolutionary Reconstruction) is an algorithm for inferring evolutionary phylogenies using multiple longitudinal bulk DNA sequencing samples from the same patient. CALDER improves upon previous methods by enforcing the evolutionary relationships that are expected between temporally ordered samples. 

Setup
------------------------
The setup process for CALDER requires the following steps:

### Download
Download CALDER.  The following command clones the current CALDER repository from GitHub:

    git clone https://github.com/raphael-group/calder.git

### Requirements
The following software is required for CALDER.

* Linux/Unix or Windows
* [Java 8 Runtime Environment (JRE)](https://www.oracle.com/technetwork/java/javase/downloads/jre8-downloads-2133155.html)
* [Gurobi optimizer version 8.0](http://www.gurobi.com/index)

### Testing
With the dependencies set up correctly, the following command will run CALDER on the provided test input:

    java -jar calder.jar -i CLL003_clustered.txt

This should take no more than a few seconds to run and the output should match sample_output.txt.

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
The main output of CALDER is a text file, CALDER_output.txt. By default, this file lists run details and, for each optimal solution, the Fhat, U, and tree, as well as the tmin and tmax values for each clone. See sample_output.txt for an example.

Stay tuned for additional scripts to visualize output.

### Command line options

    -i,--input <arg>       input file path
    -o,--output <arg>      output file path (default "CALDER_output.txt")
    -t,--threads <arg>     number of threads (default 1)
    -a,--alpha <arg>       confidence level alpha (default 0.9)
    -h,--threshold <arg>   detection threshold h (default 0.01)
    -n,--intervals         print confidence intervals (default false)
  
Additional information
----------------

### License
See `LICENSE` for license information.

### Citation
If you use CALDER in your work, please cite the following manuscript. (information to come)
