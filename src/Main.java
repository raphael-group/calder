import org.apache.commons.cli.*;

import java.io.*;
import java.time.LocalTime;
import java.util.LinkedList;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class Main {
    public static String INFILE = null;
    public static String OUTDIR = "";

    public static String OUTFILE = "CALDER_output.txt";
    public static int MAX_NUM_OPTIMA_OUTPUT = 20;

    public static boolean PRINT_ANCESTRY_GRAPH = false;
    public static boolean PRINT_CONFIDENCE_INTERVALS = false;
    public static int THREADS = 3;
    public static double CONFIDENCE = .9;
    public static double MINIMUM_USAGE_H = .01;

    public static void main(String[] args) {
        Options options = new Options();
        Option input = new Option("i", "input", true, "input file path");
        input.setRequired(true);
        options.addOption(input);

        Option output = new Option("o", "output", true, "output file path");
        output.setRequired(false);
        options.addOption(output);

        Option alpha = new Option("a", "alpha", true, "confidence level alpha (default 0.9)");
        alpha.setRequired(false);
        options.addOption(alpha);

        Option threads = new Option("t", "threads", true, "number of threads (default 1)");
        threads.setRequired(false);
        options.addOption(threads);

        Option detectionThreshold = new Option("h", "threshold", true, "detection threshold h (default 0.01)");
        detectionThreshold.setRequired(false);
        options.addOption(detectionThreshold);

        Option printCI = new Option("n", "intervals", false, "print confidence intervals (default false)");
        printCI.setRequired(false);
        options.addOption(printCI);


        // Option printGraph = new Option("g", "print-graph", false, "print confidence intervals (default false)");
        //printGraph.setRequired(false);
        //options.addOption(printGraph);


        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd = null;

        try{
            cmd = parser.parse(options, args);
        } catch (ParseException e){
            System.out.println(e.getMessage());
            formatter.printHelp("calder", options);

            System.exit(1);
        }

        try {
            INFILE = cmd.getOptionValue("input");
            if (cmd.hasOption("output")){
                OUTFILE = cmd.getOptionValue("output");
            }
            if (cmd.hasOption("threads")) {
                THREADS = Integer.parseInt(cmd.getOptionValue("threads"));
            }
            if (cmd.hasOption("threshold")) {
                MINIMUM_USAGE_H = Double.parseDouble(cmd.getOptionValue("threshold"));
                if(MINIMUM_USAGE_H <= 0 || MINIMUM_USAGE_H >= 1){
                    System.out.println("Minimum usage h must be in (0, 1)");
                    System.exit(1);
                }
            }
            if (cmd.hasOption("alpha")) {
                CONFIDENCE = Double.parseDouble(cmd.getOptionValue("alpha"));
                if(CONFIDENCE <= 0 || CONFIDENCE >= 1){
                    System.out.println("Confidence level alpha must be in (0, 1)");
                    System.exit(1);
                }
            }
            if(cmd.hasOption("print-graph")){
                PRINT_ANCESTRY_GRAPH = true;
            }
            if(cmd.hasOption("intervals")){
                PRINT_CONFIDENCE_INTERVALS = true;
            }
        } catch(NumberFormatException e){
            System.out.println(e.getMessage());
            formatter.printHelp("calder", options);
            System.exit(1);
        }

        // Set up CALDER
        Calder.init();

        // Load data from file (TSV, basically AncesTree format)
        String fname = INFILE;
        Instance I = Instance.fromFile(fname);



        PrintStream origErr = System.err;
        PrintStream origOut = System.out;

        PrintStream logOut = null;
        try {
            logOut = new PrintStream(OUTDIR + OUTFILE);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        PrintStream dummy = new PrintStream(new OutputStream() {
            @Override
            public void write(int b) throws IOException {

            }});

        System.setOut(logOut);
        if(PRINT_CONFIDENCE_INTERVALS){
            System.out.println(I);
        }


        int totalTrees;
        int discardedTrees;
        LinkedList<Tree> maximalFTrees;
        int nFeasible, biggestTrees;
        ILPResultCollector coll;

        ExecutorService pool = Executors.newWorkStealingPool(THREADS);
        LinkedList<Tree> trees;
        int doneSize = I.nMuts;

        System.setOut(dummy);
        System.setErr(dummy);

        do {
            discardedTrees = 0;
            trees = Calder.solveLVAFFP(I, doneSize);
            totalTrees = trees.size();
            do {
                // Select only those maximal trees that account for a maximal amount of frequency in F
                maximalFTrees = Calder.findBestMaximalTrees(I, trees);
                biggestTrees = maximalFTrees.size();

                CountDownLatch latch = new CountDownLatch(maximalFTrees.size());
                LinkedList<SolveILPTask> tasks = new LinkedList<>();
                coll = new ILPResultCollector(I);

                SolveILPTask task;
                for (Tree T : maximalFTrees) {
                    task = new SolveILPTask(I, T, latch, coll);
                    tasks.add(task);
                    pool.submit(task);
                    //task.run();
                }

                try {
                    latch.await();
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }


                discardedTrees += biggestTrees;
                nFeasible = coll.getNumFeasible();
            } while (nFeasible == 0 && maximalFTrees.size() > 0);
            Calder.resetCap();
            doneSize = trees.get(0).vertices.size();
            System.out.println("Found " + nFeasible + " feasible trees of size " + doneSize + " ================================");
            doneSize = trees.get(0).vertices.size() - 1;
        } while(nFeasible == 0 && doneSize > 0);

        System.setOut(logOut);
        //System.setOut(origOut);
        System.setErr(origErr);

        System.out.println("Trees of size " + trees.get(0).vertices.size() + " out of a possible " + I.normalReads[0].length);
        System.out.println("Found " + biggestTrees + " trees accounting for maximal frequency out of " + totalTrees + " total");
        System.out.println("(discarded " + (discardedTrees - biggestTrees) + " infeasible trees with higher frequencies)");
        System.out.println("Found " + coll.getNumFeasible() +  " feasible trees of " + biggestTrees);
        System.out.println("Found " + coll.getBestResultsSize() + " trees with optimal objective of " + coll.getBestScore());
        System.out.println("Printing individual solutions:");
        System.out.println(coll.printBestResults());

        //System.out.println(coll.printAllResults());

        try {
            int c = 0;
            for(ILPResult res : coll.getBestResults()){
                if (c > MAX_NUM_OPTIMA_OUTPUT){
                    break;
                }
                PrintWriter writer = new PrintWriter(new File(OUTDIR + "tree" + c + ".txt"));
                writer.write(res.printConcise());
                writer.close();
                c++;
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

    }


    public static void log(String message){
        LocalTime time = LocalTime.now();
        System.out.println(time + " " + message);
    }

}
