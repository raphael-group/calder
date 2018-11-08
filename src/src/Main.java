import sun.awt.image.ImageWatched;

import java.io.*;
import java.time.LocalTime;
import java.util.LinkedList;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class Main {
    public static String INFILE = "..\\Results\\CLL003_cl_readcounts.txt";
    public static String OUTDIR = "..\\Output\\";
    public static String LOGFILE = "CALDER_log.txt";
    public static int MAX_NUM_OPTIMA_OUTPUT = 20;
    public static boolean PRINT_CONFIDENCE_INTERVALS = false;
    public static final int THREADS = 6;
    public static double CONFIDENCE = .9;
    public static final double DETECTION_THRESHOLD_H = .01;

    public static void main(String[] args) throws FileNotFoundException {
        // Set up CALDER

        // Load data from file (TSV, basically AncesTree format)
        String fname = INFILE;
        Instance I = Instance.fromFile(fname);

        if(PRINT_CONFIDENCE_INTERVALS){
            System.out.println(I);
            System.out.println(I.intervals.length + " " + I.intervals[0].length + " " + I.intervals[0][0].length);
            System.out.println(I.intervals.length + " " + I.intervals[1].length + " " + I.intervals[1][0].length);
        }

        PrintStream origOut = System.out;
        PrintStream logOut = new PrintStream(OUTDIR + LOGFILE);

        PrintStream dummy = new PrintStream(new OutputStream() {
            @Override
            public void write(int b) throws IOException {

            }});

        int totalTrees;
        int iter = 0;
        int discardedTrees = 0;
        LinkedList<Tree> maximalFTrees;
        int nFeasible, biggestTrees;
        ILPResultCollector coll;

        ExecutorService pool = Executors.newWorkStealingPool(THREADS);
        LinkedList<Tree> trees;
        int doneSize = I.nMuts;

        System.setOut(dummy);

        do {
            trees = Calder.solveLVAFFP(I, doneSize);
            totalTrees = trees.size();
            do {
                // Select only those maximal trees that account for a maximal amount of frequency in F
                maximalFTrees = Calder.findBestMaximalTrees(I, trees);
                biggestTrees = maximalFTrees.size();

                CountDownLatch latch = new CountDownLatch(maximalFTrees.size());
                //LinkedList<SolveILPTask> tasks = new LinkedList<>();
                coll = new ILPResultCollector(I);



                Calder.init();
                SolveILPTask task;
                for (Tree T : maximalFTrees) {
                    task = new SolveILPTask(I, T, latch, coll);
                    //tasks.add(task);
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
            doneSize = trees.get(0).vertices.size() - 1;
        } while(nFeasible == 0 && doneSize > 0);

        System.setOut(logOut);

        System.out.println("Trees of size " + trees.get(0).vertices.size() + " out of a possible " + I.normalReads[0].length);
        System.out.println("Found " + biggestTrees + " trees accounting for maximal frequency out of " + totalTrees + " total");
        System.out.println("(discarded " + discardedTrees + " infeasible trees with higher frequencies)");
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
                PrintWriter writer = new PrintWriter(new File(OUTDIR + "\\tree" + c + ".txt"));
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
