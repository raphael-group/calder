import net.sf.javailp.Constraint;
import net.sf.javailp.Result;

import java.io.*;
import java.text.DecimalFormat;
import java.time.Duration;
import java.time.Instant;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.lang.management.*;

public class Main {
    static String INFILE = null;
    static String OUTDIR = "current/"; //TODO: change this back to empty by default

    static String OUTFILE = "log.txt";
    static int MAX_NUM_OPTIMA_OUTPUT = 1;

    // don't need this if using Gurobi - automatically uses all detectable processors
    //public static int THREADS = 3;

    static double CONFIDENCE = .9;
    static double MINIMUM_USAGE_H = .01;
    static int PRECISION_DIGITS = 6;

    static boolean PRINT_ANCESTRY_GRAPH = true;
    static boolean PRINT_CONFIDENCE_INTERVALS = true;
    static boolean PRINT_EFFECTIVE_CONFIDENCE = true;
    static boolean COUNT = false;

    public static void main(String[] args) {
        long start = System.currentTimeMillis();

        INFILE = "CLL003_clustered.txt";

        //INFILE = "../Results/Data/CALDER format/rz_9_readcounts.txt";
        //INFILE = "../Results/Data/CALDER format/rz_9_readcounts.txt";
        //INFILE = "../Results/Data/CALDER format/rz_9_readcounts.txt";
        //INFILE = "../Results/Data/CALDER format/rz_9_readcounts.txt";
        //INFILE = "../Results/Data/CALDER format/rz_9_readcounts.txt";
        //INFILE = "../Results/Data/CALDER format/rz_9_readcounts.txt";
        //INFILE = "../Results/Data/CALDER format/rz_9_readcounts.txt";
        //INFILE = "../Results/Data/CALDER format/rz_9_readcounts.txt";
        //INFILE = "../Results/Data/CALDER format/rz_9_readcounts.txt";

        //INFILE = "../Results/Data/CALDER format/eirew_a_readcounts.txt";
        //INFILE = "../Results/Data/CALDER format/eirew_b_readcounts.txt";
        //INFILE = "../Results/Data/CALDER format/eirew_494_readcounts.txt";
        //INFILE = "../Results/Data/CALDER format/eirew_a_cl_readcounts.txt";

        //INFILE = "../Results/Data/StJude/SJALL013787_cl_readcounts.txt";
        //INFILE = "../Results/Data/CALDER format/sim1_readcounts.txt";
        //INFILE = "../Results/Data/sim/instance0_cl_readcounts.txt";

        System.out.println(Arrays.toString(args));

        PRINT_EFFECTIVE_CONFIDENCE = false;
        PRINT_CONFIDENCE_INTERVALS = false;
        PRINT_ANCESTRY_GRAPH = false;

        if(args.length > 0){
            INFILE = args[0];
        }


        //////////////////////////////////////////////////////////
        // Set up CALDER
        Calder.init();

        // Load data from file (TSV, basically AncesTree format)
        String fname = INFILE;


        Instance I = Instance.fromFile(fname);
        int n = I.nMuts;
        int m = I.nSamples;
        ArrayList<VertexData> mutations = new ArrayList<>();

        double[][] mins = I.intervals[0];
        double[][] maxes = I.intervals[1];
        double[] myMin, myMax;
        for(int i = 0; i < n; i++){
            myMin = new double[m];
            myMax = new double[m];
            for(int t = 0; t < m; t++){
                myMin[t] = mins[t][i];
                myMax[t] = maxes[t][i];
            }
            mutations.add(new VertexData(i, myMin, myMax));
        }

        if(PRINT_CONFIDENCE_INTERVALS){
            System.out.println(I);
        }

        System.out.println("Building ancestry graph");

        Graph G = Graph.buildAncestryGraph(mutations);
        if(Main.PRINT_ANCESTRY_GRAPH){
            System.out.println("Ancestry graph:");
            System.out.println(G);

        }

        PrintStream origErr = System.err;
        PrintStream origOut = System.out;

        //TODO: maybe uncomment and implement as option
        /*
        PrintStream logOut = null;
        try {
            logOut = new PrintStream(OUTDIR + OUTFILE);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        */
        //System.setOut(logOut);
        //System.setOut(dummy);
        //System.setErr(dummy);

        // solve big ILP
        System.out.println("Starting to solve ILP");

        Result r = Calder.solve(I, G);

        // construct an ILPResult object which unpacks the ILP variables
        ILPResult result = new ILPResult(r, I, G);

        COUNT = false;

        if(!COUNT){
            //System.out.println(result);
        }

        int total = 1; // 1 edge incoming to root
        for(int i = 0; i < I.nMuts; i++){
            for(Integer j : G.outEdges.get(i)){
                total += (int) r.getPrimalValue("x_" + i + "_" + j);
                //System.out.println("x_" + i + "_" + j + " = " + r.getPrimalValue("x_" + i + "_" + j));
            }
        }
        System.out.println("Total of " + total + " vertices included");
        for(int i = 0; i < I.nMuts; i++){
            //System.out.println("d_" + i + " = " + r.getPrimalValue("d_" + i));
        }

        String[] infile_tkns = INFILE.split("/");
        String tkn = infile_tkns[infile_tkns.length - 1].split("\\.")[0].split("_")[0];

        // Write output solution to file tree0
        int maximal = result.nClones;
        int c = 0;
        try {
            PrintWriter writer = new PrintWriter(new File(OUTDIR + tkn + "_" + "tree" + c + ".txt"));
            //writer.write(result.toString() + "");
            if(!COUNT){
                writer.write(result.toStringConcise());
                writer.close();
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        c++;

        //Find additional solutions by adding a trivial constraint for each previous solution
        LinkedList<Constraint> extraConstraints = new LinkedList<>();
        if (COUNT){
            MAX_NUM_OPTIMA_OUTPUT = Integer.MAX_VALUE;
        }
        while(c < MAX_NUM_OPTIMA_OUTPUT){
            // Add trivial restraint to
            extraConstraints.add(Calder.constructDummyConstraint(result, c));

            r = Calder.solve(I, G, extraConstraints, true);
            if(r != null && (result = new ILPResult(r, I, G)).T.vertices.size() == maximal) {
                System.out.println("Solution number " + (c+1) + "------------------------------------------");
                //System.out.println(result);
                try {

                    PrintWriter writer = new PrintWriter(new File(OUTDIR + tkn + "_" + "tree" + c + ".txt"));
                    //writer.write(result.toString() + "");
                    if(!COUNT){
                        writer.write(result.toStringConcise());
                        writer.close();
                    }
                } catch (FileNotFoundException e) {
                    e.printStackTrace();
                }
                c++;
            } else {
                System.out.println("No more maximal trees");
                break;
            }
        }

        System.out.println(c);

        ThreadMXBean bean = ManagementFactory.getThreadMXBean();
        DecimalFormat df = new DecimalFormat();
        df.setMaximumFractionDigits(3);
        long stop = System.currentTimeMillis();
        if(bean.isCurrentThreadCpuTimeSupported()){
            long cpuTime = bean.getCurrentThreadCpuTime(); // system + user time
            long userTime = bean.getCurrentThreadUserTime(); // user time
            long systemTime = bean.getCurrentThreadCpuTime() - bean.getCurrentThreadUserTime();
            System.out.println("wall: " + (stop - start) / ((float) 1000)); // wall time
            System.out.println("user: " + userTime / ((float) 1000000000));
            System.out.println("system: " + systemTime / ((float) 1000000000));
            MemoryMXBean membean =  ManagementFactory.getMemoryMXBean();
            long mem = membean.getNonHeapMemoryUsage().getCommitted() + membean.getHeapMemoryUsage().getCommitted();

            try {
                PrintWriter writer = new PrintWriter(new File("current_time/" + tkn + ".txt"));
                writer.write("wall: " + (stop - start) / ((float) 1000) + "\n"); // wall time
                writer.write("user: " + userTime / ((float) 1000000000) + "\n");
                writer.write("system: " + systemTime / ((float) 1000000000) + "\n");
                writer.write("memory: " + mem / ((float) 1000) + "\n");
                writer.close();
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }
        }
    }


}
