import net.sf.javailp.Constraint;
import net.sf.javailp.Result;

import java.io.*;
import java.math.BigDecimal;
import java.time.LocalTime;
import java.util.ArrayList;
import java.util.ConcurrentModificationException;
import java.util.LinkedList;


public class Main {
    static String INFILE = null;
    static String OUTDIR = "";

    static String OUTFILE = "CALDER_output.txt";
    static int MAX_NUM_OPTIMA_OUTPUT = 3;

    // don't need this if using Gurobi - automatically uses all available processors
    //public static int THREADS = 3;
    static double CONFIDENCE = .9;
    static double MINIMUM_USAGE_H = .01;

    static boolean PRINT_ANCESTRY_GRAPH = false;
    static boolean PRINT_CONFIDENCE_INTERVALS = false;
    static boolean PRINT_EFFECTIVE_CONFIDENCE = true;

    public static void main(String[] args) {
        INFILE = "CLL003_clustered.txt";

        //INFILE = "../Results/Data/CALDER format/rz_9_readcounts.txt";


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

        PrintStream logOut = null;
        try {
            logOut = new PrintStream(OUTDIR + OUTFILE);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        //System.setOut(logOut);
        //System.setOut(dummy);
        //System.setErr(dummy);

        // solve big ILP
        System.out.println("Starting to solve ILP");

        Result r = Calder.solve(I, G);

        // construct an ILPResult object which unpacks the ILP variables
        ILPResult result = new ILPResult(r, I, G);
        System.out.println(result);
        System.out.println(result.totalPurity());

        int total = 0;
        for(int t = 0; t < I.nSamples; t++){
            for(int i = 0; i < I.nMuts; i++){
                total += (int) r.getPrimalValue("q_" + t + "_" + i);
            }
        }
        //System.out.println("Total q= " + total);

        // Write output solution to file tree0
        int maximal = result.nClones;
        int c = 0;
        try {
            PrintWriter writer = new PrintWriter(new File(OUTDIR + "tree" + c + ".txt"));
            writer.write(result + "");
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        c++;

        //Find additional solutions by adding a trivial constraint for each previous solution
        LinkedList<Constraint> extraConstraints = new LinkedList<>();
        while(c < MAX_NUM_OPTIMA_OUTPUT){
            // Add trivial restraint to
            extraConstraints.add(Calder.constructDummyConstraint(result));

            r = Calder.solve(I, G, extraConstraints);
            if(r != null && (result = new ILPResult(r, I, G)).T.vertices.size() == maximal) {
                System.out.println(result);
                try {
                    PrintWriter writer = new PrintWriter(new File(OUTDIR + "tree" + c + ".txt"));
                    writer.write(result + "");
                } catch (FileNotFoundException e) {
                    e.printStackTrace();
                }
                c++;
            } else {
                System.out.println("No more maximal trees");
                break;
            }
        }

    }


}
