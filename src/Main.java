import net.sf.javailp.Constraint;
import net.sf.javailp.Result;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.LinkedList;
import java.lang.management.*;

public class Main {
    static String INFILE = null;
    static String OUTDIR = "outdir/";

    static int MAX_NUM_OPTIMA_OUTPUT = 1;

    // don't need this if using Gurobi - automatically uses all detectable processors
    //public static int THREADS = 3;

    static double CONFIDENCE = .9;
    static double MINIMUM_USAGE_H = .01;
    static int PRECISION_DIGITS = 6;

    static boolean PRINT_ANCESTRY_GRAPH = false;
    static boolean PRINT_CONFIDENCE_INTERVALS = false;
    static boolean PRINT_EFFECTIVE_CONFIDENCE = false;
    static boolean COUNT = false;
    static boolean TIMING = false;
    static Objective OBJECTIVE = Objective.L0;
    static boolean LONGITUDINAL = true;
    static JavaILPSolver SOLVER = JavaILPSolver.GUROBI;
    static int SOLVER_VERBOSE = 0;
    static int SOLVER_TIMEOUT = 36000;
    static boolean REMOVE_CNA = false;

    enum Objective {
        L0, L1, L0center;
    }
    enum JavaILPSolver{
        LP_SOLVE, GUROBI, GLPK;
    }

    public static void main(String[] args) {

        if (INFILE == null){
            INFILE = "CLL003_clustered.txt";
            //INFILE = "../Simulation/BranchingSimData/exome2/CALDER/exome2_119.txt";
            COUNT = false;
        }

        if (!Files.exists(Paths.get(OUTDIR))) {
            new File(OUTDIR).mkdirs();
        }

        long start = System.currentTimeMillis();

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
        Result r = Calder.solve(I, G);

        // construct an ILPResult object which unpacks the ILP variables
        ILPResult result = new ILPResult(r, I, G);

        if(!COUNT){
            System.out.println(result);
        }

        int total = 1; // 1 edge incoming to root
        for(int i = 0; i < I.nMuts; i++){
            for(Integer j : G.outEdges.get(i)){
                total += (int) r.getPrimalValue("x_" + i + "_" + j);
                //System.out.println("x_" + i + "_" + j + " = " + r.getPrimalValue("x_" + i + "_" + j));
            }
        }
        System.out.println("Maximal tree contains " + total + " out of " + (G.vertices.size() - 1) + " mutations/clusters");
        /*
        // For each mutation/cluster, prints the variable indicating whether or not it was included in the solution
        for(int i = 0; i < I.nMuts; i++){
            System.out.println("d_" + i + " = " + r.getPrimalValue("d_" + i));
        }
        */

        String[] infile_tkns = INFILE.split("/");
        String tkn = infile_tkns[infile_tkns.length - 1].split("\\.")[0].split("_")[0];

        // Write output solution to file tree0
        int maximal = result.nClones;
        int c = 0;
        if(!COUNT){
            writeSolution(tkn, result, c + 1);
        }
        c++;

        //Find additional solutions by adding a trivial constraint for each previous solution
        LinkedList<Constraint> extraConstraints = new LinkedList<>();
        if (COUNT){
            MAX_NUM_OPTIMA_OUTPUT = Integer.MAX_VALUE;
        }
        while(c < MAX_NUM_OPTIMA_OUTPUT){
            // Add trivial constraint to prohibit returning exactly the same tree
            extraConstraints.add(Calder.constructDummyConstraint(result, c));

            r = Calder.solve(I, G, extraConstraints);
            if(r != null && (result = new ILPResult(r, I, G)).T.vertices.size() == maximal) {
                System.out.println("Solution number " + (c+1) + "------------------------------------------");
                writeSolution(tkn, result, c + 1);
                c++;
            } else {
                System.out.println("No more maximal trees");
                break;
            }
        }

        ThreadMXBean bean = ManagementFactory.getThreadMXBean();
        DecimalFormat df = new DecimalFormat();
        df.setMaximumFractionDigits(3);
        long stop = System.currentTimeMillis();
        if(TIMING && bean.isCurrentThreadCpuTimeSupported()){
            long cpuTime = bean.getCurrentThreadCpuTime(); // system + user time
            long userTime = bean.getCurrentThreadUserTime(); // user time
            long systemTime = bean.getCurrentThreadCpuTime() - bean.getCurrentThreadUserTime();
            System.out.println("wall: " + (stop - start) / ((float) 1000)); // wall time
            System.out.println("user: " + userTime / ((float) 1000000000));
            System.out.println("system: " + systemTime / ((float) 1000000000));
            MemoryMXBean membean =  ManagementFactory.getMemoryMXBean();
            long mem = membean.getNonHeapMemoryUsage().getCommitted() + membean.getHeapMemoryUsage().getCommitted();

            try {
                PrintWriter writer = new PrintWriter(new File(OUTDIR + File.separator + tkn + "_" + "time.txt"));
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
    public static void writeSolution(String fname, ILPResult result, int id){
        try {
<<<<<<< Updated upstream
            PrintWriter writer = new PrintWriter(new File(OUTDIR + File.separator + fname + "_" + "tree" + id + ".txt"));
            writer.write(result.treeToString(fname + "_tree" + id));
            writer.close();

            writer = new PrintWriter(new File(OUTDIR + File.separator + fname + "_" + "soln" + id + ".txt"));
=======
            PrintWriter writer = new PrintWriter(new File(OUTDIR + File.separator + fname + "_" + "tree" + id + ".dot"));
            writer.write(result.treeToString(fname + "_tree" + id));
            writer.close();

            writer = new PrintWriter(new File(OUTDIR + File.separator + fname + "_" + "soln" + id + ".csv"));
>>>>>>> Stashed changes
            writer.write(result.solnToString());
            writer.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }


}
