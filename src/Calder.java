import net.sf.javailp.*;
import net.sf.javailp.Result;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.*;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import static java.lang.System.*;

public class Calder {

    static int opp1 = 0;
    static int helpedPrune = 0;

    private static SolverFactory factory;
    private static Solver solver;
    private static double cap = Double.MAX_VALUE;

    public static void init(){
        factory = new SolverFactoryGLPK();
        factory.setParameter(Solver.VERBOSE, 1);
        factory.setParameter(Solver.TIMEOUT, 30);
        solver = factory.get();
    }

    public static LinkedList<Tree> findBestMaximalTrees(Instance i, LinkedList<Tree> trees){
        if(cap == 0){
            return new LinkedList<>();
        }
        LinkedList<Tree> result = new LinkedList<>();

        double[] columnSums = new double[i.intervals[0][0].length];
        for(int r = 0; r < i.intervals[0].length; r++){
            for(int c = 0; c < i.intervals[0][0].length; c++){
                columnSums[c] += i.intervals[0][r][c];
            }
        }

        double biggest = 0;
        double myTotalF;
        for(Tree T: trees){
            // Count the total frequency accounted for by the vertices in the tree
            myTotalF = 0;
            for(Integer v : T.vertices){
                myTotalF += columnSums[v];
            }

            if(myTotalF > biggest && myTotalF < cap){
                // If I'm the best so far, I'm the only good tree
                biggest = myTotalF;
                result = new LinkedList<>();
                result.add(T);
            } else if (myTotalF == biggest) {
                // If I'm also best, I'm another good tree
                result.add(T);
            }
        }

        cap = biggest;

        return result;
    }

    /**
     * Given a tree T, finds the longitudinal usage matrix U with minimal |U| such that the inferred frequencies are
     * within confidence bounds and every positive usage is at least H
     * @param I Instance object with confidence intervals constructed from read count data
     * @param T enumerated tree
     * @return ILPResult object encapsulating the result OR null if no feasible U was found
     */
    public static ILPResult inferU(Instance I, Tree T){
        Problem problem = new Problem();
        int nClones = T.vertices.size();
        int nSamples = T.nSamples;

        int t, i;
        Linear linear;

        // Compute fbar values (mean of each confidence interval)
        double[][] fbars = new double[nSamples][nClones];
        HashMap<Integer, Integer> idToIndex = new HashMap<>();
        HashMap<Integer, Integer> indexToId = new HashMap<>();
        for(t = 0; t < nSamples; t++){
            i = 0;
            for(Integer v : T.vertices){
                idToIndex.put(v, i);
                indexToId.put(i, v);
                fbars[t][i] = (I.intervals[0][t][v] + I.intervals[1][t][v]) / 2;
                assert fbars[t][i] >= 0;
                assert fbars[t][i] <= 1;
                i++;
            }
        }

        // Construct objective function
        linear = new Linear();
        for(t = 0; t < nSamples; t++){
            for(i = 0; i < nClones; i++){
                // using dummy variable d_tp to represent fbar_tp - fhat_tp
                //linear.add(REG_LAMBDA, "uhat_" + t + "_" + i);
                linear.add(1, "d_" + t + "_" + i);
            }
        }

        problem.setObjective(linear, OptType.MIN);

        // Add constraints for dummy objective variable

        for(t = 0; t < nSamples; t++){
            for(i = 0; i < nClones; i++){
                linear = new Linear();
                linear.add(1, "fhat_" + t + "_" + i);
                linear.add(1, "d_" + t + "_" + i);
                problem.add(new Constraint(linear, ">=", fbars[t][i]));

                linear = new Linear();
                linear.add(1, "fhat_" + t + "_" + i);
                linear.add(-1, "d_" + t + "_" + i);
                problem.add(new Constraint(linear, "<=", fbars[t][i]));
            }
        }


        // Add constraints for fhat and uhat
        for(t = 0; t < nSamples; t++){
            // for each sample, sum across clones of uhat must be <= 1
            linear = new Linear();
            for(i = 0; i < nClones; i++){
                linear.add(1, "uhat_" + t + "_" + i);
            }
            problem.add(new Constraint(linear, "<=", 1));

            for(i = 0; i < nClones; i++){
                //fhat is bounded by intervals
                linear = new Linear();
                linear.add(1, "fhat_" + t + "_" + i);
                assert I.intervals[0][t][indexToId.get(i)] < I.intervals[1][t][indexToId.get(i)];
                assert I.intervals[0][t][indexToId.get(i)] >= 0;
                problem.add(new Constraint(linear, ">=", I.intervals[0][t][indexToId.get(i)]));
                problem.add(new Constraint(linear, "<=", I.intervals[1][t][indexToId.get(i)]));
                //problem.add(new Constraint(linear, ">=", 0));
                //problem.add(new Constraint(linear, "<=", 1));

                //uhat must be non-negative
                linear = new Linear();
                linear.add(1, "uhat_" + t + "_" + i);
                problem.add(new Constraint(linear, ">=", 0));

                // Add sum condition equivalence
                linear = new Linear();
                linear.add(1, "uhat_" + t + "_" + i);
                linear.add(-1, "fhat_" + t + "_" + i);
                if(T.outEdges.containsKey(indexToId.get(i))){
                    for(Integer child : T.outEdges.get(indexToId.get(i))){
                        linear.add(1, "fhat_" + t + "_" + idToIndex.get(child));
                    }
                }
                problem.add(new Constraint(linear, "=", 0));

                // Define constraints using y
                // t < tmin -> y = 0
                linear = new Linear();
                linear.add(1, "tmin_"+ i);
                linear.add(nSamples, "y_" + t + "_" + i);
                problem.add(new Constraint(linear, "<=", t + nSamples));

                // t >= tmin -> y = 1
                linear = new Linear();
                linear.add(1, "tmin_"+ i);
                linear.add(nSamples, "y_" + t + "_" + i);
                problem.add(new Constraint(linear, ">=", t + 1));

                // y = 0 -> fhat = 0
                linear = new Linear();
                linear.add(1, "fhat_" + t + "_" + i);
                linear.add(-1, "y_" + t + "_" + i);
                problem.add(new Constraint(linear, "<=", 0));

                //Define constraints using z
                // t >= tmax -> z = 1
                linear = new Linear();
                linear.add(1, "tmax_" + i);
                linear.add(nSamples, "z_" + t + "_" + i);
                problem.add(new Constraint(linear, ">=", t + 1));

                // t < tmax -> z = 0
                linear = new Linear();
                linear.add(1, "tmax_" + i);
                linear.add(nSamples, "z_" + t + "_" + i);
                problem.add(new Constraint(linear, "<=", t + nSamples));

                // z = 1 -> uhat = 0
                linear = new Linear();
                linear.add(1, "uhat_" + t + "_" + i);
                linear.add(1, "z_" + t + "_" + i);
                problem.add(new Constraint(linear, "<=", 1));

                // uhat >= h for tmin <= t < tmax
                linear = new Linear();
                linear.add(1, "uhat_" + t + "_" + i);
                linear.add(-1, "y_" + t + "_" + i);
                linear.add(1, "z_" + t + "_" + i);
                problem.add(new Constraint(linear, ">=", Main.DETECTION_THRESHOLD_H - 1));
            }
        }

        //lineage continuity
        for(Integer parent : T.vertices){
            if(T.outEdges.containsKey(parent)){
                for(Integer child : T.outEdges.get(parent)){
                    linear = new Linear();
                    linear.add(1, "tmin_" + idToIndex.get(child));
                    linear.add(-1, "tmax_" + idToIndex.get(parent));
                    problem.add(new Constraint(linear, "<=", 0));

                    /*
                    // ancestry condition
                    for(t = 0; t < nSamples; t++){
                        linear = new Linear();
                        linear.add(1, "fhat_" + t + "_" + idToIndex.get(parent));
                        linear.add(-1, "fhat_" + t + "_" + idToIndex.get(child));
                        problem.add(new Constraint(linear, ">=", 0));
                    }
                    */
                }
            }
        }


        // Set up variable definitions (types, constraints for integers)
        for(i = 0; i < nClones; i++) {
            problem.setVarType("tmin_" + i, Integer.class);
            problem.setVarType("tmax_" + i, Integer.class);

            //tmin is in [0, m]
            linear = new Linear();
            linear.add(1, "tmin_" + i);
            problem.add(new Constraint(linear, "<=", nSamples));
            linear = new Linear();
            linear.add(1, "tmin_" + i);
            problem.add(new Constraint(linear, ">=", 0));

            // tmax is in [0, m + 1]
            linear = new Linear();
            linear.add(1, "tmax_" + i);
            problem.add(new Constraint(linear, "<=", nSamples + 1));
            linear = new Linear();
            linear.add(1, "tmax_" + i);
            problem.add(new Constraint(linear, ">=", 0));

            // tmin must be at least as small as tmax
            linear = new Linear();
            linear.add(1, "tmin_" + i);
            linear.add(-1, "tmax_" + i);
            problem.add(new Constraint(linear, "<=", 0));


            for(t = 0; t < nSamples; t++){
                problem.setVarType("d_" + t + "_" + i, Double.class);
                problem.setVarType("fhat_" + t + "_" + i, Double.class);
                problem.setVarType("uhat_" + t + "_" + i, Double.class);

                problem.setVarType("y_" + t + "_" + i, Integer.class);
                problem.setVarType("z_" + t + "_" + i, Integer.class);

                // y is binary
                linear = new Linear();
                linear.add(1, "y_" + t + "_" + i);
                problem.add(new Constraint(linear, "<=", 1));
                linear = new Linear();
                linear.add(1, "y_" + t + "_" + i);
                problem.add(new Constraint(linear, ">=", 0));

                // z is binary
                linear = new Linear();
                linear.add(1, "z_" + t + "_" + i);
                problem.add(new Constraint(linear, "<=", 1));
                linear = new Linear();
                linear.add(1, "z_" + t + "_" + i);
                problem.add(new Constraint(linear, ">=", 0));

            }
        }

        try{
            solver = factory.get();
            Result result = solver.solve(problem);


            return new ILPResult(result, nClones, nSamples, indexToId, idToIndex, T, I.rowLabels, I.colLabels);
        } catch (AssertionError e){
            System.err.println("Failed to find solution");
            //TODO: return dummy ILPResult to avoid NullPointerException
            return null;
        }

    }

    /**
     * Same as inferU except 1) infers |U| to minimize |Fhat-F| and 2) does not use minimum usage h
     * @param I instance of LVAFFP with CIs
     * @param T enumerated tree
     * @return ILPResult object encapsulating the result OR null if no feasible U was found
     */
    public static ILPResult inferUwithoutH(Instance I, Tree T){
        Problem problem = new Problem();
        int nClones = T.vertices.size();
        int nSamples = T.nSamples;

        int t, i;
        Linear linear;

        // Compute fbar values (mean of each confidence interval)
        double[][] fbars = new double[nSamples][nClones];
        HashMap<Integer, Integer> idToIndex = new HashMap<>();
        HashMap<Integer, Integer> indexToId = new HashMap<>();
        for(t = 0; t < nSamples; t++){
            i = 0;
            for(Integer v : T.vertices){
                idToIndex.put(v, i);
                indexToId.put(i, v);
                fbars[t][i] = (I.intervals[0][t][v] + I.intervals[1][t][v]) / 2;
                assert fbars[t][i] >= 0;
                assert fbars[t][i] <= 1;
                i++;
            }
        }

        // Construct objective function
        linear = new Linear();
        for(t = 0; t < nSamples; t++){
            for(i = 0; i < nClones; i++){
                // using dummy variable d_tp to represent fbar_tp - fhat_tp
                //linear.add(REG_LAMBDA, "uhat_" + t + "_" + i);
                linear.add(1, "d_" + t + "_" + i);
            }
        }

        problem.setObjective(linear, OptType.MIN);

        // Add constraints for dummy objective variable

        for(t = 0; t < nSamples; t++){
            for(i = 0; i < nClones; i++){
                linear = new Linear();
                linear.add(1, "fhat_" + t + "_" + i);
                linear.add(1, "d_" + t + "_" + i);
                problem.add(new Constraint(linear, ">=", fbars[t][i]));

                linear = new Linear();
                linear.add(1, "fhat_" + t + "_" + i);
                linear.add(-1, "d_" + t + "_" + i);
                problem.add(new Constraint(linear, "<=", fbars[t][i]));
            }
        }


        // Add constraints for fhat and uhat
        for(t = 0; t < nSamples; t++){
            // for each sample, sum across clones of uhat must be <= 1
            linear = new Linear();
            for(i = 0; i < nClones; i++){
                linear.add(1, "uhat_" + t + "_" + i);
            }
            problem.add(new Constraint(linear, "<=", 1));

            for(i = 0; i < nClones; i++){
                //fhat is bounded by intervals
                linear = new Linear();
                linear.add(1, "fhat_" + t + "_" + i);
                assert I.intervals[0][t][indexToId.get(i)] < I.intervals[1][t][indexToId.get(i)];
                assert I.intervals[0][t][indexToId.get(i)] >= 0;
                problem.add(new Constraint(linear, ">=", I.intervals[0][t][indexToId.get(i)]));
                problem.add(new Constraint(linear, "<=", I.intervals[1][t][indexToId.get(i)]));
                //problem.add(new Constraint(linear, ">=", 0));
                //problem.add(new Constraint(linear, "<=", 1));

                //uhat must be non-negative
                linear = new Linear();
                linear.add(1, "uhat_" + t + "_" + i);
                problem.add(new Constraint(linear, ">=", 0));

                // Add sum condition equivalence
                linear = new Linear();
                linear.add(1, "uhat_" + t + "_" + i);
                linear.add(-1, "fhat_" + t + "_" + i);
                if(T.outEdges.containsKey(indexToId.get(i))){
                    for(Integer child : T.outEdges.get(indexToId.get(i))){
                        linear.add(1, "fhat_" + t + "_" + idToIndex.get(child));
                    }
                }
                problem.add(new Constraint(linear, "=", 0));

                // Define constraints using y
                // t < tmin -> y = 0
                linear = new Linear();
                linear.add(1, "tmin_"+ i);
                linear.add(nSamples, "y_" + t + "_" + i);
                problem.add(new Constraint(linear, "<=", t + nSamples));

                // t >= tmin -> y = 1
                linear = new Linear();
                linear.add(1, "tmin_"+ i);
                linear.add(nSamples, "y_" + t + "_" + i);
                problem.add(new Constraint(linear, ">=", t + 1));

                // y = 0 -> uhat = 0
                linear = new Linear();
                linear.add(1, "fhat_" + t + "_" + i);
                linear.add(-1, "y_" + t + "_" + i);
                problem.add(new Constraint(linear, "<=", 0));

                //Define constraints using z
                // t >= tmax -> z = 1
                linear = new Linear();
                linear.add(1, "tmax_" + i);
                linear.add(nSamples, "z_" + t + "_" + i);
                problem.add(new Constraint(linear, ">=", t + 1));

                // t < tmax -> z = 0
                linear = new Linear();
                linear.add(1, "tmax_" + i);
                linear.add(nSamples, "z_" + t + "_" + i);
                problem.add(new Constraint(linear, "<=", t + nSamples));

                // z = 1 -> uhat = 0
                linear = new Linear();
                linear.add(1, "uhat_" + t + "_" + i);
                linear.add(1, "z_" + t + "_" + i);
                problem.add(new Constraint(linear, "<=", 1));

                // uhat >= h for tmin <= t < tmax
                linear = new Linear();
                linear.add(1, "uhat_" + t + "_" + i);
                linear.add(-1, "y_" + t + "_" + i);
                linear.add(1, "z_" + t + "_" + i);
                problem.add(new Constraint(linear, ">=", - 1 + .000001));
            }
        }

        //lineage continuity
        for(Integer parent : T.vertices){
            if(T.outEdges.containsKey(parent)){
                for(Integer child : T.outEdges.get(parent)){
                    linear = new Linear();
                    linear.add(1, "tmin_" + idToIndex.get(child));
                    linear.add(-1, "tmax_" + idToIndex.get(parent));
                    problem.add(new Constraint(linear, "<=", 0));

                    // ancestry condition
                    for(t = 0; t < nSamples; t++){
                        linear = new Linear();
                        linear.add(1, "fhat_" + t + "_" + idToIndex.get(parent));
                        linear.add(-1, "fhat_" + t + "_" + idToIndex.get(child));
                        problem.add(new Constraint(linear, ">=", 0));
                    }
                }
            }
        }


        // Set up variable definitions (types, constraints for integers)
        for(i = 0; i < nClones; i++) {
            problem.setVarType("tmin_" + i, Integer.class);
            problem.setVarType("tmax_" + i, Integer.class);

            //tmin is in [0, m]
            linear = new Linear();
            linear.add(1, "tmin_" + i);
            problem.add(new Constraint(linear, "<=", nSamples));
            linear = new Linear();
            linear.add(1, "tmin_" + i);
            problem.add(new Constraint(linear, ">=", 0));

            // tmax is in [0, m]
            linear = new Linear();
            linear.add(1, "tmax_" + i);
            problem.add(new Constraint(linear, "<=", nSamples));
            linear = new Linear();
            linear.add(1, "tmax_" + i);
            problem.add(new Constraint(linear, ">=", 0));

            // tmin must be at least as small as tmax
            linear = new Linear();
            linear.add(1, "tmin_" + i);
            linear.add(-1, "tmax_" + i);
            problem.add(new Constraint(linear, "<=", 0));


            for(t = 0; t < nSamples; t++){
                problem.setVarType("d_" + t + "_" + i, Double.class);
                problem.setVarType("fhat_" + t + "_" + i, Double.class);
                problem.setVarType("uhat_" + t + "_" + i, Double.class);

                problem.setVarType("y_" + t + "_" + i, Integer.class);
                problem.setVarType("z_" + t + "_" + i, Integer.class);

                // y is binary
                linear = new Linear();
                linear.add(1, "y_" + t + "_" + i);
                problem.add(new Constraint(linear, "<=", 1));
                linear = new Linear();
                linear.add(1, "y_" + t + "_" + i);
                problem.add(new Constraint(linear, ">=", 0));

                // z is binary
                linear = new Linear();
                linear.add(1, "z_" + t + "_" + i);
                problem.add(new Constraint(linear, "<=", 1));
                linear = new Linear();
                linear.add(1, "z_" + t + "_" + i);
                problem.add(new Constraint(linear, ">=", 0));

            }
        }
        //System.out.println(problem.getConstraints());

        try{
            solver = factory.get();
            Result result = solver.solve(problem);

            return new ILPResult(result, nClones, nSamples, indexToId, idToIndex, T, I.rowLabels, I.colLabels);
        } catch (AssertionError e){
            System.err.println("Failed to find solution");
            //TODO: return dummy ILPResult
            return null;
        }

    }

    /**
     * Infers a usage matrix within confidence interval bounds to minimize |U| (no other constraints)
     * @param I instance of LVAFFP with CIs
     * @param T enumerated treee
     * @return ILPResult object encapsulating the result OR null if no feasible U was found
     */
    public static ILPResult inferAnyU(Instance I, Tree T){
        Problem problem = new Problem();
        int nClones = T.vertices.size();
        int nSamples = T.nSamples;

        int t, i;
        Linear linear;
        double[][] fbars = new double[nSamples][nClones];
        HashMap<Integer, Integer> idToIndex = new HashMap<>();
        HashMap<Integer, Integer> indexToId = new HashMap<>();
        for(t = 0; t < nSamples; t++){
            i = 0;
            for(Integer v : T.vertices){
                idToIndex.put(v, i);
                indexToId.put(i, v);
                i++;
            }
        }

        // Construct objective function
        linear = new Linear();
        for(t = 0; t < nSamples; t++){
            for(i = 0; i < nClones; i++){
                linear.add(1, "uhat_" + t + "_" + i);
            }
        }

        problem.setObjective(linear, OptType.MAX);

        // Add constraints for fhat and uhat
        for(t = 0; t < nSamples; t++){
            // for each sample, sum across clones of uhat must be <= 1
            linear = new Linear();
            for(i = 0; i < nClones; i++){
                linear.add(1, "uhat_" + t + "_" + i);
            }
            problem.add(new Constraint(linear, "<=", 1));

            for(i = 0; i < nClones; i++){
                //fhat is bounded by intervals
                linear = new Linear();
                linear.add(1, "fhat_" + t + "_" + i);
                assert I.intervals[0][t][indexToId.get(i)] < I.intervals[1][t][indexToId.get(i)];
                assert I.intervals[0][t][indexToId.get(i)] >= 0;
                problem.add(new Constraint(linear, ">=", I.intervals[0][t][indexToId.get(i)]));
                problem.add(new Constraint(linear, "<=", I.intervals[1][t][indexToId.get(i)]));
                //problem.add(new Constraint(linear, ">=", 0));
                //problem.add(new Constraint(linear, "<=", 1));

                //uhat must be non-negative
                linear = new Linear();
                linear.add(1, "uhat_" + t + "_" + i);
                problem.add(new Constraint(linear, ">=", 0));

                // Add sum condition equivalence
                linear = new Linear();
                linear.add(1, "uhat_" + t + "_" + i);
                linear.add(-1, "fhat_" + t + "_" + i);
                if(T.outEdges.containsKey(indexToId.get(i))){
                    for(Integer child : T.outEdges.get(indexToId.get(i))){
                        linear.add(1, "fhat_" + t + "_" + idToIndex.get(child));
                    }
                }
                problem.add(new Constraint(linear, "=", 0));

                /*
                // y = 0 -> u >= h
                linear = new Linear();
                linear.add(1, "uhat_" + t + "_" + i);
                linear.add(DETECTION_THRESHOLD_H, "y_" + t + "_" + i);
                problem.add(new Constraint(linear, ">=", DETECTION_THRESHOLD_H));

                // y = 1 -> u <= 0
                linear = new Linear();
                linear.add(1, "uhat_" + t + "_" + i);
                linear.add(1, "y_" + t + "_" + i);
                problem.add(new Constraint(linear, "<=", 1));
                */
            }
        }

        // Set up variable definitions (types, constraints for integers)
        for(i = 0; i < nClones; i++) {

            for(t = 0; t < nSamples; t++){
                problem.setVarType("fhat_" + t + "_" + i, Double.class);
                problem.setVarType("uhat_" + t + "_" + i, Double.class);

                /*
                problem.setVarType("y_" + t + "_" + i, Integer.class);

                // y is binary
                linear = new Linear();
                linear.add(1, "y_" + t + "_" + i);
                problem.add(new Constraint(linear, "<=", 1));
                linear = new Linear();
                linear.add(1, "y_" + t + "_" + i);
                problem.add(new Constraint(linear, ">=", 0));
                */
            }
        }

        try{
            solver = factory.get();
            Result result = solver.solve(problem);

            return new ILPResult(result, nClones, nSamples, indexToId, idToIndex, T, I.rowLabels, I.colLabels);
        } catch (AssertionError e){
            System.err.println("Failed to find solution");
            return null;
        }

    }


    /**
     * Solve the LVAFFP represented by the input frequency matrix
     * @param I instance of the LVAFFP
     * @return list of enumerated trees
     */
    public static LinkedList<Tree> solveLVAFFP(Instance I, int doneSize){
        assert I != null;
        double[][][] F = I.intervals;
        assert F.length == 2;
        assert F[0].length > 0;
        assert F[0][0].length > 0;
        assert F[0].length == F[1].length;
        assert F[0][0].length == F[1][0].length;

        int m = F[0].length;
        int n = F[0][0].length;
        double[][] mins = F[0];
        double[][] maxes = F[1];

        ArrayList<VertexData> mutations = new ArrayList<>();
        int label = 0, tmin = -1;
        double[] myMin, myMax;
        for(int i = 0; i < n; i++){
            myMin = new double[m];
            myMax = new double[m];
            for(int t = 0; t < m; t++){
                myMin[t] = mins[t][i];
                if (tmin < 0 && myMin[t] > Main.DETECTION_THRESHOLD_H){
                    tmin = t;
                }
                myMax[t] = maxes[t][i];
            }
            mutations.add(new VertexData(label++, myMin, myMax, tmin));
        }

        Graph G = Graph.buildAncestryGraph(mutations);
        //TODO: add flag for printing ancestry graph
        //System.out.println(G);

        GabowMyersEnum gm = new GabowMyersEnum(G, I.nSamples);
        return gm.run(doneSize);
    }


    // Pretty-prints a 2D array with column and row labels
    static String print2DArray(int[][] arr, String[] colLabels, String[] rowLabels){
        assert arr != null;
        assert colLabels != null;
        assert colLabels.length == arr[0].length;
        assert rowLabels != null;
        assert rowLabels.length == arr.length;

        StringBuilder b = new StringBuilder();
        int i, j;
        b.append("\t ");
        for(i = 0; i < arr[0].length; i++){
            b.append(colLabels[i]);
            b.append(spaces(6 - colLabels[i].length()));
        }
        b.append('\n');
        for(i = 0; i < arr.length; i++){
            b.append(rowLabels[i]);
            b.append(": ");
            for(j = 0; j < arr[i].length; j++){
                b.append(arr[i][j]);
                b.append(spaces(6 - (arr[i][j] + "").length()));
            }
            b.append('\n');
        }

        return b.toString();
    }

    static String print2DArray(int[][] arr, int[] colLabels, int[] rowLabels){
        int i;
        String[] newCols = new String[colLabels.length];
        for(i = 0; i < colLabels.length; i++){
            newCols[i] = colLabels[i] + "";
        }
        String[] newRows = new String[rowLabels.length];
        for(i = 0; i < rowLabels.length; i++){
            newRows[i] = rowLabels[i] + "";
        }
        return print2DArray(arr, newCols, newRows);
    }

    static String print2DArray(int[][] arr, int[] colLabels){
        int i;
        String[] newCols = new String[colLabels.length];
        for(i = 0; i < colLabels.length; i++){
            newCols[i] = colLabels[i] + "";
        }
        String[] rowLabels = new String[arr.length];
        for(i = 0; i < arr.length; i++){
            rowLabels[i] = "(" + i + ")";
        }
        return print2DArray(arr, newCols, rowLabels);
    }

    static String print2DArray(int[][] arr){
        assert arr != null;
        assert arr.length > 0;

        String[] rowLabels = new String[arr.length];
        String[] colLabels = new String[arr[0].length];

        int i = 0;
        int max = arr.length + arr[0].length;

        while (i < max) {
            if(i < arr.length){
                rowLabels[i] = "(" + i + ")";
            }
            if(i < arr[0].length){
                colLabels[i] = "(" + i + ")";
            }
            i++;
        }


        return print2DArray(arr, colLabels, rowLabels);
    }
    static String print2DArray(double[][] arr){
        assert arr != null;
        assert arr.length > 0;
        assert arr[0].length > 0;

        String[] rowLabels = new String[arr.length];
        String[] colLabels = new String[arr[0].length];

        int i = 0;
        int max = arr.length + arr[0].length;

        while (i < max) {
            if(i < arr.length){
                rowLabels[i] = "(" + i + ")";
            }
            if(i < arr[0].length){
                colLabels[i] = "(" + i + ")";
            }
            i++;
        }


        return print2DArray(arr, colLabels, rowLabels);
    }

    static String print2DArray(double[][] arr, String[] colLabels, String[] rowLabels){
        assert arr != null;
        assert colLabels != null;
        assert colLabels.length == arr[0].length;
        assert rowLabels != null;
        assert rowLabels.length == arr.length;

        DecimalFormat df = new DecimalFormat();
        df.setMaximumFractionDigits(3);

        StringBuilder b = new StringBuilder();
        int i, j;
        b.append("\t ");
        for(i = 0; i < arr[0].length; i++){
            b.append(colLabels[i]);
            b.append("\t");
        }
        b.append('\n');
        for(i = 0; i < arr.length; i++){
            b.append(rowLabels[i]);
            b.append(": ");
            b.append("\t");
            for(j = 0; j < arr[i].length; j++){
                b.append(df.format(arr[i][j]));
                b.append("\t");
            }
            b.append('\n');
        }

        return b.toString();
    }

    // Returns a string with the designated number of spaces
    private static String spaces(int howMany){
        StringBuilder b = new StringBuilder();
        for(int i = 0; i < howMany; i++){
            b.append(' ');
        }
        return b.toString();
    }

    public static void resetCounts(){
        opp1 = 0;
        helpedPrune = 0;
    }
    public static void resetCap(){
        cap = Double.MAX_VALUE;
    }



}
