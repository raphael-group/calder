import net.sf.javailp.*;
import net.sf.javailp.Result;

import java.math.BigDecimal;
import java.util.LinkedList;
import java.util.List;

public class Calder {

    private static SolverFactory factory;
    private static Solver solver;

    public static void init(){
        factory = new SolverFactoryGurobi();
        factory.setParameter(Solver.VERBOSE, 0);
        factory.setParameter(Solver.TIMEOUT, 30);

        solver = factory.get();

    }

    public static Result solve(Instance I, Graph G){
        return solve(I, G, new LinkedList<>());
    }

    public static Result solve(Instance I, Graph G,  List<Constraint> extraConstraints){
        Problem problem = new Problem();
        int nClones = I.nMuts;
        int nSamples = I.nSamples;
        int t, i, j, k;
        Linear linear;

        for(Constraint c : extraConstraints){
            problem.add(c);
        }

        // Add dummy root vertex to G with edges to all other vertices
        G.addVertex(nClones);
        for(j = 0; j < nClones; j++){
            G.addEdge(nClones, j);
        }

        // Construct objective function
        setL0Norm(problem, I, G, nSamples, nClones);

        // Define dummy product variable a_tij = fhat_tj * x_ij
        for(t = 0; t < nSamples; t++){
            for (i = 0; i < nClones + 1; i++){
                for(Integer my_j : G.outEdges.get(i)){
                    linear = new Linear();
                    linear.add(1, "a_" + t + "_" + i + "_" + my_j);
                    linear.add(-1, "x_" + i + "_" + my_j);
                    problem.add(new Constraint(linear, "<=", 0));

                    linear = new Linear();
                    linear.add(1, "a_" + t + "_" + i + "_" + my_j);
                    linear.add(-1, "fhat_" + t + "_" + my_j);
                    problem.add(new Constraint(linear, "<=", 0));

                    linear = new Linear();
                    linear.add(1, "a_" + t + "_" + i + "_" + my_j);
                    linear.add(-1, "fhat_" + t + "_" + my_j);
                    linear.add(-1, "x_" + i + "_" + my_j);
                    problem.add(new Constraint(linear, ">=", -1));

                    linear = new Linear();
                    linear.add(1, "a_" + t + "_" + i + "_" + my_j);
                    problem.add(new Constraint(linear, ">=", 0));
                }
            }
        }


        // Define dummy product variable b_tij = u_tj * x_ij
        for(t = 0; t < nSamples; t++){
            for (i = 0; i < nClones; i++){
                for(Integer my_j : G.outEdges.get(i)){
                    linear = new Linear();
                    linear.add(1, "b_" + t + "_" + i + "_" + my_j);
                    linear.add(-1, "x_" + i + "_" + my_j);
                    problem.add(new Constraint(linear, "<=", 0));

                    linear = new Linear();
                    linear.add(1, "b_" + t + "_" + i + "_" + my_j);
                    linear.add(-1, "u_" + t + "_" + my_j);
                    problem.add(new Constraint(linear, "<=", 0));

                    linear = new Linear();
                    linear.add(1, "b_" + t + "_" + i + "_" + my_j);
                    linear.add(-1, "u_" + t + "_" + my_j);
                    linear.add(-1, "x_" + i + "_" + my_j);
                    problem.add(new Constraint(linear, ">=", -1));

                    linear = new Linear();
                    linear.add(1, "b_" + t + "_" + i + "_" + my_j);
                    problem.add(new Constraint(linear, ">=", 0));
                }
            }
        }

        // Define dummy product variable c_ij = x_ij * tmin_j
        for (i = 0; i < nClones + 1; i++){
            for(Integer my_j : G.outEdges.get(i)){
                linear = new Linear();
                linear.add(1, "c_" + i + "_" + my_j);
                linear.add(-1 * nSamples, "x_" + i + "_" + my_j);
                problem.add(new Constraint(linear, "<=", 0));

                linear = new Linear();
                linear.add(1, "c_" + i + "_" + my_j);
                linear.add(-1, "tmin_" + my_j);
                problem.add(new Constraint(linear, "<=", 0));

                linear = new Linear();
                linear.add(1, "c_" + i + "_" + my_j);
                linear.add(-1, "tmin_" + my_j);
                linear.add(-1 * nClones, "x_" + i + "_" + my_j);
                problem.add(new Constraint(linear, ">=", -1 * nClones));

                linear = new Linear();
                linear.add(1, "c_" + i + "_" + my_j);
                problem.add(new Constraint(linear, ">=", 0));
            }
        }

        // The root (index nClones) must have exactly 1 outgoing edge
        linear = new Linear();
        for(j = 0; j < nClones; j++){ // using all clones instead of delta_r because they are equivalent
            linear.add(1, "x_" + nClones + "_" + j);
        }
        problem.add(new Constraint(linear, "=", 1));


        // Each NON-ROOT vertex must have at least as many incoming edges as outgoing
        for(j = 0; j < nClones; j++){
            for(Integer my_k : G.outEdges.get(j)){
                linear = new Linear();
                linear.add(-1, "x_" + j + "_" + my_k);
                for(Integer my_i : G.inEdges.get(j)) {
                    linear.add(1, "x_" + my_i + "_" + j);
                }
                problem.add(new Constraint(linear, ">=", 0));
            }
        }

        // Each vertex can have at most 1 incoming edge
        for(j = 0; j < nClones; j++){
            linear = new Linear();
            for(Integer my_i : G.inEdges.get(j)){
                linear.add(1, "x_" + my_i + "_" + j);
            }
            problem.add(new Constraint(linear, "<=", 1));
        }


        // Add constraints for fhat and u
        for(t = 0; t < nSamples; t++){
            // for each sample, sum across clones of u must be <= 1
            linear = new Linear();
            for(i = 0; i < nClones; i++){
                //linear.add(1, "u_" + t + "_" + i);
                for(Integer my_j : G.outEdges.get(i)){
                    linear.add(1, "b_" + t + "_" + i + "_" + my_j);

                }
            }
            problem.add(new Constraint(linear, "<=", 1));

            for(i = 0; i < nClones; i++){
                //fhat is bounded by intervals
                linear = new Linear();
                linear.add(1, "fhat_" + t + "_" + i);
                assert I.intervals[0][t][i] < I.intervals[1][t][i];
                assert I.intervals[0][t][i] >= 0;
                problem.add(new Constraint(linear, ">=", new BigDecimal(I.intervals[0][t][i])));
                problem.add(new Constraint(linear, "<=", I.intervals[1][t][i]));

                //u must be non-negative
                linear = new Linear();
                linear.add(1, "u_" + t + "_" + i);
                problem.add(new Constraint(linear, ">=", 0));

                // Add sum condition equivalence
                linear = new Linear();
                linear.add(1, "u_" + t + "_" + i);
                linear.add(-1, "fhat_" + t + "_" + i);
                for(Integer my_j : G.outEdges.get(i)){
                    linear.add(1, "a_" + t + "_" + i + "_" + my_j);
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

                // z = 1 -> u = 0
                linear = new Linear();
                linear.add(1, "u_" + t + "_" + i);
                linear.add(1, "z_" + t + "_" + i);
                problem.add(new Constraint(linear, "<=", 1));

                // u >= h for tmin <= t < tmax
                linear = new Linear();
                linear.add(1, "u_" + t + "_" + i);
                linear.add(-1, "y_" + t + "_" + i);
                linear.add(1, "z_" + t + "_" + i);
                problem.add(new Constraint(linear, ">=", Main.MINIMUM_USAGE_H - 1));
            }
        }

        // lineage continuity: x_ij * tmin_j <= tmax_i
        for(i = 0; i < nClones; i++){
            for(Integer my_j : G.outEdges.get(i)){
                linear = new Linear();
                linear.add(1, "c_" + i + "_" + my_j);
                linear.add(-1, "tmax_" + i);
                problem.add(new Constraint(linear, "<=", 0));
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
                problem.setVarType("y_" + t + "_" + i, VarType.BOOL);
                problem.setVarType("z_" + t + "_" + i, VarType.BOOL);
            }


        }

        // x is binary, c is integer
        for(i = 0; i < nClones + 1; i++){
            for(Integer my_j : G.outEdges.get(i)){
                problem.setVarType("x_" + i + "_" + my_j, VarType.BOOL);
                problem.setVarType("c_" + i + "_" + my_j, VarType.INT);
            }
        }

        try{
            return solver.solve(problem);
        } catch (AssertionError e){
            System.out.println("Failed to find solution");
            //TODO: return dummy ILPResult to avoid NullPointerException
            return null;
        }

    }


    private static void setDefaultObjective(Problem p, Instance I, Graph G, int nSamples, int nClones){
        int t, i;
        // Add a term for each possible edge in the tree (most mutations)

        Linear linear = new Linear();
        for(i = 0; i < nClones + 1; i++){
            for(Integer my_j : G.outEdges.get(i)) {
                linear.add(1000000, "x_" + i + "_" + my_j);
            }
        }
        // Add a term for each frequency matrix entry present in sample (largest CCF mutations)
        for(t = 0; t < nSamples; t++){
            for(i = 0; i < nClones + 1; i++){
                for(Integer my_j : G.outEdges.get(i)) {
                    linear.add(100000.0 * I.intervals[0][t][my_j] / ((double) nClones * nSamples), "x_" + i + "_" + my_j);
                }
            }
        }

        // Add a term for each (included) usage matrix entry
        // For any tj, b_tij is nonzero for at most one incoming edge (i,j) because x specifies a tree
        for(t = 0; t < nSamples; t++) {
            for (i = 0; i < nClones; i++) {
                for(Integer my_j : G.outEdges.get(i)) {
                    linear.add(((double) -1) / ((double) nClones * nSamples), "b_" + t + "_" + i + "_" + my_j);
                }
            }
        }
        p.setObjective(linear, OptType.MAX);
    }

    private static void setL0Norm(Problem problem, Instance I, Graph G, int nSamples, int nClones){
        int t, i, j;
        Linear linear;


        // Construct dummy variables q_tj = 1 {j is in tree, u_tj = 0}
        for(t = 0; t < nSamples; t++){
            for(j = 0; j < nClones; j++){
                problem.setVarType("q_" + t + "_" + j, VarType.BOOL);

                linear = new Linear();
                linear.add(1, "q_" + t + "_" + j);
                for(Integer my_i : G.inEdges.get(j)){
                    linear.add(-1, "x_" + my_i + "_" + j);
                }
                problem.add(new Constraint(linear, "<=" , 0));

                linear = new Linear();
                linear.add(1, "q_" + t + "_" + j);
                linear.add(1, "u_" + t + "_" + j);
                problem.add(new Constraint(linear, "<=", 1));
            }
        }

        // Add a term for each possible edge in the tree (most mutations)
        linear = new Linear();
        for(i = 0; i < nClones + 1; i++){
            for(Integer my_j : G.outEdges.get(i)) {
                linear.add(1000000, "x_" + i + "_" + my_j);
            }
        }
        // Add a term for each frequency matrix entry present in sample (largest CCF mutations)
        for(t = 0; t < nSamples; t++){
            for(i = 0; i < nClones + 1; i++){
                for(Integer my_j : G.outEdges.get(i)) {
                    linear.add(100000.0 * I.intervals[0][t][my_j] / ((double) nClones * nSamples), "x_" + i + "_" + my_j);
                }
            }
        }

        /*
        // Add a term for each (included) usage matrix entry
        // For any tj, b_tij is nonzero for at most one incoming edge (i,j) because x specifies a tree
        for(t = 0; t < nSamples; t++) {
            for (i = 0; i < nClones; i++) {
                for(Integer my_j : G.outEdges.get(i)) {
                    linear.add(((double) -1) / ((double) nClones * nSamples), "b_" + t + "_" + i + "_" + my_j);
                }
            }
        }
        */

        // Add a term for each 0 value in U
        for(t = 0; t < nSamples; t++) {
            for (i = 0; i < nClones; i++) {
                linear.add(1, "q_" + t + "_" + i);
            }
        }

        problem.setObjective(linear, OptType.MAX);
    }

    /**
     * In order to not return the same tree twice, constructs a dummy constraint to avoid returning a solution with the
     * same set of edges.
     * @param r ILPResult object constructed from previous solution
     * @return constraint to prevent the edges from being returned
     */
    public static Constraint constructDummyConstraint(ILPResult r){
        System.out.println(r.T.toString());

        Linear linear = new Linear();
        int total = 0;
        for(int i : r.T.vertices){
            for(int j : r.T.outEdges.get(i)){
                linear.add(1, "x_" + i + "_" + j);
                total++;
            }
        }

        // Any other tree must be missing at least 1 edge from this set (because this solution has maximal # of edges)
        return new Constraint(linear, "<=", total - 1);
    }
}
