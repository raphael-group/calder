import net.sf.javailp.*;
import net.sf.javailp.Result;

import java.math.BigDecimal;
import java.util.LinkedList;
import java.util.List;

public class Calder {

    private static SolverFactory factory;
    private static Solver solver;

    public static void init(){
        switch(Main.SOLVER){
            case LP_SOLVE:
                factory = new SolverFactoryLpSolve();
                break;
            case CPLEX:
                factory = new SolverFactoryCPLEX();
                break;
            case GUROBI:
                factory = new SolverFactoryGurobi();
                break;
            case MOSEK:
                factory = new SolverFactoryMosek();
                break;
            case GLPK:
                factory = new SolverFactoryGLPK();
                break;
        }
        factory.setParameter(Solver.VERBOSE, 0);
        factory.setParameter(Solver.TIMEOUT, 36000);

        solver = factory.get();

    }

    public static Result solve(Instance I, Graph G){
        return solve(I, G, new LinkedList<>(), true);
    }

    public static Result solveNonLongitudinal(Instance I, Graph G){
        return solve(I, G, new LinkedList<>(), false);
    }

    public static Result solve(Instance I, Graph G,  List<Constraint> extraConstraints, boolean longitudinal){
        Problem problem = new Problem();
        int nClones = I.nMuts;
        int nSamples = I.nSamples;
        int t, i, j, counter = 0;
        Linear linear;

        // Add extra constraints (i.e. don't reproduce previous solutions)
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
        //setDefaultObjective(problem, I, G, nSamples, nClones);

        // Define dummy product variable a_tij = fhat_tj * x_ij
        for(t = 0; t < nSamples; t++){
            for (i = 0; i < nClones; i++){
                for(Integer my_j : G.outEdges.get(i)){
                    linear = new Linear();
                    linear.add(1, "a_" + t + "_" + i + "_" + my_j);
                    linear.add(-1, "x_" + i + "_" + my_j);
                    problem.add(new Constraint("constraint" + counter++, linear, "<=", 0));

                    linear = new Linear();
                    linear.add(1, "a_" + t + "_" + i + "_" + my_j);
                    linear.add(-1, "fhat_" + t + "_" + my_j);
                    problem.add(new Constraint("constraint" + counter++, linear, "<=", 0));

                    linear = new Linear();
                    linear.add(1, "a_" + t + "_" + i + "_" + my_j);
                    linear.add(-1, "fhat_" + t + "_" + my_j);
                    linear.add(-1, "x_" + i + "_" + my_j);
                    problem.add(new Constraint("constraint" + counter++, linear, ">=", -1));

                    linear = new Linear();
                    linear.add(1, "a_" + t + "_" + i + "_" + my_j);
                    problem.add(new Constraint("constraint" + counter++, linear, ">=", 0));
                }
            }
        }


        // Define dummy product variable b_ti = w_i * u_ti
        for(t = 0; t < nSamples; t++){
            for (i = 0; i < nClones; i++){
                linear = new Linear();
                linear.add(1, "b_" + t + "_" + i);
                linear.add(-1, "w_" + i);
                problem.add(new Constraint("constraint" + counter++, linear, "<=", 0));

                linear = new Linear();
                linear.add(1, "b_" + t + "_" + i);
                linear.add(-1, "u_" + t + "_" + i);
                problem.add(new Constraint("constraint" + counter++, linear, "<=", 0));

                linear = new Linear();
                linear.add(1, "b_" + t + "_" + i);
                linear.add(-1, "u_" + t + "_" + i);
                linear.add(-1, "w_" + i);
                problem.add(new Constraint("constraint" + counter++, linear, ">=", -1));

                linear = new Linear();
                linear.add(1, "b_" + t + "_" + i);
                problem.add(new Constraint("constraint" + counter++, linear, ">=", 0));
            }
        }

        // Define dummy product variable c_ij = x_ij * tmin_j
        for (i = 0; i < nClones + 1; i++){
            for(Integer my_j : G.outEdges.get(i)){
                linear = new Linear();
                linear.add(1, "c_" + i + "_" + my_j);
                linear.add(-1 * nSamples, "x_" + i + "_" + my_j);
                problem.add(new Constraint("constraint" + counter++, linear, "<=", 0));

                linear = new Linear();
                linear.add(1, "c_" + i + "_" + my_j);
                linear.add(-1, "tmin_" + my_j);
                problem.add(new Constraint("constraint" + counter++, linear, "<=", 0));

                linear = new Linear();
                linear.add(1, "c_" + i + "_" + my_j);
                linear.add(-1, "tmin_" + my_j);
                linear.add(-1 * nClones, "x_" + i + "_" + my_j);
                problem.add(new Constraint("constraint" + counter++, linear, ">=", -1 * nClones));

                linear = new Linear();
                linear.add(1, "c_" + i + "_" + my_j);
                problem.add(new Constraint("constraint" + counter++, linear, ">=", 0));
            }
        }

        // Add depth variables d_i in [1, n] for each i
        for(i = 0; i < nClones; i++){
            problem.setVarType("d_" + i, Integer.class);
            linear = new Linear();
            linear.add(1, "d_" + i);
            problem.add(new Constraint("constraint" + counter++, linear, ">=", 1));

            linear = new Linear();
            linear.add(1, "d_" + i);
            problem.add(new Constraint("constraint" + counter++, linear, "<=", nClones));
        }
        // Set up depth of root
        problem.setVarType("d_" + nClones, Integer.class);
        linear = new Linear();
        linear.add(1, "d_" + nClones);
        problem.add(new Constraint("constraint" + counter++, linear, "=", 0));

        // Add depth constraints to prevent cycles: x_ij = 1 -> d_i = d_j - 1
        for(i = 0; i < nClones + 1; i++){
            for(Integer my_j : G.outEdges.get(i)){
                linear = new Linear();
                linear.add(1, "d_" + i);
                linear.add(nClones + 1, "x_" + i + "_" + my_j);
                linear.add(-1, "d_" + my_j);
                problem.add(new Constraint("constraint" + counter++, linear, "<=", nClones));

                linear = new Linear();
                linear.add(1, "d_" + my_j);
                linear.add(-1, "d_" + i);
                linear.add(nClones + 1, "x_" + i + "_" + my_j);
                problem.add(new Constraint("constraint" + counter++, linear, "<=", nClones + 2));
            }
        }

        // Add a binary w_i variable for each vertex i where w_i indicates that i is in the solution
        for(i = 0; i < nClones; i++){
            problem.setVarType("w_" + i, VarType.BOOL);
            linear = new Linear();
            linear.add(1, "w_" + i);
            for(Integer my_j:  G.inEdges.get(i)){
                linear.add(-1, "x_" + my_j + "_" + i);
            }
            problem.add(new Constraint("constraint" + counter++, linear, "=", 0));
        }

        // The root (index nClones) must have exactly 1 outgoing edge
        linear = new Linear();
        for(j = 0; j < nClones; j++){ // using all clones instead of delta_r because they are equivalent
            linear.add(1, "x_" + nClones + "_" + j);
        }
        problem.add(new Constraint("constraint" + counter++, linear, "=", 1));


        // Each NON-ROOT vertex must have at least as many incoming edges as outgoing
        for(j = 0; j < nClones; j++){
            for(Integer my_k : G.outEdges.get(j)){
                linear = new Linear();
                linear.add(-1, "x_" + j + "_" + my_k);
                for(Integer my_i : G.inEdges.get(j)) {
                    linear.add(1, "x_" + my_i + "_" + j);
                }
                problem.add(new Constraint("constraint" + counter++, linear, ">=", 0));
            }
        }

        // Each vertex can have at most 1 incoming edge
        for(j = 0; j < nClones; j++){
            linear = new Linear();
            for(Integer my_i : G.inEdges.get(j)){
                linear.add(1, "x_" + my_i + "_" + j);
            }
            problem.add(new Constraint("constraint" + counter++, linear, "<=", 1));
        }


        // Add constraints for fhat and u
        for(t = 0; t < nSamples; t++){
            // for each sample, sum across clones of u must be <= 1
            linear = new Linear();
            for(i = 0; i < nClones; i++){
                linear.add(1, "b_" + t + "_" + i);
            }
            problem.add(new Constraint("constraint" + counter++, linear, "<=", 1));

            for(i = 0; i < nClones; i++){
                //fhat is bounded by intervals
                linear = new Linear();
                linear.add(1, "fhat_" + t + "_" + i);
                assert I.intervals[0][t][i] < I.intervals[1][t][i];
                assert I.intervals[0][t][i] >= 0;
                problem.add(new Constraint("constraint" + counter++, linear, ">=", new BigDecimal(I.intervals[0][t][i])));
                problem.add(new Constraint("constraint" + counter++, linear, "<=", I.intervals[1][t][i]));

                //u must be non-negative
                linear = new Linear();
                linear.add(1, "u_" + t + "_" + i);
                problem.add(new Constraint("constraint" + counter++, linear, ">=", 0));

                // Add sum condition equivalence
                linear = new Linear();
                linear.add(1, "u_" + t + "_" + i);
                linear.add(-1, "fhat_" + t + "_" + i);
                for(Integer my_j : G.outEdges.get(i)){
                    linear.add(1, "a_" + t + "_" + i + "_" + my_j);
                }
                problem.add(new Constraint("constraint" + counter++, linear, "=", 0));

                if(longitudinal) {
                    // Define constraints using y
                    // t < tmin -> y = 0
                    linear = new Linear();
                    linear.add(1, "tmin_" + i);
                    linear.add(nSamples, "y_" + t + "_" + i);
                    problem.add(new Constraint("constraint" + counter++, linear, "<=", t + nSamples));

                    // t >= tmin -> y = 1
                    linear = new Linear();
                    linear.add(1, "tmin_" + i);
                    linear.add(nSamples, "y_" + t + "_" + i);
                    problem.add(new Constraint("constraint" + counter++, linear, ">=", t + 1));

                    // y = 0 -> fhat = 0
                    linear = new Linear();
                    linear.add(1, "fhat_" + t + "_" + i);
                    linear.add(-1, "y_" + t + "_" + i);
                    problem.add(new Constraint("constraint" + counter++, linear, "<=", 0));

                    //Define constraints using z
                    // t >= tmax -> z = 1
                    linear = new Linear();
                    linear.add(1, "tmax_" + i);
                    linear.add(nSamples, "z_" + t + "_" + i);
                    problem.add(new Constraint("constraint" + counter++, linear, ">=", t + 1));

                    // t < tmax -> z = 0
                    linear = new Linear();
                    linear.add(1, "tmax_" + i);
                    linear.add(nSamples, "z_" + t + "_" + i);
                    problem.add(new Constraint("constraint" + counter++, linear, "<=", t + nSamples));

                    // z = 1 -> u = 0
                    linear = new Linear();
                    linear.add(1, "u_" + t + "_" + i);
                    linear.add(1, "z_" + t + "_" + i);
                    problem.add(new Constraint("constraint" + counter++, linear, "<=", 1));

                    // u >= h for tmin <= t < tmax (w_i term allows u to vary for vertices that are not in the solution)
                    linear = new Linear();
                    linear.add(1, "u_" + t + "_" + i);
                    linear.add(-1, "y_" + t + "_" + i);
                    linear.add(1, "z_" + t + "_" + i);
                    linear.add(-1, "w_" + i);
                    problem.add(new Constraint("constraint" + counter++, linear, ">=", Main.MINIMUM_USAGE_H - 2));
                    //problem.add(new Constraint(linear, ">=", .001 - 1));
                }
            }
        }

        if (longitudinal){
            // lineage continuity: x_ij * tmin_j <= tmax_i
            for(i = 0; i < nClones; i++){
                for(Integer my_j : G.outEdges.get(i)){
                    linear = new Linear();
                    linear.add(1, "c_" + i + "_" + my_j);
                    linear.add(-1, "tmax_" + i);
                    problem.add(new Constraint("constraint" + counter++, linear, "<=", 0));
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
            problem.add(new Constraint("constraint" + counter++, linear, "<=", nSamples));
            linear = new Linear();
            linear.add(1, "tmin_" + i);
            problem.add(new Constraint("constraint" + counter++, linear, ">=", 0));

            // tmax is in [0, m + 1]
            linear = new Linear();
            linear.add(1, "tmax_" + i);
            problem.add(new Constraint("constraint" + counter++, linear, "<=", nSamples + 1));
            linear = new Linear();
            linear.add(1, "tmax_" + i);
            problem.add(new Constraint("constraint" + counter++, linear, ">=", 0));

            // tmin must be at least as small as tmax
            linear = new Linear();
            linear.add(1, "tmin_" + i);
            linear.add(-1, "tmax_" + i);
            problem.add(new Constraint("constraint" + counter++, linear, "<=", 0));

        }
        for(i = 0; i < nClones; i++){
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

        //System.out.println(problem.getVariables());
        //System.out.println(problem.getConstraints());

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
                linear.add(Math.pow(10, Main.PRECISION_DIGITS + 1), "x_" + i + "_" + my_j);
            }
        }

        // Add a term for each frequency matrix entry present in sample (prioritize largest CCF mutations)
        for(t = 0; t < nSamples; t++) {
            for (i = 0; i < nClones; i++) {
                linear.add(Math.pow(10, Main.PRECISION_DIGITS) * I.intervals[0][t][i] / ((double) nClones * nSamples), "w_" + i);
            }
        }

        // Add a term for each (included) usage matrix entry
        for(t = 0; t < nSamples; t++) {
            for (i = 0; i < nClones; i++) {
                linear.add(((double) -1) / ((double) nClones * nSamples), "b_" + t + "_" + i);
            }
        }
        p.setObjective(linear, OptType.MAX);
    }

    private static void setL0Norm(Problem problem, Instance I, Graph G, int nSamples, int nClones){
        int t, i, j, counter = 0;
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
                problem.add(new Constraint("l0constraint" + counter++, linear, "<=" , 0));

                linear = new Linear();
                linear.add(1, "q_" + t + "_" + j);
                linear.add(1, "u_" + t + "_" + j);
                problem.add(new Constraint("l0constraint" + counter++,linear, "<=", 1));
            }
        }

        // Add a term for each possible edge in the tree (most mutations)
        linear = new Linear();
        for(i = 0; i < nClones + 1; i++){
            for(Integer my_j : G.outEdges.get(i)) {
                linear.add(Math.pow(10, Main.PRECISION_DIGITS + 1), "x_" + i + "_" + my_j);
            }
        }

        // Add a term for each frequency matrix entry present in sample (prioritize largest CCF mutations)
        for(t = 0; t < nSamples; t++) {
            for (i = 0; i < nClones; i++) {
                linear.add(Math.pow(10, Main.PRECISION_DIGITS) * I.intervals[0][t][i] /
                        ((double) nClones * nSamples), "w_" + i);
            }
        }

        // Add a term for each 0 value in U
        for(t = 0; t < nSamples; t++) {
            for (i = 0; i < nClones; i++) {
                linear.add(((double) 1) / ((double) nClones * nSamples), "q_" + t + "_" + i);
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
    public static Constraint constructDummyConstraint(ILPResult r, int counter){
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
        return new Constraint("extraConstraint" + counter, linear, "<=", total - 1);
    }
}
