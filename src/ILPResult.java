import net.sf.javailp.Linear;
import net.sf.javailp.Result;

import java.lang.reflect.Array;
import java.math.BigDecimal;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.HashMap;

public class ILPResult {
    final Result result;
    final double objective;
    final double[][] fhat;
    final double[][] u;
    final int nClones;
    final int nSamples;
    final String[] rowLabels;
    final String[] colLabels;
    private final HashMap<Integer,Integer> indexToId;
    private final HashMap<Integer,Integer> idToIndex;
    final Tree T;
    final double term1;
    final double term2;
    final double term3;
    final double term4;


    public ILPResult(Result res, Instance I, Graph ancestryGraph){
        assert res != null;
        assert I.nMuts > 0;
        assert I.nSamples > 0;

        this.result = res;
        this.objective = (double) res.getObjective();
        this.nSamples = I.nSamples;

        //find the root by looking at all x_nClones_i
        int i, j, val, root = -1;
        boolean foundRoot = false;
        for(i = 0; i < I.nMuts; i++){
            val = (int) res.getPrimalValue("x_" + I.nMuts + "_" + i);
            if(val == 1){
                if(foundRoot){
                    System.out.println("Found multiple root vertices");
                } else {
                    foundRoot = true;
                    root = i;
                }
            }
        }

        assert root > -1; // assert that we have a root

        Tree T = new Tree(new Graph(), root, nSamples);
        T.addVertex(root);
        int clonesInTree = 1;

        // build tree from positive x values
        boolean[] vertexPresent = new boolean[I.nMuts];
        vertexPresent[root] = true;

        for(i = 0; i < I.nMuts; i++){
            for(Integer my_j : ancestryGraph.outEdges.get(i)){
                val = (int) res.getPrimalValue("x_" + i + "_" + my_j);
                if (val == 1){
                    T.addEdge(i, my_j);
                    clonesInTree += 1;
                    vertexPresent[my_j] = true;
                }
            }
        }
        assert T.inEdges.get(root).size() == 0; // the root should have no incoming vertices

        //initialize nClones with the actual size of the tree
        nClones = clonesInTree;

        // create indexToId and idToIndex maps: id is the original value in [0, nMuts), index is mapped to [0, nClones]
        HashMap<Integer, Integer> indexToId = new HashMap<>();
        HashMap<Integer, Integer> idToIndex = new HashMap<>();
        int idx = 0;
        for(i = 0; i < I.nMuts; i++){
            if(vertexPresent[i]){
                idToIndex.put(i, idx);
                indexToId.put(idx, i);
                idx++;
            }
        }

        this.indexToId = indexToId;
        this.idToIndex = idToIndex;
        this.T = T;


        this.rowLabels = I.rowLabels;
        this.colLabels = I.colLabels;

        u = new double[nSamples][nClones];
        int t, id;

        // u has columns for only those mutations that are in the tree (there are nClones of them)
        for(t = 0; t < nSamples; t++){
            for(i = 0; i < nClones; i++){
                id = indexToId.get(i);
                u[t][i] = (double) result.getPrimalValue("u_" + t + "_" + id);
            }
        }

        // fhat has columns for only those mutations that are in the tree (there are nClones of them)
        fhat = new double[nSamples][nClones];
        for(t = 0; t < nSamples; t++){
            for(i = 0; i < nClones; i++) {
                id = indexToId.get(i);
                fhat[t][i] = (double) result.getPrimalValue("fhat_" + t + "_" + id);

            }
        }

        /*
        // check to make sure that q_tj = 1 iff u_tj = 0 and j is in the tree
        int[][] Q = new int[nSamples][nClones];
        for(t = 0; t < nSamples; t++){
            for(i = 0; i < nClones; i++) {
                id = indexToId.get(i);
                Q[t][i] = (int) result.getPrimalValue("q_" + t + "_" + id);
                System.out.println(Q[t][i] + " " + u[t][i]);
            }
        }
        */


        // print objective function components
        double term1 = 0;
        // Add a term for each possible edge in the tree (most mutations)
        for(i = 0; i < I.nMuts + 1; i++){
            for(Integer my_j : ancestryGraph.outEdges.get(i)) {
                term1 += Math.pow(10, Main.PRECISION_DIGITS + 1) * (int) result.getPrimalValue("x_" + i + "_" + my_j);
            }
        }
        this.term1 = term1;

        double term2 = 0;
        // Add a term for each frequency matrix entry present in sample (prioritize largest CCF mutations)
        for(t = 0; t < nSamples; t++) {
            for (i = 0; i < I.nMuts; i++) {
                term2 += Math.pow(10, Main.PRECISION_DIGITS) * (I.intervals[0][t][i] / ((double) nClones * nSamples)) * (int) result.getPrimalValue("w_" + i);
            }
        }
        this.term2 = term2;

        double term3 = 0;
        // Add a term for each 0 value in U
        for(t = 0; t < nSamples; t++) {
            for (i = 0; i < I.nMuts; i++) {
                term3 += (1.0 / ((double) nClones * nSamples)) * (int) result.getPrimalValue("q_" + t + "_" + i);
            }
        }
        this.term3 = term3;

        double term4 = 0;
        if(Main.OBJECTIVE == Main.Objective.L0center){
            for(t = 0; t < nSamples; t++) {
                for (i = 0; i < nClones; i++) {
                    // Add difference term to objective
                    term4 += -1.0 / ((double) nClones * nSamples * nClones * nSamples) * (double) result.getPrimalValue("d_" + t + "_" + i);
                }
            }
        }
        this.term4 = term4;
    }

    public int totalClonePresence(){
        int result = 0;
        for(int t = 0; t < nSamples; t++){
            for(int i = 0; i < nClones; i ++){
                if(u[t][i] > 0){
                    result++;
                }
            }
        }
        return result;
    }
    public double totalPurity(){
        double result = 0;
        for(int t = 0; t < nSamples; t++){
            for(int i = 0; i < nClones; i ++){
                if(u[t][i] > 0){
                    result+= u[t][i];
                }
            }
        }
        return result;
    }

    public int[] clonePrecense(){
        int[] result = new int[nSamples];
        for(int t = 0; t < nSamples; t++){
            for(int i = 0; i < nClones; i ++){
                if(u[t][i] > 0){
                    result[t]++;
                }
            }
        }
        return result;
    }

    public String toString(){
        return prettyPrint();
    }


    public String verbosePrint(){
        StringBuilder b = new StringBuilder();

        int t, i;

        b.append("Objective function value: ");
        b.append(objective);
        b.append('\n');


        String[] myColLabels = new String[nClones];
        for(i = 0; i < nClones; i++) {
            if(colLabels.length > 1){
                myColLabels[i] = colLabels[indexToId.get(i)];
            } else {
                myColLabels[i] = indexToId.get(i) + "";
            }
        }

        String[] myRowLabels = new String[nSamples];
        for(i = 0; i < nSamples; i++){
            if(rowLabels.length > 1){
                myRowLabels[i] = rowLabels[i];
            } else {
                myRowLabels[i] = i + "";
            }
        }

        b.append("Fhat:\n");
        b.append(Util.print2DArray(fhat, myColLabels, myRowLabels));
        b.append('\n');

        b.append("U:\n");
        b.append(Util.print2DArray(u, myColLabels, myRowLabels));
        b.append('\n');

        b.append("tmin:");
        int id;
        for(i = 0; i < nClones; i++) {
            b.append(' ');
            id = indexToId.get(i);
            b.append(result.getPrimalValue("tmin_" + id));
        }

        b.append('\n');

        b.append("tmax:");
        for(i = 0; i < nClones; i++) {
            b.append(' ');
            id = indexToId.get(i);
            b.append(result.getPrimalValue("tmax_" + id));
        }
        b.append('\n');
        b.append('\n');

        Graph G = new Graph();
        int root = idToIndex.get(T.root);
        for(Integer u : T.vertices){
            for (Integer v : T.outEdges.get(u)){
                G.addEdge(idToIndex.get(u), idToIndex.get(v));
            }
        }

        b.append("Tree:\n");
        b.append(new Tree(G, root, nSamples));

        return b.toString();
    }

    public String treeToString(String name){
        return T.toDot(name, colLabels);
    }

    public String treeToString(){
        return treeToString("solution_tree");
    }

    public String solnToString(){
        StringBuilder b = new StringBuilder();
        DecimalFormat df = new DecimalFormat();
        df.setMaximumFractionDigits(Main.PRECISION_DIGITS);

        assert u.length == nSamples;
        assert u[0].length == nClones;

        // header row for fhat
        b.append("Fhat");
        b.append(System.lineSeparator());
        b.append("samples");
        int t, i;
        for(i = 0; i < nClones; i++){
            b.append(',');
            b.append(colLabels[indexToId.get(i)]);
        }
        b.append(System.lineSeparator());

        for(t = 0; t < nSamples; t++){
            b.append(rowLabels[t]);
            for(i = 0; i < nClones; i++){
                b.append(',');
                b.append(df.format(fhat[t][i]));
            }
            b.append(System.lineSeparator());
        }


        b.append(System.lineSeparator());

        b.append("U");
        b.append(System.lineSeparator());
        b.append("samples");
        for(i = 0; i < nClones; i++){
            b.append(',');
            b.append(colLabels[indexToId.get(i)]);
        }
        b.append(System.lineSeparator());

        for(t = 0; t < nSamples; t++){
            b.append(rowLabels[t]);
            for(i = 0; i < nClones; i++){
                b.append(',');
                b.append(df.format(u[t][i]));
            }
            b.append(System.lineSeparator());
        }

        return b.toString();

    }

    public String prettyPrint(){
        StringBuilder b = new StringBuilder();

        int t, i;
        String[] myColLabels = new String[nClones];
        for(i = 0; i < nClones; i++) {
            if(colLabels.length > 1){
                myColLabels[i] = colLabels[indexToId.get(i)];
            } else {
                myColLabels[i] = indexToId.get(i) + "";
            }
        }

        String[] myRowLabels = new String[nSamples];
        for(i = 0; i < nSamples; i++){
            if(rowLabels.length > 1){
                myRowLabels[i] = rowLabels[i];
            } else {
                myRowLabels[i] = i + "";
            }
        }

        System.out.println("Objective: ");
        System.out.println(objective);
        System.out.println("Terms: ");
        System.out.println(term1);
        System.out.println(term2);
        System.out.println(term3);
        System.out.println(term4);

        b.append("Clone proportion matrix U");
        b.append('\n');
        b.append(Util.print2DArray(u, myColLabels, myRowLabels, Main.PRECISION_DIGITS));
        b.append('\n');

        b.append("Tree root vertex: ");
        b.append(colLabels[T.root]);
        b.append("\n");
        b.append("Tree edges:");
        b.append('\n');
        boolean first = true;
        for(Integer u : T.vertices){
            for (Integer v : T.outEdges.get(u)){
                if(first){
                    first = false;
                } else {
                    b.append('\n');
                }
                b.append(colLabels[u]);
                b.append(" -> ");
                b.append(colLabels[v]);
            }
        }
        b.append('\n');

        return b.toString();
    }
}
