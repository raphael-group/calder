import net.sf.javailp.Result;

import java.math.BigDecimal;
import java.util.HashMap;

public class ILPResult {
    final Result result;
    final Number objective;
    final double[][] fhat;
    final double[][] u;
    final int nClones;
    final int nSamples;
    final String[] rowLabels;
    final String[] colLabels;
    private final HashMap<Integer,Integer> indexToId;
    private final HashMap<Integer,Integer> idToIndex;
    final Tree T;

    public ILPResult(Result res, Instance I, Graph ancestryGraph){
        assert res != null;
        assert I.nMuts > 0;
        assert I.nSamples > 0;

        this.result = res;
        this.objective = res.getObjective();
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
        StringBuilder b = new StringBuilder();

        int t, i;
        /*
        double distance = 0;
        for(t = 0; t < nSamples; t++){
            for(i = 0; i < nClones; i++){
                distance += Math.abs(fbars[t][i] - (double) result.getPrimalValue("fhat_" + t + "_" + i));
            }
        }
        */
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
    public String printConcise(){
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


        b.append(Util.print2DArray(u, myColLabels, myRowLabels));
        b.append('\n');

        Graph G = new Graph();
        int root = idToIndex.get(T.root);
        for(Integer u : T.vertices){
            for (Integer v : T.outEdges.get(u)){
                G.addEdge(idToIndex.get(u), idToIndex.get(v));
            }
        }
        b.append(new Tree(G, root, nSamples));

        return b.toString();
    }
}
