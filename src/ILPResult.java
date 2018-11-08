import net.sf.javailp.Result;

import java.util.HashMap;

public class ILPResult {
    public final Result result;
    public final double objective;
    public final double[][] fhat;
    public final double[][] uhat;
    public final int nClones;
    public final int nSamples;
    public final String[] rowLabels;
    public final String[] colLabels;
    public final HashMap<Integer,Integer> indexToId;
    public final HashMap<Integer,Integer> idToIndex;
    public final Tree T;

    public ILPResult(Result res, int nClones, int nSamples, HashMap<Integer,Integer> indexToId,
                     HashMap<Integer,Integer> idToIndex, Tree T, String[] rowLabels, String[] colLabels){
        assert res != null;
        assert nClones > 0;
        assert nSamples > 0;

        this.result = res;
        this.objective = (double) res.getObjective();
        this.nClones = nClones;
        this.nSamples = nSamples;
        this.indexToId = indexToId;
        this.idToIndex = idToIndex;
        this.T = T;


        this.rowLabels = rowLabels;
        this.colLabels = colLabels;


        uhat = new double[nSamples][nClones];
        int i, t, id;
        for(t = 0; t < nSamples; t++){
            for(i = 0; i < nClones; i++){
                uhat[t][i] = (double) result.getPrimalValue("uhat_" + t + "_" + i);

            }
        }
        fhat = new double[nSamples][nClones];
        for(t = 0; t < nSamples; t++){
            for(i = 0; i < nClones; i++) {
                    id = indexToId.get(i);
                    fhat[t][i] = (double) result.getPrimalValue("fhat_" + t + "_" + i);

            }
        }
    }

    public int totalClonePresence(){
        int result = 0;
        for(int t = 0; t < nSamples; t++){
            for(int i = 0; i < nClones; i ++){
                if(uhat[t][i] > 0){
                    result++;
                }
            }
        }
        return result;
    }
    public int[] clonePrecense(){
        int[] result = new int[nSamples];
        for(int t = 0; t < nSamples; t++){
            for(int i = 0; i < nClones; i ++){
                if(uhat[t][i] > 0){
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
        b.append(Calder.print2DArray(fhat, myColLabels, myRowLabels));
        b.append('\n');

        b.append("U:\n");
        b.append(Calder.print2DArray(uhat, myColLabels, myRowLabels));
        b.append('\n');

        b.append("tmin:");
        int id;
        for(i = 0; i < nClones; i++) {
            b.append(' ');
            b.append(result.getPrimalValue("tmin_" + i));
        }

        b.append('\n');

        b.append("tmax:");
        for(i = 0; i < nClones; i++) {
            b.append(' ');
            b.append(result.getPrimalValue("tmax_" + i));
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


        b.append(Calder.print2DArray(uhat, myColLabels, myRowLabels));
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
