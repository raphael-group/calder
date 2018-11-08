import java.util.HashMap;
import java.util.LinkedList;

public class ILPResultCollector {
    private double bestScore = Double.MAX_VALUE;
    private LinkedList<ILPResult> bestResults = new LinkedList<>();
    private LinkedList<ILPResult> allResults = new LinkedList<>();
    private LinkedList<ILPResult> noHResults = new LinkedList<>();
    private HashMap<String, LinkedList<ILPResult>> feasible = new HashMap<>();
    private HashMap<String, LinkedList<Tree>> infeasible = new HashMap<>();
    private final Instance I;

    private int numFeasible = 0;

    public ILPResultCollector(Instance I){
        this.I = I;
    }

    public synchronized void addResult(ILPResult res){
        assert res != null;

        Tree T = res.T;
        allResults.add(res);

        numFeasible++;
        double myObj = res.objective;

        if (myObj < bestScore) {
            bestResults = new LinkedList<>();
            bestResults.add(res);
            bestScore = myObj;
        } else if (myObj == bestScore){
            bestResults.add(res);
        }
        if(feasible.containsKey(T.vertexEncoding())){
            feasible.get(T.vertexEncoding()).add(res);
        } else {
            LinkedList<ILPResult> results = new LinkedList<>();
            results.add(res);
            feasible.put(T.vertexEncoding(), results);
        }
    }

    public synchronized void addBadResult(Tree T){
        if(infeasible.containsKey(T.vertexEncoding())){
            infeasible.get(T.vertexEncoding()).add(T);
        } else {
            LinkedList<Tree> trees = new LinkedList<>();
            trees.add(T);
            infeasible.put(T.vertexEncoding(), trees);
        }
    }

    public synchronized void addNoHResult(ILPResult res){
        noHResults.add(res);
    }

    public double getBestScore(){
        return bestScore;
    }
    public int getNumFeasible(){
        return numFeasible;
    }
    public int getBestResultsSize(){
        return bestResults.size();
    }
    public String printBestResults(){
        StringBuilder b = new StringBuilder();
        for(ILPResult res : bestResults){
            b.append(res.toString());
            b.append('\n');
        }
        b.append('\n');

        return b.toString();
    }
    public String printAllResults(){
        StringBuilder b = new StringBuilder();
        for(ILPResult res : allResults){
            b.append(res.toString());
            b.append('\n');
        }
        b.append('\n');

        return b.toString();
    }
    public LinkedList<ILPResult> getBestResults(){
        return bestResults;
    }

    public void printClonePresence(){
        int[] totalClonePresence = new int[I.colLabels.length * I.rowLabels.length];
        for(ILPResult res : bestResults){
            totalClonePresence[res.totalClonePresence()]++;
        }
        int mymin = 0;
        int numAtMin = 0;
        for(int i = 0; i < totalClonePresence.length; i++){
            if(totalClonePresence[i] > 0){
                mymin = i;
                numAtMin = totalClonePresence[i];
                break;
            }
        }
        int numOptimaMin0s = 0;
        for(ILPResult res : bestResults){
            int mycount = 0;
            for(int i = 0; i < res.uhat.length; i++){
                for(int j = 0; j < res.uhat[0].length; j++){
                    if (res.uhat[i][j] > 0){
                        mycount++;
                    }
                }
            }
            if(mycount == mymin){
                numOptimaMin0s ++;
            }
        }

        System.out.println("Minimum number of non-0 entries: " + mymin);
        System.out.println("Found " + numAtMin + " solutions with this many entries, including " + numOptimaMin0s + " optima");

        //System.out.println(Arrays.toString(totalClonePresence));
    }

    public void printPairs(){
        for(String encoding : infeasible.keySet()){
            if(feasible.containsKey(encoding)){
                int nFeas = feasible.get(encoding).size();
                int nInfeas = infeasible.get(encoding).size();
                System.out.println("Found " + nFeas + " feasible trees and " + nInfeas +
                        " infeasible trees on vertices " + encoding);

                if (encoding.equals("2_3_4_5_6_7_8_13_14_")) {
                    System.out.println("--------------------------------------------------");
                    for (ILPResult res : feasible.get(encoding)) {
                        System.out.println(res.toString());
                    }
                    System.out.println("Infeasible results:----------------------------------------------------");
                    for (Tree T : infeasible.get(encoding)) {
                        /*
                        double[][] U = T.minU();
                        String[] myColLabels = new String[T.vertices.size()];
                        //System.out.println(Arrays.toString(i.colLabels));
                        //System.out.println(T.indexToId);
                        for (int p = 0; p < myColLabels.length; p++) {
                            myColLabels[p] = I.colLabels[T.indexToId.get(p)];
                        }
                        System.out.println(Calder.print2DArray(U, myColLabels, I.rowLabels));
                        System.out.println(T);
                        */
                        System.out.println(Calder.inferUwithoutH(I, T));
                    }
                }
            }
        }
    }
    public void printNoHResults(){
        System.out.println("Found results with no H: " + noHResults.size());
        for(ILPResult res : noHResults){
            System.out.println(res);
        }
    }
}
