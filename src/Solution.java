import java.util.*;

public class Solution {
    private final int nVertices;
    public final Tree T;

    private int[][] B;
    // Labels for BOTH the ROWS (taxa/clones) and COLUMNS (mutations) of B
    private int[] BLabels;

    private String[][] stringU;
    private int[][] U;
    private int[] UColSizes;
    private int isLongitudinal = 0;

    private final int factor;
    private final int nSamples;

    private HashMap<Integer, Integer> idToIdx;
    private HashMap<Integer, Integer> idxToId;

    public static boolean verbose = false;

    public Solution(Tree T, int nSamples, int factor){
        this.factor = factor;
        this.nSamples = nSamples;

        nVertices = T.vertexData.size();
        assert nVertices == T.outEdges.keySet().size();
        assert nVertices == T.vertices.size();

        this.T = copyTree(T);
        B = generateB();
        stringU = new String[nSamples][T.vertices.size()];
        //U = generateU();
    }

    public String BtoString(){
        return Calder.print2DArray(B, BLabels, BLabels);
    }

    public String UtoString(){
        return Calder.print2DArray(U, BLabels);
    }

    public String toString(){
        // TODO: implement using above 2 functions
        return "";
    }

    private Tree copyTree(Tree T){
        Graph G = new Graph(T);
        return new Tree(G, T.root, T.nSamples);
    }

    // Builds the B matrix that equivalently represents the tree T
    private int[][] generateB(){
        int[][] B = new int[nVertices][nVertices];

        int r = T.root;

        //do DFS on T to generate B
        LinkedList<Integer> next = new LinkedList<>();
        HashSet<Integer> visited = new HashSet<>();

        HashMap<Integer, Integer> childToParent = new HashMap<>();
        idToIdx = new HashMap<>();
        idxToId = new HashMap<>();
        BLabels = new int[nVertices];

        next.push(r);
        visited.add(r);

        Integer curr, parent;
        String myName;
        int i = 0, parentIdx;
        while(!next.isEmpty()){
            curr = next.pop();

            int[] myRow = new int[nVertices];
            idToIdx.put(curr, i);
            idxToId.put(i, curr);
            BLabels[i] = curr;

            // Add diagonal element to my row
            myRow[i] = 1;

            // Populate row with my parent's mutations
            parent = childToParent.get(curr);
            if (parent != null){
                parentIdx = idToIdx.get(parent);
                for (int j = 0; j < nVertices; j++){
                    myRow[j] = Math.max(myRow[j], B[parentIdx][j]);
                }
            }

            //TODO: maybe do something with the name I'm building here in the output - use as label?
            // Builds clone name by aggregating mutations up to the root
            /*
            myName = curr.toString();
            while(parent != null){
                myName = parent.toString() + myName;
                parent = childToParent.get(parent);
            }
            */

            for(Integer dest: T.outEdges.get(curr)){
                assert !visited.contains(dest);

                visited.add(dest);
                childToParent.put(dest, curr);
                next.push(dest);

            }
            B[i++] = myRow;
        }
        return B;
    }

    // Builds the U matrix determined by the tree T and frequency matrix F
    /*
    private int[][] generateU(){
        int[][] U = new int[nSamples][nVertices];

        // i indexes rows and j indexes columns in this method

        int myClone;
        UColSizes = new int[nVertices];

        int[][] columns = new int[nVertices][];

        int[] myCol;
        for(int j = 0; j < nVertices; j++){
            myClone = BLabels[j];
            myCol = Arrays.copyOf(T.vertexData.get(myClone).freq, nSamples);

            for(Integer child: T.outEdges.get(myClone)){
                for (int i = 0; i < nSamples; i++){
                    myCol[i] -= T.vertexData.get(child).freq[i];
                }
            }
            columns[j] = myCol;
        }

        int intSlack;
        for(int i = 0; i < nSamples; i++){
            for(int j = 0; j < nVertices; j++){
                intSlack = columns[j][i];
                U[i][j] = intSlack;
                stringU[i][j] = restoreDecimal(intSlack, factor);
            }
        }

        return U;
    }
    */

    // Convenience method to rescale an (internal representation) integer to its actual double value
    private String restoreDecimal(int value, int factor){
        StringBuilder sb = new StringBuilder();
        sb.append("0.");
        while(factor > 10){
            factor /= 10;
            sb.append("0");
        }
        sb.append(value);
        return sb.toString();
    }

    // Retrieves the longitudinality of this solution (memoized)
    public boolean isLongitudinal(){
        if (isLongitudinal == 0) {
            //
            /*if (checkLongitudinal()) { //TODO: fix checkLongitudinal
                isLongitudinal = 1;
            } else {
                isLongitudinal = -1;
            }*/
        }
        return isLongitudinal > 0;
    }

    /*
    // Checks the sum condition and longitudinal conditions on each edge of T
    private boolean checkLongitudinal(){
        LinkedList<Integer> q = new LinkedList<>();
        q.add(T.root);

        int[] labelToIdx = new int[BLabels.length];
        for(int i = 0; i < BLabels.length; i++){
            labelToIdx[BLabels[i]] = i;
        }

        int curr, tmax, tmin, t, childBorn, myIdx;
        while(!q.isEmpty()) {
            curr = q.pop();

            // Columns in U are rearranged according to the sorting of B, need to map index back to column
            myIdx = labelToIdx[curr];

            tmin = T.vertexData.get(curr).tmin;
            assert tmin < U.length;

            // Compute my tmax
            //tmax is the first sample after tmin in which U is 0
            tmax = tmin;
            while(tmax < U.length && U[tmax][myIdx] > 0){
                tmax++;
            }

            // Check that my extinction (at time tmax) was permanent
            for(t = tmax; t < U.length; t++){
                if(U[t][myIdx] < 0){
                    // Sum condition violated
                    return false;
                }

                if(U[t][myIdx] > 0){
                    if (verbose){
                        System.err.println("Permanent extinction violated for clone [" + + curr + "] in sample [" + t + "]; tmin=" + tmin + ", tmax=" + tmax);
                        System.err.println("U = " + U[t][myIdx]);
                    }
                    return false;
                }
            }

            // Check that my children could have been born while I was alive
            for(Integer child : T.outEdges.get(curr)){
                childBorn = T.vertexData.get(child).tmin;
                if(childBorn > tmax){
                    if(verbose){
                        System.err.println("Lineage continuity violated");
                    }
                    return false;
                } else if (childBorn < tmin){
                    System.out.println("this happened");
                }
                q.add(child);
            }

        }
        return true;
    }
    */

    // Two solutions are equal iff their corresponding trees are equal
    // NOTE: two Solutions are only comparable if they solved the same problem/frequency matrix F
    public boolean equals(Object o){
        if(o instanceof  Solution){
            Solution s = (Solution) o;
            //TODO: add some assertion that represents the assumption that two solutions are from the same problem/F
            return s.T.equals(T);
        } else {
            return false;
        }
    }

}
