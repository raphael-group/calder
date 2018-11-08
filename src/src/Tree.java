import java.util.*;

public class Tree extends Graph {
    public final int root;
    public final int nSamples;
    private LinkedList<Integer> label;
    private int valid = 0;
    HashMap<Integer, double[]> newMins;
    HashMap<Integer, Integer> idToIndex;
    HashMap<Integer, Integer> indexToId;


    public Tree(Graph T, Integer root, int nSamples){
        super(T);
        this.root = root;
        this.nSamples = nSamples;
    }
    

    public boolean equals(Object o){
        if (o instanceof  Tree){
            return ((Tree) o).getRootLabel().equals(this.getRootLabel());
        } else {
            return false;
        }
    }

    // Returns true if there is an assignment of fhat in [minf, maxf] for every mutation and every sample s.t. the sum
    // condition is met for all edges in the tree
    public boolean isValid(){
        if(valid == 0){
            if(checkValid()){
                valid = 1;
                return true;
            } else {
                valid = -1;
                return false;
            }
        } else {
            return valid > 0;
        }
    }

    private boolean checkValid(){
        newMins = new HashMap<>();
        // Compute the new min for each node (populating newMins)
        validHelper(root);

        // Verify that for each parent v, for each sample t, v.maxf[t] > newmin[t]
        int t;
        double[] vMaxes, myMins;
        for(Integer v : vertices){
            vMaxes = vertexData.get(v).maxf;
            myMins = newMins.get(v);
            for(t = 0; t < nSamples; t++){
                //System.out.println(vMaxes[t] + " " + myMins[t]);
                if(vMaxes[t] < myMins[t]){
                    //System.err.println("Vertex " + v + " failed as parent of " + outEdges.get(v) + " : " + vMaxes[t] + " < " + myMins[t]);
                    return false;
                } else {

                }
            }
        }
        /*
        System.out.println(this);
        for(Integer v : newMins.keySet()){
            System.out.println(v + ": " + Arrays.toString(newMins.get(v)));
        }
        */

        return true;
    }

    public double getSum(){
        double total = 0;
        VertexData vData;
        for (Integer v : vertices){
            vData = Graph.vertexData.get(v);
            for(int t = 0; t < nSamples; t++){
                total += vData.minf[t];
            }
        }
        return total;
    }

    // Computes newMin for each vertex in the subtree rooted at v
    private void validHelper(int v){
        for(Integer child : outEdges.get(v)){
            // if child isn't in the hashmap, call validHelper on the child
            if(!newMins.containsKey(child)){
                validHelper(child);
            }
        }
        double[] myMins = new double[nSamples];
        double[] myLb = vertexData.get(v).minf;
        double[] childMins;
        for(Integer child : outEdges.get(v)){
            childMins = newMins.get(child);
            for(int t = 0; t < nSamples; t++){
                myMins[t] += childMins[t];
            }
        }
        for(int t = 0; t < nSamples; t++){
            if(myLb[t] > myMins[t]){
                myMins[t] = myLb[t];
            }
        }
        newMins.put(v, myMins);
    }


    // Retrieves the label for this tree (memoized)
    public LinkedList<Integer> getRootLabel(){
        if(this.label != null){
            return this.label;
        } else {
            this.label = getLabel(root);
            return this.label;
        }
    }

    // Computes the label for this tree to be compared to other trees (Aho, Hopcroft and Ullman algorithm)
    private LinkedList<Integer> getLabel(Integer r) {

        HashSet<Integer> children = outEdges.get(r);
        assert children != null;

        LinkedList<Integer> label;
        LinkedList<LinkedList<Integer>> childrenLabels;
        if (children.size() == 0) {
            // If I'm a leaf, return 01 as my label
            label = new LinkedList<>();
            label.add(r);
            label.add(r);
            return label;
        } else {
            // First, recursively get the labels of all my children
            childrenLabels = new LinkedList<>();
            for (Integer child : children) {
                childrenLabels.addLast(getLabel(child));
            }
            childrenLabels.sort(new LabelComparator());

            // Prepend id for my label
            label = new LinkedList<>();
            label.addLast(r);

            // Add sorted childrens' labels in order
            for (LinkedList<Integer> l : childrenLabels) {
                label.addAll(l);
            }
            // Append my id to my label
            label.addLast(r);
            return label;
        }
    }
    public double[][] minU(){
        assert isValid();
        double[][] myU = new double[nSamples][vertices.size()];

        //TODO: fill using mins
        double[] minf;
        double[] inducedMins;
        int i = 0;
        idToIndex = new HashMap<>();
        indexToId = new HashMap<>();
        for(int v : vertices){
            idToIndex.put(v, i);
            indexToId.put(i, v);
            minf = Graph.vertexData.get(v).minf;
            inducedMins = newMins.get(v);
            for(int t = 0; t < nSamples; t++){
                if(inducedMins[t] > minf[t]){
                    // My children take all of my usage (they induce my min to be higher than my min freq)
                    myU[t][i] = 0;
                } else {
                    myU[t][i] = minf[t];
                    for(Integer child : outEdges.get(v)){
                        // Each vertex has its min allowed frequency, so use this to compute sum condition
                        myU[t][i] -= newMins.get(child)[t];
                    }
                }

            }
            i++;
        }

        return myU;
    }

    private class LabelComparator implements Comparator<LinkedList<Integer>> {
        public int compare(LinkedList<Integer> l1, LinkedList<Integer> l2) {
            if (l1 == null) {
                return -1;
            } else if (l2 == null) {
                return 1;
            } else if (l1.size() > l2.size()) {
                return 1;
            } else if (l1.size() < l2.size()) {
                return -1;
            } else {
                int a, b;
                while (!l1.isEmpty()) {
                    a = l1.poll();
                    b = l2.poll();
                    if (a != b) {
                        return a - b;
                    }
                }
                return 0;
            }
        }
    }
    public String vertexEncoding(){
        int prev = -1;
        StringBuilder sb = new StringBuilder();
        for(Integer v : vertices){
            sb.append(v);
            sb.append('_');
        }
        return sb.toString();
    }

    public String toString(){
        if(idToIndex == null){
            return super.toString();
        } else {
            StringBuilder b = new StringBuilder();
            for(Integer v: vertices){
                b.append(idToIndex.get(v));
                b.append(": ");
                for(Integer dest: outEdges.get(v)){
                    b.append('(');
                    b.append(idToIndex.get(v));
                    b.append(',');
                    b.append(idToIndex.get(dest));
                    b.append(") ");
                }
                b.append('\n');
            }
            return b.toString();
        }
    }

}
