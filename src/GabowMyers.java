import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;

public abstract class GabowMyers {
    protected LinkedList<Tree> results;
    protected Graph myG;
    protected Tree L;

    protected int nSamples;

    protected int doneSize;


    public LinkedList<Tree> run(){
        return run(myG.vertices.size());
    }

    public LinkedList<Tree> run(int doneSize){
        // Make a copy of G
        Graph G;
        Tree T;
        LinkedList<Edge> frontier;

        this.doneSize = doneSize;

        LinkedList<Tree> valid;
        Collection<Integer> rootCandidates;
        do {
            results = new LinkedList<>();

            System.out.println("Looking for trees with " + this.doneSize + " vertices");

            if(this.doneSize == Graph.vertexData.size()){
                rootCandidates = myG.findRootCandidates();
            } else {
                rootCandidates = myG.vertices;
            }
            for(Integer r : rootCandidates){
                Calder.resetCounts();

                // Make a copy of G for each iteration to avoid unexpected behavior
                G = new Graph(myG);

                // Initialize T with just the root
                T = new Tree(myG.newEmptyGraph(), r, nSamples);

                frontier = new LinkedList<>();
                for(Integer dest: G.outEdges.get(r)){
                    if(checkSumCondition(T, r, dest)){
                        frontier.push(Edge.getEdge(r, dest));
                    }
                }

                grow(G, T, frontier, r);
            }
            this.doneSize--;

            System.out.println("Found " + results.size() + " trees");
            valid = new LinkedList<>();
            for (Tree T1 : results){
                if (T1.isValid()){
                    valid.add(T1);
                }
            }
            System.out.println("Found " + valid.size() + " valid trees");

            //System.out.println("Sum condition helped prune [" + Calder.helpedPrune + "] out of [" +
            // Calder.opp1 + "] opportunities");

        } while(valid.size() == 0 && doneSize > 0);

        return valid;
    }

    protected abstract void output(Tree T, int root);

    /**
     * Finds all spanning trees in the graph G containing the tree T rooted at root
     * @param G ancestry graph
     * @param T working tree
     * @param frontier edges to explore
     * @param root root of T
     */
    protected void grow(Graph G, Tree T, LinkedList<Edge> frontier, int root){
        // if T has V vertices then we've found a spanning tree so we're done
        if(T.vertices.size() == doneSize){
            output(T, root);
            return;
        }

        LinkedList<Edge> FF = new LinkedList<>();
        boolean b = false;
        Edge e;
        int u, v;
        LinkedList<Edge> newFrontier;
        LinkedList<Edge> toRemove;

        while(!b && !frontier.isEmpty()) {
            e = frontier.pop();
            u = e.source;
            v = e.dest;

            // Gabow-Myers invariants: u is in T, v is not in T, e = (u,v) is not in T
            // May have removed the vertex u if I removed the last edge in the graph
            if(T.vertices.size() > 0){
                assert T.vertices.contains(u);
                assert !T.outEdges.get(u).contains(v);
            }
            assert !T.vertices.contains(v);

            T.addEdge(u, v);

            newFrontier = new LinkedList<>(frontier);

            // Add all new frontier outEdges from my destination
            for (Integer dest : G.outEdges.get(v)) {
                if (!T.vertices.contains(dest) && checkSumCondition(T, v, dest)) {
                    newFrontier.push(Edge.getEdge(v, dest));
                }
            }

            // Remove all outEdges from frontier that could violate tree conditions
            toRemove = new LinkedList<>();
            for (Edge newE : newFrontier){
                if (newE.dest == v) {
                    toRemove.add(newE);
                } else if(newE.source == u && !checkSumCondition(T, u, newE.dest)){
                    toRemove.add(newE);
                }
            }
            newFrontier.removeAll(toRemove);
            grow(G, T, newFrontier, root);

            T.removeEdge(u, v);
            G.removeEdge(u, v);

            FF.push(e);

            if(L != null){
                b = L.bridgeTest(v);
            }
        }

        for (Edge myE: FF){
            frontier.push(myE);
            G.addEdge(myE.source, myE.dest);
        }
    }

    private boolean checkSumCondition(Tree T, int parent, int newChild){
        int t;
        VertexData parentData = T.vertexData.get(parent);
        VertexData newChildData = T.vertexData.get(newChild);
        VertexData childData;

        Calder.opp1++;

        // Before adding the new edge to the tree, check the Sum Condition
        double[] pSlack = parentData.maxf.clone();
        if(T.outEdges.containsKey(parent)){
            for(Integer child : T.outEdges.get(parent)){
                childData = T.vertexData.get(child);
                for(t = 0; t < nSamples; t++){
                    pSlack[t] -= childData.minf[t];
                }
            }
        }
        for(t = 0; t < nSamples; t++){
            pSlack[t] -= newChildData.minf[t];
            if(pSlack[t] < 0){
                Calder.helpedPrune++;
                return false;
            }
        }

        return true;
    }

}
