import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class GabowMyersMT {
    private Graph myG;
    private CountDownLatch latch;
    private ExecutorService pool = Executors.newWorkStealingPool(Main.THREADS);

    private int nSamples;
    private HashMap<Integer, Integer> lastEdges = new HashMap<>();

    public GabowMyersMT(Graph G, int nSamples){
        this.myG = G;
    }

    public LinkedList<Tree> run() throws InterruptedException {
        // Make a copy of G
        Graph G;
        Tree T;
        LinkedList<Edge> frontier;

        int doneSize = myG.vertices.size();
        LinkedList<Tree> valid;
        Collection<Integer> rootCandidates;
        LinkedList<RunGMTask> tasks = new LinkedList<>();

        do {

            if(doneSize == Graph.vertexData.size()){
                rootCandidates = myG.findRootCandidates();
            } else {
                rootCandidates = myG.vertices;
            }
            latch = new CountDownLatch(rootCandidates.size());

            System.out.println("Looking for trees with " + doneSize + " vertices");
            System.out.println("Found " + rootCandidates.size() + " candidate roots");

            //TODO: print information on number of root candidates

            for(Integer r : rootCandidates){
                Calder.resetCounts();

                // Make a copy of G for each iteration to avoid unexpected behavior
                G = new Graph(myG);

                // Initialize initialT with just the root
                T = new Tree(myG.newEmptyGraph(), r, nSamples);

                frontier = new LinkedList<>();
                for(Integer dest: G.outEdges.get(r)){
                    if(checkSumCondition(T, r, dest)){
                        frontier.push(Edge.getEdge(r, dest));
                    }
                }
                RunGMTask task = new RunGMTask(G, T, frontier, r, doneSize, latch);
                tasks.add(task);
                pool.submit(task);
            }
            latch.await();
            LinkedList<Tree> results = new LinkedList<>();
            for(RunGMTask task: tasks){
                results.addAll(task.myResults);
            }
            doneSize--;

            System.out.println("Found " + results.size() + " trees");
            valid = new LinkedList<>();
            for (Tree T1 : results){
                if (T1.isValid()){
                    valid.add(T1);
                }
            }
            System.out.println("Found " + valid.size() + " valid trees");
            System.out.println("Sum condition helped prune [" + Calder.helpedPrune + "] out of [" +
                    Calder.opp1 + "] opportunities");
        } while(valid.size() == 0 && doneSize > 0);

        return valid;
    }

    private class RunGMTask implements Runnable{
        private final Graph initialG;
        private final Tree initialT;
        private Tree L;
        private LinkedList<Edge> initialFrontier;
        private final int root;
        private final CountDownLatch latch;

        private final int doneSize;
        final LinkedList<Tree> myResults;

        public RunGMTask(Graph initialG, Tree initialT, LinkedList<Edge> initialFrontier, int root, int doneSize, CountDownLatch latch){
            this.initialG = initialG;
            this.root = root;
            this.doneSize = doneSize;
            this.latch = latch;

            this.initialT = initialT;
            this.initialFrontier = initialFrontier;
            this.myResults = new LinkedList<>();
        }

        @Override
        public void run() {
            // enumerate trees
            grow(initialG, initialT, initialFrontier);

            System.out.println("Thread: found " + myResults.size() + " trees of size " + doneSize);
            latch.countDown();
        }

        /**
         * Finds all spanning trees in the graph G containing the tree initialT rooted at root
         * @param G working ancestry graph
         * @param T working tree
         * @param frontier edges to explore
         */
        private void grow(Graph G, Tree T, LinkedList<Edge> frontier){
            // if initialT has V vertices then we've found a spanning tree so we're done
            if(T.vertices.size() == doneSize){
                L = new Tree(T, root, nSamples);
                if(L.isValid()){
                    myResults.add(L);
                }
                return;
            }

            LinkedList<Edge> FF = new LinkedList<>();
            boolean b = false;
            Edge e;
            int u, v;
            LinkedList<Edge> newFrontier;
            LinkedList<Edge> toRemove;

            while(!b && !frontier.isEmpty()) {
                //System.out.println(initialFrontier);
                //System.out.println(lastEdges);

                //System.out.println("About to visit " + initialFrontier.peek() + ". My tree has " + initialT.vertices.size() + " nodes.");

                e = frontier.pop();
                u = e.source;
                v = e.dest;



                // Gabow-Myers invariants: u is in initialT, v is not in initialT, e = (u,v) is not in initialT
                // May have removed the vertex u if I removed the last edge in the graph
                if(T.vertices.size() > 0){
                    assert T.vertices.contains(u);
                    assert !T.outEdges.get(u).contains(v);
                }
                assert !T.vertices.contains(v);

                T.addEdge(u, v);

                newFrontier = new LinkedList<>(frontier);

                // Add all new initialFrontier outEdges from my destination
                for (Integer dest : G.outEdges.get(v)) {
                    if (!T.vertices.contains(dest) && checkSumCondition(T, v, dest)) {
                        newFrontier.push(Edge.getEdge(v, dest));
                    }
                }

                // Remove all outEdges from initialFrontier that could violate tree conditions
                toRemove = new LinkedList<>();
                for (Edge newE : newFrontier){
                    if (newE.dest == v) {
                        toRemove.add(newE);
                    } else if(newE.source == u && !checkSumCondition(T, u, newE.dest)){
                        toRemove.add(newE);
                    }
                }
                newFrontier.removeAll(toRemove);

                //System.out.println(newFrontier);
                grow(G, T, newFrontier);

                T.removeEdge(u, v);
                G.removeEdge(u, v);

                FF.push(e);

                if(L != null){
                    b = L.bridgeTest(v);
                }

                //System.out.println("Finished visiting " + e);
            }

            for (Edge myE: FF){
                frontier.push(myE);
                G.addEdge(myE.source, myE.dest);
            }
        }


    }

    private static boolean checkSumCondition(Tree T, int parent, int newChild){
        int t;
        VertexData parentData = Graph.vertexData.get(parent);
        VertexData newChildData = Graph.vertexData.get(newChild);
        VertexData childData;
        int nSamples = parentData.minf.length;

        Calder.opp1++;

        // Before adding the new edge to the tree, check the Sum Condition
        double[] pSlack = parentData.maxf.clone();
        if(T.outEdges.containsKey(parent)){
            for(Integer child : T.outEdges.get(parent)){
                childData = Graph.vertexData.get(child);
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
