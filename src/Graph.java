import java.io.Serializable;
import java.util.*;

import static java.lang.System.exit;

public class Graph implements Serializable {
    final HashSet<Integer> vertices;
    static HashMap<Integer, VertexData> vertexData = null;
    final HashMap<Integer, HashSet<Integer>> outEdges;
    private final HashMap<Integer, HashSet<Integer>> inEdges;

    static HashMap<Edge, Integer> edgeEncoding = null;
    static HashMap<Integer, Edge> edgeDecoding = null;
    static int nEdges = 0;

    public Graph(){
        vertices = new HashSet<>();
        outEdges = new HashMap<>();
        inEdges = new HashMap<>();
        if(vertexData == null){
            vertexData = new HashMap<>();
        }
    }

    public Graph(Graph G){
        if(vertexData == null){
            vertexData = new HashMap<>();
        }

        vertices = new HashSet<>();
        outEdges = new HashMap<>();
        inEdges = new HashMap<>();

        HashSet<Integer> myEdges;
        for (Integer v: G.vertices){
            vertices.add(v);

            myEdges = new HashSet<>(G.outEdges.get(v));
            outEdges.put(v, myEdges);

            myEdges = new HashSet<>(G.inEdges.get(v));
            inEdges.put(v, myEdges);

        }
    }

    public Graph newEmptyGraph(){
        return new Graph();
    }

    // This method is used to originally construct the graph - stores the data associated with this vertex
    public void addVertex(VertexData v){
        if(vertices.contains(v.id)){
            return;
        }

        vertices.add(v.id);

        vertexData.put(v.id, v);

        outEdges.put(v.id, new HashSet<>());
        inEdges.put(v.id, new HashSet<>());

    }

    // This method is used to manipulate the graph during the algorithm
    public void addVertex(int id){
        if(vertices.contains(id)){
            return;
        }

        vertices.add(id);
        outEdges.put(id, new HashSet<>());
        inEdges.put(id, new HashSet<>());
    }

    // Removes a vertex from the graph (for the purposes of structure - does not remove its data)
    public void removeVertex(int id){
        vertices.remove(id);
        outEdges.remove(id);
        inEdges.remove(id);
    }

    // Adds an edge (tuple of int ids) to the graph, adding its endpoints to the graph if they are not present
    public void addEdge(int source, int dest){
        if (!vertices.contains(source)){
            addVertex(source);
        }
        if (!vertices.contains(dest)){
            addVertex(dest);
        }
        outEdges.get(source).add(dest);
        inEdges.get(dest).add(source);
    }

    // Removes an edge (tuple of int ids) from the graph, removing its endpoints from the vertices set if they become
    // disconnected without the edge.
    public void removeEdge(int source, int dest){
        if(!outEdges.containsKey(source)){
            Main.log("ERROR: tried to remove edge but source not found: " + source);
            throw new AssertionError();
        }
        HashSet<Integer> candidates = outEdges.get(source);
        if (!candidates.contains(dest)){
            Main.log("ERROR: tried to remove edge not in graph: (" + source + "," + dest + ")");
            throw new AssertionError();
        }
        assert inEdges != null;
        assert inEdges.containsKey(dest);
        assert inEdges.get(dest).contains(source);

        outEdges.get(source).remove(dest);
        inEdges.get(dest).remove(source);

        // If the removal of this edge disconnects the source vertex from the graph, remove it
        if(outEdges.get(source).isEmpty() && inEdges.get(source).isEmpty()){
            removeVertex(source);
        }

        // If the removal of this edge disconnects the destination vertex from the graph, remove it
        if(outEdges.get(dest).isEmpty() && inEdges.get(dest).isEmpty()){
            removeVertex(dest);
        }
    }

    // Checks for an edge (w,v) where w is not a descendant of v
    public boolean bridgeTest(int v){
        HashSet<Integer> descendants = getDescendants(v);

        if(!inEdges.containsKey(v)){
            // this can happen because the last result tree may not have contained this vertex
            return false;
        }
        for(Integer w: inEdges.get(v)){
            if(!descendants.contains(w)){
                return false;
            }
        }

        return true;
    }

    // Use DFS to find all proper descendants of a vertex in this graph (assumed to be a tree)
    public HashSet<Integer> getDescendants(int v){
        LinkedList<Integer> next = new LinkedList<>();
        HashSet<Integer> result = new HashSet<>();
        HashSet<Integer> visited = new HashSet<>();

        next.push(v);
        visited.add(v);

        Integer curr;
        while(!next.isEmpty()){
            curr = next.pop();

            if(outEdges.get(curr) == null){
                //System.out.println(outEdges);
                //System.out.println(curr);

                // TODO: figure out where this actually comes from and fix it
                continue;
            }
            for(Integer dest: outEdges.get(curr)){
                if (visited.contains(dest)){
                    Main.log("ERROR: tree condition violated (DFS visited vertex " + dest + " twice");
                    break;
                }
                visited.add(dest);

                next.push(dest);
                result.add(dest);
            }
        }

        return result;
    }

    // output the graph in DOT format
    public String toDot(){
        StringBuilder b = new StringBuilder();
        b.append("digraph Graph {\n");

        for(Integer v: vertices){
            for(Integer dest: outEdges.get(v)){
                b.append("    ");
                b.append(v);
                b.append(" -> ");
                b.append(dest);
                b.append( ";\n");
            }
        }
        b.append("}");

        return b.toString();
    }

    public String toString(){
        StringBuilder b = new StringBuilder();
        for(Integer v: vertices){
            b.append(v);
            b.append(": ");
            for(Integer dest: outEdges.get(v)){
                b.append('(');
                b.append(v);
                b.append(',');
                b.append(dest);
                b.append(") ");
            }
            b.append('\n');
        }
        return b.toString();
    }

    /**
     * Uses Tarjan's strongly connected components algorithm to find all possible spanning tree roots
     * @return list of root candidates
     */
    public ArrayList<Integer> findRootCandidates(){
        ArrayList<Integer> result = new ArrayList<>();

        // run Tarjan's SCC algorithm
        StrongConnect s = new StrongConnect();
        s.run();

        // find the strongly connected component(s) with no incoming edges
        boolean foundOne = false;
        boolean valid = true;

        for(ArrayList<Integer> component : s.components){
            for(Integer v : component){
                if(!valid){
                    break;
                }
                for(Integer source : inEdges.get(v)){
                    if (s.componentOf[source] != s.componentOf[v]){
                        // Can't be the root component if it has an incoming edge from another component
                        valid = false;
                    }
                }

            }
            if(valid && !foundOne){
                // This SCC has no incoming edges, therefore has our root candidates
                result = component;
                foundOne = true;
            } else if (valid){
                // Already found an SCC with no incoming edges so there are no possible roots
                return new ArrayList<>();
            }
            valid = true;

        }
        return result;
    }

    // Run Tarjan's strongly connected components algorithm
    private class StrongConnect {
        //TODO: internally, assumes that vertices are indexed from 0 to V-1
        private int indexIter;
        private int componentIter;
        private int[] index;
        private int[] lowlink;
        private boolean[] onStack;
        private LinkedList<Integer> S;
        ArrayList<ArrayList<Integer>> components;
        int[] componentOf;

        private StrongConnect(){
            indexIter = 1;
            componentIter = 0;
            components = new ArrayList<>();
            S = new LinkedList<>();
            index = new int[vertices.size()];
            lowlink = new int[vertices.size()];
            onStack = new boolean[vertices.size()];
            componentOf = new int[vertices.size()];
        }

        private void run(){
            for(Integer v: vertices){
                if(index[v] == 0){
                    strongconnect(v);
                }
            }
        }
        private void strongconnect(int v){
            index[v] = indexIter;
            lowlink[v] = indexIter;
            indexIter++;
            S.push(v);
            onStack[v] = true;

            for(Integer w : outEdges.get(v)){
                if(index[w] == 0){
                    strongconnect(w);
                    lowlink[v] = Math.min(lowlink[v], lowlink[w]);
                } else if (onStack[w]){
                    lowlink[v] = Math.min(lowlink[v], index[w]);
                }
            }

            if(lowlink[v] == index[v]){
                ArrayList<Integer> myComponent = new ArrayList<>();
                int w;
                do{
                    w = S.pop();
                    onStack[w] = false;
                    componentOf[w] = componentIter;
                    myComponent.add(w);
                } while(w != v);
                components.add(myComponent);
                componentIter++;
            }
        }



    }

    // Builds an ancestry graph by considering all pairwise evolutionary relationships meeting the ancestry conditions
    public static Graph buildAncestryGraph(Collection<VertexData> mutations){
        edgeEncoding = new HashMap<>();
        edgeDecoding = new HashMap<>();
        int i = 0;

        Graph G = new Graph();

        for(VertexData v : mutations){
            G.addVertex(v);
        }
        Edge e;
        boolean validEdge;
        for(VertexData u : mutations){
            for(VertexData v : mutations){
                if (!u.equals(v)){
                    validEdge = true;
                    for(int t = 0; t < u.maxf.length; t++){
                        if(u.maxf[t] < v.minf[t]){
                            validEdge = false;
                            break;
                        }
                    }
                    if(validEdge){
                        G.addEdge(u.id, v.id);
                        e = Edge.getEdge(u.id, v.id);
                        edgeEncoding.put(e, i);
                        edgeDecoding.put(i++, e);
                    }

                }
            }
        }

        /*
        System.out.println(edgeEncoding);
        System.out.println(i);
        */
        nEdges = i;

        return G;
    }

    public static EncodedTree encodeTree(Tree T){
        BitSet edges = new BitSet(nEdges);
        Edge e;
        for(Integer source : T.vertices){
            for(Integer dest : T.outEdges.get(source)){
                e = Edge.getEdge(source, dest);
                edges.set(edgeEncoding.get(e));
            }
        }
        //TODO: use tree encoding to save space during enumeration
        return new EncodedTree(T.root, edges);
    }
    public static Tree decodeTree(EncodedTree enc){
        Tree T = new Tree(new Graph(), enc.root, vertexData.values().iterator().next().minf.length);
        int idx = enc.edges.nextSetBit(0);
        Edge e;
        while(idx != -1){
            e = edgeDecoding.get(idx);
            T.addEdge(e.source, e.dest);
            idx = enc.edges.nextSetBit(idx + 1);
        }
        return T;
    }

}
