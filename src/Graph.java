import java.io.Serializable;
import java.util.*;

public class Graph implements Serializable {
    final HashSet<Integer> vertices;
    static HashMap<Integer, VertexData> vertexData = null;
    final HashMap<Integer, HashSet<Integer>> outEdges;
    final HashMap<Integer, HashSet<Integer>> inEdges;

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
            Util.log("ERROR: tried to remove edge but source not found: " + source);
            throw new AssertionError();
        }
        HashSet<Integer> candidates = outEdges.get(source);
        if (!candidates.contains(dest)){
            Util.log("ERROR: tried to remove edge not in graph: (" + source + "," + dest + ")");
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

        nEdges = i;

        return G;
    }

}
