import java.util.HashMap;

public class Edge {
    private static final HashMap<Integer, HashMap<Integer, Edge>> allEdges = new HashMap<>();

    public final int source;
    public final int dest;

    private Edge(int source, int dest){
        this.source = source;
        this.dest = dest;
    }

    public static Edge getEdge(int source, int dest){
        if (allEdges.containsKey(source) && allEdges.get(source).containsKey(dest)){
            return allEdges.get(source).get(dest);
        } else {
            Edge e = new Edge(source, dest);
            if(!allEdges.containsKey(source)){
                allEdges.put(source, new HashMap<>());
            }
            allEdges.get(source).put(dest, e);

            return e;
        }
    }

    public String toString(){
        return "(" + source + "," + dest + ")";
    }

    public boolean equals(Object o){
        if (o instanceof Edge){
            Edge e = (Edge) o;
            return source == e.source && dest == e.dest;
        } else {
            return false;
        }
    }
}
