import java.util.BitSet;

public class EncodedTree {
    final int root;
    final BitSet edges;

    public EncodedTree(int root, BitSet edges){
        this.root = root;
        this.edges = edges;
    }
}
