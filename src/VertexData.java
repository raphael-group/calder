import java.lang.reflect.Array;
import java.util.Arrays;

public class VertexData {
    public final int id;
    public final double[] minf;
    public final double[] maxf;

    public VertexData(int id, double[] minf, double[] maxf){
        this.id = id;

        // Confidence intervals around the column of the frequency matrix corresponding to this mutation
        this.minf = minf;
        this.maxf = maxf;
    }

    public boolean equals(Object o){
        if (o instanceof VertexData){
            VertexData other = (VertexData) o;
            return other.id == id && Arrays.equals(minf, other.minf) && Arrays.equals(maxf, other.maxf);
        }
        return false;
    }

    public String toString(){
        StringBuilder b = new StringBuilder();
        b.append(id);
        b.append(": ");
        b.append(Arrays.toString(minf));
        b.append(", ");
        b.append(Arrays.toString(maxf));
        return b.toString();
    }
}
