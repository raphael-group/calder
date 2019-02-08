import java.math.BigDecimal;
import java.text.DecimalFormat;
import java.time.LocalTime;

public class Util {
    // Pretty-prints a 2D array with column and row labels
    static String print2DArray(int[][] arr, String[] colLabels, String[] rowLabels){
        assert arr != null;
        assert colLabels != null;
        assert colLabels.length == arr[0].length;
        assert rowLabels != null;
        assert rowLabels.length == arr.length;

        StringBuilder b = new StringBuilder();
        int i, j;
        b.append("\t ");
        for(i = 0; i < arr[0].length; i++){
            b.append(colLabels[i]);
            b.append(spaces(6 - colLabels[i].length()));
        }
        b.append('\n');
        for(i = 0; i < arr.length; i++){
            b.append(rowLabels[i]);
            b.append(": ");
            for(j = 0; j < arr[i].length; j++){
                b.append(arr[i][j]);
                b.append(spaces(6 - (arr[i][j] + "").length()));
            }
            b.append('\n');
        }

        return b.toString();
    }

    static String print2DArray(int[][] arr, int[] colLabels, int[] rowLabels){
        int i;
        String[] newCols = new String[colLabels.length];
        for(i = 0; i < colLabels.length; i++){
            newCols[i] = colLabels[i] + "";
        }
        String[] newRows = new String[rowLabels.length];
        for(i = 0; i < rowLabels.length; i++){
            newRows[i] = rowLabels[i] + "";
        }
        return print2DArray(arr, newCols, newRows);
    }


    static String print2DArray(int[][] arr, int[] colLabels){
        int i;
        String[] newCols = new String[colLabels.length];
        for(i = 0; i < colLabels.length; i++){
            newCols[i] = colLabels[i] + "";
        }
        String[] rowLabels = new String[arr.length];
        for(i = 0; i < arr.length; i++){
            rowLabels[i] = "(" + i + ")";
        }
        return Util.print2DArray(arr, newCols, rowLabels);
    }

    static String print2DArray(BigDecimal[][] arr, String[] colLabels, String[] rowLabels){
        assert arr != null;
        assert colLabels != null;
        assert colLabels.length == arr[0].length;
        assert rowLabels != null;
        assert rowLabels.length == arr.length;

        DecimalFormat df = new DecimalFormat();
        df.setMaximumFractionDigits(3);

        StringBuilder b = new StringBuilder();
        int i, j;
        b.append("\t ");
        for(i = 0; i < arr[0].length; i++){
            b.append(colLabels[i]);
            b.append("\t");
        }
        b.append('\n');
        for(i = 0; i < arr.length; i++){
            b.append(rowLabels[i]);
            b.append(": ");
            b.append("\t");
            for(j = 0; j < arr[i].length; j++){
                b.append(df.format(arr[i][j].doubleValue()));
                b.append("\t");
            }
            b.append('\n');
        }

        return b.toString();
    }

    // Returns a string with the designated number of spaces
    private static String spaces(int howMany){
        StringBuilder b = new StringBuilder();
        for(int i = 0; i < howMany; i++){
            b.append(' ');
        }
        return b.toString();
    }


    static String print2DArray(double[][] arr, String[] colLabels, String[] rowLabels){
        assert arr != null;
        assert colLabels != null;
        assert colLabels.length == arr[0].length;
        assert rowLabels != null;
        assert rowLabels.length == arr.length;

        DecimalFormat df = new DecimalFormat();
        df.setMaximumFractionDigits(3);

        StringBuilder b = new StringBuilder();
        int i, j;
        b.append("\t ");
        for(i = 0; i < arr[0].length; i++){
            b.append(colLabels[i]);
            b.append("\t");
        }
        b.append('\n');
        for(i = 0; i < arr.length; i++){
            b.append(rowLabels[i]);
            b.append(": ");
            b.append("\t");
            for(j = 0; j < arr[i].length; j++){
                b.append(df.format(arr[i][j]));
                b.append("\t");
            }
            b.append('\n');
        }

        return b.toString();
    }

    public static void log(String message){
        LocalTime time = LocalTime.now();
        System.out.println(time + " " + message);
    }
}
