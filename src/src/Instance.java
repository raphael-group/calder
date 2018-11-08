import org.apache.commons.math3.distribution.BetaDistribution;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

import static java.lang.System.exit;

public class Instance {
    final int[][] variantReads;
    final int[][] normalReads;
    final double[][][] intervals;
    final String[] colLabels;
    final String[] rowLabels;
    final int nSamples;
    final int nMuts;

    /*
    public Instance(int[][] normalReads, int[][] variantReads){
        assert normalReads != null;
        assert normalReads.length > 0;
        assert normalReads[0].length > 0;
        assert normalReads.length == variantReads.length;
        assert normalReads[0].length == variantReads[1].length;

        this.rowLabels = null;
        this.colLabels = null;

        this.normalReads = normalReads;
        this.variantReads = variantReads;

        this.nSamples = normalReads.length;
        this.nMuts = normalReads[0].length;

        double[][][] result = computeIntervals(normalReads, variantReads);
        intervals = new double[2][][];
        intervals[0] = result[0];
        intervals[1] = result[1];
    }
    */

    public Instance(int[][] normalReads, int[][] variantReads, String[] colLabels, String[] rowLabels){
        assert normalReads != null;
        assert normalReads.length > 0;
        assert normalReads[0].length > 0;
        assert normalReads.length == variantReads.length;
        assert normalReads[0].length == variantReads[1].length;

        this.normalReads = normalReads;
        this.variantReads = variantReads;

        this.nSamples = normalReads.length;
        this.nMuts = normalReads[0].length;

        double[][][] result = computeIntervals(normalReads, variantReads);
        intervals = new double[2][][];
        intervals[0] = result[0];
        intervals[1] = result[1];

        this.colLabels = colLabels;
        this.rowLabels = rowLabels;
    }

    public static Instance fromFile(String fname){
        return readInstanceFromFile(fname);
    }


    private double[][][] computeIntervals(int[][] normalReads, int[][] variantReads){
        double[][] mins = new double[normalReads.length][normalReads[0].length];
        double[][] maxes = new double[normalReads.length][normalReads[0].length];

        int i, j;
        double freq;
        double[] interval;
        for(i = 0; i < normalReads.length; i++){
            for(j = 0; j < normalReads[0].length; j++){
                interval = computeBetaInterval(variantReads[i][j] + 1, normalReads[i][j] + 1);
                mins[i][j] = interval[0];
                maxes[i][j] = interval[1];
            }
        }

        return new double[][][]{mins, maxes};

    }


    private double[] computeBetaInterval(int alpha, int beta){
        double[] result = new double[2];

        // Compute denominator for multiple hypothesis correction
        int denom = nSamples * nMuts;

        BetaDistribution b = new BetaDistribution(alpha, beta);
        double resultingConfidence = 1 - ((1 - Main.CONFIDENCE) / (double) denom);

        //TODO: add flag for printing resulting confidence
        //System.out.println(resultingConfidence);

        double side = resultingConfidence / 2;

        // Compute the 1-INTERVAL_WIDTH equal-tailed posterior probability interval

        //if the mean of the distribution is 0, only use the right tail
        if(alpha == 1){
            result[0] = 0;
            result[1] = b.inverseCumulativeProbability(side * 2);
        } else if (alpha == beta) {
            result[0] = b.inverseCumulativeProbability(1 - side * 2);
            result[1] = 1;
        } else {
            result[0] = b.inverseCumulativeProbability(0.5 - side);
            result[1] = b.inverseCumulativeProbability(0.5 + side);
        }
        assert result[0] < .5;
        result[0] *= 2;
        result[1] *= 2;
        if(result[1] > 1){
            result[1] = 1;
        }
        // If the lower end of the interval is sufficiently low, set it to 0 to allow for absence
        if (result[0] < Main.DETECTION_THRESHOLD_H){
            result[0] = 0;
        }

        return result;
    }

    /**
     * Converts a double matrix to an integer one for internal representation, returning the result and the scaling factor
     * @param mins matrix of lower CI bounds to convert to ints
     * @param maxes matrix of upper CI bounds to convert to ints
     * @return ToIntResult encapsulating the matrix and scaling factor
     */
    private static ToIntResult toInt(double[][] mins, double[][] maxes){
        assert mins != null;
        assert mins.length > 0;
        assert mins[0].length > 0;
        assert maxes != null;
        assert maxes.length > 0;
        assert maxes[0].length > 0;
        assert mins.length == maxes.length;
        assert mins[0].length == maxes[0].length;

        // Find the longest value, to use as a scaling factor for all values
        int maxLen = 0;
        int myLength;
        String[] tokens;
        double myVal;
        for (int i = 0; i < mins.length; i++){
            for(int j = 0; j < mins[0].length; j++){
                myVal = mins[i][j];

                try{
                    assert myVal >= 0;
                    assert myVal <= 1;
                } catch(AssertionError e){
                    System.err.println("Encountered value outside of [0,1] when frequency was expected");
                    exit(1);
                }

                tokens = (myVal + "").split("\\.");
                if(tokens.length < 2){
                    myLength = 0;
                } else {
                    myLength = tokens[1].length();
                }

                if (myLength > maxLen){
                    maxLen = myLength;
                }
            }
        }

        for (int i = 0; i < maxes.length; i++){
            for(int j = 0; j < maxes[0].length; j++){
                myVal = maxes[i][j];

                try{
                    assert myVal >= 0;
                    assert myVal <= 1;
                } catch(AssertionError e){
                    System.err.println("Encountered value outside of [0,1] when frequency was expected");
                    exit(1);
                }

                tokens = (myVal + "").split("\\.");
                if(tokens.length < 2){
                    myLength = 0;
                } else {
                    myLength = tokens[1].length();
                }

                if (myLength > maxLen){
                    maxLen = myLength;
                }
            }
        }

        // Truncate at 8 decimal places to avoid int overflow
        if(maxLen > 8){
            maxLen = 8;
        }

        int factor = (int) Math.pow(10, maxLen);
        int[][] newMins = new int[mins.length][mins[0].length];
        int[][] newMaxes = new int[mins.length][mins[0].length];
        for(int i = 0; i < mins.length; i++){
            for(int j = 0; j < mins[0].length; j++){
                newMins[i][j] = (int) Math.round(factor * mins[i][j]);
                newMaxes[i][j] = (int) Math.round(factor * maxes[i][j]);
            }
        }

        return new ToIntResult(newMins, newMaxes, factor);
    }

    /**
     * Reads an AncesTree-formatted TSV file with read counts and returns a corresponding Instance
     * @param filename TSV file containing integer read count values
     *                 (column and row labels, alternating normal and variant read count columns)
     * @return Instance encapsulating the input data
     */
    private static Instance readInstanceFromFile(String filename){
        ArrayList<ArrayList<Integer>> normalCounts = new ArrayList<>();
        ArrayList<ArrayList<Integer>> variantCounts = new ArrayList<>();
        String[] colLabels = new String[1];
        LinkedList<String> rowLabels = new LinkedList<>();

        // read TSV file row by row
        String line = "";
        String[] tokens = new String[1];
        ArrayList<Integer> normalRow, variantRow;
        String tkn1 = "", tkn2 = "";
        int firstLength = -1;
        int length, i, normalCount, variantCount;
        String header;
        File f = new File(filename);
        try (Scanner scanner = new Scanner(f)) {
            //Reads in header row


            header = scanner.nextLine();
            tokens = header.split("\t");

            assert tokens.length % 2 == 1;

            colLabels = new String[(tokens.length - 1) / 2];
            for (i = 0; i < colLabels.length; i++){
                colLabels[i] = tokens[2 * i + 1];
            }

            while (scanner.hasNextLine()) {
                line = scanner.nextLine();


                normalRow = new ArrayList<>();
                variantRow = new ArrayList<>();

                tokens = line.split("\t");

                length = tokens.length;
                if(firstLength < 0){
                    firstLength = length;
                }
                assert length == firstLength;
                rowLabels.add(tokens[0]);

                // Starts at 1 to skip row labels
                //TODO: actually check to see if there are row labels, or just change this back
                for(i = 1; i < tokens.length; i+= 2){
                    tkn1 = tokens[i].replaceAll("[ \n]", "");
                    tkn2 = tokens[i + 1].replaceAll("[ \n]", "");
                    normalRow.add(Integer.parseInt(tkn1));
                    variantRow.add(Integer.parseInt(tkn2));
                }
                normalCounts.add(normalRow);
                variantCounts.add(variantRow);
            }
        } catch(NumberFormatException e){
            System.err.println("Unable to parse file: found a non-integer value");
            System.err.println(line);
            System.err.println(Arrays.toString(tokens));
            System.err.println(tkn1);
            System.err.println(tkn2);
            e.printStackTrace();
            exit(1);
        } catch(AssertionError e){
            System.err.println("Unable to parse file: lines have different lengths");
            exit(1);
        } catch (FileNotFoundException e) {
            System.err.println("Unable to parse file: file not found");
            e.printStackTrace();
            exit(1);
        }

        int[][] nCounts = new int[normalCounts.size()][normalCounts.get(0).size()];
        int[][] vCounts = new int[normalCounts.size()][normalCounts.get(0).size()];
        int j;
        for(i = 0; i < nCounts.length; i++){
            for(j = 0; j < nCounts[0].length; j++){
                nCounts[i][j] = normalCounts.get(i).get(j);
                vCounts[i][j] = variantCounts.get(i).get(j);
            }
        }
        String[] actualRowLabels = new String[rowLabels.size()];
        for(i = 0; i < rowLabels.size(); i++){
            actualRowLabels[i] = rowLabels.get(i);
        }
        assert colLabels.length == nCounts[0].length;
        assert rowLabels.size() == nCounts.length;

        return new Instance(nCounts, vCounts, colLabels, actualRowLabels);
    }

    public String toString(){
        StringBuilder sb = new StringBuilder();

        sb.append(Calder.print2DArray(normalReads, colLabels, rowLabels));
        sb.append(Calder.print2DArray(variantReads, colLabels, rowLabels));
        sb.append('\n');
        sb.append(Calder.print2DArray(intervals[0], colLabels, rowLabels));
        sb.append(Calder.print2DArray(intervals[1], colLabels, rowLabels));


        return sb.toString();
    }

}
