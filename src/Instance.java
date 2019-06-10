import org.apache.commons.math3.distribution.BetaDistribution;

import java.io.File;
import java.io.FileNotFoundException;
import java.text.DecimalFormat;
import java.util.*;
import java.util.stream.Collectors;

import static java.lang.System.exit;

public class Instance {
    private final int[][] variantReads;
    private final int[][] normalReads;
    final double[][][] intervals;
    final String[] colLabels;
    final String[] rowLabels;
    final int nSamples;
    final int nMuts;
    private boolean printedEffectiveConfidence = false;

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

    private Instance(int[][] normalReads, int[][] variantReads, double[][][] intervals,
                     String[] colLabels, String[] rowLabels, int nSamples, int nMuts){
        this.normalReads = normalReads;
        this.variantReads = variantReads;
        this.intervals = intervals;
        this.colLabels = colLabels;
        this.rowLabels = rowLabels;
        this.nSamples = nSamples;
        this.nMuts = nMuts;
    }

    private static Instance constructInstance(int[][] normalReads, int[][] variantReads, String[] colLabels, String[] rowLabels){
        assert normalReads != null;
        assert normalReads.length > 0;
        assert normalReads[0].length > 0;
        assert normalReads.length == variantReads.length;
        assert normalReads[0].length == variantReads[1].length;

        int nSamples = normalReads.length;
        int nMuts = normalReads[0].length;
        double[][] mins = new double[normalReads.length][normalReads[0].length];
        double[][] maxes = new double[normalReads.length][normalReads[0].length];

        int i, j;
        double freq;
        double[] interval;
        ArrayList<Integer> toRemove = new ArrayList<>();

        // Compute denominator for multiple hypothesis correction
        int denom = nSamples * nMuts;

        double resultingConfidence = 1 - ((1 - Main.CONFIDENCE) / (double) denom);

        // Print the effective confidence level
        if(Main.PRINT_EFFECTIVE_CONFIDENCE) {
            DecimalFormat df = new DecimalFormat("00.00");
            System.out.println("Effective confidence level (after multiple hypothesis correction): " + df.format(100 * resultingConfidence) + "%");
        }

        boolean gaveWarning = false;
        for(i = 0; i < normalReads.length; i++){
            for(j = 0; j < normalReads[0].length; j++){
                interval = computeBetaInterval(variantReads[i][j] + 1, normalReads[i][j] + 1, resultingConfidence);
                if(interval.length == 3 && interval[2] == -1){
                    if (Main.REMOVE_CNA){
                        if(!toRemove.contains(j)){
                            toRemove.add(j);
                        }
                        continue;
                    } else {
                        if(!gaveWarning){
                            System.out.println("WARNING: found mutation(s) with abnormally high frequency (confidence interval lower bound > 0.5).");
                            gaveWarning = true;
                        }
                        System.out.println("Mutation " + colLabels[j] + " had original interval [" + interval[0] + ", " +
                                interval[1] + "], truncating to [0.49999, 0.5].");
                        interval[0] = 0.99998;
                        interval[1] = 1;
                    }

                }
                mins[i][j] = Math.floor(interval[0] * Math.pow(10, Main.PRECISION_DIGITS)) /  Math.pow(10, Main.PRECISION_DIGITS);
                maxes[i][j] = Math.floor(interval[1] * Math.pow(10, Main.PRECISION_DIGITS)) /  Math.pow(10, Main.PRECISION_DIGITS);
            }
        }
        if(gaveWarning){
            System.out.println("(to remove these mutations from consideration, specify the -r option)");
        }
        double[][][] result = new double[][][]{mins, maxes};

        if (Main.REMOVE_CNA && toRemove.size() > 0){
            // create a new Instance object with abnormally high-frequency mutations/clusters removed

            // Determine which mutations do not have abnormally high frequencies
            ArrayList<Integer> remainingMuts = new ArrayList<Integer>();
            HashMap<Integer, Integer> newToOld = new HashMap<>();
            int k = 0;
            for(i = 0; i < nMuts; i++){
                if(k < toRemove.size() && toRemove.get(k) == i){
                    // mutation i has abnormally high frequencies
                    k++;
                } else {
                    // mutation i is not in toRemove
                    remainingMuts.add(i);
                    newToOld.put(i - k, i);
                }
            }

            int nMutsRemaining = remainingMuts.size();
            assert nMutsRemaining < nMuts;

            // Construct new read count matrices
            String[] newCols = new String[nMutsRemaining];
            int[][] newNormalReads = new int[nSamples][nMutsRemaining];
            int[][] newVariantReads = new int[nSamples][nMutsRemaining];
            for (int t = 0; t < nSamples; t++) {
                for (i = 0; i < nMutsRemaining; i++) {
                    newNormalReads[t][i] = normalReads[t][newToOld.get(i)];
                    newVariantReads[t][i] = variantReads[t][newToOld.get(i)];
                }
            }
            // Construct new column label array
            for (i = 0; i < nMutsRemaining; i++){
                newCols[i] = colLabels[newToOld.get(i)];
            }
            System.out.println("Removed mutations with abnormally high frequencies: " +
                    toRemove.stream().map(x -> colLabels[x]).collect(Collectors.toList()));

            return constructInstance(newNormalReads, newVariantReads, newCols, rowLabels);
        } else {
            double[][][] intervals = new double[2][][];
            intervals[0] = result[0];
            intervals[1] = result[1];

            return new Instance(normalReads, variantReads, intervals, colLabels, rowLabels, nSamples, nMuts);
        }

    }


    public static Instance fromFile(String fname){
        return readInstanceFromFile(fname);
    }


    private static double[] computeBetaInterval(int alpha, int beta, double resultingConfidence){
        double[] result = new double[2];
        BetaDistribution b = new BetaDistribution(alpha, beta);
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


        if (result[0] > .5){
            return new double[]{result[0], result[1], -1};
        }


        result[0] *= 2;
        result[1] *= 2;
        if(result[1] > 1){
            result[1] = 1;
        }
        // If the lower end of the interval is sufficiently low, set it to 0 to allow for absence
        if (result[0] < Main.MINIMUM_USAGE_H){
            result[0] = 0;
        }

        return result;
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
        int length, i;
        String header;
        File f = new File(filename);
        try (Scanner scanner = new Scanner(f)) {
            //Reads in header row


            header = scanner.nextLine();
            tokens = header.split("[\t ]");
            //System.out.println(Arrays.toString(tokens));
            String[] realTokens = new String[tokens.length];
            int j = 0;

            // Remove empty entries from header row
            for(i = 0; i < tokens.length; i++){
                if(tokens[i].equals("")){
                    System.out.println("Found empty entry in column [" + (i + 1) + "] of the header row, skipping it.");
                } else {
                    realTokens[j++] = tokens[i];
                }
            }
            realTokens = Arrays.copyOfRange(realTokens, 0, j);
            if(j % 2 == 1){
                System.err.println("Unable to parse file: found an odd number of non-empty entries in the header row.");
                System.exit(1);
            }

            colLabels = new String[(realTokens.length) / 2];
            for (i = 0; i < colLabels.length; i++){
                colLabels[i] = realTokens[2 * i];
            }
            firstLength = realTokens.length + 1;
            //System.out.println(Arrays.toString(colLabels));
            //System.out.println(firstLength);

            System.out.println("Found [" + colLabels.length + "] mutations/clusters in header " +
                    "(expecting 1 + (2 * num_clusters) = " +  (2 * colLabels.length + 1) + " entries in each row).");

            int rowNum = 0;
            while (scanner.hasNextLine()) {
                rowNum++;
                line = scanner.nextLine();

                normalRow = new ArrayList<>();
                variantRow = new ArrayList<>();

                tokens = line.split("\t");
                length = tokens.length;

                if (length != firstLength){
                    System.err.println("Unable to parse file: row [" + rowNum + "] has [" + length + "] entries, " +
                            "expected [" + (colLabels.length * 2 + 1) + "].");
                    exit(1);
                }

                rowLabels.add(tokens[0]);

                // Starts at 1 to skip row labels
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

        return constructInstance(nCounts, vCounts, colLabels, actualRowLabels);
    }

    public String toString(){
        StringBuilder sb = new StringBuilder();

        //sb.append(Calder.print2DArray(normalReads, colLabels, rowLabels));
        //sb.append(Calder.print2DArray(variantReads, colLabels, rowLabels));
        sb.append("Interval lower bounds:\n");
        sb.append(Util.print2DArray(intervals[0], colLabels, rowLabels));
        sb.append("Interval upper bounds:\n");
        sb.append(Util.print2DArray(intervals[1], colLabels, rowLabels));


        return sb.toString();
    }

}
