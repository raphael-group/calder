import org.apache.commons.cli.*;

public class JarMain {
    public static void main(String[] args){
        Options options = new Options();
        Option input = new Option("i", "input", true, "input file path");
        input.setRequired(true);
        options.addOption(input);

        Option output = new Option("o", "output", true, "output file path");
        output.setRequired(false);
        options.addOption(output);

        Option alpha = new Option("a", "alpha", true, "confidence level alpha (default 0.9)");
        alpha.setRequired(false);
        options.addOption(alpha);

        // Gurobi automatically uses all available processors
        /*
        Option threads = new Option("t", "threads", true, "number of threads (default 1)");
        threads.setRequired(false);
        options.addOption(threads);
        */

        Option detectionThreshold = new Option("h", "threshold", true, "detection threshold h (default 0.01)");
        detectionThreshold.setRequired(false);
        options.addOption(detectionThreshold);

        Option printCI = new Option("n", "intervals", false, "print confidence intervals");
        printCI.setRequired(false);
        options.addOption(printCI);

        Option nSolns = new Option("s", "solutions", false, "maximum number of optimal solutions to return (default 1)");
        nSolns.setRequired(false);
        options.addOption(nSolns);

        Option printConfidence = new Option("c", "printconf", true, "print effective confidence level");
        printConfidence.setRequired(false);
        options.addOption(printConfidence);

        Option printGraph = new Option("g", "print-graph", false, "print ancestry graph");
        printGraph.setRequired(false);
        options.addOption(printGraph);


        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd = null;

        try{
            cmd = parser.parse(options, args);
        } catch (ParseException e){
            System.out.println(e.getMessage());
            formatter.printHelp("calder", options);

            System.exit(1);
        }

        try {
            Main.INFILE = cmd.getOptionValue("input");
            if (cmd.hasOption("output")){
                Main.OUTFILE = cmd.getOptionValue("output");
            }
            // Gurobi automatically uses all available processors
            /*
            if (cmd.hasOption("threads")) {
                Main.THREADS = Integer.parseInt(cmd.getOptionValue("threads"));
            }
            */
            if (cmd.hasOption("threshold")) {
                Main.MINIMUM_USAGE_H = Double.parseDouble(cmd.getOptionValue("threshold"));
                if(Main.MINIMUM_USAGE_H <= 0 || Main.MINIMUM_USAGE_H >= 1){
                    System.out.println("Minimum usage h must be in (0, 1)");
                    System.exit(1);
                }
            }
            if (cmd.hasOption("alpha")) {
                Main.CONFIDENCE = Double.parseDouble(cmd.getOptionValue("alpha"));
                if(Main.CONFIDENCE <= 0 || Main.CONFIDENCE >= 1){
                    System.out.println("Confidence level alpha must be in (0, 1)");
                    System.exit(1);
                }
            }
            if (cmd.hasOption("solutions")) {
                Main.MAX_NUM_OPTIMA_OUTPUT = Integer.parseInt(cmd.getOptionValue("solutions"));
                if(Main.MAX_NUM_OPTIMA_OUTPUT <= 0){
                    System.out.println("Maximum number of solutions must be a positive integer");
                    System.exit(1);
                }
            }

            if(cmd.hasOption("print-graph")){
                Main.PRINT_ANCESTRY_GRAPH = true;
            }
            if(cmd.hasOption("intervals")){
                Main.PRINT_CONFIDENCE_INTERVALS = true;
            }
            if(cmd.hasOption("printconf")){
                Main.PRINT_EFFECTIVE_CONFIDENCE = true;
            }
        } catch(NumberFormatException e){
            System.out.println(e.getMessage());
            formatter.printHelp("calder", options);
            System.exit(1);
        }

        Main.main(args);
    }

}
