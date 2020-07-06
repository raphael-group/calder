import org.apache.commons.cli.*;

public class JarMain {
    public static void main(String[] args){
        Options options = new Options();
        Option input = new Option("i", "input", true, "input file path");
        input.setRequired(true);
        options.addOption(input);

        Option output = new Option("o", "output", true, "output directory");
        output.setRequired(true);
        options.addOption(output);

        Option alpha = new Option("a", "alpha", true, "confidence level alpha (default 0.9)");
        alpha.setRequired(false);
        options.addOption(alpha);

        Option detectionThreshold = new Option("h", "threshold", true, "detection threshold h (default 0.01)");
        detectionThreshold.setRequired(false);
        options.addOption(detectionThreshold);

        Option printCI = new Option("n", "intervals", false, "print confidence intervals");
        printCI.setRequired(false);
        options.addOption(printCI);

        Option nSolns = new Option("s", "solutions", true, "maximum number of optimal solutions to return (default 1)");
        nSolns.setRequired(false);
        options.addOption(nSolns);

        Option printConfidence = new Option("c", "printconf", false, "print effective confidence level");
        printConfidence.setRequired(false);
        options.addOption(printConfidence);

        Option printGraph = new Option("g", "print-graph", false, "print ancestry graph");
        printGraph.setRequired(false);
        options.addOption(printGraph);

        Option time = new Option("t", "time", false, "track and output timing information");
        time.setRequired(false);
        options.addOption(time);

        Option solver = new Option("v", "solver", true, "MILP solver back-end (default gurobi)");
        solver.setRequired(false);
        options.addOption(solver);

        Option nonlong = new Option("N", "nonlongitudinal", false, "remove longitudinal constraints");
        nonlong.setRequired(false);
        options.addOption(nonlong);

        Option objective = new Option("O", "objective", true, "objective function (l0, l1, or l0center)");
        objective.setRequired(false);
        options.addOption(objective);

        Option solver_verbose = new Option("sv", "verbose", true, "verbosity setting for JavaILP solver (effect depends on solver)");
        solver_verbose.setRequired(false);
        options.addOption(solver_verbose);

        Option solver_timeout = new Option("st", "timeout", true, "timeout setting for JavaILP solver");
        solver_timeout.setRequired(false);
        options.addOption(solver_timeout);

        Option tossCNA = new Option("r", "remove-columns", false, "discard mutations/clusters with abnormally high frequencies");
        tossCNA.setRequired(false);
        options.addOption(tossCNA);

        Option enumerateTrees = new Option("e", "enumerate", false, "enumerate all maximal trees instead of just optimal solutions");
        enumerateTrees.setRequired(false);
        options.addOption(enumerateTrees);

        Option details = new Option("d", "details", false, "print detailed values of objective function terms");
        details.setRequired(false);
        options.addOption(details);

        Option min_total_usage = new Option("m", "min-total", true, "minimum total clone proportion across all samples (default 0)");
        min_total_usage.setRequired(false);
        options.addOption(min_total_usage);

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
                Main.OUTDIR = cmd.getOptionValue("output");
            }
            // Gurobi automatically uses all available processors

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
            if(cmd.hasOption("objective")){
                String s = cmd.getOptionValue("objective");
                if(s.equals("l0") || s.equals("L0")){
                    Main.OBJECTIVE = Main.Objective.L0;
                } else if(s.equals("l1") || s.equals("L1")){
                    Main.OBJECTIVE = Main.Objective.L1;

                } else if(s.toLowerCase().equals("l0center")){
                    Main.OBJECTIVE = Main.Objective.L0center;
                }
                    else {
                    System.out.println("Objective must be one of the following: L0, L1");
                    System.exit(1);
                }
            }
            if(cmd.hasOption("solver")) {
                String s = cmd.getOptionValue("solver");
                if (s.equals("lpsolve")) {
                    Main.SOLVER = Main.JavaILPSolver.LP_SOLVE;
                } else if (s.equals("gurobi")) {
                    Main.SOLVER = Main.JavaILPSolver.GUROBI;
                } else if (s.equals("glpk")) {
                    Main.SOLVER = Main.JavaILPSolver.GLPK;
                } else {
                    System.out.println("Solver must be one of the following: lpsolve, gurobi, glpk. Other solvers are " +
                            "not supported.");
                    System.exit(1);
                }
            }
            if(cmd.hasOption("verbose")){
                int verbose;
                try{
                    verbose = Integer.parseInt(cmd.getOptionValue("verbose"));
                    if (verbose < 0){
                        throw new NumberFormatException();
                    }
                    Main.SOLVER_VERBOSE = verbose;
                } catch(NumberFormatException e){
                    System.out.println("Solver verbose setting must be a positive integer.");
                    System.exit(1);
                }
            }
            if(cmd.hasOption("timeout")){
                int timeout;
                try{
                    timeout = Integer.parseInt(cmd.getOptionValue("timeout"));
                    if (timeout < 0){
                        throw new NumberFormatException();
                    }
                    Main.SOLVER_TIMEOUT = timeout;
                } catch(NumberFormatException e){
                    System.out.println("Solver timeout must be a positive integer.");
                    System.exit(1);
                }
            }

            if(cmd.hasOption("min-total")){
                double minusage;
                try{
                    minusage = Double.parseDouble(cmd.getOptionValue("min-total"));
                    if (minusage < 0){
                        throw new NumberFormatException();
                    }
                    Main.MIN_TOTAL_USAGE = minusage;
                } catch(NumberFormatException e){
                    System.out.println("Minimum total usage must be non-negative.");
                    System.exit(1);
                }
            }

            Main.PRINT_ANCESTRY_GRAPH = cmd.hasOption("print-graph");
            Main.PRINT_CONFIDENCE_INTERVALS = cmd.hasOption("intervals");
            Main.PRINT_EFFECTIVE_CONFIDENCE = cmd.hasOption("printconf");
            Main.TIMING = cmd.hasOption("time");
            Main.LONGITUDINAL = !cmd.hasOption("nonlongitudinal");
            Main.REMOVE_CNA = cmd.hasOption("remove-columns");
            Main.ENUMERATE_TREES = cmd.hasOption("enumerate");
            Main.PRINT_OBJECTIVE_DETAILS = cmd.hasOption("details");

        } catch(NumberFormatException e){
            System.out.println(e.getMessage());
            formatter.printHelp("calder", options);
            System.exit(1);
        }

        Main.main(args);
    }

}
