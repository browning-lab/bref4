/*
 * Copyright 2025 Brian L. Browning
 *
 * This file is part of the bref4 program.
 *
 * Licensed under the Apache License, Version 2.0 (the License);
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http:/www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an AS IS BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package bref4;

import blbutil.Const;
import blbutil.Utilities;
import blbutil.Validate;
import java.io.File;
import java.util.Map;

/**
 * <p>Class {@code Bref4Par} represents the command line parameters for the
 * {@code bref4} program.</p>
 *
 * <p>Instances of class {@code Bref4Par} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class Bref4Par {

    private final String[] args;

    private final String in;
    private final String out;
    private final int nthreads;

    // Undocumented parameters
    private final int bitsPerLevel;
    private final int maxNonmajor;

    private static final int DEF_NTHREADS = Runtime.getRuntime().availableProcessors();

    private static final int DEF_BITS_PER_LEVEL = 2;
    static final int DEF_MAX_NONMAJOR = Integer.MAX_VALUE; // forces value to be calculated from nHaps

    /**
     * Constructs an {@code Bref4Par} instance from the command line
     * parameters. The Java Virtual Machine will exit with an error message
     * if the command line contains an invalid parameter.
     *
     * @param args the command line parameters
     * @throws IllegalArgumentException if the command line parameter is
     * incorrectly specified
     * @throws NullPointerException if {@code (args == null)}
    */
    public Bref4Par(String[] args) {
        if (args.length==0 || args[0].equalsIgnoreCase("help")) {
            System.out.println(programInfo());
            System.out.println();
            System.out.println(usage());
            System.exit(0);
        }
        this.args = args.clone();
        int IMAX = Integer.MAX_VALUE;
        Map<String, String> argsMap = Validate.argsToMap(args, '=');

        this.in = Validate.stringArg("in", argsMap, false, null, null);
        this.out = Validate.stringArg("out", argsMap, false, null, null);
        validateInputAndOutputFiles(in, out);

        this.nthreads = Validate.intArg("nthreads", argsMap, false, DEF_NTHREADS, 1, IMAX);

        // Undocumented parameters
        this.bitsPerLevel = Validate.intArg("bits-per-level", argsMap, false, DEF_BITS_PER_LEVEL, 1, IMAX);
        this.maxNonmajor = Validate.intArg("max-nonmajor", argsMap, false, DEF_MAX_NONMAJOR, 0, IMAX);
        Validate.confirmEmptyMap(argsMap);
    }

    /**
     * Prints the program information and usage instructions to {@code stdout}.
     */
    public static void printHelp() {
        System.out.println();
        System.out.println(programInfo());
        System.out.println();
        System.out.println(Bref4Par.usage());
    }

    /**
     * Return a string representation of the program version and description.
     * The exact details of the returned string are unspecified and subject
     * to change.
     * @return a string representation of the program version and description
     */
    public static String programInfo() {
        return "bref4 version "
                + Bref4Main.VERSION
                + Const.nl
                + Const.nl
                + "The bref4 program compresses and filters phased sequence data.";
    }

    /**
     * Returns the command line.
     * @return the command line
     */
    public static String commandLine() {
        return ProcessHandle
                 .current()
                 .info()
                 .commandLine()
                 .orElse("");
    }

    /**
     * Returns a string representation of the command line and bref4 version.
     * The exact details of the representation are unspecified and subject
     * to change.
     * @return a string representation of the command line and bref4 version
     */
    public static String commandAndVersion() {
        StringBuilder sb = new StringBuilder(commandLine());
        sb.append("  # bref4.jar (version ");
        sb.append(Bref4Main.VERSION);
        sb.append(')');
        return sb.toString();
    }

    private static boolean isValidInputOrOutput(String filename) {
        File file = new File(filename);
        return file.isDirectory()==false &&
                (isVcfFormat(filename) || filename.endsWith(".bref4"));
    }

    /**
     * Returns {@code true} if the specified string equals "-"
     * or ends with ".vcf", ".vcf.gz", or ".vcf.bgz", and returns
     * {@code false} otherwise.
     * @param string an string
     * @return if the specified string equals "-" or ends with ".vcf",
     * ".vcf.gz", or ".vcf.bgz"
     */
    public static boolean isVcfFormat(String string) {
        return (string.equals("-")
                || string.endsWith(".vcf.gz")
                || string.endsWith(".vcf.bgz")
                || string.endsWith(".vcf"));
    }

    private static void appendCommandLineAndUsage(StringBuilder sb) {
        sb.append(" \"");
        sb.append(commandLine());
        sb.append('"');
        sb.append(Const.nl);
        sb.append(Const.nl);
        sb.append(usage());
    }

    /**
     * Returns a string with bref4 usage instructions.
     * The format of the returned string is unspecified and subject to change.
     * @return a string with bref4 usage instructions
     */
    public static String usage() {
        String nl = Const.nl;
        return    "Usage: "  + nl
                + "  " + Bref4Main.COMMAND + " [parameters]" + nl
                + nl
                + "Input and output file parameters:" + nl
                + "  in=[input file]                                        (required)" + nl
                + "  out=[output file]                                      (required)" + nl
                + nl
                + "  The filename suffix must indicate the file type: " + nl
                + nl
                + "    uncompressed VCF (\"*.vcf\")" + nl
                + "    gzip-compressed VCF (\"*.vcf.gz\" or \"*.vcf.bgz\")" + nl
                + "    bref4 (\"*.bref4\")" + nl
                + nl
                + "    Replace \"[input file]\" with \"-\" to read an uncompressed VCF file from stdin" + nl
                + "    Replace \"[output file]\" with \"-\" to write an uncompressed VCF file to stdout" + nl
                + nl
                + "General parameters:" + nl
                + "  nthreads=<number of threads>                           (default: all CPU cores)" + nl
                + nl;
    }

    private static void validateInputAndOutputFiles(String in, String out) {
        if (in==null) {
            exitWithMissingFileError("in");
        }
        if (out==null) {
            exitWithMissingFileError("out");
        }
        if (isValidInputOrOutput(in)==false) {
            exitWithInvalidFileError("in", in);
        }
        if (isValidInputOrOutput(out)==false) {
            exitWithInvalidFileError("out", out);
        }
        if (in!=null && in.equals("-")==false && in.equals(out)) {
            exitWithFileClobberError(in);
        }
    }

    private static void exitWithMissingFileError(String inOrOut) {
        assert inOrOut.equals("in") || inOrOut.equals("out");
        StringBuilder sb = new StringBuilder(200);
        sb.append("Error: Missing ");
        sb.append(inOrOut);
        sb.append("put file. The ");
        sb.append(inOrOut);
        sb.append("put file is specified with the \"");
        sb.append(inOrOut);
        sb.append("=\" parameter.");
        sb.append(Const.nl);
        sb.append(Const.nl);
        sb.append("Command:");
        sb.append(Const.nl);
        appendCommandLineAndUsage(sb);
        Utilities.exit(sb.toString());
    }

    private static void exitWithInvalidFileError(String inOrOut, String value) {
        assert inOrOut.equals("in") || inOrOut.equals("out");
        StringBuilder sb = new StringBuilder(200);
        sb.append("Error: Invalid ");
        sb.append(inOrOut);
        sb.append("put file \"");
        sb.append(value);
        sb.append("\". The ");
        sb.append(inOrOut);
        sb.append("put file must be'-' (for std");
        sb.append(inOrOut);
        sb.append(')');
        sb.append(Const.nl);
        sb.append("or a filename that ends in \".vcf.gz\", \".vcf.bgz\", \".vcf\", or \".bref4\"");
        sb.append(Const.nl);
        sb.append(Const.nl);
        sb.append("Command:");
        sb.append(Const.nl);
        appendCommandLineAndUsage(sb);
        Utilities.exit(sb.toString());
    }

    private static void exitWithFileClobberError(String filename) {
        StringBuilder sb = new StringBuilder(200);
        sb.append("Error: input_file and output_file are the same file: \"");
        sb.append(filename);
        sb.append('"');
        sb.append(Const.nl);
        sb.append(Const.nl);
        sb.append("Command:");
        sb.append(Const.nl);
        appendCommandLineAndUsage(sb);
        Utilities.exit(sb.toString());
    }

    /**
     * Returns the command line arguments.
     * @return the command line arguments
     */
    public String[] args() {
        return args.clone();
    }

    /**
     * Returns the input VCF filename or {@code '-'} if the input VCF file is
     * read from stdin.
     * @return the input VCF filename or {@code '-'} if the input VCF file is
     * read from stdin
     */
    public String in() {
        return in;
    }

    /**
     * Returns the output VCF filename or {@code '-'} if the output VCF file
     * is written to stdout.
     * @return the output VCF filename or {@code '-'} if the output VCF file
     * is written to stdout
     */
    public String out() {
        return out;
    }

    /**
     * Returns the {@code nthreads} parameter.
     * @return the {@code nthreads} parameter
     */
    public int nThreads() {
        return nthreads;
    }

    /**
     * Returns the {@code bits-per-level} parameter.
     * @return the {@code bits-per-level} parameter
     */
    public int bitsPerLevel() {
        return bitsPerLevel;
    }

    /**
     * Returns the {@code max-nonmajor} parameter.
     * @return the {@code max-nonmajor} parameter
     */
    public int maxNonmajor() {
        return maxNonmajor;
    }
}
