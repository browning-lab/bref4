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

import java.util.Locale;

/**
 * <p>Class {@code BrefMain} performs data filtering and converts between
 * files in VCF and bref4 format.</p>
 *
 * <p>Instances of class {@code Bref4Main} are not thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class Bref4Main {

    static final String EXECUTABLE = "bref4.jar";
    static final String VERSION = "0.1 (alpha release)";
    static final String COPYRIGHT = "Copyright (C) 2025 Brian L. Browning";
    static final String COMMAND = "java -jar " + EXECUTABLE;

    private Bref4Main() {
        // private constructor to prevent instantiation
    }

    /**
     * The {@code main()} method is the entry point to the bref4 program.
     * See the {@code Bref4Par.usage()} method for usage instructions.
     * The Java Virtual Machine will exit with an error message if there
     * is an error in the command line arguments.
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        Locale.setDefault(Locale.US);
        Bref4Par par = getPar(args);
        System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism",
                String.valueOf(par.nThreads()));

        String inFile = par.in();
        if (inFile.endsWith(".vcf.gz") || inFile.endsWith(".vcf.bgz") || inFile.equals("-")) {
            VcfInput.processVcf(par);
        }
        else if (inFile.endsWith(".bref4")) {
            Bref4Input.processBref4(par);
        }
        else {
            assert false; // Bref4Par constructor performs parameter validation
        }
    }

    private static Bref4Par getPar(String[] args) {
        if (args.length==0
                || (args.length==1 && args[0].toLowerCase().equals("help"))) {
            Bref4Par.printHelp();
            System.exit(0);
        }
        return new Bref4Par(args);
    }
}
