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

import blbutil.FileIt;
import blbutil.InputIt;
import blbutil.Utilities;
import java.io.File;
import vcf.AlleleRefGTRecIt;

/**
 * <p>Class {@code VcfInput} process VCF and gzipped VCF input files.</p>
 *
 * <p>Instances of class {@code VcfInput} are not thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class VcfInput {

    private VcfInput() {
        // private constructor to prevent instantiation
    }

    /**
     * Reads VCF data from {@code this.par().in()}, filters the data,
     * and writes the filtered data to {@code this.par().out()}.
     * The output file will be a bref4 file if
     * {@code (par.out().endsWith(".bref4") == true)}.
     * The VCF records will be read from {@code stdin} if
     * {@code (par.in().equals("-") == true)}, and filtered VCF records will
     * be written to {@code stdout} if {@code (par.out().equals("-") == true)}.
     * @param par the bref4 command line parameters
     * @throws NullPointerException if {@code (par == null)}
     */
    public static void processVcf(Bref4Par par) {
        String in  = par.in();
        if (in.endsWith(".vcf.gz")==false
                && in.endsWith(".vcf.bgz")==false
                && in.equals("-")==false) {
            throw new IllegalArgumentException(in);
        }
        if (par.out().endsWith(".bref4")) {
            writeBref4(par);
        }
        else {
            writeVCF(par);
        }
    }

    private static AlleleRefGTRecIt refIt(Bref4Par par) {
        String input = par.in();
        int bufferSize = par.nThreads() << 3;
        if (input.equals("-")) {
            return new AlleleRefGTRecIt(InputIt.fromStdIn(), bufferSize);
        }
        else {
            File inFile = new File(input);
            int nBufferedBlocks = (par.nThreads() << 3);
            FileIt<String> it = InputIt.fromBGZipFile(inFile, nBufferedBlocks);
            return new AlleleRefGTRecIt(it, bufferSize);
        }
    }

    private static void writeBref4(Bref4Par par) {
        try (AlleleRefGTRecIt it = refIt(par);
                Bref4Writer out = new Bref4Writer(par, it.vcfHeader())) {
            while (it.hasNext()) {
                out.write(it.next());
            }
        }
        catch (Throwable t) {
            Utilities.exit(t);
        }
    }

    private static void writeVCF(Bref4Par par) {
        try (AlleleRefGTRecIt refIt = refIt(par)) {
            VcfOutput.writeVCF(refIt, par);
        }
    }
}
