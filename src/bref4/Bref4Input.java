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

import blbutil.Utilities;
import blbutil.VcfFileIt;
import java.io.ByteArrayOutputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.concurrent.atomic.AtomicLong;
import vcf.RefGTRec;

/**
 * <p>Class {@code Bref4Input} converts files in bref version 4 format
 * into VCF format.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class Bref4Input {

    private Bref4Input() {
        // private constructor to prevent instantiation
    }

    /**
     * Reads bref4 data from {@code this.par().in()}, filters the data,
     * and writes the filtered data to {@code this.par().out()}.
     * @param par the bref4 command line parameters
     * @throws NullPointerException if {@code (par == null)}
     */
    public static void processBref4(Bref4Par par) {
        if (par.out().endsWith(".bref4")) {
            writeBref4(par);
        }
        else {
            try (VcfFileIt<RefGTRec> brefIt = new Bref4It(par)) {
                VcfOutput.writeVCF(brefIt, par);
            }
        }
    }

    private static void writeBref4(Bref4Par par) {
        int maxBlocks = par.nThreads() << 4;
        AtomicLong bytesWritten = new AtomicLong(); // use AtomicLong because brefOut.size() may overflow
        try (Bref4Reader reader = new Bref4Reader(par);
                DataOutputStream brefOut = Bref4Utils.dataOutputStream(par.out());
                ByteArrayOutputStream indexBytes = new ByteArrayOutputStream();
                DataOutputStream indexOut = new DataOutputStream(indexBytes)) {
            Bref4Writer.writeBref4Header(reader.header(), brefOut);
            bytesWritten.addAndGet(brefOut.size());
            byte[][] inputBlocks;
            do {
                inputBlocks = reader.readBlocks(maxBlocks);
                Arrays.stream(inputBlocks)
                    .parallel()
                    .filter(ba -> (ba.length>0))
                    .forEachOrdered(ba -> {
                        try {
                            brefOut.writeInt(ba.length);
                            brefOut.write(ba);
                            long offset = bytesWritten.getAndAdd(Integer.BYTES + ba.length);
                            Bref4Index.writeBref4BlockToIndex(offset, ba, indexOut);
                        }
                        catch (IOException ex) {
                            Utilities.exit(ex, "Error writing file");
                        }
                    });
            } while (inputBlocks[inputBlocks.length-1].length>0);
            brefOut.writeInt(0);
            long fileOffset = bytesWritten.addAndGet(Integer.BYTES);
            Bref4Index.writeBref4Index(indexBytes, fileOffset, brefOut);
        } catch (IOException ex) {
            Utilities.exit(ex, "Error writing file");
        }
    }
}
