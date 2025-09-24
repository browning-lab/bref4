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

import blbutil.BGZIPOutputStream;
import blbutil.BooleanArray;
import blbutil.FileUtil;
import blbutil.Utilities;
import blbutil.VcfFileIt;
import ints.UnsignedByteArray;
import java.io.BufferedOutputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.function.IntFunction;
import java.util.stream.IntStream;
import vcf.RefGTRec;
import vcf.VcfHeader;
import vcf.VcfWriter;

/**
 * <p>Class {@code VcfInput} process VCF and gzipped VCF input files.</p>
 *
 * <p>Instances of class {@code VcfInput} are not thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class VcfOutput {

    private VcfOutput() {
        // private constructor to prevent instantiation
    }

    /**
     * Writes the specified VCF records to {@code par.out()}.  The VCF
     * records will be written to {@code stdout} if
     * {@code (par.out().equals("-") == true)}.
     * @param it an iterator that returns VCF records
     * @param par the bref4 command line parameters
     * @throws NullPointerException if {@code ((it == null) || (par == null))}
     */
    public static void writeVCF(VcfFileIt<RefGTRec> it, Bref4Par par) {
        int nThreads = Math.max(par.nThreads(), 1);
        String out = par.out();
        try (OutputStream os = outputStream(out)) {
            VcfHeader vcfHeader = it.vcfHeader();
            int nHaps = vcfHeader.nSamples() << 1;
            int nBufferedRecs = nBufferedRecs(nHaps, nThreads);
            boolean bgzip = out.endsWith(".gz") || out.endsWith(".bgz");
            writeVcfHeader(par, vcfHeader, bgzip, os);
            while (it.hasNext()) {
                RefGTRec[] recs = readRecArrays(it, nBufferedRecs);
                int nSteps = Math.min(recs.length, (nThreads<<2));
                int[] stepEnds = startsAndEnds(0, recs.length, nSteps);
                IntFunction<UnsignedByteArray> stepToByteArray =
                        i -> recordsToBytes(recs, null, stepEnds[i-1], stepEnds[i], bgzip);
                UnsignedByteArray[] bytes = IntStream.range(1, stepEnds.length)
                    .parallel()
                    .mapToObj(stepToByteArray)
                    .toArray(UnsignedByteArray[]::new);
                writeBytes(out, bytes, os);
            }
            if (bgzip) {
                BGZIPOutputStream.writeEmptyBlock(os);
            }
        }
        catch (Throwable t) {
            Utilities.exit(t);
        }
    }

    private static PrintWriter baosPrintWriter(ByteArrayOutputStream baos,
            boolean bgzip) {
        if (bgzip) {
            return new PrintWriter(new BGZIPOutputStream(baos, false));
        }
        else {
            return new PrintWriter(baos);
        }
    }

    private static OutputStream outputStream(String outFile) {
        if (outFile.equals("-")) {
            return new BufferedOutputStream(System.out);
        }
        else {
            return FileUtil.bufferedOutputStream(new File(outFile));
        }
    }

    private static int nBufferedRecs(int nHaps, int nThreads) {
        long maxAlleles = 1L<<32;
        long nBufferedRecs = nThreads<<10;
        if (nBufferedRecs*nHaps < maxAlleles) {
            return (int) nBufferedRecs;
        }
        else {
            nBufferedRecs = (long) Math.ceil(maxAlleles/nHaps);
            return (int) Math.min(nBufferedRecs, Integer.MAX_VALUE);
        }
    }

    private static void writeVcfHeader(Bref4Par par, VcfHeader vcfHeader,
             boolean bgzip, OutputStream os) {
        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        try (PrintWriter baosOut = baosPrintWriter(baos, bgzip)) {
            VcfWriter.copyMetaInfoLines(vcfHeader, baosOut);
            boolean quoteValue = true;
            baosOut.println(VcfHeader.metaInfoLine("bref4Command",
                    Bref4Par.commandAndVersion(), quoteValue));
            VcfWriter.writeHeaderLine(vcfHeader.samples().ids(), baosOut);
        }   // closes baosOut to flush buffer (if there is one) to baos
        try {
            baos.writeTo(os);
        } catch (IOException e) {
            writeError(par.out(), e);
        }
    }

    private static RefGTRec[] readRecArrays(VcfFileIt<RefGTRec> it, int bufferSize) {
        ArrayList<RefGTRec> buffer = new ArrayList<>(bufferSize);
        for (int j=0; j<bufferSize && it.hasNext(); ++j) {
            buffer.add(it.next());
        }
        return buffer.toArray(RefGTRec[]::new);
    }

    private static int[] startsAndEnds(int start, int end, int nSteps) {
        int step = ((end - start) + (nSteps - 1))/nSteps;
        int[] startsAndEnds = new int[nSteps + 1];
        int index = 0;
        for (int m=start; m<end; m+=step) {
            startsAndEnds[index++] = m;
        }
        startsAndEnds[index] = end;
        return startsAndEnds;
    }

    private static UnsignedByteArray recordsToBytes(RefGTRec[] recs,
            BooleanArray isHaploid, int start, int end, boolean bgzip) {
        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        try (PrintWriter baosOut = baosPrintWriter(baos, bgzip)) {
            IntStream.range(start, end)
                    .forEach(j -> baosOut.println(recs[j].toVcfRecord(isHaploid)));
        }
        return new UnsignedByteArray(baos);
    }

    private static void writeBytes(String out, UnsignedByteArray[] bytes,
            OutputStream os) {
        try {
            for (UnsignedByteArray uba : bytes) {
                uba.write(os);
            }
        } catch (IOException e) {
            writeError(out, e);
        }
    }

    private static void writeError(String out, Exception e) {
        StringBuilder sb = new StringBuilder();
        sb.append("Error writing to \"");
        sb.append(out);
        sb.append('"');
        if (out.equals("-")) {
            sb.append(" (stdout)");
        }
        Utilities.exit(e, sb.toString());
    }
}
