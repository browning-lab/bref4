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
package vcf;

import blbutil.BlockLineReader;
import blbutil.FileIt;
import blbutil.Utilities;
import blbutil.VcfFileIt;
import java.util.ArrayDeque;
import java.util.Arrays;
import java.util.Deque;
import java.util.NoSuchElementException;
import java.util.function.Function;

/**
 * <p>Class {@code AlleleRefGTRecIt} represents  an iterator whose
 * {@code next()} method returns an allele VCF record having
 * phased, non-missing genotypes.
 * </p>
 * <p>Instances of class {@code AlleleRefGTRecIt} are not thread-safe.
 * </p>
 * <p>Methods of this class will terminate the Java Virtual Machine with
 * an error message if an I/O error occurs.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class AlleleRefGTRecIt implements VcfFileIt<RefGTRec> {

    private final VcfHeader vcfHeader;
    private final Function<String, RefGTRec> mapper;

    private final BlockLineReader reader;
    private final Deque<RefGTRec> recBuffer;

    /**
     * Create a new {@code AlleleCodedRefIt} instance from the
     * specified data.
     * @param it an iterator that returns lines of a VCF file
     * @param bufferSize the number of VCF records stored in a buffer
     * @throws IllegalArgumentException if a line of the VCF file returned by
     * {@code it} cannot be parsed
     * @throws IllegalArgumentException if {@code (bufferSize < 1)}
     * @throws NullPointerException if {@code (it == null)}
     */
    public AlleleRefGTRecIt(FileIt<String> it, int bufferSize) {
        if (bufferSize < 1) {
            throw new IllegalArgumentException(String.valueOf(bufferSize));
        }
        String[] vcfHeaderLines = VcfHeader.readVcfHeader(it, null, null)
                .toArray(String[]::new);
        String firstDataLine = it.hasNext() ? it.next() : null;
        if (firstDataLine==null) {
            throw new IllegalArgumentException("Missing VCF data lines (" + it.source() + ")");
        }
        boolean[] isDiploid = VcfHeader.isDiploid(firstDataLine);
        this.vcfHeader = new VcfHeader(it.source(), vcfHeaderLines, isDiploid, null);
        this.mapper = (String s) -> {
            return RefGTRec.alleleRefGTRec(new VcfRecGTParser(vcfHeader, s, false));
        } ;
        this.recBuffer = new ArrayDeque<>(bufferSize);
        int nBlocks = 1;
        this.reader = BlockLineReader.create(it, bufferSize, nBlocks);
        fillRecBuffer(firstDataLine);
    }

    private void fillRecBuffer() {
        fillRecBuffer(null);
    }

    private void fillRecBuffer(String firstDataLine) {
        while (recBuffer.isEmpty()) {
            String[] lines = reader.next();
            if (firstDataLine!=null) {
                lines = combine(firstDataLine, lines);
                firstDataLine = null;
            }
            if (lines==BlockLineReader.SENTINAL) {
                return;
            }
            else {
                RefGTRec[] recs = parseAndFilterLines(lines);
                recBuffer.addAll(Arrays.asList(recs));
            }
        }
    }

    private String[] combine(String firstDataLine, String[] lines) {
        String[] modLines = new String[lines.length + 1];
        modLines[0] = firstDataLine;
        System.arraycopy(lines, 0, modLines, 1, lines.length);
        return modLines;
    }

    private RefGTRec[] parseAndFilterLines(String[] lines) {
        RefGTRec[] recs = null;
        try {
            recs = Arrays.stream(lines)
                    .parallel()
                    .map(mapper)
                    .toArray(RefGTRec[]::new);
        }
        catch (Throwable t) {
            Utilities.exit(t);
        }
        return recs;
    }

    @Override
    public void close() {
        reader.close();
        recBuffer.clear();
    }

    /**
     * Returns {@code true} if the iteration has more elements, and returns
     * {@code false} otherwise.
     * @return {@code true} if the iteration has more elements
     */
    @Override
    public boolean hasNext() {
        return !recBuffer.isEmpty();
    }

    /**
     * Returns the next element in the iteration.
     * @return the next element in the iteration
     * @throws NoSuchElementException if the iteration has no more elements
     */
    @Override
    public RefGTRec next() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }
        RefGTRec first = recBuffer.removeFirst();
        if (recBuffer.isEmpty()) {
            fillRecBuffer();
        }
        return first;
    }

    /**
     * The {@code remove} method is not supported by this iterator.
     * @throws UnsupportedOperationException if this method is invoked
     */
    @Override
    public void remove() {
        throw new UnsupportedOperationException(this.getClass().toString());
    }

    @Override
    public String source() {
        return reader.source();
    }

    @Override
    public Samples samples() {
        return vcfHeader.samples();
    }

    @Override
    public VcfHeader vcfHeader() {
        return vcfHeader;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(80);
        sb.append(this.getClass().toString());
        sb.append(" : ");
        sb.append(reader.source());
        return sb.toString();
    }
}
