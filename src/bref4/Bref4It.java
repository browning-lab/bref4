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
import java.io.File;
import java.io.IOException;
import java.util.ArrayDeque;
import java.util.Arrays;
import java.util.NoSuchElementException;
import vcf.RefGTRec;
import vcf.Samples;
import vcf.VcfHeader;

/**
 * <p>Class {@code Bref4It} is an iterator that returns VCF records containing
 * phased, non-missing genotypes from a bref4 file.</p>
 *
 * <p>The Java Virtual Machine will terminate with an error message
 * if an I/O error or file format error is detected when reading the
 * bref4 file.</p>
 *
 * <p>Instances of class {@code Bref4It} are not thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class Bref4It implements VcfFileIt<RefGTRec> {

    private final Bref4Reader bref4Reader;
    private final VcfHeader vcfHeader;
    private final Bref4BlockInflater blockInflater;
    private final int maxBlocks;
    private final ArrayDeque<ArrayDeque<RefGTRec>> buffers;

    private ArrayDeque<RefGTRec> buffer;

    /**
     * Constructs a new {@code Bref4It} instance from the specified data.
     * @param par the bref4 analysis parameters
     * @throws IllegalArgumentException if
     * {@code ( (new File(par.in())).exists() == false )}
     * @throws NullPointerException if {@code (par == null)}
     */
    public Bref4It(Bref4Par par) {
        this(new File(par.in()), par.nThreads());
    }

    /**
     * Constructs a new {@code Bref4It} instance from the specified data.
     * @param bref4File the bref4 file
     * @param nThreads the number of computational threads
     * @throws IllegalArgumentException if {@code (bref4File.exists() == false)}
     * @throws IllegalArgumentException if {@code (nThreads < 1)}
     * @throws NullPointerException if {@code (bref4File == null)}
     */
    public Bref4It(File bref4File, int nThreads) {
        if (bref4File.exists()==false) {
            String s = "bref4 file does not exist: " + bref4File;
            throw new IllegalArgumentException(s);
        }
        if (nThreads<1) {
            throw new IllegalArgumentException(String.valueOf(nThreads));
        }
        this.bref4Reader = new Bref4Reader(bref4File);
        Bref4Header bref4Header = bref4Reader.header();
        this.vcfHeader = bref4Header.vcfHeader();
        this.blockInflater = new Bref4BlockInflater(bref4Header);
        this.maxBlocks = nThreads<<4;
        this.buffers = new ArrayDeque<>(maxBlocks);
        readBlocksIntoBuffer();
        this.buffer = buffers.isEmpty() ? new ArrayDeque<>(0) : buffers.removeFirst();
    }

    private void readBlocksIntoBuffer() {
        byte[][] bref4Blocks = bref4Reader.readBlocks(maxBlocks);
        Arrays.stream(bref4Blocks)
                .parallel()
                .map(blockInflater)
                .map(recArray -> new ArrayDeque<>(Arrays.asList(recArray)))
                .forEachOrdered(buffers::add);
    }

    /**
     * Returns {@code true} if the iteration has more elements, and returns
     * {@code false} otherwise.
     * @return {@code true} if the iteration has more elements
     */
    @Override
    public boolean hasNext() {
        // An empty empty buffer should only occur after returning the final bref4 block
        return !buffer.isEmpty();
    }

    /**
     * Returns the next element in the iteration.
     * @return the next element in the iteration
     * @throws NoSuchElementException if the iteration has no more elements
     */
    @Override
    public RefGTRec next() {
        if (hasNext()==false) {
            throw new NoSuchElementException();
        }
        RefGTRec rec = buffer.removeFirst();
        if (buffer.isEmpty()) {
            if (buffers.isEmpty()) {
                readBlocksIntoBuffer();
            }
            assert buffers.isEmpty()==false;
            buffer = buffers.removeFirst();
        }
        return rec;
    }

    @Override
    public void close() {
        try {
            bref4Reader.close();
        } catch (IOException ex) {
            Utilities.exit(ex, "Error closing file");
        }
        buffer.clear();
    }

    @Override
    public String source() {
        return bref4Reader.header().source();
    }

    @Override
    public VcfHeader vcfHeader() {
        return vcfHeader;
    }

    @Override
    public Samples samples() {
        return bref4Reader.header().samples();
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(80);
        sb.append(this.getClass().toString());
        sb.append(' ');
        sb.append('[');
        sb.append(bref4Reader.header().source());
        sb.append(']');
        return sb.toString();
    }
}
