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
import java.io.DataInputStream;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;

/**
 * <p>Class {@code Bref4Reader} contains methods for reading data in a bref4 file.</p>
 *
 * <p>Instances of class {@code Bref4Reader} are not thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class Bref4Reader implements AutoCloseable {

    private static final String READ_ERR = "Error reading file";
    private static final byte[] EMPTY_BLOCK = new byte[0];

    private final Bref4Header bref4Header;
    private DataInputStream din;
    private int currentBlock = 0;
    private final int endBlock = Integer.MAX_VALUE;   // exclusive end

    /**
     * Constructs a new {@code Bref4Reader} instance from the specified data.
     * The constructor will print an error message and terminate the
     * Java Virtual Machine if an I/O error occurs or a file format error
     * is detected.
     * @param par the bref4 command line parameters
     * @throws NullPointerException if {@code (par == null)}
     */
    public Bref4Reader(Bref4Par par) {
        this(new File(par.in()));
    }

    /**
     * Constructs a new {@code Bref4Reader} instance from the specified data.
     * The constructor will print an error message and terminate the
     * Java Virtual Machine if an I/O error occurs or a file format error
     * is detected.
     * @param bref4File the bref4File that will be read
     * @throws NullPointerException if {@code (bref4File == null)}
     */
    public Bref4Reader(File bref4File) {
        String source = bref4File.toString();
        this.currentBlock = 0;
        this.din = Bref4Utils.dataInputStream(bref4File);
        this.bref4Header = new Bref4Header(source, din);
    }

    /**
     * Reads and returns the next bref4 block.
     * This method first reads and discards the integer storing the number of
     * bytes in the bref4 block and then reads and returns the bref4 block.
     * The method returns an array of length 0 if {@code this.close()} has
     * previously been invoked or if all bref4 blocks have been read.
     *
     * @return the sequence of bytes in the next bref4 block
     */
    public byte[] readBlock() {
        if (currentBlock<endBlock) {
            try {
                int nBytes = din.readInt();
                if (nBytes==0) {
                    currentBlock = endBlock;
                    return EMPTY_BLOCK;
                }
                else {
                    byte[] bytes = new byte[nBytes];
                    din.readFully(bytes);
                    return bytes;
                }
            } catch (IOException ex) {
                Utilities.exit(ex, READ_ERR);
            }
        }
        return EMPTY_BLOCK;
    }

    /**
     * Returns an array of bref4 blocks obtained by invoking
     * the {@code this.readBlock()} method {@code n} times. If there are
     * fewer than {@code n} bref4 blocks remaining, then the
     * returned array will contain the remaining bref4 blocks followed by
     * a {@code byte[]} array of length 0.
     *
     * @param n the number of bref4 blocks to read
     * @return an array of bref4 blocks
     * @throws NegativeArraySizeException if {@code nBlocks < 0}
     */
    public byte[][] readBlocks(int n) {
        byte[][] blocks = new byte[n][];
        for (int j=0; j<n; ++j) {
            blocks[j] = readBlock();
            if (blocks[j].length==0) {
                blocks = Arrays.copyOf(blocks, j+1);
                break;
            }
        }
        return blocks;
    }

    /**
     * Returns the bref4 file header.
     * @return the bref4 file header
     */
    public Bref4Header header() {
        return bref4Header;
    }

    /**
     * Closes the data input stream associated with this instance.
     * @throws IOException if an I/O exception occurs
     */
    @Override
    public void close() throws IOException {
        currentBlock = endBlock;
        din.close();
    }

    /**
     * Returns a string description of {@code this} instance.  The exact
     * details ofvthe description are unspecified and subject to change.
     * @return a string description of {@code this} instance
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(80);
        sb.append(this.getClass().toString());
        return sb.toString();
    }
}