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

import beagleutil.ChromIds;
import blbutil.Const;
import blbutil.Utilities;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.DataInputStream;
import java.io.DataOutput;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * <p>Class {@code Bref4Index} stores a binary reference format version 4
 * (bref4) index file.</p>
 *
 * <p>Instances of class {@code Bref4Index} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class Bref4Index {

    // row (first dimension) is chromosome index
    private final int[][] startPos;
    private final int[][] endPos;
    private final long[][] offset;
    private final int nBlocks;

    /**
     * Constructs a new {@code Bref4Index} instance from the specified
     * bref4 index file. The constructor will print an error
     * message and terminate the Java Virtual Machine if an I/O error occurs or
     * a file format error is detected.
     * @param bref4File a bref4 index file
     *
     * @throws NullPointerException if {@code (bref4File == null)}
     */
    public Bref4Index(File bref4File) {
        if (bref4File==null) {
            throw new NullPointerException("bref4File==null");
        }
        ArrayList<ArrayList<IndexRecord>> index = readIndex(bref4File);
        int nChrom = index.size();
        int[][] startPos0 = new int[nChrom][];
        int[][] endPos0 = new int[nChrom][];
        long[][] offsets0 = new long[nChrom][];
        int blockCnt = 0;
        for (int j=0; j<nChrom; ++j) {
            ArrayList<IndexRecord> recList = index.get(j);
            int n = recList.size();
            blockCnt += n;
            int[] starts = new int[n];
            int[] ends = new int[n];
            long[] offsets = new long[n];
            for (int k=0; k<n; ++k) {
                IndexRecord indexData = recList.get(k);
                starts[k] = indexData.startPos;
                ends[k] = indexData.endPos;
                offsets[k] = indexData.offset;
            }
            startPos0[j] = starts;
            endPos0[j] = ends;
            offsets0[j] = offsets;
        }
        this.startPos = startPos0;
        this.endPos = endPos0;
        this.offset = offsets0;
        this.nBlocks = blockCnt;
    }

    private static ArrayList<ArrayList<IndexRecord>> readIndex(File bref4File) {
        ChromIds chromIds = ChromIds.instance();
        ArrayList<ArrayList<IndexRecord>> indexData = new ArrayList<>();
        long lastByteOffset = bref4File.length() - Long.BYTES;
        try (RandomAccessFile raf = new RandomAccessFile(bref4File, "r")) {
            raf.seek(lastByteOffset);
            long indexOffset = raf.readLong();
            raf.seek(indexOffset);
            long offset0;
            while ((offset0 = raf.readLong()) != -1) {
                String chrId = raf.readUTF();
                int startPos0 = raf.readInt();
                int endPos0 = raf.readInt();
                IndexRecord indexRecord = new IndexRecord(startPos0, endPos0, offset0);
                int chrIndex = chromIds.getIndex(chrId);
                while (chrIndex>=indexData.size()) {
                    indexData.add(new ArrayList<>());
                }
                indexData.get(chrIndex).add(indexRecord);
            }
        }
        catch (IOException ex) {
            Utilities.exit(ex, "Error reading bref file index, " + bref4File
                    + ". Exiting program");
        }
        return indexData;
    }

    private static class IndexRecord {

        private final int startPos;
        private final int endPos;
        private final long offset;

        public IndexRecord(int startPos, int endPos, long offset) {
            this.startPos = startPos;
            this.endPos = endPos;
            this.offset = offset;
        }
    }

    /**
     * Return the number of bref4 blocks.
     * @return the number of bref4 blocks
     */
    public int nBlocks() {
        return nBlocks;
    }

    /**
     * Return the number of bref4 blocks on the specified chromosome.
     * @param chromIndex the chromosome index
     * @return the number of bref4 blocks on the specified chromosome
     * @throws IndexOutOfBoundsException if {@code (chromIndex < 0)}
     */
    public int nBlocks(int chromIndex) {
        if (chromIndex>=startPos.length) {
            return 0;
        }
        else {
            return startPos[chromIndex].length;
        }
    }

    /**
     * Return the position of the first record in the specified bref4 block
     * @param chromIndex the chromosome index
     * @param block the bref4 block index
     * @return the position of the first record in the specified bref4 block
     * @throws IndexOutOfBoundsException if {@code (chromIndex < 0)}
     * @throws IndexOutOfBoundsException if
     * {@code (block < 0) || (block >= this.nBlocks(chromIndex))}
     */
    public int startPos(int chromIndex, int block) {
        return startPos[chromIndex][block];
    }

    /**
     * Return the position of the last record in the specified bref4 block
     * @param chromIndex the chromosome index
     * @param block the bref4 block index
     * @return the position of the last record in the specified bref4 block
     * @throws IndexOutOfBoundsException if {@code (chromIndex < 0)}
     * @throws IndexOutOfBoundsException if
     * {@code (block < 0) || (block >= this.nBlocks(chromIndex))}
     */
    public int endPos(int chromIndex, int block) {
        return endPos[chromIndex][block];
    }

    /**
     * Return the byte offset of the specified bref4 block
     * @param chromIndex the chromosome index
     * @param block the bref4 block index
     * @return the byte offset of the specified bref4 block
     * @throws IndexOutOfBoundsException if {@code (chromIndex < 0)}
     * @throws IndexOutOfBoundsException if
     * {@code (block < 0) || (block >= this.nBlocks(chromIndex))}
     */
    public long offset(int chromIndex, int block) {
        return offset[chromIndex][block];
    }

    /**
     * Returns the index of the first bref4 block whose genomic interval
     * includes the specified genomic position. If no bref4 block
     * includes the position, the index of the first bref4 block for
     * the chromosome whose genomic interval begins after the specified
     * position is returned. If there are no bref4 blocks for the chromosome
     * that contain or begin after the specified
     * position, {@code this.nBlocks(chrom)} is returned.
     * @param chromIndex a chromosome index
     * @param pos a position
     * @return the index of the of a bref4 block whose genomic interval
     * includes the specified genomic position
     * @throws IndexOutOfBoundsException if {@code (chromIndex < 0)}
     */
    public int block(int chromIndex, int pos) {
        if (chromIndex >= endPos.length || endPos[chromIndex].length==0) {
            return 0;
        }
        int x = Arrays.binarySearch(endPos[chromIndex], pos);
        if (x>=0) {
            while (x>0 && endPos[chromIndex][x-1]==pos) {
                --x;
            }
            return x;
        }
        else {
            return -x-1;
        }
    }

    /**
     * Returns the index of the first bref4 block whose start position is
     * greater than the specified genomic position. If there are no bref4
     * blocks for the chromosome whose start position is greater than
     * the specified position, {@code this.nBlocks(chrom)} is returned.
     * @param chromIndex a chromosome index
     * @param pos a position
     * @return the index of the first bref4 block after the block that contains
     * the specified position, or, if no such bref4 block exists, the first
     * bref4 block whose genomic interval begins after specified genomic
     * position
     * @throws IndexOutOfBoundsException if {@code (chromIndex < 0)}
     */
    public int nextBlock(int chromIndex, int pos) {
        if (chromIndex >= startPos.length || startPos[chromIndex].length==0) {
            return 0;
        }
        int x = Arrays.binarySearch(startPos[chromIndex], pos);
        if (x>=0) {
            ++x;
            while ((x<startPos[chromIndex].length) &&
                    (startPos[chromIndex][x]==pos)) {
                ++x;
            }
            return x;
        }
        else {
            return -x-1;
        }
    }

    /**
     * Write the specified block data to a bref4 index file. The Java virtual
     * machine will exit with an error if an {@code IOException} is thrown.
     * No bytes are written to the index file if {@code (bytes.length == 0)}
     * The contract for this method is undefined if {@code bytes} is not a
     * valid bref4 block.
     * @param offset the number of bytes in the bref4 file preceding the
     * bref4 block data
     * @param bytes the bytes in the bref4 block
     * @param dos the interface to the output index file
     * @throws NullPointerException if
     * {@code (bytes == null) || (indexOut == null)}
     */
    public static void writeBref4BlockToIndex(long offset, byte[] bytes,
            DataOutputStream dos) {
        if (bytes.length>0) {
            try {
                DataInputStream dis = new DataInputStream(new ByteArrayInputStream(bytes));
                dis.readInt(); // nRecords
                int endPos = dis.readInt();
                dis.readByte(); // nMaps
                String chrom = dis.readUTF();
                int startPos = dis.readInt();
                writeBref4BlockToIndex(offset, chrom, startPos, endPos, dos);
            } catch (IOException ex) {
                Utilities.exit(ex, "Error parsing bref4 block");
            }
        }
    }

    /**
     * Write the specified block data to a bref4 index file.  The Java virtual
     * machine will exit with an error if an {@code IOException} is thrown.
     * @param offset the number of bytes in the bref4 file preceding the
     * bref4 block data
     * @param chrom the chromosome
     * @param startPos the position of the first VCF record in the block
     * @param endPos the position of the last VCF record in the block
     * @param indexOut the interface to the output index file
     * @throws NullPointerException if
     * {@code (chrom == null) || (indexOut == null)}
     */
    public static void writeBref4BlockToIndex(long offset, String chrom,
            int startPos, int endPos, DataOutput indexOut) {
        try {
            indexOut.writeLong(offset);
            indexOut.writeUTF(chrom);
            indexOut.writeInt(startPos);
            indexOut.writeInt(endPos);
        } catch (IOException ex) {
            Utilities.exit(ex, "Error writing block to bref4 index");
        }
    }

    /**
     * Write the bref4 index to the specified output stream. The Java virtual
     * machine will exit with an error if an {@code IOException} is thrown.
     * @param indexBytes a byte array output stream storing the bref4 index
     * @param bytesWritten the number of bytes previously written to
     * {@code bref4Out}
     * @param bref4Out the output stream to the bref4 file
     * @throws NullPointerException if
     * {@code (indexBytes == null) || (bref4Out == null)}
     */
    public static void writeBref4Index(ByteArrayOutputStream indexBytes,
            long bytesWritten, DataOutputStream bref4Out) {
        try (DataOutputStream indexOut = new DataOutputStream(indexBytes)) {
            indexOut.writeLong(-1);  // sentinal signalling end of bref4 index
            indexOut.writeLong(bytesWritten);   // file offset of index
            indexBytes.writeTo(bref4Out);
        } catch (IOException ex) {
            Utilities.exit(ex, "Error appending bref4 index to bref4 file");
        }
    }

    /**
     * Returns a string description of {@code this}.  The exact details of
     * the description are unspecified and subject to change.
     * @return a string description of {@code this}
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(80);
        String[] chromIds = ChromIds.instance().ids();
        sb.append("BLOCK\tCHROM\tSTART\tEND\tOFFSET\tDIFF");
        sb.append(Const.nl);
        int blockCnt = 0;
        for (int j=0; j<offset.length; ++j) {
            for (int k=0, n=offset[j].length; k<n; ++k) {
                sb.append(++blockCnt);
                sb.append(Const.tab);
                sb.append(chromIds[j]);
                sb.append(Const.tab);
                sb.append(startPos(j,k));
                sb.append(Const.tab);
                sb.append(endPos(j,k));
                sb.append(Const.tab);
                sb.append(offset(j,k));
                sb.append(Const.tab);
                sb.append(offset(j,k) - offset(j, Math.max(0, k-1)));
                sb.append(Const.nl);
            }
        }
        return sb.toString();
    }
}