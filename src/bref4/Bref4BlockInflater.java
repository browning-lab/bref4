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
import blbutil.Utilities;
import ints.IndexArray;
import ints.IntArray;
import java.io.ByteArrayInputStream;
import java.io.DataInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.function.Function;
import vcf.Marker;
import vcf.RefGTRec;
import vcf.Samples;
import vcf.MapRefGTRec;
import vcf.IntArrayRefGTRec;

/**
 * <p>Class {@code Bref4BlockInflater} converts a bref4 block into
 * a {@code RefGTRec[]} array.
 *
 * <p>Instances of class {@code Bref4BlockInflater} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class Bref4BlockInflater implements Function<byte[], RefGTRec[]> {

    private static final String READ_ERR = "Error: invalid bref4 block";

    private final Bref4Header bref4Header;
    private final Samples samples;
    private final int nHaps;

    /**
     * Constructs a new {@code Bref4BlockInflater} instance from the
     * specified data.
     * @param bref4Header the bref4 file header
     * @throws NullPointerException if {@code (bref4Header == null)}
     */
    public Bref4BlockInflater(Bref4Header bref4Header) {
        this.bref4Header = bref4Header;
        this.samples = bref4Header.samples();
        this.nHaps = bref4Header.filteredHapIndices().size();
    }

    /**
     * Returns the bref4 file header.
     * @return the bref4 file header
     */
    public Bref4Header bref4Header() {
        return bref4Header;
    }

    /**
     * Converts the specified bref4 block byte sequence into a
     * {@code RefGTRec[]} array and returns the
     * {@code RefGTRec[]} array. Returns {@code RefGTRec[0]} if
     * {@code (bytes.length == 0)}.  The contract
     * for this method is undefined if {@code (bytes.length > 0)} and the
     * number of haplotypes in the bref4 block is not equal to
     * {@code this.bref4Header().filteredHapIndices().size()}.  The Java
     * virtual machine will exit with an error message if the specified
     * {@code bref4} block does not conform to the {@code bref4} block
     * specification.
     * @param bytes a bref4 block byte sequence
     * @return the VCF records in the block as a {@code RefGTRec[]} array
     * @throws NullPointerException if {@code (bytes == null)}
     */
    @Override
    public RefGTRec[] apply(byte[] bytes) {
        if (bytes.length==0) {
            return new RefGTRec[0];
        }
        // Permit maps[0].size()<maps[0].valueSize() due to sample filtering
        int[] mapBuffer = new int[bref4Header.invFilteredHapIndices().size()];
        ArrayList<RefGTRec> records = new ArrayList<>();
        DataInputStream dis = new DataInputStream(new ByteArrayInputStream(bytes));
        try {
            int nRecs = dis.readInt();
            int lastRecPos = dis.readInt();
            int nMaps  = dis.readByte();
            String chrom = dis.readUTF();
            short chromIndex = (short) ChromIds.instance().getIndex(chrom);
            IndexArray[] maps = new IndexArray[nMaps];
            int nHapToSeqMaps = ((nMaps + 1) >> 1);
            int lastPos = 0;
            IntArray hapToSeq = null;
            for (int j=0; j<nRecs; ++j) {
                int pos = lastPos + Bref4Utils.readRestrictedInt(dis);
                Marker marker = Marker.readNonPosFields(chromIndex, pos, dis);
                lastPos = marker.pos();
                byte startMapIndex = dis.readByte();
                if (startMapIndex>=0) {
                    Bref4Utils.readMaps(maps, nHaps, startMapIndex, dis);
                    if (nHapToSeqMaps==nMaps) {
                        records.add(new IntArrayRefGTRec(marker, samples, maps[0].intArray()));
                    }
                    else {
                        if (startMapIndex < nHapToSeqMaps) {
                            hapToSeq = compose(maps, 0, nHapToSeqMaps, mapBuffer);
                        }
                        IntArray seqToAllele = compose(maps, nHapToSeqMaps, maps.length, mapBuffer);
                        records.add(new MapRefGTRec(marker, samples, hapToSeq, seqToAllele));
                    }
                }
                else if (startMapIndex==Bref4Utils.ALLELE_REC) {
                    int[][] alToHaps = Bref4Utils.readAllelesArray(dis, marker.nAlleles());
                    records.add(RefGTRec.alleleRefGTRec(marker, samples, alToHaps));
                }
                else {
                    Utilities.exit(READ_ERR);
                }
            }
        } catch (Throwable ex) {
            Utilities.exit(ex, READ_ERR);
        }
        return records.toArray(RefGTRec[]::new);
    }

    private static IntArray compose(IndexArray[] maps, int from, int to,
            int[] buffer) {
        IntArray initMap = maps[from].intArray();
        int size = initMap.size();
        for (int k=0; k<size; ++k) {
            buffer[k] = initMap.get(k);
        }
        for (int j=(from+1); j<to; ++j) {
            IntArray map = maps[j].intArray();
            for (int k=0; k<size; ++k) {
                buffer[k] = map.get(buffer[k]);
            }
        }
        return IntArray.create(buffer, 0, initMap.size(), maps[to-1].valueSize());
    }

    /**
     * Returns a string description of {@code this}.  The exact details of
     * the description are unspecified and subject to change.
     * @return a string description of {@code this}
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(80);
        sb.append(this.getClass().toString());
        return sb.toString();
    }
}