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

import ints.IndexArray;
import ints.IntArray;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.stream.IntStream;

/**
 * <p>Class {@code SeqCoder4} stores a map from haplotype to allele
 * sequence. The map is created from a list of {@code Bref4Rec} objects.</p>
 *
 * <p>Instances of class {@code SeqCoder4} are not thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class SeqCoder4 {

    /**
     * The maximum number permitted alleles in a {@code Bref4Rec} object
     * passed to the {@code SeqCoder4.add()} method.
     */
    public static final int MAX_N_ALLELES = 256;

    private static final int NOT_ASSIGNED = -1;
    private static final int ASSIGNED = Integer.MAX_VALUE;

    private final int nHaps;
    private final int maxNSeq;
    private final int[] hap2Seq;
    private final int[] seq2Cnt;

    private int nSeq;
    private int updateMapSize;
    private final int[][] updateMap;
    private final ArrayList<Bref4Rec> bref4RecList;

    /**
     * Constructs a new {@code SeqCoder4} instance.
     * @param nHaps the number of haplotypes
     * @param maxNSeq the maximum permitted number of distinct allele
     * sequences.
     *
     * @throws IllegalArgumentException if {@code nHaps < 1 || maxNSeq < 1}
     */
    public SeqCoder4(int nHaps, int maxNSeq) {
        if (nHaps < 1) {
            throw new IllegalArgumentException(String.valueOf(nHaps));
        }
        if (maxNSeq < 1) {
            throw new IllegalArgumentException(String.valueOf(maxNSeq));
        }
        this.nHaps = nHaps;
        this.maxNSeq = maxNSeq;

        this.hap2Seq = new int[nHaps];
        this.seq2Cnt = new int[maxNSeq];
        this.seq2Cnt[0] = nHaps;
        this.nSeq = 1;
        this.updateMapSize = 8;
        this.updateMap = initializeUpdateMap(updateMapSize, maxNSeq);
        this.bref4RecList = new ArrayList<>();
    }

    private static int[][] initializeUpdateMap(int nMapAlleles, int maxNSeq) {
        int[][] updateMap = new int[MAX_N_ALLELES][];
        for (int j=0; j<nMapAlleles; ++j) {
            updateMap[j] = new int[maxNSeq];
        }
        return updateMap;
    }

    /**
     * Attempts to add the specified record to the map from haplotype
     * to allele sequence, and returns {@code true} if the record
     * was successfully added. If {@code false} is returned,
     * {@code this.clear()} should be invoked before repeating the
     * attempt to add the record.
     *
     * @param rec a list of haplotype alleles
     * @return {@code true} if the specified record was successfully added
     *
     * @throws IllegalArgumentException if {@code rec.size() != this.nHaps()}
     * @throws IllegalArgumentException if
     * {@code rec.marker().nAlleles() > SeqCoder4.MAX_N_ALLELES}
     * @throws NullPointerException if {@code rec == null}
     */
    public boolean add(Bref4Rec rec) {
        if (rec.size()!=nHaps) {
            throw new IllegalArgumentException(String.valueOf(rec.size()));
        }
        int nAlleles = rec.marker().nAlleles();
        if (nAlleles>MAX_N_ALLELES) {
            throw new IllegalArgumentException(String.valueOf(nAlleles));
        }
        if (nAlleles>updateMapSize) {
            addToUpdateMap(nAlleles);
        }
        int[][] alleleToHaps = rec.alleleToHaps();
        int nullAllele = rec.nullRow();
        boolean success = setUpdateMap(alleleToHaps, nullAllele);
        if (success) {
            updateHap2Seq(alleleToHaps);
            bref4RecList.add(rec);
        }
        return success;
    }

    private void addToUpdateMap(int nAlleles) {
        for (int j=updateMapSize; j<nAlleles; ++j) {
            updateMap[j] = new int[maxNSeq];
        }
        updateMapSize = nAlleles;
    }

    private boolean setUpdateMap(int[][] alleleToHap, int nullAllele) {
        int nSeqAtStart = nSeq;
        int nAlleles = alleleToHap.length;
        resetUpdateMap(alleleToHap, nullAllele);
        int[] nullAlleleMap = updateMap[nullAllele];
        for (int al=0; al<nAlleles; ++al) {
            if (al!=nullAllele) {
                int[] alleleMap = updateMap[al];
                for (int h : alleleToHap[al]) {
                    int seq = hap2Seq[h];
                    if (alleleMap[seq]==NOT_ASSIGNED) {
                        if (nullAlleleMap[seq]==NOT_ASSIGNED) {
                            nullAlleleMap[seq] = ASSIGNED;
                            alleleMap[seq] = seq;
                        }
                        else {
                            alleleMap[seq] = nSeq++;
                        }
                    }
                }
            }
        }
        boolean success = nSeq <= maxNSeq;
        if (success==false) {
            nSeq = nSeqAtStart;
        }
        else {
            Arrays.fill(seq2Cnt, nSeqAtStart, nSeq, 0);
        }
        return success;
    }

    private void resetUpdateMap(int[][] alleleToHaps, int nullAllele) {
        int[] seqToNullAlleleCnt = seqToNullAlleleCnt(alleleToHaps);
        IntStream.range(0, alleleToHaps.length)
                .parallel()
                .forEach(i -> Arrays.fill(updateMap[i], 0, nSeq, NOT_ASSIGNED));

        for (int seq=0; seq<nSeq; ++seq) {
            if (seqToNullAlleleCnt[seq]>0) {
                updateMap[nullAllele][seq] = seq;
            }
        }
    }

    private int[] seqToNullAlleleCnt(int[][] alleleToHaps) {
        int[] nullAlleleCnt = Arrays.copyOf(seq2Cnt, nSeq);
        IntStream.range(0, alleleToHaps.length)
                .flatMap(i -> alleleToHaps[i]==null ? IntStream.empty() : Arrays.stream(alleleToHaps[i]))
                .forEach(i -> --nullAlleleCnt[hap2Seq[i]]);
        return nullAlleleCnt;
    }

    private void updateHap2Seq(int[][] alleleToHaps) {
        for (int al=0; al<alleleToHaps.length; ++al) {
            int[] alMap = updateMap[al];
            if (alleleToHaps[al]!=null) {
                for (int h : alleleToHaps[al]) {
                    int oldSeq = hap2Seq[h];
                    int newSeq = alMap[oldSeq];
                    if (newSeq != oldSeq) {
                        hap2Seq[h] = newSeq;
                        --seq2Cnt[oldSeq];
                        ++seq2Cnt[newSeq];
                    }
                }
            }
        }
    }

    /**
     * Returns the number of allele sequences.
     * @return the number of allele sequences
     */
    public int nSeq() {
        return nSeq;
    }

    /**
     * Return the map from haplotype to allele sequence.
     * @return the map from haplotype to allele sequence
     */
    public IndexArray hapToSeq() {
        IntArray intArray = IntArray.packedCreate(hap2Seq, nSeq);
        return new IndexArray(intArray, nSeq);
    }

    /**
     * Returns {@code this.bref4Recs().length}.
     * @return {@code this.bref4Recs().length}
     */
    public int nBref4Recs() {
        return bref4RecList.size();
    }

    /**
     * Returns the list of {@code Bref4Rec} objects that have been
     * added to the map from haplotype to allele sequence.
     * @return the list of {@code Bref4Rec} objects that have been
     * added to the map from haplotype to allele sequence
     */
    public Bref4Rec[] bref4Recs() {
        return bref4RecList.toArray(Bref4Rec[]::new);
    }

    /**
     * Returns the list of {@code Bref4Rec} objects obtained by
     * applying the specified map to the elements of {@code this.brefRecs()}.
     * @param hapToSeq a map from haplotype to sequence
     * @return the list of {@code Bref4Rec} objects obtained by
     * applying the specified map to the elements of {@code this.brefRecs()}
     * @throws IllegalArgumentException if {@code hapToSeq.size() != this.nHaps()}
     * @throws NullPointerException if {@code hapToSeq == null}
     */
    public Bref4Rec[] mappedBref4Recs(IndexArray hapToSeq) {
        if (hapToSeq.size() != nHaps) {
            throw new IllegalArgumentException(String.valueOf(hapToSeq.size()));
        }
        return bref4RecList.stream()
                .parallel()
                .map(rec -> rec.applyMap(hapToSeq))
                .toArray(Bref4Rec[]::new);
    }

    /**
     * Clears the stored list of {@code Bref4Rec} objects and the
     * associated map from haplotype to allele sequence. When this method
     * returns, all haplotypes will map to the empty allele sequence.
     */
    public void clear() {
        // initialize with empty sequence (seq index is 0)
        Arrays.fill(hap2Seq, 0);
        nSeq = 1;
        seq2Cnt[0] = nHaps;
        bref4RecList.clear();
    }

    /**
     * Returns the number of haplotypes.
     * @return the number of haplotypes
     */
    public int nHaps() {
        return nHaps;
    }

    /**
     * Returns the maximum permitted number of distinct allele sequences.
     * @return the maximum permitted number of distinct allele sequences
     */
    public int maxNSeq() {
        return maxNSeq;
    }
}