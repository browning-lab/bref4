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
import java.util.Arrays;
import vcf.Marker;
import vcf.RefGTRec;
import vcf.VcfRecGTParser;

/**
 * <p>Class {@code Bref4AlleleRec} is a {@code Bref4Rec} instance.</p>
 *
 * <p>Instances of class {@code Bref4AlleleRec} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class Bref4AlleleRec implements Bref4Rec {

    private final Marker marker;
    private final int size;
    private final int[][] alleleToSeqs;
    private final int nullRow;

    /**
     * Constructs a new {@code AllelicBref4Rec} instance from the specified
     * data.
     * @param rec a VCF record with phased, non-missing alleles.
     * @throws NullPointerException if {@code (rec == null)}
     */
    public Bref4AlleleRec(RefGTRec rec) {
        this.marker = rec.marker();
        this.size = rec.size();
        this.alleleToSeqs = rec.alleleToHaps();
        this.nullRow = RefGTRec.nullRow(alleleToSeqs);
    }

    /**
     * Constructs a new {@code AllelicBref4Rec} instance from the specified
     * data.
     * @param gtp a VCF record with phased, non-missing alleles
     * @throws IllegalArgumentException if a format error is detected in the
     * VCF record
     * @throws IllegalArgumentException if the VCF record contains an
     * unphased or missing genotype
     * @throws NullPointerException if {@code (gtp == null)}
     */
    public Bref4AlleleRec(VcfRecGTParser gtp) {
        this.marker = gtp.marker();
        this.size = gtp.nSamples()<<1;
        this.alleleToSeqs = gtp.nonMajAlleleIndices();
        this.nullRow = RefGTRec.nullRow(alleleToSeqs);
    }

    /**
     * Constructs a new {@code AllelicBref4Rec} instance from the specified
     * data. The contract for this constructor is unspecified if there exists
     * haplotypes {@code h1} and {@code h2} such that
     * {@code (0 <= h1 && h1 < h2 && h2 < this.size())} and
     * {@code (rec.get(h1) != rec.get(h2))} and
     * {@code (map.get(h1) == map.get(h2))}
     * @param map a map from haplotype to sequence
     * @param rec the {@code AllelicBref4Rec} whose haplotype indices will be
     * mapped
     * @throws IllegalArgumentException if {@code (map.size() != rec.size())}
     * @throws NullPointerException if
     * {@code ((map == null) || (rec == null))}
     */
    private Bref4AlleleRec(IndexArray map, Bref4AlleleRec rec) {
        if (map.size() != rec.size()) {
            throw new IllegalArgumentException(String.valueOf(map.size()));
        }
        int[][] newAlleleToHaps = new int[rec.alleleToSeqs.length][];
        for (int j=0; j<rec.alleleToSeqs.length; ++j) {
            if (rec.alleleToSeqs[j]!=null) {
                newAlleleToHaps[j] = Arrays.stream(rec.alleleToSeqs[j])
                        .map(i -> map.get(i))
                        .distinct()
                        .sorted()
                        .toArray();
            }
        }
        this.marker = rec.marker();
        this.size = map.valueSize();
        this.nullRow = rec.nullRow;
        this.alleleToSeqs = newAlleleToHaps;
    }

    @Override
    public Marker marker() {
        return marker;
    }

    @Override
    public int size() {
        return size;
    }

    @Override
    public int get(int hap) {
        if (hap < 0 || hap >= size) {
            throw new IndexOutOfBoundsException(String.valueOf(hap));
        }
        for (int j=0; j<alleleToSeqs.length; ++j) {
            if (j != nullRow) {
                if (Arrays.binarySearch(alleleToSeqs[j], hap) >= 0) {
                    return j;
                }
            }
        }
        return nullRow;
    }

    @Override
    public int[][] alleleToHaps() {
        int[][] ia = new int[alleleToSeqs.length][];
        for (int j=0; j<ia.length; ++j) {
            if (j!=nullRow) {
                ia[j] = alleleToSeqs[j].clone();
            }
        }
        return ia;
    }

    @Override
    public IndexArray hapToAllele() {
        int[] ia = new int[size];
        Arrays.fill(ia, nullRow);
        for (int j=0; j<alleleToSeqs.length; ++j) {
            if (j!=nullRow) {
                for (int k : alleleToSeqs[j]) {
                    ia[k] = j;
                }
            }
        }
        int nAlleles = marker.nAlleles();
        IntArray intArray = IntArray.packedCreate(ia, nAlleles);
        return new IndexArray(intArray, nAlleles);
    }

    @Override
    public int nullRow() {
        return nullRow;
    }

    @Override
    public Bref4Rec applyMap(IndexArray map) {
        return new Bref4AlleleRec(map, this);
    }
}
