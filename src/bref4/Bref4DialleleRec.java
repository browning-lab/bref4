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
 * <p>Class {@code Bref4DialleleRec} is a {@code Bref4Rec} instance
 * with two alleles.</p>
 *
 * <p>Instances of class {@code Bref4DialleleRec} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class Bref4DialleleRec implements Bref4Rec {

    private final Marker marker;
    private final int size;
    private final int nullAllele;
    private final int nonNullAllele;
    private final int[] nonNullSeqIndices;

    /**
     * Constructs a new {@code DiallelicAlleleCodedBref4Rec} instance from
     * the specified data.
     * @param rec a VCF record with phased, non-missing alleles
     * @throws IllegalArgumentException if {@code (rec.marker().nAlleles() != 2)}
     * @throws NullPointerException if {@code (rec == null)}
     */
    public Bref4DialleleRec(RefGTRec rec) {
        if (rec.marker().nAlleles()!=2) {
            throw new IllegalArgumentException(String.valueOf(rec.marker()));
        }
        this.marker = rec.marker();
        this.size = rec.size();
        int[][] alleleToHaps = rec.alleleToHaps();
        this.nullAllele = alleleToHaps[0]==null ? 0 : 1;
        this.nonNullAllele = (1 - nullAllele);
        this.nonNullSeqIndices = alleleToHaps[nonNullAllele];
    }

    /**
     * Constructs a new {@code DiallelicAlleleCodedBref4Rec} instance from
     * the specified data.
     * @param gtp a VCF record with phased, non-missing alleles
     * @throws IllegalArgumentException if {@code gtp.nAlleles() != 2}
     * @throws IllegalArgumentException if a format error is detected in the
     * VCF record
     * @throws IllegalArgumentException if the VCF record contains an
     * unphased or missing genotype
     * @throws NullPointerException if {@code (gtp == null)}
     */
    public Bref4DialleleRec(VcfRecGTParser gtp) {
        if (gtp.nAlleles()!=2) {
            throw new IllegalArgumentException(String.valueOf(gtp.nAlleles()));
        }
        this.marker = gtp.marker();
        this.size = 2*gtp.nSamples();
        int[][] alleleToHaps = gtp.nonMajAlleleIndices();
        this.nullAllele = alleleToHaps[0]==null ? 0 : 1;
        this.nonNullAllele = (1 - nullAllele);
        this.nonNullSeqIndices = alleleToHaps[nonNullAllele];
    }

    /**
     * Constructs a new {@code DiallelicBref4Rec} instance from the specified
     * data. The contract for this constructor is unspecified if there exists
     * haplotypes {@code h1} and {@code h2} such that
     * {@code (0 <= h1 && h1 < h2 && h2 < this.size())} and
     * {@code (rec.get(h1) != rec.get(h2))} and
     * {@code (map.get(h1) == map.get(h2))}.
     * @param map a map from haplotype to sequence
     * @param rec the {@code DiallelicBref4Rec} whose haplotype indices
     * will be mapped
     * @throws IllegalArgumentException if {@code (map.size() != rec.size())}
     * @throws NullPointerException if {@code ((map == null) || (rec == null))}
     */
    private Bref4DialleleRec(IndexArray map, Bref4DialleleRec rec) {
        if (map.size() != rec.size()) {
            throw new IllegalArgumentException(String.valueOf(map.size()));
        }
        this.marker = rec.marker();
        this.size = map.valueSize();
        this.nullAllele = rec.nullAllele;
        this.nonNullAllele = rec.nonNullAllele;
        this.nonNullSeqIndices = Arrays.stream(rec.nonNullSeqIndices)
                        .map(i -> map.get(i))
                        .distinct()
                        .sorted()
                        .toArray();
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
        if (Arrays.binarySearch(nonNullSeqIndices, hap) >= 0) {
            return nonNullAllele;
        }
        else {
            return nullAllele;
        }
    }

    @Override
    public int[][] alleleToHaps() {
        int[][] ia = new int[2][];
        ia[nonNullAllele] = nonNullSeqIndices.clone();
        return ia;
    }

    @Override
    public IndexArray hapToAllele() {
        int[] ia = new int[size];
        Arrays.fill(ia, nullAllele);
        for (int j : nonNullSeqIndices) {
            ia[j] = nonNullAllele;
        }
        int nAlleles = marker.nAlleles();
        IntArray intArray = IntArray.packedCreate(ia, nAlleles);
        return new IndexArray(intArray, nAlleles);
    }

    @Override
    public int nullRow() {
        return nullAllele;
    }

    @Override
    public Bref4Rec applyMap(IndexArray map) {
        return new Bref4DialleleRec(map, this);
    }
}
