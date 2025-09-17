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
import vcf.Marker;
import vcf.RefGTRec;
import vcf.VcfRecGTParser;

/**
 * <p>Interface {@code Bref4Rec} stores a marker and the sequence
 * indices that carry each marker allele.  If the {@code Bref4Rec} instance
 * was returned by the {@code Bref4Rec.applyMap()} method, the number of
 * sequence indices may be less than the number of haplotypes in an analysis.</p>
 *
 * <p>Classes that implement {@code Bref4Rec} are required to be immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface Bref4Rec extends IntArray {

    /**
     * Returns a new {@code Bref4Rec} instance constructed from the specified
     * data.
     * @param rec a VCF record with phased, non-missing alleles.
     * @return a new {@code Bref4Rec} instance constructed from the specified
     * data
     * @throws NullPointerException if {@code (rec == null)}
     */
    public static Bref4Rec from(RefGTRec rec) {
        if (rec.marker().nAlleles()==2) {
            return new Bref4DialleleRec(rec);
        }
        else {
            return new Bref4AlleleRec(rec);
        }
    }

    /**
     * Returns a new {@code Bref4Rec} instance constructed from the specified
     * data.
     * @param gtp a VCF record with phased, non-missing alleles.
     * @return a new {@code Bref4Rec} instance constructed from the specified
     * data
     * @throws NullPointerException if {@code (gtp == null)}
     */
    public static Bref4Rec from(VcfRecGTParser gtp) {
        if (gtp.marker().nAlleles()==2) {
            return new Bref4DialleleRec(gtp);
        }
        else {
            return new Bref4AlleleRec(gtp);
        }
    }

    /**
     * Returns the marker.
     * @return the marker
     */
    Marker marker();

    /**
     * Returns the number of sequences.
     * @return the number of sequences
     */
    @Override
    int size();

    /**
     * Returns an array of length {@code this.marker().nAlleles()}
     * whose {@code j}-th element is either {@code null} or is an increasing
     * list of sequence indices that carry allele {@code j}. Exactly
     * one element of the returned array will be {@code null}.
     *
     * @return an array of length {@code this.marker().nAlleles()}
     * whose {@code j}-th element is either {@code null} or is an increasing
     * list of sequence indices that carry allele {@code j}
     */
    int[][] alleleToHaps();

    /**
     * Returns an {@code IndexArray} with {@code this.size()} elements
     * that maps sequence to allele.
     *
     * @return an {@code IndexArray} with {@code this.size()} elements
     * that maps sequence to allele
     */
    IndexArray hapToAllele();

    /**
     * Returns the index of the {@code null} element of
     * {@code this.alleleToHaps()}.
     *
     * @return the index of the {@code null} element of
     * {@code this.alleleToHaps()}
     */
    int nullRow();

    /**
     * Returns an {@code Bref4Rec} that is obtained from {@code this} by
     * applying the specified map to the haplotype indices.  The
     * contract for this method is unspecified if there exists
     * haplotype indices {@code h1} and {@code h2} such that {@code (h1 != h2)},
     * {@code (0 < h1 && h1 < this.size())},
     * {@code (0 < h2 && h2 < this.size())},
     * {@code (map.get(h1) == map.get(h2))}, and
     * {@code (this.hapToAllele().get(h1) != this.hapToAllelle.get(h2))}.
     *
     * @param map an array that maps the sequence indices onto a new set
     * of sequence indices
     * @return an {@code Bref4Rec} that is obtained from {@code this} by
     * applying the specified map to the sequence indices
     * @throws IllegalArgumentException if
     * {@code (this.size() != map.size())}
     * @throws NullPointerException if {@code (map == null)}
     */
    Bref4Rec applyMap(IndexArray map);
}
