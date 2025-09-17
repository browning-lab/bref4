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

import blbutil.Const;
import blbutil.Utilities;
import ints.IntArray;
import ints.WrappedIntArray;
import java.io.ByteArrayOutputStream;
import java.io.DataInput;
import java.io.DataOutputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.stream.IntStream;
import vcf.Samples;
import vcf.VcfHeader;

/**
 * <p>Class {@code Bref4Header} represents the header of a bref4 file (binary
 * reference format version 4).</p>
 *
 * <p>Instances of class {@code Bref4Header} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class Bref4Header {

    private final String source;
    private final String[] metaInfoLines;
    private final IntArray filteredHapIndices;
    private final IntArray invFilteredHapIndices;
    private final Samples filteredSamples;

    /**
     * Reads and constructs a {@code Bref4Header} instance from the specified
     * data input stream. The constructor reads an {@code int} that identifies
     * the data as bref4 data, an {@code int} that contains the number of
     * bytes used to store the meta-information lines and header lines in
     * bref4 format, a string array containing VCF meta-information lines,
     * and a string array containing sample identifiers from the data
     * input stream.  All strings are coded in modified UTF-8 format.
     * The Java Virtual Machine will print an error message and exit if an
     * I/O error or a file format error is detected.
     *
     * @param source a string describing the bref4 data source
     * @param dataIn a data input stream
     * @throws NullPointerException if
     * {@code ((source == null) || (dataIn == null))}
     */
    public Bref4Header(String source, DataInput dataIn) {
        if (source==null) {
            throw new NullPointerException("source==null");
        }
        String[] metaInfoLines0 = null;
        String[] sampleIds0 = null;
        try {
            Bref4Utils.checkMagicNumber(dataIn.readInt());
            int nBytes0 = dataIn.readInt();
            metaInfoLines0 = Bref4Utils.readStringArray(dataIn);
            sampleIds0 = Bref4Utils.readStringArray(dataIn);
        } catch (IOException ex) {
            Utilities.exit(ex, "Error reading file");
        }
        int[] filteredSampleIndices = filteredSampleIndices(sampleIds0);
        if (filteredSampleIndices.length==0) {
            noSamplesError(source);
        }
        this.source = source;
        this.metaInfoLines = metaInfoLines0;
        this.filteredHapIndices = hapIndices(filteredSampleIndices);
        this.invFilteredHapIndices = invArray(filteredHapIndices, (sampleIds0.length<<1));
        this.filteredSamples = samples(sampleIds0, filteredSampleIndices);
    }

    private static int[] filteredSampleIndices(String[] ids) {
        return IntStream.range(0, ids.length)
                .parallel()
                .toArray();
    }

    private static void noSamplesError(String source) {
        String err = "No samples remaining after sample exclusions";
        String info = Const.nl + "Error      :  " + err
                + Const.nl     + "Bref4 file :  " + source;
        Utilities.exit(new Throwable(err), info);
    }

    /**
     * Returns the array of haplotype indices corresponding to the
     * specified array of sample indices. The returned array
     * will have length {@code (2 * sampIndices.length)}.  The
     * {@code (2 * k)}-th element of the returned array is
     * {@code (2 * sampIndices[k])}. The {@code ((2 * k) + 1)}-th element
     * of the returned array is {@code ((2 * sampIndices[k]) + 1)}.
     * @param sampIndices a list of sample indices
     * @return an array of haplotype indices
     *
     * @throws IllegalArgumentException if {@code (sampIndices.length >= (1<<30))}
     * @throws NullPointerException if {@code (sampIndices == null)}
     */
    private static IntArray hapIndices(int[] sampIndices) {
        if (sampIndices.length >= (1<<30)) {
            // 2*sampIndices.length will wrap around to a negative value
            throw new IllegalArgumentException(String.valueOf(sampIndices.length));
        }
        int[] hapIndices = new int[sampIndices.length<<1];
        for (int j=0, index=0; j<sampIndices.length; ++j) {
            int hap1 = sampIndices[j]<<1;
            hapIndices[index++] = hap1;
            hapIndices[index++] = (hap1 | 0b1);
        }
        return new WrappedIntArray(hapIndices);
    }

    /**
     * Returns the list of filtered samples
     * @param sampleIds A sorted list of sample identifiers
     * @param sampleIndices a list of {@code sampleIds} array indices
     * @return the specified list of samples
     * @throws IllegalArgumentException if any two elements of
     * {@code sampleIndices} are equal
     * @throws IndexOutOfBoundsException if there is a {@code j} such that
     * {@code (0 &le j) && (j < sampleIndices.length)} and
     * {@code (sampleIndices[j] < 0) || (sampleIndices[j] >= sampleIds.length)}
     * @throws NullPointerException if
     * {@code ((sampleIds = null) || (sampleIndices == null))}
     */
    private static Samples samples(String[] sampleIds, int[] sampleIndices) {
        String[] ids = new String[sampleIndices.length];
        boolean[] isDiploid = new boolean[sampleIndices.length];
        for (int j=0; j<ids.length; ++j) {
            ids[j] = sampleIds[sampleIndices[j]];
            isDiploid[j] = true;
        }
        return new Samples(ids, isDiploid);
    }

    /**
     * Returns an array with the specified length that maps
     * {@code hapIndices[j]} to {@code j}. All indices of the returned array
     * that are not {@code hapIndices} array values are assigned the value
     * {@code -1}.
     * @param size the length of the returned array
     * @param hapIndices an array with nonnegative values
     * @return an inverse array
     *
     * @throws IllegalArgumentException if {@code size < 0}
     * @throws IllegalArgumentException if
     * {@code (hapIndices[j] == hapIndices[k])} for some
     * {@code ((0 <= j) && (j < k) && (k < hapIndices.size()))}
     * @throws IndexOutOfBoundsException if
     * {@code (hapIndices[j] < 0) || (hapIndices[j] >= size)} for some
     * {@code ((0 <= j) && (j < hapIndices.size()))}
     * @throws NullPointerException if {@code (hapIndices == null)}
     */
    private static IntArray invArray(IntArray hapIndices, int size) {
        if (size < 0) {
            throw new IllegalArgumentException(String.valueOf(size));
        }
        int[] inverseArray = IntStream.range(0, size)
                .map(j -> -1)
                .toArray();
        for (int j=0, n=hapIndices.size(); j<n; ++j) {
            if (inverseArray[hapIndices.get(j)] >= 0) {
                String s = "duplicate array value: " + hapIndices.get(j);
                throw new IllegalArgumentException(s);
            }
            inverseArray[hapIndices.get(j)] = j;
        }
        return new WrappedIntArray(inverseArray);
    }

    /**
     * Return a string description of the bref4 data source.
     * @return a string description of the bref4 data source
     */
    public String source() {
        return source;
    }

    /**
     * Returns the list of filtered haplotype indices.
     * @return the list of filtered haplotype indices
     */
    public IntArray filteredHapIndices() {
        return filteredHapIndices;
    }

    /**
     * Returns an array with length equals the to number of haplotypes
     * before sample filtering that maps {@code this.includedHapIndices()[j]}
     * to {@code j}. All other elements of the returned array will have
     * the value {@code -1}.
     * @return an array with length equals the to number of haplotypes
     * before sample filtering that maps {@code this.includedHapIndices()[j]}
     * to {@code j}.
     */
    public IntArray invFilteredHapIndices() {
        return invFilteredHapIndices;
    }

    /**
     * Returns the number of VCF meta-information lines. VCF meta-information
     * lines are lines that precede the VCF header line. A VCF meta-information
     * line must begin with "##".
     *
     * @return the number of VCF meta-information lines
     */
     public int nMetaInfoLines() {
         return metaInfoLines.length;
     }

    /**
      * Returns the specified VCF meta-information line.

      * @param index a VCF meta-information line index
      * @return the specified VCF meta-information line
      *
      * @throws IndexOutOfBoundsException if
      * {@code ((index < 0) || (index >= this.nMetaInfoLines()))}
      */
     public String metaInfoLine(int index) {
         return metaInfoLines[index];
     }

    /**
     * Returns the VCF meta-information lines.
     * @return the VCF meta-information lines
     */
    public String[] metaInfoLines() {
        return metaInfoLines.clone();
    }

     /**
      * Returns the number of samples before sample filtering.
      * @return the number of samples before sample filtering
      */
    public int nUnfilteredSamples() {
        return invFilteredHapIndices.size()>>1;
    }

     /**
      * Returns the number of samples after sample filtering.
      * @return the number of samples after sample filtering
      */
     public int nFilteredSamples() {
         return filteredSamples.size();
     }

    /**
     * Returns the list of filtered samples.
     * @return the list of filtered samples
     */
    public Samples samples() {
        return filteredSamples;
    }

    /**
     * Writes a bref4 header for the specified data to the specified
     * {@code DataOutputStream}. It is the caller's responsibility
     * to ensure that the VCF meta-information lines are correctly
     * formatted. Each string meta-information line must not terminate
     * with an the end-of-line separator.
     *
     * @param metaInfoLines the VCF meta-information lines
     * @param sampleIds the list of sample identifiers
     * @param brefOut the {@code DataOutputStream} to which the
     * bref4 header will be written
     * @throws NullPointerException if
     * {@code ((metaInfoLines == null) || (sampleIds == null) || (brefOut == null))}
     * @throws NullPointerException if any element of the {@code metaInfoLines}
     * array or {@code sampleIds} array is {@code null}
     */
    static void writeHeader(String[] metaInfoLines, String[] sampleIds,
            DataOutputStream brefOut) {
        ByteArrayOutputStream headerData = new ByteArrayOutputStream(1<<16);
        DataOutputStream out = new DataOutputStream(headerData);
        try {
            Bref4Utils.writeStringArray(metaInfoLines, out);
            Bref4Utils.writeStringArray(sampleIds, out);

            brefOut.writeInt(Bref4Utils.MAGIC_NUMBER_V4);
            brefOut.writeInt(out.size());
            headerData.writeTo(brefOut);
        } catch (IOException ex) {
            Utilities.exit(ex, "Error writing file");
        }
    }

    /**
     * Returns a {@code VcfHeader} that has the same meta-information lines
     * and the same list of samples as this {@code Bref4Header}.
     * @return a {@code VcfHeader} that has the same meta-information lines
     * and the same list of samples as this {@code Bref4Header}
     */
    public VcfHeader vcfHeader() {
        String[] sampleIds = filteredSamples.ids();
        boolean[] isDiploid = new boolean[sampleIds.length];
        Arrays.fill(isDiploid, true);
        StringBuilder sb = new StringBuilder(80 + 8*sampleIds.length);
        sb.append(VcfHeader.HEADER_PREFIX);
        for (String id : sampleIds) {
            sb.append('\t');
            sb.append(id);
        }
        int len = metaInfoLines.length + 1;
        String[] lines = Arrays.copyOf(metaInfoLines, len);
        lines[metaInfoLines.length] = sb.toString();
        return new VcfHeader(source, lines, isDiploid);
    }
}
