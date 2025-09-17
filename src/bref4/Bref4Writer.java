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
import ints.IntList;
import java.io.ByteArrayOutputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicLong;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import vcf.Marker;
import vcf.RefGTRec;
import vcf.VcfHeader;

/**
 * <p>Class {@code Bref4Writer} writes phased, non-missing genotypes
 * to a binary reference format version 4 (bref4) file.</p>
 *
 * <p>Instances of class {@code Bref4Writer} are not thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class Bref4Writer implements AutoCloseable {

    /**
     * The maximum number of samples whose data can be compressed.
     */
    public static final int MAX_SAMPLES = (1 << 30) - 1;

    private int lastChromIndex = Integer.MIN_VALUE;
    private final HashSet<Integer> chromIndices = new HashSet<>();

    private final Bref4Par par;
    private final VcfHeader vcfHeader;
    private final int maxAlleleRecNonNullCount;
    private final int[] levelToMaxNSeq;
    private final int maxMaps;
    private final int maxMapRecNAlleles;
    private final DataOutputStream bref4Out;
    private final ByteArrayOutputStream indexBytes;
    private final DataOutputStream indexOut;
    private final SeqCoder4 coder;
    private final ArrayList<RefGTRec> refGTRecs;
    private final int nHaps;
    private final ExecutorService singleThreadExecutor;
    private final AtomicLong bytesWritten; // necessary since brefOut.size() may overflow

    /**
     * Constructs a new {@code Bref4fWriter} for the specified data.The
     * Java virtual machine will exit with an error message if an I/O
     * error occurs.
     *
     * @param par the bref4 analysis parameters
     * @param vcfHeader the VCF meta-information lines and header line
     * @throws IllegalArgumentException if
     * {@code (vcfHeader.samples().size() > Bref4Writer.MAX_SAMPLES)}
     * @throws NullPointerException if
     * {@code ((par == null) || (vcfHeader == null))}
     */
    public Bref4Writer(Bref4Par par, VcfHeader vcfHeader) {
        int nSamples = vcfHeader.samples().size();
        if (nSamples > MAX_SAMPLES) {
            throw new IllegalArgumentException(String.valueOf(nSamples));
        }
        this.par = par;
        this.vcfHeader = vcfHeader;
        this.nHaps = nSamples<<1;
        this.levelToMaxNSeq = levelToMaxNSeq(nHaps, par.bitsPerLevel());
        this.maxAlleleRecNonNullCount = maxNonmajor(par, nHaps);
        this.maxMaps = levelToMaxNSeq.length + 1;
        int largestMaxNSeq = levelToMaxNSeq.length==0 ? nHaps : levelToMaxNSeq[0];
        int smallestMaxNSeq = levelToMaxNSeq.length==0 ? nHaps : levelToMaxNSeq[levelToMaxNSeq.length-1];
        this.maxMapRecNAlleles = Math.min(smallestMaxNSeq, SeqCoder4.MAX_N_ALLELES);
        this.bref4Out = Bref4Utils.dataOutputStream(par.out());
        this.indexBytes = new ByteArrayOutputStream(1<<10);
        this.indexOut = new DataOutputStream(this.indexBytes);

        this.coder = new SeqCoder4(nHaps, largestMaxNSeq);
        this.refGTRecs = new ArrayList<>();
        this.singleThreadExecutor = Executors.newSingleThreadExecutor();
        writeBref4Header(vcfHeader, bref4Out);
        this.bytesWritten = new AtomicLong(bref4Out.size());
    }

    private static int[] levelToMaxNSeq(int nHaps, int extraBitsPerLevel) {
        assert extraBitsPerLevel > 0;
        int nSeq = 16;
        int halfNHaps = nHaps>>1;
        IntList list = new IntList(10);
        while (nSeq <= halfNHaps) {
            list.add(nSeq);
            nSeq <<= extraBitsPerLevel;
        }
        return reverseArray(list);
    }

    private static int maxNonmajor(Bref4Par par, int nHaps) {
        if (par.maxNonmajor()==Bref4Par.DEF_MAX_NONMAJOR) {
            int floorLogBase2 = 31 - Integer.numberOfLeadingZeros(nHaps);
            //  5,000 samples = (1<<12) samples = (1<<13) haps  -> 13 - 11 = 2
            // 38,079 samples = (1<<15) samples = (1<<16) haps  -> 16 - 11 = 5
            return Math.max(4, 4*(floorLogBase2 - 11));
        }
        else {
            return par.maxNonmajor();
        }
    }

    private static int[] reverseArray(IntList list) {
        int[] ia = list.toArray();
        int n = ia.length >> 1;
        for (int i=0, j=ia.length-1; i<n; ++i, --j) {
            int tmp = ia[i];
            ia[i] = ia[j];
            ia[j] = tmp;
        }
        return ia;
    }

    /**
     * Write a bref4 header to the specified {@code DataOutputStream}
     * @param bref4Header the input bref4 header
     * @param bref4Out the object to which the bref4 header will be written
     * @throws NullPointerException if
     * {@code ((bref4Header == null) || (bref4Out == null))}
     */
    static void writeBref4Header(Bref4Header bref4Header,
            DataOutputStream bref4Out) {
        writeBref4Header(bref4Header.metaInfoLines(), bref4Header.samples().ids(),
                bref4Out);
    }

    /**
     * Write a bref4 header to the specified {@code DataOutputStream}
     * @param bref4Header the input VCF header
     * @param bref4Out the object to which the bref4 header will be written
     * @throws NullPointerException if
     * {@code ((vcfHeader == null) || (bref4Out == null))}
     */
    static void writeBref4Header(VcfHeader vcfHeader, DataOutputStream bref4Out) {
        writeBref4Header(vcfHeader.metaInfoLines(), vcfHeader.samples().ids(),
                bref4Out);
    }

    private static void writeBref4Header(String[] metaInfoLines, String[] ids,
            DataOutputStream bref4Out) {
        boolean quoteValue = true;
        metaInfoLines = VcfHeader.addMetaInfoLine(metaInfoLines,
                "bref4Command", Bref4Par.commandAndVersion(), quoteValue);
        Bref4Header.writeHeader(metaInfoLines, ids, bref4Out);
    }

    /**
     * Returns the bref4 analysis parameters.
     * @return the bref4 analysis parameters
     */
    public Bref4Par parameters() {
        return par;
    }

    /**
     * Returns the VCF header that is written the bref4 output file.
     * @return the VCF header that is written the bref4 output file
     */
    public VcfHeader vcfHeader() {
        return vcfHeader;
    }

    /**
     * Writes the specified phased VCF record in bref4 format.
     * The Java virtual machine will exit with an error message if an I/O
     * error occurs during method execution.
     * @param refGTRec a VCF record with phased, non-missing genotypes
     * @throws IllegalArgumentException if
     * {@code (rec.size() != 2*this.vcfHeader().samples().size())}
     * @throws NullPointerException if {@code (refGTRec == null)}
     */
    public void write(RefGTRec refGTRec) {
        if (refGTRec.size() != nHaps) {
            throw new IllegalArgumentException(String.valueOf(refGTRec.size()));
        }
        if (lastChromIndex == Integer.MIN_VALUE) {
            lastChromIndex = refGTRec.marker().chromIndex();
        }
        if (refGTRecs.size()==Integer.MAX_VALUE
                || (refGTRec.marker().chromIndex()!=lastChromIndex)) {
            writeBlockAndClearBuffer();
            lastChromIndex = refGTRec.marker().chromIndex();
            checkChromContiguity(lastChromIndex);
        }
        if (storeAsHapCodedRec(refGTRec)) {
            Bref4Rec brefRec = Bref4Rec.from(refGTRec);
            boolean success = coder.add(brefRec);
            if (success==false) {
                writeBlockAndClearBuffer();
                success = coder.add(brefRec);
                assert success;
            }
        }
        refGTRecs.add(refGTRec);
    }

    private void checkChromContiguity(int chromIndex) {
        if (chromIndices.add(chromIndex)==false) {
            String s = "ERROR: VCF records for chromosome "
                    + ChromIds.instance().id(chromIndex)
                    + " are not contiguous";
            Utilities.exit(s);
        }
    }

    private boolean storeAsHapCodedRec(RefGTRec rec) {
        return (rec.nullRowAlleleCnt() > maxAlleleRecNonNullCount)
                && (rec.marker().nAlleles() <= maxMapRecNAlleles);
    }

    private void writeBlockAndClearBuffer() {
        RefGTRec[] refGTRecsArray = refGTRecs.toArray(new RefGTRec[0]);
        ArrayList<ArrayList<IndexArray>> maps = allocateMapLists(coder.nBref4Recs());
        Bref4Rec[] mappedBrefRecs;
        if (levelToMaxNSeq.length>0 && maps.size()>0) {
            IndexArray hap2Seq = coder.hapToSeq();
            maps.get(0).add(hap2Seq);
            mappedBrefRecs = coder.mappedBref4Recs(hap2Seq);
        }
        else {
            mappedBrefRecs = coder.bref4Recs();
        }
        coder.clear();
        refGTRecs.clear();
        singleThreadExecutor.submit(() -> {
            try {
                storeMaps(0, maps, mappedBrefRecs);
                writeBref4Block(refGTRecsArray, maps);
            }
            catch (Throwable t) {
                Utilities.exit(t);
            }
        });
    }

    private static ArrayList<ArrayList<IndexArray>> allocateMapLists(int size) {
        return IntStream.range(0, size)
                .parallel()
                .mapToObj(j -> new ArrayList<IndexArray>())
                .collect(Collectors.toCollection(ArrayList<ArrayList<IndexArray>>::new));
    }

    private void storeMaps(int level, List<ArrayList<IndexArray>> maps,
            Bref4Rec[] brefRecs) {
        if (++level < levelToMaxNSeq.length && maps.size()>0) {
            SeqCoder4 subCoder = new SeqCoder4(brefRecs[0].size(), levelToMaxNSeq[level]);
            int lastStart = 0;
            for (int j=0; j<brefRecs.length; ++j) {
                boolean success = subCoder.add(brefRecs[j]);
                if (success==false) {
                    IndexArray hapToSeq = subCoder.hapToSeq();
                    maps.get(lastStart).add(hapToSeq);
                    storeMaps(level, maps.subList(lastStart, j),
                            subCoder.mappedBref4Recs(hapToSeq));
                    subCoder.clear();
                    lastStart = j;
                    success = subCoder.add(brefRecs[j]);
                    assert success;
                }
            }
            IndexArray hapToSeq = subCoder.hapToSeq();
            maps.get(lastStart).add(hapToSeq);
            storeMaps(level, maps.subList(lastStart, brefRecs.length),
                    subCoder.mappedBref4Recs(hapToSeq));
        }
        else {
            for (int j=0; j<brefRecs.length; ++j) {
                maps.get(j).add(brefRecs[j].hapToAllele());
            }
        }
    }

    private void writeBref4Block(RefGTRec[] allRecs,
            List<ArrayList<IndexArray>> maps) {
        try {
            if (allRecs.length>0) {
                Marker firstMarker = allRecs[0].marker();
                int lastPos = allRecs[allRecs.length-1].marker().pos();
                Bref4Index.writeBref4BlockToIndex(bytesWritten.get(), firstMarker.chromID(),
                        firstMarker.pos(), lastPos, indexOut);
                int mapIndex = 0;
                int lastPosOffset = 0;
                ByteArrayOutputStream buffer = new ByteArrayOutputStream(1<<16);
                DataOutputStream out = new DataOutputStream(buffer);
                out.writeInt(allRecs.length);
                out.writeInt(lastPos);
                out.writeByte(maps.isEmpty() ? 0 : maps.get(0).size()); // number of maps
                out.writeUTF(firstMarker.chromID());
                for (int j=0; j<allRecs.length; ++j) {
                    RefGTRec rec = allRecs[j];
                    lastPosOffset = writeMarker(rec.marker(), lastPosOffset, out);
                    if (storeAsHapCodedRec(rec)) {
                        writeMapRec(maps.get(mapIndex++), out);
                    }
                    else {
                        writeAlleleRec(rec.alleleToHaps(), out);
                    }
                }
                assert mapIndex==maps.size();
                bref4Out.writeInt(out.size());
                buffer.writeTo(bref4Out);
                int nBytes = Integer.BYTES + buffer.size();
                bytesWritten.addAndGet(nBytes);
            }
        } catch (Throwable ex) {
            Utilities.exit(ex, "Error writing file");
        }
    }

    private int writeMarker(Marker marker, int pos0, DataOutputStream out) throws IOException {
        int pos = marker.pos();
        Bref4Utils.writeRestrictedInt((pos - pos0), out);
        marker.writeNonPosFields(out);
        return pos;
    }

    private void writeMapRec(ArrayList<IndexArray> maps, DataOutputStream out)
            throws IOException {
        int nMaps = maps.size();
        int offset = maxMaps - nMaps;
        out.write(offset);
        for (int i=0; i<nMaps; ++i) {
            Bref4Utils.writePackedArray(maps.get(i),  out);
        }
    }

    private void writeAlleleRec(int[][] indices, DataOutputStream dos)
            throws IOException {
        dos.writeByte(Bref4Utils.ALLELE_REC);
        Bref4Utils.writeAllelesArray(indices, dos);
    }

    /**
     * Closes this {@code Bref4Writer} and releases all I/O resources.
     * @throws IOException if an I/O error occurs.
     */
    @Override
    public void close() throws IOException {
        if (refGTRecs.isEmpty()==false) {
            writeBlockAndClearBuffer();
        }
        try {
            singleThreadExecutor.shutdown();
            singleThreadExecutor.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
        }
        catch (InterruptedException e) {
            Utilities.exit(e);
        }
        try {
            bref4Out.writeInt(0);    // write sentinal
            long bref4FileOffset = bytesWritten.addAndGet(Integer.BYTES);
            Bref4Index.writeBref4Index(indexBytes, bref4FileOffset, bref4Out);
            bref4Out.close();
        } catch (IOException ex) {
            Utilities.exit(ex, "Error writing file");
        }
    }
}
