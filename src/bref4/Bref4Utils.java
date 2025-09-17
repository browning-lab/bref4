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
import blbutil.FileUtil;
import blbutil.Utilities;
import ints.IndexArray;
import ints.IntArray;
import java.io.ByteArrayOutputStream;
import java.io.DataInput;
import java.io.DataInputStream;
import java.io.DataOutput;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.io.File;
import java.io.UTFDataFormatException;
import vcf.Marker;

/**
 * <p>Class {@code Bref4Utils} contains static utility methods for bref4
 * format.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class Bref4Utils {

    /**
     * The initial integer in a bref4 file.
     */
    public static final int MAGIC_NUMBER_V4=25597034;

    /**
     * A byte-value indicating that the indices of haplotypes carrying
     * a VCF record's non-major alleles are stored.
     */
    public static final byte ALLELE_REC = -1;

    private static final String READ_ERR = "Error reading file";

    private Bref4Utils() {
        // private constructor to prevent instantiation
    }

    static void checkMagicNumber(long magicNumber) {
        if (magicNumber!=Bref4Utils.MAGIC_NUMBER_V4) {
            String s = "ERROR: Unrecognized input file.  Was the input file created "
                    + Const.nl + "with a different version of the bref program?";
            Utilities.exit(s);
        }
    }

    /**
     * Writes a representation of the specified integer or to the specified
     * output stream using a variable number of bytes. The specified
     * integer must satisfy {@code ((-1 <= value) && (value < 0x3f000000))}.
     * The output bytes are defined by:
     * <pre>{@code
        if (value == -1) {
            out.write(0xff);
        }
        else {
            if (value &ge; (1 << 22)) {
                value |= (0b11 << 30); // encode number of additional bytes
                out.write(value >> 24);
                out.write(value >> 16);
                out.write(value >> 8);
            } else if (value &ge; (1 << 14)) {
                value |= (0b10 << 22); // encode number of additional bytes
                out.write(value >> 16);
                out.write(value >> 8);
            } else if (value &ge; (1 << 6)) {
                value |= (0b01 << 14); // encode number of additional bytes
                out.write(value >> 8);
            }
            out.write(value);
        }
     * }</pre>
     * @param value the value to be written
     * @param out the output destination
     * @throws IllegalArgumentException if
     * {@code ((value < -1) || (value >= 0x3f000000))}
     * @throws IOException if an I/O error occurs
     */
    public static void writeRestrictedInt(int value, DataOutput out)
            throws IOException {
        if (value == -1) {
            out.write(0xff);
        }
        else {
            if (value < 0 || value >= 0x3f000000) { // ensure that first stored byte is not 0xff
                throw new IllegalArgumentException(String.valueOf(value));
            }
            if (value >= (1 << 22)) {
                value |= (0b11 << 30); // encode number of additional bytes
                out.write(value >> 24);
                out.write(value >> 16);
                out.write(value >> 8);
            } else if (value >= (1 << 14)) {
                value |= (0b10 << 22); // encode number of additional bytes
                out.write(value >> 16);
                out.write(value >> 8);
            } else if (value >= (1 << 6)) {
                value |= (0b01 << 14); // encode number of additional bytes
                out.write(value >> 8);
            }
            out.write(value);
        }
    }

    /**
     * Reads and returns an integer {@code x} in the range
     * {@code ((-1 <= x) && (x < 0x3f000000))} that was written
     * with the {@code MarkerUtils.writeRestrictedInt()} method.
     * @param in the input source
     * @return an integer {@code x} in the range
     * {@code ((-1 <= x) && (x < 0x3f000000))}
     * @throws IOException if an I/O error occurs
     * @throws NullPointerException if {@code (in == null)}
     */
    public static int readRestrictedInt(DataInput in) throws IOException {
        int value = in.readUnsignedByte();
        if (value==0xff) {
            return -1;
        }
        else {
            int nAdditionalBytes = value >> 6;
            value &= 0x3f; // clear first two bits storing nAdditionalBytes
            for (int j=0; j<nAdditionalBytes; ++j) {
                value <<= 8;
                value |= in.readUnsignedByte();
            }
            return value;
        }
    }


    /**
     * Writes the specified string array to the specified output stream.  The
     * output stream is written with the following code:
    * <pre>
        out.writeInt(sa.length);
        for (String s : sa) {
            out.writeUTF(s);
        }
    * </pre>
     * @param sa a string array
     * @param out the data output stream to which the marker will be written
     * @throws IOException if an I/O error occurs
     * @throws NullPointerException if {@code (sa == null)} or if there is
     * an array fieldsIndex {@code j} such that
     * {@code ((0 <= j) && (j < sa.length) && (sa[j] == null))}
     * @throws NullPointerException if {@code (out == null)}
     */
    public static void writeStringArray(String[] sa, DataOutput out)
            throws IOException {
        out.writeInt(sa.length);
        for (String s : sa) {
            out.writeUTF(s);
        }
    }

    /**
     * Reads and returns a string array from the specified input stream.
     * A string array must be encoded in the input stream as a
     * 4 byte integer specifying the array length followed by the string
     * array elements in order of increasing index.  Each string array
     * element must be encoded in modified UTF-8 format.
     * @param dataIn a data input stream
     * @return a string array or {@code null}
     * @throws EOFException if specified stream ends before reading all the
     * bytes
     * @throws IOException if an I/O error occurs
     * @throws IllegalArgumentException if the integer specifying the array
     * length is negative
     * @throws NullPointerException if {@code (dataIn ==  null)}
     * @throws UTFDataFormatException  if the bytes in a string read from
     * {@code dataIn} are not a valid modified UTF-8 encoding of a string.
     */
    public static String[] readStringArray(DataInput dataIn) throws IOException {
        int length = dataIn.readInt();
        if (length<0) {
            throw new IllegalArgumentException(String.valueOf(length));
        } else if (length==0) {
            return new String[0];
        } else {
            String[] sa = new String[length];
            for (int j=0; j<sa.length; ++j) {
                sa[j] = dataIn.readUTF();
            }
            return sa;
        }
    }

    /**
     * Returns a buffered {@code DataInputStream} that reads from the
     * specified file. If the input stream cannot be opened, an error
     * message will be printed and the Java interpreter will exit.
     * @param file an input file
     * @return a buffered {@code DataInputStream}
     * @throws NullPointerException if {@code (file == null)}
     */
    public static DataInputStream dataInputStream(File file) {
        return new DataInputStream(FileUtil.bufferedInputStream(file));
    }

    /**
     * Returns a buffered {@code DataOutputStream} that writes to the
     * specified file.
     * @param pathname the file pathname
     * @return a buffered {@code DataOutputStream} that writes to the
     * specified file
     * @throws NullPointerException if {@code (pathname == null)}
     */
    public static DataOutputStream dataOutputStream(String pathname) {
        return new DataOutputStream(FileUtil.bufferedOutputStream(new File(pathname)));
    }

    /**
     * Writes the specified index array to the specified output stream.
     * The index array is encoded as a byte specifying the number of bits
     * encoding each array element followed by a sequence of bytes containing
     * the bits that encode each array element. If an array element is
     * encoded using {@code nBitsPerValue} bits, then the array elements
     * will be stored in {@code (((size*nBitsPerValue) + 7) / 8)} bytes.
     * The number of bits used to encode each elements is
     * {@code (Integer.SIZE - Integer.numberOfLeadingZeros(ia.valueSize() - 1))}.
     * @param ia an index array
     * @param out a data output stream
     * @throws IllegalArgumentException if there exists {@code j} such that
     * {@code ((0 <= j < ia.size()) && (ia.get(j) >= ia.valueSize()))}
     * @throws IOException if an I/O error occurs
     * @throws NullPointerException if {@code ((ia == null) || (out ==  null))}
     */
    public static void writePackedArray(IndexArray ia, DataOutput out)
            throws IOException {
        int bitsPerValue =
                (Integer.SIZE - Integer.numberOfLeadingZeros(ia.valueSize() - 1));
        int size = ia.size();
        long nBits = size*bitsPerValue;
        int valueSize = ia.valueSize();
        int divideBy64Shift = 6;
        long nWords = (size*bitsPerValue + (Long.SIZE - 1)) >> divideBy64Shift;
        long[] words = new long[(int) nWords];
        long bitIndex = 0;
        for (int j=0; j<size; ++j) {
            int value = ia.get(j);
            if (value<0 || value>=valueSize) {
                throw new IllegalArgumentException(String.valueOf(value));
            }
            for (int k=0; k<bitsPerValue; ++k) {
                if ((value & 1)==1) {
                    words[(int) (bitIndex >> divideBy64Shift)] |= (1L << bitIndex);
                }
                value >>= 1;
                ++bitIndex;
            }
        }
        writeRestrictedInt(valueSize, out);
        writePackedData(words, nBits, out);
    }

    private static void writePackedData(long[] words, long nBits, DataOutput out)
            throws IOException{
        int nLeftOverBits = (int) (nBits & 63);
        if (nLeftOverBits==0 || nLeftOverBits>56) {
            for (int j=0; j<words.length; ++j) {
                out.writeLong(words[j]);
            }
        }
        else {
            int lengthM1 = words.length-1;
            for (int j=0; j<lengthM1; ++j) {
                out.writeLong(words[j]);
            }
            long remainingBytes = words[lengthM1];
            for (int j=0; j<nLeftOverBits; j+=8) {
                out.write((byte) remainingBytes);
                remainingBytes >>= 8;
            }
        }
    }

    /**
     * Reads and returns an index array from the specified input stream.
     * The index array must be encoded as a byte specifying the number of bits
     * encoding each array element followed by a sequence of bytes
     * containing the bits that encode each array element.  If an array element
     * is encoded using {@code nBitsPerValue} bits, then the elements in the
     * array should be stored in {@code (((size*nBitsPerValue) + 7) / 8)} bytes.
     * @param dataIn a data input stream
     * @param size the length of the array that will be read
     * @return an fieldsIndex array
     * @throws EOFException if specified stream ends before reading all the
     * bytes
     * @throws IllegalArgumentException if {@code (size <= 0)}
     * @throws IllegalArgumentException if the first byte read is negative or
     * greater than or equal to {@code Integer.SIZE}
     * @throws IOException if an I/O error occurs
     * @throws NullPointerException if {@code (dataIn ==  null)}
     */
    public static IndexArray readPackedIndexArray(DataInput dataIn, int size)
            throws IOException {
        if (size<=0) {
            throw new IllegalArgumentException(String.valueOf(size));
        }
        int valueSize = readRestrictedInt(dataIn);
        byte nBitsPerValue =
                (byte) (Integer.SIZE - Integer.numberOfLeadingZeros(valueSize - 1));
        long[] packedArray = readPackedArray(size, nBitsPerValue, dataIn);
        return toIndexArray(packedArray, nBitsPerValue, size, valueSize);
    }

    private static long[] readPackedArray(int size, byte nBitsPerValue,
            DataInput dataIn) throws IOException {
        final int divideBy64Shift = 6;
        long nBits = ((long) size)*nBitsPerValue;
        int nLeftOverBits = (int) (nBits & 63);
        long nWords = (nBits + (Long.SIZE - 1)) >> divideBy64Shift;
        long[] words = new long[(int) nWords];
        if (nLeftOverBits==0 || nLeftOverBits>56) {
            for (int j=0; j<words.length; ++j) {
                words[j] = dataIn.readLong();
            }
        }
        else {
            int lengthM1 = words.length-1;
            for (int j=0; j<lengthM1; ++j) {
                words[j] = dataIn.readLong();
            }
            for (int offset=0; offset<nLeftOverBits; offset+=8) {
                words[lengthM1] |= ((dataIn.readByte() & 0xffL) << offset);
            }
        }
        return words;
    }

    private static IndexArray toIndexArray(long[] packedArray, int nBitsPerValue,
            int size, int valueSize) {
        assert (0<nBitsPerValue) && (nBitsPerValue<(Integer.SIZE-1));
        final int divideBy64Shift = 6;
        int[] values = new int[size];
        long valueMask = (1<<nBitsPerValue) - 1;
        long startBit = 0;
        for (int j=0; j<values.length; ++j) {
            int wordIndex = (int) (startBit >> divideBy64Shift);
            long bitOffset = startBit & 63;  // equals (startBit % Long.SIZE)
            long bits = packedArray[wordIndex] >>> bitOffset;
            if (bitOffset + nBitsPerValue > Long.SIZE) {
                bits |= (packedArray[wordIndex+1] << (-bitOffset));
            }
            int value = (int) (bits & valueMask);
            if (value<0 || value>=valueSize) {
                throw new IllegalArgumentException(String.valueOf(value));
            }
            values[j] = value;
            startBit += nBitsPerValue;
        }
        IntArray intArray = IntArray.create(values, valueSize);
        return new IndexArray(intArray, valueSize);
    }

    /**
     * Writes the specified marker to the specified output stream.
     * @param marker the marker to be written
     * @param lastPos the position of the most recent marker written, or 0
     * @param out the data output stream
     * @throws NullPointerException if {@code (marker == null) || (out == null)}
     * @throws IOException if an {@code IOException} is encountered
     */
    public static void writeMarker(Marker marker, int lastPos, DataOutputStream out)
            throws IOException {
        Bref4Utils.writeRestrictedInt((marker.pos() - lastPos), out);
        marker.writeNonPosFields(out);
    }

    /**
     * Writes the specified allele record data.
     * @param indices the array to be written
     * @param out the data output stream
     * @throws IOException if an {@code IOException} is encountered
     * @throws NullPointerException if
     * {@code (indices == null) || (out == null)}
     */
    public static void writeAllelesArray(int[][] indices, DataOutputStream out)
            throws IOException {
        for (int[] ia : indices) {
            if (ia==null) {
                writeRestrictedInt(-1, out);
            } else {
                writeRestrictedInt(ia.length, out);
                for (int hap : ia) {
                    out.writeInt(hap);
                }
            }
        }
    }

    /**
     * Reads the allele record data written with the {@code writeAllelesArray()}
     * method from the specified input stream and returns the allele
     * record data.
     * @param dis the input stream from which the data will be read
     * @param nAlleles the number of {@code int[]} elements in the returned
     * array
     * @return an integer array
     * @throws EOFException if the input stream reaches the end before
     * reading all of the bytes
     * @throws IOException if an {@code IOException} occurs
     * @throws NullPointerException if {@code (dis == null)}
     */
    public static int[][] readAllelesArray(DataInputStream dis, int nAlleles)
            throws IOException {
        int[][] alToHaps = new int[nAlleles][];
        for (int j = 0; j < nAlleles; ++j) {
            int length = Bref4Utils.readRestrictedInt(dis);
            if (length == -1) {
                alToHaps[j] = null;
            } else {
                int[] ia = new int[length];
                for (int k=0; k<ia.length; ++k) {
                    ia[k] = dis.readInt();
                }
                alToHaps[j] = ia;
            }
        }
        return alToHaps;
    }

    /**
     * Writes the specified map record to the specified output stream
     * @param maps the haplotype maps
     * @param startIndex the index of the first haplotype map to be written
     * @param dos the output stream
     * @throws IOException if an {@code IOException} is encountered
     * @throws NullPointerException if
     * {@code ((maps == null) || (dos == null))}
     * @throws NullPointerException if
     * {@code ((startIndex <= j) && (j < maps.length) && (maps[j] == null))}
     */
    public static void writeMaps(IndexArray[] maps, int startIndex,
            DataOutputStream dos) throws IOException {
        dos.write(startIndex);
        for (int j=startIndex; j<maps.length; ++j) {
            Bref4Utils.writePackedArray(maps[j],  dos);
        }
    }

    /**
     * Reads the specified map records from the specified input stream.
     * @param maps the haplotype maps
     * @param startIndex the index of the first haplotype map to be read
     * @param nHaps the number of haplotypes
     * @param dis an input stream
     * @throws IOException if an {@code IOException} is encountered
     * @throws NullPointerException if
     * {@code ((maps == null) || (dis == null))}
     * @throws NullPointerException if
     * {@code ((startIndex <= j) && (j < maps.length) && (maps[j] == null))}
     */
    public static void readMaps(IndexArray[] maps, int nHaps,
            int startIndex, DataInputStream dis) throws IOException {
        int size = startIndex==0 ? nHaps : maps[startIndex-1].valueSize();
        for (int k=startIndex; k<maps.length; ++k) {
            maps[k] = Bref4Utils.readPackedIndexArray(dis, size);
            size = maps[k].valueSize();
        }
    }

    /**
     * Returns a byte array containing the specified 4-byte {@code nRecs}
     * integer, followed by the specified 4-byte {@code lastOutPos} integer,
     * followed by the bytes in the specified {@code ByteArrayOutputStream}.
     * The {@code nRecs} and {@code lastOutPos} integers are written in
     * big-endian format.
     * @param nRecs the number of records
     * @param lastOutPos the last VCF record POS
     * @param buffer a sequence of bytes
     * @return a byte array containing the specified data
     * @throws IllegalArgumentException if {@code (nRecs < 1)}
     * @throws NullPointerException if {@code (buffer == null)}
     */
    public static byte[] prependData(int nRecs, int lastOutPos,
            ByteArrayOutputStream buffer) {
        if (nRecs < 1) {
            throw new IllegalArgumentException(String.valueOf(nRecs));
        }
        try {
            int nBytes = (2*Integer.BYTES) + buffer.size();
            ByteArrayOutputStream baos = new ByteArrayOutputStream(nBytes);
            baos.write(nRecs >>> 24);
            baos.write(nRecs >>> 16);
            baos.write(nRecs >>> 8);
            baos.write(nRecs);
            baos.write(lastOutPos >>> 24);
            baos.write(lastOutPos >>> 16);
            baos.write(lastOutPos >>> 8);
            baos.write(lastOutPos);
            buffer.writeTo(baos);
            return baos.toByteArray();
        } catch (IOException ex) {
            assert false;
            Utilities.exit(ex, READ_ERR);
        }
        assert false;
        return null;
    }

    /**
     * Returns a byte array containing the specified 4-byte {@code nRecs}
     * integer, followed by the specified 4-byte {@code lastOutPos} integer,
     * followed by the least significant byte of the specified {@code nMaps}
     * integer, followed by the bytes in the specified
     * {@code ByteArrayOutputStream}. The {@code nRecs} and {@code lastOutPos}
     * integers are written in big-endian format.
     * @param nRecs the number of records
     * @param lastOutPos the last VCF record POS
     * @param nMaps the number of maps in {@code MapRefGTRec} objects
     * @param buffer a sequence of bytes
     * @return a byte array containing the specified data
     * @throws IllegalArgumentException if {@code (nRecs < 1)}
     * @throws IllegalArgumentException if {@code (nMaps > 255)}
     * @throws NullPointerException if {@code (buffer == null)}
     */
    public static byte[] prependData(int nRecs, int lastOutPos, int nMaps,
            ByteArrayOutputStream buffer) {
        if (nRecs < 1) {
            throw new IllegalArgumentException(String.valueOf(nRecs));
        }
        if (nMaps > 255) {
            throw new IllegalArgumentException(String.valueOf(nMaps));
        }
        try {
            int nBytes = (2*Integer.BYTES + Byte.BYTES) + buffer.size();
            ByteArrayOutputStream baos = new ByteArrayOutputStream(nBytes);
            baos.write(nRecs >>> 24);
            baos.write(nRecs >>> 16);
            baos.write(nRecs >>> 8);
            baos.write(nRecs);
            baos.write(lastOutPos >>> 24);
            baos.write(lastOutPos >>> 16);
            baos.write(lastOutPos >>> 8);
            baos.write(lastOutPos);
            baos.write(nMaps);
            buffer.writeTo(baos);
            return baos.toByteArray();
        } catch (IOException ex) {
            assert false;
            Utilities.exit(ex, READ_ERR);
        }
        assert false;
        return null;
    }
}
