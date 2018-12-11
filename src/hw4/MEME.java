package hw4;

import java.io.IOException;
import java.util.*;

public class MEME {
    public static void main(String[] args) throws IOException {
        // Read sequence inputs (FASTA) from console
        Map<String, String> proteins = new HashMap<>(); // Maps names to seqs
        List<String> listNames = new ArrayList<>(); // List of all protein names
        Scanner console = new Scanner(System.in);
        String line = console.nextLine();
        String proteinName = line;
        String proteinContent = "";
        // Get the first protein
        while (console.hasNextLine()) {
            line = console.nextLine();
            if (line.length() > 0 && line.charAt(0) == '>') {
                // Beginning of a new protein entry
                proteins.put(proteinName, proteinContent);
                listNames.add(proteinName);
                proteinName = line;
                proteinContent = "";
            } else {
                // Continuation of a protein
                proteinContent = proteinContent + filterDNASequence(line);

            }
        }
        // Close the last protein
        proteins.put(proteinName, proteinContent);
        listNames.add(proteinName);

        // TODO
        // Create original count matrix
        FourMatrix mx = makeCountMatrix(listNames, proteins);
        FourTuple ps = new FourTuple(3);
        mx = addPseudo(mx, ps);
    }

    /**
     * Generates a count matrix given the list of protein names (and a name->sequence map)
     *
     * @param listNames names of proteins to include
     * @param proteins  mapping of name to sequence
     * @return FourMatrix of given proteins
     */
    public static FourMatrix makeCountMatrix(List<String> listNames, Map<String, String> proteins) {
        FourMatrix mx = new FourMatrix();
        for (String seqName : listNames) {
            mx.addSequence(seqName, proteins.get(seqName));
        }
        return mx;
    }

    /**
     * Returns the sum of mx and pseudo (pseudo added to each line of mx)
     *
     * @param mx     count matrix to add to
     * @param pseudo pseudocount vector to add
     * @return sum of mx and pseudo
     */
    public static FourMatrix addPseudo(FourMatrix mx, FourTuple pseudo) {
        return mx.addVector(pseudo);
    }

    /**
     * Generate a frequency matrix from a count matrix.
     *
     * @param countMatrix count matrix to use
     * @return frequency matrix (each coord = fraction of sequences that this nucleotide is here at this position)
     */
    public static FourMatrix makeFrequencyMatrix(FourMatrix countMatrix) {
        // TODO
        return null;
    }

    /**
     * Returns a String with all nonstandard bases translated to standard bases
     * (ACGT all caps, with letters other than that turned to T)
     *
     * @param sequence DNA sequence
     * @return normalized DNA sequence
     */
    private static String filterDNASequence(String sequence) {
        StringBuilder filtered = new StringBuilder();
        for (int i = 0; i < sequence.length(); i++) {
            if (sequence.charAt(i) == 'A' || sequence.charAt(i) == 'a') {
                filtered.append('A');
            } else if (sequence.charAt(i) == 'C' || sequence.charAt(i) == 'c') {
                filtered.append('C');
            } else if (sequence.charAt(i) == 'G' || sequence.charAt(i) == 'g') {
                filtered.append('G');
            } else {
                // T is the catchall letter
                filtered.append('T');
            }
        }
        return filtered.toString();
    }


    public static class FourMatrix {
        /**
         * Maps a sequence name to its base oldCounts
         */
        public Map<String, FourTuple> oldCounts;

        /**
         * List of nucleotide counts for each position where
         * counts.get(i) = counts over all included sequences for position i
         */
        public List<FourTuple> counts;

        /**
         * List of proteins included in the count
         */
        public List<String> proteins;

        public FourMatrix() {
            this.counts = new ArrayList<>();
            this.proteins = new ArrayList<>();
        }

        /**
         * Adds a DNA sequence to the count matrix. If already there, replaces previous value.
         *
         * @param seqName name of the sequence to add
         * @param seq     sequence to add
         */
        public void addSequence(String seqName, String seq) {
            this.proteins.add(seqName);

            for (int i = 0; i < seq.length(); i++) {
                if (counts.size() <= i) {
                    // Haven't already added this position
                    FourTuple count = new FourTuple();
                    count.addBase(seq.charAt(i));
                    counts.add(count);
                } else {
                    // Have already had this position before: increment
                    counts.get(i).addBase(seq.charAt(i));
                }
            }
        }

        /**
         * Returns the sum of this and vector (vector added to each line of this)
         *
         * @param vector pseudocount vector to add
         * @return sum of this and vector
         */
        public FourMatrix addVector(FourTuple vector) {
            // TODO fix
            FourMatrix result = new FourMatrix();
            for (String line : this.oldCounts.keySet()) {
                result.oldCounts.put(line, this.oldCounts.get(line).add(vector));
            }
            return result;
        }

    }


    public static class FourTuple {
        /**
         * [0]=A
         * [1]=C
         * [2]=G
         * [3]=T
         */
        public double[] a = new double[4];

        /**
         * Generates a 0 FourTuple
         */
        public FourTuple() {
        }

        /**
         * Generates a num FourTuple
         *
         * @param num default number
         */
        public FourTuple(double num) {
            this.a[0] = num;
            this.a[1] = num;
            this.a[2] = num;
            this.a[3] = num;
        }

        /**
         * Generate a FourTuple containing oldCounts of all DNA bases
         *
         * @param seq DNA sequencne to add
         */
        public FourTuple(String seq) {
            for (int i = 0; i < seq.length(); i++) {
                char letter = seq.charAt(i);
                if (letter == 'A') {
                    a[0]++;
                } else if (letter == 'C') {
                    a[1]++;
                } else if (letter == 'G') {
                    a[2]++;
                } else if (letter == 'T') {
                    a[3]++;
                } else {
                    throw new RuntimeException("seq " + seq + " has nonstandard DNA base " + letter);
                }
            }
        }

        /**
         * Returns the vector sum of this and other
         *
         * @param other other vector to add
         * @return this+other vector
         */
        public FourTuple add(FourTuple other) {
            FourTuple result = new FourTuple();
            result.a[0] = this.a[0] + other.a[0];
            result.a[1] = this.a[1] + other.a[1];
            result.a[2] = this.a[2] + other.a[2];
            result.a[3] = this.a[3] + other.a[3];
            return result;
        }

        /**
         * Increments the position associated with the given base
         * (if base is nonstandard, does nothing)
         *
         * @param base DNA base in standardized form
         */
        public void addBase(char base) {
            if (base == 'A') {
                this.a[0]++;
            } else if (base == 'C') {
                this.a[1]++;
            } else if (base == 'G') {
                this.a[2]++;
            } else if (base == 'T') {
                this.a[3]++;
            }
        }
    }
}
