package hw4;

import java.io.BufferedReader;
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
        CountMatrix mx = makeCountMatrix(listNames, proteins);
        FourTuple ps = new FourTuple(3);
        mx = addPseudo(mx, ps);
    }

    /**
     * Generates a count matrix given the list of protein names (and a name->sequence map)
     *
     * @param listNames names of proteins to include
     * @param proteins  mapping of name to sequence
     * @return CountMatrix of given proteins
     */
    public static CountMatrix makeCountMatrix(List<String> listNames, Map<String, String> proteins) {
        CountMatrix mx = new CountMatrix();
        for (String seqName : listNames) {
            mx.addSequence(seqName, proteins.get(seqName));
        }
        return mx;
    }

    /**
     * Returns the sum of mx and pseudo (pseudo added to each line of mx)
     *
     * @param mx     CountMatrix to add to
     * @param pseudo pseudocount vector to add
     * @return sum of mx and pseudo
     */
    public static CountMatrix addPseudo(CountMatrix mx, FourTuple pseudo) {
        return mx.addVector(pseudo);
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


    public static class CountMatrix {
        /**
         * Maps a sequence name to its base counts
         */
        public Map<String, FourTuple> counts;

        public CountMatrix() {
            this.counts = new HashMap<>();
        }

        /**
         * Adds a DNA sequence to the count matrix. If already there, replaces previous value.
         *
         * @param seqName name of the sequence to add
         * @param seq     sequence to add
         */
        public void addSequence(String seqName, String seq) {
            FourTuple tup = new FourTuple(seq);
            this.counts.put(seqName, tup);
        }

        /**
         * Returns the sum of this and vector (vector added to each line of this)
         *
         * @param vector pseudocount vector to add
         * @return sum of this and vector
         */
        public CountMatrix addVector(FourTuple vector) {
            CountMatrix result = new CountMatrix();
            for (String line : this.counts.keySet()) {
                result.counts.put(line, this.counts.get(line).add(vector));
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
         * Generate a FourTuple containing counts of all DNA bases
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
    }
}
