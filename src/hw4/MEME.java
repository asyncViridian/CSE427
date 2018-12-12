package hw4;

import java.io.IOException;
import java.util.*;

public class MEME {
    public static void main(String[] args) throws IOException {
        // Read sequence inputs (FASTA, DNA) from console
        List<String> DNASequences = new ArrayList<>(); // Maps names to seqs
        Scanner console = new Scanner(System.in);
        String line = console.nextLine();
        String DNAcontent = "";
        // Get the first protein
        while (console.hasNextLine()) {
            line = console.nextLine();
            if (line.length() > 0 && line.charAt(0) == '>') {
                // Beginning of a new protein entry
                DNASequences.add(DNAcontent);
                DNAcontent = "";
            } else {
                // Continuation of a protein
                DNAcontent = DNAcontent + filterDNASequence(line);

            }
        }
        // Close the last protein
        DNASequences.add(DNAcontent);

        // Generate background frequency
        FourTuple background = new FourTuple(0.25);
        // Generate pseudocount frequency
        FourTuple pseudocounts = new FourTuple(0.25);
        // Get seed pseudocount frequency (this is a special case!
        FourTuple seedPseudocount = new FourTuple(0.15 / ((4 * 0.85) - 1));
        // Initialize motif width
        int motifWidth = 10;

        // Get seed subsequences (EM Initialization!)
        List<String> seeds = new ArrayList<>();
        String firstSeq = DNASequences.get(0);
        for (int i = 0; i + motifWidth < firstSeq.length(); i += motifWidth / 2) {
            seeds.add(firstSeq.substring(i, i + motifWidth));
        }
        seeds.add(firstSeq.substring(firstSeq.length() - motifWidth));


        // Collect entropy data
        List<List<Double>> entropies = new ArrayList<>();

        // Initialize WMMs for each of the S seed WMMs
        List<FourMatrix> WMMSeeds = new ArrayList<>();
        for (int i = 0; i < seeds.size(); i++) {
            entropies.add(new ArrayList<Double>());

            String seed = seeds.get(i);

            // Initialize some stuff
            List<String> subseqs = new ArrayList<>();
            subseqs.add(seed);
            List<Double> weights = new ArrayList<>();
            weights.add(1.0);
            // Create the WMMs
            FourMatrix countMatrix = makeCountMatrix(subseqs, weights);
            FourMatrix totalCounts = addPseudo(countMatrix, seedPseudocount);
            FourMatrix frequencyMatrix = makeFrequencyMatrix(totalCounts);
            FourMatrix resultWMM = makeWMM(frequencyMatrix, background);

            WMMSeeds.add(resultWMM);
            entropies.get(i).add(resultWMM.entropy);
        }

        // For each of the S seed WMMs, run three E/M step pairs
        List<FourMatrix> WMMResults = new ArrayList<>();
        for (int i = 0; i < WMMSeeds.size(); i++) {
            FourMatrix WMMSeed = WMMSeeds.get(i);
            // First pair
            List<List<Double>> z = EStep(WMMSeed, DNASequences);
            FourMatrix newWMM = MStep(z, DNASequences, motifWidth, background, pseudocounts);
            entropies.get(i).add(newWMM.entropy);
            // Second pair
            z = EStep(newWMM, DNASequences);
            newWMM = MStep(z, DNASequences, motifWidth, background, pseudocounts);
            entropies.get(i).add(newWMM.entropy);
            // Third pair
            z = EStep(newWMM, DNASequences);
            newWMM = MStep(z, DNASequences, motifWidth, background, pseudocounts);
            entropies.get(i).add(newWMM.entropy);

            WMMResults.add(newWMM);
        }

        // Sort WMMs by relative entropy, lowest to highest
        WMMResults.sort(new Comparator<FourMatrix>() {
            @Override
            public int compare(FourMatrix a, FourMatrix b) {
                return Double.compare(a.entropy, b.entropy);
            }
        });

        // Pick highest, median, lowest entropy WMM
        FourMatrix A = WMMResults.get(WMMResults.size() - 1);
        FourMatrix B = WMMResults.get(WMMResults.size() / 2);
        FourMatrix C = WMMResults.get(0);

        // Run 7 more E/M step pairs
        // Incidentally, this code style is atrocious
        List<FourMatrix> MoreWMMResults = new ArrayList<>();
        for (int i = 0; i < WMMResults.size(); i++) {
            FourMatrix WMM = WMMResults.get(i);

            // First pair
            List<List<Double>> z = EStep(WMM, DNASequences);
            FourMatrix newWMM = MStep(z, DNASequences, motifWidth, background, pseudocounts);
            entropies.get(i).add(newWMM.entropy);
            // Second pair
            z = EStep(newWMM, DNASequences);
            newWMM = MStep(z, DNASequences, motifWidth, background, pseudocounts);
            entropies.get(i).add(newWMM.entropy);
            // Third pair
            z = EStep(newWMM, DNASequences);
            newWMM = MStep(z, DNASequences, motifWidth, background, pseudocounts);
            entropies.get(i).add(newWMM.entropy);
            // Fourth pair
            z = EStep(newWMM, DNASequences);
            newWMM = MStep(z, DNASequences, motifWidth, background, pseudocounts);
            entropies.get(i).add(newWMM.entropy);
            // Fifth pair
            z = EStep(newWMM, DNASequences);
            newWMM = MStep(z, DNASequences, motifWidth, background, pseudocounts);
            entropies.get(i).add(newWMM.entropy);
            // Sixth pair
            z = EStep(newWMM, DNASequences);
            newWMM = MStep(z, DNASequences, motifWidth, background, pseudocounts);
            entropies.get(i).add(newWMM.entropy);
            // Seventh pair
            z = EStep(newWMM, DNASequences);
            newWMM = MStep(z, DNASequences, motifWidth, background, pseudocounts);
            entropies.get(i).add(newWMM.entropy);

            MoreWMMResults.add(newWMM);
        }

        // Sort second-round WMMs by relative entropy, lowest to highest
        MoreWMMResults.sort(new Comparator<FourMatrix>() {
            @Override
            public int compare(FourMatrix a, FourMatrix b) {
                return Double.compare(a.entropy, b.entropy);
            }
        });

        // Get the best of the second-round WMMs
        FourMatrix S = MoreWMMResults.get(WMMResults.size() - 1);

        // Print out entropies
        System.out.println("Entropy matrix:");
        for (int i = 0; i < entropies.size(); i++) {
            String entropyLine = "";
            for (int j = 0; j < entropies.get(i).size(); j++) {
                entropyLine += String.format("%-8.3f ", entropies.get(i).get(j));
            }
            System.out.println(entropyLine);
        }
        System.out.println();

        // Print out the top picks
        System.out.println("Frequency matrix for A:\n" + A.freqMatrix);
        System.out.println("Frequency matrix for B:\n" + B.freqMatrix);
        System.out.println("Frequency matrix for C:\n" + C.freqMatrix);
        System.out.println("Frequency matrix for S:\n" + S.freqMatrix);

        // Generate histograms
        // TODO
    }

    /**
     * Generates a count matrix given the list of sequences (and a name->sequence map)
     *
     * @param seqs    list of DNA sequences to include
     * @param weights amount to weigh each given sequence by
     * @return count matrix for given DNA sequences
     */
    public static FourMatrix makeCountMatrix(List<String> seqs, List<Double> weights) {
        if (seqs.size() != weights.size()) {
            throw new IllegalArgumentException("seqs and weights have different lengths");
        }

        FourMatrix mx = new FourMatrix();
        for (int i = 0; i < seqs.size(); i++) {
            mx.addSequence(seqs.get(i), weights.get(i));
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
        return countMatrix.convertToFrequency();
    }

    /**
     * Computes a weight matrix (and associated relative entropy) given a frequency matrix and background distribution
     *
     * @param frequencyMatrix frequency matrix to use
     * @param background      background distribution
     * @return WM for given params with relative entropy value
     */
    public static FourMatrix makeWMM(FourMatrix frequencyMatrix, FourTuple background) {
        FourMatrix result = frequencyMatrix.wmmScale(background);
        result.freqMatrix = frequencyMatrix;
        return result;
    }

    /**
     * Scores a list of given sequences using the given WMM
     *
     * @param WMMMatrix WMM to score with
     * @param seqs      list of sequences to score along
     * @return a map from DNA sequence to list of scores
     * (the score at index i = score of sequence of length k beginning at index i of dna sequence)
     */
    public static Map<String, List<Double>> scanWMM(FourMatrix WMMMatrix, List<String> seqs) {
        Map<String, List<Double>> scores = new HashMap<>();
        for (String seq : seqs) {
            scores.put(seq, scanSequence(WMMMatrix, seq));
        }
        return scores;
    }

    /**
     * Returns a list of scores for a single sequence of equal or longer length than the WMM
     *
     * @param WMMMatrix WMM to score on
     * @param seq       sequence to score
     * @return list of scores where get(i) = score of subsequence beginning at index i
     */
    private static List<Double> scanSequence(FourMatrix WMMMatrix, String seq) {
        List<Double> scores = new ArrayList<Double>();
        for (int i = 0; i < seq.length() - WMMMatrix.length() + 1; i++) {
            scores.add(WMMMatrix.score(seq.substring(i, i + WMMMatrix.length())));
        }
        return scores;
    }

    /**
     * MEME E-Step:
     * Calculate Z_ij (variable 1=motif in seq # 1 starts at position j, 0 = not)
     *
     * @param WMMMatrix WMM to use
     * @param seqs      sequences to calculate off of
     * @return List where list.get(i) gets variables for ith sequence, list.get(i).get(j) gets Z_ij
     */
    public static List<List<Double>> EStep(FourMatrix WMMMatrix, List<String> seqs) {
        Map<String, List<Double>> scannedScores = scanWMM(WMMMatrix, seqs);

        List<List<Double>> z = new ArrayList<>();
        for (String seq : seqs) {
            // Add all the probabilities (unscaled)
            double sum = 0;
            List<Double> seqList = new ArrayList<>();
            for (Double score : scannedScores.get(seq)) {
                Double newValue = Math.pow(2, score);
                seqList.add(newValue);
                sum += newValue;
            }

            // Scale the probabilities so they sum to 1 for each sequence
            for (int i = 0; i < seqList.size(); i++) {
                seqList.set(i, seqList.get(i) / sum);
            }

            z.add(seqList);
        }

        return z;
    }

    /**
     * MEME M-Step
     * Calculate best WMM given weighting of subsequences
     *
     * @param z    Weight of each subsequence
     * @param seqs All input DNA sequences
     */
    public static FourMatrix MStep(List<List<Double>> z, List<String> seqs, int motifWidth, FourTuple background, FourTuple pseudocounts) {
        List<String> subsequences = new ArrayList<>();
        List<Double> weights = new ArrayList<>();
        // Build up the list of subsequences and weights for each subsequence to input into makeCountMatrix
        for (int i = 0; i < seqs.size(); i++) {
            String seq = seqs.get(i);
            for (int j = 0; j <= seq.length() - motifWidth; j++) {
                subsequences.add(seq.substring(j, j + motifWidth));
                weights.add(z.get(i).get(j));
            }
        }

        FourMatrix countMatrix = makeCountMatrix(subsequences, weights);
        FourMatrix totalCounts = addPseudo(countMatrix, pseudocounts);
        FourMatrix frequencyMatrix = makeFrequencyMatrix(totalCounts);
        return makeWMM(frequencyMatrix, background);
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
         * List of nucleotide counts for each position where
         * counts.get(i) = counts over all included sequences for position i
         */
        private List<FourTuple> counts;

        /**
         * -1 if this matrix has no appropriate entropy value associated,
         * positive value if it does (and if so, equals the relative entropy value)
         */
        public double entropy = -1;

        /**
         * Null if this matrix has no appropriate frequency matrix associated,
         * nonnull if it does (and if so, equals the frequency matrix)
         */
        public FourMatrix freqMatrix = null;

        public FourMatrix() {
            this.counts = new ArrayList<>();
        }

        /**
         * Adds a DNA sequence to the count matrix. If already there, replaces previous value.
         *
         * @param seq    sequence to add
         * @param weight amount to weigh this sequence by
         */
        public void addSequence(String seq, Double weight) {
            for (int i = 0; i < seq.length(); i++) {
                if (counts.size() <= i) {
                    // Haven't already added this position
                    FourTuple count = new FourTuple();
                    count.addBase(seq.charAt(i), weight);
                    counts.add(count);
                } else {
                    // Have already had this position before: increment
                    counts.get(i).addBase(seq.charAt(i), weight);
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
            FourMatrix result = new FourMatrix();
            for (int i = 0; i < this.counts.size(); i++) {
                result.counts.add(new FourTuple());
                result.counts.set(i, result.counts.get(i).add(this.counts.get(i)));
                result.counts.set(i, result.counts.get(i).add(vector));
            }
            return result;
        }

        /**
         * Return a frequency(this) matrix
         *
         * @return conversion of this (count-style) to a (frequency-style) matrix
         */
        public FourMatrix convertToFrequency() {
            FourMatrix result = new FourMatrix();
            for (int i = 0; i < this.counts.size(); i++) {
                result.counts.add(this.counts.get(i).frequencyNormalize());
            }
            return result;
        }

        /**
         * Return a WMM(this frequency matrix) given background frequency and pseudocounts.
         * Also sets entropy value :)
         *
         * @param background background frequency
         * @return WMM(this)
         */
        public FourMatrix wmmScale(FourTuple background) {
            // calculate WM values
            FourMatrix result = new FourMatrix();
            for (int i = 0; i < this.counts.size(); i++) {
                result.counts.add(this.counts.get(i).wmmScale(background));
            }

            // calculate relative entropy
            result.entropy = 0;
            for (int i = 0; i < result.counts.size(); i++) {
                result.entropy += result.counts.get(i).indivEntropy(this.counts.get(i));
            }

            return result;
        }

        /**
         * Returns a score for this WMM on a single sequence (matching length)
         *
         * @param sequence DNA sequence with the exact same length as this WMM
         * @return score of that sequence on this WMM
         */
        public double score(String sequence) {
            if (this.length() != sequence.length()) {
                throw new RuntimeException("Sequence length != WMM length");
            }

            double score = 0;
            for (int i = 0; i < sequence.length(); i++) {
                score += this.counts.get(i).getBase(sequence.charAt(i));
            }
            return score;
        }

        /**
         * @return the length of the sequences this matrix encodes
         */
        public int length() {
            return this.counts.size();
        }

        @Override
        public String toString() {
            String result = String.format("%-8s %-8s %-8s %-8s\n",
                    "A", "C", "G", "T");
            for (int i = 0; i < this.length(); i++) {
                result += String.format("%-8.3f %-8.3f %-8.3f %-8.3f\n",
                        this.counts.get(i).a[0],
                        this.counts.get(i).a[1],
                        this.counts.get(i).a[2],
                        this.counts.get(i).a[3]);
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
        public void addBase(char base, double weight) {
            if (base == 'A') {
                this.a[0] += weight;
            } else if (base == 'C') {
                this.a[1] += weight;
            } else if (base == 'G') {
                this.a[2] += weight;
            } else if (base == 'T') {
                this.a[3] += weight;
            }
        }

        /**
         * Returns the number associated with a specific base
         *
         * @param base DNA base in standard form
         * @return double associated with that base in this vector
         */
        public double getBase(char base) {
            if (base == 'A') {
                return this.a[0];
            } else if (base == 'C') {
                return this.a[1];
            } else if (base == 'G') {
                return this.a[2];
            } else { // if (base == 'T')
                return this.a[3];
            }
        }

        /**
         * Return a frequency-normalized version of this vector
         * (divide all elements by the sum of all elements)
         *
         * @return normalized vector
         */
        public FourTuple frequencyNormalize() {
            double sum = this.a[0] + this.a[1] + this.a[2] + this.a[3];
            FourTuple result = new FourTuple();
            for (int i = 0; i < 4; i++) {
                result.a[i] = this.a[i] / sum;
            }
            return result;
        }

        /**
         * Return a vector scaled wmm-style
         *
         * @param background background frequencies
         * @return WMM-scaled vector
         */
        public FourTuple wmmScale(FourTuple background) {
            FourTuple result = new FourTuple();
            for (int i = 0; i < 4; i++) {
                // Calculate log2((actual+pseudo)/background) in each cell
                result.a[i] = Math.log(this.a[i] / background.a[i]) / Math.log(2);
            }
            return result;
        }

        /**
         * Calculate relative entropy for an individual column
         *
         * @param frequency frequency of each base
         * @return relative entropy of this individual base coord
         */
        public double indivEntropy(FourTuple frequency) {
            double result = 0;
            for (int i = 0; i < 4; i++) {
                result += this.a[i] * frequency.a[i];
            }
            return result;
        }

        @Override
        public String toString() {
            return "(" + this.a[0] + "," + this.a[1] + "," + this.a[2] + "," + this.a[3] + ")";
        }
    }
}
