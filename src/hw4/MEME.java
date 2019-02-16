package hw4;

// Joyce Zhou

import javafx.util.Pair;

import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.*;
import java.util.*;
import java.util.List;

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

        // TODO there is a bug in the E-M algos somewhere, but not sure where exactly

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
        int mA = generateHistogram(A, "A", DNASequences);
        int mB = generateHistogram(B, "B", DNASequences);
        int mC = generateHistogram(C, "C", DNASequences);
        int mS = generateHistogram(S, "S", DNASequences);
        System.out.println("Location of max c(j) for A: " + mA);
        System.out.println("Location of max c(j) for B: " + mB);
        System.out.println("Location of max c(j) for C: " + mC);
        System.out.println("Location of max c(j) for S: " + mS);
        System.out.println();

        // Do ROC stuff
        FourMatrix[] wmms = {A, B, C, S};
        generateROC(wmms, DNASequences);
    }

    /**
     * Creates a histogram of top hits for the given WMM. Also returns the most likely motif starting position.
     *
     * @param WMMMatrix WMM to use
     * @param filename  filename to save the histogram under (filename.png)
     * @param seqs      sequences to use
     * @return the most likely motif starting index
     */
    public static int generateHistogram(FourMatrix WMMMatrix, String filename, List<String> seqs) {
        filename = "histogram_" + filename + ".png";

        // Generate scores data
        Map<String, List<Double>> scores = scanWMM(WMMMatrix, seqs);
        // Initialize histogram data to all 0
        List<Integer> counts = new ArrayList<>();
        for (int i = 0; i < scores.get(seqs.get(0)).size(); i++) {
            counts.add(0);
        }
        // Find the maximum score in each
        for (String seq : scores.keySet()) {
            List<Double> s = scores.get(seq);
            int maxIndex = 0;
            Double max = s.get(0);
            for (int i = 0; i < s.size(); i++) {
                if (s.get(i) - max > 0) {
                    maxIndex = i;
                    max = s.get(i);
                }
            }

            // And increment histogram
            counts.set(maxIndex, counts.get(maxIndex) + 1);
        }

        // Find the most likely motif point
        int maxIndex = 0;
        Integer max = counts.get(0);
        for (int i = 0; i < counts.size(); i++) {
            if (counts.get(i) - max > 0) {
                maxIndex = i;
                max = counts.get(i);
            }
        }

        // draw image:
        BufferedImage image = new BufferedImage(800, 600, BufferedImage.TYPE_INT_ARGB);
        Graphics g = image.getGraphics();
        // Draw bars
        for (int i = 0; i < counts.size(); i++) {
            // Draw fill
            g.setColor(Color.CYAN);
            g.fillRect(scale(0, counts.size(), i, 20, 780),
                    550 - scale(0, max, counts.get(i), 0, 500),
                    (780 - 20) / counts.size(),
                    scale(0, max, counts.get(i), 0, 500));
            // Draw outline
            g.setColor(Color.GRAY);
            g.drawRect(scale(0, counts.size(), i, 20, 780),
                    550 - scale(0, max, counts.get(i), 0, 500),
                    (780 - 20) / counts.size(),
                    scale(0, max, counts.get(i), 0, 500));
            // Draw label if appropriate
            if (counts.get(i) > 0) {
                // Draw position #
                g.drawString("" + i,
                        scale(0, counts.size(), i, 20, 780),
                        565);
                // Draw # hits
                g.drawString(counts.get(i) + "",
                        scale(0, counts.size(), i, 20, 780),
                        547 - scale(0, max, counts.get(i), 0, 500));
            }
        }
        // Draw bottom axis
        g.setColor(Color.BLACK);
        g.drawLine(20, 550, 780, 550);
        g.drawLine(780, 555, 780, 550);
        g.drawString("" + counts.size(), 780, 575);
        for (int i = 0; i < counts.size(); i += 10) {
            g.drawLine(scale(0, counts.size(), i, 20, 780), 555,
                    scale(0, counts.size(), i, 20, 780), 550);
            g.drawString("" + i, scale(0, counts.size(), i, 20, 780), 575);
        }
        File outputFile = new File(filename);

        try {
            ImageIO.write(image, "png", outputFile);
        } catch (IOException e) {
            System.err.println("Unable to generate file " + filename);
        }

        return maxIndex;
    }

    public static void generateROC(FourMatrix[] WMMMatrix, List<String> seqs) {
        if (WMMMatrix.length != 4) {
            throw new IllegalArgumentException("array of WMMs doesnt contain exactly 4");
        }

        String filename = "ROC.png";

        Double[] fullTau = new Double[4]; // {A,B,C,S}
        Color[] colors = {Color.BLUE, Color.ORANGE, Color.GREEN, Color.RED};

        // draw image:
        BufferedImage image = new BufferedImage(800, 600, BufferedImage.TYPE_INT_ARGB);
        Graphics g = image.getGraphics();
        // draw axes and boring background stuff
        g.setColor(Color.BLACK);
        g.drawLine(15, 20, 780, 20);
        g.drawString("1", 10, 25);
        g.drawLine(15, 580, 780, 580);
        g.drawString("0", 10, 585);
        g.drawLine(20, 20, 20, 585);
        g.drawString("0", 20, 595);
        g.drawLine(780, 20, 780, 585);
        g.drawString("1", 780, 595);
        g.setColor(Color.LIGHT_GRAY);
        g.drawLine(20, 580, 780, 20);
        g.drawString("TPR", 5, 300);
        g.drawString("FPR", 400, 595);
        // draw curves and associated AUC
        fullTau[0] = drawROC(g, colors[0], "A", WMMMatrix[0], seqs, 600, 520, false);
        fullTau[1] = drawROC(g, colors[1], "B", WMMMatrix[1], seqs, 600, 535, false);
        fullTau[2] = drawROC(g, colors[2], "C", WMMMatrix[2], seqs, 600, 550, true);
        fullTau[3] = drawROC(g, colors[3], "S", WMMMatrix[3], seqs, 600, 565, false);
        File outputFile = new File(filename);

        try {
            ImageIO.write(image, "png", outputFile);
        } catch (IOException e) {
            System.err.println("Unable to generate file " + filename);
        }
    }

    private static Double drawROC(Graphics g, Color c, String name, FourMatrix WMMMatrix, List<String> seqs,
                                  int aucX, int aucY, boolean printStats) {
        // Generate scores data
        Map<String, List<Double>> scores = scanWMM(WMMMatrix, seqs);
        // Initialize list of "true" motif indexes
        // (motifs.get(i) = index of motif in seqs.get(i))
        List<Integer> motifs = new ArrayList<>();
        // Find the maximum score in each
        for (String seq : seqs) {
            List<Double> s = scores.get(seq);
            int maxIndex = 0;
            Double max = s.get(0);
            for (int i = 0; i < s.size(); i++) {
                if (s.get(i) - max > 0) {
                    maxIndex = i;
                    max = s.get(i);
                }
            }

            // And increment histogram
            motifs.add(maxIndex);
        }

        // Generate mapping from score to true-falseness
        List<Pair<Double, Boolean>> annotated = new ArrayList<>();
        int totalPos = 0;
        int totalNeg = 0;
        for (int i = 0; i < seqs.size(); i++) {
            String seq = seqs.get(i);
            for (int j = 0; j < scores.get(seq).size(); j++) {
                if (motifs.get(i).equals(j)) {
                    // If the motif in this sequence = current reading index
                    annotated.add(new Pair<>(scores.get(seq).get(j), Boolean.TRUE));
                    totalPos++;
                } else {
                    annotated.add(new Pair<>(scores.get(seq).get(j), Boolean.FALSE));
                    totalNeg++;
                }
            }
        }
        annotated.sort(new Comparator<Pair<Double, Boolean>>() {
            @Override
            public int compare(Pair<Double, Boolean> a, Pair<Double, Boolean> b) {
                return Double.compare(a.getKey(), b.getKey());
            }
        });

        // Find the minimum tau that hits at least all true positives?
        int minTauIndex = 0;
        while (annotated.get(minTauIndex).getValue().equals(Boolean.FALSE)) {
            minTauIndex++;
        }
        Double minTau = annotated.get(minTauIndex).getKey();

        // draw ROC with graphics object
        g.setColor(c);
        double auc = 0;
        int tp = 0;
        int fp = 0;
        int tn = totalNeg;
        int fn = totalPos;
        Pair<Float, Float> prev = new Pair<>(0f, 0f);
        float leftAnchor = 0;
        for (int i = annotated.size() - 1; i >= 0; i--) {
            Pair<Double, Boolean> atau = annotated.get(i);

            // Adjust our counts of each true/false pos/neg etc.
            if (atau.getValue()) {
                tp++;
                fn--;
            } else {
                fp++;
                tn--;
            }
            float tpr = 1f * tp / totalPos;
            float fpr = 1f * fp / totalNeg;

            // Also calculate AUC
            if (!atau.getValue()) {
                // if fpr increased, add to the area & set anchor
                auc += (fpr - leftAnchor) * tpr;
                leftAnchor = fpr;

            }

            // Draw the line!
            g.drawLine(scale(0, 1, prev.getKey(), 20, 780),
                    600 - scale(0, 1, prev.getValue(), 20, 580),
                    scale(0, 1, fpr, 20, 780),
                    600 - scale(0, 1, tpr, 20, 580));
            prev = new Pair<>(fpr, tpr);


            // Print stuff that we want to print
            if (printStats && atau.getKey().equals(minTau)) {
                System.out.println("For WMM C:");
                System.out.println("tau: " + atau.getKey());
                System.out.println(" TP: " + tp);
                System.out.println(" FP: " + fp);
                System.out.println(" TN: " + tn);
                System.out.println(" FN: " + fn);
                System.out.println("TPR: " + tpr);
                System.out.println("FPR: " + fpr);
                printStats = false;
            }
        }
        // Fenceposting
        g.drawLine(scale(0, 1, prev.getKey(), 20, 580),
                scale(0, 1, prev.getValue(), 20, 780),
                scale(0, 1, 1, 20, 580),
                scale(0, 1, 1, 20, 780));

        // calculate and print AC
        g.drawString(name + " WMM AUC: " + auc, aucX, aucY);

        return minTau;
    }

    private static int scale(float vL, float vH, float v, float L, float H) {
        return Math.round((1f * (v - vL) / (vH - vL)) * (H - L) + L);
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
