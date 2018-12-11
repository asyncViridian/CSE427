package hw3;

import general.Pair;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Iterator;
import java.util.Stack;

public class GCPatchFinder {
    public static final double BEGIN_STATE1 = 0.9999;
    public static final double BEGIN_STATE2 = 0.0001;

    public static void main(String[] args) throws IOException {
        if (args.length != 1) {
            System.err.println("Usage: GCPatchFinder (filename)");
            return;
        }

        BufferedReader in = new BufferedReader(new FileReader(args[0]));
        // Get FASTA format input: assume properly formatted
        // Load the sequence into input iterator
        InputIterator seqIterator = new InputIterator(in);

        // Initialize our HMM
        GCPatchHMM currentHMM = new GCPatchHMM(
                new GCPatchHMM.EmitProbabilities(0.25, 0.25, 0.25, 0.25),
                new GCPatchHMM.EmitProbabilities(0.2, 0.3, 0.3, 0.2),
                0.9999, 0.0001, 0.01, 0.99);

        for (int i = 0; i < 9; i++) {
            // Perform Viterbi traceback on input
            HMMTrainer trainer = new HMMTrainer();
            viterbiTraceback(currentHMM, seqIterator, 5, trainer);
            System.out.println();
            in.close();
            in = new BufferedReader(new FileReader(args[0]));
            seqIterator = new InputIterator(in);

            // And complete Viterbi training (iterate!)
            currentHMM = trainer.getModel();
        }

        // Run our final analysis
        viterbiTraceback(currentHMM, seqIterator, -1, new HMMTrainer());
    }

    private static void viterbiTraceback(GCPatchHMM hmm,
                                                Iterator<Character> inputSeq,
                                                int numResultsToPrint,
                                                HMMTrainer trainer) {
        System.out.println(hmm);
        double startTime = System.nanoTime();
        if (!inputSeq.hasNext()) {
            // No input at all. >:[
            System.err.println("No properly formatted input found");
            return;
        }

        // Initialize states with beginning information
        // currentStep[0] is the state for background, [1] for GC
        Character first = inputSeq.next();
        GCPatchHMM.HMMState[] states = {GCPatchHMM.HMMState.BACKGROUND, GCPatchHMM.HMMState.GC};
        ViterbiNode[] currentStep = {
                new ViterbiNode(null, states[0],
                        Math.log(GCPatchFinder.BEGIN_STATE1
                                * hmm.emitProb(first, states[0])),
                        1, first),
                new ViterbiNode(null, states[1],
                        Math.log(GCPatchFinder.BEGIN_STATE2
                                * hmm.emitProb(first, states[1])),
                        1, first)};

        // Step through input seq
        while (inputSeq.hasNext()) {
            Character c = inputSeq.next();
            ViterbiNode[] nextStep = {
                    calculateNext(currentStep, hmm, states[0], c),
                    calculateNext(currentStep, hmm, states[1], c)};

            // And update steps
            currentStep = nextStep;
        }

        // Get the best result
        ViterbiNode result;
        if (currentStep[0].logProb > currentStep[1].logProb) {
            result = currentStep[0];
        } else {
            result = currentStep[1];
        }

        // Print out the result information
        printResultInfo(result, numResultsToPrint, trainer);
        System.out.println("Time taken for this iteration: " +
                (System.nanoTime() - startTime) + " nanoseconds");
    }

    /**
     * Return the next Viterbi node based off of a set of previous nodes, a HMM,
     * the type of state the next node will be, and an emission.
     *
     * @param currentStep set of previous nodes
     * @param hmm         HMM to calculate based off of
     * @param targetState type of node
     * @param c           emission
     * @return new ViterbiNode with the correct log_e(probability)
     */
    private static ViterbiNode calculateNext(ViterbiNode[] currentStep,
                                             GCPatchHMM hmm,
                                             GCPatchHMM.HMMState targetState,
                                             Character c) {
        ViterbiNode bestNode = null;
        double bestLogProb = 0;

        // Check each input state
        for (ViterbiNode testNode : currentStep) {
            double transitionProb = hmm.transProb(testNode.state, targetState);
            double emissionProb = hmm.emitProb(c, targetState);
            double testLogProb = testNode.logProb + Math.log(transitionProb * emissionProb);

            // If not yet initialized or it beats the previous winner, set new winner
            if (bestNode == null || testLogProb > bestLogProb) {
                bestNode = testNode;
                bestLogProb = testLogProb;
            }
        }

        // Return result!
        return new ViterbiNode(bestNode, targetState, bestLogProb, bestNode.index + 1, c);
    }

    private static void printResultInfo(ViterbiNode n, int k, HMMTrainer trainer) {
        System.out.println("log(probability) = \t" + n.logProb);

        // Read through and get hits
        Stack<Pair<Integer>> hits = new Stack<>();
        boolean inHit = false;
        int lastEnding = 1;
        while (n != null) {
            trainer.addToTrainer(n);
            if (inHit && n.state == GCPatchHMM.HMMState.BACKGROUND) {
                // A hit just ended!
                hits.push(new Pair<>(n.index + 1, lastEnding));
                lastEnding = -1;
                inHit = false;
            } else if (!inHit && n.state == GCPatchHMM.HMMState.GC) {
                // A hit begins!
                lastEnding = n.index;
                inHit = true;
            }
            n = n.prev;
        }
        // Fencepost the last case
        if (inHit) {
            hits.push(new Pair<>(1, lastEnding));
        }

        // Print out all the hits we got
        System.out.println("Total hits = \t" + hits.size());
        if (k == -1 || k > hits.size()) {
            // Print all the hits if (k == -1) or if there are fewer than k hits
            k = hits.size();
        }
        System.out.printf("|%-8s|%-8s|%-8s|\n", "Start", "End", "Length");
        System.out.println("|--------|--------|--------|");
        for (int i = 0; i < k; i++) {
            Pair<Integer> hit = hits.pop();
            System.out.printf("|%-8d|%-8d|%-8d|\n",
                    hit.a, hit.b, hit.b - hit.a + 1);
        }
    }

    private static class HMMTrainer {
        // Note: 1 because we use pseudocounts!

        // Counts of each emission in a state (A,C,G,T respectively)
        private int[] backgroundEmit = {1, 1, 1, 1};
        private int[] GCEmit = {1, 1, 1, 1};
        // Counts of each transition (B>B,B>G,G>B,G>G respectively)
        private int[] numTransitions = {1, 1, 1, 1};

        public void addToTrainer(ViterbiNode n) {
            // increment emit
            if (n.state == GCPatchHMM.HMMState.BACKGROUND) {
                if (n.emitted == 'A') {
                    backgroundEmit[0]++;
                } else if (n.emitted == 'C') {
                    backgroundEmit[1]++;
                } else if (n.emitted == 'G') {
                    backgroundEmit[2]++;
                } else if (n.emitted == 'T') {
                    backgroundEmit[3]++;
                }
            } else { //n.state== GCPatchHMM.HMMState.GC
                if (n.emitted == 'A') {
                    GCEmit[0]++;
                } else if (n.emitted == 'C') {
                    GCEmit[1]++;
                } else if (n.emitted == 'G') {
                    GCEmit[2]++;
                } else if (n.emitted == 'T') {
                    GCEmit[3]++;
                }
            }

            // increment transition
            if (n.state == GCPatchHMM.HMMState.BACKGROUND) {
                if (n.prev != null) {
                    if (n.prev.state == GCPatchHMM.HMMState.BACKGROUND) {
                        numTransitions[0]++;
                    } else if (n.prev.state == GCPatchHMM.HMMState.GC) {
                        numTransitions[2]++;
                    }
                }
            } else { //n.state== GCPatchHMM.HMMState.GC
                if (n.prev != null) {
                    if (n.prev.state == GCPatchHMM.HMMState.BACKGROUND) {
                        numTransitions[1]++;
                    } else if (n.prev.state == GCPatchHMM.HMMState.GC) {
                        numTransitions[3]++;
                    }
                }
            }
        }

        public GCPatchHMM getModel() {
            double totalBEmit = backgroundEmit[0] + backgroundEmit[1] + backgroundEmit[2] + backgroundEmit[3];
            double totalGEmit = GCEmit[0] + GCEmit[1] + GCEmit[2] + GCEmit[3];
            double totalFromB = numTransitions[0] + numTransitions[1];
            double totalFromG = numTransitions[2] + numTransitions[3];
            return new GCPatchHMM(
                    new GCPatchHMM.EmitProbabilities(
                            (double) backgroundEmit[0] / totalBEmit,
                            (double) backgroundEmit[1] / totalBEmit,
                            (double) backgroundEmit[2] / totalBEmit,
                            (double) backgroundEmit[3] / totalBEmit),
                    new GCPatchHMM.EmitProbabilities(
                            (double) GCEmit[0] / totalGEmit,
                            (double) GCEmit[1] / totalGEmit,
                            (double) GCEmit[2] / totalGEmit,
                            (double) GCEmit[3] / totalGEmit),
                    numTransitions[0] / totalFromB,
                    numTransitions[1] / totalFromB,
                    numTransitions[2] / totalFromG,
                    numTransitions[3] / totalFromG);
        }
    }

    private static class ViterbiNode {
        public final ViterbiNode prev;
        public final GCPatchHMM.HMMState state;
        public final double logProb;
        public final int index;
        public final char emitted;

        public ViterbiNode(ViterbiNode prev, GCPatchHMM.HMMState state, double logProb, int index, char emitted) {
            this.prev = prev;
            this.state = state;
            this.logProb = logProb;
            this.index = index;
            this.emitted = emitted;
        }
    }

    private static class InputIterator implements Iterator<Character> {

        private final BufferedReader seq;
        private int nextChar;

        public InputIterator(BufferedReader seq) throws IOException {
            this.seq = seq;

            // Skip first line (comment). Only read the first actual gene sequence.
            this.seq.readLine();
            this.nextChar = seq.read();
        }

        @Override
        public boolean hasNext() {
            return !(nextChar == -1) && !(nextChar == '>');
        }

        @Override
        public Character next() {
            // Uppercasify
            char c = Character.toUpperCase((char) nextChar);

            // Move on the "next character"
            try {
                do {
                    nextChar = seq.read();
                } while (nextChar == '\n' || nextChar == '\r');
            } catch (IOException e) {
                e.printStackTrace();
            }

            // Default to T if it's not ACGT
            if (c != 'A' && c != 'C' && c != 'G' && c != 'T') {
                c = 'T';
            }
            return c;
        }
    }
}
