package hw3;

/**
 * A HMM to model G-C patches in a genome sequence.
 */
public class GCPatchHMM {
    private final EmitProbabilities state1;
    private final EmitProbabilities state2;
    private final double p1to1;
    private final double p1to2;
    private final double p2to1;
    private final double p2to2;

    /**
     * Create a new HMM specifically designed to model GC patches.
     *
     * @param state1 the background emission probabilities
     * @param state2 the GC patch emission probabilities
     * @param p1to1  the probability that state 1 will remain in state 1
     * @param p1to2  the probability that state 1 will transition to state 2
     * @param p2to1  the probability that state 2 will transition to state 1
     * @param p2to2  the probability that state 2 will remain in state 2
     */
    public GCPatchHMM(EmitProbabilities state1, EmitProbabilities state2,
                      double p1to1, double p1to2, double p2to1, double p2to2) {
        this.state1 = state1;
        this.state2 = state2;
        this.p1to1 = p1to1;
        this.p1to2 = p1to2;
        this.p2to1 = p2to1;
        this.p2to2 = p2to2;
    }

    /**
     * Return the probability that a character was emitted from a particular state
     *
     * @param c     the character
     * @param state the state that it was emitted in
     * @return probability P(emit | c, state)
     */
    public double emitProb(Character c, HMMState state) {
        switch (state) {
            case BACKGROUND:
                return this.state1.emitProb(c);
            case GC:
                return this.state2.emitProb(c);
            default: // This should never happen
                return -1;
        }
    }

    /**
     * Return the probability of a state transition
     *
     * @param from the state originating from
     * @param to   the state transitioning to
     * @return p(from transitions to to)
     */
    public double transProb(HMMState from, HMMState to) {
        if (from == HMMState.BACKGROUND) {
            if (to == HMMState.BACKGROUND) {
                return this.p1to1;
            } else { //to == HMMState.GC
                return this.p1to2;
            }
        } else { // from == HMMState.GC
            if (to == HMMState.BACKGROUND) {
                return this.p2to1;
            } else { //to == HMMState.GC
                return this.p2to2;
            }
        }
    }

    public String toString() {
        return String.format("HMM parameters:\n" +
                        "%-8s%-8s%-8s%-8s%-8s\n" +
                        "%-8s" + this.state1 + "\n" +
                        "%-8s" + this.state2 + "\n" +
                        "Trans\tBG\tCG\n" +
                        "%-8s%-8.5f%-8.5f\n" +
                        "%-8s%-8.5f%-8.5f\n" +
                        "%-8s%-8.5f%-8.5f\n",
                "Emit", "A", "C", "G", "T",
                "BG",
                "CG",
                "Begin", GCPatchFinder.BEGIN_STATE1, GCPatchFinder.BEGIN_STATE2,
                "BG", this.p1to1, this.p1to2,
                "CG", this.p2to1, this.p2to2);
    }

    /**
     * The possible states of the HMM
     */
    public enum HMMState {
        /**
         * State 1 of the GCPatchHMM; represents background model (no patch)
         */
        BACKGROUND,
        /**
         * State 2 of the GCPatchHMM; represents a GC patch
         */
        GC
    }

    /**
     * Emission probabilities for a particular HMM state
     */
    public static class EmitProbabilities {
        public final double pA;
        public final double pC;
        public final double pG;
        public final double pT;

        /**
         * Create a new set of emission probabilities
         *
         * @param emitA the probability of emitting A
         * @param emitC the probability of emitting C
         * @param emitG the probability of emitting G
         * @param emitT the probability of emitting T
         */
        public EmitProbabilities(double emitA, double emitC, double emitG, double emitT) {
            this.pA = emitA;
            this.pC = emitC;
            this.pG = emitG;
            this.pT = emitT;
        }

        /**
         * Return the probability of emitting character C in this set of emit probabilities
         *
         * @param character the character being emitted
         * @return p(emit | c, this state)
         */
        public double emitProb(Character character) {
            char c = Character.toUpperCase(character.charValue());
            if (c == 'A') {
                return pA;
            } else if (c == 'C') {
                return pC;
            } else if (c == 'G') {
                return pG;
            } else if (c == 'T') {
                return pT;
            } else { // User put in bad input >:[
                return 0;
            }
        }

        public String toString() {
            return String.format("%-8.5f%-8.5f%-8.5f%-8.5f", this.pA, this.pC, this.pG, this.pT);
        }
    }
}
