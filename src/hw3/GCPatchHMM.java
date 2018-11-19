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
    }
}
