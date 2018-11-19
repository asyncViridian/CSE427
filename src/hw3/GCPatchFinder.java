package hw3;

public class GCPatchFinder {
    private static final double BEGIN_STATE1 = 0.9999;
    private static final double BEGIN_STATE2 = 0.0001;
    private static GCPatchHMM currentHMM = new GCPatchHMM(
            new GCPatchHMM.EmitProbabilities(0.25, 0.25, 0.25, 0.25),
            new GCPatchHMM.EmitProbabilities(0.2, 0.3, 0.3, 0.2),
            0.9999, 0.0001, 0.01, 0.99);

    public static void main(String[] args) {

    }
}
