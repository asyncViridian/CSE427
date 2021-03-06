package hw1;

import general.Pair;

import java.util.*;

/**
 * @author Joyce Zhou
 */
public class LocalSeqAlign {

    private static final int GAP_PENALTY = -4;
    private static final char GAP_CHARACTER = '-';
    private static final int CONSOLE_WIDTH = 60;

    // Returns a linked list representing trace-back of the highest substring score
    private static LinkedNode align(BLOSUM62Table table, String seq1, String seq2, boolean printTable) {
        Map<Pair<Integer>, LinkedNode> substringScoreTable = new HashMap<>();

        // Build the baseline max substring values
        for (int i = 0; i <= seq1.length(); i++) {
            Pair<Integer> coords = new Pair<Integer>(i, 0);
            substringScoreTable.put(coords, new LinkedNode(coords, 0, null));
        }
        for (int j = 0; j <= seq2.length(); j++) {
            Pair<Integer> coords = new Pair<Integer>(0, j);
            substringScoreTable.put(coords, new LinkedNode(coords, 0, null));
        }

        // Keep track of maximum value
        LinkedNode maxScoreFound = new LinkedNode(null, 0, null);

        // Iterate through all of the cells as they gain 3 top-left neighbors
        for (int i = 1; i <= seq1.length(); i++) {
            for (int j = 1; j <= seq2.length(); j++) {
                char seq1char = seq1.charAt(i - 1);
                char seq2char = seq2.charAt(j - 1);

                // Calculate the scores of each possible option
                LinkedNode nodeAlignBoth = substringScoreTable.get(new Pair<Integer>(i - 1, j - 1));
                int scoreAlignBoth = nodeAlignBoth.value + table.lookup(seq1char, seq2char);
                LinkedNode nodeSkipSeq1 = substringScoreTable.get(new Pair<Integer>(i, j - 1));
                int scoreSkipSeq1 = nodeSkipSeq1.value + GAP_PENALTY;
                LinkedNode nodeSkipSeq2 = substringScoreTable.get(new Pair<Integer>(i - 1, j));
                int scoreSkipSeq2 = substringScoreTable.get(new Pair<Integer>(i - 1, j)).value + GAP_PENALTY;
                int scoreWorstCase = 0;

                // Pick the best score and link together the terms
                int bestScore = Math.max(Math.max(scoreAlignBoth, scoreSkipSeq1),
                        Math.max(scoreSkipSeq2, scoreWorstCase));
                LinkedNode bestTerm;
                if (bestScore == scoreAlignBoth) {
                    bestTerm = nodeAlignBoth;
                } else if (bestScore == scoreSkipSeq1) {
                    bestTerm = nodeSkipSeq1;
                } else if (bestScore == scoreSkipSeq2) {
                    bestTerm = nodeSkipSeq2;
                } else {
                    bestTerm = null;
                }
                Pair<Integer> coords = new Pair<Integer>(i, j);
                LinkedNode newnode = new LinkedNode(coords, bestScore, bestTerm);
                substringScoreTable.put(coords, newnode);

                // Check if we have a new best score
                if (bestScore > maxScoreFound.value) {
                    maxScoreFound = newnode;
                }
            }
        }

        // Print the score table (if applicable)
        if (printTable) {
            printScoreTable(substringScoreTable, seq1, seq2, Math.round(Math.ceil(Math.log10(maxScoreFound.value))));
        }
        return maxScoreFound;
    }

    private static void printScoreTable(Map<Pair<Integer>, LinkedNode> scoreTable, String seq1, String seq2,
                                        long spacing) {
        if (seq1.length() > 15) {
            return;
        }
        if (seq2.length() > 15) {
            return;
        }

        seq1 = " " + seq1;
        seq2 = " " + seq2;

        // Print header row
        System.out.printf("  |");
        for (int i = 0; i < seq2.length(); i++) {
            System.out.printf(" %" + spacing + "C", seq2.charAt(i));
        }
        System.out.println();

        // Print separator row
        System.out.printf("--+");
        for (int i = 0; i < (spacing + 1) * seq2.length(); i++) {
            System.out.print("-");
        }
        System.out.println();

        // Print each score row
        for (int i = 0; i < seq1.length(); i++) {
            System.out.printf("%C |", seq1.charAt(i));
            for (int j = 0; j < seq2.length(); j++) {
                System.out.printf(" %" + spacing + "d", scoreTable.get(new Pair<Integer>(i, j)).value);
            }
            System.out.println();
        }

        System.out.println();
    }

    private static void printAlignment(BLOSUM62Table table, LinkedNode node, String name1, String name2, String seq1,
                                       String seq2) {
        StringBuilder subseq1 = new StringBuilder();
        int index1 = node.coords.a;
        StringBuilder subseq2 = new StringBuilder();
        int index2 = node.coords.b;

        // Build aligned form of strings in reverse order
        while (node.prev != null) {
            // Check seq1
            if (node.coords.a == node.prev.coords.a) {
                // Gap in seq1
                subseq1.append(GAP_CHARACTER);
            } else {
                // Match in seq1
                subseq1.append(seq1.charAt(node.coords.a - 1));
                index1 = node.coords.a;
            }
            // Check seq2
            if (node.coords.b == node.prev.coords.b) {
                // Gap in seq2
                subseq2.append(GAP_CHARACTER);
            } else {
                // Match in seq2
                subseq2.append(seq2.charAt(node.coords.b - 1));
                index2 = node.coords.b;
            }
            node = node.prev;
        }

        // Reverse stringbuilders
        subseq1 = subseq1.reverse();
        subseq2 = subseq2.reverse();

        // Get alignment similarities
        StringBuilder similar = new StringBuilder();
        for (int i = 0; i < subseq1.length(); i++) {
            if (subseq1.charAt(i) == subseq2.charAt(i)) {
                // Exact match
                similar.append(subseq1.charAt(i));
            } else if (subseq1.charAt(i) != '-' && subseq2.charAt(i) != '-'
                    && table.lookup(subseq1.charAt(i), subseq2.charAt(i)) > 0) {
                // Positive substitution score (no gap)
                similar.append('+');
            } else {
                // Nothing we want to mark
                similar.append(' ');
            }
        }

        // Print this stuff to the console!
        for (int i = 0; i < subseq1.length(); i += CONSOLE_WIDTH) {
            System.out.printf("%8s:%4d %s%n", name1, index1,
                    subseq1.substring(i, Math.min(i + CONSOLE_WIDTH, subseq1.length())));
            System.out.printf("%13s %s%n", "", similar.substring(i, Math.min(i + CONSOLE_WIDTH, similar.length())));
            System.out.printf("%8s:%4d %s%n", name2, index2,
                    subseq2.substring(i, Math.min(i + CONSOLE_WIDTH, subseq2.length())));
            System.out.println();
        }
    }

    private static double calculatePValue(BLOSUM62Table table, String seq1, String seq2, int foundScore, Random random, int numPermutations) {
        int numHigher = 0;
        int numTrials = 0;
        for (int i = 0; i < numPermutations; i++) {
            // Generate permutation
            String permutation = seq2;
            for (int j = seq2.length() - 1; j > 0; j--) {
                int swap = random.nextInt(j + 1);
                // Swap positions of characters at indexes j, swap
                // Note that swap < j
                if (swap != j) {
                    permutation = permutation.substring(0, swap)
                            + permutation.charAt(j)
                            + permutation.substring(swap + 1, j)
                            + permutation.charAt(swap)
                            + permutation.substring(j + 1);
                }
            }

            // Calculate score for permutation
            int permuteScore = align(table, seq1, permutation, false).value;

            // Compare permutation permuteScore against actual permuteScore; maintain tally
            if (permuteScore >= foundScore) {
                numHigher++;
            }
            numTrials++;
        }
        return ((double) numHigher + 1) / (numTrials + 1);
    }

    private static int performComparisonAnalysis(BLOSUM62Table table, String name1, String name2, String seq1,
                                                 String seq2, Random random, int numPermutations) {
        LinkedNode result = LocalSeqAlign.align(table, seq1, seq2, true);
        LocalSeqAlign.printAlignment(table, result, name1, name2, seq1, seq2);
        System.out.printf("Total score: " + result.value + "%n");
        if (numPermutations > 0) {
            System.out.printf("Estimated P-value: %e", calculatePValue(table, seq1, seq2, result.value, random, numPermutations));
        }
        System.out.println();
        return result.value;
    }

    private static String filterProteinSequence(String sequence) {
        StringBuilder filtered = new StringBuilder();
        for (int i = 0; i < sequence.length(); i++) {
            for (int j = 0; j < BLOSUM62Table.INPUT_HEADER.length; j++) {
                if (BLOSUM62Table.INPUT_HEADER[j] == sequence.charAt(i)) {
                    filtered.append(sequence.charAt(i));
                    break;
                }
            }
        }
        return filtered.toString();
    }

    /**
     * Locally aligns all FASTA-formatted sequences entered via System.in (command line) with the given parameters.
     *
     * @param args --permutations [number of permutations to run to calculate P-value]
     *             --seed [seed for RNG]
     */
    public static void main(String[] args) {
        // Usage message
        if (args.length % 2 != 0) {
            System.err.println("Usage: java LocalSeqAlign [--permutations PPP] [--seed NNN]");
            return;
        }

        BLOSUM62Table table = new BLOSUM62Table();
        Random random = new Random();

        // Read arguments parameters
        int numPermutations = 0;
        if (args.length != 0) {
            for (int i = 0; i < args.length; i += 2) {
                if (args[i].equals("--permutations")) {
                    // Set permutation count
                    numPermutations = Integer.valueOf(args[i + 1]);
                } else if (args[i].equals("--seed")) {
                    // Set random seed
                    random = new Random(Integer.valueOf(args[i + 1]));
                }
            }
        }

        // Read sequence inputs (FASTA)
        Map<String, String> proteins = new HashMap<>();
        List<String> listNames = new ArrayList<>();
        Scanner console = new Scanner(System.in);
        String line = console.nextLine();
        String proteinName = line.split("\\|")[1];
        String proteinContent = "";
        // Get the first protein
        while (console.hasNextLine()) {
            line = console.nextLine();
            if (line.length() > 0 && line.charAt(0) == '>') {
                // Beginning of a new protein entry
                proteins.put(proteinName, proteinContent);
                listNames.add(proteinName);
                proteinName = line.split("\\|")[1];
                proteinContent = "";
            } else {
                // Continuation of a protein
                proteinContent = proteinContent + filterProteinSequence(line);

            }
        }
        // Close the last protein
        proteins.put(proteinName, proteinContent);
        listNames.add(proteinName);

        // Do all alignments
        List<Integer> scores = new ArrayList<>();
        for (int i = 0; i < listNames.size(); i++) {
            for (int j = i + 1; j < listNames.size(); j++) {
                scores.add(performComparisonAnalysis(table, listNames.get(i), listNames.get(j),
                        proteins.get(listNames.get(i)), proteins.get(listNames.get(j)), random, numPermutations));
            }
        }
        System.out.println();

        // Print score matrix
        int current = 0;
        for (int i = 0; i < listNames.size(); i++) {
            for (int j = 0; j < listNames.size(); j++) {
                if (i >= j) {
                    // 0 in lower triangle and main diagonal
                    System.out.printf("%5d ", 0);
                } else {
                    System.out.printf("%5d ", scores.get(current));
                    current++;
                }
            }
            System.out.println();
        }
    }

}
