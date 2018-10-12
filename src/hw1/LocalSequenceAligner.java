package hw1;

import java.util.*;

public class LocalSequenceAligner {

	private static final int GAP_PENALTY = -4;
	private static final char GAP_CHARACTER = '-';
	private static final int CONSOLE_WIDTH = 60;

	// Returns a linked list representing trace-back of the highest substring score
	private static LinkedNode align(BLOSUM62Table table, String seq1, String seq2) {
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
		for (int sum = 2; sum <= seq1.length() + seq2.length(); sum++) {
			for (int i = Math.max(1, sum - seq1.length() - 1); i < sum && i <= seq1.length(); i++) {
				int j = sum - i;
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

		return maxScoreFound;
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
			} else if (table.lookup(subseq1.charAt(i), subseq2.charAt(i)) > 0) {
				// Positive substitution score
				similar.append('+');
			} else {
				// Nothing we want to mark
				similar.append(' ');
			}
		}

		// Print this stuff to the console!
		for (int i = 0; i < subseq1.length(); i += CONSOLE_WIDTH) {
			System.out.printf("%8s:%4d %s\n", name1, index1,
					subseq1.substring(i, Math.min(i + CONSOLE_WIDTH, subseq1.length())));
			System.out.printf("%13s %s\n", "", similar.substring(i, Math.min(i + CONSOLE_WIDTH, similar.length())));
			System.out.printf("%8s:%4d %s\n", name2, index2,
					subseq2.substring(i, Math.min(i + CONSOLE_WIDTH, subseq2.length())));
			System.out.println();
		}
	}

	public static void main(String[] args) {
		BLOSUM62Table table = new BLOSUM62Table();
		LinkedNode result = LocalSequenceAligner.align(table, "KEVLAR", "KNIEVIL");
		LocalSequenceAligner.printAlignment(table, result, "name1", "name2", "KEVLAR", "KNIEVIL");
		// TODO empirical p-value calculation using sequence permutation
		// TODO actual main method control with permutation count and seeds
		// TODO
	}

}
