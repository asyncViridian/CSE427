package hw1;

import java.util.*;

public class LocalSequenceAligner {

	private static final int GAP_PENALTY = -4;

	public static int align(BLOSUM62Table table, String seq1, String seq2) {
		Map<Pair<Integer>, Integer> substringScoreTable = new HashMap<>();

		// Build the baseline max substring values
		for (int i = 0; i <= seq1.length(); i++) {
			substringScoreTable.put(new Pair<Integer>(i, 0), 0);
		}
		for (int j = 0; j <= seq2.length(); j++) {
			substringScoreTable.put(new Pair<Integer>(0, j), 0);
		}

		// Keep track of maximum value
		int maxScoreFound = 0;
		int maxScoreRow = 0;
		int maxScoreCol = 0;

		// Iterate through all of the cells as they gain 3 top-left neighbors
		for (int sum = 2; sum <= seq1.length() + seq2.length(); sum++) {
			for (int i = Math.max(1, sum - seq1.length() - 1); i < sum && i <= seq1.length(); i++) {
				int j = sum - i;
				char seq1char = seq1.charAt(i - 1);
				char seq2char = seq2.charAt(j - 1);

				// Calculate the scores of each possible option
				int scoreAlignBoth = substringScoreTable.get(new Pair<Integer>(i - 1, j - 1))
						+ table.lookup(seq1char, seq2char);
				int scoreSkipSeq1 = substringScoreTable.get(new Pair<Integer>(i, j - 1)) + GAP_PENALTY;
				int scoreSkipSeq2 = substringScoreTable.get(new Pair<Integer>(i - 1, j)) + GAP_PENALTY;
				int scoreWorstCase = 0;
				int bestScore = Math.max(Math.max(scoreAlignBoth, scoreSkipSeq1),
						Math.max(scoreSkipSeq2, scoreWorstCase));
				substringScoreTable.put(new Pair<Integer>(i, j), bestScore);

				// Check if we have a new best score
				if (bestScore > maxScoreFound) {
					maxScoreFound = bestScore;
					maxScoreRow = i;
					maxScoreCol = j;
				}
			}
		}

		// TODO return max score location
		return maxScoreFound;
	}

	public static void main(String[] args) {
		BLOSUM62Table table = new BLOSUM62Table();
		LocalSequenceAligner.align(table, "KEVLAR", "KNIEVIL");
		// TODO
	}

}
