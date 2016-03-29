
import java.util.*;

public class PlayerSkeleton {


	/** We need to go through the field and fill it with new 1's at the 
	    positions where the piece was dropped, then return the new field */
	public int[][] dropPiece(int orient, int slot, int piece, int[][] currField){

		return currField;

	}

	/** Calculate the new fields score based on our awesome 
		heuristics */
	public int calcScore(int[][] newField){

		// Our magic numbers :)
		double a = -0.510066, b = 0.760666, c = -0.35663, d = -0.184473;
		int score = 0;

		//int height = calcHeight(newField); - Ryan
		//int linesCompleted = calcLinesCompleted(newField) -- James
		//int numHoles = calcNumHoles(newField) -- Michael
		//int bumpiness = calcBumpiness(newField) -- Andres

		// Multiply weights with scores and add together...
		return score;

	}


	//implement this function to have a working system
	public int pickMove(State s, int[][] legalMoves) {

		int orient = 0, slot = 0, rowCounter = 0,piece;
		TreeMap<Integer, Integer> scores = new TreeMap<Integer, Integer>();

		//System.out.println("** Piece " + s.getNextPiece());

		// Go through every possible move
		for(int i = 0; i < legalMoves.length; i++){

			for(int j = 0; j < legalMoves[i].length; j++){

				orient = legalMoves[i][j];
				slot = legalMoves[i][j];

			}

			// Get the representation of the field when the piece is dropped and collisions calculated
			int[][] newField = dropPiece(orient, slot, s.getNextPiece(), s.getField());
			// How well does this state perform?
			int score = calcScore(newField);
			// We map this score to the current move (which is defined by row index)
			scores.put(score, rowCounter++);

		}

		return scores.firstKey(); // We want the move associated to the highest score
	
	}
	
	public static void main(String[] args) {
		State s = new State();
		new TFrame(s);
		PlayerSkeleton p = new PlayerSkeleton();
		while(!s.hasLost()) {
			s.makeMove(p.pickMove(s,s.legalMoves()));
			s.draw();
			s.drawNext(0,0);
			try {
				Thread.sleep(300);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		System.out.println("You have completed "+s.getRowsCleared()+" rows.");
	}
	
}
