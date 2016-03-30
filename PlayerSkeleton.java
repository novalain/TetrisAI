
import java.util.*;

public class PlayerSkeleton {

	public int calcHeight(int[][] newField){
		int height = 0;
		int maxRow = newField.length;
		int maxCol = newField[0].length;

		//start in the top left of the board
		//go down each colom untill you find the first peice
		//add that to the height
		for(int col = maxCol - 1; col >= 0; col--){
			for(int row = maxRow - 1; row >=0; row--){
				if(newField[row][col] == 1){
					System.out.println("row, col "+ row +" "+ col);
					height += (row + 1);
					break;
				}
			}
		}
		//System.out.println(height);
		return height;
	}


	/** We need to go through the field and fill it with new 1's at the 
	    positions where the piece was dropped, this is the state we want to evaluate later */
	public int[][] dropPiece(State state, int orient, int slot){

		// Problem: can't copy the old state. want to theoretically make a move

		//State copiedState = new State();
		//copiedState = state;
		//copiedState.makeMove(orient, slot);
		return state.getField();

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
		int height = calcHeight(newField);
		// Multiply weights with scores and add together...
		return score;

	}


	//implement this function to have a working system
	public int pickMove(State s, int[][] legalMoves) {

		int orient, slot, piece, rowCounter = 0;
		TreeMap<Integer, Integer> scores = new TreeMap<Integer, Integer>();

		//System.out.println("** Piece " + s.getNextPiece());

		// Go through every possible move
		for(int i = 0; i < legalMoves.length; i++){

			orient = legalMoves[i][0];
			slot = legalMoves[i][1];

			// Get the representation of the field when the piece is dropped and collisions calculated
			int[][] newField = dropPiece(s, orient, slot);
			// How well does this state perform?
			int score = calcScore(newField);
			// We map this score to the current move (which is defined by row index)
			scores.put(score, rowCounter++);

			//System.out.println("****************");

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
