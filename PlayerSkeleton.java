
import java.util.*;
import java.io.*;


public class PlayerSkeleton {

	/** Calculate number of filled lines in our theoretical field */
	public int calcLinesCleared(int[][] newField){

		int linesCleared = 0;

		for(int i = 0; i < State.ROWS; i++ ){

			boolean isCompleted = true;

			for(int j = 0; j < State.COLS; j++){
				if(newField[i][j] == 0)
					isCompleted = false;
			}
			if(isCompleted)
				linesCleared++;
		}
		return linesCleared;
	}


	public int calcHeight(int[][] newField){
		int height = 0;
		int maxRow = newField.length;
		int maxCol = newField[0].length;

		//start in the top left of the board
		//go down each colom untill you find the first peice
		//add that to the height
		for(int col = maxCol - 1; col >= 0; col--){
			for(int row = maxRow - 1; row >=0; row--){
				if(newField[row][col] != 0){
					//System.out.println("row, col "+ row +" "+ col);
					height += (row + 1);
					break;
				}
			}
		}
		//System.out.println(height);
		return height;
	}

	/** Copied from class State. Update the fields 0's and 1's based on current theoretical move 

	Would rather copy a state and use function copiedState.makeMove() but can't make deep copy
	since we cannot modify class State */
	public void makeTheoreticalMove(final int orient, final int slot, int[][] field, final int[][] pWidth, final int[][] pHeight, 
									final int[][][] pTop, final int[][][] pBottom, final int[] top, final int nextPiece, int turn) {
		
		turn++;

		//height if the first column makes contact
		int height = top[slot]-pBottom[nextPiece][orient][0];

		//for each column beyond the first in the piece
		for(int c = 1; c < pWidth[nextPiece][orient]; c++) {
			height = Math.max(height,top[slot+c]-pBottom[nextPiece][orient][c]);
		}

		//check if game ended
		if(height+pHeight[nextPiece][orient] >= State.ROWS) {
			System.out.println("We lost... This state sucks. Crash right now, better fix inc");
		}
		
		//for each column in the piece - fill in the appropriate blocks
		for(int i = 0; i < pWidth[nextPiece][orient]; i++) {
			//from bottom to top of brick
			for(int h = height+pBottom[nextPiece][orient][i]; h < height+pTop[nextPiece][orient][i]; h++) {
				field[h][i+slot] = turn;
			}
		}

		//adjust top
		for(int c = 0; c < pWidth[nextPiece][orient]; c++) {
			top[slot+c] = height+pTop[nextPiece][orient][c];
		}

	}

	/** Calculate the new fields score based on our awesome 
		heuristics */
	public double calcScore(final int[][] newField){

		// Our magic numbers :)
		double a = -0.510066, b = 0.760666, c = -0.35663, d = -0.184473;
		int score = 0;

		int height = calcHeight(newField); //- Ryan
		int linesCompleted = calcLinesCleared(newField);
		//int numHoles = calcNumHoles(newField); //-- Michael
		//int bumpiness = calcBumpiness(newField) -- Andres
		//int height = calcHeight(newField);
		//Multiply weights with scores and add together...
		//System.out.println(" lines completed of new field " + linesCompleted);

		return height*a + linesCompleted*b;
	
	}

	/**  Performs a deep copy of our old state to new state */
	public static int[][] deepCopy(int[][] original){

		int[][] result = new int[original.length][];

		for(int i = 0; i < original.length; i++){

			result[i] = Arrays.copyOf(original[i], original[i].length);

		}

		return result;

	}

	//implement this function to have a working system
	public int pickMove(State s, int[][] legalMoves) {

		int orient, slot, piece, rowCounter = 0;
		TreeMap<Double, Integer> scores = new TreeMap<Double, Integer>();

		System.out.println("** Curr piece " + s.getNextPiece());

		// Go through every possible move
		for(int i = 0; i < legalMoves.length; i++){

			orient = legalMoves[i][0];
			slot = legalMoves[i][1];
			
			// Perform deep copies, we don't want to mess up our original board
			int[][] newField = deepCopy(s.getField());
			int[] newTop = Arrays.copyOf(s.getTop(), s.getTop().length);

			// Get the representation of the field when the piece is dropped and collisions calculated (newField will change and passed as a pointer)
			makeTheoreticalMove(orient, slot, newField, s.getpWidth(), s.getpHeight(), s.getpTop(), s.getpBottom(), newTop, s.getNextPiece(), s.getTurnNumber());
			// How well does this state perform?
			double score = calcScore(newField);
			// We map this score to the current move (which is defined by row index)
			scores.put(score, rowCounter++);

		}
			
		/*
		for(Map.Entry<Double,Integer> entry : scores.entrySet()) {
		  Double key = entry.getKey();
		  Integer value = entry.getValue();

		  System.out.println(key + " => " + value);
		}

		System.out.println("Making move " + scores.get(scores.lastKey()));*/

		return scores.get(scores.lastKey()); // We want the move associated to the highest score
	
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
