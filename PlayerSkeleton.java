
import java.util.*;
import java.io.*;


public class PlayerSkeleton {

	/** Calculate number of filled lines in our theoretical field */
	public int calcLinesCleared(int[][] newField){

		int linesCleared = 0;

		for(int i = 0; i < State.ROWS; i++ ){

			int isCompleted = 1;

			for(int j = 0; j < State.COLS; j++){
				if(newField[i][j] == 0){
					isCompleted = 0;
					break;
				}
			}
			linesCleared += isCompleted;
		}

		return linesCleared;
	}

	/** Calculate number of holes in our theoretical field */
	public int calcNumHoles(int[][] newField) {
		
		int holes = 0;

		for (int j = 0; j < State.COLS; j++){

			// Boolean representation of if we have found the columns top yet
			int top = 0;

			// Do we need to check top row? If we're there we lost right?
			for (int i = State.ROWS - 1; i >= 0; i--) {

				// If the top filled square of the column is not yet found and we found a filled square 
				if (newField[i][j] != 0 && top == 0) {
					// Set top to found
					top = 1;
				}

				// If the column's top has been found and we have found an empty square
				if (top == 1 && newField[i][j] == 0) {
					//It must be a hole, so increment
					holes += 1;
				}
			}
		}

		return holes;
	}

	/** Calculate bumpiness in our theoretical field */
	public int calcBumpiness(int[][] newField){

		int bumpiness = 0;
		int lRowHeight = 0, rRowHeight = 0;
		
		for(int i = 0; i < State.COLS; i++ ){
			for(int j = State.ROWS-1; j >= 0; j--){
				if(newField[j][i] != 0 || j == 0) {
					rRowHeight = j+(newField[j][i] != 0 ? 1:0);
					if(i!= 0) {
						bumpiness += Math.abs(lRowHeight-rRowHeight);
					}
					lRowHeight = rRowHeight;
					// System.out.printf("%d, ", rRowHeight);
					break;
				}
			}
			
			
		}
		// System.out.println();
		return bumpiness;
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
		try{
		//for each column in the piece - fill in the appropriate blocks
		for(int i = 0; i < pWidth[nextPiece][orient]; i++) {
			//from bottom to top of brick
			for(int h = height+pBottom[nextPiece][orient][i]; h < height+pTop[nextPiece][orient][i]; h++) {
				field[h][i+slot] = turn;
			}
		}
	}
		catch(Exception e) {

		}

		//adjust top
		for(int c = 0; c < pWidth[nextPiece][orient]; c++) {
			top[slot+c] = height+pTop[nextPiece][orient][c];
		}

		// System.out.println("Move");
		// 	for(int i = 0; i < State.ROWS; i++ ){
		// 		for(int j = 0; j < State.COLS; j++){
		// 			System.out.printf("%d ", field[i][j]);
		// 		}
		// 		System.out.println();
		// 	}

	}

	/** Calculate the new fields score based on our awesome 
		heuristics */
	public double calcScore(final int[][] newField){

		// Our magic numbers :)
		double a = -0.510066, b = 0.760666, c = -0.35663, d = -0.184473;

		double height = calcHeight(newField); //- Ryan
		double linesCompleted = calcLinesCleared(newField);
		double numHoles = calcNumHoles(newField); //-- Jimmay!!!!
		double bumpiness = calcBumpiness(newField); //-- Andres
		// System.out.println("Height " + height + " linesCompleted: " + linesCompleted + " bumpiness: " + bumpiness + " HOLES: " + numHoles);
		//int height = calcHeight(newField);
		//Multiply weights with scores and add together...
		//System.out.println(" lines completed of new field " + linesCompleted);

		return height*a + linesCompleted*b + numHoles*c + bumpiness*d;
	
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

		int orient, slot;
		// TreeMap<Double, Integer> scores = new TreeMap<Double, Integer>();
		int moveKey = 0;
		double highestScore = -10000;

		// System.out.println("** Curr piece " + s.getNextPiece());
		int pWidth[][] = s.getpWidth();
		int pHeight[][] = s.getpHeight();
		int pTop[][][] = s.getpTop();
		int pBottom[][][] = s.getpBottom();
		int top[] = s.getTop();
		int field[][] = s.getField();
		int nextPiece = s.getNextPiece();
		int turnNumber = s.getTurnNumber();
		
		// Go through every possible move
		for(int i = 0; i < legalMoves.length; i++){

			orient = legalMoves[i][0];
			slot = legalMoves[i][1];
			
			// Perform deep copies, we don't want to mess up our original board
			int[][] newField = deepCopy(field);
			int[] newTop = Arrays.copyOf(top, top.length);

			// Get the representation of the field when the piece is dropped and collisions calculated (newField will change and passed as a pointer)
			makeTheoreticalMove(orient, slot, newField, pWidth, pHeight, pTop, pBottom, newTop, nextPiece, turnNumber);
			// How well does this state perform?
			double score = calcScore(newField);
			// We map this score to the current move (which is defined by row index)
			
			// scores.put(score, i);
			
			// System.out.printf("%f ", score);
			if(score > highestScore){
				moveKey = i;
				highestScore = score;
			}

		}
			
		/*
		for(Map.Entry<Double,Integer> entry : scores.entrySet()) {
		  Double key = entry.getKey();
		  Integer value = entry.getValue();

		  System.out.println(key + " => " + value);

		}
*/
		// System.out.println("Making move " + scores.get(scores.lastKey()) + " " + scores.lastKey() + " " + moveKey + " " + highestScore);

		return moveKey; // We want the move associated to the highest score
	
	}
	
	public static void main(String[] args) {
		State s = new State();
		new TFrame(s);
		PlayerSkeleton p = new PlayerSkeleton();
		while(!s.hasLost()) {
			s.makeMove(p.pickMove(s,s.legalMoves()));
			s.draw();
			// System.out.println()
			s.drawNext(0,0);
			try {
				Thread.sleep(30);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		System.out.println("You have completed "+s.getRowsCleared()+" rows.");
	}
	
}
