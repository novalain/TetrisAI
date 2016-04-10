
import java.util.*;
import java.io.*;
import java.util.Random;


public class PlayerSkeleton {
	private static Random randnum;

	/**
	 * Copied from class State. Update the fields 0's and 1's based on current theoretical move 
	 * Would rather copy a state and use function copiedState.makeMove() but can't make deep copy
	 * since we cannot modify class State
	 * @param
	 * @param
	 * @param
	 * @param
	 * @param
	 * @param
	 * @param
	 * @param
	 * @param
	 * @param
	 * @return
	 */
	public Boolean makeTheoreticalMove(final int orient, final int slot, int[][] field, final int[][] pWidth, final int[][] pHeight, 
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
			return false;
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
		return true;
	}

	/**
	 * @param newField The board state used to calculate the fitness score
	 * @param newTop An array containing the top indices of the new board
	 * @param weightFactors The weight factors for the heuristic values.
	 * @return Returns the fitness score of this board state
	 */
	public double fitnessFunction(final int[][] newField, final int[] newTop, double[] weightFactors){

		// Our magic numbers :)
		int maxRow = newField.length;
		int maxCol = newField[0].length;

		double height = 0;
		double completeLines = 0;
		double numHoles = 0;
		double bumpiness = 0;
		
		int lowestTop = 100;

		// Calculate Height and Bumpiness
		for(int col = 0; col<maxCol; col++){
			height+= newTop[col];
			if(newTop[col] < lowestTop) lowestTop = newTop[col];
			if(col!= 0) {
				bumpiness += Math.abs(newTop[col]-newTop[col-1]);
			}	
		}
		// Calculate Height and Bumpiness
		for (int j = 0; j < maxCol; j++){
			// Boolean representation of if we have found the columns top yet
			// Do we need to check top row? If we're there we lost right?
			for (int i = newTop[j]; i >= 0; i--) {
				// If the top filled square of the column is not yet found and we found a filled square 
				if(newField[i][j] == 0) {
					//It must be a hole, so increment
					numHoles += 1;
				}
			}
		}
		
		for(int i = 0; i < lowestTop; i++ ){
			int isCompleted = 1;
			for(int j = 0; j < State.COLS; j++){
				if(newField[i][j] == 0){
					isCompleted = 0;
					break;
				}
			}
			completeLines += isCompleted;
		}

		return height*weightFactors[0] + completeLines*weightFactors[1] + numHoles*weightFactors[2] + bumpiness*weightFactors[3];
	
	}

	/**
	 * Performs a deep copy of our old state to new state
	 * @param original original state
	 * @return Deep copy of the original state
	 */
	public static int[][] deepCopy(int[][] original){

		int[][] result = new int[original.length][];

		for(int i = 0; i < original.length; i++){

			result[i] = Arrays.copyOf(original[i], original[i].length);

		}

		return result;

	}

	/**
	 * Method that obtains the best moe by iterating through the list of legalMoves and calculating the 
	 * best score based on a fitness function
	 * @param s The present state of the game
	 * @param legalMoves Array representing the available legal moves
	 * @param weightFactors The parameter list representing the weight of each heuristic value on the final function
	 * @return The best move calculated by the fitness function
	 */
	public int pickMove(State s, int[][] legalMoves, double[] weightFactors) {


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
			if(makeTheoreticalMove(orient, slot, newField, pWidth, pHeight, pTop, pBottom, newTop, nextPiece, turnNumber)){
				
				// How well does this state perform?
				double score = fitnessFunction(newField, newTop, weightFactors);
				// We map this score to the current move (which is defined by row index)
				
				// scores.put(score, i);
				
				// System.out.printf("%f ", score);
				if(score > highestScore){
					moveKey = i;
					highestScore = score;
				}
			}

		}

		return moveKey; // We want the move associated to the highest score
	
	}
	/**
	 * Commented out but might used for parallel
	 */
	// public List<Output> processInputs(List<Input> inputs)
	//         throws InterruptedException, ExecutionException {

	//     int threads = Runtime.getRuntime().availableProcessors();
	//     ExecutorService service = Executors.newFixedThreadPool(threads);

	//     List<Future<Output>> futures = new ArrayList<Future<Output>>();
	//     for (final Input input : inputs) {
	//         Callable<Output> callable = new Callable<Output>() {
	//             public Output call() throws Exception {
	//                 Output output = new Output();
	//                 // process your input here and compute the output
	//                 return output;
	//             }
	//         };
	//         futures.add(service.submit(callable));
	//     }

	//     service.shutdown();

	//     List<Output> outputs = new ArrayList<Output>();
	//     for (Future<Output> future : futures) {
	//         outputs.add(future.get());
	//     }
	//     return outputs;
	// }

	/**
	 * @param fitnessP The weight factors passed into the fitness function
	 * @param sleepTime The time given to the thread to sleep between drawing each new move
	 * @param maxScore Maximum permited score for this game
	 * @return Returns the score obtained by the AI
	 */
	public static int RunAI(double[] fitnessP, int sleepTime, int maxScore){
		State s = new State();
		new TFrame(s);
		PlayerSkeleton p = new PlayerSkeleton();
		
		while(!s.hasLost() && s.getRowsCleared() < maxScore) {
			s.makeMove(p.pickMove(s,s.legalMoves(), fitnessP));
			s.draw();
			s.drawNext(0,0);
			try {
				Thread.sleep(sleepTime);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		return s.getRowsCleared();
	}
	/**
	 * Creates an initial population of randomly generated unit vectors in the 4D unit sphere for the 
	 * heuristic values
	 * @param populationCount The size of the desired population
	 * @param weightsCount Number of desired weight factors in each individual
	 * @return Returns the initial population
	 */
	public static double[][] getIntialPopulation(int populationCount, int weightsCount){
		double pEvolve[][] = new double[populationCount][weightsCount];
		int fitnessP[] = new int[populationCount];
		
		for(int i = 0; i< populationCount; i++) {
			fitnessP[i] = 0;
			double magnitude = 0.0;
			for (int j = 0; j<weightsCount; j++) {
				pEvolve[i][j] = randnum.nextDouble()*2.0 - 1.0;
				magnitude += Math.abs(pEvolve[i][j]);
			}
			for (int j = 0; j<weightsCount; j++) {
				pEvolve[i][j] /= magnitude;
				
			}
		}
		return pEvolve;
	}
	
	/**
	 * Runs the evolution with the use of genetic algorithms
	 * @param heuristicsCount Number of heurstic weights to search for
	 * @return The result of the evolution, namely the heuristic value weights
	 */
	public static double[] RunEvolution(int heuristicsCount){
		System.out.println("Starting Evolution");
		int numGames= 3;
		int populationCount = 100;
		double finalParameters[] = new double[heuristicsCount];

		double pEvolve[][] = getIntialPopulation(populationCount, heuristicsCount);
		int fitnessP[] = new int[populationCount];

		// This will be run in parallel

		for(int i = 0; i< populationCount; i++){
			System.out.println("Evolving "+ pEvolve[i][0] + " " + pEvolve[i][1] + " " + pEvolve[i][2] + " " + pEvolve[i][3]);
			for (int j = 0; j<numGames; j++) {
				int rowsCleared = RunAI(pEvolve[i], 30, 200);
				System.out.println("Game " + (j+1)+ ": "+ rowsCleared);
				fitnessP[i] += rowsCleared;
			}
			System.out.println("Total " + fitnessP[i]);
		}

		// End parallel
		Random rand = new Random(); 
		LinkedList crossOver = new LinkedList();
		while(crossOver.size() <= populationCount/3) {
			// select 10%
			// select fittest individuals and add them to the crossOver population
		}
		// Continue with evolution
		return finalParameters;
	}

	/**
	 * Main method for running either the game or the evolution
	 *
	 * @param args If an argument is passed(Whatever it is) it will run the evolution, otherwise it will run the game
	 */
	public static void main(String[] args) {
		randnum = new Random();
		double weightFactors[] = {-0.510066, 0.760666, -0.35663, -0.184473};
		if(args.length > 0){
			weightFactors = RunEvolution(4);
			// Returns empty, not yet finished
			System.out.println("Completed evolution with: " + weightFactors[0] + " " + weightFactors[1] + " " + weightFactors[2] + " " + weightFactors[3]);
		}
		System.out.println("You have completed "+ RunAI(weightFactors, 300, Integer.MAX_VALUE)+" rows.");
	}
	
}
