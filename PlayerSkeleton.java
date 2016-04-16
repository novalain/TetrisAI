
import java.util.*;
import java.io.*;
import java.util.Random;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

public class PlayerSkeleton {
	private static Random randnum;
	public static boolean runEvoLookahead = false;

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
	public double newFitnessFunction(final int[][] newField, final int[] newTop, double[] weights){
		int maxRow = newField.length;
		int maxCol = newField[0].length;

		double landingHeight = 0; // Done
		double rowsCleared = 0; // Done
		double rowTransitions = 0; // Done
		double columnTransitions = 0; // Done
		double numHoles = 0; // Done
		double wellSums = 0;
		int moveNumber = -1;

		for(int i = 0; i<maxCol; i++) {
			for (int j  = newTop[i]-1; j >=0; j--) {
				if(newField[j][i] == 0) numHoles++;
			}
			if(newField[Math.max(newTop[i]-1, 0)][i] > moveNumber) {
				moveNumber = newField[Math.max(newTop[i]-1, 0)][i];
				
				landingHeight = newTop[i];
			}
		}

		for(int i = 0; i<maxRow; i++) {
			boolean lastCell = false;
			boolean currentCell = false;
			int rowIsClear = 1;
			for (int j = 0; j<maxCol; j++) {
				currentCell = false;
				if(newField[i][j] == 0) {
					rowIsClear = 0;
					currentCell = true;
				}
				
				if(lastCell != currentCell) {
					rowTransitions++;
				}
				lastCell = currentCell;
			}
			rowsCleared+=rowIsClear;
			if(currentCell) rowTransitions++;
		}

		for(int i = 0; i<maxCol; i++) {
			boolean lastCell = true;
			boolean currentCell = false;
			for (int j = 0; j<maxRow-1; j++) {
				currentCell = (newField[j][i] != 0);
				if(lastCell != currentCell) {
					columnTransitions++;
				}
				lastCell = currentCell;
			}
		}

		for(int i = 1; i<maxCol-1; i++) {
			for(int j = 0; j < maxRow; j++) {
				if(newField[j][i] == 0 && newField[j][i-1] != 0 && newField[j][i+1] != 0) {
					wellSums++;
					for (int k = j -1; k >=0; k--) {
						if(newField[k][i] == 0) wellSums++;
						else break;
					}
				}
			}
		}


		for(int j = 0; j < maxRow; j++) {
			if(newField[j][0] == 0 && newField[j][1] != 0) {
				wellSums++;
				for (int k = j -1; k >=0; k--) {
					if(newField[k][0] == 0) wellSums++;
					else break;
				}
			}
			if(newField[j][maxCol-1] == 0 && newField[j][maxCol-2] != 0) {
				wellSums++;
				for (int k = j -1; k >=0; k--) {
					if(newField[k][maxCol-1] == 0) wellSums++;
					else break;
				}
			}
		}

		return landingHeight*weights[0] + rowsCleared*weights[1] + rowTransitions*weights[2] + 
		columnTransitions*weights[3] + numHoles*weights[4] + wellSums*weights[5];
		
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
	/**
	 * Method that obtains the best moe by iterating through the list of legalMoves and calculating the 
	 * best score based on a fitness function
	 * @param s The present state of the game
	 * @param legalMoves Array representing the available legal moves
	 * @param weightFactors The parameter list representing the weight of each heuristic value on the final function
	 * @return The best move calculated by the fitness function
	 */
	public int pickMove(State s, int[][] legalMoves, double[] weightFactors) {

		ArrayList<Integer> possibleMoves = new ArrayList<Integer>();

		int orient, slot;
		int moveKey = 0;
		double highestScore = -10000;

		int pWidth[][] = s.getpWidth();
		int pHeight[][] = s.getpHeight();
		int pTop[][][] = s.getpTop();
		int pBottom[][][] = s.getpBottom();
		int top[] = s.getTop();
		int field[][] = s.getField();
		int nextPiece = s.getNextPiece();
		int turnNumber = s.getTurnNumber();
		for(int i = 0; i < legalMoves.length; i++){

			orient = legalMoves[i][0];
			slot = legalMoves[i][1];
			
			// Perform deep copies, we don't want to mess up our original board
			int[][] newField = deepCopy(field);
			int[] newTop = Arrays.copyOf(top, top.length);

			// Get the representation of the field when the piece is dropped and collisions calculated (newField will change and passed as a pointer)
			if(makeTheoreticalMove(orient, slot, newField, pWidth, pHeight, pTop, pBottom, newTop, nextPiece, turnNumber)){
				double score = newFitnessFunction(newField, newTop, weightFactors);
				// Equal to highest
				if(Math.abs(score - highestScore) < 0.000000001){
					//System.out.println("Warning, scores are equal!");
					possibleMoves.add(i);
				}
				// New highest score
				else if(score > highestScore){
					possibleMoves.clear();
					possibleMoves.add(i);
					highestScore = score;
				}
			}
		}

		return possibleMoves.size() == 0 ? 0 : possibleMoves.get(randnum.nextInt(possibleMoves.size()));
	
	}


	/**
	 * Method that obtains the best moe by iterating through the list of legalMoves and calculating the 
	 * best score based on a fitness function
	 * @param s The present state of the game
	 * @param legalMoves Array representing the available legal moves
	 * @param weightFactors The parameter list representing the weight of each heuristic value on the final function
	 * @return The best move calculated by the fitness function
	 */
	public int pickMoveWithLookahead(State s, int[][] legalMoves, double[] weightFactors) {
		ArrayList<Integer> possibleMoves = new ArrayList<Integer>();
		ArrayList<Integer> possibleScores = new ArrayList<Integer>();
		IntegerDouble[] fitnessScores = new IntegerDouble[legalMoves.length];

		// init
		for(int i = 0 ; i < legalMoves.length; i++){
			fitnessScores[i] = new IntegerDouble(-10000,-10000);
		}

		int orient, slot;
		int moveKey = 0;
		double highestScore = -10000;

		int pWidth[][] = s.getpWidth();
		int pHeight[][] = s.getpHeight();
		int pTop[][][] = s.getpTop();
		int pBottom[][][] = s.getpBottom();
		int top[] = s.getTop();
		int field[][] = s.getField();
		int nextPiece = s.getNextPiece();
		int turnNumber = s.getTurnNumber();

		int[][] newField = new int[field.length][];
		int[] newTop = new int[top.length];

		// Go through every possible move
		for(int i = 0; i < legalMoves.length; i++){

			orient = legalMoves[i][0];
			slot = legalMoves[i][1];
			
			// Perform deep copies, we don't want to mess up our original board
			newField = deepCopy(field);
			newTop = Arrays.copyOf(top, top.length);

			// Get the representation of the field when the piece is dropped and collisions calculated (newField will change and passed as a pointer)
			if(makeTheoreticalMove(orient, slot, newField, pWidth, pHeight, pTop, pBottom, newTop, nextPiece, turnNumber)){
				// How well does this state perform?
				double score = newFitnessFunction(newField, newTop, weightFactors);
				IntegerDouble fitnessScore = new IntegerDouble(score, i);
				fitnessScores[i] = fitnessScore;

			} else {
				return 0;
			}

		}

		// Need this stuff copied from state classto get legal moves from a arbitrary pieace ..
		int[][][] legalMovesForPiece = new int[State.N_PIECES][][];
		for(int i = 0; i < State.N_PIECES; i++) {
			//figure number of legal moves
			int n = 0;
			for(int j = 0; j < s.getpOrients()[i]; j++) {
				//number of locations in this orientation
				n += State.COLS+1-s.getpWidth()[i][j];
			}
			//allocate space
			legalMovesForPiece[i] = new int[n][2];
			//for each orientation
			n = 0;
			for(int j = 0; j < s.getpOrients()[i]; j++) {
				//for each slot
				for(int k = 0; k < State.COLS+1-pWidth[i][j];k++) {
					legalMovesForPiece[i][n][State.ORIENT] = j;
					legalMovesForPiece[i][n][State.SLOT] = k;
					n++;
				}
			}
		}

		Arrays.sort(fitnessScores);
		int numCandidatesToLookahead = 5;
		IntegerDouble[] bestSpots = new IntegerDouble[numCandidatesToLookahead];

		for(int i = 0; i < numCandidatesToLookahead; i++){
			bestSpots[i] = fitnessScores[fitnessScores.length - 1 - i];
		}

		for(int i = 0; i < numCandidatesToLookahead; i++){
			double[] bestScoresForFutureTiles = new double[7];
			// Iterate through all pieces and store the max score of the lookahead
			for(int j = 0; j < 7; j++){
				// get legal moves for this piece
				int[][] legalMovesForLookahead = legalMovesForPiece[j];

				IntegerDouble highestScoreForLookaheadPiece = new IntegerDouble(0,0);
				highestScore = -10000;

				// Determine the best spot where these pieces can land on
				for(int k = 0; k < legalMovesForLookahead.length; k++){

					orient = legalMovesForLookahead[k][0];
					slot = legalMovesForLookahead[k][1];
					
					// Perform new copies of our already theoretical board , theroretical x2 ish
					int[][] newFieldLookahead = deepCopy(newField);
					int[] newTopLookahead = Arrays.copyOf(newTop, newTop.length);

					// Get the representation of the field when the piece is dropped and collisions calculated (newField will change and passed as a pointer)
					if(makeTheoreticalMove(orient, slot, newFieldLookahead, pWidth, pHeight, pTop, pBottom, newTopLookahead, j, ++turnNumber)){
						double score = newFitnessFunction(newFieldLookahead, newTopLookahead, weightFactors);
						if(score > highestScore){
							highestScore = score;

						}
					}
				}
				bestScoresForFutureTiles[j] = highestScore;
			}

			// Take the worst of these score and add to out best spots candidate
			Arrays.sort(bestScoresForFutureTiles);
			bestSpots[i].first += bestScoresForFutureTiles[0]; // add worst of the set

		}

		Arrays.sort(bestSpots);
		return (int)bestSpots[bestSpots.length - 1].second;	
	}


	//Start of Parallelization code
	public static class gameThread implements Callable<Integer> {
		private final double[][] pEvolve;
		private final int guy;
		private final int numGames;
		private int rowsCleared = 0;

		public gameThread(double[][] probEvolve, int i, int games) {
			this.pEvolve = probEvolve;
			this.guy = i;
			this.numGames = games;
		}

		public Integer call() {
			try {
				for (int i = 0; i < numGames; i++) {
					rowsCleared += RunAI(pEvolve[guy], 30, Integer.MAX_VALUE, false, runEvoLookahead);
				}
			} catch (Exception e) {
				e.printStackTrace();
			}
		return rowsCleared;
		}
	}

	public static IntegerPair[] parallelizePerson(double[][] pEvolve, int numGames, int populationCount){

	    IntegerPair fitnessP[] = new IntegerPair[populationCount];

	    int threads = Runtime.getRuntime().availableProcessors();
	    ExecutorService executor = Executors.newFixedThreadPool(threads);
	    ArrayList<Future<Integer>> futures = new ArrayList<Future<Integer>>();

	    for(int i = 0; i < populationCount; i++){
			fitnessP[i] = new IntegerPair(0, i);
	    	Callable<Integer> callable = new gameThread(pEvolve, i, numGames);
		    Future<Integer> future = executor.submit(callable);
		    futures.add(future);
	    }

	    executor.shutdown();

	    int j = 0;
	    int totalRows = 0;
	    for (Future<Integer> future : futures) {
	    	try {
	    		fitnessP[j].first += future.get();
	    		j += 1;
	    	} catch (InterruptedException e) {
	    		e.printStackTrace();
	    	} catch (ExecutionException e) {
	    		System.out.println(e);
	    	}
	    }
	    return fitnessP;
	}
	//End of Parallelization code


	/**
	 * @param fitnessP The weight factors passed into the fitness function
	 * @param sleepTime The time given to the thread to sleep between drawing each new move
	 * @param maxPieces Maximum permited score for this game
	 * @return Returns the score obtained by the AI
	 */
	public static int RunAI(double[] fitnessP, int sleepTime, int maxPieces, Boolean playingMode, Boolean runLookahead){
		State s = new State();
		TFrame frame = null;
		if(playingMode){
			frame = new TFrame(s);
		}

		PlayerSkeleton p = new PlayerSkeleton();
		int piecesDrawn = 0;
		while(!s.hasLost() && piecesDrawn < maxPieces) {
			if(runLookahead)
				s.makeMove(p.pickMoveWithLookahead(s,s.legalMoves(), fitnessP));
			else
				s.makeMove(p.pickMove(s,s.legalMoves(), fitnessP));
			if(playingMode){
				s.draw();
				s.drawNext(0,0);
				
				try {
					Thread.sleep(sleepTime);
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			
			piecesDrawn++;
		}
		if(playingMode) {
			frame.dispose();
		}
		return s.getRowsCleared();
	}
	/**
	 * Creates an initial population of randomly generated unit vectors in the 4D unit sphere for the 
	 * heuristic values
	 * @param populationCount The size of the desired population
	 * @param heuristicsCount Number of desired weight factors in each individual
	 * @return Returns the initial population
	 */
	public static double[][] getIntialPopulation(int populationCount, int heuristicsCount){
		double pEvolve[][] = new double[populationCount][heuristicsCount];
		//int fitnessP[] = new int[populationCount];
		
		for(int i = 0; i< populationCount; i++) {
			//fitnessP[i].first = 0;
			double magnitude = 0.0;

			for (int j = 0; j<heuristicsCount; j++) {
				pEvolve[i][j] = randnum.nextDouble()*20.0 - 10.0;
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
	
		int numGames = 5; // 5 games is enough, actually
		int populationCount = 200; // Biology says we need at least 160
		double finalParameters[] = new double[heuristicsCount];
		double pEvolve[][] = getIntialPopulation(populationCount, heuristicsCount);
		int averageScore = -1;
		int generationCount = 0;
		// When to stop?
		
		while(averageScore < Integer.MAX_VALUE - 1 && generationCount < 50) {
			generationCount++;
			System.out.printf("Generation: %d\n", generationCount);


			int maxScore = -1;
			System.out.println("Starting generation playing");

			//New parallelization method
			IntegerPair fitnessP[] = parallelizePerson(pEvolve, numGames, populationCount);
 
			int offspringCount = 0;
			int tenPercent = (int)(populationCount*0.1);
			int thirtyPercent = (int)(populationCount*0.3);
			double[][] offsprings = new double[populationCount][heuristicsCount];
			System.out.println("Selecting parents");
			// We want create offsprings up to 30% of the population
			while(offspringCount < thirtyPercent){

					// Select 10% of the population at random
					double[][] randomSelection = new double[tenPercent][heuristicsCount];
					int[] randomSelectionFitness = new int[tenPercent];
					int count = 0;
					for(int i = 0; i < populationCount; i += tenPercent) {
						int randIndex = i + randnum.nextInt(tenPercent);
						randomSelection[count] = pEvolve[randIndex];
						randomSelectionFitness[count++] = fitnessP[randIndex].first;
					}

					// Get the two most fittest from this random selection  
					IntegerPair maxFitness = new IntegerPair(-1, -1); // Value, Index
					IntegerPair maxFitnessSecond = new IntegerPair(-1, -1);

					for(int i = 0; i < tenPercent; i++){

						if(randomSelectionFitness[i] > maxFitness.first){
							maxFitnessSecond.first = maxFitness.first;
							maxFitnessSecond.second = maxFitness.second;
							maxFitness.first = randomSelectionFitness[i];
							maxFitness.second = i;

						} else if (randomSelectionFitness[i] > maxFitnessSecond.first){
							maxFitnessSecond.first = randomSelectionFitness[i];
							maxFitnessSecond.second = i;
						}
					}

					double[] p1 = randomSelection[maxFitness.second];
					double[] p2 = randomSelection[maxFitnessSecond.second];
					double[] offspring = new double[heuristicsCount];

					// Weighted average crossover
					for(int i = 0; i < heuristicsCount; i++){
						offspring[i] = ((randnum.nextDouble() > 0.5) ? p1[i]:p2[i]);
					}
		
					// Mutation
					if(randnum.nextDouble() < 0.10){
						int heuristicMutation = randnum.nextInt(heuristicsCount);
						offspring[heuristicMutation] += (offspring[heuristicMutation]*-0.1 + offspring[heuristicMutation]*0.2*randnum.nextDouble());
					}

					offsprings[offspringCount++] = offspring;

			}

			System.out.println("Replacing offspring");
			// NOW we need to get the weakest 30%, delete them and replace with new offsprings array
			// Then we can run the evolution again
			Arrays.sort(fitnessP);

			for (int i = 0; i < thirtyPercent; i++ ) {
				pEvolve[fitnessP[i].second] = offsprings[i];
			}
			finalParameters = pEvolve[fitnessP[populationCount -1].second];
			System.out.println("Playing game with bestVector so far");
			System.out.println("BestVector "+ finalParameters[0] + " " + finalParameters[1] + " " + finalParameters[2] + " " + finalParameters[3]
				+ " " + finalParameters[4] + " " + finalParameters[5]);
			averageScore = 0;
			for (int j = 0; j < 3; j++) {
				averageScore += RunAI(finalParameters, 20, Integer.MAX_VALUE, false, runEvoLookahead);
			}
			averageScore /= 3;
			System.out.printf("Average score of bestVector in generation: %d\n", averageScore);

		}

		// Run next evolution
		return finalParameters;
	}

	/**
	 * Main method for running either the game or the evolution
	 *
	 * @param args If an argument is passed(Whatever it is) it will run the evolution, otherwise it will run the game
	 */
	public static void main(String[] args) {
		int sumOfScore =0;
		randnum = new Random();
		// Our optimal weights
		double weightFactors[] = {-3.1472553592987946, 2.46883837144299,-2.2945510371937452,-7.521200744605782, -6.510648376902374,-2.1908554239402918}; // My best found values
		if(args.length > 0){
			
			if(args[0].equals("-e")){
				System.out.println("Starting Evolution");
				weightFactors = RunEvolution(6);
				// Returns empty, not yet finished
				System.out.println("Completed evolution with: " + weightFactors[0] + " " + weightFactors[1] + " " + weightFactors[2] + " " + weightFactors[3]  + " " + weightFactors[4]);
			} else if (args[0].equals("-l")){
				System.out.println("Running AI with lookahead");
				System.out.println("You have completed "+ RunAI(weightFactors, 30, Integer.MAX_VALUE, false, true) +" rows with lookahead.");
				return;
			} else if(args[0].equals("-el")){
				runEvoLookahead = true;
				System.out.println("Starting Evolution with lookahead");
				weightFactors = RunEvolution(6);
				// Returns empty, not yet finished
				System.out.println("Completed lookahead evolution with: " + weightFactors[0] + " " + weightFactors[1] + " " + weightFactors[2] + " " + weightFactors[3]  + " " + weightFactors[4]);
				System.out.println("You have completed "+ RunAI(weightFactors, 30, Integer.MAX_VALUE, false, true) +" rows with lookahead.");
				return;
			} 
		}

		System.out.println("You have completed "+ RunAI(weightFactors, 30, Integer.MAX_VALUE, false, false) +" rows.");
	}
	
}

class IntegerPair implements Comparable<IntegerPair>{

	public int first = 0;
	public int second = 0;

	public IntegerPair(int first, int second){
		this.first = first;
		this.second = second;
	}
	public IntegerPair()
    {

    }

    @Override
    public int compareTo(IntegerPair o) {
        return first < o.first ? -1 : first > o.first ? 1 : 0;
    }

}



class IntegerDouble implements Comparable<IntegerDouble>{

	public double first = 0;
	public double second = 0;

	public IntegerDouble(double first, double second){
		this.first = first;
		this.second = second;
	}
	public IntegerDouble()
    {

    }
    @Override
    public int compareTo(IntegerDouble o) {
       
        return first < o.first ? -1 : first > o.first ? 1 : 0;

    }

}


