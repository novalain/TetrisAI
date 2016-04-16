
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


		// int lowestTop = -1;
		
			// for(int i = 0; i<maxCol; i++) {
			// 	System.out.printf("%d ", newTop[i]);
			// }
			// System.out.println();
		for(int i = 0; i<maxCol; i++) {
			for (int j  = newTop[i]-1; j >=0; j--) {
				if(newField[j][i] == 0) numHoles++;
			}
			// System.out.println(Math.max(newTop[i]-1, 0));
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
				// if(!currentCell && newField[j+1][i] !=0) numHoles++;
				if(lastCell != currentCell) {
					columnTransitions++;
				}
				lastCell = currentCell;
			}
			// if(!currentCell) columnTransitions++;
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
		
		// System.out.println(landingHeight); // Done
		// System.out.println(rowsCleared); // Done
		// System.out.println(rowTransitions); // Done
		// System.out.println(columnTransitions); // Done
		// System.out.println(numHoles); // Done
		// System.out.println(wellSums);

		return landingHeight*weights[0] + rowsCleared*weights[1] + rowTransitions*weights[2] + 
		columnTransitions*weights[3] + numHoles*weights[4] + wellSums*weights[5];
		
	}

	// public double fitnessFunction(final int[][] newField, final int[] newTop, double[] weights){

	// 	// Our magic numbers :)
	// 	int maxRow = newField.length;
	// 	int maxCol = newField[0].length;

		
	// 	double completeLines = 0; // Done
	// 	double numHoles = 0; // Done // Squared
	// 	// double numHolesSq = 0; // Squared // Done
	// 	double maxHoleHeight = 0; // Squared // Done
	// 	// double maxHoleHeightSq = 0; // Done
	// 	double maxColumnHeight = 0; //Squares // Done
	// 	// double maxColumnHeightSq = 0; // Done
	// 	double columnWithHoles = 0; // Done
	// 	double rowsWithHoles = 0; // Done
	// 	double totalHeight = 0; // Done
	// 	double lowestPlayableRow = 0; // Don't know what this is //Have no idea how to do
	// 	double bumpiness = 0; // Roughness in the paper // Done
	// 	double maxPitDepth = 0; // Squared // Have no idea how to do
	// 	// double maxPitDepthSq = 0;
	// 	double slope = 0; // Done
	// 	double convexity = 0; // Have no idea how to do
		

	// 	double bumpinessNew = 0;
		
	// 	int lowestTop = 100;

	// 	// Calculate Height and Bumpiness
	// 	for(int col = 0; col<maxCol; col++){
	// 		if(newTop[col] > maxColumnHeight) maxColumnHeight = newTop[col] -1;
	// 		totalHeight+= newTop[col];
	// 		if(newTop[col] < lowestTop) lowestTop = newTop[col];
	// 		if(col!= 0) {
	// 			bumpiness += Math.abs(newTop[col]-newTop[col-1]);
	// 			slope += newTop[col]-newTop[col-1];
	// 		}	
	// 	}
	// 	slope = Math.abs(slope);

	// 	boolean[] rowIsHoleFree = new boolean[maxRow];

	// 	// Calculate totalHeight and Bumpiness
	// 	for (int j = 0; j < maxCol; j++){
	// 		Boolean columnHoleFree = false;
	// 		// Boolean representation of if we have found the columns top yet
	// 		// Do we need to check top row? If we're there we lost right?
	// 		for (int i = newTop[j]-1; i >= 0; i--) {
	// 			// If the top filled square of the column is not yet found and we found a filled square 
	// 			if(newField[i][j] == 0) {

	// 				if(i > maxHoleHeight) maxHoleHeight = i;
	// 				if(columnHoleFree) {
	// 					columnWithHoles++;
	// 					columnHoleFree = true;
	// 				}

	// 				if(rowIsHoleFree[i]) {
	// 					rowsWithHoles++;
	// 					rowIsHoleFree[i] = true;
	// 				}
	// 				//It must be a hole, so increment
	// 				numHoles += 1;

	// 			}
	// 		}
	// 	}
		
	// 	for(int i = 0; i < lowestTop; i++ ){
	// 		int isCompleted = 1;
	// 		for(int j = 0; j < State.COLS; j++){
	// 			if(newField[i][j] == 0){
	// 				isCompleted = 0;
	// 				break;
	// 			}
	// 		}
	// 		completeLines += isCompleted;
	// 	}
	// 	// System.out.println(completeLines); // Done
	// 	// System.out.println(numHoles); // Done // Squared
	// 	// System.out.println(numHolesSq); // Squared // Done
	// 	// System.out.println(maxHoleHeight); // Squared // Done
	// 	// System.out.println(maxHoleHeightSq); // Done
	// 	// System.out.println(maxColumnHeight); //Squares // Done
	// 	// System.out.println(maxColumnHeightSq); // Done
	// 	// System.out.println(columnWithHoles); // Done
	// 	// System.out.println(rowsWithHoles); // Done
	// 	// System.out.println(totalHeight); // Done
	// 	// System.out.println(lowestPlayableRow); // Don't know what this is //Have no idea how to do
	// 	// System.out.println(bumpiness); // Roughness in the paper // Done
	// 	// System.out.println(maxPitDepth); // Squared // Have no idea how to do
	// 	// System.out.println(maxPitDepthSq);
	// 	// System.out.println(slope); // Done
	// 	// System.out.println(convexity); // Have no idea how to do
	// 	double score = completeLines*weights[0] + numHoles*weights[1] + Math.pow(numHoles, 2)*weights[2] + maxHoleHeight*weights[3] +
	// 	Math.pow(maxHoleHeight, 2)*weights[4] + maxColumnHeight*weights[5] + Math.pow(maxColumnHeight, 2)*weights[6] + columnWithHoles*weights[7] +
	// 	rowsWithHoles*weights[8] + totalHeight*weights[9] + lowestPlayableRow*weights[10] + bumpiness*weights[11] + maxPitDepth*weights[12] + 
	// 	Math.pow(maxPitDepth,2)*weights[13] + slope*weights[14] + convexity*weights[15];
	// 	// System.out.println(score);
	// 	return score;
	// 	//return totalHeight*weights[0] + completeLines*weights[1] + numHoles*weights[2] + bumpiness*weights[3];
	
	// }

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

		//TreeMap<Double, Integer> scores = new TreeMap<Double, Integer>();
		//Map<Double, List<Integer>> scores = new Map<Double, int[]>(); 
		ArrayList<Integer> possibleMoves = new ArrayList<Integer>();

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
		// System.out.println(legalMoves.length);
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
				// for(int k = 0; k  < newField.length; k ++) { 
				// 	for(int j = 0; j<newField[0].length; j++) {
				// 		System.out.printf("%03d  ", newField[k][j]);
				// 	}
				// 	System.out.printf("\n");
				// }
				// System.out.printf("\n");
				double score = newFitnessFunction(newField, newTop, weightFactors);
				//assert score != brute_force_score;

				/*if(Math.abs(score - brute_force_score) > 0.000000001)
					System.out.println("** DIFF IN BRUTE AND FITNESS");*/
					
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


	//Start of Parallelization code
	public static class gameThread implements Callable<Integer> {
		private final double[][] pEvolve;
		private final int guy;
		private int rowsCleared;

		public gameThread(double[][] probEvolve, int i) {
			this.pEvolve = probEvolve;
			this.guy = i;
		}

		public Integer call() {
			//System.out.println("Thread " + guy + " has BEGUN");
			try {
				rowsCleared = RunAI(pEvolve[guy], 300, Integer.MAX_VALUE, false);
			} catch (Exception e) {
				e.printStackTrace();
			}
		//System.out.println("Thread " + guy + " has FINISHED");
		//System.out.println("Guy " + guy + " cleared: " + rowsCleared + " rows");
		return rowsCleared;
		}
	}

	public static int parallelizePerson(double[][] pEvolve, int guy, int numGames){

	    int threads = Runtime.getRuntime().availableProcessors();
	    //System.out.println("We have: " + threads + " threads to play on");
	    ExecutorService executor = Executors.newFixedThreadPool(threads);
	    ArrayList<Future<Integer>> futures = new ArrayList<Future<Integer>>();
	    // Now that we don't play that many games, we don't use the whole CPU. We were creating 100 background threads
	    // Which is not efficient due to communication time.
	    for (int i = 0; i < numGames; i++) {
	    	Callable<Integer> callable = new gameThread(pEvolve, guy);
		    Future<Integer> future = executor.submit(callable);
		    futures.add(future);
	    }

	    executor.shutdown();

	    int totalRows = 0;
	    for (Future<Integer> future : futures) {
	    	try {
	    		totalRows += future.get();
	    		//System.out.println("TotalRows = " + totalRows);
	    	} catch (InterruptedException e) {
	    		e.printStackTrace();
	    	} catch (ExecutionException e) {
	    		System.out.println(e);
	    	}
	    }
	    //System.out.println("Guy " + guy + " cleared: " + totalRows + " in total");
	    return totalRows;
	}
	//End of Parallelization code


	/**
	 * @param fitnessP The weight factors passed into the fitness function
	 * @param sleepTime The time given to the thread to sleep between drawing each new move
	 * @param maxPieces Maximum permited score for this game
	 * @return Returns the score obtained by the AI
	 */
	public static int RunAI(double[] fitnessP, int sleepTime, int maxPieces, Boolean playingMode){
		State s = new State();
		TFrame frame = null;
		if(playingMode){
			frame = new TFrame(s);
		}

		PlayerSkeleton p = new PlayerSkeleton();
		int piecesDrawn = 0;
		while(!s.hasLost() && piecesDrawn < maxPieces) {
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
				// magnitude += Math.pow(pEvolve[i][j], 2);
			}
			// magnitude = Math.sqrt(magnitude);
			// for (int j = 0; j<heuristicsCount; j++) {
			// 	pEvolve[i][j] /= magnitude;
			// }
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
		int numGames = 5; // 5 games is enough, actually
		int populationCount = 200; // Biology says we need at least 160
		double finalParameters[] = new double[heuristicsCount];
		double pEvolve[][] = getIntialPopulation(populationCount, heuristicsCount);
		IntegerPair fitnessP[] = new IntegerPair[populationCount];
		int averageScore = -1;
		int generationCount = 0;
		// When to stop?
		
		while(averageScore < Integer.MAX_VALUE - 1 && generationCount < 50) {
			generationCount++;
			System.out.printf("Generation: %d\n", generationCount);


			int maxScore = -1;
			System.out.println("Starting generation playing");
			for(int i = 0; i < populationCount; i++){
				fitnessP[i] = new IntegerPair(0, i);

				// System.out.println("Playing subject: "+ i + "\n" + pEvolve[i][0] + " " + pEvolve[i][1] + " " + pEvolve[i][2] + " " + pEvolve[i][3]);
				// for (int j = 0; j<numGames; j++) {
				// 	//Thread singleGame = new gameThread(pEvolve, i);
				// 	int rowsCleared = RunAI(pEvolve[i], 30, 700, false);
				// 	// System.out.println("Game " + (j+1)+ ": "+ rowsCleared);
				// 	// System.out.printf("%d", fitnessP.first);
				// 	// breakpoint;
				// 	fitnessP[i].first += rowsCleared;

				// }
				int totalCleared = parallelizePerson(pEvolve, i, numGames);
				fitnessP[i].first += totalCleared;

				// if(fitnessP[i].first > maxScore) {
				// 	maxScore = fitnessP[i].first;
				// 	finalParameters = pEvolve[i];
				// }

			}
			// System.out.printf("Max TotalScore %d\n", maxScore);
			// System.out.println("BestVector of current generation");
			// System.out.println(finalParameters[0] + " " + finalParameters[1] + " " + finalParameters[2] + " " + finalParameters[3]);
			// Select parents and produce offsprings part 
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

						} else if (randomSelectionFitness[i] > maxFitnessSecond.first && randomSelectionFitness[i] < maxFitness.first){
							maxFitnessSecond.first = randomSelectionFitness[i];
							maxFitnessSecond.second = i;
						}
					}

					// int f_p1 = maxFitness.first;
					// int f_p2 = maxFitnessSecond.first;
					double[] p1 = randomSelection[maxFitness.second];
					double[] p2 = randomSelection[maxFitnessSecond.second];
					double[] offspring = new double[heuristicsCount];

					// System.out.println("Elite parent1: " + p1[0] + " " + p1[1] + " " + p1[2] + " " + p1[3] + " with fitness " + f_p1);
					// System.out.println("Elite parent2: " + p2[0] + " " + p2[1] + " " + p2[2] + " " + p1[3] + " with fitness " + f_p2);

					// double[] test = new double[heuristicsCount];
					// Weighted average crossover
					for(int i = 0; i < heuristicsCount; i++){
						offspring[i] = ((randnum.nextDouble() > 0.5) ? p1[i]:p2[i]);
						// offspring[i] = p1[i] * f_p1 + p2[i] * f_p2;
						// test[i] = p1[i] * f_p1 + p2[i] * f_p2;
					}
		
					// Mutation
					if(randnum.nextDouble() < 0.10){
						int heuristicMutation = randnum.nextInt(heuristicsCount);
						offspring[heuristicMutation] += (offspring[heuristicMutation]*-0.1 + offspring[heuristicMutation]*0.2*randnum.nextDouble());
					}

					// // Normalize 
					// double magnitude = 0;
					// for (int i = 0; i<heuristicsCount; i++) {
					// 	magnitude += Math.pow(offspring[i], 2);
					// }
					// magnitude = Math.sqrt(magnitude);
					// for (int i = 0; i<heuristicsCount; i++) {
					// 	offspring[i] /= magnitude;
					// }

					offsprings[offspringCount++] = offspring;

			}

			System.out.println("Replacing offspring");
			// NOW we need to get the weakest 30%, delete them and replace with new offsprings array
			// Then we can run the evolution again
			// System.out.println("Sorting");
			Arrays.sort(fitnessP);
			// System.out.println("Sorted");
			// for (int i = 0; i < populationCount; i++) {
			// 	if(fitnessP[i].second == 0)
			// 		System.out.println(fitnessP[i].first + " " + fitnessP[i].second);
			// }
			// 
			for (int i = 0; i < thirtyPercent; i++ ) {
				// System.out.println("Replacing: " + fitnessP[i].first);
				pEvolve[fitnessP[i].second] = offsprings[i];
				// System.out.println("Replacing");
			}
			finalParameters = pEvolve[fitnessP[populationCount -1].second];
			System.out.println("Playing game with bestVector so far");
			System.out.println("BestVector "+ finalParameters[0] + " " + finalParameters[1] + " " + finalParameters[2] + " " + finalParameters[3]
				+ " " + finalParameters[4] + " " + finalParameters[5]);
			averageScore = 0;
			for (int j = 0; j < 3; j++) {
				averageScore += RunAI(finalParameters, 20, Integer.MAX_VALUE, false);
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
		// double weightFactors[] = {-0.510066, 0.760666, -0.35663, -0.184483};// 0};
		// Dude's weights
		// double weightFactors[] =  {-4.500158825082766, 3.4181268101392694, -3.2178882868487753, -9.348695305445199, -7.899265427351652, -3.3855972247263626};
		double weightFactors[] = {-3.1472553592987946, 2.46883837144299,-2.2945510371937452,-7.521200744605782, -6.510648376902374,-2.1908554239402918}; // My best found values
		if(args.length > 0){
			weightFactors = RunEvolution(6);
			// Returns empty, not yet finished
			System.out.println("Completed evolution with: " + weightFactors[0] + " " + weightFactors[1] + " " + weightFactors[2] + " " + weightFactors[3]  + " " + weightFactors[4]);
		}

		System.out.println("You have completed "+ RunAI(weightFactors, 30, Integer.MAX_VALUE, false) +" rows.");
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

