
import java.util.*;
import java.io.*;
import java.util.Random;


//Why did I create a new class? Needs to implement runnable to be able to call the threads
//It was not working when I did this within class PlayerSkeleton
//Going to optimize (basically remove main and just have PlayerSkeleton implements runnable)
//Main will thus have to summon parallel threads, but atm easier for me to keep this code seperated
/*public class gameThread extends PlayerSkeleton implements runnable {
	int threadNum;
	int numGames;
	int populationCount;
	double[][] pEvolve;
	int[] fitnessP;

	public void run() {
		for(int i = 0; i< populationCount; i++){
			System.out.println("Evolving "+ pEvolve[i][0] + " " + pEvolve[i][1] + " " + pEvolve[i][2] + " " + pEvolve[i][3]);
			for (int j = 0; j<numGames; j++) {
				int rowsCleared = RunAI(pEvolve[i], 30, 200, evolutionMode);
				System.out.println("Game " + (j+1)+ ": "+ rowsCleared);
				fitnessP[i].first += rowsCleared;
			}
			System.out.println("Total " + fitnessP[i].first);
		}
	}

	public gameThread(int i, int games, int popCount, double probEvolve, int fitnessProb) {
		threadNum = i;
		numGames = games;
		populationCount = popCount;
		pEvolve = probEvolve;
		fitnessP = fitnessProb;
	}
}*/

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


		double bumpinessNew = 0;
		
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
			for (int i = newTop[j]-1; i >= 0; i--) {
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

	public double bruteForceFitness(final int[][] newField){

		double a = -0.510066, b = 0.760666, c = -0.35663, d = -0.184473;
		int score = 0;

		// calc height
		int height = 0;
		int maxRow = newField.length;
		int maxCol = newField[0].length;
		
		// Height
		for(int col = maxCol - 1; col >= 0; col--){
			for(int row = maxRow - 1; row >=0; row--){
				if(newField[row][col] != 0){
					//System.out.println("row, col "+ row +" "+ col);
					height += (row + 1);
					break;
				}
			}
		}

		// Lines cleared 
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


		// Holes
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

		// Bumbiness
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

		return height*a + linesCleared*b + holes*c + bumpiness*d;
		
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
	 * @param weightsCount Number of desired weight factors in each individual
	 * @return Returns the initial population
	 */
	public static double[][] getIntialPopulation(int populationCount, int weightsCount){
		double pEvolve[][] = new double[populationCount][weightsCount];
		//int fitnessP[] = new int[populationCount];
		
		for(int i = 0; i< populationCount; i++) {
			//fitnessP[i].first = 0;
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
		int numGames= 100;
		int populationCount = 1000;
		double finalParameters[] = new double[heuristicsCount];
		double pEvolve[][] = getIntialPopulation(populationCount, heuristicsCount);
		IntegerPair fitnessP[] = new IntegerPair[populationCount];
		int highestScore = -1;
		int generationCount = 0;
		// When to stop?
		while(highestScore < Integer.MAX_VALUE) {
			generationCount++;
			System.out.printf("Generation: %d", generationCount);
			//int numThreads = 4;
			//ArrayList<Thread> gameThreads = new ArrayList<Thread>();

			//Create an array containing as many threads as we can efficiently run on our system
			//Having issues passing pEvolve and fitnessP
			/*for (int i = 0; i < numThreads; i++) {
				gameThreads.add(new Thread(new gameThread(i, numGames, populationCount, pEvolve, fitnessP)));
			}
			//Start the threads
			for (int i = 0; i < gameThreads.size(); i++) {
				gameThreads.get(i).start();
			}*/

			// End parallel

			//YET TO IMPLEMENT: Adding results to an array and spitting it back or writing to file
			//so that we can run evolution processing using the parallel outputs

			//This is now computed within the thread (not complete)
			int gameNumber = 0;
			int maxScore = -1;
			
			for(int i = 0; i < populationCount; i++){
				fitnessP[i] = new IntegerPair(0, i);

				// System.out.println("Playing subject: "+ i + "\n" + pEvolve[i][0] + " " + pEvolve[i][1] + " " + pEvolve[i][2] + " " + pEvolve[i][3]);
				for (int j = 0; j<numGames; j++) {
					int rowsCleared = RunAI(pEvolve[i], 30, 700, false);
					// System.out.println("Game " + (j+1)+ ": "+ rowsCleared);
					// System.out.printf("%d", fitnessP.first);
					// breakpoint;
					fitnessP[i].first += rowsCleared;
					
					gameNumber++;
				}
				if(fitnessP[i].first > maxScore) {
					maxScore = fitnessP[i].first;
					finalParameters = pEvolve[i];
				}

			}
			System.out.printf("Game: %d TotalScore %d\n", gameNumber, maxScore);
			System.out.println("BestVector "+ finalParameters[0] + " " + finalParameters[1] + " " + finalParameters[2] + " " + finalParameters[3]);
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

					int f_p1 = maxFitness.first;
					int f_p2 = maxFitnessSecond.first;
					double[] p1 = randomSelection[maxFitness.second];
					double[] p2 = randomSelection[maxFitnessSecond.second];
					double[] offspring = new double[heuristicsCount];

					// System.out.println("Elite parent1: " + p1[0] + " " + p1[1] + " " + p1[2] + " " + p1[3] + " with fitness " + f_p1);
					// System.out.println("Elite parent2: " + p2[0] + " " + p2[1] + " " + p2[2] + " " + p1[3] + " with fitness " + f_p2);

					// Weighted average crossover
					for(int i = 0; i < heuristicsCount; i++)
						offspring[i] = p1[i] * f_p1 + p2[i] * f_p2;
					
					// Mutation
					if(randnum.nextDouble() < 0.05){
						offspring[randnum.nextInt(4)] += (-2 + 4 * randnum.nextDouble());
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
			highestScore = RunAI(pEvolve[fitnessP[populationCount -1].second], 20, Integer.MAX_VALUE, true);
			// if(rowsCleared > highestScore) { 
			finalParameters = pEvolve[fitnessP[populationCount -1].second];
			System.out.println("BestVector "+ finalParameters[0] + " " + finalParameters[1] + " " + finalParameters[2] + " " + finalParameters[3]);
			// }
			System.out.printf("TotalScore of best new population: %d\n", highestScore);
			

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
		randnum = new Random();
		double weightFactors[] = {-0.510066, 0.760666, -0.35663, -0.184483};
		
		if(args.length > 0){
			weightFactors = RunEvolution(4);
			// Returns empty, not yet finished
			System.out.println("Completed evolution with: " + weightFactors[0] + " " + weightFactors[1] + " " + weightFactors[2] + " " + weightFactors[3]);
		}
		System.out.println("You have completed "+ RunAI(weightFactors, 20, Integer.MAX_VALUE, true)+" rows.");
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

