##TetrisAI
This is a utility-based agent that for every move evaluates the boards state using a utility function based on a linear weighted sum of heuristics.

<p><img src="/img/tetris_res.gif" width = "128"/></p>

### Heuristics 
 The boards state consists of six heuristics which are: 

* Number of holes
* Landing height
* Row Transitions
* Column transitions
* Well Sums
* Completed Lines 

The optimal parameters were obtained by running a genetic algorithm on multiple threads. The optimal values can be found in PlayerSkeleton.java

###To Run:
- Compile the project using ```javac```
- Run a new evolution with ```java PlayerSkeleton.java -e``` **or** run with our best found parameters with ```java PlayerSkeleton.java```

The GUI is currently disabled to optimize the calculations, to activate the GUI see comments in code where to toggle some boolean values 

###Files
* State.java - tetris simulation
* TFrame.java - frame that draws the board
* TLabel.java - drawing library
* PlayerSkeleton.java - setup for implementing a player
	
	
#####State.java
This is the tetris simulation.  It keeps track of the state and allows you to 
make moves.  The board state is stored in field (a double array of integers) and
is accessed by getField().  Zeros denote an empty square.  Other values denote
the turn on which that square was placed.  NextPiece (accessed by getNextPiece)
contains the ID (0-6) of the piece you are about to play.

Moves are defined by two numbers: the SLOT, the leftmost column of the piece and
the ORIENT, the orientation of the piece.  Legalmoves gives an nx2 int array
containing the n legal moves.  A move can be made by specifying the two
parameters as either 2 ints, an int array of length 2, or a single int
specifying the row in the legalMoves array corresponding to the appropriate move.

It also keeps track of the number of lines cleared - accessed by getRowsCleared().

draw() draws the board.
drawNext() draws the next piece above the board
clearNext() clears the drawing of the next piece so it can be drawn in a different
	slot/orientation

#####TFrame.java
This extends JFrame and is instantiated to draw a state.
It can save the current drawing to a .png file.
The main function allows you to play a game manually using the arrow keys.

#####TLabel.java
This is a drawing library.

#####PlayerSkeleton.java
An example of how to implement a player.
The main function plays a game automatically (with visualization).
