# CSCI 208, Spring 2021, Programming Project 4

## Write a Robust Bowling Scorer
In this assignment, you will start your code *almost* from scratch, developing an application that reads in a file of bowling results and reports the scores for all the players.  Your code should deal with *all* potential runtime errors, and in addition it *must* make use of exceptions for at least one such error.  That is:  You must deal with *all* errors in some way (with or without exceptions), but you have to use exceptions at least once.  You *can* create your own exceptions, but you do not have to (you can just handle existing exceptions, if you like).

Your program should prompt the user for the input file, then read the input file, estimate the scores for all players represented in the file, and output those scores.  I have written a *brittle* routine to score (and potentially print) games for you.  It does not handle errors very well, so you'll need to manage that part.

## Input File Format
The file will represent a single game from several players.  Each player's game will be on a line by itself.  The line contains the player's name, followed by a list of the number of pins they got for each ball they threw, all separated by spaces.  For example:

```    1     2   3    5    5    6     7   8    9    10
Paul  9 1  10   5 1  7 1  8 2  7 2  10   5 3  9 0  8 2 7
```   

This is a spare in the first frame, a strike in the second frame, etc.  The game above is worth 128 points.  The frame-byframe count looks like this:

*  20, 36, 42, 50, 67, 76, 94, 102, 111, 128  

## How Bowling Is Scored
A game of bowling is comprised of 10 frames, and in each frame there are 10 pins the knock down with a ball (with one exception).  The player gets up to two attempts each frame to knock pins down.  There are three things that can happen in most frames that affect scoring:

1.  Some pins are knocked down -- the player's score is increased by the total number of pins knocked down.  For instance, if the player knocks down 5 pins, then a remaining 2 pins, a 7 is added to her or his score.
2.  The player needs two balls to knock all the pins down.  This is called a *spare*.  In this case, the player adds 10 points to her or his score, *plus* however many pins she or he knocks down in the first roll of the *next* frame.
3.  The player requires only *one* ball to know all the pins down.  This is called a *strike*.  In this case, the player adds 10 points to her or his score, *plus* however many pins she or he knocks down in her or his next *two* rolls.

The one exception is the 10th frame.  If the player gets a spare, she or he gets to throw one more ball.  If the player gets a strike, she or he gets to throw two more balls.  

The highest possible game is obtained by throwing 12 strikes in a row.  This is worth 300 points.

## Be Careful About Assumptions
Note, because of strikes and the 10th frame rule, not every game results in the same number of throws.  A player that doesn't get any strikes and does not spare in the 10th frame will roll 20 balls.  A player who gets *all* strikes will roll only 12.  Likewise, not every frame has the same number of balls thrown.

Spares and strikes have a *compounding* affect on the score, because those pins are actually counted more than once.  A spare in the first frame and a 5 on the first throw of the second, gets 15 points in the first frame (10+5), but those five pins in the second frame *also* count in the second frame.




## Using GitHub Classroom
You have obtained this code by accepting my GitHub Classroom invitation link.  I will grade the materials I find in your repository on GitHub.  So you are responsible for making sure you have committed and pushed your code to the source repository *before* the due date.
