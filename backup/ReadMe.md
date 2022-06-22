# CSCI 208, Project 3, Fall 2021

Compete the program provided that builds a sorted linked list of movies, sorted by rating.  Your program must
* Query the user as to a file name
* Read in that file containing a list of movie titles and their ratings (format below)
* Create a linked list, inserting the movies in the correct order into the list
* Output that list to the screen in sorted order when you are complete.


## Movie File Format:

Each line of the file will contain information about a movie:  the average rating followed by the name of the movie. The rating will be a double, and the movie name a string that *can contain spaces*.  There are several ways to handle this, the easiest of which is to use the extraction operator to get the double then *getline()* to get the rest of the line.  Don't forget to check edge cases (empty lines, the last line, etc.).  Below is an example.  I took these from IMDB's "Top Rated Movies", so please do not infer anything about my movie tastes from it.

```
8.9 The Dark Knight
9.2 The Shawshank Redemption
9.0 The Gadfather: Part II
9.1 The Godfather
```

## Output Format

The output should use the same style of format, but be in *descending order* by rating, like this:

```
9.2 The Shawshank Redemption
9.1 The Godfather
9.0 The Gadfather: Part II
8.9 The Dark Knight
```



## How and What to Submit
Submission will happen via this GitHub respository.  So you'll need to make sure you *add* all necessary files , and you will also need to make sure your code is committed and pushed to GitHub *before* the submission due date.  There's no need to do anything in BlackBoard other than accept the initial invitation link.

Here are the files that should be present in your repository:
* Any necessary C++ source files
* A *Makefile* that builds your program
* At least one example input files
* Any special instructions for running your program (if needed).


## Hints For Completing This

* You have now seen many examples of a linked list, so look over your examples.
* The hardest part of this assignment is inserting *in order*, which requires traversing to the correct position and inserting a node there.
* Start simple and gradually add complexity (see notes below).
* **SKETCH THE INSERTION PROCESS** before you write your code.  
* Keep in mind that there are several cases -- think them through:
    1.  What is the list is empty when you are ready to insert?  
    2.  What if there is a list but the movie you are inserting must go in front?
    3.  Same question, but at the end?
    4.  What if the movie goes somewhere in the middle?


## Possible Stages of Development

This is the hardest programming assignment.  Start early.  Don't try to build the whole program at once.  Here is one plausible set of stages:

1. Start by just reading in the file, parsing the score and the name, and outputting those fields without a linked list.
2. Once that words, make your program build the link list by inserting to the front.  Use labs and my example code to help you.
3. Once you have insertion at the front working, write a routine that traverses that list and outputs it in the correct format.
4. Once you can read the file, insert into the front of the link list, and output in the correct format, you are almost there!  Modify your program to perform the insertion in order instead of always at the front.


