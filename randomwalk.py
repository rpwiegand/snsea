## This is a silly little script to help demonstrate that random walks in d>2 stay bounded 
## relatively "close" to the starting point.

import numpy as np
import pandas as pd
import sys


def walk(d=5, steps=100, sigma=1):
    """
    Take `step` number of steps, starting at the origin, where each step is a random
    normal deviate from the previous step, with standard deviation `sigma`.  The dimensionality
    of the space is specified using the `d` parameter.  Return the distance of the final point from the origin.
    """
    X = np.random.normal(size=d*steps, scale=sigma).reshape((d,steps))
    Y = np.cumsum(X, axis=0)
    return np.sum(np.square(Y), axis=0)



def main(numTrials=30, d=5, steps=100, sigma=1):
    """
    Take `numTrials` random walks and return the distances obtained.
    """
    df = pd.DataFrame()

    for trial in range(numTrials):
        trialLabel = "Trial" + str(trial)
        df[trialLabel] = walk(d, steps, sigma)

    return df



if __name__ == '__main__':
    numTrials = 3
    d = 5
    steps = 100
    sigma = 1

    try:
        numTrials = int(sys.argv[1])
        d = int(sys.argv[2])
        steps = int(sys.argv[3])
        sigma = float(sys.argv[4])
    except:
        pass


    df = main(numTrials, d, steps, sigma)
    outputName = "randomwalk-" + str(d) + "-" + str(steps) + "-" + str(sigma) + ".csv"
    print("Writing output to: ", outputName)
    df.to_csv(outputName)


