import sys
import numpy as np
import configReader
import scipy.stats as stats
import os, shutil

def initialize(n, allzero=False, sigma=0.0):
  """
  Initialize the single-population EA.  If the allzero flag is on, then 
  the individual is initialized at the 0^n position.  If sigma is 0, then
  it is assumed the individual is a binary string.  Otherwise, it is assumed
  the individual is a real valued vector.  When not all zeros, binary strings
  are initialized uniformly at random from {0,1}^n, and real vectors are
  initialized from [0,1]^n uniformly at random.
  """
  # First, setup an all-zero vector
  x = np.array([0]*n)

  # If not allzero and we're binary, init in {0,1}^n i.i.d.
  if (not allzero and sigma <= 0.0):
    x = np.random.random_integers(low=0, high=1, size=n)

  # If not allzero and we're numeric, init in [0,1]^n
  elif (not allzero and sigma > 0.0):
    x = np.random.uniform(low=0.0, high=1.0, size=n)

  # Return the initialized vector
  return (x)


def isInsideEuclideanBoundary(x, low=0, high=1):
  """
  Determine whether the selected point is inside or outside
  the square boundary given.  By default, this boundary is
  the unit square.
  """
  inside = True
  for xitem in x:
    if (xitem < low) or (xitem > high):
      inside =False
  return (inside)


def mutate(x, pm, sigma=0.0, bound=[0,1]):
  """
  Mutate the individual, where pm is the probability of mutating for 
  binary representation (ignored for numeric), and sigma is the Gaussian
  spread for each real-valued gene.  If sigma is 0 or lower, we assume
  a binary representation.  If the bound variable is set, repeat 
  mutations until child is inside the bound.  To turn this off,
  set bound=None.
  """
  # Get the length and initialize a child vector
  n = len(x)
  child = np.array([0]*n)

  # If we're using a binary representation, then flip bits
  # at random according to pm, i.i.d.
  if (sigma <= 0.0):
    mask = np.random.choice([0, 1], size=n, p=[1-pm,pm])
    y = x.copy()
    child = y ^ mask

  # Otherwise, jitter the gene by sigma according to
  # a normal distribution, i.i.d. for each gene.  If there
  # is a bound, then repeat mutation until a valid child
  # is produced.
  else:
    inside = False
    stuckCount = 0
    while (not inside):
      # Check if we're stuck
      if stuckCount > 100:
        raise Exception("The re-mutate loop is stuck ... should not take this many iterations to find a good point.")
      stuckCount += 1

      # Generate a mutation
      offsets = np.random.normal(scale=sigma, size=n)
      child = x + offsets

      # Check if the mutation is inside the bounds
      if (bound==None):
        inside=True
      else:
        inside = isInsideEuclideanBoundary(child, bound[0], bound[1])

  # Return the child
  return (child)


def getL2Norm(x1, x2):
  """
  Return the Euclidean distance betweeen two points.
  """
  return (np.sqrt( sum( (x1-x2)**2 ) ))


def computeSparseness(x, archive, k):
  """
  Give a potential candidate for the archive, the archive, and k, 
  estimate the sparseness contribution of x.  This is the
  average distances over the k closest neighbors in the archive.
  """
  # Compute the distance to x from every point in the archive
  # and put them in a vector.
  distances = []
  for xi in archive:
    distances.append( getL2Norm(x,xi) )

  # Arrange things so that we can get either
  # the k closest values or the size of the 
  # sarchive, whichever is smaller
  numElements = min(len(distances), k)
  distances.sort()

  # Estimate sparseness as the average over that subset, then return
  sparseness = float(sum(distances[0:(numElements+1)]))/float(numElements)
  return (sparseness)


def getArchiveDistances(archive, k, n):
  """
  Get minimum archive distances and minimum distance to all 1 string.
  """
  # Setup the arbitrary "piton" measure (the all 1 binary string)
  minDistTo1N = n
  OneN = np.array([1.0]*n)

  # Initialize the distances
  distances = []

  # Compute pair-wise distances of all points in the archive.
  # Also, keep track of the distance to the all 1 string, 1^n.
  for idx in range(len(archive)):
    minDistTo1N = min( minDistTo1N, getL2Norm(archive[idx], OneN) )
    tmpDistances = []
    for jdx in range(len(archive)):
      if (not idx == jdx):
        xi = archive[idx]
        xj = archive[jdx]
        tmpDistances.append( getL2Norm(xi,xj) )

    # Use these pairwise distances to estimate the sparseness of
    # each point in the archive to the rest of the archive
    # (excluding itself).
    numElements = min(len(tmpDistances), k)
    
    if numElements > 0:
      tmpDistances.sort()
      dist = float(sum(tmpDistances[0:(numElements+1)]))/float(numElements)
      distances.append( dist)

  # Return the sparseness distances for all points in the archive, as well
  # as the minimum distances of any point to the "piton" point at the all 1
  # string.
  return (distances, minDistTo1N)


def estimatePackingEpsilon(archive, sampleSize):
  """
  Compute estimates for the best epsilon estimate for this 
  archive's epsilon-packing.  The result is a tuple of stats
  for the distances of points within the archive:
    (min, Q1, median, mean, Q3, max)
  """
  archiveSize = len(archive)
  sampleSize = min(sampleSize, archiveSize-1)
  sampleDistances = [1E100]
    
  # Compute epsilon metrics for every point in the archive
  for idx in range(archiveSize):
    # Shuffle all the indexes to other points in
    # the archive, except this point
    shuffledIndexes = range(archiveSize)
    shuffledIndexes.remove(idx)
    np.random.shuffle(shuffledIndexes)

    # Compute all distances from idx to sampled
    # points from the archive
    for jdx in shuffledIndexes[0:sampleSize]:
      if (not idx == jdx):
        xi = archive[idx]
        xj = archive[jdx]
        sampleDistances.append( getL2Norm(xi,xj) )

  # Estimate the min, Q1, median, Q3, and max for these
  ##epsilonMetrics = np.percentile(sampleDistances, [0,25,50,75,100])
  ##epsilonMetrics.append( np.mean(sampleDistances) )
  ##epsilonMetrics.extend( np.percentile(sampleDistances, [75,100]) )

  # Return the estimation metrics for the epsilon-Packing
  return (min(sampleDistances))
  #return (tuple(epsilonMetrics))


def estimateCoverEpsilon(archive, sampleSize, n, sigma=0.0):
  """
  Compute estimates for the best epsilon estimate for this 
  archive's epsilon-cover.  The result is a tuple of stats
  for the distances of points within the archive:
    (min, Q1, median, mean, Q3, max)
  """
  archiveSize = len(archive)
  sampleSize = min(sampleSize, 2**n)
  sampleDistances = []

  for sampleIdx in range(sampleSize):
    archiveDistances = []
    for idx in range(archiveSize):
      xi = archive[idx]  
      xj = initialize(n, False, sigma)
      archiveDistances.append( getL2Norm(xi,xj) )
    sampleDistances.append( max(archiveDistances) )

  meanEpsilon = 0.0
  if (len(sampleDistances) > 1):
    meanEpsilon = np.mean(sampleDistances)
    
  # Estimate the min, Q1, median, Q3, and max for these
  ##epsilonMetrics = np.percentile(sampleDistances, [0,25,50,75,100])
  ##epsilonMetrics.append( np.mean(sampleDistances) )
  ##epsilonMetrics.extend( np.percentile(sampleDistances, [75,100]) )

  # Return the estimation metrics for the epsilon-Packing
  return (meanEpsilon)
  #return (tuple(epsilonMetrics))
    

def altArchiveReportHeader():
  print "Trial \t Generation \t CoverEpsilon \t PackingEpsilon \t ArchiveSize"

def altArchiveReport(archive, n, gen, trial, sampleSize, sigma):
  """
  Print output for the trial and various epsilon metrics
  """
  packingEps = estimatePackingEpsilon(archive, sampleSize)#[-1]
  coverEps   = estimateCoverEpsilon(archive, sampleSize, n, sigma)#[0]

  print int(trial), "\t", int(gen), "\t",\
        coverEps, "\t", packingEps, "\t",\
        len(archive)


def archiveReportHeader():
  print "Trial \t Generation \t MinSparseness \t AvgSparseness \t MaxSparseness \t ArchiveSize \t MinDistToN"

def archiveReport(archive, k, gen, trial):
  """
  Print output for the trial and various sparseness measures.
  """
  # Get sparseness measures
  n = len(archive[0])
  sparseness, minDistTo1N = getArchiveDistances(archive, k, n)

  # If we have a positive sparseness, report those results
  if len(sparseness) > 0:
    print int(trial), "\t", int(gen), "\t",\
          min(sparseness), "\t", np.mean(sparseness), "\t", max(sparseness), "\t",\
          len(archive), "\t", int(minDistTo1N)

  # Otherwise, report just the minimum distance to the all 1 string
  else:
    print int(trial), "\t", int(gen), "\t",\
          0.0, "\t", 0.0, "\t", float(n), "\t",\
          len(archive), int(minDistTo1N)

  # Flush standard out so we see the output in a timely fashion
  sys.stdout.flush()


def clearVisualizationDir(vizDirName):
  try:
    shutil.rmtree(vizDirName)
  except:
    print "Could not remove directory:", vizDirName

  try:
    os.mkdir(vizDirName)
  except:
    print "Could not make directory:", vizDirName

    
def writeVisualizationFile(vizDirName, gen, archive):
  """
  Write a file for reading and visualizing in Paraview
  """
  filename = "test" + str(gen) + ".csv"
  dirname = vizDirName.strip()
  fullPathFilename = os.path.join(dirname, filename)
  f = open(fullPathFilename, "w")
  f.write("x, y, z, idx\n")
  for idx in range(len(archive)):
    lineStr = ""
    #print "YY: ", gen, idx, py, rhoMin,
    for arg in archive[idx]:
      lineStr += str(arg) + ", "
    lineStr += str(idx) + '\n'
    f.write(lineStr)
  f.close()
  

        
def snsea(n, rhoMin, k, trial, pm=0.0, sigma=0.0, maxGenerations=100, allzero=True, vizDirName="visualizationData"):
  """
  This is the main routine for the program.  It takes the mutation probability information,
  the size of the string, the sparseness criteral and runs until maxGenerations is hit.
  We give it which trial it is working on so that the reporting method knows; however, the
  EA itself doesn't know or care.
  """
  # Initialize the individual and the archive
  x = initialize(n, allzero, sigma)
  archive = [x]

  # If we've not specified the probability of mutation, assume it is 1/n
  if (pm <= 0.0):
    pm = 1.0/float(n)

  # If we're writing for paraview visualization, clear the directory
  if (not vizDirName == "NOVIZ"):
    clearVisualizationDir(vizDirName)

  # Loop through generation counter
  for gen in range(maxGenerations):
    y = mutate(x, pm, sigma)
    py = computeSparseness(y, archive, k)
    if (py > rhoMin):
      archive.append( y )
      if (not vizDirName == "NOVIZ"):
        writeVisualizationFile(vizDirName, gen, archive)

    # Report results ever 100 generations
    if ( (gen % 100) == 0):
      altArchiveReport(archive, n, gen, trial, 10000, sigma)

    # Select an individual at random from the archive to serve
    # as a parent
    randSelectIdx = np.random.random_integers(low=0, high=(len(archive)-1))
    x = archive[randSelectIdx]

  # Return the archive, which is the solution in this case
  return (archive)


if __name__ == '__main__':
  configFileName = ""
  if (len(sys.argv) > 1):
    configFileName = sys.argv[1].strip()
  
  # Configuration parameters for command line and INI file
  configDefaults = {"n":32,\
                    "k":3,\
                    "rhoMin":0.25,\
                    "pm":0.0,\
                    "maxGenerations":50000,\
                    "numTrials":30,\
                    "sigma":0.0,\
                    "archiveDistancesFile":'archiveDistances.out',\
                    "vizDirName":"visualizationData"}
  configObj = configReader.buildArgObject(configFileName,'snsea',configDefaults,False)

  print
  print "Running SNS-EA ..."
  altArchiveReportHeader()
  for trial in range(configObj.numTrials):
    archive = snsea(configObj.n,\
                    configObj.rhoMin,\
                    configObj.k,\
                    trial,\
                    configObj.pm,\
                    configObj.sigma,\
                    maxGenerations=configObj.maxGenerations,\
                    allzero=True,
                    vizDirName=configObj.vizDirName)

  # Provide all the distance measures for the archive of the final
  # trial.
  if (False):
    print
    print "Collecting intra-archive sparseness distances ..."
    archiveDistances = getArchiveDistances(archive, configObj.k)
    archiveDistanceStrings = []
    for dist in archiveDistances:
      archiveDistanceStrings.append( 'XX: ' + str(dist) + '\n' )

  # Provide some statistics and output
  if False:
    if len(archiveDistances) > 0:
      output = ['\n']
      output.append("Statistics:\n")
      output.append("  Genome length:               " + str(configObj.n) + '\n')
      output.append("  Size of genome space:        " + str(2**configObj.n) + '\n')
      output.append("  rhoMin:                      " + str(configObj.rhoMin) + '\n')
      output.append("  Final archive size:          " + str(len(archive)) + '\n')
      output.append("  Min archive sparseness:      " + str(stats.tmin(archiveDistances)) + '\n')
      output.append("  Average distance in archive: " + str(stats.tmean(archiveDistances)) + '\n')
      output.append("  Median distance in archive:  " + str(np.median(archiveDistances)) + '\n')
      output.append("  Min archive sparseness:      " + str(stats.tmax(archiveDistances)) + '\n')
      output.append("  StdDev distances in archive: " + str(stats.tstd(archiveDistances)) + '\n')
      output.append("  IQR distances in archive:    " + str(stats.iqr(archiveDistances)) + '\n') 
      output.append('\n')
    else:
      output = ['\n', 'WARNING:  Archive contained only one candidate\n','\n']

    for myLine in output:
      sys.stdout.write( myLine )

    # Write the archive distances to a file
    print "Writing final archive distances to", configObj.archiveDistancesFile
    print
  
    f = open(configObj.archiveDistancesFile, 'w')
    f.writelines(archiveDistanceStrings)
    f.writelines(output)
    f.close()
  
