import sys
import numpy as np
import configReader
import scipy.stats as stats
import os, shutil

def initialize(n, allzero=False, sigma=0.0, bounds=(0,1)):
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

  # If not allzero and we're numeric, init in [lp,ub]^n
  elif (not allzero and sigma > 0.0):
    x = np.random.uniform(low=bounds[0], high=bounds[1], size=n)

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
      inside = False
  return (inside)


def inSpecialArea(parent, criteria, variance):
  special=False
  if (getDistance(parent, np.array([criteria] * len(parent))) < variance) or\
     (parent[1] < 0) or (parent[2] < 0) or (parent[0] < 0):
    special=True
  return (special)
    
  
def mutate(x, pm, sigma=0.0, bound=[0,1], useEscapeSphere=False):
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
      if stuckCount > 50:
        #raise Exception("The re-mutate loop is stuck ... should not take this many iterations to find a good point.")
        inside=True
      stuckCount += 1

      # Generate a mutation
      offsets = np.random.normal(scale=sigma, size=n)
      child = x + offsets

      # Check if the mutation is inside the bounds
      if (bound==None):
        inside=True
      elif (useEscapeSphere) and (inSpecialArea(x, 0.25, .2)):
        inside=True
      else:
        inside = isInsideEuclideanBoundary(child, bound[0], bound[1])

  # Return the child
  return (child)


def isBinary(x):
  """
  Determine if the array is made up of just 0's and 1's.
  """
  binary = True
  for arg in x:
    if (arg > 0) and (arg < 1):
      binary = False
  return (binary)


def getDistance(x1, x2):
  """
  Return the distance betweeen two points.  Use L2 norm (Euclidean) distance
  for real values and Hamming for binary spaces.
  """
  distance = sum( (x1-x2)**2 )

  # If this is a real-value, use L2, otherwise using Hamming distance
  if not (isBinary(x1) and isBinary(x2)):
    distance = np.sqrt(distance)
    
  return (distance)


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
    distances.append( getDistance(x,xi) )

  # Arrange things so that we can get either
  # the k closest values or the size of the 
  # sarchive, whichever is smaller
  numElements = min(len(distances), k)
  distances.sort()

  # Estimate sparseness as the average over that subset, then return
  sparseness = float(sum(distances[0:(numElements+1)]))/float(numElements)
  return (sparseness)


def getPairwiseSparsenessMetrics(archive, k, n):
  """
  Get minimum archive distances according to the sparseness metric,
  as well as minimum distance to 1^n.
  """
  # Setup the arbitrary "piton" measure (the all 1 binary string)
  minDistTo1N = n
  maxDistInArchive = 0
  OneN = np.array([1.0]*n)

  # Initialize the distances
  distances = []

  # Compute pair-wise distances of all points in the archive.
  # Also, keep track of the distance to the all 1 string, 1^n.
  for idx in range(len(archive)):
    minDistTo1N = min( minDistTo1N, getDistance(archive[idx], OneN) )
    tmpDistances = []
    for jdx in range(len(archive)):
      if (not idx == jdx):
        xi = archive[idx]
        xj = archive[jdx]
        dist = getDistance(xi,xj)
        tmpDistances.append( dist )
        maxDistInArchive = max( [maxDistInArchive, dist] )
          
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
  return (distances, minDistTo1N, maxDistInArchive)


def estimatePackingEpsilon(archive, sampleSize, maxPacking=sys.float_info.max):
  """
  Compute the best epsilon estimate for this archive's epsilon-packing.
  The result is the maximum distance between any two points in the archive
  divided by 2.  See epsilon-Packing definition in KNN literature.
    --> Larger epsilon mean that the archive is more efficiently spread out
  """
  archiveSize = len(archive)
  sampleSize = min(sampleSize, archiveSize-1)
  sampleDistances = [0]#[maxPacking]
    
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
        sampleDistances.append( getDistance(xi,xj) )

  # Return the estimation metrics for the epsilon-Packing
  #   Subset Y of U is a eps-Packing iff D(x,y) > 2eps for all x,y \in Y
  return (max(sampleDistances)/2)


def estimateCoverEpsilon(archive, sampleSize, n, sigma=0.0, bounds=(0,1)):
  """
  Compute the best epsilon estimate for this archive's epsilon-cover.
  We do this by sampling the whole space and finding the the point
  with the smallest distance to that point to *any* point in the archive.  
  Our estimate is the maximum of all such distances.  We also include the
  1^n string as one of those points as a kind of "piton" measure since our
  algorithms tend to start at 0^n.  See epsilon-cover definition in KNN 
  literature.
    --> Smaller epsilon means fewer points are "close" to whole space
  """
  archiveSize = len(archive)
  sampleSize = min(sampleSize, 2**n)
  sampleDistances = []
  xi, xj, xk = (None, None, None)

  for sampleIdx in range(sampleSize):
    # Sample a point from the space
    xj = initialize(n, False, sigma, bounds)
      
    # Find the closest point in archive to a random point
    archiveDistances = []
    for idx in range(archiveSize):
      xi = archive[idx]  
      archiveDistances.append( getDistance(xi,xj) )
    sampleDistances.append( min(archiveDistances) )

    # Find the closest points in archive to 1^n
    xk = np.array([1]*len(xi))    
    archiveDistances = []
    for idx in range(archiveSize):
      xi = archive[idx]  
      archiveDistances.append( getDistance(xi,xk) )
    sampleDistances.append( min(archiveDistances) )

  # Our epsilon estimate is the *maximum* of those closest points
  epsilon = 0.0
  if (len(sampleDistances) > 1):
    epsilon = max(sampleDistances)
    
  # Return the estimation metrics for the epsilon-Packing
  #   Subet Y of U is an eps-cover if for every x \in U,
  #   there is some y \in Y where D(x,y) < eps
  return (epsilon)


#################################################################
# Example: There exists a 2-net for U={0,1}^4 with archive of  #
#           size 6.                                            #
#                                                              #
#      A = {0000, 0011, 1100, 0110, 1001, 1111}                #
#    |A| = 6                                                   #
#                                                              #
#  1) No point in {0,1}^4 is more than 2 away from some point  #
#     in A.  So the archive is a 2-Cover                       #
#                                                              #
#  2) The largest distance between any two points in A is 4,   #
#     and 4/2 = 2.  So the archive is a 2-Packing.             #
#                                                              #
# Therefore this archive is an 2-*optimal* archive.            #
#################################################################


def altArchiveReportHeader():
  print "XX: Trial \t Generation \t CoverEpsilon \t PackingEpsilon \t MinArchiveSparseness \t ArchiveSize"

def altArchiveReport(archive, n, gen, trial, sampleSize, sigma, k, bounds):
  """
  Print output for the trial and various epsilon metrics
  """
  maxPackingDistance = getDistance(np.array([0]*n), np.array([1]*n))
  packingEps = estimatePackingEpsilon(archive, sampleSize, maxPackingDistance)
  coverEps   = estimateCoverEpsilon(archive, sampleSize, n, sigma)
  maxDistInArchive = 0
  
  minSparse = -1
  if (len(archive) > 1):
    sparsenessVals, minAllOne, maxDistInArchive  = getPairwiseSparsenessMetrics(archive, k, n)
    minSparse = min(sparsenessVals)
    
  print "XX:", int(trial), '\t', int(gen), '\t',\
        coverEps, '\t', packingEps, '\t', \
        minSparse, '\t', \
        len(archive)

  # Flush standard out so we see the output in a timely fashion
  sys.stdout.flush()


def archiveReportHeader():
  print "Trial \t Generation \t MinSparseness \t AvgSparseness \t MaxSparseness \t ArchiveSize \t MinDistToN"

def archiveReport(archive, k, gen, trial):
  """
  Print output for the trial and various sparseness measures.
  """
  # Get sparseness measures
  n = len(archive[0])
  sparseness, minDistTo1N, maxDistInArchive = getPairwiseSparsenessMetrics(archive, k, n)

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
  filename = "archive" + str(gen) + ".csv"
  dirname = vizDirName.strip()
  fullPathFilename = os.path.join(dirname, filename)
  f = open(fullPathFilename, "w")
  f.write("x, y, z, idx\n")
  for idx in range(len(archive)):
    lineStr = ""
    for arg in archive[idx]:
      lineStr += str(arg) + ", "
    lineStr += str(idx) + '\n'
    f.write(lineStr)
  f.close()
  

def isAlreadyInArchive(archive, x):
  """
  Check to see if the candidate is already in the archive.  
  """
  already = False
  for xi in archive:
    if tuple(xi) == tuple(x):
      already=True
  return (already)


def isArchiveOutOfBounds(archive, bounds):
  """
  Check to see if there are any points n the archive that are
  outside the specified bounds.
  """
  outOfBounds = False
  for x in archive:
    for arg in x:
      if (arg > max(bounds)):
        outOfBounds = True
      elif (arg < min(bounds)):
        outOfBounds = True
  return(outOfBounds)

        
def snsea(n, rhoMin, k, trial, pm=0.0, sigma=0.0, maxGenerations=100, allzero=True, \
          vizDirName="NOVIZ", reportFreq=100, boundMutation=True, useEscapeSphere=False):
  """
  This is the main routine for the program.  It takes the mutation probability information,
  the size of the string, the sparseness criteral and runs until maxGenerations is hit.
  We give it which trial it is working on so that the reporting method knows; however, the
  EA itself doesn't know or care.
  """
  # Initialize the individual and the archive
  x = initialize(n, allzero, sigma)
  y = initialize(n, allzero, sigma)
  archive = [x]

  # If we've not specified the probability of mutation, assume it is 1/n
  if (pm <= 0.0):
    pm = 1.0/float(n)

  # If we're writing for paraview visualization, clear the directory
  if (not vizDirName == "NOVIZ"):
    clearVisualizationDir(vizDirName)

  # Loop through generation counter
  for gen in range(maxGenerations):
    if boundMutation:
      y = mutate(x, pm, sigma, (0,1), useEscapeSphere)
    else:
      y = mutate(x, pm, sigma, None, useEscapeSphere)

    # Only compute sparseness and consider for admissio
    # if we've generated a new point.
    py = 0
    if not isAlreadyInArchive(archive, y):
      py = computeSparseness(y, archive, k)

    if (py >= rhoMin):
      archive.append( y )
      if (not vizDirName == "NOVIZ"):
        writeVisualizationFile(vizDirName, gen, archive)

    # Report results ever 100 generations
    if ( (gen % reportFreq) == 0) and (boundMutation):
      altArchiveReport(archive, n, gen, trial, 10000, sigma, k, None)
    else:
      upperBound = sigma * maxGenerations
      altArchiveReport(archive, n, gen, trial, 10000, sigma, k, (-upperBound, upperBound))

    # Select an individual at random from the archive to serve
    # as a parent
    randSelectIdx = np.random.random_integers(low=0, high=(len(archive)-1))
    x = archive[randSelectIdx]

  # Return the archive, which is the solution in this case
  return (archive)


def writeArchive(archive, archiveFilename):
  """
  Write the archive out to the specified file.
  """
  archiveStrings = []
  for x in archive:
    outStr = ''
    for idx in range(len(x)-1):
      outStr += str(x[idx]) + ','
    outStr += str(x[-1]) + '\n'
    archiveStrings.append( outStr )
    
  f = open(archiveFilename, 'w')
  f.writelines(archiveStrings)
  f.close()
    
  
if __name__ == '__main__':
  configFileName = ""
  if (len(sys.argv) > 1):
    configFileName = sys.argv[1].strip()
  
  # Configuration parameters for command line and INI file
  configDefaults = {"n":32,\
                    "k":3,\
                    "rhoMin":2.0,\
                    "pm":0.0,\
                    "maxGenerations":50000,\
                    "numTrials":1,\
                    "sigma":0.0,\
                    "archiveFilename":'NOARCHIVEWRITE',\
                    "vizDirName":"NOVIZ",\
                    "reportFrequency":1,\
                    "boundMutation":True,
                    "useEscapeSphere":False}
  configObj = configReader.buildArgObject(configFileName,'snsea',configDefaults,False)
  
  # Flush std I/O so that it prints early during long runs
  sys.stdout.flush()

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
                    vizDirName=configObj.vizDirName,
                    reportFreq=configObj.reportFrequency,\
                    boundMutation=configObj.boundMutation,
                    useEscapeSphere=configObj.useEscapeSphere)
    if (isArchiveOutOfBounds(archive, (0,1))):
      print "Trial: ", trial, " archive contains points OUTOFBOUNDS"

  if (not configObj.archiveFilename == 'NOARCHIVEWRITE'):
    writeArchive(archive, configObj.archiveFilename)
   
