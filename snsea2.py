import sys
import numpy as np
import configReader
import scipy.stats as stats


def initialize(n, allzero=False, sigma=0.0):
  x = np.array([0]*n)

  if (not allzero and sigma <= 0.0):
    x = np.random.random_integers(low=0, high=1, size=n)
  elif (not allzero and sigma > 0.0):
    x = np.random.uniform(low=0.0, high=1.0, size=n)
    
  return (x)



def mutate(x, pm, sigma=0.0):
  n = len(x)
  child = np.array([0]*n)
  
  if (sigma <= 0.0):
    mask = np.random.choice([0, 1], size=n, p=[1-pm,pm])
    y = x.copy()
    child = y ^ mask
  else:
    offsets = np.random.normal(scale=sigma, size=n)
    child = x + offsets
    
  return (child)


def computeSparseness(x, archive, k):
  distances = []
  for xi in archive:
    distances.append( sum( np.sqrt((x-xi)**2) ) )

  numElements = min(len(distances), k)
  distances.sort()
  sparseness = float(sum(distances[0:(numElements+1)]))/float(numElements)

  return (sparseness)


def getArchiveDistances(archive, k, n):
  minDistTo1N = n
  OneN = np.array([1.0]*n)
  distances = []
  for idx in range(len(archive)):
    minDistTo1N = min(minDistTo1N, sum( (archive[idx]-OneN)**2 ) )
    tmpDistances = []
    for jdx in range(len(archive)):
      if (not idx == jdx):
        xi = archive[idx]
        xj = archive[jdx]
        tmpDistances.append( sum( (xi-xj)**2 ) )

    numElements = min(len(tmpDistances), k)

    if numElements > 0:
      tmpDistances.sort()
      dist = float(sum(tmpDistances[0:(numElements+1)]))/float(numElements)
      distances.append( dist)
    
  return (distances, minDistTo1N)
      

def archiveReport(archive, k, gen, trial):
  n = len(archive[0])
  sparseness, minDistTo1N = getArchiveDistances(archive, k, n)

  if len(sparseness) > 0:
    print int(trial), "\t", int(gen), "\t",\
          min(sparseness), "\t", np.mean(sparseness), "\t", max(sparseness), "\t",\
          len(archive), "\t", int(minDistTo1N)
  else:
    print int(trial), "\t", int(gen), "\t",\
          0.0, "\t", 0.0, "\t", float(n), "\t",\
          len(archive), int(minDistTo1N)

  sys.stdout.flush()


        
def snsea(n, rhoMin, k, trial, pm=0.0, sigma=0.0, maxGenerations=100, allzero=True):
  x = initialize(n, allzero, sigma)
  archive = [x]
  #print "  Generation:", 0, "   ArchiveSize:", len(archive)

  if (pm <= 0.0):
    pm = 1.0/float(n)

  for gen in range(maxGenerations):
    y = mutate(x, pm, sigma)
    py = computeSparseness(y, archive, k)
    if (py > rhoMin):
      archive.append( y )
      ## --- vvv --- Temporary --- vvv -----
      f = open("test" + str(gen) + ".csv", "w")
      f.write("x, y, z, idx\n")
      for idx in range(len(archive)):
        lineStr = ""
        #print "YY: ", gen, idx, py, rhoMin,
        for arg in archive[idx]:
          lineStr += str(arg) + ", "
        lineStr += str(idx) + '\n'
        f.write(lineStr)
      f.close()
        #  print arg,
        #print
      ## --- ^^^ --- Temporary --- ^^^ -----
      #print "  Generation:", gen, "   ArchiveSize:", len(archive)

    if ( (gen % 100) == 0):
      archiveReport(archive, k, gen, trial)
    
    randSelectIdx = np.random.random_integers(low=0, high=(len(archive)-1))
    x = archive[randSelectIdx]

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
                    "archiveDistancesFile":'archiveDistances.out'}
  configObj = configReader.buildArgObject(configFileName,'snsea',configDefaults,False)

  print
  print "Running SNS-EA ..."
  print "Trial \t Generation \t MinSparseness \t AvgSparseness \t MaxSparseness \t ArchiveSize \t MinDistToN"
  for trial in range(configObj.numTrials):
    archive = snsea(configObj.n,\
                    configObj.rhoMin,\
                    configObj.k,\
                    trial,\
                    configObj.pm,\
                    configObj.sigma,\
                    maxGenerations=configObj.maxGenerations,\
                    allzero=True)

  if (False):
    print
    print "Collecting intra-archive sparseness distances ..."
    archiveDistances = getArchiveDistances(archive, configObj.k)
    archiveDistanceStrings = []
    for dist in archiveDistances:
      archiveDistanceStrings.append( 'XX: ' + str(dist) + '\n' )

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
    
    print "Writing final archive distances to", configObj.archiveDistancesFile
    print
  
    f = open(configObj.archiveDistancesFile, 'w')
    f.writelines(archiveDistanceStrings)
    f.writelines(output)
    f.close()
  
