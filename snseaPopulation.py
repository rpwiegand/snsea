import sys
import numpy as np
import configReader
import scipy.stats as stats
import os, shutil, os.path
import snseaBase as sb


def initializePopulation(n, popSize, sigma, bounds=(0,1)):
  """
  Create an initial population of size popSize 
  """
  population = []
  for idx  in range(popSize):
    x = sb.initializeIndividual(n, False, sigma, bounds)
    population.append( x )

  return population


def updateArchive(population, archive, k, rhoMin):
  """
  Given a set of individuals, add any to the archive that meet the sparseness
  criterion
  """
  for y in population:
    if not sb.isAlreadyInArchive(archive, y):
      py = sb.computeSparseness(y, archive, k)

      # Is the sparseness at least as large as rhoMin?
      if (py >= rhoMin):
        archive.append( y )
  

def generateChildren(parents, llambda, pm, sigma, boundMutation, archive, archiveSelectProb):
  children = []
  mu = len(parents)
  parentCount = 0

  for idx in range(llambda):    
    x = None

    if (np.random.uniform() < archiveSelectProb):
      # If we're pulling from the archive, then grab a random parent from there
      archiveIdx = np.random.randint(0, len(archive))
      x = archive[archiveIdx]
      #print(".", archiveIdx, x)
    
    else:
      # Keeping looping through the parent population to produce new children
      # until we've filled the child population
      parentIdx = parentCount % mu
      x = parents[parentIdx]
      parentCount += 1
      #print("o", parentIdx, x, parents)

    if boundMutation:
      y = sb.mutateIndividual(x, pm, sigma, (0,1), False)
      children.append( y )
    else:
      y = sb.mutateIndividual(x, pm, sigma, None, False)
      children.append( y )

    #print("DBGGG:  x=", x, "     y=", y)
  return children


def selectNewParents(population, k, mu, compareSet=None):
  """
  Return the mu best new parents.
  """
  sortablePopulation = []

  # compareSet is the set of individuals
  # over which will will compute sparseness
  # for *fitness purposes only*.  By default,
  # this is the population itself.  But we
  # can choose a different baseline ... the archive,
  # for example
  if compareSet==None:
    compareSet = population

  # Compute sparsness of every individual relative to the others
  for y in population:
    py = sb.computeSparseness(y, compareSet, k)
    sortablePopulation.append( (py, y) )

  # Sort these, then strip out the mu best individuals
  sortablePopulation.sort(key=lambda tup: tup[0])
  newParents = list(zip(*sortablePopulation))[1][0:mu]
  
  return list(newParents)

      
def snsea(n, mu, llambda, rhoMin, k, trial, pm=0.0, sigma=0.0, maxGenerations=100, allzero=True, \
          reportFreq=100, boundMutation=True, plusStrategy=False, vizDirName="NONE", msrChildren=False,
          fitnessByArchive=False, convergenceTest=False, archiveSelectProb=0.0, archiveFilename="NOARCHIVEWRITE"):
  """
  This is the main routine for the program.  It takes the mutation probability information,
  the size of the string, the sparseness criteral and runs until maxGenerations is hit.
  We give it which trial it is working on so that the reporting method knows; however, the
  EA itself doesn't know or care.
  """
  # Initialize the population
  parents = initializePopulation(n, mu, sigma)

  # Initialize the archive with {0}^n and any parents who might qualify
  archive = [ sb.initializeIndividual(n, allzero, sigma) ]
  updateArchive(parents, archive, k, rhoMin)

  # If we've not specified the probability of mutation, assume it is 1/n
  if (pm <= 0.0):
    pm = 1.0/float(n)

  # Estimate the upper limit for continuous spaces to get to in any dimension
  upperBound = sigma * maxGenerations

  # These are only used if we're testing for convergence
  lastMinCoverGen = 0
  minCover = sys.float_info.max

  # Loop through generation counter
  for gen in range(maxGenerations):
    # Parents have some kids!
    children = generateChildren(parents, llambda, pm, sigma, boundMutation, archive, archiveSelectProb)

    # Update the archive with any children who meet the criteria
    updateArchive(children, archive, k, rhoMin)

    # Report results every reportFreq generations
    coverEst = None
    if ( (gen % reportFreq) == 0) and (boundMutation):
      coverEst = sb.archiveReport(archive, n, gen, trial, 10000, sigma, k, None)
    elif (msrChildren):
      coverEst = sb.archiveReport(archive, n, gen, trial, 10000, sigma, k, (-upperBound, upperBound), children)
    else:
      coverEst = sb.archiveReport(archive, n, gen, trial, 10000, sigma, k, (-upperBound, upperBound))

    # If we're testing for convergence, see if we've seen the lowest cover so far
    if convergenceTest and (coverEst == None):  # No need to estimate twice, if we've already reported
      coverEst = sb.estimateCoverEpsilon(archive, 10000, n, sigma)
    if coverEst < minCover: # Keep track of the min, including which gen
      minCover = coverEst
      lastMinCoverGen = gen

    if (not vizDirName == "NONE") and (not vizDirName == "NOVIZ"):
       writeVisualizationFilePop(vizDirName, trial, gen, archive, children, parents)

    # Allow the user to specify whether selection is based on sparseness
    # of population members compared to the population itself (default) or
    # the archive
    compareSet = None
    if fitnessByArchive:
      compareSet = archive

    # Use plus or comma strategy.  Comma is the default.
    selectSet = children
    if plusStrategy:
      selectSet= children + parents

    # Select individuals to go into the next generation
    parents = selectNewParents(selectSet, k, mu, compareSet)

    # Ugly bail out to shorten runs that have converged
    if convergenceTest and (gen > 100) and (lastMinCoverGen < gen/2):
      break


  # If we're testing for convergence, then report that
  if convergenceTest and (lastMinCoverGen < maxGenerations/2):
    print("YY: ", trial, '\t', minCover, '\t', lastMinCoverGen, '\t', maxGenerations, '\tCONVERGED')
  elif convergenceTest:
    print("YY: ", trial, '\t', minCover, '\t', lastMinCoverGen, '\t', maxGenerations, '\tNOT CONVERGED')

  # Maybe write out the archive?
  if (not vizDirName == "NONE") and (not vizDirName == "NOVIZ") and (not archiveFilename == "NOARCHIVEWRITE"):
    sb.writeVisualizationFile(vizDirName, trial, archive, archiveFilename)

  # Return the archive, which is the solution in this case
  return (archive)


def writeVisualizationFilePop(vizDirName, trial, gen, archive, children, parents):
  """
  Write a file for reading and visualizing in Paraview
  """
  # Get file handle read for writing
  filename = "viz-archive-and-pop.csv"
  dirname = vizDirName.strip()
  fullPathFilename = os.path.join(dirname, filename)
  
  # Open file and write the header
  if (not os.path.exists(fullPathFilename)):
    f = open(fullPathFilename, "w")
    f.write("trial,generation,whichPop,x,y,idx\n")
  else:
    f = open(fullPathFilename, "a+")

  # Write the archive
  for idx in range(len(archive)):
    lineStr = str(trial) + ","
    lineStr += str(gen) + ",archive,"
    lineStr += str(archive[idx][0]) + ','
    lineStr += str(archive[idx][1]) + ','    
    lineStr += str(idx) + '\n'
    f.write(lineStr)

  # Write the child population
  for idx in range(len(children)):
    lineStr = str(trial) + ","
    lineStr += str(gen) + ",children,"
    lineStr += str(children[idx][0]) + ','
    lineStr += str(children[idx][1]) + ','    
    lineStr += str(idx) + '\n'
    f.write(lineStr)

  # Write the parent population
  for idx in range(len(parents)):
    lineStr = str(trial) + ","
    lineStr += str(gen) + ",parents,"
    lineStr += str(parents[idx][0]) + ','
    lineStr += str(parents[idx][1]) + ','    
    lineStr += str(idx) + '\n'
    f.write(lineStr)
    
  f.close()
  
  
if __name__ == '__main__':
  configFileName = ""
  if (len(sys.argv) > 1):
    configFileName = sys.argv[1].strip()
  
  # Configuration parameters for command line and INI file
  configDefaults = {"n":32,\
                    "k":3,\
                    "mu":5,\
                    "llambda":35,\
                    "rhoMin":2.0,\
                    "pm":0.0,\
                    "maxGenerations":500,\
                    "numTrials":1,\
                    "startTrialNum":0,\
                    "sigma":0.0,\
                    "reportFrequency":1,\
                    "boundMutation":True,\
                    "usePlusStrategy":False,\
                    "vizDirName":"NONE",\
                    "measureChildren":False,\
                    "fitnessByArchive":False,\
                    "convergenceTest":False,\
                    "archiveSelectProb":0.0}                      
  configObj = configReader.buildArgObject(configFileName, 'snsea',configDefaults,False)
  
  # Flush std I/O so that it prints early during long runs
  sys.stdout.flush()

  print()
  print("Running population-based SNS-EA ...")
  sb.archiveReportHeader()
  for trial in range(configObj.startTrialNum, configObj.startTrialNum+configObj.numTrials):
    archive = snsea(configObj.n,\
                    configObj.mu,\
                    configObj.llambda,\
                    configObj.rhoMin,\
                    configObj.k,\
                    trial,\
                    configObj.pm,\
                    configObj.sigma,\
                    maxGenerations=configObj.maxGenerations,\
                    allzero=True,\
                    reportFreq=configObj.reportFrequency,\
                    boundMutation=configObj.boundMutation,\
                    plusStrategy=configObj.usePlusStrategy,
                    vizDirName=configObj.vizDirName,\
                    msrChildren=configObj.measureChildren,\
                    fitnessByArchive=configObj.fitnessByArchive,\
                    convergenceTest=configObj.convergenceTest,
                    archiveSelectProb=configObj.archiveSelectProb,
                    archiveFilename=configObj.archiveFilename)
   
