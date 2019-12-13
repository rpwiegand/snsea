import sys
import numpy as np
import configReader
import scipy.stats as stats
import os, shutil
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
  

def generateChildren(parents, llambda, pm, sigma, boundMutation):
  children = []
  mu = len(parents)

  for idx in range(llambda):
    # Keeping looping through the parent population to produce new children
    # until we've filled the child population
    parentIdx = idx % mu

    # Get a parent and produce a child via mutation
    x = parents[parentIdx]
    if boundMutation:
      y = sb.mutateIndividual(x, pm, sigma, (0,1), False)
      children.append( y )
    else:
      y = sb.mutateIndividual(x, pm, sigma, None, False)
      children.append( y )

  return children


def selectNewParents(population, k, mu):
  """
  Return the mu best new parents.
  """
  sortablePopulation = []

  # Compute sparsness of every individual relative to the others
  for y in population:
    py = sb.computeSparseness(y, population, k)
    sortablePopulation.append( (py, y) )

  # Sort these, then strip out the mu best individuals
  sortablePopulation.sort(key=lambda tup: tup[0])
  newParents = list(zip(*sortablePopulation))[1][0:mu]
  
  return list(newParents)

      
def snsea(n, mu, llambda, rhoMin, k, trial, pm=0.0, sigma=0.0, maxGenerations=100, allzero=True, \
          reportFreq=100, boundMutation=True, plusStrategy=True):
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

  # Loop through generation counter
  for gen in range(maxGenerations):
    # Parents have some kids!
    children = generateChildren(parents, llambda, pm, sigma, boundMutation)

    # Update the archive with any children who meet the criteria
    updateArchive(children, archive, k, rhoMin)

    # Report results ever 100 generations
    if ( (gen % reportFreq) == 0) and (boundMutation):
      sb.archiveReport(archive, n, gen, trial, 10000, sigma, k, None)
    else:
      sb.archiveReport(archive, n, gen, trial, 10000, sigma, k, (-upperBound, upperBound))

    parents = selectNewParents(children + parents, k, mu)

  # Return the archive, which is the solution in this case
  return (archive)

  
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
                    "usePlusStrategy":False}
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
                    plusStrategy=configObj.usePlusStrategy)
   
