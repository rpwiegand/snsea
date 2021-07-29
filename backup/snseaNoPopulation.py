import sys
import numpy as np
import configReader
import scipy.stats as stats
import os, shutil
import snseaBase as sb


        
def snsea(n, rhoMin, k, trial, pm=0.0, sigma=0.0, maxGenerations=100, allzero=True, \
          vizDirName="NOVIZ", reportFreq=100, boundMutation=True, useEscapeSphere=False):
  """
  This is the main routine for the program.  It takes the mutation probability information,
  the size of the string, the sparseness criteral and runs until maxGenerations is hit.
  We give it which trial it is working on so that the reporting method knows; however, the
  EA itself doesn't know or care.
  """
  # Initialize the individual and the archive
  x = sb.initializeIndividual(n, allzero, sigma)
  y = sb.initializeIndividual(n, allzero, sigma)
  archive = [x]

  # If we've not specified the probability of mutation, assume it is 1/n
  if (pm <= 0.0):
    pm = 1.0/float(n)

  # If we're writing for paraview visualization, clear the directory
  if (not vizDirName == "NOVIZ"):
    sb.clearVisualizationDir(vizDirName)

  # Loop through generation counter
  for gen in range(maxGenerations):
    if boundMutation:
      y = sb.mutateIndividual(x, pm, sigma, (0,1), useEscapeSphere)
    else:
      y = sb.mutateIndividual(x, pm, sigma, None, useEscapeSphere)

    # Only compute sparseness and consider for admissio
    # if we've generated a new point.
    py = 0
    if not sb.isAlreadyInArchive(archive, y):
      py = sb.computeSparseness(y, archive, k)

    if (py >= rhoMin):
      archive.append( y )
      if (not vizDirName == "NOVIZ"):
        sb.writeVisualizationFile(vizDirName, gen, archive)

    # Report results ever 100 generations
    if ( (gen % reportFreq) == 0) and (boundMutation):
      sb.archiveReport(archive, n, gen, trial, 10000, sigma, k, None)
    else:
      upperBound = sigma * maxGenerations
      archiveReport(archive, n, gen, trial, 10000, sigma, k, (-upperBound, upperBound))

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
                    "rhoMin":2.0,\
                    "pm":0.0,\
                    "maxGenerations":500,\
                    "numTrials":1,\
                    "startTrialNum":0,\
                    "sigma":0.0,\
                    "archiveFilename":'NOARCHIVEWRITE',\
                    "vizDirName":"NOVIZ",\
                    "reportFrequency":1,\
                    "boundMutation":True,\
                    "useEscapeSphere":False}
  configObj = configReader.buildArgObject(configFileName,'snsea',configDefaults,False)
  
  # Flush std I/O so that it prints early during long runs
  sys.stdout.flush()

  print
  print "Running SNS-EA ..."
  sb.archiveReportHeader()
  for trial in range(configObj.startTrialNum, configObj.startTrialNum+configObj.numTrials):
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
                    boundMutation=configObj.boundMutation,\
                    useEscapeSphere=configObj.useEscapeSphere)
    if (isArchiveOutOfBounds(archive, (0,1))):
      print "Trial: ", trial, " archive contains points OUTOFBOUNDS"

  if (not configObj.archiveFilename == 'NOARCHIVEWRITE'):
    sb.writeArchive(archive, configObj.archiveFilename)
   
