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
  upperBound = sigma * maxGenerations
  boundReportArg = None
  if not boundMutation:
    boundReportArg = (-upperBound, upperBound)

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

    # Report results every reportFreq generations
    if ((gen % reportFreq) == 0) :
      sb.archiveReport(archive, n, gen, trial, 10000, sigma, k, boundReportArg)

    # Select an individual at random from the archive to serve
    # as a parent
    #randSelectIdx = np.random.random_integers(low=0, high=(len(archive)-1))
    if (len(archive) <= 1):
      randSelectIdx = 0
    else:
      randSelectIdx = np.random.randint(low=0, high=(len(archive)-1))
    x = archive[randSelectIdx]

  # Return the archive, which is the solution in this case
  return (archive)


def snseaConvergenceTester(n, rhoMin, k, trial, pm=0.0, sigma=0.0, minGenerations=50, maxGenerations=2000,\
          allzero=True, reportFreq=100, boundMutation=True, convergenceFactor=2):
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
  upperBound = sigma * maxGenerations
  boundReportArg = None
  if not boundMutation:
    boundReportArg = (-upperBound, upperBound)

  # If we've not specified the probability of mutation, assume it is 1/n
  if (pm <= 0.0):
    pm = 1.0/float(n)

  # Loop through generation counter
  lastMinCoverGen = maxGenerations
  minCover = sys.float_info.max
  gen = 0
  # Keep looping as long as all of the following are true:
  #   1) We've not hit the max generations yet
  #   2) We've not "converged"
  #   3) It's still possible we might converge
  # --> "Convergence" is defined as having not updated the lowest cover
  #     in as many generations as it took to find it the first time.
  #     but we require the algorithm to go at least min generations.
  #     If we've passed the halfway mark to the max generations and
  #     we are still updating the cover, then we're not going to make it.
  while (gen < maxGenerations) and\
        ((gen < minGenerations) or (gen < convergenceFactor*lastMinCoverGen)) and\
        (lastMinCoverGen < maxGenerations/convergenceFactor):
    if boundMutation:
      y = sb.mutateIndividual(x, pm, sigma, (0,1), False)
    else:
      y = sb.mutateIndividual(x, pm, sigma, None, False)

    # Only compute sparseness and consider for admissio
    # if we've generated a new point.
    py = 0
    if not sb.isAlreadyInArchive(archive, y):
      py = sb.computeSparseness(y, archive, k)

    if (py >= rhoMin):
      archive.append( y )

    # Report results every reportFreq generations
    if ((gen % reportFreq) == 0) :
      sb.archiveReport(archive, n, gen, trial, 10000, sigma, k, boundReportArg)

    # Record the last time we saw the smallest cover
    coverEstimate = sb.estimateCoverEpsilon(archive, 10000, n, sigma)
    if coverEstimate < minCover:
      minCover = coverEstimate
      lastMinCoverGen = gen

    # Select an individual at random from the archive to serve
    # as a parent
    #randSelectIdx = np.random.random_integers(low=0, high=(len(archive)-1))
    if (len(archive) <= 1):
      randSelectIdx = 0
    else:
      randSelectIdx = np.random.randint(low=0, high=(len(archive)-1))
    x = archive[randSelectIdx]

    # Increment the generation counter
    gen += 1

  # Report the last generation where the minimum cover was found
  convStr = "CONVERGED"
  if (gen > maxGenerations/convergenceFactor):
    convStr = "NOT CONVERGED"
  print("YY: ", trial, '\t', minCover, '\t', lastMinCoverGen, '\t', maxGenerations, '\t', convStr)

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
                    "minGenerations":50,\
                    "maxGenerations":500,\
                    "numTrials":1,\
                    "startTrialNum":0,\
                    "sigma":0.0,\
                    "archiveFilename":'NOARCHIVEWRITE',\
                    "vizDirName":"NOVIZ",\
                    "reportFrequency":1,\
                    "boundMutation":True,\
                    "useEscapeSphere":False,\
                    "convergenceTest":False}
  configObj = configReader.buildArgObject(configFileName,'snsea',configDefaults,False)

  if (configObj.sigma <= 0.0) and (configObj.rhoMin >= configObj.n):
    configObj.rhoMin = np.ceil(0.75*np.log2(configObj.n))
    print("Setting rhoMin to", configObj.rhoMin)

  # Flush std I/O so that it prints early during long runs
  sys.stdout.flush()

  print()
  print("Running SNS-EA ...")

  # Print headers
  sb.archiveReportHeader()
  if configObj.convergenceTest:
    print("YY: Trial\tMinCover\tLastMinCoverGen\tMaxGen\tConvergeFlag")

  for trial in range(configObj.startTrialNum, configObj.startTrialNum+configObj.numTrials):
    if not configObj.convergenceTest:
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
      if (sb.isArchiveOutOfBounds(archive, (0,1))):
        print("Trial: ", trial, " archive contains points OUTOFBOUNDS")

    else:
      archive = snseaConvergenceTester(configObj.n,\
                      configObj.rhoMin,\
                      configObj.k,\
                      trial,\
                      configObj.pm,\
                      configObj.sigma,\
                      minGenerations=configObj.minGenerations,\
                      maxGenerations=configObj.maxGenerations,\
                      allzero=True,
                      reportFreq=configObj.reportFrequency,\
                      boundMutation=configObj.boundMutation)
  if (not configObj.archiveFilename == 'NOARCHIVEWRITE'):
    sb.writeArchive(archive, configObj.archiveFilename)

