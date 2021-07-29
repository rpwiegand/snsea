import sys, re, os

def printMergedDataset(mergedDataset):
  """
  This routine takes a mergedDataset and prints it to the
  screen in order, with a header.
  """
  # Print the header
  header = "Trial\tGeneration\tCoverEpsilon\tPackingEpsilon\tMinArchiveSparseness\tArchiveSize"
  print(header)

  for itemKey in mergedDataset:
    trial,generation,cover,packing,sparseness,size = mergedDataset[itemKey]
    print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(trial,generation,cover,packing,sparseness,size))

  
def readFileIntoDictionary(filename):
  """
  This routine takes a file name and reads the file into a dictionary
  structured as follows.  The key of the dictionary is a tuple containing
  (trial, generation).  The value is a tuple contaiing:
  (trial,generation,cover,packing,sparseness,size).
  The routine returns this dictionary, as well as the maximum generation
  without performing any cleaning.
  """
  # Read all contents of the file into an array of strings.
  dataFileHandler = open(filename)
  dataTextLines = dataFileHandler.readlines()
  dataFileHandler.close()

  # Initialize the data dictionary and max values
  dataDict = {}
  maxGeneration = 0

  # 
  for textLine in dataTextLines:
    lineItems = textLine.strip().split()

    try:
      # Read the line items into their constituent variables
      #  the lines start with "XX: ", so ignore lineItems[0]
      trial = int(lineItems[1])
      generation = int(lineItems[2])
      cover = float(lineItems[3])
      packing = float(lineItems[4])
      sparseness = float(lineItems[5])
      size = int(lineItems[6])

      # Updated the data dictionary, keyed on the trial and generation
      dataDict[(trial,generation)] = (trial,generation,cover,packing,sparseness,size) 

      # Update the max values
      if (generation > maxGeneration):
        maxGeneration = generation

    # Ignore any lines that cannot be properly parsed
    except:
      pass

  return dataDict, maxGeneration


def isRecordBad(recordTuple, maxAcceptableGen=350):
  # Initilizations
  trial,generation,cover,packing,sparseness,size = recordTuple
  badRecord = (generation == maxAcceptableGen) and (sparseness < 0)
  return badRecord

  
def cleanDataSetOfFailedRuns(dataset, maxGeneration):
  """
  This routne looks to see if the last generation of each trial
  has bad measure data, then excludes all the data for that trial.
  It returns a cleaned dataset, as well as a list of the excluded
  trials.
  """
  # Look through each record to find the bad trials
  excludedTrials = set()
  for itemKey in dataset:
    trial, generation = itemKey
    recordTuple = dataset[itemKey]
    if isRecordBad(recordTuple):
      excludedTrials.add(trial)

  # Rebuild a new dataset excluding bad trials
  cleanedDataset = dict()
  for itemKey in dataset:
    trial, generation = itemKey
    if not (trial in excludedTrials):
      cleanedDataset[itemKey] = dataset[itemKey]

  return cleanedDataset, excludedTrials


def readAndCleanAllDataFiles(filenameList):
  """
  This routine takes a list of filenames and produces a list of
  datasets that have been read and cleaned of bad runs.
  """
  
  # Initialize the list of datafiles
  datasetList = []

  #Spin through all the filenames in the list
  for filename in filenameList:
    # Try to read and clean the data, then put it in the list
    try:
      df, mg = readFileIntoDictionary(filename)
      ndf, et = cleanDataSetOfFailedRuns(df, mg)
      datasetList.append(ndf)

    # Report any failures
    except:
      print("Could not read or parse file: " + filename)
      pass
    
  return  datasetList


def mergeAllDatasets(datasetList):
  """
  Take a list of datasets and merge them so that the trial counters make sense
  and are unique.  Return the merged dictionary.
  """
  revisedTrialCounter = -1
  mergedDataset = dict()
  previousTrial = -1

  for dataset in datasetList:
    for itemKey in dataset:
      if revisedTrialCounter >= 50:
        continue
      
      # Get the old record
      trial,generation,cover,packing,sparseness,size = dataset[itemKey]
      
      # Check to see if we must increment the trial counter
      if not (trial == previousTrial):
        revisedTrialCounter += 1
        previousTrial = trial
        
      # Change the trial value to the current counter
      trial = revisedTrialCounter

      # Add the new record to the merged dataset
      mergedDataset[(trial,generation)] = (trial,generation,cover,packing,sparseness,size)


  return mergedDataset
      

      
if __name__== "__main__":
   filenameList = sys.argv[1:]
   dfList = readAndCleanAllDataFiles(filenameList)
   mdf = mergeAllDatasets(dfList)
   printMergedDataset(mdf)
