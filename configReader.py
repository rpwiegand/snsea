"""
This module provides a simple utility to help get parameters from ini
files and the comand line for running the various programs in this
repository.  The main function is buildArgObject, which takes config
file and command line parameter information, as well as a defaults
dictionary, and produces an object with fields matching the requested
fields in the defaults dictionary.  Command-line parameters are
preferred over config file entries, which are preferred over the defaults.
"""

import os, sys
from ConfigParser import SafeConfigParser
import re

class empty(object):
  pass

GLOBAL_QUIET = False

def cleanUpKeyValuePair(rawKey, rawValueStr):
  """
  Given a raw key and a raw value string, clean them up and
  remove a leading "-" if necessary.
  """
  paramStr = rawValueStr.strip()
  key = rawKey.strip().strip('-')
  return key, paramStr


def valToBool(rawVal):
  """
  Convert string or integer value from the config file or command line to a boolean.
  """
  retVal = False

  if type(rawVal) == str:
    retVal = (rawVal.strip().strip('"').lower()[0] in ['y','t','1'])

  else:
    try:
      retVal = bool(rawVal)
    except:
      retVal = False

  return retVal


def strToList(rawStr, defaultList):
  """
  Convert string value form the config file or command line to a list of
  items the same type as the default list.
  """
  # Get the type of the items in the default list.  If the list is empty,
  # assume the type is a string.
  listItemType = str
  if (len(defaultList) >= 1):
    listItemType = type(defaultList[0])
    
  strippedString = rawStr.strip().strip("]").strip("[").strip("}").strip("{")
  stringItems = re.split("[,;]",strippedString)

  returnList = []
  for item in stringItems:
    try:
      returnList.append( listItemType(item.strip()) )
    except:
      print "Warning:  Config Reader could not parse list item", item, \
            "because it is not of expected type", listItemType.__name__, \
            ".  Using default list instead."
      returnList = defaultList
      break

  return returnList

  
def getKeyedCMDLineParams(defaultDict):
  """
  Given a defaults dictionary, return a new dictionary of all
  command-line key=value pairs where the key matches a key in the
  defaults dictionary.  The fields are type cast to the same types as
  in the defaults dictionary.
  """
  cmdLineParams = {}
  
  for param in sys.argv:
    items = param.split("=")
    
    if (len(items)==2):
      key, paramStr = cleanUpKeyValuePair(items[0],items[1])
      if (key in defaultDict):
        defaultVal = defaultDict[key]
        try:
          if (type(defaultVal) == bool):
            cmdLineParams[key] = valToBool(paramStr)
          elif (type(defaultVal) == list):
            cmdLineParams[key] = strToList(paramStr,defaultVal)
          else:
            cmdLineParams[key] = type(defaultVal)(paramStr)
        except:
          pass
      
  return cmdLineParams


def getKeyedIniParams(filename, sectionname, defaultDict):
  """
  Given a defaults dictionary, return a new dictionary of all
  key=value pairs from the specified ini file and section where the
  key matches a key in the defaults dictionary.  The fields are type
  cast to the same types as in the defaults dictionary.
  """
  parser = None
  sections = []
  iniParams = {}

  if not os.path.isfile(filename):
    if (not filename[0] == '-') and (not GLOBAL_QUIET):
      print "Warning:  File", filename, "does not appear to exist."
      print
    return iniParams
  
  try:
    parser = SafeConfigParser()
    parser.read(filename)
    sections = parser.sections()
  except:
    print "Warning: Could not read config file", filename
    print
    return iniParams

  if not sectionname in sections:
    print "Warning:  Config file", filename, "contains no section called", sectionname
    print "Sections:", sections
    print
    return iniParams

  for key in defaultDict:
    lcKey = key.lower().strip()
    defaultVal = defaultDict[key]
    try:
      iniVal = parser.get(sectionname,lcKey)
      if (type(defaultVal) == bool):
        iniParams[key] = valToBool(iniVal)
      elif (type(defaultVal) == list):
        iniParams[key] = strToList(iniVal,defaultVal)
      else:
        iniParams[key] = type(defaultVal)(iniVal)
      
    except:
      pass
    
  return iniParams
  

def checkForHelpParam(defaultDict):
  helpFlag = False
  
  for cmdParam in sys.argv:
    if (cmdParam.strip().strip('-').lower()[0:4] == "help"):
      helpFlag =True
      print
      print "Below is a list of the parameters for this program and their default values."
      print "Parameters can be set in a configuration file or at the command line."
      for key in defaultDict:
        print "  ",key,"=", defaultDict[key]
      print
      
  return helpFlag

    
def buildArgObject(filename, sectionname, defaultDict, quiet=False):
  """
  Return an object with a member variable corresponding to every key in the
  defaultDict.  If there is a command-line parameter override, then the value
  comes from the command-line.  Otherwise, if there is a config-file override,
  then the value comes from the config file.  Otherwise the value of the member
  variable is the default value passed into the routine. NOTE:  The default
  values provide all type information, so DO NOT USE None as a default value.
  """
  requiredFailureFlag = False
  
  if checkForHelpParam(defaultDict):
    sys.exit(1)
  
  cmdLineParamsDict = getKeyedCMDLineParams(defaultDict)
  iniParamsDict     = {}
  if (not filename.strip() == ""):
    iniParamsDict = getKeyedIniParams(filename,sectionname,defaultDict)

  if (not quiet) and (len(iniParamsDict) > 0):
    print "Loaded some config params from config file", filename, "from section", sectionname

  argObject = empty()
  outputLines = []
  for key in defaultDict:
    val = defaultDict[key]
    suffix = "(defaulted)"
   
    if (key in cmdLineParamsDict):
      val = cmdLineParamsDict[key]
      suffix = "(from command line)"
      outputLines.insert(0, "  " + (key + ":").ljust(20) + " " + str(val) + " " + suffix)

    elif (key in iniParamsDict):
      val = iniParamsDict[key]
      suffix = "(from config file)"
      outputLines.insert(0, "  " + (key + ":").ljust(20) + " " + str(val) + " " + suffix)
 
    elif (type(val) == str) and (val.strip()[0:13] == "PROMPT_STRING"):
      promptStr = val.split(":")[1]
      sys.stdout.write(promptStr + ": ")
      val = raw_input().strip()
      suffix = "(prompt)"
      if val.strip() == '':
        requiredFailureFalg = True
      outputLines.append ("  " + (key + ":").ljust(20) + " " + str(val) + " " + suffix)

    else:
      outputLines.insert(0, "  " + (key + ":").ljust(20) + " " + str(val) + " " + suffix)
     
    setattr(argObject,key,val)
    if (type(val) == str) and (val.strip()[0:8] == "REQUIRED"):
      requiredFailureFlag = True

  if (not quiet):
    print
    print "Configuration parameters:"
    for line in outputLines:
      print line
    print

  if requiredFailureFlag:
    print
    print "ERROR:  At least one required parameter was not set."
    print
    sys.exit(1)
    
  return argObject

 
if __name__ == '__main__':
  foo = buildArgObject('rpwtest.ini','foobar',{"MyString":"this is a string", "MyInt":2, "MyFloat":-1.0, "MyOtherFloat":-10.1, "MyBool":False, "MyList":[1,2,3]})
