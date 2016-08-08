##############################################################################################
# Modified csv2json.py in order to parse csv files avoiding a certain prescale = n HLT paths #
##############################################################################################
#!/usr/bin/env python

import sys
import optparse
import re
from FWCore.PythonUtilities.LumiList import LumiList


if __name__ == '__main__':
    
    parser = optparse.OptionParser ("Usage: %prog input.csv")
    parser.add_option ('--output', dest='output', type='string',
                       help='Save output to file OUTPUT')
    parser.add_option ('--runIndex', dest='runIndex', type='int',
                       default = 0,
                       help='column to be converted to run number (default %default)')
    parser.add_option ('--lumiIndex', dest='lumiIndex', type='int',
                       default = 2,
                       help='column to be converted to lumi section number (default %default)')
    parser.add_option ('--HLTprescIndex', dest='HLTprescIndex', type='int',
                       default = 5,
                       help='column to be converted to HLT prescale number (default %default)')
    # required parameters
    (options, args) = parser.parse_args()
    if len (args) != 1:
        raise RuntimeError, "Must provide exactly one input file"

    sepRE = re.compile (r'[\s,;:]+')
    runLumiDict = {}    
    events = open (args[0], 'r')
    runIndex, lumiIndex, HLTprescIndex = options.runIndex, options.lumiIndex, options.HLTprescIndex
    minPieces = max (runIndex, lumiIndex, HLTprescIndex) + 1
    for line in events:
        pieces = sepRE.split (line.strip())
        if len (pieces) < minPieces:
            continue
        try:
            run, lumi, HLTpresc = int( pieces[runIndex] ), int( pieces[lumiIndex] ), int( pieces[HLTprescIndex] )
        except:
            continue
        if HLTpresc == 1: # Set the prescale you want to avoid
            continue
        runLumiDict.setdefault (run, []).append (lumi)
    jsonList = LumiList (runsAndLumis = runLumiDict) 
    if options.output:
        jsonList.writeJSON (options.output)
    else:
        print jsonList
