#!/usr/bin/env python

#################################################################################$$
# Copyright (c) 2011-2014, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted (subject to the limitations in the
# disclaimer below) provided that the following conditions are met:
#
#  * Redistributions of source code must retain the above copyright
#  notice, this list of conditions and the following disclaimer.
#
#  * Redistributions in binary form must reproduce the above
#  copyright notice, this list of conditions and the following
#  disclaimer in the documentation and/or other materials provided
#  with the distribution.
#
#  * Neither the name of Pacific Biosciences nor the names of its
#  contributors may be used to endorse or promote products derived
#  from this software without specific prior written permission.
#
# NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
# GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
# BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
# USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
# OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
# SUCH DAMAGE.
#################################################################################$$

# from smrtanalysis-2.1.1.128549

""" Description of the functionality of this script """

import sys, os
import optparse, logging
from pbpy.io.FastaIO import SimpleFastqReader as sfqr, qvCharToInt
import numpy

class ScriptToRun:
    def __init__( self ):
        self.__parseArgs( )
        self.__initLog( )

    def __parseArgs( self ):
        """Handle command line argument parsing"""
        
        usage = "%prog [--help] [options] ARG1"
        parser = optparse.OptionParser( usage=usage, description=__doc__ )

        parser.add_option( "-l", "--logFile", help="Specify a file to log to. Defaults to stderr." )
        parser.add_option( "-d", "--debug", action="store_true", help="Increases verbosity of logging" )
        parser.add_option( "-i", "--info", action="store_true", help="Display informative log entries" )
        parser.add_option( "-p", "--profile", action="store_true", help="Profile this script, dumping to <scriptname>.profile" )
        parser.add_option( "--qvCut",  help="quality value cutoff")
        parser.add_option( "--trimFront",  help="beginning base pairs to remove(gmin)")
        parser.add_option( "--out",  help="out file(gdefaults to stdout)")
        parser.add_option( "--fastaOut",  help="out fasta file(gdefaults to None)")
        parser.add_option( "--minSeqLen",  help="min seq len to keep")

        parser.set_defaults(logFile=None, debug=False, info=False, profile=False, 
                            qvCut=10.0, out=None, trimFront=0, fastaOut=None, minSeqLen=50)
        
        self.opts, args = parser.parse_args( )

        if len(args) != 1:
            parser.error("Expected a single argument.")

        self.inFN = args[0]
        self.qvCut = float(self.opts.qvCut)
        self.trimFront = int(self.opts.trimFront)
        self.out = sys.stdout
        if self.opts.out:
            self.out = open(self.opts.out, 'w')
        if self.opts.fastaOut:
            self.fastaOut = open(self.opts.fastaOut, 'w')
        self.minSeqLen = int(self.opts.minSeqLen)

    def _trimFastq(self, qvCut):
        for entry in sfqr(self.inFN):
            qArray = map(lambda x: ord(x) -33 - qvCut, entry.quality)
            
            maxSum = 0
            maxStartIndex = 0
            maxEndIndex = 0
            currSum = qArray[0]
            currStartIndex = 0

            for i in range(1, len(qArray)):
                if currSum > 0:
                    currSum = currSum + qArray[i]
                else:
                    currSum = qArray[i]
                    currStartIndex = i
                if currSum > maxSum:                    
                    maxSum = currSum
                    maxStartIndex = currStartIndex
                    maxEndIndex = i

            if maxEndIndex-maxStartIndex < self.minSeqLen:
                continue
            maxStartIndex = max(maxStartIndex, self.trimFront)
            newName = "%s/%i_%i" %(entry.name, maxStartIndex, maxEndIndex)
            self.out.write("@%s\n%s\n+%s\n%s\n" 
                   % (newName, entry.sequence[maxStartIndex:maxEndIndex+1], 
                      newName, entry.quality[maxStartIndex:maxEndIndex+1]))
            if self.opts.fastaOut:
                self.fastaOut.write(">%s\n%s\n" 
                   % (newName, entry.sequence[maxStartIndex:maxEndIndex+1]))

    def __initLog( self ):
        """Sets up logging based on command line arguments. Allows for three levels of logging:
        logging.error( ): always emitted
        logging.info( ): emitted with --info or --debug
        logging.debug( ): only with --debug"""

        logLevel = logging.DEBUG if self.opts.debug else logging.INFO if self.opts.info else logging.ERROR
        logFormat = "%(asctime)s [%(levelname)s] %(message)s"
        if self.opts.logFile != None:
            logging.basicConfig( filename=self.opts.logFile, level=logLevel, format=logFormat )
        else:
            logging.basicConfig( stream=sys.stderr, level=logLevel, format=logFormat )
                                                                 
    def run( self ):
        """Executes the body of the script."""
    
        logging.info("Log level set to INFO")
        logging.debug("Log Level set to DEBUG")

        self._trimFastq(self.qvCut)

        return 0

if __name__ == "__main__":
    app = ScriptToRun()
    if app.opts.profile:
        import cProfile
        cProfile.run( 'app.run()', '%s.profile' % sys.argv[0] )
    sys.exit( app.run() )
