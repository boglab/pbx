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
# updated imports to find Amos-related things in pbaha instead of pbpy

import sys, logging, optparse 
from pbaha.io.ReadMatcherIO import parseRm4
from pbaha.io.AmosBank import AmosBank
from pbaha.model.AlignmentHit import AlignmentHit
from itertools import ifilter, groupby
from copy import deepcopy
from pbaha.model.Range import Range, Ranges

"""Converts a blasr4 alignment file to an AMOS LAY message file containing
   query reads ordered by the alignments to the target reads."""

class AlignToLayouts:

    def __init__( self ):
        self.__parseArgs( )

    def __parseArgs( self ):
        usage = "Usage: %prog in.alignb4 in.bank > out.lay"
        parser = optparse.OptionParser(usage=usage, description=__doc__)

        parser.add_option("--debug", action="store_true", 
                help="Outputs a log to stderr with helpful debug info.")
        parser.add_option("--logFile", help="Log file")
        parser.add_option("--overlapTolerance", 
                help="Amount of overlap to put reads into same layout.")
        parser.add_option("--mappingFile", 
                help="File to store mapping of corrected reads to contig ids.")
        parser.add_option("--includeTarget", action="store_true",
                help="Include targets in the layout. Note must be in bank.")
        parser.add_option("--mandatoryLayoutIdsFile", 
                help="Include each of these ids (bank EIDs) in a layout even if"
                     " they do not show up in any alignment.")
        parser.add_option("--restrictedRangesFile", 
                help="A file of ranges to restrict considered hits.  "
                     "Each line is a range with format: <target> <start1> <end1> . . .")
        parser.add_option("--readCoords", action="store_true",
                help="Input file has coords reported on the read rather than "
                     "subread (for blasr run on {pls|bas}.h5 files only).")
        parser.add_option("--allowPartialAlignments", action="store_true",
                help="Don't require reads to be fully aligned or full "
                     "overlap to participate in a layout.")

        parser.add_option("--trimHit", 
                help="Trim each hit by this amount.")

        parser.set_defaults(overlapTolerance=50, 
                            mappingFile=None,
                            includeTarget=False, 
                            mandatoryLayoutIdsFile=False, 
                            restrictedRangesFile=None,
                            trimHit=False,
                            readCoords=False,
                            allowPartialAlignments=False)
                            
        (self.opts, args) = parser.parse_args()
        self._configureLogging()
        
        if len(args) != 2: parser.error("Expected two arguments.")

        self.inAlignFile = args[0]
        self.inBank      = args[1]
        self.endTolerance = 20
        self.outMappingFile = self.opts.mappingFile if self.opts.mappingFile \
                                else "%s.mapping" % self.inAlignFile
        self.trimHit = int(self.opts.trimHit) if self.opts.trimHit else False
                
        self.overlapTolerance = int(self.opts.overlapTolerance)
        if self.opts.includeTarget:
            logging.info(
                "Including target in layouts. Overlap tol will be ignored")
        self.layoutCount = 0
        self.readCoords = self.opts.readCoords

    def _configureLogging( self ):
        """Does basic logging config if debug is True."""

        level = logging.DEBUG if self.opts.debug else logging.INFO
        format = "%(asctime)s [%(levelname)s] %(message)s" 
        if self.opts.logFile is not None:
            logging.basicConfig( 
                filename=self.opts.logFile, level=level, format=format)
        else: 
            logging.basicConfig( 
                stream=sys.stderr, level=level, format=format)

    def groupOverlappingHits(self, hits):
        """Breaks a list of hits into sublists that overlap one another.
        Returns: A list of lists of AlignmentHit objects"""

        hits.sort(key=lambda x: x.target_start)
        lastEnd = hits[0].target_start + len(hits[0])
        hitGroup = [ hits[0] ]
        hitGroups = []
        for h in hits[1:]:
            if (self.opts.includeTarget or 
                h.target_start  < lastEnd - self.overlapTolerance): 
                hitGroup.append(h)
            else:
                hitGroups.append(hitGroup)
                hitGroup = [h]
            # TODO perhaps len should be replaced with target_end - target_start ?
            lastEnd = max(lastEnd, h.target_start + len(h)) 
        if len(hitGroup): hitGroups.append(hitGroup)

        return hitGroups

    def _makeTargetSelfHit(self, targetId, bank):
        """Adds a hit between the target and itself to hitGroup, so that the target
        will ultimately be included in the layout."""

        hit = AlignmentHit()
        hit.query_id = hit.target_id = targetId
        hit.query_strand = hit.target_strand = "+"
        length = bank.getRead(targetId).getLength() 
        hit.query_end = hit.query_length = hit.target_end = hit.target_length = length
        hit.query_start = hit.target_start = 0 
        return hit

    def groupToLayout(self, hitGroup, eid2iid, mappingFile):
        """Converts a list of AlignmentHit objects in hitGroup to an AMOS layout 
        message, using the EID to IID mapping in eid2iid. The mappingFile keeps 
        track of a mapping between layout IIDs and hit targets. 

        Returns: A string containing the AMOS layout message."""

        # TODO take mappingFile out of this routine?
        self.layoutCount += 1
        lastEnd = hitGroup[0].target_start + len(hitGroup[0])
        for h in hitGroup:
            lastEnd = max(lastEnd, h.target_start + len(h))

        message = ["{LAY\niid:%i\n" % self.layoutCount]
        targetStart = hitGroup[0].target_start
        for h in hitGroup:
            iid = eid2iid[h.query_id]
            message.append( "{TLE\n" )
            if h.query_strand == "+":
                message.append( "clr:%i,%i\n" % (h.query_start, h.query_end) )
            else:
                message.append( "clr:%i,%i\n" % (h.query_end, h.query_start) )
            message.append( "off:%i\n" % (h.target_start - targetStart) )
            message.append( "src:%i\n" % iid )
            message.append( "}\n" )
        message.append( "}\n" )
        mappingFile.write("contig%i,%s\n" \
                          % (self.layoutCount, hitGroup[0].target_id))
        return "".join( message )

    def _satisfiesHitCriteria(self, hit, targetRanges=None):
        # check if the hit overlaps restricted ranges
        if targetRanges is not None:
            hitRange = Range(hit.target_start, hit.target_end)
            if hit.target_id in targetRanges:
                restrictedRanges = targetRanges[hit.target_id]
                if not restrictedRanges.intersects(hitRange): 
                    return False
            else:
                return False
        # check if hit is fully aligned
        return self.opts.allowPartialAlignments or hit.fullyAligned() or hit.hasFullOverlap()

    def _restrictedRanges(self):
        """Returns a dict mapping target id to a Ranges() object. 
           Returns None if no restricted Ranges."""
        if self.opts.restrictedRangesFile is None: return None
        targetRanges = {}
        rf = open (self.opts.restrictedRangesFile)
        for line in rf.xreadlines():
            lineList = line.strip().split()
            target = lineList[0]
            start, end = map(int, lineList[1:])
            targetRanges.setdefault(target, Ranges()).addRange(Range(start, end))
        return targetRanges

    def _mandatoryLayoutIds(self):
        """Returns: A set of ids that must be put into a layout."""
        if self.opts.mandatoryLayoutIdsFile: 
            return set([ l.rstrip() for l in open(
                            self.opts.mandatoryLayoutIdsFile).readlines() ])
        else:
            return set([])

    def _trimHits(self, hits, trimHit, minHitLength=30):
        newHits = []
        for h in hits:
            h.query_start  += trimHit
            h.target_start += trimHit
            h.query_end    -= trimHit
            h.target_end   -= trimHit
            if (h.query_end - h.query_start > minHitLength and 
                h.target_end - h.target_start > minHitLength):
                newHits.append(h)
        return newHits

    def run(self):
        layoutsOut = sys.stdout
        targetRanges = self._restrictedRanges()
        hits = ifilter(lambda x: self._satisfiesHitCriteria(x, targetRanges), 
                                 parseRm4(self.inAlignFile, 
                                 readCoords=self.readCoords))
        bank = AmosBank( self.inBank )
        hitsInBank = ifilter(lambda x: inBank(x, bank), hits)
        
        didOutputGroup = False

        mandatoryIds = self._mandatoryLayoutIds() 
        eid2iid = {r.getEID():r.getIID() for r in bank.readIterator()}
        mappingFile = open(self.outMappingFile, "w")
        for targetId, hitIterator in groupby(hitsInBank, lambda x: x.target_id):
            hitList = list(hitIterator)
            if self.trimHit:
                hitList = self._trimHits(hitList, self.trimHit)
            if not hitList: continue
            if self.opts.includeTarget: 
                hitList.insert(0, self._makeTargetSelfHit(targetId, bank))
            hitGroups = self.groupOverlappingHits(hitList)
            for g in hitGroups: 
                print >>layoutsOut, self.groupToLayout(g, eid2iid, mappingFile),
                didOutputGroup = True
            if targetId in mandatoryIds: mandatoryIds.remove(targetId)

        for targetId in mandatoryIds:
            hitGroup = [ self._makeTargetSelfHit(targetId, bank) ]
            print >>layoutsOut, self.groupToLayout(hitGroup, eid2iid, mappingFile),

        layoutsOut.close()
        mappingFile.close()

        if not didOutputGroup: 
            logging.fatal("Failed to output any layouts!")
            return 1
        return 0

def inBank(hit, bank):
    """Is the hit query_id in the bank"""
    retVal = True
    try:
        bank.getReadIIDfromEID(hit.query_id)
    except:
        retVal = False
    return retVal

if __name__ == "__main__":
    
    sys.exit( AlignToLayouts().run() )
