#!/usr/bin/env python
"""
From https://raw.github.com/PacificBiosciences/HBAR-DTK/master/src/CA_best_edge_to_GML.py

Copyright (c) 2011-2013, Pacific Biosciences of California, Inc.

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted (subject to the limitations in the
disclaimer below) provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above
   copyright notice, this list of conditions and the following
   disclaimer in the documentation and/or other materials provided
   with the distribution.

 * Neither the name of Pacific Biosciences nor the names of its
   contributors may be used to endorse or promote products derived
   from this software without specific prior written permission.

NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
SUCH DAMAGE.

"""

"""
A simple script to convert Celera(R) Assembler's "best.edges" to a GML which can be used to feed into Gephi to
check the topology of the best overlapping graph.

Usage:
    python CA_best_edge_to_GML.py asm.gkp_store asm.tigStore best.edge output.gml 
"""

import networkx as nx
import os
import shlex
import sys
import subprocess

if "check_output" not in dir( subprocess ): # duck punch it in!
    def f(*popenargs, **kwargs):
        if 'stdout' in kwargs:
            raise ValueError('stdout argument not allowed, it will be overridden.')
        process = subprocess.Popen(stdout=subprocess.PIPE, *popenargs, **kwargs)
        output, unused_err = process.communicate()
        retcode = process.poll()
        if retcode:
            cmd = kwargs.get("args")
            if cmd is None:
                cmd = popenargs[0]
            raise subprocess.CalledProcessError(retcode, cmd)
        return output
    subprocess.check_output = f

gkp_store = sys.argv[1]
tig_store = sys.argv[2]
best_edge_file = sys.argv[3]
output_file = sys.argv[4]



G=nx.DiGraph()
frg_to_tig = {}
args = shlex.split("tigStore -g %s -t %s 1 -D unitiglist" % (gkp_store, tig_store ))
out = subprocess.check_output(args)
out = out.split("\n")
for l in out:
    l = l.strip().split()
    if len(l) == 0: continue
    if l[0] == "maID": continue
    unitig_id = int(l[0])

    os.system("tigStore -g %s -t %s 1 -d frags -u %d > frag_list" % ( gkp_store, tig_store, unitig_id) )

    args = shlex.split( "tigStore -g %s -t %s 1 -d frags -u %d" % ( gkp_store, tig_store, unitig_id) )
    f_out = subprocess.check_output(args)
    f_out = f_out.split("\n")
    for l in f_out:
        """FRG    1453 179419,182165"""
        l = l.replace(",", " ")
        l = l.strip().split()
        if len(l) == 0: continue
        frg_id = l[1]
        frg_to_tig[frg_id] = unitig_id

with open(best_edge_file) as f:
    for l in f:
        if l[0] == "#": continue
        l = l.strip().split()
        id1, lib_id, best5, o1, best3, o3 = l[:6]
        G.add_node(id1, unitig="utg%s" % frg_to_tig[id1])
        if best5 != "0":
            G.add_edge(best5, id1)
        if best3 != "0":
            G.add_edge(id1, best3)
        #G[id1]["unitig"] = frg_to_tig[id1]

nx.write_gml(G, output_file)