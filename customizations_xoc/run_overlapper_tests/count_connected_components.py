import networkx as nx
import sys

G = nx.read_gml(sys.argv[1])

concomps = nx.weakly_connected_components(G)
concompslen = ",".join([str(len(concomp)) for concomp in concomps])
print("%d: %s" % (len(concomps), concompslen))
