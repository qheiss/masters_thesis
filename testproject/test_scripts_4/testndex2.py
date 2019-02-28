from ndex2.nice_cx_network import NiceCXNetwork
import ndex2.client as nc
import ndex2
import networkx as nx
import pandas as pd
import os



#nice_cx_network = ndex2.create_nice_cx_from_server(server='public.ndexbio.org', uuid='becec556-86d4-11e7-a10d-0ac135e8bacf')
#nice_cx_network = ndex2.create_nice_cx_from_server(server='public.ndexbio.org', uuid='d14db454-3d18-11e8-a935-0ac135e8bacf')
nice_cx_network = ndex2.create_nice_cx_from_server(server='public.ndexbio.org', uuid='9c38ce6e-c564-11e8-aaa6-0ac135e8bacf')
#nice_cx_network = ndex2.create_nice_cx_from_server(server='public.ndexbio.org', uuid='1093e665-86da-11e7-a10d-0ac135e8bacf')

for node_id, node in nice_cx_network.get_nodes():
	print(node_id)
	print(node.get('n'))
	print(node.get('r'))
	print(nice_cx_network.get_node_attribute_value(node_id,'GeneName_A'))
#print("edges")

#for edge_id, edge in nice_cx_network.get_edges():
	#print(edge.get('s'))
	#print(edge.get('t'))
