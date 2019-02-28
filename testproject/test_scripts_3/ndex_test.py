import json
from pprint import pprint
from os import listdir, makedirs, path
from os.path import isfile, isdir, join, abspath, dirname, exists, basename, splitext
import string
from multiprocessing import Pool
import time
import pandas as pd
import numpy as np
import itertools
import networkx as nx
import matplotlib.pyplot as plt
import json
import imp
import mygene
import os
#import datetime
import numpy as np
import matplotlib.pyplot as plt
import mpld3

import seaborn as sns
import pandas as pd
from numpy import array

from pybiomart import Dataset
import ndex.client as nc

from pybiomart import Dataset

fh = open("biogrid1.cx",encoding='utf-8')
fn = fh.read()
lines = fn.split("\n")
lines[3].replace("@","")
lines[3].replace("uniprot:","uniprot")
lines[3].replace("signor:","signor")
lines[3].replace(" ","")
lines[3].replace("ncbigene:","")
lines[3].replace("\\n","")
#lines[3].replace("
node_line = lines[3][11:-3].replace("ncbigene:","")
print(node_line)
nodelinesplit = node_line.split(", ")
#print(nodelinesplit)
dictlist = []
node_dict = {}
node_line_2 = "[" + node_line + "]"
#print(nodelinesplit[len(nodelinesplit)-1])
tmp2 = json.loads(node_line_2)
print(tmp2)
for item in tmp2:
	print(item)
	#tmp = json.loads(item)
	dictlist.append(item)
	if not(any(c.islower() for c in item['r'])):
		node_dict[item['@id']] = item['r']
		#print(tmp['@id'])
#print(dictlist)
lines[4].replace("@","")
lines[4].replace("uniprot:","uniprot")
lines[4].replace("signor:","signor")
lines[4].replace(" ","")
lines[4].replace("\\n","")
edge_line = lines[4][11:-3].rstrip()
edge_line_2 = "[" + edge_line + "]"
print(edge_line_2)
edgelinesplit = edge_line.split(", ")	
#print(edgelinesplit)
edgelist = []
tmp4 = json.loads(edge_line_2)
#dataset = Dataset(name='hsapiens_gene_ensembl',
#	                  host='http://www.ensembl.org')
#conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name','entrezgene'])
ret = []
for item in tmp4:
	print(item)
	#dictlist.append(tmp)
	if(item['s'] in node_dict and item['t'] in node_dict):
		source = node_dict[item['s']]
		target = node_dict[item['t']]
		if(source != target):
			baz = [str(source),str(target)]
			ret.append("\t".join(baz))
		#print(source + "\t" + target)
		#prot1_nbr = conv.index[conv['Gene name'] == source]
		#prot2_nbr = conv.index[conv['Gene name'] == target]
		#if(source in conv['Gene name'].tolist() and target in conv['Gene name'].tolist()):
#			prot1 = conv.loc[prot1_nbr,'NCBI gene ID'].values[0]
#			prot2 = conv.loc[prot2_nbr,'NCBI gene ID'].values[0]
#			if(prot1 != prot2):
#				baz = [str(int(prot1)), str(int(prot2))]
#				ret.append("\t".join(baz))
print("\n".join(ret))
