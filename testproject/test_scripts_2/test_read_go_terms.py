import mygene
from pybiomart import Dataset
mg = mygene.MyGeneInfo()
import math

#fh1 = open("test_gene_list_2.txt")
#lines = fh1.readlines()
#for line in lines:
#	print(line.split("\t")[1])
#	print(mg.getgene(line.split("\t")[1], fields='all'))

dataset = Dataset(name='hsapiens_gene_ensembl',
		                  host='http://www.ensembl.org')
conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name','entrezgene'])

convlist = conv['Gene name'].tolist()


text_file = open("test_go_1.txt","w")

fh1 = open("test_go_list.txt")
lines = fh1.readlines()
term_dict = {}
ctr = 0
for line in lines:
	if(ctr == 50000):
		print(ctr)
		print(term_dict)
		ctr = 0
		for key in term_dict:
			if(key !=float('nan')):
				text_file.write(str(int(key)) + "\t" + ";".join(term_dict[key]))
				text_file.write("\n")
	geneid=line.split("\t")[2]
	if(geneid in convlist):
		tmp1 = conv.index[conv['Gene name'] == geneid]
		tmp2 = conv.loc[tmp1,'NCBI gene ID'].values[0]
		goterm=line.split("\t")[4]
		#print(goterm)
		if(geneid in term_dict):
			term_dict[tmp2].append(goterm)
		else:
			term_dict[tmp2] = []
			term_dict[tmp2].append(goterm)
	ctr = ctr + 1

print(term_dict)


for key in term_dict:
	text_file.write(str(int(key)) + "\t" + ";".join(term_dict[key]))
	text_file.write("\n")
text_file.close()
