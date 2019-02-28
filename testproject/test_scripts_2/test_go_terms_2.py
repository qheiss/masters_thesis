import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt

glist = ['IGKV4-1', 'CD55', 'IGKC', 'PPFIBP1', 'ABHD4', 'PCSK6', 'PGD', 'ARHGDIB', 'ITGB2', 'CARD6']

enr = gp.enrichr(gene_list=glist,
                 description='test_name',
                 # gene_sets='KEGG_2016',
                 # or gene_sets='KEGG_2016,KEGG_2013',
                 gene_sets=['KEGG_2016','KEGG_2013'],
                 outdir='test/enrichr_kegg',
                 cutoff=0.5 # test dataset, use lower value of range(0,1)
                )

