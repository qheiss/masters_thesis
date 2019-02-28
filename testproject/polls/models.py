import datetime
from django.db import models
from django.utils import timezone
# Create your models here.
import os, sys
from django.shortcuts import render

# Create your views here.
from django.http import HttpResponse
from django.template import loader
from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from django.shortcuts import get_object_or_404, render
from django.urls import reverse
from django.views import generic
from scipy.stats.stats import pearsonr
import plotly
import plotly.plotly as py
import plotly.graph_objs as go
import plotly.io as pio
import plotly.offline
#from .models import Choice, Question
from datetime import datetime
from networkx.readwrite import json_graph
import json
from django.shortcuts import render_to_response,render
from django.template import RequestContext
from django.http import HttpResponseRedirect
from django.urls import reverse
import polls
#from polls.models import Document
from polls.forms import DocumentForm
#from polls.models import Upload,UploadForm
import numpy as np
import matplotlib.pyplot as plt
import mpld3

import seaborn as sns
import pandas as pd
from numpy import array

import matplotlib.patches as mpatches


import networkx as nx
from bokeh.io import show, output_notebook, output_file, save
from bokeh.plotting import figure
from bokeh.models import Circle, HoverTool, TapTool, BoxSelectTool
from bokeh.models.graphs import from_networkx
from bokeh.transform import linear_cmap
from bokeh.models import ColumnDataSource, LabelSet
from bokeh.models.graphs import NodesAndLinkedEdges, EdgesAndLinkedNodes
from biomart import BiomartServer
from bokeh.embed import components
from bokeh.palettes import Spectral4
from bokeh.models import Plot, Range1d, MultiLine, Circle, HoverTool, TapTool, BoxSelectTool

from pybiomart import Dataset

class Question(models.Model):
    question_text = models.CharField(max_length=200)
    pub_date = models.DateTimeField('date published')
    def __str__(self):
        return self.question_text
    def was_published_recently(self):
        return self.pub_date >= timezone.now() - datetime.timedelta(days=1)



class Choice(models.Model):
    question = models.ForeignKey(Question, on_delete=models.CASCADE)
    choice_text = models.CharField(max_length=200)
    votes = models.IntegerField(default=0)
    def __str__(self):
        return self.choice_text


class Document(models.Model):
#	docfile = models.FileField(upload_to='documents/%Y/%m/%d')
	docfile = models.FileField(upload_to='/home/quirin/bla')
	class Meta:
		app_label = 'Document'

class DocumentForm(models.Model):
	docfile = models.FileField(
	#label='Select a file',
	help_text='max. 42 megabytes'
	)
	class Meta:
		app_label = 'DocumentForm'

from django.forms import ModelForm

class Upload(models.Model):
        pic = models.FileField(upload_to="images/")
        upload_date=models.DateTimeField(auto_now_add =True)

# FileUpload form class.
class UploadForm(ModelForm):
        class Meta:
                model = Upload
                fields = ('pic',)

class GraphForm(models.Model):
	def handle_upload_3_4(fn,prot_fn):
		patients = []
		
		patients1 = []
		patients2 = []
		genes = []
		geneNames = []
		data = {}
		data1 = {}
		data2 = {}
		group1 = []
		group2 = []
		group_labels1 = []
		group_labels2 = []
		group1_data = []
		group2_data = []
		
		patients.append([1,2,3,4])
		patients1.append(['1','2'])
		#patients.append(['3','4'])
		patients2.append(['3','4'])
		group_labels1.append([1,1])
		group_labels2.append([2,2])
		logstring = "Creating Plots for given input files... \n\n"
		logstring = "Reading gene expression data... \n"
		ctr_bar = 0
		df2 = pd.read_csv(fn,delim_whitespace=True,header=None,index_col=0)
		#diff = {}
		ctr2 = 0
		#group1.append(group_labels1)
		#group2.append(group_labels2)
		survival_info =  df2.loc['SURVIVAL',:].apply(pd.to_numeric).tolist()
		#print(survival_info)
		for k in range(0,len(survival_info)):
			if(survival_info[k] == 1):
				group1.append(k)
			else:
				group2.append(k)
		#group1.append(df2.index[df2['SURVIVAL']==1].tolist())
		#group1.append(df2.index[df2['SURVIVAL']==0].tolist())
		genes = ()
		conv = pd.read_csv("Output4.txt",delim_whitespace=True,header=0,index_col=0)	
		genes_3 = {}
		genes_4 = {}
		logstring = logstring + "\n\n Matching gene and protein IDs... \n"
		#foobar_2 = list(df2.iloc[4:,0])
		foobar_2 = list(df2.index)[4:]
		#print(foobar_2)
		ensg_id_list = conv.iloc[:,0].tolist()
		conv_dict = conv.to_dict('index')
		conv_dict_2 = {}
		#print(ensg_id_list)
		df2shape = df2.shape[1]-1
		for i in foobar_2:
			#print(data1[1000])
			#here could be check whether level of expression is >100
			#if(float(df2.loc[i].tolist()[2]) > 100.0):
			ensg_id = i.split('.')[0]
			if(ensg_id in conv_dict):
				prot_id = conv_dict[ensg_id]['Gene_name']
				conv_dict_2[prot_id]=i
				genes_4.update({prot_id:0})
		G = nx.Graph()
		nodes = []
		#dataset = Dataset(name='hsapiens_gene_ensembl',
	        #          host='http://www.ensembl.org')
		#conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
		#print(conv.iloc[conv['Gene stable ID'] == 'ENSG00000223972']['Gene name'])
		#print(conv.index[conv['Gene stable ID'] == 'ENSG00000284332'].tolist())
		#dataset.list_attributes()
		genes = {}
		genes_3 = {}
		logstring = logstring + "\n\nReading PPI file...\n"
		for line in prot_fn.readlines():
			lineSplit = line.decode('utf8').split()
			prot1 = lineSplit[0]
			prot2 = lineSplit[1]
			#logstring = logstring + prot1 + prot2
			if(prot1 not in nodes and prot1 in genes_4):
				nodes.append(prot1)
				#logstring = logstring + str(genes_3[prot1])
				diff = 0
				diff_curr = 0
				gene_tmp = conv_dict_2[prot1]
				for j in range(0,df2shape):
						if(j in group1):
							#print(df2.loc[i])
							diff_curr = diff_curr + float(df2.loc[gene_tmp].tolist()[j])
						else:
							#print(df2.loc[i])
							diff_curr = diff_curr - float(df2.loc[gene_tmp].tolist()[j])
				diff_curr_2 = diff_curr / 2000.0
				#diff.update({i:diff_curr_2})
				G.add_node(prot1, Name=prot1, d=float(diff_curr_2))
				genes_3.update({prot_id:diff_curr_2})
				#G.add_node(prot1, Name=prot1, d=float(genes_3[prot1]))
				genes.update({prot1:diff_curr_2})
			if(prot2 not in nodes and prot2 in genes_4):
				nodes.append(prot2)
				#G.add_node(prot2, Name=prot2, d=float(genes_3[prot2]))
				#genes.update({prot2:genes_3[prot2]})
				gene_tmp = conv_dict_2[prot2]
				diff = 0
				diff_curr = 0
				for j in range(0,df2shape):
						if(j in group1):
							#print(df2.loc[i])
							diff_curr = diff_curr + float(df2.loc[gene_tmp].tolist()[j])
						else:
							#print(df2.loc[i])
							diff_curr = diff_curr - float(df2.loc[gene_tmp].tolist()[j])
				diff_curr_2 = diff_curr / 2000.0
				#diff.update({i:diff_curr_2})
				G.add_node(prot2, Name=prot2, d=float(diff_curr_2))
				genes_3.update({prot_id:diff_curr_2})
				#G.add_node(prot1, Name=prot1, d=float(genes_3[prot1]))
				genes.update({prot2:diff_curr_2})
			#logstring = logstring + str(genes_3[prot1]) + str(genes_3[prot2])
			#if(prot2 in genes_3 and prot1 in genes_3):
			G.add_edge(prot1,prot2)
		output_notebook()
		plot = figure(x_range=(-1.1, 1.1), y_range=(-1.1, 1.1))
		# add tools to the plot
		plot.add_tools(HoverTool(tooltips=[("Name", "@d")]), 
	               TapTool(), 
	               BoxSelectTool())
		# create bokeh graph
		graph = from_networkx(G, nx.spring_layout, iterations=1000, scale=1, center=(0,0))
		
		# add name to node data
		#graph.node_renderer.data_source.data['d'] = list(G.nodes())
		
		
		graph.node_renderer.data_source.data['d'] = [genes[i] for i in G.nodes]
		
		# add name to node data
		graph.node_renderer.data_source.data['Name'] = list(G.nodes())
		x,y = zip(*graph.layout_provider.graph_layout.values())
		#x = [i for i in G.nodes]
		#y = [genes[i] for i in G.nodes]
		node_labels = nx.get_node_attributes(G, 'Name')
		#logstring = logstring + str(node_labels)
		#source = ColumnDataSource({'x': x, 'y': y,
		 #                          'Name': [node_labels[i] for i in range(len(x))]})
		z = tuple([(bar-0.1) for bar in x])
		#print(x)
		#print(z)
		source = ColumnDataSource(data=dict(x=x, y=y, 
	                           Name=[node_labels[i] for i in genes.keys()]))
		#labels2 = LabelSet(x='z', y='y', text='Name', source=source,
                #  background_fill_color='white',level='glyph',render_mode='canvas')
		#labels2 = LabelSet(x='x', y='y', text='Name',x_offset=-30, source=source,
                #  background_fill_color='white',level='glyph',background_fill_alpha=0.0,render_mode='canvas')
		#print(labels2)
		#logstring = logstring + str(graph.layout_provider.graph_layout.values())
		# add club to node data
		#graph.node_renderer.data_source.data['club'] = [i[1]['club'] for i in G.nodes(data=True)]
		
		# set node size
		#graph.node_renderer.glyph = Circle(size=10, fill_color=linear_cmap('d', 'Spectral8', min(G.nodes()), max(G.nodes())))
		#graph.node_renderer.glyph = Circle(size=100, fill_color=linear_cmap('d', 'Spectral8', -2.5, 2.5))
		labels2 = LabelSet(x='x', y='y', text='Name',text_font_size="8pt",x_offset=-20, source=source,
                  background_fill_color='white',background_fill_alpha=0.0,level='glyph',render_mode='canvas')
		graph.node_renderer.glyph = Circle(size=60, fill_color=linear_cmap('d', 'Spectral8', -2.5, 2.5))
		
		plot.renderers.append(graph)
		plot.renderers.append(labels2)
		#graph.inspection_policy = EdgesAndLinkedNodes()
		plot.add_layout(labels2)
		#logstring = logstring + "\n\nPPI Graph created..."
		#output_file("polls/interactive_graphs.html")
		save(plot)
		#dataArray1 = np.array(group1)
		#dataArray2 = np.array(group2)
		#endData = np.concatenate((dataArray1, dataArray2), axis=1)
		#logstring = logstring + str(endData)
		red_patch = mpatches.Patch(color='red', label='Group1')
		blue_patch = mpatches.Patch(color='blue', label='Group2')
		#lut = dict(zip(set(endData[0]), sns.hls_palette(len(set(endData[0])), l=0.5, s=0.8)))
		#col_colors = pd.DataFrame(endData[0])[0].map(lut)
		#print(col_colors)
		colordict={0:'#BB0000',1:'#0000BB'}
		#colors = np.array(['#BB0000','#BB0000','#0000BB','#0000BB'])
		#df = pd.DataFrame(data=endData[1:,0:],index=geneNames,columns=patients)
		#df2 = pd.DataFrame(data=endData[0,0:], index='',columns=patients)
		#my_palette = dict(zip(df[.unique(), ["orange","yellow","brown"]))
		#row_colors = df2.cyl.map(my_palette)
		#fig, (ax1, ax2) = plt.subplots(1,2,sharex=True,sharey=True)
		#logstring = logstring + str(df)
		#df2 = pd.read_csv('testgenes13.txt',delim_whitespace=True,header=None,index_col=0)
		print(df2)
		#print(foobar)
		#print(colors)
		df5 = df2.sort_values(by='SURVIVAL',axis=1)
		foobar = df5.loc['SURVIVAL',:].apply(pd.to_numeric)
		colors = [colordict[i] for i in foobar]
		#df6 = df2.iloc[5:50000,1:].apply(pd.to_numeric)
		print(df5)
		colors = [colordict[i] for i in foobar]
		#print(colors)
		df6 = df5.iloc[5:50000,:].apply(pd.to_numeric)
		#print(df6)
		df7 = df6[df6.iloc[:,2]>100.0]
		#print(df7)
		print(df7)
		sns.clustermap(df7,method="average", figsize=(13, 13),col_colors=colors,col_cluster=False)
		#sns.clustermap(dataArray,method="average", figsize=(13, 13),ax=ax2)	
		#plt.imshow(dataArray, cmap='hot', interpolation='nearest')
		#plt.yticks(range(data_present_length), geneNames, size='small')
		#plt.xticks(range(len(patients)), patients, size='small')
		plt.legend(handles=[red_patch,blue_patch],loc=2,mode="expand",bbox_to_anchor=(3, 1))
		#plt.switch_backend('SVG') 
		plt.savefig("/home/quirin/testproject/polls/static/test.png")
		#print(logstring)
		script, div = components(plot)
		plot_1=plt.gcf()
		plt.close()
		return(script,div,plot_1)
	def handle_upload_3_survival_3(fn,prot_fn):
		patients = []
		
		patients1 = []
		patients2 = []
		genes = []
		geneNames = []
		data = {}
		data1 = {}
		data2 = {}
		group1 = []
		group2 = []
		group_labels1 = []
		group_labels2 = []
		group1_data = []
		group2_data = []
		
		patients.append([1,2,3,4])
		patients1.append(['1','2'])
		#patients.append(['3','4'])
		patients2.append(['3','4'])
		group_labels1.append([1,1])
		group_labels2.append([2,2])
		logstring = "Creating Plots for given input files... \n\n"
		logstring = "Reading gene expression data... \n"
		ctr_bar = 0		
		df2 = pd.read_csv(fn,delim_whitespace=True,header=None,index_col=0)
		diff = {}
		ctr2 = 0
		#group1.append(group_labels1)
		#group2.append(group_labels2)
		group1.append([1,1])
		group2.append([2,2])		
		survival_info =  df2.loc['SURVIVAL',:].apply(pd.to_numeric).tolist()
		survival_yrs =  df2.loc['SURVIVAL_YRS',:].apply(pd.to_numeric).tolist()
		print(survival_info)
		for k in range(0,len(survival_info)):
			if(survival_info[k] == 1):
				group1.append(k)
			else:
				group2.append(k)
		#group1.append(df2.index[df2['SURVIVAL']==1].tolist())
		#group1.append(df2.index[df2['SURVIVAL']==0].tolist())
		genes = ()
		dataset = Dataset(name='hsapiens_gene_ensembl',
	                  host='http://www.ensembl.org')
		conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
		genes_3 = {}
		genes_corr = {}
		logstring = logstring + "\n\n Matching gene and protein IDs... \n"
		#foobar_2 = list(df2.iloc[4:,0])
		foobar_2 = list(df2.index)[4:]
		print(foobar_2)
		ensg_id_list = conv['Gene stable ID'].tolist()
		for i in foobar_2:
			#print(data1[1000])
			#here could be check whether level of expression is >100
			if(float(df2.loc[i].tolist()[2]) > 100.0):
				ensg_id = i.split('.')[0]
				if(ensg_id in ensg_id_list):
				#if(len(data1[i]) == len(patients1) and len(data2[i]) == len(patients2)):
					#print(data1[1000])
					ctr2 = ctr2 + 1
					#print(i)
					#if(data1[i][0] > 10000.0):
					#	print(data1[i])
					geneNames.append(ctr2)
					#group1.append(data1[i])
					#group2.append(data2[i])
					#print(conv.iloc[conv['Gene stable ID'] == 'ENSG00000223972']['Gene name'])
					
					print(ensg_id)
					prot_nbr = conv.index[conv['Gene stable ID'] == ensg_id].tolist()
					print(prot_nbr)
					prot_id = conv.loc[prot_nbr,'Gene name'].values[0]
					#print(prot_id)
					diff_curr = 0
					for j in range(0,(df2.shape[1]-1)):
						if(j in group1):
							print(df2.loc[i])
							diff_curr = diff_curr + float(df2.loc[i].tolist()[j])
						else:
							print(df2.loc[i])
							diff_curr = diff_curr - float(df2.loc[i].tolist()[j])
					if(1==1):
						#group1.append(data1[i])
						#group2.append(data2[i])
						logstring = logstring + str(prot_id)
						diff_curr_2 = diff_curr / 2000.0
						diff.update({i:diff_curr_2})
						genes_3.update({prot_id:diff_curr_2})
						tmp_3 = df2.loc[i].tolist()
						tmp_4 = [float(x) for x in tmp_3]
						tmp_corr = pearsonr(tmp_4,survival_yrs)
						print(tmp_corr)
						genes_corr.update({prot_id:tmp_corr})
						#print(prot_id)
		
		G = nx.Graph()
		nodes = []
		dataset = Dataset(name='hsapiens_gene_ensembl',
	                  host='http://www.ensembl.org')
		conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
		#print(conv.iloc[conv['Gene stable ID'] == 'ENSG00000223972']['Gene name'])
		#print(conv.index[conv['Gene stable ID'] == 'ENSG00000284332'].tolist())
		dataset.list_attributes()
		genes = {}
		logstring = logstring + "\n\nReading PPI file...\n"
		for line in prot_fn.readlines():
			lineSplit = line.decode('utf8').split()
			prot1 = lineSplit[0]
			prot2 = lineSplit[1]
			logstring = logstring + prot1 + prot2
			if(prot1 not in nodes and prot1 in genes_3):
				nodes.append(prot1)
				logstring = logstring + str(genes_3[prot1])
				G.add_node(prot1, Name=prot1, d=float(genes_3[prot1]))
				genes.update({prot1:genes_3[prot1]})
			if(prot2 not in nodes and prot2 in genes_3):
				nodes.append(prot2)
				G.add_node(prot2, Name=prot2, d=float(genes_3[prot2]))
				genes.update({prot2:genes_3[prot2]})
			logstring = logstring + str(genes_3[prot1]) + str(genes_3[prot2])
			G.add_edge(prot1,prot2)
		#diff_frame = pd.DataFrame(list(genes.items()), columns=['Gene', 'Difference'])
		genelist_diff = []
		#diff_frame_2 = diff_frame.sort_values('Difference')
		#for i in range(0,len(diff_frame_2)):
		#	curr_tmp = diff_frame_2.iloc[i,0]
		#	print(curr_tmp)
		#	genelist_diff.append({'gene_name':curr_tmp,'difference':diff_frame_2.loc[i,'Difference']})
		genelist_ret =[]
		for k,v in list(genes.items()):
			#print(k)
			#print(v)
			genelist_diff.append({'gene_name':k,'difference':v})
			genelist_ret.append({'gene_name':k,'difference':v,'correlation':genes_corr[k]})
		output_notebook()
		plot = figure(x_range=(-1.1, 1.1), y_range=(-1.1, 1.1))
		# add tools to the plot
		plot.add_tools(HoverTool(tooltips=[("Name", "@d")]), 
	               TapTool(), 
	               BoxSelectTool())
		# create bokeh graph
		graph = from_networkx(G, nx.spring_layout, iterations=1000, scale=1, center=(0,0))
		
		# add name to node data
		#graph.node_renderer.data_source.data['d'] = list(G.nodes())
		
		
		graph.node_renderer.data_source.data['d'] = [genes[i] for i in G.nodes]
		
		# add name to node data
		graph.node_renderer.data_source.data['Name'] = list(G.nodes())
		x,y = zip(*graph.layout_provider.graph_layout.values())
		#x = [i for i in G.nodes]
		#y = [genes[i] for i in G.nodes]
		node_labels = nx.get_node_attributes(G, 'Name')
		logstring = logstring + str(node_labels)
		#source = ColumnDataSource({'x': x, 'y': y,
		 #                          'Name': [node_labels[i] for i in range(len(x))]})
		
		source = ColumnDataSource(data=dict(x=x, y=y, 
	                           Name=[node_labels[i] for i in genes.keys()]))
		#labels2 = LabelSet(x='x', y='y', text='Name', source=source,
                #  background_fill_color='white',background_fill_alpha=0.0,level='glyph',render_mode='canvas')
		#print(labels2)
		#logstring = logstring + str(graph.layout_provider.graph_layout.values())
		# add club to node data
		#graph.node_renderer.data_source.data['club'] = [i[1]['club'] for i in G.nodes(data=True)]
		
		# set node size
		#graph.node_renderer.glyph = Circle(size=10, fill_color=linear_cmap('d', 'Spectral8', min(G.nodes()), max(G.nodes())))
		#graph.node_renderer.glyph = Circle(size=100, fill_color=linear_cmap('d', 'Spectral8', -2.5, 2.5))
		labels2 = LabelSet(x='x', y='y', text='Name',text_font_size="8pt",x_offset=-20, source=source,
                  background_fill_color='white',background_fill_alpha=0.0,level='glyph',render_mode='canvas')
		graph.node_renderer.glyph = Circle(size=60, fill_color=linear_cmap('d', 'Spectral8', -2.5, 2.5))
		plot.renderers.append(graph)
		plot.renderers.append(labels2)
		#graph.inspection_policy = EdgesAndLinkedNodes()
		plot.add_layout(labels2)
		logstring = logstring + "\n\nPPI Graph created..."
		output_file("polls/interactive_graphs.html")
		save(plot)
		#dataArray1 = np.array(group1)
		#dataArray2 = np.array(group2)
		#endData = np.concatenate((dataArray1, dataArray2), axis=1)
		#logstring = logstring + str(endData)
		red_patch = mpatches.Patch(color='red', label='Group1')
		blue_patch = mpatches.Patch(color='blue', label='Group2')
		#lut = dict(zip(set(endData[0]), sns.hls_palette(len(set(endData[0])), l=0.5, s=0.8)))
		#col_colors = pd.DataFrame(endData[0])[0].map(lut)
		#print(col_colors)
		colordict={0:'#BB0000',1:'#0000BB'}
		#colors = np.array(['#BB0000','#BB0000','#0000BB','#0000BB'])
		#df = pd.DataFrame(data=endData[1:,0:],index=geneNames,columns=patients)
		#df2 = pd.DataFrame(data=endData[0,0:], index='',columns=patients)
		#my_palette = dict(zip(df[.unique(), ["orange","yellow","brown"]))
		#row_colors = df2.cyl.map(my_palette)
		#fig, (ax1, ax2) = plt.subplots(1,2,sharex=True,sharey=True)
		#logstring = logstring + str(df)
		#df2 = pd.read_csv(fn,delim_whitespace=True,header=None,index_col=0)
		#foobar = df2.loc['SURVIVAL',:].apply(pd.to_numeric)
		#print(foobar)
		#print(colors)			
		#df22 = pd.read_csv(fn,delim_whitespace=True,header=None,index_col=0)
		df6 = df2.iloc[5:50000,1:].apply(pd.to_numeric)
		df5 = df2.sort_values(by='SURVIVAL',axis=1)
		foobar = df5.loc['SURVIVAL',:].apply(pd.to_numeric)
		#colors = [colordict[i] for i in foobar]
		print(df5)
		colors = [colordict[i] for i in foobar]
		#print(colors)
		df6 = df5.iloc[5:50000,:].apply(pd.to_numeric)
		print(df6)
		df7 = df6[df6.iloc[:,2]>100.0]
		print(df7)
		sns.clustermap(df7,method="average", figsize=(13, 13),col_colors=colors,col_cluster=False)
		#sns.clustermap(dataArray,method="average", figsize=(13, 13),ax=ax2)	
		#plt.imshow(dataArray, cmap='hot', interpolation='nearest')
		#plt.yticks(range(data_present_length), geneNames, size='small')
		#plt.xticks(range(len(patients)), patients, size='small')
		plt.legend(handles=[red_patch,blue_patch],loc=2,mode="expand",bbox_to_anchor=(3, 1))
		#plt.switch_backend('SVG') 
		plt.savefig("/home/quirin/testproject/polls/static/test3.png")
		print(logstring)
		script, div = components(plot)
		plot_1=plt.gcf()
		survival_yrs = df2.loc['SURVIVAL_YRS',:].apply(pd.to_numeric)
		survival_1 = {0:1.0}
		survival_2 = {0:1.0}
		survivor_nbr = sum(foobar)
		death_nbr = len(foobar) - sum(foobar)
		count1 = 1.0
		count2 = 1.0
		for j in range(1,len(survival_yrs)):
			if(foobar[j] == 0):
				survival_1[j] = count1 - (1/death_nbr)
				count1 = count1 - (1/death_nbr)
			if(foobar[j] == 1):
				survival_2[j] = count2 - (1/survivor_nbr)
				count2 = count2 - (1/survivor_nbr)
		trace1 = go.Scatter(
		x=list(survival_1.keys()),
		y=list(survival_1.values()),
		mode='lines+markers',
		name="'Group 1'",
		hoverinfo='name',
		line=dict(
		shape='hv'))
		trace2 = go.Scatter(
		x=list(survival_2.keys()),
		y=list(survival_2.values()),
		mode='lines+markers',
		name="'Group 2'",
		hoverinfo='name',
		line=dict(
		shape='hv'))
		data99 = [trace1,trace2]
		layout = dict(
		legend=dict(
		y=0.5,
 		traceorder='reversed',
		font=dict(size=16)))
		fig = dict(data=data99, layout=layout)
		#plot_3 = go.Figure(data=data99,layout=layout)
		#plot_div = plot(plot_3, output_type='div', include_plotlyjs=False)
		#plotly.offline.plot(fig, filename='/home/quirin/testproject/polls/templates/polls/line-shapes.html')
		plot_div=plotly.offline.plot(fig, auto_open=False, output_type='div') 
		print(genelist_diff)
		plt.close()
		return(script,div,plot_1,plot_div,genelist_ret)
	def handle_upload_3_survival_4(fn,prot_fn):
		patients = []
		
		patients1 = []
		patients2 = []
		genes = []
		geneNames = []
		data = {}
		data1 = {}
		data2 = {}
		group1 = []
		group2 = []
		group_labels1 = []
		group_labels2 = []
		group1_data = []
		group2_data = []
		
		patients.append([1,2,3,4])
		patients1.append(['1','2'])
		#patients.append(['3','4'])
		patients2.append(['3','4'])
		group_labels1.append([1,1])
		group_labels2.append([2,2])
		logstring = "Creating Plots for given input files... \n\n"
		logstring = "Reading gene expression data... \n"
		ctr_bar = 0		
		df2 = pd.read_csv(fn,delim_whitespace=True,header=None,index_col=0)
		diff = {}
		ctr2 = 0
		#group1.append(group_labels1)
		#group2.append(group_labels2)
		group1.append([1,1])
		group2.append([2,2])		
		survival_info =  df2.loc['SURVIVAL',:].apply(pd.to_numeric).tolist()
		survival_yrs =  df2.loc['SURVIVAL_YRS',:].apply(pd.to_numeric).tolist()
		print(survival_info)
		for k in range(0,len(survival_info)):
			if(survival_info[k] == 1):
				group1.append(k)
			else:
				group2.append(k)
		#group1.append(df2.index[df2['SURVIVAL']==1].tolist())
		#group1.append(df2.index[df2['SURVIVAL']==0].tolist())
		genes = ()
		dataset = Dataset(name='hsapiens_gene_ensembl',
	                  host='http://www.ensembl.org')
		conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
		genes_3 = {}
		genes_corr = {}
		logstring = logstring + "\n\n Matching gene and protein IDs... \n"
		#foobar_2 = list(df2.iloc[4:,0])
		foobar_2 = list(df2.index)[4:]
		print(foobar_2)
		ensg_id_list = conv['Gene stable ID'].tolist()
		for i in foobar_2:
			#print(data1[1000])
			#here could be check whether level of expression is >100
			if(float(df2.loc[i].tolist()[2]) > 100.0):
				ensg_id = i.split('.')[0]
				if(ensg_id in ensg_id_list):
				#if(len(data1[i]) == len(patients1) and len(data2[i]) == len(patients2)):
					#print(data1[1000])
					ctr2 = ctr2 + 1
					#print(i)
					#if(data1[i][0] > 10000.0):
					#	print(data1[i])
					geneNames.append(ctr2)
					#group1.append(data1[i])
					#group2.append(data2[i])
					#print(conv.iloc[conv['Gene stable ID'] == 'ENSG00000223972']['Gene name'])
					
					print(ensg_id)
					prot_nbr = conv.index[conv['Gene stable ID'] == ensg_id].tolist()
					print(prot_nbr)
					prot_id = conv.loc[prot_nbr,'Gene name'].values[0]
					#print(prot_id)
					diff_curr = 0
					for j in range(0,(df2.shape[1]-1)):
						if(j in group1):
							print(df2.loc[i])
							diff_curr = diff_curr + float(df2.loc[i].tolist()[j])
						else:
							print(df2.loc[i])
							diff_curr = diff_curr - float(df2.loc[i].tolist()[j])
					if(1==1):
						#group1.append(data1[i])
						#group2.append(data2[i])
						logstring = logstring + str(prot_id)
						diff_curr_2 = diff_curr / 2000.0
						diff.update({i:diff_curr_2})
						genes_3.update({prot_id:diff_curr_2})
						tmp_3 = df2.loc[i].tolist()
						tmp_4 = [float(x) for x in tmp_3]
						tmp_corr = pearsonr(tmp_4,survival_yrs)
						print(tmp_corr)
						genes_corr.update({prot_id:tmp_corr})
						#print(prot_id)
		
		G = nx.Graph()
		nodes = []
		dataset = Dataset(name='hsapiens_gene_ensembl',
	                  host='http://www.ensembl.org')
		conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
		#print(conv.iloc[conv['Gene stable ID'] == 'ENSG00000223972']['Gene name'])
		#print(conv.index[conv['Gene stable ID'] == 'ENSG00000284332'].tolist())
		dataset.list_attributes()
		genes = {}
		logstring = logstring + "\n\nReading PPI file...\n"
		for line in prot_fn.readlines():
			lineSplit = line.decode('utf8').split()
			prot1 = lineSplit[0]
			prot2 = lineSplit[1]
			logstring = logstring + prot1 + prot2
			if(prot1 not in nodes and prot1 in genes_3):
				nodes.append(prot1)
				logstring = logstring + str(genes_3[prot1])
				G.add_node(prot1, Name=prot1, d=float(genes_3[prot1]))
				genes.update({prot1:genes_3[prot1]})
			if(prot2 not in nodes and prot2 in genes_3):
				nodes.append(prot2)
				G.add_node(prot2, Name=prot2, d=float(genes_3[prot2]))
				genes.update({prot2:genes_3[prot2]})
			logstring = logstring + str(genes_3[prot1]) + str(genes_3[prot2])
			if(prot1 in genes and prot2 in genes):		
				G.add_edge(prot1,prot2)
		#diff_frame = pd.DataFrame(list(genes.items()), columns=['Gene', 'Difference'])
		genelist_diff = []
		#diff_frame_2 = diff_frame.sort_values('Difference')
		#for i in range(0,len(diff_frame_2)):
		#	curr_tmp = diff_frame_2.iloc[i,0]
		#	print(curr_tmp)
		#	genelist_diff.append({'gene_name':curr_tmp,'difference':diff_frame_2.loc[i,'Difference']})
		genelist_ret =[]
		for k,v in list(genes.items()):
			#print(k)
			#print(v)
			genelist_diff.append({'gene_name':k,'difference':v})
			genelist_ret.append({'gene_name':k,'difference':v,'correlation':genes_corr[k]})
		output_notebook()
		plot = figure(x_range=(-1.5, 1.5), y_range=(-1.5, 1.5))
		# add tools to the plot
		plot.add_tools(HoverTool(tooltips=[("Diff", "@d")]), 
	               TapTool(), 
	               BoxSelectTool())
		# create bokeh graph
		graph = from_networkx(G, nx.spring_layout, iterations=1000, scale=1, center=(0,0))
		
		# add name to node data
		#graph.node_renderer.data_source.data['d'] = list(G.nodes())
		
		
		graph.node_renderer.data_source.data['d'] = [genes[i] for i in G.nodes]
		
		# add name to node data
		graph.node_renderer.data_source.data['Name'] = list(G.nodes())
		
		graph.edge_renderer.selection_glyph = MultiLine(line_color="#000000", line_width=4)
		graph.selection_policy = NodesAndLinkedEdges()
		graph.inspection_policy = EdgesAndLinkedNodes()
		graph.node_renderer.data_source.data['Name'] = list(G.nodes())
		x,y = zip(*graph.layout_provider.graph_layout.values())
		node_labels = nx.get_node_attributes(G, 'Name')
		z = tuple([(bar-0.1) for bar in x])
		#make bokeh graph with transparent labels
		source = ColumnDataSource(data=dict(x=x, y=y, 
	                           Name=[node_labels[i] for i in genes.keys()]))
		#labels2 = LabelSet(x='x', y='y', text='Name',text_font_size="8pt",x_offset=-20, source=source,
                #  background_fill_color='white',background_fill_alpha=0.0,level='glyph',render_mode='canvas')
		#graph.node_renderer.glyph = Circle(size=60, fill_color=linear_cmap('d', 'Spectral8', -2.5, 2.5))
		labels2 = LabelSet(x='x', y='y', text='Name',text_font_size="8pt",x_offset=-20, source=source,
                  background_fill_color='white',background_fill_alpha=0.0,level='glyph',render_mode='canvas')
		graph.node_renderer.glyph = Circle(size=60, fill_color=linear_cmap('d', 'Spectral8', -2.5, 2.5))
		plot.renderers.append(graph)
		plot.renderers.append(labels2)
		#graph.inspection_policy = EdgesAndLinkedNodes()
		plot.add_layout(labels2)
		logstring = logstring + "\n\nPPI Graph created..."
		output_file("polls/interactive_graphs.html")
		save(plot)
		#dataArray1 = np.array(group1)
		#dataArray2 = np.array(group2)
		#endData = np.concatenate((dataArray1, dataArray2), axis=1)
		#logstring = logstring + str(endData)
		red_patch = mpatches.Patch(color='red', label='Group1')
		blue_patch = mpatches.Patch(color='blue', label='Group2')
		#lut = dict(zip(set(endData[0]), sns.hls_palette(len(set(endData[0])), l=0.5, s=0.8)))
		#col_colors = pd.DataFrame(endData[0])[0].map(lut)
		#print(col_colors)
		colordict={0:'#BB0000',1:'#0000BB'}
		#colors = np.array(['#BB0000','#BB0000','#0000BB','#0000BB'])
		#df = pd.DataFrame(data=endData[1:,0:],index=geneNames,columns=patients)
		#df2 = pd.DataFrame(data=endData[0,0:], index='',columns=patients)
		#my_palette = dict(zip(df[.unique(), ["orange","yellow","brown"]))
		#row_colors = df2.cyl.map(my_palette)
		#fig, (ax1, ax2) = plt.subplots(1,2,sharex=True,sharey=True)
		#logstring = logstring + str(df)
		#df2 = pd.read_csv(fn,delim_whitespace=True,header=None,index_col=0)
		#foobar = df2.loc['SURVIVAL',:].apply(pd.to_numeric)
		#print(foobar)
		#print(colors)			
		#df22 = pd.read_csv(fn,delim_whitespace=True,header=None,index_col=0)
		df6 = df2.iloc[5:50000,1:].apply(pd.to_numeric)
		df5 = df2.sort_values(by='SURVIVAL',axis=1)
		foobar = df5.loc['SURVIVAL',:].apply(pd.to_numeric)
		#colors = [colordict[i] for i in foobar]
		print(df5)
		colors = [colordict[i] for i in foobar]
		#print(colors)
		df6 = df5.iloc[5:50000,:].apply(pd.to_numeric)
		print(df6)
		df7 = df6[df6.iloc[:,2]>100.0]
		print(df7)
		sns.clustermap(df7,method="average", figsize=(13, 13),col_colors=colors,col_cluster=False)
		#sns.clustermap(dataArray,method="average", figsize=(13, 13),ax=ax2)	
		#plt.imshow(dataArray, cmap='hot', interpolation='nearest')
		#plt.yticks(range(data_present_length), geneNames, size='small')
		#plt.xticks(range(len(patients)), patients, size='small')
		plt.legend(handles=[red_patch,blue_patch],loc=2,mode="expand",bbox_to_anchor=(3, 1))
		#plt.switch_backend('SVG') 
		plt.savefig("/home/quirin/testproject/polls/static/test3.png")
		print(logstring)
		script, div = components(plot)
		plot_1=plt.gcf()
		survival_yrs = df2.loc['SURVIVAL_YRS',:].apply(pd.to_numeric)
		survival_1 = {0:1.0}
		survival_2 = {0:1.0}
		survivor_nbr = sum(foobar)
		death_nbr = len(foobar) - sum(foobar)
		count1 = 1.0
		count2 = 1.0
		for j in range(1,len(survival_yrs)):
			if(foobar[j] == 0):
				survival_1[j] = count1 - (1/death_nbr)
				count1 = count1 - (1/death_nbr)
			if(foobar[j] == 1):
				survival_2[j] = count2 - (1/survivor_nbr)
				count2 = count2 - (1/survivor_nbr)
		trace1 = go.Scatter(
		x=list(survival_1.keys()),
		y=list(survival_1.values()),
		mode='lines+markers',
		name="'Group 1'",
		hoverinfo='name',
		line=dict(
		shape='hv'))
		trace2 = go.Scatter(
		x=list(survival_2.keys()),
		y=list(survival_2.values()),
		mode='lines+markers',
		name="'Group 2'",
		hoverinfo='name',
		line=dict(
		shape='hv'))
		data99 = [trace1,trace2]
		layout = dict(
		legend=dict(
		y=0.5,
 		traceorder='reversed',
		font=dict(size=16)))
		fig = dict(data=data99, layout=layout)
		#plot_3 = go.Figure(data=data99,layout=layout)
		#plot_div = plot(plot_3, output_type='div', include_plotlyjs=False)
		#plotly.offline.plot(fig, filename='/home/quirin/testproject/polls/templates/polls/line-shapes.html')
		plot_div=plotly.offline.plot(fig, auto_open=False, output_type='div') 
		print(genelist_diff)
		return(script,div,plot_1,plot_div,genelist_ret)
	def handle_upload_8(fn,prot_fn):
		patients = []
		patients1 = []
		patients2 = []
		genes = []
		geneNames = []
		data = {}
		group1 = []
		group2 = []
		group1_data = []
		group2_data = []
		logstring = "Creating Plots for given input files... \n\n"
		logstring = "Reading gene expression data... \n"
		ctr_bar = 0
		df2 = pd.read_csv(fn,delim_whitespace=True,header=None,index_col=0)
		diff = {}
		ctr2 = 0
		survival_info =  df2.loc['SURVIVAL',:].apply(pd.to_numeric).tolist()
		#assign patients to groups based on 'survival' (can also use 'cluster') column
		for k in range(0,len(survival_info)):
			if(survival_info[k] == 1):
				group1.append(k)
			else:
				group2.append(k)
		genes = ()
		#conv is the conversion dataframe from ENSGuvw.xy.z to PROT-ID
		conv = pd.read_csv("Output4.txt",delim_whitespace=True,header=0,index_col=0)	
		genes_3 = {}
		genes_4 = {}
		logstring = logstring + "\n\n Matching gene and protein IDs... \n"
		#SOME OLD CODE		
		#foobar_2 = list(df2.iloc[4:,0])	
		#print(df2)
		#df16 = df2.iloc[5:,:]
		#df17 = df16.loc[df2[1]>100.0]
		#print(df17)
		foobar_2 = list(df2.index)[4:]
		ensg_id_list = conv.iloc[:,0].tolist()
		conv_dict = conv.to_dict('index')
		conv_dict_2 = {}
		conv_dict_3 = {}
		df2shape = df2.shape[1]-1
		#fill dictionaries (prot -> gene and gene -> prot) and confirm gene occurence in expression data (dict->faster than check whole DF)
		for i in foobar_2:
			#here could be check whether level of expression is >100
			#if(float(df2.loc[i].tolist()[2]) > 100.0):
			ensg_id = i.split('.')[0]
			if(ensg_id in conv_dict):
				prot_id = conv_dict[ensg_id]['Gene_name']
				conv_dict_2[prot_id]=i	
				genes_4.update({prot_id:0})
		G = nx.Graph()
		G_2 = nx.Graph()
		nodes = []
		genes = {}
		genes_5 = {}
		genes_3 = {}
		#read file with PPI, calculate expression difference of respective genes between groups
		logstring = logstring + "\n\nReading PPI file...\n"
		for line in prot_fn.readlines():
			lineSplit = line.split()
			prot1 = lineSplit[0]
			prot2 = lineSplit[1]
			#logstring = logstring + prot1 + prot2
			#insert nodes for new genes, insert edge for each PPI
			if(prot1 not in nodes and prot1 in genes_4):
				nodes.append(prot1)
				diff = 0
				diff_curr = 0
				gene_tmp = conv_dict_2[prot1]
				for j in range(0,df2shape):
						if(j in group1):
							diff_curr = diff_curr + float(df2.loc[gene_tmp].tolist()[j])
						else:
							diff_curr = diff_curr - float(df2.loc[gene_tmp].tolist()[j])
				diff_curr_2 = diff_curr / 2000.0
				G.add_node(prot1, Name=prot1, d=float(diff_curr_2))
				genes_3.update({prot_id:diff_curr_2})
				genes.update({prot1:diff_curr_2})
			if(prot2 not in nodes and prot2 in genes_4):
				nodes.append(prot2)
				gene_tmp = conv_dict_2[prot2]
				diff_curr = 0
				for j in range(0,df2shape):
						if(j in group1):
							diff_curr = diff_curr + float(df2.loc[gene_tmp].tolist()[j])
						else:
							diff_curr = diff_curr - float(df2.loc[gene_tmp].tolist()[j])
				diff_curr_2 = diff_curr / 2000.0
				G.add_node(prot2, Name=prot2, d=float(diff_curr_2))
				genes_3.update({prot_id:diff_curr_2})
				genes.update({prot2:diff_curr_2})
				#logstring = logstring + str(genes_3[prot1]) + str(genes_3[prot2])
			if(prot2 in genes_4 and prot1 in genes_4):
				G.add_edge(prot1,prot2)
				if(prot2=="MT-ND1" or prot2=="MT-ND5"):
					G_2.add_node(prot1, Name=prot1, d=float(diff_curr_2))
					G_2.add_node(prot2, Name=prot2, d=float(diff_curr_2))
					G_2.add_edge(prot1,prot2)
					genes_5.update({prot1:diff_curr_2})
					genes_5.update({prot2:diff_curr_2})
		output_notebook()
		plot = figure(x_range=(-1.5, 1.5), y_range=(-1.5, 1.5))
		# add tools to the plot
		plot.add_tools(HoverTool(tooltips=None),
	               TapTool(), 
	               BoxSelectTool())
		# create bokeh graph
		#graph = from_networkx(G, nx.spring_layout, iterations=1000, scale=1, center=(0,0))
		graph = from_networkx(G, nx.circular_layout, scale=1, center=(0,0))
		graph_2 = from_networkx(G_2, nx.circular_layout, scale=1, center=(0,0))
		# add name to node data
		
		
		graph.node_renderer.data_source.data['d'] = [genes[i] for i in G.nodes]
		
		# add name to node data
		graph.edge_renderer.selection_glyph = MultiLine(line_color="#000000", line_width=4)
		graph.selection_policy = NodesAndLinkedEdges()
		graph.inspection_policy = EdgesAndLinkedNodes()
		graph.node_renderer.data_source.data['Name'] = list(G.nodes())
		x,y = zip(*graph.layout_provider.graph_layout.values())
		node_labels = nx.get_node_attributes(G, 'Name')
		z = tuple([(bar-0.1) for bar in x])
		#make bokeh graph with transparent labels
		source = ColumnDataSource(data=dict(x=x, y=y, 
	                           Name=[node_labels[i] for i in genes.keys()]))
		labels2 = LabelSet(x='x', y='y', text='Name',text_font_size="8pt",x_offset=-20, source=source,
                  background_fill_color='white',background_fill_alpha=0.0,level='glyph',render_mode='canvas')
		graph.node_renderer.glyph = Circle(size=60, fill_color=linear_cmap('d', 'Spectral8', -2.5, 2.5))

		graph_2.node_renderer.data_source.data['d'] = [genes[i] for i in G_2.nodes]
		
		# add name to node data
		graph_2.edge_renderer.selection_glyph = MultiLine(line_color="#000000", line_width=4)
		graph_2.selection_policy = NodesAndLinkedEdges()
		graph_2.inspection_policy = EdgesAndLinkedNodes()
		print(list(G_2.nodes()))
		graph_2.node_renderer.data_source.data['Name'] = list(G_2.nodes())
		x_2,y_2 = zip(*graph_2.layout_provider.graph_layout.values())
		node_labels_2 = nx.get_node_attributes(G_2, 'Name')
		z_2 = tuple([(bar-0.2) for bar in x_2])
		#make bokeh graph with transparent labels
		source_2 = ColumnDataSource(data=dict(x=x_2, y=y_2, 
	                           Name=[node_labels_2[i] for i in genes_5.keys()]))
		labels2_2 = LabelSet(x='x', y='y', text='Name',text_font_size="8pt",x_offset=-20, source=source_2,
                  background_fill_color='white',background_fill_alpha=0.0,level='glyph',render_mode='canvas')
		graph_2.node_renderer.glyph = Circle(size=100, fill_color=linear_cmap('d', 'Spectral8', -2.5, 2.5))
		plot.renderers.append(graph)
		plot.renderers.append(graph_2)
		plot.renderers.append(labels2_2)
		plot.renderers.append(labels2)
		plot.add_layout(labels2_2)
		plot.add_layout(labels2)
		#begin making heatmap here
		red_patch = mpatches.Patch(color='red', label='Group1')
		blue_patch = mpatches.Patch(color='blue', label='Group2')
		colordict={0:'#BB0000',1:'#0000BB'}
		df5 = df2.sort_values(by='SURVIVAL',axis=1)
		foobar = df5.loc['SURVIVAL',:].apply(pd.to_numeric)
		colors = [colordict[i] for i in foobar]
		df6 = df5.iloc[5:,:].apply(pd.to_numeric)
		#arbitrary choice of some genes to make heatmap smaller
		df7 = df6[df6.iloc[:,3]>100.0]
		print(df7)
		sns.clustermap(df7,method="average", figsize=(13, 13),robust=True,col_colors=colors,col_cluster=False)
		plt.legend(handles=[red_patch,blue_patch],loc=2,mode="expand",bbox_to_anchor=(3, 1))
		plt.savefig("/home/quirin/testproject/polls/static/test.png")
		print(logstring)
		script, div = components(plot)	
		plot_1=plt.gcf()
		#make survival info + plot based on surivival years in data. data w/o survival years are dealt in seperate methods
		survival_yrs = df2.loc['SURVIVAL_YRS',:].apply(pd.to_numeric)
		survival_1 = {0:1.0}
		survival_2 = {0:1.0}
		survivor_nbr = sum(foobar)
		death_nbr = len(foobar) - sum(foobar)
		count1 = 1.0
		count2 = 1.0
		for j in range(1,len(survival_yrs)):
			if(foobar[j] == 0):
				survival_1[j] = count1 - (1/death_nbr)
				count1 = count1 - (1/death_nbr)
			if(foobar[j] == 1):
				survival_2[j] = count2 - (1/survivor_nbr)
				count2 = count2 - (1/survivor_nbr)
		trace1 = go.Scatter(
		x=list(survival_1.keys()),
		y=list(survival_1.values()),
		mode='lines+markers',
		name="'Group 1'",
		hoverinfo='name',
		line=dict(
		shape='hv'))
		trace2 = go.Scatter(
		x=list(survival_2.keys()),
		y=list(survival_2.values()),
		mode='lines+markers',
		name="'Group 2'",
		hoverinfo='name',
		line=dict(
		shape='hv'))
		data99 = [trace1,trace2]
		layout = dict(
		legend=dict(
		y=0.5,
 		traceorder='reversed',
		font=dict(size=16)))
		fig = dict(data=data99, layout=layout)
		plot_div=plotly.offline.plot(fig, auto_open=False, output_type='div') 
		plt.close()
		return(script,div,plot_1,plot_div)
	def handle_upload_7(fn,prot_fn):
		patients = []
		patients1 = []
		patients2 = []
		genes = []
		geneNames = []
		data = {}
		group1 = []
		group2 = []
		group1_data = []
		group2_data = []
		logstring = "Creating Plots for given input files... \n\n"
		logstring = "Reading gene expression data... \n"
		ctr_bar = 0
		df2 = pd.read_csv(fn,delim_whitespace=True,header=None,index_col=0)
		diff = {}
		ctr2 = 0
		survival_info =  df2.loc['SURVIVAL',:].apply(pd.to_numeric).tolist()
		#assign patients to groups based on 'survival' (can also use 'cluster') column
		for k in range(0,len(survival_info)):
			if(survival_info[k] == 1):
				group1.append(k)
			else:
				group2.append(k)
		genes = ()
		#conv is the conversion dataframe from ENSGuvw.xy.z to PROT-ID
		conv = pd.read_csv("Output4.txt",delim_whitespace=True,header=0,index_col=0)	
		genes_3 = {}
		genes_4 = {}
		logstring = logstring + "\n\n Matching gene and protein IDs... \n"
		#SOME OLD CODE		
		#foobar_2 = list(df2.iloc[4:,0])	
		#print(df2)
		#df16 = df2.iloc[5:,:]
		#df17 = df16.loc[df2[1]>100.0]
		#print(df17)
		foobar_2 = list(df2.index)[4:]
		ensg_id_list = conv.iloc[:,0].tolist()
		conv_dict = conv.to_dict('index')
		conv_dict_2 = {}
		conv_dict_3 = {}
		df2shape = df2.shape[1]-1
		#fill dictionaries (prot -> gene and gene -> prot) and confirm gene occurence in expression data (dict->faster than check whole DF)
		for i in foobar_2:
			#here could be check whether level of expression is >100
			#if(float(df2.loc[i].tolist()[2]) > 100.0):
			ensg_id = i.split('.')[0]
			if(ensg_id in conv_dict):
				prot_id = conv_dict[ensg_id]['Gene_name']
				conv_dict_2[prot_id]=i	
				genes_4.update({prot_id:0})
		G = nx.Graph()
		nodes = []
		genes = {}
		genes_3 = {}
		#read file with PPI, calculate expression difference of respective genes between groups
		logstring = logstring + "\n\nReading PPI file...\n"
		for line in prot_fn.readlines():
			lineSplit = line.split()
			prot1 = lineSplit[0]
			prot2 = lineSplit[1]
			#logstring = logstring + prot1 + prot2
			#insert nodes for new genes, insert edge for each PPI
			if(prot1 not in nodes and prot1 in genes_4):
				nodes.append(prot1)
				diff = 0
				diff_curr = 0
				gene_tmp = conv_dict_2[prot1]
				for j in range(0,df2shape):
						if(j in group1):
							diff_curr = diff_curr + float(df2.loc[gene_tmp].tolist()[j])
						else:
							diff_curr = diff_curr - float(df2.loc[gene_tmp].tolist()[j])
				diff_curr_2 = diff_curr / 2000.0
				G.add_node(prot1, Name=prot1, d=float(diff_curr_2))
				genes_3.update({prot_id:diff_curr_2})
				genes.update({prot1:diff_curr_2})
			if(prot2 not in nodes and prot2 in genes_4):
				nodes.append(prot2)
				gene_tmp = conv_dict_2[prot2]
				diff_curr = 0
				for j in range(0,df2shape):
						if(j in group1):
							diff_curr = diff_curr + float(df2.loc[gene_tmp].tolist()[j])
						else:
							diff_curr = diff_curr - float(df2.loc[gene_tmp].tolist()[j])
				diff_curr_2 = diff_curr / 2000.0
				G.add_node(prot2, Name=prot2, d=float(diff_curr_2))
				genes_3.update({prot_id:diff_curr_2})
				genes.update({prot2:diff_curr_2})
				#logstring = logstring + str(genes_3[prot1]) + str(genes_3[prot2])
			if(prot2 in genes_4 and prot1 in genes_4):
				G.add_edge(prot1,prot2)
		output_notebook()
		plot = figure(x_range=(-1.5, 1.5), y_range=(-1.5, 1.5))
		# add tools to the plot
		plot.add_tools(HoverTool(tooltips=None),
	               TapTool(), 
	               BoxSelectTool())
		# create bokeh graph
		#graph = from_networkx(G, nx.spring_layout, iterations=1000, scale=1, center=(0,0))
		graph = from_networkx(G, nx.circular_layout, scale=1, center=(0,0))
		
		# add name to node data
		
		
		graph.node_renderer.data_source.data['d'] = [genes[i] for i in G.nodes]
		
		# add name to node data
		graph.edge_renderer.selection_glyph = MultiLine(line_color="#000000", line_width=4)
		graph.selection_policy = NodesAndLinkedEdges()
		graph.inspection_policy = EdgesAndLinkedNodes()
		graph.node_renderer.data_source.data['Name'] = list(G.nodes())
		x,y = zip(*graph.layout_provider.graph_layout.values())
		node_labels = nx.get_node_attributes(G, 'Name')
		z = tuple([(bar-0.1) for bar in x])
		#make bokeh graph with transparent labels
		source = ColumnDataSource(data=dict(x=x, y=y, 
	                           Name=[node_labels[i] for i in genes.keys()]))
		labels2 = LabelSet(x='x', y='y', text='Name',text_font_size="8pt",x_offset=-20, source=source,
                  background_fill_color='white',background_fill_alpha=0.0,level='glyph',render_mode='canvas')
		graph.node_renderer.glyph = Circle(size=60, fill_color=linear_cmap('d', 'Spectral8', -2.5, 2.5))
		plot.renderers.append(graph)
		plot.renderers.append(labels2)
		plot.add_layout(labels2)
		#begin making heatmap here
		red_patch = mpatches.Patch(color='red', label='Group1')
		blue_patch = mpatches.Patch(color='blue', label='Group2')
		colordict={0:'#BB0000',1:'#0000BB'}
		df5 = df2.sort_values(by='SURVIVAL',axis=1)
		foobar = df5.loc['SURVIVAL',:].apply(pd.to_numeric)
		colors = [colordict[i] for i in foobar]
		df6 = df5.iloc[5:,:].apply(pd.to_numeric)
		#arbitrary choice of some genes to make heatmap smaller
		df7 = df6[df6.iloc[:,3]>100.0]
		print(df7)
		sns.clustermap(df7,method="average", figsize=(13, 13),robust=True,col_colors=colors,col_cluster=False)
		plt.legend(handles=[red_patch,blue_patch],loc=2,mode="expand",bbox_to_anchor=(3, 1))
		plt.savefig("/home/quirin/testproject/polls/static/test.png")
		print(logstring)
		script, div = components(plot)	
		plot_1=plt.gcf()
		#make survival info + plot based on surivival years in data. data w/o survival years are dealt in seperate methods
		survival_yrs = df2.loc['SURVIVAL_YRS',:].apply(pd.to_numeric)
		survival_1 = {0:1.0}
		survival_2 = {0:1.0}
		survivor_nbr = sum(foobar)
		death_nbr = len(foobar) - sum(foobar)
		count1 = 1.0
		count2 = 1.0
		for j in range(1,len(survival_yrs)):
			if(foobar[j] == 0):
				survival_1[j] = count1 - (1/death_nbr)
				count1 = count1 - (1/death_nbr)
			if(foobar[j] == 1):
				survival_2[j] = count2 - (1/survivor_nbr)
				count2 = count2 - (1/survivor_nbr)
		trace1 = go.Scatter(
		x=list(survival_1.keys()),
		y=list(survival_1.values()),
		mode='lines+markers',
		name="'Group 1'",
		hoverinfo='name',
		line=dict(
		shape='hv'))
		trace2 = go.Scatter(
		x=list(survival_2.keys()),
		y=list(survival_2.values()),
		mode='lines+markers',
		name="'Group 2'",
		hoverinfo='name',
		line=dict(
		shape='hv'))
		data99 = [trace1,trace2]
		layout = dict(
		legend=dict(
		y=0.5,
 		traceorder='reversed',
		font=dict(size=16)))
		fig = dict(data=data99, layout=layout)
		plot_div=plotly.offline.plot(fig, auto_open=False, output_type='div') 
		plt.close()
		return(script,div,plot_1,plot_div)
	def handle_script_output(T,row_colors1,col_colors1,G2,means,genes_all,adjlist,genes1):
		#G = nx.Graph()
		#G_2 = nx.Graph()
		def color_for_graph(v):
			cmap_custom = {-4:'rgb(255, 0, 0)',-3:'rgb(255, 153, 51)',-2:'rgb(255, 204, 0)',-1:'rgb(255, 255, 0)',0:'rgb(204, 255, 51)',1:'rgb(153, 255, 51)',2:'rgb(102, 255, 51)',3:'rgb(51, 204, 51)'}
			v = v*2
			#tmp98 = int(v + (0.5 if v > 0 else -0.5))
			tmp98 = int(v)
			if(v < -4):
				tmp98 = -4
			if(v > 3):
				tmp98 = 3
			return(cmap_custom[tmp98])
		
		nodes = []
		nodecolors = []
		names = []
		genes = {}
		genes_5 = {}
		genes_3 = {}
		#read file with PPI, calculate expression difference of respective genes between groups
		G_list = list(G2.nodes())
		ctr = 0
		G = nx.Graph()
		print(G_list)
		print(genes1)
		#for G_tmp in G_list:
		#	genes.update({G_tmp:0})	
		#	tp = "circle"
		#	if(G_tmp in genes1):
		#		tp = "square"
		#	G.add_node(G_tmp, Name=G_tmp, d=float(means[ctr]),color=color_for_graph(means[ctr]),type=tp,label=G_tmp)
		#	#G.add_node(G_tmp, Name=G_tmp, d=float(means[ctr]),color=color_for_graph(means[ctr]),type=tp,label=G_tmp,borderColor="#000")
		#	#G.add_node(G_tmp, Name=G_tmp, d=float(means[ctr]))
		#	nodecolors.append(color_for_graph(means[ctr]))
		#	ctr = ctr + 1
		#	names.append(G_tmp)
		genelist_ret = []
		for G_tmp in genes_all:
			genes.update({G_tmp:0})	
			tp = "circle"
			if(G_tmp in genes1):
				tp = "square"
			G.add_node(G_tmp, Name=G_tmp, d=float(means[ctr]),color=color_for_graph(means[ctr]),type=tp,label=G_tmp)
			#G.add_node(G_tmp, Name=G_tmp, d=float(means[ctr]),color=color_for_graph(means[ctr]),type=tp,label=G_tmp,borderColor="#000")
			#G.add_node(G_tmp, Name=G_tmp, d=float(means[ctr]))
			nodecolors.append(color_for_graph(means[ctr]))
			#genelist_ret.append({'gene_name':G_tmp,'diff':float(means[ctr])})
			ctr = ctr + 1
			names.append(G_tmp)
		ctr = 0
		for edg in adjlist:
			G.add_edge(edg[0],edg[1],id=ctr,color="rgb(0,0,0)")
			ctr = ctr + 1
		pos = nx.spring_layout(G)	
		#print(pos)
		x_pos = {}
		y_pos = {}
		for k in pos:	
			x_pos[k] = pos[k][0]
			y_pos[k] = pos[k][1]
		#print(x_pos)
		edgl = {}
		ctr = 0
		
		#nodecolors = [color_for_graph(i) for i in means]
		#for edg in G.edges():
		#	edgl[edg] = ctr
		#	ctr = ctr + 1
		nx.set_node_attributes(G,x_pos,'x')
		nx.set_node_attributes(G,x_pos,'x')
		nx.set_node_attributes(G,y_pos,'y')
		#nx.set_edge_attributes(G,edgl,'id')
		nx.set_node_attributes(G,10,'size')
		#nx.set_node_attributes(G,nodecolors,'color')
		#print(json_graph.node_link_data(G))
		jsn = json_graph.node_link_data(G)
		jsn2 = str(json.dumps(jsn))
		jsn33 = jsn2.replace('links','edges')
		jsn44 = jsn33.replace('Label','label')
		jsn55 = jsn44.replace('bels','bel')
		jsn3 = jsn55.replace('\"directed\": false, \"multigraph\": false, \"graph\": {},','') 
		#print(jsn3)
		with open("polls/static/test15.json", "w") as text_file:
   			text_file.write(jsn3)		
		#nx.write_gexf(G,"test.gexf")
		output_notebook()
		plot = figure(x_range=(-1.5, 1.5), y_range=(-1.5, 1.5))
		# add tools to the plot
		plot.add_tools(HoverTool(tooltips=None),
	               TapTool(), 
	               BoxSelectTool())
		# create bokeh graph
		#graph = from_networkx(G, nx.spring_layout, iterations=1000, scale=1, center=(0,0))
		graph = from_networkx(G, nx.spring_layout, scale=1, center=(0,0))		
		#graph = from_networkx(G, nx.circular_layout, scale=1, center=(0,0))
		# add name to node data
		
		
		graph.node_renderer.data_source.data['d'] = [genes[i] for i in G.nodes]
		
		# add name to node data
		graph.edge_renderer.selection_glyph = MultiLine(line_color="#000000", line_width=4)
		graph.selection_policy = NodesAndLinkedEdges()
		graph.inspection_policy = EdgesAndLinkedNodes()
		graph.node_renderer.data_source.data['Name'] = list(G.nodes())
		x,y = zip(*graph.layout_provider.graph_layout.values())
		node_labels = nx.get_node_attributes(G, 'Name')
		z = tuple([(bar-0.1) for bar in x])
		#make bokeh graph with transparent labels
		source = ColumnDataSource(data=dict(x=x, y=y, 
	                           Name=[node_labels[i] for i in genes.keys()]))
		labels2 = LabelSet(x='x', y='y', text='Name',text_font_size="8pt",x_offset=-20, source=source,
                  background_fill_color='white',background_fill_alpha=0.0,level='glyph',render_mode='canvas')
		graph.node_renderer.glyph = Circle(size=60, fill_color=linear_cmap('d', 'Spectral8', -2.5, 2.5))

		plot.renderers.append(graph)
		plot.renderers.append(labels2)
		plot.add_layout(labels2)
		#begin making heatmap here
		red_patch = mpatches.Patch(color='#4FB6D3', label='SCC')
		blue_patch = mpatches.Patch(color='#22863E', label='ADK')
		colordict={0:'#BB0000',1:'#0000BB'}
		#sns.clustermap(T,method="average", figsize=(13, 13),robust=True,col_colors=col_colors1,row_colors=row_colors1,col_cluster=False)
		#sns.clustermap(T, figsize=(13, 13),col_colors=col_colors1,row_colors=row_colors1)		
		#plt.xlabel('Genes')
		#plt.ylabel('Patients')
		g = sns.clustermap(T, figsize=(13, 13),col_colors=col_colors1,row_colors=row_colors1)			
		ax = g.ax_heatmap
		ax.set_xlabel("Genes")
		ax.set_ylabel("Patients")
		#plt.legend(handles=[red_patch,blue_patch],loc=2,mode="expand",bbox_to_anchor=(3, 1))
		plt.savefig("/home/quirin/testproject/polls/static/test.png")
		script, div = components(plot)	
		plot_1=plt.gcf()
		plt.clf()
		patch_1 = mpatches.Patch(color='#ff0000', label='-2')
		patch_2 = mpatches.Patch(color='#ff9933', label='-1.5')
		patch_3 = mpatches.Patch(color='#ffcc00', label='-1')
		patch_4 = mpatches.Patch(color='#ffff00', label='-0.5')
		patch_5 = mpatches.Patch(color='#ccff33', label='0.5')
		patch_6 = mpatches.Patch(color='#99ff33', label='1')
		patch_7 = mpatches.Patch(color='#66ff33', label='1.5')
		patch_8 = mpatches.Patch(color='#33cc33', label='2')
		ax = plt.gca()
		ax.set_facecolor('white')
		plt.legend(handles=[patch_1,patch_2,patch_3,patch_4,patch_5,patch_6,patch_7,patch_8],loc=2,mode="expand",bbox_to_anchor=(0.5, 0.5))
		plt.savefig("/home/quirin/testproject/polls/static/legend_test.png")
		
		with open("polls/static/output_console.txt", "w") as text_file:
        		text_file.write("Done!")
		return(script,div,plot_1)



	def handle_script_output_ownfile(T,row_colors1,col_colors1,G2,means,adjlist,genes1):
		#G = nx.Graph()
		#G_2 = nx.Graph()
		def color_for_graph(v):
			cmap_custom = {-4:'rgb(255, 0, 0)',-3:'rgb(255, 153, 51)',-2:'rgb(255, 204, 0)',-1:'rgb(255, 255, 0)',0:'rgb(204, 255, 51)',1:'rgb(153, 255, 51)',2:'rgb(102, 255, 51)',3:'rgb(51, 204, 51)'}
			v = v*2
			#tmp98 = int(v + (0.5 if v > 0 else -0.5))
			tmp98 = int(v)
			if(v < -4):
				tmp98 = -4
			if(v > 3):
				tmp98 = 3
			return(cmap_custom[tmp98])
		
		nodes = []
		nodecolors = []
		names = []
		genes = {}
		genes_5 = {}
		genes_3 = {}
		#read file with PPI, calculate expression difference of respective genes between groups
		G_list = list(G2.nodes())
		ctr = 0
		G = nx.Graph()
		print(G_list)
		print(genes1)
		for G_tmp in G_list:
			genes.update({G_tmp:0})	
			tp = "circle"
			if(G_tmp in genes1):
				tp = "square"
			G.add_node(G_tmp, Name=G_tmp, d=float(means[ctr]),color=color_for_graph(means[ctr]),type=tp,label=G_tmp)
			#G.add_node(G_tmp, Name=G_tmp, d=float(means[ctr]),color=color_for_graph(means[ctr]),type=tp,label=G_tmp,borderColor="#000")
			#G.add_node(G_tmp, Name=G_tmp, d=float(means[ctr]))
			nodecolors.append(color_for_graph(means[ctr]))
			ctr = ctr + 1
			names.append(G_tmp)
		ctr = 0
		for edg in adjlist:
			G.add_edge(edg[0],edg[1],id=ctr)
			ctr = ctr + 1
		pos = nx.spring_layout(G)	
		#print(pos)
		x_pos = {}
		y_pos = {}
		for k in pos:	
			x_pos[k] = pos[k][0]
			y_pos[k] = pos[k][1]
		#print(x_pos)
		edgl = {}
		ctr = 0
		
		#nodecolors = [color_for_graph(i) for i in means]
		#for edg in G.edges():
		#	edgl[edg] = ctr
		#	ctr = ctr + 1
		nx.set_node_attributes(G,x_pos,'x')
		nx.set_node_attributes(G,x_pos,'x')
		nx.set_node_attributes(G,y_pos,'y')
		#nx.set_edge_attributes(G,edgl,'id')
		nx.set_node_attributes(G,10,'size')
		#nx.set_node_attributes(G,nodecolors,'color')
		#print(json_graph.node_link_data(G))
		jsn = json_graph.node_link_data(G)
		jsn2 = str(json.dumps(jsn))
		jsn33 = jsn2.replace('links','edges')
		jsn44 = jsn33.replace('Label','label')
		jsn55 = jsn44.replace('bels','bel')
		jsn3 = jsn55.replace('\"directed\": false, \"multigraph\": false, \"graph\": {},','') 
		#print(jsn3)
		with open("polls/static/test15.json", "w") as text_file:
   			text_file.write(jsn3)		
		#nx.write_gexf(G,"test.gexf")
		output_notebook()
		plot = figure(x_range=(-1.5, 1.5), y_range=(-1.5, 1.5))
		# add tools to the plot
		plot.add_tools(HoverTool(tooltips=None),
	               TapTool(), 
	               BoxSelectTool())
		# create bokeh graph
		#graph = from_networkx(G, nx.spring_layout, iterations=1000, scale=1, center=(0,0))
		graph = from_networkx(G, nx.spring_layout, scale=1, center=(0,0))		
		#graph = from_networkx(G, nx.circular_layout, scale=1, center=(0,0))
		# add name to node data
		
		
		graph.node_renderer.data_source.data['d'] = [genes[i] for i in G.nodes]
		
		# add name to node data
		graph.edge_renderer.selection_glyph = MultiLine(line_color="#000000", line_width=4)
		graph.selection_policy = NodesAndLinkedEdges()
		graph.inspection_policy = EdgesAndLinkedNodes()
		graph.node_renderer.data_source.data['Name'] = list(G.nodes())
		x,y = zip(*graph.layout_provider.graph_layout.values())
		node_labels = nx.get_node_attributes(G, 'Name')
		z = tuple([(bar-0.1) for bar in x])
		#make bokeh graph with transparent labels
		source = ColumnDataSource(data=dict(x=x, y=y, 
	                           Name=[node_labels[i] for i in genes.keys()]))
		labels2 = LabelSet(x='x', y='y', text='Name',text_font_size="8pt",x_offset=-20, source=source,
                  background_fill_color='white',background_fill_alpha=0.0,level='glyph',render_mode='canvas')
		graph.node_renderer.glyph = Circle(size=60, fill_color=linear_cmap('d', 'Spectral8', -2.5, 2.5))

		plot.renderers.append(graph)
		plot.renderers.append(labels2)
		plot.add_layout(labels2)
		#begin making heatmap here
		red_patch = mpatches.Patch(color='red', label='Group1')
		blue_patch = mpatches.Patch(color='blue', label='Group2')
		colordict={0:'#BB0000',1:'#0000BB'}
		sns.clustermap(T,method="average", figsize=(13, 13),robust=True,col_colors=col_colors1,row_colors=row_colors1,col_cluster=False)
		plt.legend(handles=[red_patch,blue_patch],loc=2,mode="expand",bbox_to_anchor=(3, 1))
		plt.savefig("/home/quirin/testproject/polls/static/test.png")
		script, div = components(plot)	
		plot_1=plt.gcf()
		return(script,div,plot_1)

	def handle_script_output_2(T,row_colors1,col_colors1,G2,means,adjlist):
		#G = nx.Graph()
		#G_2 = nx.Graph()
		nodes = []
		genes = {}
		genes_5 = {}
		genes_3 = {}
		#read file with PPI, calculate expression difference of respective genes between groups
		G_list = list(G2.nodes())
		ctr = 0
		G = nx.Graph()
		for G_tmp in G_list:
			genes.update({G_tmp:means[ctr]})	
			G.add_node(G_tmp, Name=G_tmp, d=float(means[ctr]))
			G.add_node(G_tmp, Name=G_tmp, d=float(means[ctr]))
			ctr = ctr + 1
		ctr = 0
		for edg in adjlist:
			G.add_edge(edg[0],edg[1],id=ctr)
			ctr = ctr + 1	
		#print(pos)
		#print(x_pos)
		edgl = {}
		ctr = 0
		#for edg in G.edges():
		#	edgl[edg] = ctr
		#	ctr = ctr + 1	
		#nx.write_gexf(G,"test.gexf")
		output_notebook()
		plot = figure(x_range=(-2.5, 2.5), y_range=(-2.5, 2.5),plot_width=1000,plot_height=1000)
		# add tools to the plot
		plot.add_tools(HoverTool(tooltips=None),
	               TapTool(), 
	               BoxSelectTool())
		# create bokeh graph
		#graph = from_networkx(G, nx.spring_layout, iterations=1000, scale=1, center=(0,0))
		graph = from_networkx(G, nx.circular_layout, scale=2, center=(0,0))		
		#graph = from_networkx(G, nx.circular_layout, scale=1, center=(0,0))
		# add name to node data
		
		
		graph.node_renderer.data_source.data['d'] = [genes[i] for i in G.nodes]
		
		# add name to node data
		graph.edge_renderer.selection_glyph = MultiLine(line_color="#000000", line_width=4)
		graph.selection_policy = NodesAndLinkedEdges()
		graph.inspection_policy = EdgesAndLinkedNodes()
		graph.node_renderer.data_source.data['Name'] = list(G.nodes())
		x,y = zip(*graph.layout_provider.graph_layout.values())
		node_labels = nx.get_node_attributes(G, 'Name')
		z = tuple([(bar-0.1) for bar in x])
		#make bokeh graph with transparent labels
		source = ColumnDataSource(data=dict(x=x, y=y, 
	                           Name=[node_labels[i] for i in genes.keys()]))
		labels2 = LabelSet(x='x', y='y', text='Name',text_font_size="8pt",x_offset=-20, source=source,
                  background_fill_color='white',background_fill_alpha=0.0,level='glyph',render_mode='canvas')
		graph.node_renderer.glyph = Circle(size=40, fill_color=linear_cmap('d', 'Spectral8', -2.5, 2.5))

		plot.renderers.append(graph)
		plot.renderers.append(labels2)
		plot.add_layout(labels2)
		#begin making heatmap here
		red_patch = mpatches.Patch(color='green', label='ADK')
		blue_patch = mpatches.Patch(color='blue', label='SDC')
		colordict={0:'#BB0000',1:'#0000BB'}
		sns.clustermap(T,method="average", figsize=(13, 13),robust=True,col_colors=col_colors1,row_colors=row_colors1,col_cluster=False)
		plt.legend(handles=[red_patch,blue_patch],loc=2,mode="expand",bbox_to_anchor=(3, 1))
		plt.savefig("/home/quirin/testproject/polls/static/test.png")
		script, div = components(plot)	
		plot_1=plt.gcf()
		return(script,div,plot_1)
	def makehref(term):
		ret = term + ""
		ret = "<a href=\"https://www.ncbi.nlm.nih.gov/gene/?term="+ret+"\">"+ret+"</a>"
		return(ret)
	def is_logged_in(username,password):
		return(1)
	
	def save_user_data(fn,prot_fn,username):
		foobar = "user_uploaded_files/" + username
		if not(os.path.isdir(foobar)):
			os.mkdir(foobar)
		fn.seek(0)
		prot_fn.seek(0)
		str1 = fn.read().decode('utf-8')
		str2 = prot_fn.read().decode('utf-8')
		#print(str1)
		filename_1 = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S_%f')[:-3]
		filename_2 = foobar + "/" + filename_1 + "_expr.txt"
		filename_3 = foobar + "/" + filename_1 + "_prot.txt"
		outfile1 = open(filename_2, "w")
		outfile1.write(str1)
		outfile1.close()
		outfile2 = open(filename_3, "w")
		outfile2.write(str2)
		outfile2.close()
	def save_user_data_2(fn,prot_fn,clinical_fn,username):
		foobar = "user_uploaded_files/" + username
		if not(os.path.isdir(foobar)):
			os.mkdir(foobar)
		fn.seek(0)
		prot_fn.seek(0)
		clinical_fn.seek(0)
		str1 = fn.read().decode('utf-8')
		str2 = prot_fn.read().decode('utf-8')
		str3 = clinical_fn.read().decode('utf-8')
		#print(str1)
		filename_1 = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S_%f')[:-3]
		filename_2 = foobar + "/" + filename_1 + "_expr.txt"
		filename_3 = foobar + "/" + filename_1 + "_prot.txt"
		filename_4 = foobar + "/" + filename_1 + "_clin.txt"
		outfile1 = open(filename_2, "w")
		outfile1.write(str1)
		outfile1.close()
		outfile2 = open(filename_3, "w")
		outfile2.write(str2)
		outfile2.close()
		outfile3 = open(filename_4, "w")
		outfile3.write(str3)
		outfile3.close()


	def save_results(username):
		foobar = "user_uploaded_files/" + username
		if not(os.path.isdir(foobar)):
			os.mkdir(foobar)
		fn.seek(0)
		prot_fn.seek(0)
		str1 = fn.read().decode('utf-8')
		str2 = prot_fn.read().decode('utf-8')
		#print(str1)
		filename_1 = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S_%f')[:-3]
		filename_2 = foobar + "/" + filename_1 + "_json.json"
		filename_3 = foobar + "/" + filename_1 + "_heatmap.png"
		outfile1 = open(filename_2, "w")
		outfile1.write(str1)
		outfile1.close()
		outfile2 = open(filename_3, "w")
		outfile2.write(str2)
		outfile2.close()
	
	def list_user_data(username):
		foobar = "user_uploaded_files/" + username
		fileslist = os.listdir(foobar)
		print(fileslist)
		bar = []
		for f in fileslist:
			if "_expr.txt" in f:
				print(f)
				ret1 = "user_uploaded_files/" + username + "/" + f
				fn_temp = f.split("_expr.txt")[0] + "_prot.txt"
				ret2 = "user_uploaded_files/" + username + "/" + fn_temp
				bar.append({'f1':ret1,'f2':ret2})
				#time_temp = f.split("_")
				#time_actual = time_temp[0] + "-" + time_temp[1] + "-" + time_temp[2] + " at " + time_temp[3] + ":" + time_temp[4]
		return bar
	def list_user_data_2(username):
		foobar = "user_uploaded_files/" + username
		fileslist = os.listdir(foobar)
		bar = []
		for f in fileslist:
			if "_json.json" in f:
				ret1 = "user_uploaded_files/" + username + "/" + f
				fn_temp = f.split("_json.json")[0] + "_heatmap.png"
				ret2 = "user_uploaded_files/" + username + "/" + fn_temp
				bar.append({'f1':ret1,'f2':ret2})
				#time_temp = f.split("_")
				#time_actual = time_temp[0] + "-" + time_temp[1] + "-" + time_temp[2] + " at " + time_temp[3] + ":" + time_temp[4]
		return bar
	
