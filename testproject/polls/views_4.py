from django.shortcuts import render




from io import StringIO
from django.http import HttpResponse
from django.template import loader
from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from django.shortcuts import get_object_or_404, render, redirect
from django.urls import reverse
from django.views import generic
import plotly
import plotly.plotly as py
import plotly.graph_objs as go
import plotly.io as pio
import plotly.offline
from .models import Choice, Question
from plotly.offline import plot_mpl
from sklearn.ensemble import RandomForestClassifier
from django.contrib.auth.models import User
from django.shortcuts import render_to_response,render
from django.template import RequestContext
from django.http import HttpResponseRedirect
from django.urls import reverse
import polls
from polls.models import Document
from polls.forms import DocumentForm
from polls.models import Upload,UploadForm,GraphForm
from django.contrib.auth import authenticate, login, logout
import numpy as np
import matplotlib.pyplot as plt
import mpld3
from polls.script3 import run_algo,algo_output,algo_output_ownfile,algo_output_ownfile_2,algo_output_ownfile_3
import seaborn as sns
import pandas as pd
from numpy import array
#from polls.tasks import run_analysis
import matplotlib.patches as mpatches
from polls.tasks import make_empty_figure,algo_output_task,script_output_task,empty_log_file,write_to_file_1,add_loading_image,remove_loading_image,script_output_task_2,show_old_data,script_output_task_3,write_metadata_to_file,list_metadata,metadata_to_string,script_output_task_4,read_ndex_file,read_ndex_file_2,read_ndex_file_3,list_metadata_2,list_metadata,script_output_task_7,script_output_task_8,script_output_task_9,list_metadata_3,run_enrichment,read_kegg_enrichment,run_go_enrichment,read_ndex_file_4,run_enrichment_2,run_go_enrichment_2,run_reac_enrichment,import_ndex,read_kegg_enrichment_2,convert_gene_list,check_input_files,script_output_task_10,list_metadata_4
import os.path
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
from django.contrib.auth import authenticate, login
from pybiomart import Dataset
from django.core.cache import cache
from datetime import datetime
def logout_2(request):
	if(request.method == "POST"):
		logout(request)
		return redirect('polls/logout.html')
	else:
		return render(request,'polls/logout.html')

def login_2(request):
	if('username' in request.POST and 'password' in request.POST):
		username = request.POST['username']
		password = request.POST['password']
		user = authenticate(request, username=username, password=password)
		if user is not None:
			login(request, user)
			#return render(request,'polls/login.html')
			return redirect('polls/clustering_6.html')
		else:
			return render(request,'polls/login.html')
			#return redirect('polls/clustering.html')
	else:
		return render(request,'polls/login.html')
		#return redirect('polls/clustering.html')
def clustering(request):
	ctr = 0
	div = ""
	script = ""
	plot2 = ""
	plot3 = ""
	genelist_diff = ""
	list_of_files = ""
	if(ctr == 0):
	        if request.method=="POST":
	                img = UploadForm(request.POST, request.FILES)
	                ctr = 1
	                filenames = request.FILES.items()
			#check if survival information is used
	                save_data = request.POST.get("save_data", None)
	                display_type = request.POST.get("include_survival", None)
	                if('pre_def_file' in request.POST):
	                	if(request.POST['pre_def_file']):
	                		fh1 = open("testgenes13.txt")
	                		fh2 = open("ppi_test_2.txt")
	                		#(script,div,plot2) = handle_upload_5(fh1,fh2)
	                		(script,div,plot2,plot3) = GraphForm.handle_upload_8(fh1,fh2)
	                		#(script,div,plot2,plot3) = handle_upload_5(fh1,fh2)
	                if('own_file' in request.POST):
	                	if(request.POST['own_file']):
	                		filename1 = request.POST.get("own_file")
	                		fh1 = open(filename1)
	                		filename2 = filename1.split("_expr.txt")[0] + "_prot.txt"
	                		fh2 = open(filename2)	                		
	                		#(script,div,plot2) = handle_upload_5(fh1,fh2)
	                		(script,div,plot2,plot3) = GraphForm.handle_upload_8(fh1,fh2)
	                		#(script,div,plot2,plot3) = handle_upload_5(fh1,fh2)
	                if('myfile' in request.FILES and 'prot_file' in request.FILES):
	                	if(request.FILES['myfile'] and request.FILES['prot_file']):
	                		if display_type in ["include_survival"]:
	                			#(script,div,plot2,plot3) = GraphForm.handle_upload_3_survival(request.FILES['myfile'],request.FILES['prot_file'])
	                			(script,div,plot2,plot3,genelist_diff) = GraphForm.handle_upload_3_survival_3(request.FILES['myfile'],request.FILES['prot_file'])
	                		else:
	                			(script,div,plot2) = GraphForm.handle_upload_3_4(request.FILES['myfile'],request.FILES['prot_file'])
	                	if save_data in ["save_data"]:
	                			if request.user.is_authenticated:
	                				username = str(request.user)
	                				GraphForm.save_user_data(request.FILES['myfile'],request.FILES['prot_file'],username)
	                if img.is_valid():
	                        img.save()
	        else:
	                img=UploadForm()
	                ctr = 1
	        if request.user.is_authenticated:
	        	username = str(request.user)
	        	list_of_files = GraphForm.list_user_data(username)
	        images=Upload.objects.all()
		#as of now, heatmap is stored statically
	        plot2 = "test.png"
		# Create user and save to the database
	       # user = User.objects.create_user('myusername', 'myemail@crazymail.com', 'mypassword')
		# Update fields and then save again
	       # user.first_name = 'John'
	       # user.last_name = 'Citizen'
	       # user.save()
	        #return render(request,'polls/clustering.html',{'form':img,'images':images, 'div':div, 'script':script,'plot2':plot2, 'plot_2':plot3, 'genelist_diff':genelist_diff})
	        return render(request,'polls/clustering.html',{'form':img,'images':images, 'div':div, 'script':script,'plot2':plot2, 'plot_2':plot3, 'genelist_diff':genelist_diff, 'list_of_files':list_of_files})


def clustering_2(request):
	if('L_g_min' in request.POST and 'L_g_max' in request.POST and 'cancer_type' in request.POST):
		lgmin = int(request.POST['L_g_min'])
		lgmax = int(request.POST['L_g_max'])
		cancer_type = request.POST.get("cancer_type")
		if(cancer_type == 1):
			pic = run_algo(1,lgmin,lgmax)
			#return render(request,'polls/login.html')
			return render(request, 'polls/clustering_2.html', {'pic':pic})
		else:
			pic = run_algo(2,lgmin,lgmax)
			return render(request,'polls/clustering_2.html', {'pic':pic})
			#return redirect('polls/clustering.html')
	else:
		return render(request,'polls/clustering_2.html')

def clustering_3(request):
	if('L_g_min' in request.POST and 'L_g_max' in request.POST and 'cancer_type' in request.POST):
		lgmin = int(request.POST['L_g_min'])
		lgmax = int(request.POST['L_g_max'])
		cancer_type = request.POST.get("cancer_type")
		if(cancer_type == 1):
			pic = run_algo(1,lgmin,lgmax)
			(T,row_colors,col_colors,G2,means,adjlist) = algo_output(1,lgmin,lgmax)
			(div,script,plot1) = GraphForm.handle_script_output_2(T,row_colors,col_colors,G2,means,adjlist)
			#return render(request,'polls/login.html')				
			plot2 = "test.png"
			return render(request, 'polls/clustering_3.html', {'form':"",'images':"",'div':div,'script':script,'plot2':plot2})
		else:
			(T,row_colors,col_colors,G2,means,adjlist) = algo_output(2,lgmin,lgmax)
			(div,script,plot1) = GraphForm.handle_script_output_2(T,row_colors,col_colors,G2,means,adjlist)
			#return render(request,'polls/login.html')				
			plot2 = "test.png"
			return render(request, 'polls/clustering_3.html', {'form':"",'images':"",'div':div,'script':script,'plot2':plot2})
		#return redirect('polls/clustering.html')
	else:
		return render(request,'polls/clustering_3.html')

def clustering_4(request):
	if('L_g_min' in request.POST and 'L_g_max' in request.POST and 'cancer_type' in request.POST):
		lgmin = int(request.POST['L_g_min'])
		lgmax = int(request.POST['L_g_max'])
		cancer_type = request.POST.get("cancer_type")
		if(cancer_type == 1):
			pic = run_algo(1,lgmin,lgmax)
			(T,row_colors,col_colors,G2,means,adjlist,genes1) = algo_output(1,lgmin,lgmax)
			(div,script,plot1) = GraphForm.handle_script_output(T,row_colors,col_colors,G2,means,adjlist,genes1)
			#return render(request,'polls/login.html')				
			plot2 = "test.png"
			return render(request, 'polls/clustering_4.html', {'form':"",'images':"",'div':div,'script':script,'plot2':plot2})
		else:
			(T,row_colors,col_colors,G2,means,adjlist,genes1) = algo_output(2,lgmin,lgmax)
			(div,script,plot1) = GraphForm.handle_script_output(T,row_colors,col_colors,G2,means,adjlist,genes1)
			#return render(request,'polls/login.html')				
			plot2 = "test.png"
			return render(request, 'polls/clustering_4.html', {'form':"",'images':"",'div':div,'script':script,'plot2':plot2})
		#return redirect('polls/clustering.html')
	else:
		return render(request,'polls/clustering_4.html')


def clustering_5(request):
	list_of_files = ""
	save_data = request.POST.get("save_data", None)
	nbr_iter = request.POST.get("nbr_iter",45)
	if('myfile' in request.FILES and 'protfile' in request.FILES):
		if(request.FILES['myfile'] and request.FILES['protfile']):
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1) = algo_output_ownfile(1,lgmin,lgmax,request.FILES['myfile'],request.FILES['protfile'],nbr_iter)
				#(T,row_colors,col_colors,G2,means,adjlist,genes1) = algo_output_ownfile(1,lgmin,lgmax,request.FILES['myfile'],request.FILES['protfile'])
				(div,script,plot1) = GraphForm.handle_script_output(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1)
				plot2 = "test.png"
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			username = str(request.user)
	                			GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
				#else:
					#cache.clear()
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data(username)	        
				#else:
				cache.clear()
				return render(request, 'polls/clustering_5.html', {'form':"",'images':"",'div':div,'script':script,'plot2':plot2, 'list_of_files':list_of_files})
	if('predef_file' in request.POST):
		if('L_g_min' in request.POST and 'L_g_max' in request.POST):
			if(request.POST['predef_file']):
				cancer_type = request.POST.get("cancer_type")
				fh2 = open("polls/data/biogrid.human.entrez.tsv")
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])				
				if(cancer_type == 1):
					fh1 = open("polls/data/lung_cancer_expr.csv")
				else:
					fh1 = open("polls/data/breast_cancer_expr.csv")
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1) = algo_output_ownfile(1,lgmin,lgmax,fh1,fh2,nbr_iter)
				#(T,row_colors,col_colors,G2,means,adjlist,genes1) = algo_output_ownfile(1,lgmin,lgmax,fh1,fh2)
				(div,script,plot1) = GraphForm.handle_script_output(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1)
				plot2 = "test.png"
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data(username)	        				
				#else:			
				cache.clear()
				return render(request, 'polls/clustering_5.html', {'form':"",'images':"",'div':div,'script':script,'plot2':plot2, 'list_of_files':list_of_files})
	if('own_file' in request.POST and request.user.is_authenticated):
		if(request.POST['own_file']):
			filename1 = request.POST.get("input_own_file")
			fh1 = open(filename1)
			filename2 = filename1.split("_expr.txt")[0] + "_prot.txt"
			fh2 = open(filename2)
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])				
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1) = algo_output_ownfile(1,lgmin,lgmax,fh1,fh2,nbr_iter)
				#(T,row_colors,col_colors,G2,means,adjlist,genes1) = algo_output_ownfile(1,lgmin,lgmax,fh1,fh2)
				(div,script,plot1) = GraphForm.handle_script_output(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1)
				plot2 = "test.png"
				cache.clear()
				return render(request, 'polls/clustering_5.html', {'form':"",'images':"",'div':div,'script':script,'plot2':plot2, 'list_of_files':list_of_files})       		
	                
	#if not(('myfile' in request.FILES and 'protfile' in request.FILES) or ('predef_file' in request.POST)):
	else:	
		if not (request.user.is_authenticated):
			cache.clear()
		else:
			username = str(request.user)
			list_of_files = GraphForm.list_user_data(username)
			cache.clear()	        				
		return render(request,'polls/clustering_5.html',{'list_of_files':list_of_files})





def clustering_6(request):
	ret_metadata1 = {}
	ret_metadata2 = {}
	ret_metadata3 = {}
	metadata_dict = []
	enrichment_dict = []
	pval_enr = 0.5
	list_of_files = ""
	list_of_files_2 = ""
	save_data = request.POST.get("save_data", None)
	nbr_iter = request.POST.get("nbr_iter",45)
	nbr_ants = request.POST.get("nbr_ants",30)
	evap = request.POST.get("evap",0.3)
	epsilon = request.POST.get("stopcr",0.02) 
	hi_sig = request.POST.get("hisig",1)
	pher_sig = request.POST.get("pher",1)
	savedata_param = "false"
	#fig = plt.figure(figsize=(10,8))
	#plt.savefig("polls/static/progress.png")
	#plt.close(fig)
	if request.user.is_authenticated:
		username = str(request.user)
		list_of_files = GraphForm.list_user_data_2(username)
	if('redo_analysis' in request.POST and request.user.is_authenticated):
		if(request.POST['redo_analysis']):
			with open("polls/static/output_console.txt", "w") as text_file:
   				text_file.write("Your request is being processed...")
			add_loading_image.delay()
			filename1 = request.POST.get("input_own_file_redo")
			fh1 = open(filename1)
			exprstr = fh1.read()
			filename2 = filename1.split("_expr.txt")[0] + "_prot.txt"
			fh2 = open(filename2)
			filename3 = filename1.split("_expr.txt")[0] + "_clin.txt"
			has_clin_data = "false"
			if(os.path.isfile(filename3)):
				fh3 = open(filename3)
				has_clin_data = "true"
				clinicalstr = fh3.read()
			#path_json = filename1
			#path_heatmap = filename1.split("_json.json")[0] + "_heatmap.png"
			ppistr = fh2.read()
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				make_empty_figure.delay()
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])	
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids) =result1.get()				
				if(has_clin_data == "true"):
					result2 = script_output_task_3.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr)
					(div,script,plot1,plot_div,ret_metadata) = result2.get()
				else:
					result2 = script_output_task.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1)
					(div,script,plot1) = result2.get()
					ret_metadata = ""
				plot2 = "test.png"
				cache.clear()			
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				list_of_files = GraphForm.list_user_data_2(username)	
				list_of_files_2 = GraphForm.list_user_data(username)      
				remove_loading_image.delay()
				return render(request, 'polls/clustering_6.html', {'form':"",'images':"",'div':"",'script':"",'plot2':plot2, 'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'ret_dat':ret_metadata})       		
	elif('myfile' in request.FILES and 'ndexfile' in request.FILES and 'protfile' in request.FILES):
		if(request.FILES['myfile'] and request.FILES['ndexfile'] and request.FILES['protfile']):
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				add_loading_image.delay()
				with open("polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				make_empty_figure.delay()
				#fh2 = open(request.FILES['protfile'])
				result3 = read_ndex_file_2.delay(request.FILES['ndexfile'].read().decode('utf-8'))
				ppistr_1 = result3.get()
				ppistr_2 = request.FILES['protfile'].read().decode('utf-8')
				ppistr = "\n".join([ppistr_2,ppistr_1])
				#fh1 = open(request.FILES['myfile'])
				exprstr = request.FILES['myfile'].read().decode('utf-8')
				savedata_param = "false"
				plot_div = ""
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			savedata_param = "true"
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				clinicalstr = "empty"
				if(request.user.is_authenticated and savedata_param == "true"):
					result2 = script_output_task_2.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,savedata_param,str(request.user))
					(div,script,plot1) = result2.get()			
				else:				
					result2 = script_output_task_4.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2)
					(div,script,plot1,plot_div,ret_metadata) = result2.get()
				plot2 = "test.png"
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			username = str(request.user)
	                			foobar = "user_uploaded_files/" + username
	                			if not(os.path.isdir(foobar)):
	                				os.mkdir(foobar)
	                			filename_1 = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S_%f')[:-3]
	                			filename_2 = foobar + "/" + username + filename_1 + "expr.txt"
	                			filename_3 = foobar + "/" + username + filename_1 + "prot.txt"
	                			with open(filename_2, "w") as text_file:
	                				text_file.write(exprstr)
	                				text_file.close()
	                			with open(filename_3, "w") as text_file_2:
	                				text_file_2.write(ppistr)
	                				text_file_2.close()
	                			#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
				#else:
					#cache.clear()
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)      	        
				#else:
				remove_loading_image.delay()
				cache.clear()				
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				return render(request, 'polls/clustering_6.html', {'list_of_files':list_of_files,'ret_metadata':ret_metadata,'result_99':ret_metadata,'list_of_files_2':list_of_files_2,'json_path':"test15.json",'heatmap_path':"test.png",'survival_graph':plot_div})
	elif('myfile' in request.FILES and 'ndex_file' in request.POST and 'ndex_name' in request.POST):
		print(request.POST.get("ndex_name"))
		if(request.FILES['myfile'] and request.POST['ndex_file'] and request.POST['ndex_name']):
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				add_loading_image.delay()
				with open("polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				make_empty_figure.delay()
				ndex_file_id = request.POST.get("ndex_name")
				if(ndex_file_id == "1"):
					fh2 = open("polls/data/APID.txt",encoding='utf-8')
				elif(ndex_file_id == "2"):
					fh2 = open("polls/data/STRING.txt",encoding='utf-8')
				elif(ndex_file_id == "3"):
					fh2 = open("polls/data/biogrid.human.entrez.tsv",encoding='utf-8')
				elif(ndex_file_id == "4"):
					fh2 = open("polls/data/HPRD.txt",encoding='utf-8')
				elif(ndex_file_id == "5"):
					fh2 = open("polls/data/ulitsky.txt",encoding='utf-8')
				elif(ndex_file_id == "6"):
					fh2 = open("polls/data/I2D.txt",encoding='utf-8')
				elif(ndex_file_id == "7"):
					fh2 = open("polls/data/biogrid_cel.txt",encoding='utf-8')
				elif(ndex_file_id == "8"):
					fh2 = open("polls/data/biogrid_ath.txt",encoding='utf-8')
				elif(ndex_file_id == "9"):
					fh2 = open("polls/data/biogrid_mus.txt",encoding='utf-8')
				elif(ndex_file_id == "10"):
					fh2 = open("polls/data/biogrid_dme.txt",encoding='utf-8')
				elif(ndex_file_id == "11"):
					fh2 = open("polls/data/biogrid_eco.txt",encoding='utf-8')
				elif(ndex_file_id == "12"):
					fh2 = open("polls/data/biogrid_sce.txt",encoding='utf-8')
				elif(ndex_file_id == "13"):
					fh2 = open("polls/data/biogrid_dan.txt",encoding='utf-8')
				ppistr = fh2.read()
				exprstr = request.FILES['myfile'].read().decode('utf-8')
				savedata_param = "false"
				plot_div = ""
				#if save_data in ["save_data"]:
	                	#	if request.user.is_authenticated:
	                	#		savedata_param = "true"
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				clinicalstr = "empty"
				if(request.user.is_authenticated and savedata_param == "true"):
					result2 = script_output_task_2.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,savedata_param,str(request.user))
					(div,script,plot1) = result2.get()			
				else:				
					result2 = script_output_task_4.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2)
					(div,script,plot1,plot_div,ret_metadata) = result2.get()
				plot2 = "test.png"
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			username = str(request.user)
	                			foobar = "user_uploaded_files/" + username
	                			if not(os.path.isdir(foobar)):
	                				os.mkdir(foobar)
	                			filename_1 = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S_%f')[:-3]
	                			filename_2 = foobar + "/" + username + filename_1 + "expr.txt"
	                			filename_3 = foobar + "/" + username + filename_1 + "prot.txt"
	                			with open(filename_2, "w") as text_file:
	                				text_file.write(exprstr)
	                				text_file.close()
	                			with open(filename_3, "w") as text_file_2:
	                				text_file_2.write(ppistr)
	                				text_file_2.close()
	                			#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
						#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
				#else:
					#cache.clear()
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)      	        
				#else:
				remove_loading_image.delay()
				cache.clear()				
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				return render(request, 'polls/clustering_6.html', {'list_of_files':list_of_files,'ret_metadata':ret_metadata,'result_99':ret_metadata,'list_of_files_2':list_of_files_2,'json_path':"test15.json",'heatmap_path':"test.png",'survival_graph':plot_div})
	elif('myfile' in request.FILES and 'parse_ndex_file' in request.POST and 'ndex_name_2' in request.POST):
		#print(request.POST.get("ndex_name"))
		if(request.FILES['myfile'] and request.POST['parse_ndex_file']):
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				add_loading_image.delay()
				with open("polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				make_empty_figure.delay()
				ndex_file_id = request.POST.get("ndex_name_2")
				if(ndex_file_id == "1"):
					result_ndex = import_ndex.delay("9c38ce6e-c564-11e8-aaa6-0ac135e8bacf")
					ppistr = result_ndex.get()
				elif(ndex_file_id == "2"):
					result_ndex = import_ndex.delay("275bd84e-3d18-11e8-a935-0ac135e8bacf")
					ppistr = result_ndex.get()
				elif(ndex_file_id == "3"):
					result_ndex = import_ndex.delay("becec556-86d4-11e7-a10d-0ac135e8bacf")
					ppistr = result_ndex.get()
				elif(ndex_file_id == "4"):
					result_ndex = import_ndex.delay("1093e665-86da-11e7-a10d-0ac135e8bacf")
					ppistr = result_ndex.get()
				#if(ndex_file_id == "1"):
				#	fh2 = open("test_scripts_3/APID.cx",encoding='utf-8')
				#elif(ndex_file_id == "2"):
				#	fh2 = open("test_scripts_3/STRING.cx",encoding='utf-8')
				#elif(ndex_file_id == "3"):
				#	fh2 = open("polls/data/biogrid.human.entrez.tsv",encoding='utf-8')
				#	ppistr = fh2.read()
				#elif(ndex_file_id == "4"):
				#	fh2 = open("test_scripts_3/HPRD.cx",encoding='utf-8')
				elif(ndex_file_id == "5"):
					fh2 = open("polls/data/ulitsky.txt",encoding='utf-8')
					ppistr = fh2.read()
				elif(ndex_file_id == "6"):
					fh2 = open("polls/data/I2D.txt",encoding='utf-8')
					ppistr = fh2.read()
				elif(ndex_file_id == "7"):
					fh2 = open("test_scripts_3/biogrid_cel.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "8"):
					fh2 = open("test_scripts_3/biogrid_ath.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "9"):
					fh2 = open("test_scripts_3/biogrid_mus.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "10"):
					fh2 = open("test_scripts_3/biogrid_dro.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "11"):
					fh2 = open("test_scripts_3/biogrid_eco.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "12"):
					fh2 = open("test_scripts_3/biogrid_sce.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "13"):
					fh2 = open("test_scripts_3/biogrid_dan.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				#if not(ndex_file_id == "1" or ndex_file_id == "2" or ndex_file_id == "3" or ndex_file_id=="6" or ndex_file_id == "5"):
				exprstr = request.FILES['myfile'].read().decode('utf-8')
				result4 = check_input_files.delay(ppistr,exprstr)
				errstr = result4.get()
				if(errstr != ""):
					request.session['errors'] = errstr
					return render(request,'polls/errorpage.html',{'errors':errstr})
				savedata_param = "false"
				plot_div = ""
				#if save_data in ["save_data"]:
	                	#	if request.user.is_authenticated:
	                	#		savedata_param = "true"
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				clinicalstr = "empty"
				if(request.user.is_authenticated and savedata_param == "true"):
					result2 = script_output_task_2.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,savedata_param,str(request.user))
					(div,script,plot1) = result2.get()			
				else:				
					clinicalstr = "empty"
					clinicaldf = ""
					survival_col_name = ""
					result2 = script_output_task_9.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf)
					(div,script,plot1,plot_div,ret_metadata) = result2.get()	
				
					#result2 = script_output_task_4.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2)
					#(div,script,plot1,plot_div,ret_metadata) = result2.get()
				plot2 = "test.png"
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			username = str(request.user)
	                			foobar = "user_uploaded_files/" + username
	                			if not(os.path.isdir(foobar)):
	                				os.mkdir(foobar)
	                			filename_1 = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S_%f')[:-3]
	                			filename_2 = foobar + "/" + username + filename_1 + "expr.txt"
	                			filename_3 = foobar + "/" + username + filename_1 + "prot.txt"
	                			with open(filename_2, "w") as text_file:
	                				text_file.write(exprstr)
	                				text_file.close()
	                			with open(filename_3, "w") as text_file_2:
	                				text_file_2.write(ppistr)
	                				text_file_2.close()
	                			#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
						#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
				#else:
					#cache.clear()
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)      	        
				#else:
				remove_loading_image.delay()
				cache.clear()				
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				return render(request, 'polls/clustering_6.html', {'list_of_files':list_of_files,'ret_metadata':ret_metadata,'result_99':ret_metadata,'list_of_files_2':list_of_files_2,'json_path':"test15.json",'heatmap_path':"test.png",'survival_graph':plot_div})
	elif('predef_file' in request.FILES and 'parse_ndex_file' in request.POST and 'ndex_name_2' in request.POST):
		#print(request.POST.get("ndex_name"))
		if(request.FILES['myfile'] and request.POST['parse_ndex_file'] and request.POST['predef_file']):
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				add_loading_image.delay()
				cancer_type = request.POST.get("cancer_type")
				if(cancer_type == "1"):
					print("babababababa")
					fh1 = open("polls/data/lung_cancer_expr.csv")
					exprstr = fh1.read()
					clinicaldf = pd.read_csv("polls/data/lung_cancer_clinical.csv")
					fh4 = open("polls/data/lung_cancer_clinical.csv")
					#fh4 = open("polls/data/breast_cancer_clinical.csv")
					is_lungc = "true"
					clinicalstr = fh4.read()
					fh4.flush()
					fh4.close()
					survival_col_name = "disease free survival in months:ch1"
				else:
					fh1 = open("polls/data/breast_cancer_expr.csv")
					exprstr = fh1.read()
					clinicaldf = pd.read_csv("polls/data/breast_cancer_clinical.csv")
					fh4 = open("polls/data/breast_cancer_clinical.csv")
					clinicalstr = fh4.read()
					fh4.flush()
					fh4.close()
					survival_col_name = "mfs (yr):ch1"
				
				with open("polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				make_empty_figure.delay()
				ndex_file_id = request.POST.get("ndex_name_2")
				if(ndex_file_id == "1"):
					result_ndex = import_ndex.delay("9c38ce6e-c564-11e8-aaa6-0ac135e8bacf")
					ppistr = result_ndex.get()
				elif(ndex_file_id == "2"):
					result_ndex = import_ndex.delay("275bd84e-3d18-11e8-a935-0ac135e8bacf")
					ppistr = result_ndex.get()
				elif(ndex_file_id == "3"):
					result_ndex = import_ndex.delay("becec556-86d4-11e7-a10d-0ac135e8bacf")
					ppistr = result_ndex.get()
				elif(ndex_file_id == "4"):
					result_ndex = import_ndex.delay("1093e665-86da-11e7-a10d-0ac135e8bacf")
					ppistr = result_ndex.get()
				#if(ndex_file_id == "1"):
				#	fh2 = open("test_scripts_3/APID.cx",encoding='utf-8')
				#elif(ndex_file_id == "2"):
				#	fh2 = open("test_scripts_3/STRING.cx",encoding='utf-8')
				#elif(ndex_file_id == "3"):
				#	fh2 = open("polls/data/biogrid.human.entrez.tsv",encoding='utf-8')
				#	ppistr = fh2.read()
				#elif(ndex_file_id == "4"):
				#	fh2 = open("test_scripts_3/HPRD.cx",encoding='utf-8')
				elif(ndex_file_id == "5"):
					fh2 = open("polls/data/ulitsky.txt",encoding='utf-8')
					ppistr = fh2.read()
				elif(ndex_file_id == "6"):
					fh2 = open("polls/data/I2D.txt",encoding='utf-8')
					ppistr = fh2.read()
				elif(ndex_file_id == "7"):
					fh2 = open("test_scripts_3/biogrid_cel.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "8"):
					fh2 = open("test_scripts_3/biogrid_ath.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "9"):
					fh2 = open("test_scripts_3/biogrid_mus.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "10"):
					fh2 = open("test_scripts_3/biogrid_dro.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "11"):
					fh2 = open("test_scripts_3/biogrid_eco.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "12"):
					fh2 = open("test_scripts_3/biogrid_sce.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "13"):
					fh2 = open("test_scripts_3/biogrid_dan.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				#if not(ndex_file_id == "1" or ndex_file_id == "2" or ndex_file_id == "3" or ndex_file_id=="6" or ndex_file_id == "5"):
				exprstr = request.FILES['myfile'].read().decode('utf-8')
				savedata_param = "false"
				plot_div = ""
				#if save_data in ["save_data"]:
	                	#	if request.user.is_authenticated:
	                	#		savedata_param = "true"
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				clinicalstr = "empty"
				if(request.user.is_authenticated and savedata_param == "true"):
					result2 = script_output_task_2.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,savedata_param,str(request.user))
					(div,script,plot1) = result2.get()			
				else:				
					clinicalstr = "empty"
					clinicaldf = ""
					survival_col_name = ""
					result2 = script_output_task_9.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf)
					(div,script,plot1,plot_div,ret_metadata) = result2.get()	
				
					#result2 = script_output_task_4.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2)
					#(div,script,plot1,plot_div,ret_metadata) = result2.get()
				plot2 = "test.png"
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			username = str(request.user)
	                			foobar = "user_uploaded_files/" + username
	                			if not(os.path.isdir(foobar)):
	                				os.mkdir(foobar)
	                			filename_1 = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S_%f')[:-3]
	                			filename_2 = foobar + "/" + username + filename_1 + "expr.txt"
	                			filename_3 = foobar + "/" + username + filename_1 + "prot.txt"
	                			with open(filename_2, "w") as text_file:
	                				text_file.write(exprstr)
	                				text_file.close()
	                			with open(filename_3, "w") as text_file_2:
	                				text_file_2.write(ppistr)
	                				text_file_2.close()
	                			#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
						#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
				#else:
					#cache.clear()
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)      	        
				#else:
				remove_loading_image.delay()
				cache.clear()				
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				return render(request, 'polls/clustering_6.html', {'list_of_files':list_of_files,'ret_metadata':ret_metadata,'result_99':ret_metadata,'list_of_files_2':list_of_files_2,'json_path':"test15.json",'heatmap_path':"test.png",'survival_graph':plot_div})

	elif('myfile' in request.FILES and 'parse_ndex_file' in request.POST and 'ndex_name_2' in request.POST and 'analyze_metadata' in request.POST):
		#print(request.POST.get("ndex_name"))
		if(request.FILES['myfile'] and request.POST['parse_ndex_file'] and request.FILES['patientdata']):
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				add_loading_image.delay()
				with open("polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				make_empty_figure.delay()
				ndex_file_id = request.POST.get("ndex_name_2")
				if(ndex_file_id == "1"):
					result_ndex = import_ndex.delay("9c38ce6e-c564-11e8-aaa6-0ac135e8bacf")
					ppistr = result_ndex.get()
				elif(ndex_file_id == "2"):
					result_ndex = import_ndex.delay("275bd84e-3d18-11e8-a935-0ac135e8bacf")
					ppistr = result_ndex.get()
				elif(ndex_file_id == "3"):
					result_ndex = import_ndex.delay("becec556-86d4-11e7-a10d-0ac135e8bacf")
					ppistr = result_ndex.get()
				elif(ndex_file_id == "4"):
					result_ndex = import_ndex.delay("1093e665-86da-11e7-a10d-0ac135e8bacf")
					ppistr = result_ndex.get()
				#if(ndex_file_id == "1"):
				#	fh2 = open("test_scripts_3/APID.cx",encoding='utf-8')
				#elif(ndex_file_id == "2"):
				#	fh2 = open("test_scripts_3/STRING.cx",encoding='utf-8')
				#elif(ndex_file_id == "3"):
				#	fh2 = open("polls/data/biogrid.human.entrez.tsv",encoding='utf-8')
				#	ppistr = fh2.read()
				#elif(ndex_file_id == "4"):
				#	fh2 = open("test_scripts_3/HPRD.cx",encoding='utf-8')
				elif(ndex_file_id == "5"):
					fh2 = open("polls/data/ulitsky.txt",encoding='utf-8')
					ppistr = fh2.read()
				elif(ndex_file_id == "6"):
					fh2 = open("polls/data/I2D.txt",encoding='utf-8')
					ppistr = fh2.read()
				elif(ndex_file_id == "7"):
					fh2 = open("test_scripts_3/biogrid_cel.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "8"):
					fh2 = open("test_scripts_3/biogrid_ath.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "9"):
					fh2 = open("test_scripts_3/biogrid_mus.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "10"):
					fh2 = open("test_scripts_3/biogrid_dro.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "11"):
					fh2 = open("test_scripts_3/biogrid_eco.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "12"):
					fh2 = open("test_scripts_3/biogrid_sce.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "13"):
					fh2 = open("test_scripts_3/biogrid_dan.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				#if not(ndex_file_id == "1" or ndex_file_id == "2" or ndex_file_id == "3" or ndex_file_id=="6" or ndex_file_id == "5"):
				exprstr = request.FILES['myfile'].read().decode('utf-8')
				result4 = check_input_files.delay(ppistr,exprstr)
				errstr = result4.get()
				if(errstr != ""):
					request.session['errors'] = errstr
					return render(request,'polls/errorpage.html',{'errors':errstr})
				savedata_param = "false"
				plot_div = ""
				clinicalstr = request.FILES['patientdata'].read().decode('utf-8')
				clinical_stringio = StringIO(clinicalstr)
				clinicaldf = pd.read_csv(clinical_stringio)
				survival_col_name = ""
				if('survival_col' in request.POST):
					print("barabsfrbasdb")
					if(request.POST['survival_col']):
						survival_col_name = request.POST['survival_col']
				
				#if save_data in ["save_data"]:
	                	#	if request.user.is_authenticated:
	                	#		savedata_param = "true"
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				clinicalstr = "empty"
				if(request.user.is_authenticated and savedata_param == "true"):
					result2 = script_output_task_2.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,savedata_param,str(request.user))
					(div,script,plot1) = result2.get()			
				else:				
					clinicalstr = "empty"
					clinicaldf = ""
					survival_col_name = ""
					result2 = script_output_task_9.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf)
					(div,script,plot1,plot_div,ret_metadata) = result2.get()	
				
					#result2 = script_output_task_4.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2)
					#(div,script,plot1,plot_div,ret_metadata) = result2.get()
				plot2 = "test.png"
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			username = str(request.user)
	                			foobar = "user_uploaded_files/" + username
	                			if not(os.path.isdir(foobar)):
	                				os.mkdir(foobar)
	                			filename_1 = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S_%f')[:-3]
	                			filename_2 = foobar + "/" + username + filename_1 + "expr.txt"
	                			filename_3 = foobar + "/" + username + filename_1 + "prot.txt"
	                			with open(filename_2, "w") as text_file:
	                				text_file.write(exprstr)
	                				text_file.close()
	                			with open(filename_3, "w") as text_file_2:
	                				text_file_2.write(ppistr)
	                				text_file_2.close()
	                			#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
						#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
				#else:
					#cache.clear()
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)      	        
				#else:
				remove_loading_image.delay()
				cache.clear()				
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				return render(request, 'polls/clustering_6.html', {'list_of_files':list_of_files,'ret_metadata':ret_metadata,'result_99':ret_metadata,'list_of_files_2':list_of_files_2,'json_path':"test15.json",'heatmap_path':"test.png",'survival_graph':plot_div})

	elif('myfile' in request.FILES and 'parse_ndex_file' in request.POST and 'ndex_file_input' in request.FILES):
		#print(request.POST.get("ndex_name"))
		if(request.FILES['myfile'] and request.POST['parse_ndex_file'] and request.FILES['ndex_file_input']):
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				add_loading_image.delay()
				with open("polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				make_empty_figure.delay()
				ndex_file_str = request.FILES['ndex_file_input'].read().decode('utf-8')				
				result_ndex = read_ndex_file_4.delay(ndex_file_str)
				ppistr = result_ndex.get()
				exprstr = request.FILES['myfile'].read().decode('utf-8')
				savedata_param = "false"
				plot_div = ""
				#if save_data in ["save_data"]:
	                	#	if request.user.is_authenticated:
	                	#		savedata_param = "true"
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				clinicalstr = "empty"
				if(request.user.is_authenticated and savedata_param == "true"):
					result2 = script_output_task_2.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,savedata_param,str(request.user))
					(div,script,plot1) = result2.get()			
				else:				
					result2 = script_output_task_4.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2)
					(div,script,plot1,plot_div,ret_metadata) = result2.get()
				plot2 = "test.png"
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			username = str(request.user)
	                			foobar = "user_uploaded_files/" + username
	                			if not(os.path.isdir(foobar)):
	                				os.mkdir(foobar)
	                			filename_1 = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S_%f')[:-3]
	                			filename_2 = foobar + "/" + username + filename_1 + "expr.txt"
	                			filename_3 = foobar + "/" + username + filename_1 + "prot.txt"
	                			with open(filename_2, "w") as text_file:
	                				text_file.write(exprstr)
	                				text_file.close()
	                			with open(filename_3, "w") as text_file_2:
	                				text_file_2.write(ppistr)
	                				text_file_2.close()
	                			#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
						#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
				#else:
					#cache.clear()
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)      	        
				#else:
				remove_loading_image.delay()
				cache.clear()				
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				return render(request, 'polls/clustering_6.html', {'list_of_files':list_of_files,'ret_metadata':ret_metadata,'result_99':ret_metadata,'list_of_files_2':list_of_files_2,'json_path':"test15.json",'heatmap_path':"test.png",'survival_graph':plot_div})

	elif('myfile' in request.FILES and 'ndex_file' in request.POST and 'ndex_name' in request.POST and 'analyze_metadata' in request.POST):
		print(request.POST.get("ndex_name"))
		if(request.FILES['myfile'] and request.POST['ndex_file'] and request.POST['ndex_name'] and request.FILES['patientdata']):
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				add_loading_image.delay()
				with open("polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				make_empty_figure.delay()
				ndex_file_id = request.POST.get("ndex_name")
				if(ndex_file_id == "1"):
					fh2 = open("polls/data/APID.txt",encoding='utf-8')
				elif(ndex_file_id == "2"):
					fh2 = open("polls/data/STRING.txt",encoding='utf-8')
				elif(ndex_file_id == "3"):
					fh2 = open("polls/data/biogrid.human.entrez.tsv",encoding='utf-8')
				elif(ndex_file_id == "4"):
					fh2 = open("polls/data/HPRD.txt",encoding='utf-8')
				elif(ndex_file_id == "5"):
					fh2 = open("polls/data/ulitsky.txt",encoding='utf-8')
				elif(ndex_file_id == "6"):
					fh2 = open("polls/data/I2D.txt",encoding='utf-8')
				elif(ndex_file_id == "7"):
					fh2 = open("polls/data/biogrid_cel.txt",encoding='utf-8')
				elif(ndex_file_id == "8"):
					fh2 = open("polls/data/biogrid_ath.txt",encoding='utf-8')
				elif(ndex_file_id == "9"):
					fh2 = open("polls/data/biogrid_mus.txt",encoding='utf-8')
				elif(ndex_file_id == "10"):
					fh2 = open("polls/data/biogrid_dme.txt",encoding='utf-8')
				elif(ndex_file_id == "11"):
					fh2 = open("polls/data/biogrid_eco.txt",encoding='utf-8')
				elif(ndex_file_id == "12"):
					fh2 = open("polls/data/biogrid_sce.txt",encoding='utf-8')
				elif(ndex_file_id == "13"):
					fh2 = open("polls/data/biogrid_dan.txt",encoding='utf-8')
				ppistr = fh2.read()
				exprstr = request.FILES['myfile'].read().decode('utf-8')
				result4 = check_input_files.delay(ppistr,exprstr)
				errstr = result4.get()
				if(errstr != ""):
					request.session['errors'] = errstr
					return render(request,'polls/errorpage.html',{'errors':errstr})
				savedata_param = "false"
				plot_div = ""
				clinicalstr = request.FILES['patientdata'].read().decode('utf-8')
				clinical_stringio = StringIO(clinicalstr)
				clinicaldf = pd.read_csv(clinical_stringio)
				survival_col_name = ""
				if('survival_col' in request.POST):
					print("barabsfrbasdb")
					if(request.POST['survival_col']):
						survival_col_name = request.POST['survival_col']
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			savedata_param = "true"
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				#clinicalstr = "empty"
				#if(request.user.is_authenticated and savedata_param == "true"):
				#	result2 = script_output_task_2.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,savedata_param,str(request.user))
				#	(div,script,plot1) = result2.get()			
				#else:				
				#result2 = script_output_task_4.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2)
				#(div,script,plot1,plot_div,ret_metadata) = result2.get()
				result2 = script_output_task_9.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf)
				(div,script,plot1,plot_div,ret_metadata) = result2.get()	
				plot2 = "test.png"
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			username = str(request.user)
	                			foobar = "user_uploaded_files/" + username
	                			if not(os.path.isdir(foobar)):
	                				os.mkdir(foobar)
	                			filename_1 = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S_%f')[:-3]
	                			filename_2 = foobar + "/" + username + filename_1 + "expr.txt"
	                			filename_3 = foobar + "/" + username + filename_1 + "prot.txt"
	                			with open(filename_2, "w") as text_file:
	                				text_file.write(exprstr)
	                				text_file.close()
	                			with open(filename_3, "w") as text_file_2:
	                				text_file_2.write(ppistr)
	                				text_file_2.close()
	                			#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
						#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
				#else:
					#cache.clear()
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)      	        
				#else:
				remove_loading_image.delay()
				cache.clear()				
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				return render(request, 'polls/clustering_6.html', {'list_of_files':list_of_files,'ret_metadata':ret_metadata,'result_99':ret_metadata,'list_of_files_2':list_of_files_2,'json_path':"test15.json",'heatmap_path':"test.png",'survival_graph':plot_div})
	
	
	
	
	elif('own_file' in request.POST and request.user.is_authenticated):
		if(request.POST['own_file']):
			make_empty_figure.delay()
			with open("polls/static/output_console.txt", "w") as text_file:
   				text_file.write("Your request is being processed...")
			filename1 = request.POST.get("input_own_file")
			#fh1 = open(filename1)
			#exprstr = fh1.read()
			#filename2 = filename1.split("_expr.txt")[0] + "_prot.txt"
			#fh2 = open(filename2)
			path_json = filename1
			path_heatmap = filename1.split("_json.json")[0] + "_heatmap.png"
			#ppistr = fh2.read()
			list_of_files = GraphForm.list_user_data_2(username)	
			list_of_files_2 = GraphForm.list_user_data(username)      
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])	
				show_old_data(path_json,path_heatmap)             	
				plot2 = "test.png"
				cache.clear()			
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				return render(request, 'polls/clustering_6.html', {'form':"",'images':"",'div':"",'script':"",'plot2':plot2, 'list_of_files':list_of_files,'list_of_files_2':list_of_files_2})       		
	elif('analyze_metadata' in request.POST and 'myfile' in request.FILES and 'protfile' in request.FILES and 'patientdata' in request.FILES):
		if(request.FILES['myfile'] and request.FILES['protfile'] and request.FILES['patientdata']):
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				add_loading_image.delay()
				with open("polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				make_empty_figure.delay()
				#fh2 = open(request.FILES['protfile'])
				ppistr = request.FILES['protfile'].read().decode('utf-8')
				#fh1 = open(request.FILES['myfile'])
				exprstr = request.FILES['myfile'].read().decode('utf-8')
				clinicalstr = request.FILES['patientdata'].read().decode('utf-8')
				clinical_stringio = StringIO(clinicalstr)
				clinicaldf = pd.read_csv(clinical_stringio)
				print(clinicaldf)
				result4 = check_input_files.delay(ppistr,exprstr)
				errstr = result4.get()
				if(errstr != ""):
					request.session['errors'] = errstr
					return render(request,'polls/errorpage.html',{'errors':errstr})
				survival_col_name = ""
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			savedata_param = "true"
				if('survival_col' in request.POST):
					print("barabsfrbasdb")
					if(request.POST['survival_col']):
						survival_col_name = request.POST['survival_col']
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				#result2 = script_output_task.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1)
				#result2 = script_output_task_4.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2)
				#result2 = script_output_task_8.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name)
				#(div,script,plot1,plot_div,ret_metadata) = result2.get()
				#result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				#(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				#print("ikkrngasndgnn")
				#result2 = script_output_task_4.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2)
				result2 = script_output_task_9.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf)
				(div,script,plot1,plot_div,ret_metadata) = result2.get()
				metd = list_metadata_3.apply_async(countdown=0)
				(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get()
				print(adjlist)
				convert_gene_list.delay(adjlist,"polls/static/genelist_temp.txt")
				print("iasdfasdfsf")
				print(ret_metadata1) 	
				#(div,script,plot1) = result2.get()			
				plot2 = "test.png"
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			savedata_param = "true"
	                			username = str(request.user)
	                			GraphForm.save_user_data_2(request.FILES['myfile'],request.FILES['protfile'],request.FILES['patientdata'],username)
						#GraphForm.save_results(username)
				#else:
					#cache.clear()
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)              
				#else:
				remove_loading_image.delay()
				cache.clear()				
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				#list_metadata.apply_async(countdown=0)
				with open("polls/static/plotly_output_2.html", "w") as text_file_3:
   					text_file_3.write(plot_div)
				return render(request, 'polls/clustering_6.html', {'form':"",'images':"",'plot_div':plot_div,'script':script,'plot2':plot2, 'list_of_files':list_of_files,'ret_dat':ret_metadata,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'list_of_files_2':list_of_files_2})

	elif('myfile' in request.FILES and 'protfile' in request.FILES):
		if(request.FILES['myfile'] and request.FILES['protfile']):
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				add_loading_image.delay()
				with open("polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				make_empty_figure.delay()
				#fh2 = open(request.FILES['protfile'])
				ppistr = request.FILES['protfile'].read().decode('utf-8')
				#fh1 = open(request.FILES['myfile'])
				exprstr = request.FILES['myfile'].read().decode('utf-8')
				savedata_param = "false"
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			savedata_param = "true"
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				#result2 = script_output_task.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1)
				clinicalstr = "empty"
				result4 = check_input_files.delay(ppistr,exprstr)
				errstr = result4.get()
				if(errstr != ""):
					request.session['errors'] = errstr
					return render(request,'polls/errorpage.html',{'errors':errstr})
				if(request.user.is_authenticated and savedata_param == "true"):
					result2 = script_output_task_2.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,savedata_param,str(request.user))
					(div,script,plot1) = result2.get()			
				else:				
					clinicalstr = "empty"
					clinicaldf = ""
					survival_col_name = ""
					result2 = script_output_task_9.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf)
					(div,script,plot1,plot_div,ret_metadata) = result2.get()	
				
					#result2 = script_output_task_4.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2)
					#(div,script,plot1,plot_div,ret_metadata) = result2.get()
				plot2 = "test.png"
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			savedata_param = "true"
	                			username = str(request.user)
	                			foobar = "user_uploaded_files/" + username
	                			if not(os.path.isdir(foobar)):
	                				os.mkdir(foobar)
	                			filename_1 = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S_%f')[:-3]
	                			filename_2 = foobar + "/" + username + filename_1 + "expr.txt"
	                			filename_3 = foobar + "/" + username + filename_1 + "prot.txt"
	                			with open(filename_2, "w") as text_file:
	                				text_file.write(exprstr)
	                				text_file.close()
	                			with open(filename_3, "w") as text_file_2:
	                				text_file_2.write(ppistr)
	                				text_file_2.close()
	                		
	                			#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
				#else:
					#cache.clear()
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)      	        
				#else:
				remove_loading_image.delay()
				cache.clear()				
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				return render(request, 'polls/clustering_6.html', {'form':"",'images':"",'div':div,'script':script,'plot2':plot2, 'list_of_files':list_of_files,'list_of_files_2':list_of_files_2})
	elif('predef_file' in request.POST and 'protfile' in request.POST):
		if('L_g_min' in request.POST and 'L_g_max' in request.POST):
			if(request.POST['predef_file'] and request.POST['protfile']):
				make_empty_figure.delay()
				with open("polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				add_loading_image.delay()
				cancer_type = request.POST.get("cancer_type")
				print("Cancer Type:" + str(cancer_type))
				ppistr = request.FILES['protfile'].read().decode('utf-8')
				#fh2 = open("polls/data/biogrid.human.entrez.tsv")
				#ppistr = fh2.read()
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				is_lungc = "false"					
				if(cancer_type == "1"):
					print("babababababa")
					fh1 = open("polls/data/lung_cancer_expr.csv")
					exprstr = fh1.read()
					clinicaldf = pd.read_csv("polls/data/lung_cancer_clinical.csv")
					fh4 = open("polls/data/lung_cancer_clinical.csv")
					#fh4 = open("polls/data/breast_cancer_clinical.csv")
					is_lungc = "true"
					clinicalstr = fh4.read()
					fh4.flush()
					fh4.close()
					survival_col_name = "disease free survival in months:ch1"
				else:
					fh1 = open("polls/data/breast_cancer_expr.csv")
					exprstr = fh1.read()
					clinicaldf = pd.read_csv("polls/data/breast_cancer_clinical.csv")
					fh4 = open("polls/data/breast_cancer_clinical.csv")
					clinicalstr = fh4.read()
					fh4.flush()
					fh4.close()
					survival_col_name = "mfs (yr):ch1"
				print(clinicaldf)
				if('survival_col' in request.POST):
					print("barabsfrbasdb")
					if(request.POST['survival_col']):
						survival_col_name = request.POST['survival_col']
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				print("ikkrngasndgnn")
				#result2 = script_output_task_4.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2)
				result2 = script_output_task_9.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf)
				(div,script,plot1,plot_div,ret_metadata) = result2.get()
				print(ret_metadata)
				metd = list_metadata_3.apply_async(countdown=0)
				(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get()
				plot2 = "test.png"
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)        				
				#else:
				remove_loading_image.delay()
				with open("polls/static/metadata_99.html", "w") as text_file_4:
   					text_file_4.write(str(ret_metadata[0]))
   					text_file_4.write(str(ret_metadata[1]))
   					text_file_4.write(str(ret_metadata[2]))
				cache.clear()
				with open("polls/static/plotly_output_2.html", "w") as text_file_3:
   					text_file_3.write(plot_div)
				#write_metadata_to_file.apply_async((ret_metadata), countdown=1)
				#list_metadata.apply_async(countdown=0)
				string_99 = metadata_to_string.apply_async(ret_metadata, countdown=0)
				result_99 = string_99.get()
				return render(request, 'polls/clustering_6.html', {'plot2':plot2, 'list_of_files':list_of_files,'ret_metadata':ret_metadata,'result_99':ret_metadata,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'list_of_files_2':list_of_files_2})
	
	elif('predef_file' in request.POST):
		if('L_g_min' in request.POST and 'L_g_max' in request.POST):
			if(request.POST['predef_file']):
				make_empty_figure.delay()
				with open("polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				add_loading_image.delay()
				cancer_type = request.POST.get("cancer_type")
				print("Cancer Type:" + str(cancer_type))
				fh2 = open("polls/data/biogrid.human.entrez.tsv")
				ppistr = fh2.read()
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				is_lungc = "false"					
				if(cancer_type == "1"):
					print("babababababa")
					fh1 = open("polls/data/lung_cancer_expr.csv")
					exprstr = fh1.read()
					clinicaldf = pd.read_csv("polls/data/lung_cancer_clinical.csv")
					fh4 = open("polls/data/lung_cancer_clinical.csv")
					#fh4 = open("polls/data/breast_cancer_clinical.csv")
					is_lungc = "true"
					clinicalstr = fh4.read()
					fh4.flush()
					fh4.close()
					survival_col_name = "disease free survival in months:ch1"
				else:
					fh1 = open("polls/data/breast_cancer_expr.csv")
					exprstr = fh1.read()
					clinicaldf = pd.read_csv("polls/data/breast_cancer_clinical.csv")
					fh4 = open("polls/data/breast_cancer_clinical.csv")
					clinicalstr = fh4.read()
					fh4.flush()
					fh4.close()
					survival_col_name = "mfs (yr):ch1"
				print(clinicaldf)
				if('survival_col' in request.POST):
					print("barabsfrbasdb")
					if(request.POST['survival_col']):
						survival_col_name = request.POST['survival_col']
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				print("ikkrngasndgnn")
				#result2 = script_output_task_4.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2)
				result2 = script_output_task_9.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf)
				(div,script,plot1,plot_div,ret_metadata) = result2.get()
				print(ret_metadata)
				metd = list_metadata_3.apply_async(countdown=0)
				(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get()
				plot2 = "test.png"
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)        				
				#else:
				remove_loading_image.delay()
				with open("polls/static/metadata_99.html", "w") as text_file_4:
   					text_file_4.write(str(ret_metadata[0]))
   					text_file_4.write(str(ret_metadata[1]))
   					text_file_4.write(str(ret_metadata[2]))
				cache.clear()
				with open("polls/static/plotly_output_2.html", "w") as text_file_3:
   					text_file_3.write(plot_div)
				#write_metadata_to_file.apply_async((ret_metadata), countdown=1)
				#list_metadata.apply_async(countdown=0)
				string_99 = metadata_to_string.apply_async(ret_metadata, countdown=0)
				result_99 = string_99.get()
				return render(request, 'polls/clustering_6.html', {'plot2':plot2, 'list_of_files':list_of_files,'ret_metadata':ret_metadata,'result_99':ret_metadata,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'list_of_files_2':list_of_files_2})
	elif('kegg_enrichment' in request.POST):
		group_for_enr = "both"
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		if('group_for_enr' in request.POST):
			group_for_enr = request.POST.get('group_for_enr')
			#print(pval_enr)
		result1 = run_enrichment_2.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_kegg")
		result2 = run_enrichment_2.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_kegg")
		result3 = run_enrichment_2.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_kegg")
		enr_results = result1.get()
		enr_results_2 = result2.get()
		enr_results_3 = result3.get()
		print("enr")
		result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
		result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
		result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
		enrichment_dict = result4.get()
		enrichment_dict_2 = result5.get()
		enrichment_dict_3 = result6.get()
		return render(request,'polls/clustering_6.html',{'list_of_files':list_of_files,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3})
	elif('enrichment_type' in request.POST):
		enr_type = request.POST.get("enrichment_type")
		group_for_enr = "both"
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		if('group_for_enr' in request.POST):
			group_for_enr = request.POST.get('group_for_enr')
			#print(pval_enr)
		enrichment_dict = {}
		enrichment_dict_2 = {}
		enrichment_dict_3 = {}
		enrichment_dict_4 = {}
		enrichment_dict_5 = {}
		if(enr_type == "kegg_enrichment"):
			result1 = run_enrichment_2.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_kegg")
			result2 = run_enrichment_2.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_kegg")
			result3 = run_enrichment_2.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_kegg")
			enr_results = result1.get()
			enr_results_2 = result2.get()
			enr_results_3 = result3.get()
			print("enr")
			result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			result7 = read_kegg_enrichment_2.delay("polls/data/test2/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt","polls/data/test3/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			enrichment_dict = result4.get()
			enrichment_dict_2 = result5.get()
			enrichment_dict_3 = result6.get()
			(enrichment_dict_4,enrichment_dict_5) = result7.get()
		elif(enr_type == "go_enrichment"):	
			result1 = run_go_enrichment_2.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_go")
			result2 = run_go_enrichment_2.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_go")
			result3 = run_go_enrichment_2.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_go")
			#result1 = run_go_enrichment.delay("genelist.txt",pval_enr)
			enr_results = result1.get()
			enr_results_2 = result2.get()
			enr_results_3 = result3.get()
			print("enr")
			result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			result7 = read_kegg_enrichment_2.delay("polls/data/test2/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt","polls/data/test3/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			enrichment_dict = result4.get()
			enrichment_dict_2 = result5.get()
			enrichment_dict_3 = result6.get()	
			(enrichment_dict_4,enrichment_dict_5) = result7.get()
		elif(enr_type == "go_molecular"):
			result1 = run_go_enrichment_2.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_go")
			result2 = run_go_enrichment_2.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_go")
			result3 = run_go_enrichment_2.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_go")
			enr_results = result1.get()
			enr_results_2 = result2.get()
			enr_results_3 = result3.get()
			print("enr")
			result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result7 = read_kegg_enrichment_2.delay("polls/data/test2/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt","polls/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			enrichment_dict = result4.get()
			enrichment_dict_2 = result5.get()
			enrichment_dict_3 = result6.get()
			(enrichment_dict_4,enrichment_dict_5) = result7.get()
		elif(enr_type == "reactome_enrichment"):
			#result1 = run_reac_enrichment.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_reactome")
			result2 = run_reac_enrichment.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_reactome")
			#result3 = run_reac_enrichment.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_reactome")
			#enr_results = result1.get()
			enr_results_2 = result2.get()
			#enr_results_3 = result3.get()
			print("enr")
			#result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_reactome/Reactome_2016.test_name.enrichr.reports.txt",pval_enr)
			#result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			#enrichment_dict = result4.get()
			enrichment_dict = {}
			enrichment_dict_2 = result5.get()
			#enrichment_dict_3 = result6.get()
			enrichment_dict_3 = {}		
		return render(request,'polls/clustering_6.html',{'list_of_files':list_of_files,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3,'enrichment_dict_4':enrichment_dict_4,'enrichment_dict_5':enrichment_dict_5})

	elif('go_enrichment' in request.POST):
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		group_for_enr = "both"
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		if('group_for_enr' in request.POST):
			group_for_enr = request.POST.get('group_for_enr')
			#print(pval_enr)
		result1 = run_go_enrichment_2.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_go")
		result2 = run_go_enrichment_2.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_go")
		result3 = run_go_enrichment_2.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_go")
		#result1 = run_go_enrichment.delay("genelist.txt",pval_enr)
		enr_results = result1.get()
		enr_results_2 = result2.get()
		enr_results_3 = result3.get()
		print("enr")
		result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
		result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
		result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
		enrichment_dict = result4.get()
		enrichment_dict_2 = result5.get()
		enrichment_dict_3 = result6.get()
		return render(request,'polls/clustering_6.html',{'list_of_files':list_of_files,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3})
	elif('go_molecular' in request.POST):
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		#result1 = run_go_enrichment.delay("genelist.txt",pval_enr)
		group_for_enr = "both"
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		if('group_for_enr' in request.POST):
			group_for_enr = request.POST.get('group_for_enr')
			#print(pval_enr)
		result1 = run_go_enrichment_2.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_go")
		result2 = run_go_enrichment_2.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_go")
		result3 = run_go_enrichment_2.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_go")
		enr_results = result1.get()
		enr_results_2 = result2.get()
		enr_results_3 = result3.get()
		print("enr")
		result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
		result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
		result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
		enrichment_dict = result4.get()
		enrichment_dict_2 = result5.get()
		enrichment_dict_3 = result6.get()
		return render(request,'polls/clustering_6.html',{'list_of_files':list_of_files,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3})
	elif('reac_enrichment' in request.POST):
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		#result1 = run_go_enrichment.delay("genelist.txt",pval_enr)
		group_for_enr = "both"
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		if('group_for_enr' in request.POST):
			group_for_enr = request.POST.get('group_for_enr')
			#print(pval_enr)
		#result1 = run_reac_enrichment.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_reactome")
		result2 = run_reac_enrichment.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_reactome")
		#result3 = run_reac_enrichment.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_reactome")
		#enr_results = result1.get()
		enr_results_2 = result2.get()
		#enr_results_3 = result3.get()
		print("enr")
		#result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
		result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_reactome/Reactome_2016.test_name.enrichr.reports.txt",pval_enr)
		#result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
		#enrichment_dict = result4.get()
		enrichment_dict = {}
		enrichment_dict_2 = result5.get()
		#enrichment_dict_3 = result6.get()
		enrichment_dict_3 = {}
		return render(request,'polls/clustering_6.html',{'list_of_files':list_of_files,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3})

	#if not(('myfile' in request.FILES and 'protfile' in request.FILES) or ('predef_file' in request.POST)):
	else:		
		#remove_loading_image.delay()
		ret_metadata = ""
		if not (request.user.is_authenticated):
			cache.clear()
		else:
			username = str(request.user)
			list_of_files = GraphForm.list_user_data_2(username)
			list_of_files_2 = GraphForm.list_user_data(username)
			cache.clear()
			#metd = list_metadata.apply_async(countdown=0)
			metd = list_metadata_3.apply_async(countdown=0)
			(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get() 
			print(ret_metadata1) 
			print("iubaerb")
			metadata_dict = [ret_metadata1,ret_metadata2,ret_metadata3]
			result2 = read_kegg_enrichment.delay("polls/data/test/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			#result2 = read_kegg_enrichment.delay("data/test/enrichr_kegg")	
			enrichment_dict = result2.get()
			#(ret_metadata,is_lungc) = metd.get()  
			#return render(request,'polls/clustering_6.html',{'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'ret_metadata':ret_metadata,'is_lungc':is_lungc})

		#make_empty_figure.apply_async(countdown=2)
		#empty_log_file.apply_async(countdown=2)	 				
		return render(request,'polls/clustering_6.html',{'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'ret_metadata':ret_metadata,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'metadata_dict':metadata_dict,'enrichment_dict':enrichment_dict})


def clustering_6_part_1(request):
	if('newAnalysis' in request.POST):
		print(request.POST.get("newAnalysis"))
		print(request.POST['newAnalysis'])
	#	request.POST._mutable = True
	#	if(request.POST['newAnalysis'] != "false"):
	#		if('done' in request.session):
	#			if(request.session['done'] == "true"):
	#				request.session['done'] = "False"
	#		request.POST['newAnalysis'] = "false"
	#if(request.session['done'] == "true"):
	if('done' in request.session):
		#return HttpResponseRedirect('polls/clustering_6_part_3.html')
		print("done")
	ret_metadata1 = {}
	ret_metadata2 = {}
	ret_metadata3 = {}
	metadata_dict = []
	enrichment_dict = []
	pval_enr = 0.5
	list_of_files = ""
	list_of_files_2 = ""
	save_data = request.POST.get("save_data", None)
	nbr_iter = request.POST.get("nbr_iter",45)
	nbr_ants = request.POST.get("nbr_ants",30)
	evap = request.POST.get("evap",0.3)
	epsilon = request.POST.get("stopcr",0.02) 
	hi_sig = request.POST.get("hisig",1)
	pher_sig = request.POST.get("pher",1)
	savedata_param = "false"
	#fig = plt.figure(figsize=(10,8))
	#plt.savefig("polls/static/progress.png")
	#plt.close(fig)
	if request.user.is_authenticated:
		username = str(request.user)
		list_of_files = GraphForm.list_user_data_2(username)
	if('redo_analysis' in request.POST and request.user.is_authenticated):
		if(request.POST['redo_analysis']):
			with open("polls/static/output_console.txt", "w") as text_file:
   				text_file.write("Your request is being processed...")
			add_loading_image.delay()
			filename1 = request.POST.get("input_own_file_redo")
			fh1 = open(filename1)
			exprstr = fh1.read()
			filename2 = filename1.split("_expr.txt")[0] + "_prot.txt"
			fh2 = open(filename2)
			filename3 = filename1.split("_expr.txt")[0] + "_clin.txt"
			has_clin_data = "false"
			if(os.path.isfile(filename3)):
				fh3 = open(filename3)
				has_clin_data = "true"
				clinicalstr = fh3.read()
			#path_json = filename1
			#path_heatmap = filename1.split("_json.json")[0] + "_heatmap.png"
			ppistr = fh2.read()
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				make_empty_figure.delay()
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])	
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids) =result1.get()				
				if(has_clin_data == "true"):
					result2 = script_output_task_3.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr)
					(div,script,plot1,plot_div,ret_metadata) = result2.get()
				else:
					result2 = script_output_task.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1)
					(div,script,plot1) = result2.get()
					ret_metadata = ""
				plot2 = "test.png"
				cache.clear()			
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				list_of_files = GraphForm.list_user_data_2(username)	
				list_of_files_2 = GraphForm.list_user_data(username)      
				remove_loading_image.delay()
				return render(request, 'polls/clustering_6_part_1.html', {'form':"",'images':"",'div':"",'script':"",'plot2':plot2, 'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'ret_dat':ret_metadata})       		
	elif('myfile' in request.FILES and 'ndexfile' in request.FILES and 'protfile' in request.FILES):
		if(request.FILES['myfile'] and request.FILES['ndexfile'] and request.FILES['protfile']):
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				add_loading_image.delay()
				with open("polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				make_empty_figure.delay()
				#fh2 = open(request.FILES['protfile'])
				result3 = read_ndex_file_2.delay(request.FILES['ndexfile'].read().decode('utf-8'))
				ppistr_1 = result3.get()
				ppistr_2 = request.FILES['protfile'].read().decode('utf-8')
				ppistr = "\n".join([ppistr_2,ppistr_1])
				#fh1 = open(request.FILES['myfile'])
				exprstr = request.FILES['myfile'].read().decode('utf-8')
				savedata_param = "false"
				plot_div = ""
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			savedata_param = "true"
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				clinicalstr = "empty"
				if(request.user.is_authenticated and savedata_param == "true"):
					result2 = script_output_task_2.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,savedata_param,str(request.user))
					(div,script,plot1) = result2.get()			
				else:				
					result2 = script_output_task_4.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2)
					(div,script,plot1,plot_div,ret_metadata) = result2.get()
				plot2 = "test.png"
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			username = str(request.user)
	                			foobar = "user_uploaded_files/" + username
	                			if not(os.path.isdir(foobar)):
	                				os.mkdir(foobar)
	                			filename_1 = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S_%f')[:-3]
	                			filename_2 = foobar + "/" + username + filename_1 + "expr.txt"
	                			filename_3 = foobar + "/" + username + filename_1 + "prot.txt"
	                			with open(filename_2, "w") as text_file:
	                				text_file.write(exprstr)
	                				text_file.close()
	                			with open(filename_3, "w") as text_file_2:
	                				text_file_2.write(ppistr)
	                				text_file_2.close()
	                			#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
				#else:
					#cache.clear()
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)      	        
				#else:
				remove_loading_image.delay()
				cache.clear()				
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				return render(request, 'polls/clustering_6_part_1.html', {'list_of_files':list_of_files,'ret_metadata':ret_metadata,'result_99':ret_metadata,'list_of_files_2':list_of_files_2,'json_path':"test15.json",'heatmap_path':"test.png",'survival_graph':plot_div})
	elif('myfile' in request.FILES and 'ndex_file' in request.POST and 'ndex_name' in request.POST):
		print(request.POST.get("ndex_name"))
		if(request.FILES['myfile'] and request.POST['ndex_file'] and request.POST['ndex_name']):
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				add_loading_image.delay()
				with open("polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				make_empty_figure.delay()
				ndex_file_id = request.POST.get("ndex_name")
				if(ndex_file_id == "1"):
					fh2 = open("polls/data/APID.txt",encoding='utf-8')
				elif(ndex_file_id == "2"):
					fh2 = open("polls/data/STRING.txt",encoding='utf-8')
				elif(ndex_file_id == "3"):
					fh2 = open("polls/data/biogrid.human.entrez.tsv",encoding='utf-8')
				elif(ndex_file_id == "4"):
					fh2 = open("polls/data/HPRD.txt",encoding='utf-8')
				elif(ndex_file_id == "5"):
					fh2 = open("polls/data/ulitsky.txt",encoding='utf-8')
				elif(ndex_file_id == "6"):
					fh2 = open("polls/data/I2D.txt",encoding='utf-8')
				elif(ndex_file_id == "7"):
					fh2 = open("polls/data/biogrid_cel.txt",encoding='utf-8')
				elif(ndex_file_id == "8"):
					fh2 = open("polls/data/biogrid_ath.txt",encoding='utf-8')
				elif(ndex_file_id == "9"):
					fh2 = open("polls/data/biogrid_mus.txt",encoding='utf-8')
				elif(ndex_file_id == "10"):
					fh2 = open("polls/data/biogrid_dme.txt",encoding='utf-8')
				elif(ndex_file_id == "11"):
					fh2 = open("polls/data/biogrid_eco.txt",encoding='utf-8')
				elif(ndex_file_id == "12"):
					fh2 = open("polls/data/biogrid_sce.txt",encoding='utf-8')
				elif(ndex_file_id == "13"):
					fh2 = open("polls/data/biogrid_dan.txt",encoding='utf-8')
				ppistr = fh2.read()
				exprstr = request.FILES['myfile'].read().decode('utf-8')
				savedata_param = "false"
				plot_div = ""
				#if save_data in ["save_data"]:
	                	#	if request.user.is_authenticated:
	                	#		savedata_param = "true"
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				clinicalstr = "empty"
				if(request.user.is_authenticated and savedata_param == "true"):
					result2 = script_output_task_2.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,savedata_param,str(request.user))
					(div,script,plot1) = result2.get()			
				else:				
					result2 = script_output_task_4.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2)
					(div,script,plot1,plot_div,ret_metadata) = result2.get()
				plot2 = "test.png"
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			username = str(request.user)
	                			foobar = "user_uploaded_files/" + username
	                			if not(os.path.isdir(foobar)):
	                				os.mkdir(foobar)
	                			filename_1 = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S_%f')[:-3]
	                			filename_2 = foobar + "/" + username + filename_1 + "expr.txt"
	                			filename_3 = foobar + "/" + username + filename_1 + "prot.txt"
	                			with open(filename_2, "w") as text_file:
	                				text_file.write(exprstr)
	                				text_file.close()
	                			with open(filename_3, "w") as text_file_2:
	                				text_file_2.write(ppistr)
	                				text_file_2.close()
	                			#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
						#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
				#else:
					#cache.clear()
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)      	        
				#else:
				remove_loading_image.delay()
				cache.clear()				
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				return render(request, 'polls/clustering_6_part_1.html', {'list_of_files':list_of_files,'ret_metadata':ret_metadata,'result_99':ret_metadata,'list_of_files_2':list_of_files_2,'json_path':"test15.json",'heatmap_path':"test.png",'survival_graph':plot_div})
	elif('myfile' in request.FILES and 'parse_ndex_file' in request.POST and 'ndex_name_2' in request.POST):
		#print(request.POST.get("ndex_name"))
		if(request.FILES['myfile'] and request.POST['parse_ndex_file']):
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				add_loading_image.delay()
				with open("polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				make_empty_figure.delay()
				ndex_file_id = request.POST.get("ndex_name_2")
				if(ndex_file_id == "1"):
					result_ndex = import_ndex.delay("9c38ce6e-c564-11e8-aaa6-0ac135e8bacf")
					ppistr = result_ndex.get()
				elif(ndex_file_id == "2"):
					result_ndex = import_ndex.delay("275bd84e-3d18-11e8-a935-0ac135e8bacf")
					ppistr = result_ndex.get()
				elif(ndex_file_id == "3"):
					result_ndex = import_ndex.delay("becec556-86d4-11e7-a10d-0ac135e8bacf")
					ppistr = result_ndex.get()
				elif(ndex_file_id == "4"):
					result_ndex = import_ndex.delay("1093e665-86da-11e7-a10d-0ac135e8bacf")
					ppistr = result_ndex.get()
				#if(ndex_file_id == "1"):
				#	fh2 = open("test_scripts_3/APID.cx",encoding='utf-8')
				#elif(ndex_file_id == "2"):
				#	fh2 = open("test_scripts_3/STRING.cx",encoding='utf-8')
				#elif(ndex_file_id == "3"):
				#	fh2 = open("polls/data/biogrid.human.entrez.tsv",encoding='utf-8')
				#	ppistr = fh2.read()
				#elif(ndex_file_id == "4"):
				#	fh2 = open("test_scripts_3/HPRD.cx",encoding='utf-8')
				elif(ndex_file_id == "5"):
					fh2 = open("polls/data/ulitsky.txt",encoding='utf-8')
					ppistr = fh2.read()
				elif(ndex_file_id == "6"):
					fh2 = open("polls/data/I2D.txt",encoding='utf-8')
					ppistr = fh2.read()
				elif(ndex_file_id == "7"):
					fh2 = open("test_scripts_3/biogrid_cel.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "8"):
					fh2 = open("test_scripts_3/biogrid_ath.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "9"):
					fh2 = open("test_scripts_3/biogrid_mus.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "10"):
					fh2 = open("test_scripts_3/biogrid_dro.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "11"):
					fh2 = open("test_scripts_3/biogrid_eco.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "12"):
					fh2 = open("test_scripts_3/biogrid_sce.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "13"):
					fh2 = open("test_scripts_3/biogrid_dan.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				#if not(ndex_file_id == "1" or ndex_file_id == "2" or ndex_file_id == "3" or ndex_file_id=="6" or ndex_file_id == "5"):
				exprstr = request.FILES['myfile'].read().decode('utf-8')
				savedata_param = "false"
				plot_div = ""
				#if save_data in ["save_data"]:
	                	#	if request.user.is_authenticated:
	                	#		savedata_param = "true"
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				clinicalstr = "empty"
				if(request.user.is_authenticated and savedata_param == "true"):
					result2 = script_output_task_2.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,savedata_param,str(request.user))
					(div,script,plot1) = result2.get()			
				else:				
					clinicalstr = "empty"
					clinicaldf = ""
					survival_col_name = ""
					result2 = script_output_task_9.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf)
					(div,script,plot1,plot_div,ret_metadata) = result2.get()	
				
					#result2 = script_output_task_4.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2)
					#(div,script,plot1,plot_div,ret_metadata) = result2.get()
				plot2 = "test.png"
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			username = str(request.user)
	                			foobar = "user_uploaded_files/" + username
	                			if not(os.path.isdir(foobar)):
	                				os.mkdir(foobar)
	                			filename_1 = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S_%f')[:-3]
	                			filename_2 = foobar + "/" + username + filename_1 + "expr.txt"
	                			filename_3 = foobar + "/" + username + filename_1 + "prot.txt"
	                			with open(filename_2, "w") as text_file:
	                				text_file.write(exprstr)
	                				text_file.close()
	                			with open(filename_3, "w") as text_file_2:
	                				text_file_2.write(ppistr)
	                				text_file_2.close()
	                			#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
						#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
				#else:
					#cache.clear()
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)      	        
				#else:
				remove_loading_image.delay()
				cache.clear()				
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				return render(request, 'polls/clustering_6_part_1.html', {'list_of_files':list_of_files,'ret_metadata':ret_metadata,'result_99':ret_metadata,'list_of_files_2':list_of_files_2,'json_path':"test15.json",'heatmap_path':"test.png",'survival_graph':plot_div})
	elif('predef_file' in request.FILES and 'parse_ndex_file' in request.POST and 'ndex_name_2' in request.POST):
		#print(request.POST.get("ndex_name"))
		if(request.FILES['myfile'] and request.POST['parse_ndex_file'] and request.POST['predef_file']):
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				add_loading_image.delay()
				cancer_type = request.POST.get("cancer_type")
				if(cancer_type == "1"):
					print("babababababa")
					fh1 = open("polls/data/lung_cancer_expr.csv")
					exprstr = fh1.read()
					clinicaldf = pd.read_csv("polls/data/lung_cancer_clinical.csv")
					fh4 = open("polls/data/lung_cancer_clinical.csv")
					#fh4 = open("polls/data/breast_cancer_clinical.csv")
					is_lungc = "true"
					clinicalstr = fh4.read()
					fh4.flush()
					fh4.close()
					survival_col_name = "disease free survival in months:ch1"
				else:
					fh1 = open("polls/data/breast_cancer_expr.csv")
					exprstr = fh1.read()
					clinicaldf = pd.read_csv("polls/data/breast_cancer_clinical.csv")
					fh4 = open("polls/data/breast_cancer_clinical.csv")
					clinicalstr = fh4.read()
					fh4.flush()
					fh4.close()
					survival_col_name = "mfs (yr):ch1"
				
				with open("polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				make_empty_figure.delay()
				ndex_file_id = request.POST.get("ndex_name_2")
				if(ndex_file_id == "1"):
					result_ndex = import_ndex.delay("9c38ce6e-c564-11e8-aaa6-0ac135e8bacf")
					ppistr = result_ndex.get()
				elif(ndex_file_id == "2"):
					result_ndex = import_ndex.delay("275bd84e-3d18-11e8-a935-0ac135e8bacf")
					ppistr = result_ndex.get()
				elif(ndex_file_id == "3"):
					result_ndex = import_ndex.delay("becec556-86d4-11e7-a10d-0ac135e8bacf")
					ppistr = result_ndex.get()
				elif(ndex_file_id == "4"):
					result_ndex = import_ndex.delay("1093e665-86da-11e7-a10d-0ac135e8bacf")
					ppistr = result_ndex.get()
				#if(ndex_file_id == "1"):
				#	fh2 = open("test_scripts_3/APID.cx",encoding='utf-8')
				#elif(ndex_file_id == "2"):
				#	fh2 = open("test_scripts_3/STRING.cx",encoding='utf-8')
				#elif(ndex_file_id == "3"):
				#	fh2 = open("polls/data/biogrid.human.entrez.tsv",encoding='utf-8')
				#	ppistr = fh2.read()
				#elif(ndex_file_id == "4"):
				#	fh2 = open("test_scripts_3/HPRD.cx",encoding='utf-8')
				elif(ndex_file_id == "5"):
					fh2 = open("polls/data/ulitsky.txt",encoding='utf-8')
					ppistr = fh2.read()
				elif(ndex_file_id == "6"):
					fh2 = open("polls/data/I2D.txt",encoding='utf-8')
					ppistr = fh2.read()
				elif(ndex_file_id == "7"):
					fh2 = open("test_scripts_3/biogrid_cel.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "8"):
					fh2 = open("test_scripts_3/biogrid_ath.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "9"):
					fh2 = open("test_scripts_3/biogrid_mus.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "10"):
					fh2 = open("test_scripts_3/biogrid_dro.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "11"):
					fh2 = open("test_scripts_3/biogrid_eco.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "12"):
					fh2 = open("test_scripts_3/biogrid_sce.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "13"):
					fh2 = open("test_scripts_3/biogrid_dan.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				#if not(ndex_file_id == "1" or ndex_file_id == "2" or ndex_file_id == "3" or ndex_file_id=="6" or ndex_file_id == "5"):
				exprstr = request.FILES['myfile'].read().decode('utf-8')
				savedata_param = "false"
				plot_div = ""
				#if save_data in ["save_data"]:
	                	#	if request.user.is_authenticated:
	                	#		savedata_param = "true"
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				clinicalstr = "empty"
				if(request.user.is_authenticated and savedata_param == "true"):
					result2 = script_output_task_2.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,savedata_param,str(request.user))
					(div,script,plot1) = result2.get()			
				else:				
					clinicalstr = "empty"
					clinicaldf = ""
					survival_col_name = ""
					result2 = script_output_task_9.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf)
					(div,script,plot1,plot_div,ret_metadata) = result2.get()	
				
					#result2 = script_output_task_4.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2)
					#(div,script,plot1,plot_div,ret_metadata) = result2.get()
				plot2 = "test.png"
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			username = str(request.user)
	                			foobar = "user_uploaded_files/" + username
	                			if not(os.path.isdir(foobar)):
	                				os.mkdir(foobar)
	                			filename_1 = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S_%f')[:-3]
	                			filename_2 = foobar + "/" + username + filename_1 + "expr.txt"
	                			filename_3 = foobar + "/" + username + filename_1 + "prot.txt"
	                			with open(filename_2, "w") as text_file:
	                				text_file.write(exprstr)
	                				text_file.close()
	                			with open(filename_3, "w") as text_file_2:
	                				text_file_2.write(ppistr)
	                				text_file_2.close()
	                			#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
						#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
				#else:
					#cache.clear()
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)      	        
				#else:
				remove_loading_image.delay()
				cache.clear()				
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				return render(request, 'polls/clustering_6_part_1.html', {'list_of_files':list_of_files,'ret_metadata':ret_metadata,'result_99':ret_metadata,'list_of_files_2':list_of_files_2,'json_path':"test15.json",'heatmap_path':"test.png",'survival_graph':plot_div})

	elif('myfile' in request.FILES and 'parse_ndex_file' in request.POST and 'ndex_name_2' in request.POST and 'analyze_metadata' in request.POST):
		#print(request.POST.get("ndex_name"))
		if(request.FILES['myfile'] and request.POST['parse_ndex_file'] and request.FILES['patientdata']):
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				add_loading_image.delay()
				with open("polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				make_empty_figure.delay()
				ndex_file_id = request.POST.get("ndex_name_2")
				if(ndex_file_id == "1"):
					result_ndex = import_ndex.delay("9c38ce6e-c564-11e8-aaa6-0ac135e8bacf")
					ppistr = result_ndex.get()
				elif(ndex_file_id == "2"):
					result_ndex = import_ndex.delay("275bd84e-3d18-11e8-a935-0ac135e8bacf")
					ppistr = result_ndex.get()
				elif(ndex_file_id == "3"):
					result_ndex = import_ndex.delay("becec556-86d4-11e7-a10d-0ac135e8bacf")
					ppistr = result_ndex.get()
				elif(ndex_file_id == "4"):
					result_ndex = import_ndex.delay("1093e665-86da-11e7-a10d-0ac135e8bacf")
					ppistr = result_ndex.get()
				#if(ndex_file_id == "1"):
				#	fh2 = open("test_scripts_3/APID.cx",encoding='utf-8')
				#elif(ndex_file_id == "2"):
				#	fh2 = open("test_scripts_3/STRING.cx",encoding='utf-8')
				#elif(ndex_file_id == "3"):
				#	fh2 = open("polls/data/biogrid.human.entrez.tsv",encoding='utf-8')
				#	ppistr = fh2.read()
				#elif(ndex_file_id == "4"):
				#	fh2 = open("test_scripts_3/HPRD.cx",encoding='utf-8')
				elif(ndex_file_id == "5"):
					fh2 = open("polls/data/ulitsky.txt",encoding='utf-8')
					ppistr = fh2.read()
				elif(ndex_file_id == "6"):
					fh2 = open("polls/data/I2D.txt",encoding='utf-8')
					ppistr = fh2.read()
				elif(ndex_file_id == "7"):
					fh2 = open("test_scripts_3/biogrid_cel.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "8"):
					fh2 = open("test_scripts_3/biogrid_ath.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "9"):
					fh2 = open("test_scripts_3/biogrid_mus.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "10"):
					fh2 = open("test_scripts_3/biogrid_dro.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "11"):
					fh2 = open("test_scripts_3/biogrid_eco.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "12"):
					fh2 = open("test_scripts_3/biogrid_sce.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				elif(ndex_file_id == "13"):
					fh2 = open("test_scripts_3/biogrid_dan.cx",encoding='utf-8')
					result_ndex = read_ndex_file_4.delay(fh2.read())
					ppistr = result_ndex.get()
				#if not(ndex_file_id == "1" or ndex_file_id == "2" or ndex_file_id == "3" or ndex_file_id=="6" or ndex_file_id == "5"):
				exprstr = request.FILES['myfile'].read().decode('utf-8')
				savedata_param = "false"
				plot_div = ""
				clinicalstr = request.FILES['patientdata'].read().decode('utf-8')
				clinical_stringio = StringIO(clinicalstr)
				clinicaldf = pd.read_csv(clinical_stringio)
				survival_col_name = ""
				if('survival_col' in request.POST):
					print("barabsfrbasdb")
					if(request.POST['survival_col']):
						survival_col_name = request.POST['survival_col']
				
				#if save_data in ["save_data"]:
	                	#	if request.user.is_authenticated:
	                	#		savedata_param = "true"
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				clinicalstr = "empty"
				if(request.user.is_authenticated and savedata_param == "true"):
					result2 = script_output_task_2.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,savedata_param,str(request.user))
					(div,script,plot1) = result2.get()			
				else:				
					clinicalstr = "empty"
					clinicaldf = ""
					survival_col_name = ""
					result2 = script_output_task_9.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf)
					(div,script,plot1,plot_div,ret_metadata) = result2.get()	
				
					#result2 = script_output_task_4.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2)
					#(div,script,plot1,plot_div,ret_metadata) = result2.get()
				plot2 = "test.png"
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			username = str(request.user)
	                			foobar = "user_uploaded_files/" + username
	                			if not(os.path.isdir(foobar)):
	                				os.mkdir(foobar)
	                			filename_1 = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S_%f')[:-3]
	                			filename_2 = foobar + "/" + username + filename_1 + "expr.txt"
	                			filename_3 = foobar + "/" + username + filename_1 + "prot.txt"
	                			with open(filename_2, "w") as text_file:
	                				text_file.write(exprstr)
	                				text_file.close()
	                			with open(filename_3, "w") as text_file_2:
	                				text_file_2.write(ppistr)
	                				text_file_2.close()
	                			#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
						#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
				#else:
					#cache.clear()
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)      	        
				#else:
				remove_loading_image.delay()
				cache.clear()				
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				return render(request, 'polls/clustering_6_part_3.html', {'list_of_files':list_of_files,'ret_metadata':ret_metadata,'result_99':ret_metadata,'list_of_files_2':list_of_files_2,'json_path':"test15.json",'heatmap_path':"test.png",'survival_graph':plot_div})

	elif('myfile' in request.FILES and 'parse_ndex_file' in request.POST and 'ndex_file_input' in request.FILES):
		#print(request.POST.get("ndex_name"))
		if(request.FILES['myfile'] and request.POST['parse_ndex_file'] and request.FILES['ndex_file_input']):
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				add_loading_image.delay()
				with open("polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				make_empty_figure.delay()
				ndex_file_str = request.FILES['ndex_file_input'].read().decode('utf-8')				
				result_ndex = read_ndex_file_4.delay(ndex_file_str)
				ppistr = result_ndex.get()
				exprstr = request.FILES['myfile'].read().decode('utf-8')
				savedata_param = "false"
				plot_div = ""
				#if save_data in ["save_data"]:
	                	#	if request.user.is_authenticated:
	                	#		savedata_param = "true"
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				clinicalstr = "empty"
				if(request.user.is_authenticated and savedata_param == "true"):
					result2 = script_output_task_2.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,savedata_param,str(request.user))
					(div,script,plot1) = result2.get()			
				else:				
					result2 = script_output_task_4.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2)
					(div,script,plot1,plot_div,ret_metadata) = result2.get()
				plot2 = "test.png"
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			username = str(request.user)
	                			foobar = "user_uploaded_files/" + username
	                			if not(os.path.isdir(foobar)):
	                				os.mkdir(foobar)
	                			filename_1 = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S_%f')[:-3]
	                			filename_2 = foobar + "/" + username + filename_1 + "expr.txt"
	                			filename_3 = foobar + "/" + username + filename_1 + "prot.txt"
	                			with open(filename_2, "w") as text_file:
	                				text_file.write(exprstr)
	                				text_file.close()
	                			with open(filename_3, "w") as text_file_2:
	                				text_file_2.write(ppistr)
	                				text_file_2.close()
	                			#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
						#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
				#else:
					#cache.clear()
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)      	        
				#else:
				remove_loading_image.delay()
				cache.clear()				
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				return render(request, 'polls/clustering_6_part_1.html', {'list_of_files':list_of_files,'ret_metadata':ret_metadata,'result_99':ret_metadata,'list_of_files_2':list_of_files_2,'json_path':"test15.json",'heatmap_path':"test.png",'survival_graph':plot_div})

	elif('myfile' in request.FILES and 'ndex_file' in request.POST and 'ndex_name' in request.POST and 'analyze_metadata' in request.POST):
		print(request.POST.get("ndex_name"))
		if(request.FILES['myfile'] and request.POST['ndex_file'] and request.POST['ndex_name'] and request.FILES['patientdata']):
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				add_loading_image.delay()
				with open("polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				make_empty_figure.delay()
				ndex_file_id = request.POST.get("ndex_name")
				if(ndex_file_id == "1"):
					fh2 = open("polls/data/APID.txt",encoding='utf-8')
				elif(ndex_file_id == "2"):
					fh2 = open("polls/data/STRING.txt",encoding='utf-8')
				elif(ndex_file_id == "3"):
					fh2 = open("polls/data/biogrid.human.entrez.tsv",encoding='utf-8')
				elif(ndex_file_id == "4"):
					fh2 = open("polls/data/HPRD.txt",encoding='utf-8')
				elif(ndex_file_id == "5"):
					fh2 = open("polls/data/ulitsky.txt",encoding='utf-8')
				elif(ndex_file_id == "6"):
					fh2 = open("polls/data/I2D.txt",encoding='utf-8')
				elif(ndex_file_id == "7"):
					fh2 = open("polls/data/biogrid_cel.txt",encoding='utf-8')
				elif(ndex_file_id == "8"):
					fh2 = open("polls/data/biogrid_ath.txt",encoding='utf-8')
				elif(ndex_file_id == "9"):
					fh2 = open("polls/data/biogrid_mus.txt",encoding='utf-8')
				elif(ndex_file_id == "10"):
					fh2 = open("polls/data/biogrid_dme.txt",encoding='utf-8')
				elif(ndex_file_id == "11"):
					fh2 = open("polls/data/biogrid_eco.txt",encoding='utf-8')
				elif(ndex_file_id == "12"):
					fh2 = open("polls/data/biogrid_sce.txt",encoding='utf-8')
				elif(ndex_file_id == "13"):
					fh2 = open("polls/data/biogrid_dan.txt",encoding='utf-8')
				ppistr = fh2.read()
				exprstr = request.FILES['myfile'].read().decode('utf-8')
				savedata_param = "false"
				plot_div = ""
				clinicalstr = request.FILES['patientdata'].read().decode('utf-8')
				clinical_stringio = StringIO(clinicalstr)
				clinicaldf = pd.read_csv(clinical_stringio)
				survival_col_name = ""
				if('survival_col' in request.POST):
					print("barabsfrbasdb")
					if(request.POST['survival_col']):
						survival_col_name = request.POST['survival_col']
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			savedata_param = "true"
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				#clinicalstr = "empty"
				#if(request.user.is_authenticated and savedata_param == "true"):
				#	result2 = script_output_task_2.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,savedata_param,str(request.user))
				#	(div,script,plot1) = result2.get()			
				#else:				
				#result2 = script_output_task_4.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2)
				#(div,script,plot1,plot_div,ret_metadata) = result2.get()
				result2 = script_output_task_9.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf)
				(div,script,plot1,plot_div,ret_metadata) = result2.get()	
				plot2 = "test.png"
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			username = str(request.user)
	                			foobar = "user_uploaded_files/" + username
	                			if not(os.path.isdir(foobar)):
	                				os.mkdir(foobar)
	                			filename_1 = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S_%f')[:-3]
	                			filename_2 = foobar + "/" + username + filename_1 + "expr.txt"
	                			filename_3 = foobar + "/" + username + filename_1 + "prot.txt"
	                			with open(filename_2, "w") as text_file:
	                				text_file.write(exprstr)
	                				text_file.close()
	                			with open(filename_3, "w") as text_file_2:
	                				text_file_2.write(ppistr)
	                				text_file_2.close()
	                			#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
						#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
				#else:
					#cache.clear()
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)      	        
				#else:
				remove_loading_image.delay()
				cache.clear()				
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				return render(request, 'polls/clustering_6_part_1.html', {'list_of_files':list_of_files,'ret_metadata':ret_metadata,'result_99':ret_metadata,'list_of_files_2':list_of_files_2,'json_path':"test15.json",'heatmap_path':"test.png",'survival_graph':plot_div})
	
	
	
	
	elif('own_file' in request.POST and request.user.is_authenticated):
		if(request.POST['own_file']):
			make_empty_figure.delay()
			with open("polls/static/output_console.txt", "w") as text_file:
   				text_file.write("Your request is being processed...")
			filename1 = request.POST.get("input_own_file")
			#fh1 = open(filename1)
			#exprstr = fh1.read()
			#filename2 = filename1.split("_expr.txt")[0] + "_prot.txt"
			#fh2 = open(filename2)
			path_json = filename1
			path_heatmap = filename1.split("_json.json")[0] + "_heatmap.png"
			#ppistr = fh2.read()
			list_of_files = GraphForm.list_user_data_2(username)	
			list_of_files_2 = GraphForm.list_user_data(username)      
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])	
				show_old_data(path_json,path_heatmap)             	
				plot2 = "test.png"
				cache.clear()			
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				return render(request, 'polls/clustering_6_part_1.html', {'form':"",'images':"",'div':"",'script':"",'plot2':plot2, 'list_of_files':list_of_files,'list_of_files_2':list_of_files_2})       		
	elif('analyze_metadata' in request.POST and 'myfile' in request.FILES and 'protfile' in request.FILES and 'patientdata' in request.FILES):
		if(request.FILES['myfile'] and request.FILES['protfile'] and request.FILES['patientdata']):
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				#return HttpResponseRedirect(reverse('clustering_6_part_2'))
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				add_loading_image.delay()
				with open("polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				make_empty_figure.delay()
				#fh2 = open(request.FILES['protfile'])
				ppistr = request.FILES['protfile'].read().decode('utf-8')
				#fh1 = open(request.FILES['myfile'])
				exprstr = request.FILES['myfile'].read().decode('utf-8')
				clinicalstr = request.FILES['patientdata'].read().decode('utf-8')
				clinical_stringio = StringIO(clinicalstr)
				clinicaldf = pd.read_csv(clinical_stringio)
				print(clinicaldf)
				survival_col_name = ""
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			savedata_param = "true"
				if('survival_col' in request.POST):
					print("barabsfrbasdb")
					if(request.POST['survival_col']):
						survival_col_name = request.POST['survival_col']
				result4 = check_input_files.delay(ppistr,exprstr)
				errstr = result4.get()
				if(errstr != ""):
					request.session['errors'] = errstr
					return render(request,'polls/errorpage.html',{'errors':errstr})
				
				#return HttpResponseRedirect('clustering_6_part_2_2',exprstr_par=exprstr,ppistr_par=ppistr)
				#return HttpResponseRedirect('polls/clustering_6_part_2_2.html',exprstr,ppistr)
				#return HttpResponseRedirect('polls/clustering_6_part_2_2.html',args={'exprstr_par':exprstr,'ppistr_par':ppistr})
				#return clustering_6_part_2_2(request,exprstr,ppistr)
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				#result2 = script_output_task.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1)
				#result2 = script_output_task_4.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2)
				#result2 = script_output_task_8.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name)
				#(div,script,plot1,plot_div,ret_metadata) = result2.get()
				#result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				#(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				#print("ikkrngasndgnn")
				#result2 = script_output_task_4.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2)
				#result2 = script_output_task_9.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf)
				session_id = request.session._get_or_create_session_key()
				output_plot_path = "polls/static/output_plotly_" + session_id + ".html"
				json_path = "polls/static/test15_" + session_id + ".json"
				result2 = script_output_task_10.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf,session_id)
				(div,script,plot1,plot_div,ret_metadata,path99,path_metadata,output_plot_path,json_path) = result2.get()
				#metd = list_metadata_3.apply_async(countdown=0)
				path_metadata = "polls/static/metadata_" + session_id + ".txt"
				metd = list_metadata_4.apply_async(args=[path_metadata],countdown=0)
				(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get()
				print("iasdfasdfsf")
				print(ret_metadata1) 	
				request.session['ppistr'] = ppistr
				#return HttpResponseRedirect(reverse('polls/clustering_6_part_2_2/'),args=(exprstr,ppistr))
				#return clustering_6_part_2_2(request,exprstr,ppistr)
				#(div,script,plot1) = result2.get()			
				plot2 = "test.png"
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			savedata_param = "true"
	                			username = str(request.user)
	                			GraphForm.save_user_data_2(request.FILES['myfile'],request.FILES['protfile'],request.FILES['patientdata'],username)
						#GraphForm.save_results(username)
				#else:
					#cache.clear()
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)              
				#else:
				remove_loading_image.delay()
				cache.clear()				
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				#session_id = request.session._get_or_create_session_key()
				print(session_id)
				request.session['done'] = "true"
				#request.session.modified = True
				print(request.session['done'])
				#request.session['plot_div'] = plot_div
				#list_metadata.apply_async(countdown=0)
				#return clustering_6_part_2_2(request)
				#request.session['plot_div'] = plot_div
				#request.session['script'] = script
				#request.session['plot2'] = plot2
				#request.session['list_of_files'] = list_of_files
				#request.session['ret_dat'] = ret_metadata
				#request.session['ret_metadata1'] = ret_metadata1
				#request.session['ret_metadata2'] = ret_metadata2
				#request.session['ret_metadata3'] = ret_metadata3
				#request.session['list_of_files_2'] = list_of_files_2
				#return render(request,'polls/clustering_6_part_2_2.html')
				#return clustering_6_part_2_2(request,plot_div,script,plot2,list_of_files,ret_metadata,ret_metadata1,ret_metadata2,ret_metadata3,list_of_files_2)
				with open("polls/static/plotly_output_2.html", "w") as text_file_3:
   					text_file_3.write(plot_div)
				return render(request, 'polls/clustering_6_part_3.html', {'form':"",'images':"",'plot_div':plot_div,'script':script,'plot2':plot2,'plot4':path99,'output_plot_path':output_plot_path,'json_path':json_path, 'list_of_files':list_of_files,'ret_dat':ret_metadata,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'list_of_files_2':list_of_files_2})
				#return render(request, 'polls/clustering_6_part_1.html')
	elif('myfile' in request.FILES and 'protfile' in request.FILES):
		if(request.FILES['myfile'] and request.FILES['protfile']):
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				add_loading_image.delay()
				with open("polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				make_empty_figure.delay()
				#fh2 = open(request.FILES['protfile'])
				ppistr = request.FILES['protfile'].read().decode('utf-8')
				#fh1 = open(request.FILES['myfile'])
				exprstr = request.FILES['myfile'].read().decode('utf-8')
				savedata_param = "false"
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			savedata_param = "true"
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				#result2 = script_output_task.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1)
				clinicalstr = "empty"
				if(request.user.is_authenticated and savedata_param == "true"):
					result2 = script_output_task_2.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,savedata_param,str(request.user))
					(div,script,plot1) = result2.get()			
				else:				
					clinicalstr = "empty"
					clinicaldf = ""
					survival_col_name = ""
					result2 = script_output_task_9.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf)
					(div,script,plot1,plot_div,ret_metadata) = result2.get()	
				
					#result2 = script_output_task_4.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2)
					#(div,script,plot1,plot_div,ret_metadata) = result2.get()
				plot2 = "test.png"
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			savedata_param = "true"
	                			username = str(request.user)
	                			foobar = "user_uploaded_files/" + username
	                			if not(os.path.isdir(foobar)):
	                				os.mkdir(foobar)
	                			filename_1 = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S_%f')[:-3]
	                			filename_2 = foobar + "/" + username + filename_1 + "expr.txt"
	                			filename_3 = foobar + "/" + username + filename_1 + "prot.txt"
	                			with open(filename_2, "w") as text_file:
	                				text_file.write(exprstr)
	                				text_file.close()
	                			with open(filename_3, "w") as text_file_2:
	                				text_file_2.write(ppistr)
	                				text_file_2.close()
	                		
	                			#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
				#else:
					#cache.clear()
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)      	        
				#else:
				remove_loading_image.delay()
				cache.clear()				
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				return render(request, 'polls/clustering_6_part_1.html', {'form':"",'images':"",'div':div,'script':script,'plot2':plot2, 'list_of_files':list_of_files,'list_of_files_2':list_of_files_2})
	elif('predef_file' in request.POST):
		if('L_g_min' in request.POST and 'L_g_max' in request.POST):
			if(request.POST['predef_file']):
				make_empty_figure.delay()
				with open("polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				add_loading_image.delay()
				cancer_type = request.POST.get("cancer_type")
				print("Cancer Type:" + str(cancer_type))
				fh2 = open("polls/data/biogrid.human.entrez.tsv")
				ppistr = fh2.read()
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				is_lungc = "false"					
				if(cancer_type == "1"):
					print("babababababa")
					fh1 = open("polls/data/lung_cancer_expr.csv")
					exprstr = fh1.read()
					clinicaldf = pd.read_csv("polls/data/lung_cancer_clinical.csv")
					fh4 = open("polls/data/lung_cancer_clinical.csv")
					#fh4 = open("polls/data/breast_cancer_clinical.csv")
					is_lungc = "true"
					clinicalstr = fh4.read()
					fh4.flush()
					fh4.close()
					survival_col_name = "disease free survival in months:ch1"
				else:
					fh1 = open("polls/data/breast_cancer_expr.csv")
					exprstr = fh1.read()
					clinicaldf = pd.read_csv("polls/data/breast_cancer_clinical.csv")
					fh4 = open("polls/data/breast_cancer_clinical.csv")
					clinicalstr = fh4.read()
					fh4.flush()
					fh4.close()
					survival_col_name = "mfs (yr):ch1"
				print(clinicaldf)
				if('survival_col' in request.POST):
					print("barabsfrbasdb")
					if(request.POST['survival_col']):
						survival_col_name = request.POST['survival_col']
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				print("ikkrngasndgnn")
				#result2 = script_output_task_4.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2)
				result2 = script_output_task_9.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf)
				(div,script,plot1,plot_div,ret_metadata) = result2.get()
				print(ret_metadata)
				metd = list_metadata_3.apply_async(countdown=0)
				(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get()
				plot2 = "test.png"
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)        				
				#else:
				remove_loading_image.delay()
				with open("polls/static/metadata_99.html", "w") as text_file_4:
   					text_file_4.write(str(ret_metadata[0]))
   					text_file_4.write(str(ret_metadata[1]))
   					text_file_4.write(str(ret_metadata[2]))
				cache.clear()
				with open("polls/static/plotly_output_2.html", "w") as text_file_3:
   					text_file_3.write(plot_div)
				#write_metadata_to_file.apply_async((ret_metadata), countdown=1)
				#list_metadata.apply_async(countdown=0)
				string_99 = metadata_to_string.apply_async(ret_metadata, countdown=0)
				result_99 = string_99.get()
				return render(request, 'polls/clustering_6_part_1.html', {'plot2':plot2, 'list_of_files':list_of_files,'ret_metadata':ret_metadata,'result_99':ret_metadata,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'list_of_files_2':list_of_files_2})
	elif('kegg_enrichment' in request.POST):
		group_for_enr = "both"
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		if('group_for_enr' in request.POST):
			group_for_enr = request.POST.get('group_for_enr')
			#print(pval_enr)
		result1 = run_enrichment_2.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_kegg")
		result2 = run_enrichment_2.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_kegg")
		result3 = run_enrichment_2.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_kegg")
		enr_results = result1.get()
		enr_results_2 = result2.get()
		enr_results_3 = result3.get()
		print("enr")
		result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
		result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
		result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
		enrichment_dict = result4.get()
		enrichment_dict_2 = result5.get()
		enrichment_dict_3 = result6.get()
		return render(request,'polls/clustering_6.html',{'list_of_files':list_of_files,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3})
	elif('enrichment_type' in request.POST):
		enr_type = request.POST.get("enrichment_type")
		group_for_enr = "both"
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		if('group_for_enr' in request.POST):
			group_for_enr = request.POST.get('group_for_enr')
			#print(pval_enr)
		enrichment_dict = {}
		enrichment_dict_2 = {}
		enrichment_dict_3 = {}
		enrichment_dict_4 = {}
		enrichment_dict_5 = {}
		if(enr_type == "kegg_enrichment"):
			result1 = run_enrichment_2.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_kegg")
			result2 = run_enrichment_2.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_kegg")
			result3 = run_enrichment_2.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_kegg")
			enr_results = result1.get()
			enr_results_2 = result2.get()
			enr_results_3 = result3.get()
			print("enr")
			result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			result7 = read_kegg_enrichment_2.delay("polls/data/test2/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt","polls/data/test3/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			enrichment_dict = result4.get()
			enrichment_dict_2 = result5.get()
			enrichment_dict_3 = result6.get()
			(enrichment_dict_4,enrichment_dict_5) = result7.get()
		elif(enr_type == "go_enrichment"):	
			result1 = run_go_enrichment_2.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_go")
			result2 = run_go_enrichment_2.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_go")
			result3 = run_go_enrichment_2.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_go")
			#result1 = run_go_enrichment.delay("genelist.txt",pval_enr)
			enr_results = result1.get()
			enr_results_2 = result2.get()
			enr_results_3 = result3.get()
			print("enr")
			result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			result7 = read_kegg_enrichment_2.delay("polls/data/test2/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt","polls/data/test3/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			enrichment_dict = result4.get()
			enrichment_dict_2 = result5.get()
			enrichment_dict_3 = result6.get()	
			(enrichment_dict_4,enrichment_dict_5) = result7.get()
		elif(enr_type == "go_molecular"):
			result1 = run_go_enrichment_2.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_go")
			result2 = run_go_enrichment_2.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_go")
			result3 = run_go_enrichment_2.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_go")
			enr_results = result1.get()
			enr_results_2 = result2.get()
			enr_results_3 = result3.get()
			print("enr")
			result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result7 = read_kegg_enrichment_2.delay("polls/data/test2/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt","polls/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			enrichment_dict = result4.get()
			enrichment_dict_2 = result5.get()
			enrichment_dict_3 = result6.get()
			(enrichment_dict_4,enrichment_dict_5) = result7.get()
		elif(enr_type == "reactome_enrichment"):
			#result1 = run_reac_enrichment.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_reactome")
			result2 = run_reac_enrichment.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_reactome")
			#result3 = run_reac_enrichment.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_reactome")
			#enr_results = result1.get()
			enr_results_2 = result2.get()
			#enr_results_3 = result3.get()
			print("enr")
			#result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_reactome/Reactome_2016.test_name.enrichr.reports.txt",pval_enr)
			#result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			#enrichment_dict = result4.get()
			enrichment_dict = {}
			enrichment_dict_2 = result5.get()
			#enrichment_dict_3 = result6.get()
			enrichment_dict_3 = {}		
		return render(request,'polls/clustering_6.html',{'list_of_files':list_of_files,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3,'enrichment_dict_4':enrichment_dict_4,'enrichment_dict_5':enrichment_dict_5})

	elif('go_enrichment' in request.POST):
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		group_for_enr = "both"
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		if('group_for_enr' in request.POST):
			group_for_enr = request.POST.get('group_for_enr')
			#print(pval_enr)
		result1 = run_go_enrichment_2.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_go")
		result2 = run_go_enrichment_2.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_go")
		result3 = run_go_enrichment_2.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_go")
		#result1 = run_go_enrichment.delay("genelist.txt",pval_enr)
		enr_results = result1.get()
		enr_results_2 = result2.get()
		enr_results_3 = result3.get()
		print("enr")
		result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
		result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
		result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
		enrichment_dict = result4.get()
		enrichment_dict_2 = result5.get()
		enrichment_dict_3 = result6.get()
		return render(request,'polls/clustering_6.html',{'list_of_files':list_of_files,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3})
	elif('go_molecular' in request.POST):
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		#result1 = run_go_enrichment.delay("genelist.txt",pval_enr)
		group_for_enr = "both"
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		if('group_for_enr' in request.POST):
			group_for_enr = request.POST.get('group_for_enr')
			#print(pval_enr)
		result1 = run_go_enrichment_2.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_go")
		result2 = run_go_enrichment_2.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_go")
		result3 = run_go_enrichment_2.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_go")
		enr_results = result1.get()
		enr_results_2 = result2.get()
		enr_results_3 = result3.get()
		print("enr")
		result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
		result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
		result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
		enrichment_dict = result4.get()
		enrichment_dict_2 = result5.get()
		enrichment_dict_3 = result6.get()
		return render(request,'polls/clustering_6.html',{'list_of_files':list_of_files,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3})
	elif('reac_enrichment' in request.POST):
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		#result1 = run_go_enrichment.delay("genelist.txt",pval_enr)
		group_for_enr = "both"
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		if('group_for_enr' in request.POST):
			group_for_enr = request.POST.get('group_for_enr')
			#print(pval_enr)
		#result1 = run_reac_enrichment.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_reactome")
		result2 = run_reac_enrichment.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_reactome")
		#result3 = run_reac_enrichment.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_reactome")
		#enr_results = result1.get()
		enr_results_2 = result2.get()
		#enr_results_3 = result3.get()
		print("enr")
		#result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
		result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_reactome/Reactome_2016.test_name.enrichr.reports.txt",pval_enr)
		#result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
		#enrichment_dict = result4.get()
		enrichment_dict = {}
		enrichment_dict_2 = result5.get()
		#enrichment_dict_3 = result6.get()
		enrichment_dict_3 = {}
		return render(request,'polls/clustering_6.html',{'list_of_files':list_of_files,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3})

	#if not(('myfile' in request.FILES and 'protfile' in request.FILES) or ('predef_file' in request.POST)):
	else:		
		#remove_loading_image.delay()
		ret_metadata = ""
		if not (request.user.is_authenticated):
			cache.clear()
		else:
			username = str(request.user)
			list_of_files = GraphForm.list_user_data_2(username)
			list_of_files_2 = GraphForm.list_user_data(username)
			cache.clear()
			#metd = list_metadata.apply_async(countdown=0)			
			#metd = list_metadata_3.apply_async(countdown=0)
			session_id = request.session._get_or_create_session_key()
			path_metadata = "polls/static/metadata_" + session_id + ".txt"
			metd = list_metadata_4.apply_async(args=[path_metadata],countdown=0)
			(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get() 
			print(ret_metadata1) 
			print("iubaerb")
			metadata_dict = [ret_metadata1,ret_metadata2,ret_metadata3]
			result2 = read_kegg_enrichment.delay("polls/data/test/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			#result2 = read_kegg_enrichment.delay("data/test/enrichr_kegg")	
			enrichment_dict = result2.get()
			#(ret_metadata,is_lungc) = metd.get()  
			#return render(request,'polls/clustering_6_part_1.html',{'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'ret_metadata':ret_metadata,'is_lungc':is_lungc})

		#make_empty_figure.apply_async(countdown=2)
		#empty_log_file.apply_async(countdown=2)
		#print(request.session['done'])
		#print(request.session['done'])
		#if(request.session['done'] == "true"):
		#	print("done")
		if('done' in request.session):
			print("clustering 6 part 3")
			print(request.session['done'])
			if(request.session['done'] == "true"): 
				session_id = request.session._get_or_create_session_key()
				path99 = "/home/quirin/testproject/polls/static/test_" + session_id + ".png"
				path_metadata = "polls/static/metadata_" + session_id + ".txt"
				output_plot_path = "polls/static/output_plotly_" + session_id + ".html"
				return render(request,'polls/clustering_6_part_3.html',{'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'ret_metadata':ret_metadata,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'path4':path99,'output_plot_path':output_plot_path,'metadata_dict':metadata_dict,'enrichment_dict':enrichment_dict})

		return render(request,'polls/clustering_6_part_1.html',{'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'ret_metadata':ret_metadata,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'metadata_dict':metadata_dict,'enrichment_dict':enrichment_dict})


def clustering_6_part_2(request):
	ret_metadata1 = {}
	ret_metadata2 = {}
	ret_metadata3 = {}
	metadata_dict = []
	enrichment_dict = []
	pval_enr = 0.5
	list_of_files = ""
	list_of_files_2 = ""
	save_data = request.POST.get("save_data", None)
	nbr_iter = request.POST.get("nbr_iter",45)
	nbr_ants = request.POST.get("nbr_ants",30)
	evap = request.POST.get("evap",0.3)
	epsilon = request.POST.get("stopcr",0.02) 
	hi_sig = request.POST.get("hisig",1)
	pher_sig = request.POST.get("pher",1)
	savedata_param = "false"
	#fig = plt.figure(figsize=(10,8))
	#plt.savefig("polls/static/progress.png")
	#plt.close(fig)
	if request.user.is_authenticated:
		username = str(request.user)
		list_of_files = GraphForm.list_user_data_2(username)
	if('analyze_metadata' in request.POST and 'myfile' in request.FILES and 'protfile' in request.FILES and 'patientdata' in request.FILES):
		if(request.FILES['myfile'] and request.FILES['protfile'] and request.FILES['patientdata']):
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				print("basfbasdbasdbasdb")
				#return HttpResponseRedirect(reverse('clustering_6_step_2'))
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				add_loading_image.delay()
				with open("polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				make_empty_figure.delay()
				#fh2 = open(request.FILES['protfile'])
				ppistr = request.FILES['protfile'].read().decode('utf-8')
				#fh1 = open(request.FILES['myfile'])
				exprstr = request.FILES['myfile'].read().decode('utf-8')
				clinicalstr = request.FILES['patientdata'].read().decode('utf-8')
				clinical_stringio = StringIO(clinicalstr)
				clinicaldf = pd.read_csv(clinical_stringio)
				print(clinicaldf)
				survival_col_name = ""
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			savedata_param = "true"
				if('survival_col' in request.POST):
					print("barabsfrbasdb")
					if(request.POST['survival_col']):
						survival_col_name = request.POST['survival_col']
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				#result2 = script_output_task.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1)
				#result2 = script_output_task_4.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2)
				#result2 = script_output_task_8.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name)
				#(div,script,plot1,plot_div,ret_metadata) = result2.get()
				#result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				#(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				#print("ikkrngasndgnn")
				#result2 = script_output_task_4.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2)
				result2 = script_output_task_9.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf)
				(div,script,plot1,plot_div,ret_metadata) = result2.get()
				metd = list_metadata_3.apply_async(countdown=0)
				(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get()
				print("iasdfasdfsf")
				print(ret_metadata1) 	
				#(div,script,plot1) = result2.get()			
				plot2 = "test.png"
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			savedata_param = "true"
	                			username = str(request.user)
	                			GraphForm.save_user_data_2(request.FILES['myfile'],request.FILES['protfile'],request.FILES['patientdata'],username)
						#GraphForm.save_results(username)
				#else:
					#cache.clear()
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)              
				#else:
				remove_loading_image.delay()
				cache.clear()				
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				#list_metadata.apply_async(countdown=0)
				with open("polls/static/plotly_output_2.html", "w") as text_file_3:
   					text_file_3.write(plot_div)
				return render(request, 'polls/clustering_6_part_3.html', {'form':"",'images':"",'plot_div':plot_div,'script':script,'plot2':plot2, 'list_of_files':list_of_files,'ret_dat':ret_metadata,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'list_of_files_2':list_of_files_2})
	elif('kegg_enrichment' in request.POST):
		group_for_enr = "both"
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		if('group_for_enr' in request.POST):
			group_for_enr = request.POST.get('group_for_enr')
			#print(pval_enr)
		result1 = run_enrichment_2.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_kegg")
		result2 = run_enrichment_2.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_kegg")
		result3 = run_enrichment_2.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_kegg")
		enr_results = result1.get()
		enr_results_2 = result2.get()
		enr_results_3 = result3.get()
		print("enr")
		result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
		result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
		result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
		enrichment_dict = result4.get()
		enrichment_dict_2 = result5.get()
		enrichment_dict_3 = result6.get()
		return render(request,'polls/clustering_6.html',{'list_of_files':list_of_files,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3})
	elif('enrichment_type' in request.POST):
		enr_type = request.POST.get("enrichment_type")
		group_for_enr = "both"
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		if('group_for_enr' in request.POST):
			group_for_enr = request.POST.get('group_for_enr')
			#print(pval_enr)
		enrichment_dict = {}
		enrichment_dict_2 = {}
		enrichment_dict_3 = {}
		enrichment_dict_4 = {}
		enrichment_dict_5 = {}
		if(enr_type == "kegg_enrichment"):
			result1 = run_enrichment_2.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_kegg")
			result2 = run_enrichment_2.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_kegg")
			result3 = run_enrichment_2.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_kegg")
			enr_results = result1.get()
			enr_results_2 = result2.get()
			enr_results_3 = result3.get()
			print("enr")
			result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			result7 = read_kegg_enrichment_2.delay("polls/data/test2/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt","polls/data/test3/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			enrichment_dict = result4.get()
			enrichment_dict_2 = result5.get()
			enrichment_dict_3 = result6.get()
			(enrichment_dict_4,enrichment_dict_5) = result7.get()
		elif(enr_type == "go_enrichment"):	
			result1 = run_go_enrichment_2.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_go")
			result2 = run_go_enrichment_2.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_go")
			result3 = run_go_enrichment_2.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_go")
			#result1 = run_go_enrichment.delay("genelist.txt",pval_enr)
			enr_results = result1.get()
			enr_results_2 = result2.get()
			enr_results_3 = result3.get()
			print("enr")
			result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			result7 = read_kegg_enrichment_2.delay("polls/data/test2/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt","polls/data/test3/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			enrichment_dict = result4.get()
			enrichment_dict_2 = result5.get()
			enrichment_dict_3 = result6.get()	
			(enrichment_dict_4,enrichment_dict_5) = result7.get()
		elif(enr_type == "go_molecular"):
			result1 = run_go_enrichment_2.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_go")
			result2 = run_go_enrichment_2.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_go")
			result3 = run_go_enrichment_2.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_go")
			enr_results = result1.get()
			enr_results_2 = result2.get()
			enr_results_3 = result3.get()
			print("enr")
			result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result7 = read_kegg_enrichment_2.delay("polls/data/test2/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt","polls/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			enrichment_dict = result4.get()
			enrichment_dict_2 = result5.get()
			enrichment_dict_3 = result6.get()
			(enrichment_dict_4,enrichment_dict_5) = result7.get()
		elif(enr_type == "reactome_enrichment"):
			#result1 = run_reac_enrichment.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_reactome")
			result2 = run_reac_enrichment.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_reactome")
			#result3 = run_reac_enrichment.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_reactome")
			#enr_results = result1.get()
			enr_results_2 = result2.get()
			#enr_results_3 = result3.get()
			print("enr")
			#result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_reactome/Reactome_2016.test_name.enrichr.reports.txt",pval_enr)
			#result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			#enrichment_dict = result4.get()
			enrichment_dict = {}
			enrichment_dict_2 = result5.get()
			#enrichment_dict_3 = result6.get()
			enrichment_dict_3 = {}		
		return render(request,'polls/clustering_6.html',{'list_of_files':list_of_files,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3,'enrichment_dict_4':enrichment_dict_4,'enrichment_dict_5':enrichment_dict_5})

	elif('go_enrichment' in request.POST):
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		group_for_enr = "both"
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		if('group_for_enr' in request.POST):
			group_for_enr = request.POST.get('group_for_enr')
			#print(pval_enr)
		result1 = run_go_enrichment_2.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_go")
		result2 = run_go_enrichment_2.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_go")
		result3 = run_go_enrichment_2.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_go")
		#result1 = run_go_enrichment.delay("genelist.txt",pval_enr)
		enr_results = result1.get()
		enr_results_2 = result2.get()
		enr_results_3 = result3.get()
		print("enr")
		result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
		result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
		result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
		enrichment_dict = result4.get()
		enrichment_dict_2 = result5.get()
		enrichment_dict_3 = result6.get()
		return render(request,'polls/clustering_6.html',{'list_of_files':list_of_files,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3})
	elif('go_molecular' in request.POST):
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		#result1 = run_go_enrichment.delay("genelist.txt",pval_enr)
		group_for_enr = "both"
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		if('group_for_enr' in request.POST):
			group_for_enr = request.POST.get('group_for_enr')
			#print(pval_enr)
		result1 = run_go_enrichment_2.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_go")
		result2 = run_go_enrichment_2.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_go")
		result3 = run_go_enrichment_2.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_go")
		enr_results = result1.get()
		enr_results_2 = result2.get()
		enr_results_3 = result3.get()
		print("enr")
		result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
		result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
		result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
		enrichment_dict = result4.get()
		enrichment_dict_2 = result5.get()
		enrichment_dict_3 = result6.get()
		return render(request,'polls/clustering_6.html',{'list_of_files':list_of_files,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3})
	elif('reac_enrichment' in request.POST):
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		#result1 = run_go_enrichment.delay("genelist.txt",pval_enr)
		group_for_enr = "both"
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		if('group_for_enr' in request.POST):
			group_for_enr = request.POST.get('group_for_enr')
			#print(pval_enr)
		#result1 = run_reac_enrichment.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_reactome")
		result2 = run_reac_enrichment.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_reactome")
		#result3 = run_reac_enrichment.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_reactome")
		#enr_results = result1.get()
		enr_results_2 = result2.get()
		#enr_results_3 = result3.get()
		print("enr")
		#result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
		result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_reactome/Reactome_2016.test_name.enrichr.reports.txt",pval_enr)
		#result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
		#enrichment_dict = result4.get()
		enrichment_dict = {}
		enrichment_dict_2 = result5.get()
		#enrichment_dict_3 = result6.get()
		enrichment_dict_3 = {}
		return render(request,'polls/clustering_6.html',{'list_of_files':list_of_files,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3})

	#if not(('myfile' in request.FILES and 'protfile' in request.FILES) or ('predef_file' in request.POST)):
	else:		
		#remove_loading_image.delay()
		ret_metadata = ""
		if not (request.user.is_authenticated):
			cache.clear()
		else:
			username = str(request.user)
			list_of_files = GraphForm.list_user_data_2(username)
			list_of_files_2 = GraphForm.list_user_data(username)
			cache.clear()
			#metd = list_metadata.apply_async(countdown=0)
			metd = list_metadata_3.apply_async(countdown=0)
			(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get() 
			print(ret_metadata1) 
			print("iubaerb")
			metadata_dict = [ret_metadata1,ret_metadata2,ret_metadata3]
			result2 = read_kegg_enrichment.delay("polls/data/test/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			#result2 = read_kegg_enrichment.delay("data/test/enrichr_kegg")	
			enrichment_dict = result2.get()
			#(ret_metadata,is_lungc) = metd.get()  
			#return render(request,'polls/clustering_6.html',{'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'ret_metadata':ret_metadata,'is_lungc':is_lungc})

		#make_empty_figure.apply_async(countdown=2)
		#empty_log_file.apply_async(countdown=2) 				
		return render(request,'polls/clustering_6_part_2.html',{'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'ret_metadata':ret_metadata,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'metadata_dict':metadata_dict,'enrichment_dict':enrichment_dict})





def clustering_6_part_2_2(request):
		#print("badfbasdbasdbasdb")
		#return HttpResponseRedirect('polls/clustering_6_part_3.html')
		#print(exprstr_par)
		print(request.session['ppistr'])
		if('done' in request.session):
			print(request.session['done'])
		#if(request.session.get('done') == "True"):
		#	print("done")
		return HttpResponseRedirect('polls/clustering_6_part_3.html')
		#return render(request,'polls/clustering_6_part_3.html')
		#return(clustering_6_part_3_2(request))




def errorpage(request):
		errors = ""
		if('errors' in request.session):
			errors = request.session['errors']
		#print("badfbasdbasdbasdb")
		#return HttpResponseRedirect('polls/clustering_6_part_3.html')
		#print(exprstr_par)
		#print(request.session['ppistr'])
		#if('done' in request.session):
		#	print(request.session['done'])
		#if(request.session.get('done') == "True"):
		#	print("done")
		#return HttpResponseRedirect('polls/clustering_6_part_3.html')
		return render(request,'polls/errorpage.html',{'errors':errors})
		#return(clustering_6_part_3_2(request))



def clustering_6_part_3(request):
	print("basdbasdbasdbasdbs")
	ret_metadata1 = {}
	ret_metadata2 = {}
	ret_metadata3 = {}
	metadata_dict = []
	enrichment_dict = []
	pval_enr = 0.5
	list_of_files = ""
	list_of_files_2 = ""
	save_data = request.POST.get("save_data", None)
	nbr_iter = request.POST.get("nbr_iter",45)
	nbr_ants = request.POST.get("nbr_ants",30)
	evap = request.POST.get("evap",0.3)
	epsilon = request.POST.get("stopcr",0.02) 
	hi_sig = request.POST.get("hisig",1)
	pher_sig = request.POST.get("pher",1)
	savedata_param = "false"
	#fig = plt.figure(figsize=(10,8))
	#plt.savefig("polls/static/progress.png")
	#plt.close(fig)
	if request.user.is_authenticated:
		username = str(request.user)
		list_of_files = GraphForm.list_user_data_2(username)
	if('analyze_metadata' in request.POST and 'myfile' in request.FILES and 'protfile' in request.FILES and 'patientdata' in request.FILES):
		if(request.FILES['myfile'] and request.FILES['protfile'] and request.FILES['patientdata']):
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				#return HttpResponseRedirect(reverse('clustering_6_step_2'))
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				add_loading_image.delay()
				with open("polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				make_empty_figure.delay()
				#fh2 = open(request.FILES['protfile'])
				ppistr = request.FILES['protfile'].read().decode('utf-8')
				#fh1 = open(request.FILES['myfile'])
				exprstr = request.FILES['myfile'].read().decode('utf-8')
				clinicalstr = request.FILES['patientdata'].read().decode('utf-8')
				clinical_stringio = StringIO(clinicalstr)
				clinicaldf = pd.read_csv(clinical_stringio)
				print(clinicaldf)
				survival_col_name = ""
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			savedata_param = "true"
				if('survival_col' in request.POST):
					print("barabsfrbasdb")
					if(request.POST['survival_col']):
						survival_col_name = request.POST['survival_col']
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				#result2 = script_output_task.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1)
				#result2 = script_output_task_4.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2)
				#result2 = script_output_task_8.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name)
				#(div,script,plot1,plot_div,ret_metadata) = result2.get()
				#result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				#(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				#print("ikkrngasndgnn")
				#result2 = script_output_task_4.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2)
				result2 = script_output_task_9.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf)
				(div,script,plot1,plot_div,ret_metadata) = result2.get()
				metd = list_metadata_3.apply_async(countdown=0)
				(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get()
				print("iasdfasdfsf")
				print(ret_metadata1) 	
				#(div,script,plot1) = result2.get()			
				plot2 = "test.png"
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			savedata_param = "true"
	                			username = str(request.user)
	                			GraphForm.save_user_data_2(request.FILES['myfile'],request.FILES['protfile'],request.FILES['patientdata'],username)
						#GraphForm.save_results(username)
				#else:
					#cache.clear()
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)              
				#else:
				remove_loading_image.delay()
				cache.clear()				
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				#list_metadata.apply_async(countdown=0)
				with open("polls/static/plotly_output_2.html", "w") as text_file_3:
   					text_file_3.write(plot_div)
				return render(request, 'polls/clustering_6_part_3.html', {'form':"",'images':"",'plot_div':plot_div,'script':script,'plot2':plot2, 'list_of_files':list_of_files,'ret_dat':ret_metadata,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'list_of_files_2':list_of_files_2})
	elif('kegg_enrichment' in request.POST):
		group_for_enr = "both"
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		if('group_for_enr' in request.POST):
			group_for_enr = request.POST.get('group_for_enr')
			#print(pval_enr)
		result1 = run_enrichment_2.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_kegg")
		result2 = run_enrichment_2.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_kegg")
		result3 = run_enrichment_2.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_kegg")
		enr_results = result1.get()
		enr_results_2 = result2.get()
		enr_results_3 = result3.get()
		print("enr")
		result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
		result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
		result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
		enrichment_dict = result4.get()
		enrichment_dict_2 = result5.get()
		enrichment_dict_3 = result6.get()
		return render(request,'polls/clustering_6_part_3.html',{'list_of_files':list_of_files,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3})
	elif('enrichment_type' in request.POST):
		enr_type = request.POST.get("enrichment_type")
		group_for_enr = "both"
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		if('group_for_enr' in request.POST):
			group_for_enr = request.POST.get('group_for_enr')
			#print(pval_enr)
		enrichment_dict = {}
		enrichment_dict_2 = {}
		enrichment_dict_3 = {}
		enrichment_dict_4 = {}
		enrichment_dict_5 = {}
		if(enr_type == "kegg_enrichment"):
			result1 = run_enrichment_2.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_kegg")
			result2 = run_enrichment_2.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_kegg")
			result3 = run_enrichment_2.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_kegg")
			enr_results = result1.get()
			enr_results_2 = result2.get()
			enr_results_3 = result3.get()
			print("enr")
			result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			result7 = read_kegg_enrichment_2.delay("polls/data/test2/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt","polls/data/test3/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			enrichment_dict = result4.get()
			enrichment_dict_2 = result5.get()
			enrichment_dict_3 = result6.get()
			(enrichment_dict_4,enrichment_dict_5) = result7.get()
		elif(enr_type == "go_enrichment"):	
			result1 = run_go_enrichment_2.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_go")
			result2 = run_go_enrichment_2.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_go")
			result3 = run_go_enrichment_2.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_go")
			#result1 = run_go_enrichment.delay("genelist.txt",pval_enr)
			enr_results = result1.get()
			enr_results_2 = result2.get()
			enr_results_3 = result3.get()
			print("enr")
			result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			result7 = read_kegg_enrichment_2.delay("polls/data/test2/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt","polls/data/test3/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			enrichment_dict = result4.get()
			enrichment_dict_2 = result5.get()
			enrichment_dict_3 = result6.get()	
			(enrichment_dict_4,enrichment_dict_5) = result7.get()
		elif(enr_type == "go_molecular"):
			result1 = run_go_enrichment_2.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_go")
			result2 = run_go_enrichment_2.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_go")
			result3 = run_go_enrichment_2.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_go")
			enr_results = result1.get()
			enr_results_2 = result2.get()
			enr_results_3 = result3.get()
			print("enr")
			result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result7 = read_kegg_enrichment_2.delay("polls/data/test2/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt","polls/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			enrichment_dict = result4.get()
			enrichment_dict_2 = result5.get()
			enrichment_dict_3 = result6.get()
			(enrichment_dict_4,enrichment_dict_5) = result7.get()
		elif(enr_type == "reactome_enrichment"):
			#result1 = run_reac_enrichment.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_reactome")
			result2 = run_reac_enrichment.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_reactome")
			#result3 = run_reac_enrichment.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_reactome")
			#enr_results = result1.get()
			enr_results_2 = result2.get()
			#enr_results_3 = result3.get()
			print("enr")
			#result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_reactome/Reactome_2016.test_name.enrichr.reports.txt",pval_enr)
			#result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			#enrichment_dict = result4.get()
			enrichment_dict = {}
			enrichment_dict_2 = result5.get()
			#enrichment_dict_3 = result6.get()
			enrichment_dict_3 = {}		
		return render(request,'polls/clustering_6_part_3.html',{'list_of_files':list_of_files,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3,'enrichment_dict_4':enrichment_dict_4,'enrichment_dict_5':enrichment_dict_5})

	elif('go_enrichment' in request.POST):
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		group_for_enr = "both"
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		if('group_for_enr' in request.POST):
			group_for_enr = request.POST.get('group_for_enr')
			#print(pval_enr)
		result1 = run_go_enrichment_2.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_go")
		result2 = run_go_enrichment_2.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_go")
		result3 = run_go_enrichment_2.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_go")
		#result1 = run_go_enrichment.delay("genelist.txt",pval_enr)
		enr_results = result1.get()
		enr_results_2 = result2.get()
		enr_results_3 = result3.get()
		print("enr")
		result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
		result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
		result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
		enrichment_dict = result4.get()
		enrichment_dict_2 = result5.get()
		enrichment_dict_3 = result6.get()
		return render(request,'polls/clustering_6.html',{'list_of_files':list_of_files,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3})
	elif('go_molecular' in request.POST):
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		#result1 = run_go_enrichment.delay("genelist.txt",pval_enr)
		group_for_enr = "both"
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		if('group_for_enr' in request.POST):
			group_for_enr = request.POST.get('group_for_enr')
			#print(pval_enr)
		result1 = run_go_enrichment_2.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_go")
		result2 = run_go_enrichment_2.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_go")
		result3 = run_go_enrichment_2.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_go")
		enr_results = result1.get()
		enr_results_2 = result2.get()
		enr_results_3 = result3.get()
		print("enr")
		result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
		result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
		result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
		enrichment_dict = result4.get()
		enrichment_dict_2 = result5.get()
		enrichment_dict_3 = result6.get()
		return render(request,'polls/clustering_6.html',{'list_of_files':list_of_files,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3})
	elif('reac_enrichment' in request.POST):
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		#result1 = run_go_enrichment.delay("genelist.txt",pval_enr)
		group_for_enr = "both"
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		if('group_for_enr' in request.POST):
			group_for_enr = request.POST.get('group_for_enr')
			#print(pval_enr)
		#result1 = run_reac_enrichment.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_reactome")
		result2 = run_reac_enrichment.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_reactome")
		#result3 = run_reac_enrichment.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_reactome")
		#enr_results = result1.get()
		enr_results_2 = result2.get()
		#enr_results_3 = result3.get()
		print("enr")
		#result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
		result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_reactome/Reactome_2016.test_name.enrichr.reports.txt",pval_enr)
		#result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
		#enrichment_dict = result4.get()
		enrichment_dict = {}
		enrichment_dict_2 = result5.get()
		#enrichment_dict_3 = result6.get()
		enrichment_dict_3 = {}
		return render(request,'polls/clustering_6_part_3.html',{'list_of_files':list_of_files,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3})

	#if not(('myfile' in request.FILES and 'protfile' in request.FILES) or ('predef_file' in request.POST)):
	else:		
		#remove_loading_image.delay()
		ret_metadata = ""
		if not (request.user.is_authenticated):
			cache.clear()
		else:
			username = str(request.user)
			list_of_files = GraphForm.list_user_data_2(username)
			list_of_files_2 = GraphForm.list_user_data(username)
			cache.clear()
			#metd = list_metadata.apply_async(countdown=0)
			metd = list_metadata_3.apply_async(countdown=0)
			(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get() 
			print(ret_metadata1) 
			print("iubaerb")
			metadata_dict = [ret_metadata1,ret_metadata2,ret_metadata3]
			result2 = read_kegg_enrichment.delay("polls/data/test/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			#result2 = read_kegg_enrichment.delay("data/test/enrichr_kegg")	
			enrichment_dict = result2.get()
			#(ret_metadata,is_lungc) = metd.get()  
			#return render(request,'polls/clustering_6.html',{'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'ret_metadata':ret_metadata,'is_lungc':is_lungc})

		#make_empty_figure.apply_async(countdown=2)
		#empty_log_file.apply_async(countdown=2)
		return render(request,'polls/clustering_6_part_3.html',{'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'ret_metadata':ret_metadata,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'metadata_dict':metadata_dict,'enrichment_dict':enrichment_dict})



def clustering_6_part_3_2(request):
	print("basdbasdbasdbasdbs")
	ret_metadata1 = {}
	ret_metadata2 = {}
	ret_metadata3 = {}
	metadata_dict = []
	enrichment_dict = []
	pval_enr = 0.5
	list_of_files = ""
	list_of_files_2 = ""
	save_data = request.POST.get("save_data", None)
	nbr_iter = request.POST.get("nbr_iter",45)
	nbr_ants = request.POST.get("nbr_ants",30)
	evap = request.POST.get("evap",0.3)
	epsilon = request.POST.get("stopcr",0.02) 
	hi_sig = request.POST.get("hisig",1)
	pher_sig = request.POST.get("pher",1)
	savedata_param = "false"
	#fig = plt.figure(figsize=(10,8))
	#plt.savefig("polls/static/progress.png")
	#plt.close(fig)
	if request.user.is_authenticated:
		username = str(request.user)
		list_of_files = GraphForm.list_user_data_2(username)
	if('analyze_metadata' in request.POST and 'myfile' in request.FILES and 'protfile' in request.FILES and 'patientdata' in request.FILES and (1 == 0)):
		if(request.FILES['myfile'] and request.FILES['protfile'] and request.FILES['patientdata']):
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				#return HttpResponseRedirect(reverse('clustering_6_step_2'))
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				add_loading_image.delay()
				with open("polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				make_empty_figure.delay()
				#fh2 = open(request.FILES['protfile'])
				ppistr = request.FILES['protfile'].read().decode('utf-8')
				#fh1 = open(request.FILES['myfile'])
				exprstr = request.FILES['myfile'].read().decode('utf-8')
				clinicalstr = request.FILES['patientdata'].read().decode('utf-8')
				clinical_stringio = StringIO(clinicalstr)
				clinicaldf = pd.read_csv(clinical_stringio)
				print(clinicaldf)
				survival_col_name = ""
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			savedata_param = "true"
				if('survival_col' in request.POST):
					print("barabsfrbasdb")
					if(request.POST['survival_col']):
						survival_col_name = request.POST['survival_col']
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				#result2 = script_output_task.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1)
				#result2 = script_output_task_4.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2)
				#result2 = script_output_task_8.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name)
				#(div,script,plot1,plot_div,ret_metadata) = result2.get()
				#result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				#(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				#print("ikkrngasndgnn")
				#result2 = script_output_task_4.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2)
				result2 = script_output_task_9.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,survival_col_name,clinicaldf)
				(div,script,plot1,plot_div,ret_metadata) = result2.get()
				metd = list_metadata_3.apply_async(countdown=0)
				(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get()
				print("iasdfasdfsf")
				print(ret_metadata1) 	
				#(div,script,plot1) = result2.get()			
				plot2 = "test.png"
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			savedata_param = "true"
	                			username = str(request.user)
	                			GraphForm.save_user_data_2(request.FILES['myfile'],request.FILES['protfile'],request.FILES['patientdata'],username)
						#GraphForm.save_results(username)
				#else:
					#cache.clear()
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)              
				#else:
				remove_loading_image.delay()
				cache.clear()				
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				session_id = request.session._get_or_create_session_key()
				path99 = "/home/quirin/testproject/polls/static/test_" + session_id + ".png"
				path_metadata = "polls/static/metadata_" + session_id + ".txt"
				output_plot_path = "polls/static/output_plotly_" + session_id + ".html"
				#list_metadata.apply_async(countdown=0)
				with open("polls/static/plotly_output_2.html", "w") as text_file_3:
   					text_file_3.write(plot_div)
				return render(request, 'polls/clustering_6_part_3.html', {'form':"",'images':"",'plot_div':plot_div,'script':script,'plot2':plot2, 'list_of_files':list_of_files,'ret_dat':ret_metadata,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'path4':path99,'output_plot_path':output_plot_path,'list_of_files_2':list_of_files_2})
	elif('kegg_enrichment' in request.POST):
		group_for_enr = "both"
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		if('group_for_enr' in request.POST):
			group_for_enr = request.POST.get('group_for_enr')
			#print(pval_enr)
		result1 = run_enrichment_2.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_kegg")
		result2 = run_enrichment_2.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_kegg")
		result3 = run_enrichment_2.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_kegg")
		enr_results = result1.get()
		enr_results_2 = result2.get()
		enr_results_3 = result3.get()
		print("enr")
		result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
		result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
		result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
		enrichment_dict = result4.get()
		enrichment_dict_2 = result5.get()
		enrichment_dict_3 = result6.get()
		return render(request,'polls/clustering_6_part_3.html',{'list_of_files':list_of_files,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3})
	elif('enrichment_type' in request.POST):
		enr_type = request.POST.get("enrichment_type")
		group_for_enr = "both"
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		if('group_for_enr' in request.POST):
			group_for_enr = request.POST.get('group_for_enr')
			#print(pval_enr)
		enrichment_dict = {}
		enrichment_dict_2 = {}
		enrichment_dict_3 = {}
		enrichment_dict_4 = {}
		enrichment_dict_5 = {}
		if(enr_type == "kegg_enrichment"):
			result1 = run_enrichment_2.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_kegg")
			result2 = run_enrichment_2.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_kegg")
			result3 = run_enrichment_2.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_kegg")
			enr_results = result1.get()
			enr_results_2 = result2.get()
			enr_results_3 = result3.get()
			print("enr")
			result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			result7 = read_kegg_enrichment_2.delay("polls/data/test2/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt","polls/data/test3/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			enrichment_dict = result4.get()
			enrichment_dict_2 = result5.get()
			enrichment_dict_3 = result6.get()
			(enrichment_dict_4,enrichment_dict_5) = result7.get()
		elif(enr_type == "go_enrichment"):	
			result1 = run_go_enrichment_2.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_go")
			result2 = run_go_enrichment_2.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_go")
			result3 = run_go_enrichment_2.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_go")
			#result1 = run_go_enrichment.delay("genelist.txt",pval_enr)
			enr_results = result1.get()
			enr_results_2 = result2.get()
			enr_results_3 = result3.get()
			print("enr")
			result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			result7 = read_kegg_enrichment_2.delay("polls/data/test2/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt","polls/data/test3/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
			enrichment_dict = result4.get()
			enrichment_dict_2 = result5.get()
			enrichment_dict_3 = result6.get()	
			(enrichment_dict_4,enrichment_dict_5) = result7.get()
		elif(enr_type == "go_molecular"):
			result1 = run_go_enrichment_2.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_go")
			result2 = run_go_enrichment_2.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_go")
			result3 = run_go_enrichment_2.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_go")
			enr_results = result1.get()
			enr_results_2 = result2.get()
			enr_results_3 = result3.get()
			print("enr")
			result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result7 = read_kegg_enrichment_2.delay("polls/data/test2/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt","polls/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			enrichment_dict = result4.get()
			enrichment_dict_2 = result5.get()
			enrichment_dict_3 = result6.get()
			(enrichment_dict_4,enrichment_dict_5) = result7.get()
		elif(enr_type == "reactome_enrichment"):
			#result1 = run_reac_enrichment.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_reactome")
			result2 = run_reac_enrichment.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_reactome")
			#result3 = run_reac_enrichment.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_reactome")
			#enr_results = result1.get()
			enr_results_2 = result2.get()
			#enr_results_3 = result3.get()
			print("enr")
			#result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_reactome/Reactome_2016.test_name.enrichr.reports.txt",pval_enr)
			#result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
			#enrichment_dict = result4.get()
			enrichment_dict = {}
			enrichment_dict_2 = result5.get()
			#enrichment_dict_3 = result6.get()
			enrichment_dict_3 = {}		
		return render(request,'polls/clustering_6_part_3.html',{'list_of_files':list_of_files,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3,'enrichment_dict_4':enrichment_dict_4,'enrichment_dict_5':enrichment_dict_5})

	elif('go_enrichment' in request.POST):
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		group_for_enr = "both"
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		if('group_for_enr' in request.POST):
			group_for_enr = request.POST.get('group_for_enr')
			#print(pval_enr)
		result1 = run_go_enrichment_2.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_go")
		result2 = run_go_enrichment_2.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_go")
		result3 = run_go_enrichment_2.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_go")
		#result1 = run_go_enrichment.delay("genelist.txt",pval_enr)
		enr_results = result1.get()
		enr_results_2 = result2.get()
		enr_results_3 = result3.get()
		print("enr")
		result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
		result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
		result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_go/GO_Biological_Process_2018.test_name.enrichr.reports.txt",pval_enr)
		enrichment_dict = result4.get()
		enrichment_dict_2 = result5.get()
		enrichment_dict_3 = result6.get()
		return render(request,'polls/clustering_6_part_3.html',{'list_of_files':list_of_files,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3})
	elif('go_molecular' in request.POST):
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		#result1 = run_go_enrichment.delay("genelist.txt",pval_enr)
		group_for_enr = "both"
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		if('group_for_enr' in request.POST):
			group_for_enr = request.POST.get('group_for_enr')
			#print(pval_enr)
		result1 = run_go_enrichment_2.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_go")
		result2 = run_go_enrichment_2.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_go")
		result3 = run_go_enrichment_2.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_go")
		enr_results = result1.get()
		enr_results_2 = result2.get()
		enr_results_3 = result3.get()
		print("enr")
		result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
		result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
		result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
		enrichment_dict = result4.get()
		enrichment_dict_2 = result5.get()
		enrichment_dict_3 = result6.get()
		return render(request,'polls/clustering_6_part_3.html',{'list_of_files':list_of_files,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3})
	elif('reac_enrichment' in request.POST):
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		#result1 = run_go_enrichment.delay("genelist.txt",pval_enr)
		group_for_enr = "both"
		if('pval_enr' in request.POST):
			pval_enr = request.POST.get('pval_enr')
			print(pval_enr)
		if('group_for_enr' in request.POST):
			group_for_enr = request.POST.get('group_for_enr')
			#print(pval_enr)
		#result1 = run_reac_enrichment.delay("genelist_1.txt",pval_enr,"polls/data/test/enrichr_reactome")
		result2 = run_reac_enrichment.delay("genelist_2.txt",pval_enr,"polls/data/test2/enrichr_reactome")
		#result3 = run_reac_enrichment.delay("genelist.txt",pval_enr,"polls/data/test3/enrichr_reactome")
		#enr_results = result1.get()
		enr_results_2 = result2.get()
		#enr_results_3 = result3.get()
		print("enr")
		#result4 = read_kegg_enrichment.delay("polls/data/test/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
		result5 = read_kegg_enrichment.delay("polls/data/test2/enrichr_reactome/Reactome_2016.test_name.enrichr.reports.txt",pval_enr)
		#result6 = read_kegg_enrichment.delay("polls/data/test3/enrichr_go/GO_Molecular_Function_2018.test_name.enrichr.reports.txt",pval_enr)
		#enrichment_dict = result4.get()
		enrichment_dict = {}
		enrichment_dict_2 = result5.get()
		#enrichment_dict_3 = result6.get()
		enrichment_dict_3 = {}
		return render(request,'polls/clustering_6_part_3.html',{'list_of_files':list_of_files,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'enrichment_dict':enrichment_dict,'enrichment_dict_2':enrichment_dict_2,'enrichment_dict_3':enrichment_dict_3})

	#if not(('myfile' in request.FILES and 'protfile' in request.FILES) or ('predef_file' in request.POST)):
	else:		
		#remove_loading_image.delay()
		ret_metadata = ""
		if not (request.user.is_authenticated):
			cache.clear()
		else:
			username = str(request.user)
			list_of_files = GraphForm.list_user_data_2(username)
			list_of_files_2 = GraphForm.list_user_data(username)
			cache.clear()
			#metd = list_metadata.apply_async(countdown=0)
			#metd = list_metadata_3.apply_async(countdown=0)
			session_id = request.session._get_or_create_session_key()
			path_metadata = "polls/static/metadata_" + session_id + ".txt"
			metd = list_metadata_4.apply_async(args=[path_metadata],countdown=0)			
			(ret_metadata1,ret_metadata2,ret_metadata3) = metd.get() 
			print(ret_metadata1) 
			print("iubaerb")
			metadata_dict = [ret_metadata1,ret_metadata2,ret_metadata3]
			result2 = read_kegg_enrichment.delay("polls/data/test/enrichr_kegg/KEGG_2013.test_name.enrichr.reports.txt",pval_enr)
			#result2 = read_kegg_enrichment.delay("data/test/enrichr_kegg")	
			enrichment_dict = result2.get()
			#(ret_metadata,is_lungc) = metd.get()  
			#return render(request,'polls/clustering_6.html',{'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'ret_metadata':ret_metadata,'is_lungc':is_lungc})
			print("clustering_6_part_3_2")
		session_id = request.session._get_or_create_session_key()
		path99 = "/home/quirin/testproject/polls/static/test_" + session_id + ".png"
		path_metadata = "polls/static/metadata_" + session_id + ".txt"
		output_plot_path = "polls/static/output_plotly_" + session_id + ".html"
		return render(request,'polls/clustering_6_part_3.html',{'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'ret_metadata':ret_metadata,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'path4':path99,'output_plot_path':output_plot_path,'metadata_dict':metadata_dict,'enrichment_dict':enrichment_dict})

			#return render(request,'polls/clustering_6_part_3.html',{'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'ret_metadata':ret_metadata,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'metadata_dict':metadata_dict,'enrichment_dict':enrichment_dict})

		#make_empty_figure.apply_async(countdown=2)
		#empty_log_file.apply_async(countdown=2) 				
		return render(request,'polls/clustering_6_part_3.html',{'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'ret_metadata':ret_metadata,'ret_metadata1':ret_metadata1,'ret_metadata2':ret_metadata2,'ret_metadata3':ret_metadata3,'metadata_dict':metadata_dict,'enrichment_dict':enrichment_dict})




def clustering_7(request):
	list_of_files = ""
	list_of_files_2 = ""
	save_data = request.POST.get("save_data", None)
	nbr_iter = request.POST.get("nbr_iter",45)
	nbr_ants = request.POST.get("nbr_ants",30)
	evap = request.POST.get("evap",0.3)
	epsilon = request.POST.get("stopcr",0.02) 
	hi_sig = request.POST.get("hisig",1)
	pher_sig = request.POST.get("pher",1)
	savedata_param = "false"
	#fig = plt.figure(figsize=(10,8))
	#plt.savefig("polls/static/progress.png")
	#plt.close(fig)
	if request.user.is_authenticated:
		username = str(request.user)
		list_of_files = GraphForm.list_user_data_2(username)
	if('redo_analysis' in request.POST and request.user.is_authenticated):
		if(request.POST['redo_analysis']):
			with open("polls/static/output_console.txt", "w") as text_file:
   				text_file.write("Your request is being processed...")
			add_loading_image.delay()
			filename1 = request.POST.get("input_own_file_redo")
			fh1 = open(filename1)
			exprstr = fh1.read()
			filename2 = filename1.split("_expr.txt")[0] + "_prot.txt"
			fh2 = open(filename2)
			filename3 = filename1.split("_expr.txt")[0] + "_clin.txt"
			has_clin_data = "false"
			if(os.path.isfile(filename3)):
				fh3 = open(filename3)
				has_clin_data = "true"
				clinicalstr = fh3.read()
			#path_json = filename1
			#path_heatmap = filename1.split("_json.json")[0] + "_heatmap.png"
			ppistr = fh2.read()
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				make_empty_figure.delay()
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])	
				plot_div = ""
				##(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1) = algo_output_ownfile_2(1,lgmin,lgmax,fh1,fh2,nbr_iter,nbr_ants,evap)
				##(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1) = algo_output_ownfile_3(1,lgmin,lgmax,fh1,fh2,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				##(T,row_colors,col_colors,G2,means,adjlist,genes1) = algo_output_ownfile(1,lgmin,lgmax,fh1,fh2)
				##(div,script,plot1) = GraphForm.handle_script_output(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1)
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids) =result1.get()				
				if(has_clin_data == "true"):
					result2 = script_output_task_3.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr)
					(div,script,plot1,plot_div,ret_metadata) = result2.get()
				else:
					result2 = script_output_task.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1)
					(div,script,plot1) = result2.get()
					ret_metadata = ""
				plot2 = "test.png"
				cache.clear()			
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				list_of_files = GraphForm.list_user_data_2(username)	
				list_of_files_2 = GraphForm.list_user_data(username)      
				remove_loading_image.delay()
				#return render(request, 'polls/clustering_6.html', {'form':"",'images':"",'div':div,'script':script,'plot2':plot2, 'list_of_files':list_of_files})       		
				#show_old_data(path_json,path_heatmap)             	
				#plot2 = "test.png"
				#cache.clear()			
				#make_empty_figure.apply_async(countdown=10)
				#empty_log_file.apply_async(countdown=10)
				return render(request, 'polls/clustering_7.html', {'list_of_files':list_of_files,'ret_metadata':ret_metadata,'result_99':ret_metadata,'list_of_files_2':list_of_files_2,'json_path':"test15.json",'heatmap_path':"test.png",'survival_graph':plot_div})
   	
	elif('myfile' in request.FILES and 'ndexfile' in request.FILES):
		if(request.FILES['myfile'] and request.FILES['ndexfile']):
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				add_loading_image.delay()
				with open("polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				make_empty_figure.delay()
				#fh2 = open(request.FILES['protfile'])
				ppistr = read_ndex_file(request.FILES['ndexfile'])
				#fh1 = open(request.FILES['myfile'])
				exprstr = request.FILES['myfile'].read().decode('utf-8')
				#(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1) = algo_output_ownfile_3(1,lgmin,lgmax,request.FILES['myfile'],request.FILES['protfile'],nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				#(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1) = algo_output_ownfile_2(1,lgmin,lgmax,request.FILES['myfile'],request.FILES['protfile'],nbr_iter,nbr_ants,evap)
				#(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1) = run_analysis.delay(1,lgmin,lgmax,request.FILES['myfile'],request.FILES['protfile'],nbr_iter,nbr_ants,evap)
				#(T,row_colors,col_colors,G2,means,adjlist,genes1) = algo_output_ownfile(1,lgmin,lgmax,request.FILES['myfile'],request.FILES['protfile'])
				#(div,script,plot1) = GraphForm.handle_script_output(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1)
				savedata_param = "false"
				plot_div = ""
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			savedata_param = "true"
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				#result2 = script_output_task.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1)
				clinicalstr = "empty"
				if(request.user.is_authenticated and savedata_param == "true"):
					result2 = script_output_task_2.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,savedata_param,str(request.user))
					(div,script,plot1) = result2.get()			
				else:				
					result2 = script_output_task_4.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2)
					(div,script,plot1,plot_div,ret_metadata) = result2.get()
				plot2 = "test.png"
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			savedata_param = "true"
	                			username = str(request.user)
	                			#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
				#else:
					#cache.clear()
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)      	        
				#else:
				remove_loading_image.delay()
				cache.clear()				
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				return render(request, 'polls/clustering_7.html', {'list_of_files':list_of_files,'ret_metadata':ret_metadata,'result_99':ret_metadata,'list_of_files_2':list_of_files_2,'json_path':"test15.json",'heatmap_path':"test.png",'survival_graph':plot_div})
	
	
	elif('own_file' in request.POST and request.user.is_authenticated):
		if(request.POST['own_file']):
			make_empty_figure.delay()
			with open("polls/static/output_console.txt", "w") as text_file:
   				text_file.write("Your request is being processed...")
			filename1 = request.POST.get("input_own_file")
			#fh1 = open(filename1)
			#exprstr = fh1.read()
			#filename2 = filename1.split("_expr.txt")[0] + "_prot.txt"
			#fh2 = open(filename2)
			path_json = filename1
			path_heatmap = filename1.split("_json.json")[0] + "_heatmap.png"
			#ppistr = fh2.read()
			list_of_files = GraphForm.list_user_data_2(username)	
			list_of_files_2 = GraphForm.list_user_data(username)      
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])	
				##(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1) = algo_output_ownfile_2(1,lgmin,lgmax,fh1,fh2,nbr_iter,nbr_ants,evap)
				##(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1) = algo_output_ownfile_3(1,lgmin,lgmax,fh1,fh2,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				##(T,row_colors,col_colors,G2,means,adjlist,genes1) = algo_output_ownfile(1,lgmin,lgmax,fh1,fh2)
				##(div,script,plot1) = GraphForm.handle_script_output(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1)
				#result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				#(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1) =result1.get()				
				#result2 = script_output_task.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1)
				#(div,script,plot1) = result2.get()
				#plot2 = "test.png"
				#cache.clear()			
				#make_empty_figure.apply_async(countdown=10)
				#empty_log_file.apply_async(countdown=10)
				#return render(request, 'polls/clustering_6.html', {'form':"",'images':"",'div':div,'script':script,'plot2':plot2, 'list_of_files':list_of_files})       		
				show_old_data(path_json,path_heatmap)             	
				plot2 = "test.png"
				cache.clear()			
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				return render(request, 'polls/clustering_6.html', {'form':"",'images':"",'div':"",'script':"",'plot2':plot2, 'list_of_files':list_of_files,'list_of_files_2':list_of_files_2})       		
	elif('analyze_metadata' in request.POST and 'myfile' in request.FILES and 'protfile' in request.FILES and 'patientdata' in request.FILES):
		if(request.FILES['myfile'] and request.FILES['protfile'] and request.FILES['patientdata']):
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				add_loading_image.delay()
				with open("polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				make_empty_figure.delay()
				#fh2 = open(request.FILES['protfile'])
				ppistr = request.FILES['protfile'].read().decode('utf-8')
				#fh1 = open(request.FILES['myfile'])
				exprstr = request.FILES['myfile'].read().decode('utf-8')
				clinicalstr = request.FILES['patientdata'].read().decode('utf-8')
				#(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1) = algo_output_ownfile_3(1,lgmin,lgmax,request.FILES['myfile'],request.FILES['protfile'],nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				#(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1) = algo_output_ownfile_2(1,lgmin,lgmax,request.FILES['myfile'],request.FILES['protfile'],nbr_iter,nbr_ants,evap)
				#(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1) = run_analysis.delay(1,lgmin,lgmax,request.FILES['myfile'],request.FILES['protfile'],nbr_iter,nbr_ants,evap)
				#(T,row_colors,col_colors,G2,means,adjlist,genes1) = algo_output_ownfile(1,lgmin,lgmax,request.FILES['myfile'],request.FILES['protfile'])
				#(div,script,plot1) = GraphForm.handle_script_output(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1)
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			savedata_param = "true"
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				#result2 = script_output_task.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1)
				result2 = script_output_task_4.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2)
				(div,script,plot1,plot_div,ret_metadata) = result2.get()
				#(div,script,plot1) = result2.get()			
				plot2 = "test.png"
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			savedata_param = "true"
	                			username = str(request.user)
	                			GraphForm.save_user_data_2(request.FILES['myfile'],request.FILES['protfile'],request.FILES['patientdata'],username)
						#GraphForm.save_results(username)
				#else:
					#cache.clear()
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)              
				#else:
				remove_loading_image.delay()
				cache.clear()				
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				#list_metadata.apply_async(countdown=0)
				with open("polls/static/plotly_output_2.html", "w") as text_file_3:
   					text_file_3.write(plot_div)
				return render(request, 'polls/clustering_7.html', {'list_of_files':list_of_files,'ret_metadata':ret_metadata,'result_99':ret_metadata,'list_of_files_2':list_of_files_2,'json_path':"test15.json",'heatmap_path':"test.png",'survival_graph':plot_div})


	elif('myfile' in request.FILES and 'protfile' in request.FILES):
		if(request.FILES['myfile'] and request.FILES['protfile']):
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				add_loading_image.delay()
				with open("polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				make_empty_figure.delay()
				#fh2 = open(request.FILES['protfile'])
				ppistr = request.FILES['protfile'].read().decode('utf-8')
				#fh1 = open(request.FILES['myfile'])
				exprstr = request.FILES['myfile'].read().decode('utf-8')
				#(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1) = algo_output_ownfile_3(1,lgmin,lgmax,request.FILES['myfile'],request.FILES['protfile'],nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				#(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1) = algo_output_ownfile_2(1,lgmin,lgmax,request.FILES['myfile'],request.FILES['protfile'],nbr_iter,nbr_ants,evap)
				#(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1) = run_analysis.delay(1,lgmin,lgmax,request.FILES['myfile'],request.FILES['protfile'],nbr_iter,nbr_ants,evap)
				#(T,row_colors,col_colors,G2,means,adjlist,genes1) = algo_output_ownfile(1,lgmin,lgmax,request.FILES['myfile'],request.FILES['protfile'])
				#(div,script,plot1) = GraphForm.handle_script_output(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1)
				savedata_param = "false"
				plot_div = ""
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			savedata_param = "true"
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				#result2 = script_output_task.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1)
				clinicalstr = "empty"
				if(request.user.is_authenticated and savedata_param == "true"):
					result2 = script_output_task_2.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,savedata_param,str(request.user))
					(div,script,plot1) = result2.get()			
				else:				
					result2 = script_output_task_4.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2)
					(div,script,plot1,plot_div,ret_metadata) = result2.get()
				plot2 = "test.png"
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			savedata_param = "true"
	                			username = str(request.user)
	                			#GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
						#GraphForm.save_results(username)
				#else:
					#cache.clear()
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)      	        
				#else:
				remove_loading_image.delay()
				cache.clear()				
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				return render(request, 'polls/clustering_7.html', {'list_of_files':list_of_files,'ret_metadata':ret_metadata,'result_99':ret_metadata,'list_of_files_2':list_of_files_2,'json_path':"test15.json",'heatmap_path':"test.png",'survival_graph':plot_div})
	elif('predef_file' in request.POST):
		if('L_g_min' in request.POST and 'L_g_max' in request.POST):
			if(request.POST['predef_file']):
				make_empty_figure.delay()
				with open("polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				add_loading_image.delay()
				cancer_type = request.POST.get("cancer_type")
				print("Cancer Type:" + str(cancer_type))
				fh2 = open("polls/data/biogrid.human.entrez.tsv")
				ppistr = fh2.read()
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])	
				is_lungc = "false"		
				if(cancer_type == "1"):
					print("babababababa")
					fh1 = open("polls/data/lung_cancer_expr.csv")
					exprstr = fh1.read()
					#fh4 = open("polls/data/breast_cancer_clinical.csv")
					fh4 = open("polls/data/lung_cancer_clinical.csv")
					clinicalstr = fh4.read()
					is_lungc = "true"
					#clinicalLines = clinicalstr.split("\n")
					#print(clinicalLines[0])
					#bababa = clinicalLines[50].replace('"','')
					#bababa2 = bababa.split(",")
					#print(bababa2[64])
					fh4.flush()
					fh4.close()
				else:
					fh1 = open("polls/data/breast_cancer_expr.csv")
					exprstr = fh1.read()
					fh4 = open("polls/data/breast_cancer_clinical.csv")
					clinicalstr = fh4.read()
					#clinicalLines = clinicalstr.split("\n")
					#print(clinicalLines[0])
					#bababa = clinicalLines[50].replace('"','')
					#bababa2 = bababa.split(",")
					#print(bababa2[64])
					fh4.flush()
					fh4.close()
				#(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1) = algo_output_ownfile_2(1,lgmin,lgmax,fh1,fh2,nbr_iter,nbr_ants,evap)
				#(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1) = algo_output_ownfile_3(1,lgmin,lgmax,fh1,fh2,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				#(T,row_colors,col_colors,G2,means,adjlist,genes1) = algo_output_ownfile(1,lgmin,lgmax,fh1,fh2)
				#(div,script,plot1) = GraphForm.handle_script_output(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1)
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,jac_1,jac_2) =result1.get()				
				#result2 = script_output_task.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1)
				#result2 = script_output_task_3.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr)
				#(div,script,plot1,plot_div,ret_metadata) = result2.get()
				#result2 = script_output_task_4.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2)
				result2 = script_output_task_7.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1,jac_2,is_lungc)
				(div,script,plot1,plot_div,ret_metadata) = result2.get()
				print(ret_metadata)
				plot2 = "test.png"
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data_2(username)	
			        	list_of_files_2 = GraphForm.list_user_data(username)        				
				#else:
				remove_loading_image.delay()
				#with open("polls/static/metadata_99.html", "w") as text_file_4:
   				#	text_file_4.write(str(ret_metadata[0]))
   				#	text_file_4.write(str(ret_metadata[1]))
   				#	text_file_4.write(str(ret_metadata[2]))
				cache.clear()
				with open("polls/static/plotly_output_2.html", "w") as text_file_3:
   					text_file_3.write(plot_div)
				#write_metadata_to_file.apply_async((ret_metadata), countdown=1)
				#list_metadata.apply_async(countdown=0)
				#string_99 = metadata_to_string.apply_async(ret_metadata, countdown=0)
				#result_99 = string_99.get()
				return render(request, 'polls/clustering_6.html', {'list_of_files':list_of_files,'ret_metadata':ret_metadata,'result_99':ret_metadata,'list_of_files_2':list_of_files_2,'json_path':"test15.json",'heatmap_path':"test.png",'survival_graph':plot_div})
	                
	#if not(('myfile' in request.FILES and 'protfile' in request.FILES) or ('predef_file' in request.POST)):
	else:		
		#remove_loading_image.delay()
		ret_metadata = ""
		if not (request.user.is_authenticated):
			cache.clear()
		else:
			username = str(request.user)
			list_of_files = GraphForm.list_user_data_2(username)
			list_of_files_2 = GraphForm.list_user_data(username)
			cache.clear()
			metd = list_metadata.apply_async(countdown=0)
			print("foo")
			ret_metadata = metd.get()  
		#make_empty_figure.apply_async(countdown=2)
		#empty_log_file.apply_async(countdown=2)	 				
		return render(request,'polls/clustering_7.html',{'list_of_files':list_of_files,'list_of_files_2':list_of_files_2,'ret_metadata':""})


def infopage(request):			
	return render(request,'polls/infopage.html')

def list_2(request):
	ctr = 0
	div = ""
	script = ""
	plot2 = ""
	logstring = ""
	plot3 = ""
	genelist_diff = ""
	if(ctr == 0):
	        if request.method=="POST":
	                img = UploadForm(request.POST, request.FILES)
	                ctr = 1
	                filenames = request.FILES.items()
	                logstring = logstring + str(list(filenames))	
			#check if option "include survival" was crossed
	                display_type = request.POST.get("include_survival", None)
			#handle pre-defined file
	                if('pre_def_file' in request.POST):
	                	if(request.POST['pre_def_file']):
	                		fh1 = open("testgenes13.txt")
	                		fh2 = open("ppi_test_2.txt")
	                		#(script,div,plot2) = handle_upload_5(fh1,fh2)
	                		(script,div,plot2,plot3) = GraphForm.handle_upload_7(fh1,fh2)
	                		#(script,div,plot2,plot3) = handle_upload_5(fh1,fh2)
			#handle manually uploaded files
	                if('myfile' in request.FILES and 'prot_file' in request.FILES):
	                	if(request.FILES['myfile'] and request.FILES['prot_file']):
	                		if display_type in ["include_survival"]:
	                			#(script,div,plot2,plot3) = GraphForm.handle_upload_3_survival(request.FILES['myfile'],request.FILES['prot_file'])
	                			(script,div,plot2,plot3,genelist_diff) = GraphForm.handle_upload_3_survival_3(request.FILES['myfile'],request.FILES['prot_file'])
	                		else:
	                			(script,div,plot2) = GraphForm.handle_upload_3_4(request.FILES['myfile'],request.FILES['prot_file'])

	                if img.is_valid():
	                        img.save()
	        else:
	                img=UploadForm()
	                ctr = 1
	        images=Upload.objects.all()
	        plot2 = "test.png"
	        return render(request,'polls/list.html',{'form':img,'images':images, 'div':div, 'script':script,'plot2':plot2, 'plot_2':plot3, 'genelist_diff':genelist_diff})

#those are all old methods, inherited from the django template. i leave them here because i used to play around with question and polls.
def index(request):
    latest_question_list = Question.objects.order_by('-pub_date')[:5]
    context = {'latest_question_list': latest_question_list}
    return render(request, 'polls/index.html', context)

def detail(request, question_id):
    try:
        question = Question.objects.get(pk=question_id)
    except Question.DoesNotExist:
        raise Http404("Question does not exist")
    return render(request, 'polls/detail.html', {'question': question})
    #return HttpResponse("You're looking at question %s." % question_id)
def results(request, question_id):
    question = get_object_or_404(Question, pk=question_id)
    return render(request, 'polls/results.html', {'question': question})
def vote(request, question_id):
    question = get_object_or_404(Question, pk=question_id)
    try:
        selected_choice = question.choice_set.get(pk=request.POST['choice'])
    except (KeyError, Choice.DoesNotExist):
        # Redisplay the question voting form.
        return render(request, 'polls/detail.html', {
            'question': question,
            'error_message': "You didn't select a choice.",
        })
    else:
        selected_choice.votes += 1
        selected_choice.save()
        # Always return an HttpResponseRedirect after successfully dealing
        # with POST data. This prevents data from being posted twice if a
        # user hits the Back button.
        return HttpResponseRedirect(reverse('polls:results', args=(question.id,)))

class IndexView(generic.ListView):
    template_name = 'polls/index.html'
    context_object_name = 'latest_question_list'

    def get_queryset(self):
        """Return the last five published questions."""
        return Question.objects.order_by('-pub_date')[:5]


class DetailView(generic.DetailView):
    model = Question
    template_name = 'polls/detail.html'


class ResultsView(generic.DetailView):
    model = Question
    template_name = 'polls/results.html'





def clustering_6_OLD(request):
	list_of_files = ""
	save_data = request.POST.get("save_data", None)
	nbr_iter = request.POST.get("nbr_iter",45)
	nbr_ants = request.POST.get("nbr_ants",30)
	evap = request.POST.get("evap",0.3)
	epsilon = request.POST.get("stopcr",0.02) 
	hi_sig = request.POST.get("hisig",1)
	pher_sig = request.POST.get("pher",1)
	#fig = plt.figure(figsize=(10,8))
	#plt.savefig("polls/static/progress.png")
	#plt.close(fig)
	if('myfile' in request.FILES and 'protfile' in request.FILES):
		if(request.FILES['myfile'] and request.FILES['protfile']):
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])
				with open("polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				make_empty_figure.delay()
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1) = algo_output_ownfile_3(1,lgmin,lgmax,request.FILES['myfile'],request.FILES['protfile'],nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				#(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1) = algo_output_ownfile_2(1,lgmin,lgmax,request.FILES['myfile'],request.FILES['protfile'],nbr_iter,nbr_ants,evap)
				#(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1) = run_analysis.delay(1,lgmin,lgmax,request.FILES['myfile'],request.FILES['protfile'],nbr_iter,nbr_ants,evap)
				#(T,row_colors,col_colors,G2,means,adjlist,genes1) = algo_output_ownfile(1,lgmin,lgmax,request.FILES['myfile'],request.FILES['protfile'])
				(div,script,plot1) = GraphForm.handle_script_output(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1)
				plot2 = "test.png"
				if save_data in ["save_data"]:
	                		if request.user.is_authenticated:
	                			username = str(request.user)
	                			GraphForm.save_user_data(request.FILES['myfile'],request.FILES['protfile'],username)
				#else:
					#cache.clear()
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data(username)	        
				#else:
				cache.clear()				
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				return render(request, 'polls/clustering_6.html', {'form':"",'images':"",'div':div,'script':script,'plot2':plot2, 'list_of_files':list_of_files})
	elif('predef_file' in request.POST):
		if('L_g_min' in request.POST and 'L_g_max' in request.POST):
			if(request.POST['predef_file']):
				make_empty_figure.delay()
				with open("polls/static/output_console.txt", "w") as text_file:
   					text_file.write("Your request is being processed...")
   					text_file.close()
				cancer_type = request.POST.get("cancer_type")
				fh2 = open("polls/data/biogrid.human.entrez.tsv")
				ppistr = fh2.read()
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])				
				if(cancer_type == 1):
					fh1 = open("polls/data/lung_cancer_expr.csv")
					exprstr = fh1.read()
				else:
					fh1 = open("polls/data/breast_cancer_expr.csv")
					exprstr = fh1.read()
				#(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1) = algo_output_ownfile_2(1,lgmin,lgmax,fh1,fh2,nbr_iter,nbr_ants,evap)
				#(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1) = algo_output_ownfile_3(1,lgmin,lgmax,fh1,fh2,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				#(T,row_colors,col_colors,G2,means,adjlist,genes1) = algo_output_ownfile(1,lgmin,lgmax,fh1,fh2)
				#(div,script,plot1) = GraphForm.handle_script_output(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1)
				result1 = algo_output_task.delay(1,lgmin,lgmax,exprstr,ppistr,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1) =result1.get()				
				result2 = script_output_task.delay(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1)
				(div,script,plot1) = result2.get()
				plot2 = "test.png"
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				if request.user.is_authenticated:
			        	username = str(request.user)
			        	list_of_files = GraphForm.list_user_data(username)	        				
				#else:			
				cache.clear()
				return render(request, 'polls/clustering_6.html', {'form':"",'images':"",'div':div,'script':script,'plot2':plot2, 'list_of_files':list_of_files})
	elif('own_file' in request.POST and request.user.is_authenticated):
		if(request.POST['own_file']):
			make_empty_figure.delay()
			with open("polls/static/output_console.txt", "w") as text_file:
   				text_file.write("Your request is being processed...")
			filename1 = request.POST.get("input_own_file")
			fh1 = open(filename1)
			filename2 = filename1.split("_expr.txt")[0] + "_prot.txt"
			fh2 = open(filename2)
			if('L_g_min' in request.POST and 'L_g_max' in request.POST):
				lgmin = int(request.POST['L_g_min'])
				lgmax = int(request.POST['L_g_max'])	
				#(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1) = algo_output_ownfile_2(1,lgmin,lgmax,fh1,fh2,nbr_iter,nbr_ants,evap)
				(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1) = algo_output_ownfile_3(1,lgmin,lgmax,fh1,fh2,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig)
				#(T,row_colors,col_colors,G2,means,adjlist,genes1) = algo_output_ownfile(1,lgmin,lgmax,fh1,fh2)
				(div,script,plot1) = GraphForm.handle_script_output(T,row_colors,col_colors,G2,means,genes_all,adjlist,genes1)
				plot2 = "test.png"
				cache.clear()			
				make_empty_figure.apply_async(countdown=10)
				empty_log_file.apply_async(countdown=10)
				return render(request, 'polls/clustering_6.html', {'form':"",'images':"",'div':div,'script':script,'plot2':plot2, 'list_of_files':list_of_files})       		
	                
	#if not(('myfile' in request.FILES and 'protfile' in request.FILES) or ('predef_file' in request.POST)):
	else:	
		if not (request.user.is_authenticated):
			cache.clear()
		else:
			username = str(request.user)
			list_of_files = GraphForm.list_user_data(username)
			cache.clear()	        				
		return render(request,'polls/clustering_6.html',{'list_of_files':list_of_files})







