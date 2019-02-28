from django.shortcuts import render

# Create your views here.
from django.http import HttpResponse
from django.template import loader
from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from django.shortcuts import get_object_or_404, render
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

from sklearn import datasets, linear_model
from sklearn.metrics import mean_squared_error, r2_score
from django.shortcuts import render_to_response,render
from django.template import RequestContext
from django.http import HttpResponseRedirect
from django.urls import reverse
import polls
from polls.models import Document
from polls.forms import DocumentForm
from polls.models import Upload,UploadForm,GraphForm
import numpy as np
import matplotlib.pyplot as plt
import mpld3

import seaborn as sns
import pandas as pd
from numpy import array

import matplotlib.patches as mpatches
###NOTE: ZEIGT NOCH NICHT DAS RICHTIGE BILD AN. ADRESSE DES BILDES ÄNDERN.

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

def homepage(request):
    null = ""
    context = {'latest_question_list': null}
    return render(request, 'polls/homepage.html', context)

def prediction(request):
	ctr = 0
	div = ""
	script = ""
	plot2 = ""
	logstring = ""
	plot3 = ""
	results = ""
	if(ctr == 0):
	        if request.method=="POST":
	                img = UploadForm(request.POST, request.FILES)
	                ctr = 1
	                filenames = request.FILES.items()
	                logstring = logstring + str(list(filenames))
	                display_type = request.POST.get("linear_prediction", None)
	                if('pre_def_file' in request.POST):
	                	if(request.POST['pre_def_file']):
	                		fh1 = open("testgenes13.txt")
	                		fh2 = open("ppi_test_2.txt")
	                		#(script,div,plot2) = handle_upload_5(fh1,fh2)
	                		results = handle_upload_2(fh1,fh2)
	                if('myfile' in request.FILES):
	                	if(request.FILES['myfile']):
	                		if display_type in ["linear_prediction"]:
	                			results = handle_upload_2(request.FILES['myfile'])
	                		else:
	                			results = handle_upload(request.FILES['myfile'])

	                #for filename, file in request.FILES.items():
	                        #handle_uploaded_file(request.FILES[filename])
	                        #if(ctr == 0):
	                        #handle_upload_1(request.FILES[filename])
	                #        (script,div,plot2) = handle_upload_3(request.FILES[filename])
	                        #return render(request, 'home.html', {'div': div, 'script': script})
	                        #ctr = 1		
                        #handle_upload_2(request.FILES[filename])
                        #makeplot()
	                if img.is_valid():
	                        #handle_uploaded_file(request.FILES['file'])
	                        img.save()
	                        #return HttpResponseRedirect(reverse('imageupload'))
	        else:
	                img=UploadForm()
	                ctr = 1
	        images=Upload.objects.all()
	        #plot2 = "test.png"
	        return render(request,'polls/predictions_3.html',{'results':results})



def prediction_linear(request):
	ctr = 0
	div = ""
	script = ""
	plot2 = ""
	logstring = ""
	plot3 = ""
	results = ""
	if(ctr == 0):
	        if request.method=="POST":
	                img = UploadForm(request.POST, request.FILES)
	                ctr = 1
	                filenames = request.FILES.items()
	                logstring = logstring + str(list(filenames))
	                display_type = request.POST.get("include_survival", None)
	                if('pre_def_file' in request.POST):
	                	if(request.POST['pre_def_file']):
	                		fh1 = open("testgenes13.txt")
	                		fh2 = open("ppi_test_2.txt")
	                		#(script,div,plot2) = handle_upload_5(fh1,fh2)
	                		results = handle_upload_2(fh1,fh2)
	                if('myfile' in request.FILES and 'prot_file' in request.FILES):
	                	if(request.FILES['myfile'] and request.FILES['prot_file']):
	                		if display_type in ["include_survival"]:
	                			results = handle_upload(request.FILES['myfile'],request.FILES['prot_file'])
	                		else:
	                			results = handle_upload(request.FILES['myfile'],request.FILES['prot_file'])

	                #for filename, file in request.FILES.items():
	                        #handle_uploaded_file(request.FILES[filename])
	                        #if(ctr == 0):
	                        #handle_upload_1(request.FILES[filename])
	                #        (script,div,plot2) = handle_upload_3(request.FILES[filename])
	                        #return render(request, 'home.html', {'div': div, 'script': script})
	                        #ctr = 1		
                        #handle_upload_2(request.FILES[filename])
                        #makeplot()
	                if img.is_valid():
	                        #handle_uploaded_file(request.FILES['file'])
	                        img.save()
	                        #return HttpResponseRedirect(reverse('imageupload'))
	        else:
	                img=UploadForm()
	                ctr = 1
	        images=Upload.objects.all()
	        #plot2 = "test.png"
	        return render(request,'polls/prediction.html',{'results':results})


def prediction_linear_corr(request):
	ctr = 0
	div = ""
	script = ""
	plot2 = ""
	logstring = ""
	plot3 = ""
	results = ""
	if(ctr == 0):
	        if request.method=="POST":
	                img = UploadForm(request.POST, request.FILES)
	                ctr = 1
	                filenames = request.FILES.items()
	                logstring = logstring + str(list(filenames))
	                display_type = request.POST.get("include_survival", None)
	                if('pre_def_file' in request.POST):
	                	if(request.POST['pre_def_file']):
	                		fh1 = open("testgenes13.txt")
	                		fh2 = open("ppi_test_2.txt")
	                		#(script,div,plot2) = handle_upload_5(fh1,fh2)
	                		results = handle_upload_3(fh1,fh2)
	                if('myfile' in request.FILES and 'prot_file' in request.FILES):
	                	if(request.FILES['myfile'] and request.FILES['prot_file']):
	                		if display_type in ["include_survival"]:
	                			results = handle_upload_3(request.FILES['myfile'],request.FILES['prot_file'])
	                		else:
	                			results = handle_upload_3(request.FILES['myfile'],request.FILES['prot_file'])

	                #for filename, file in request.FILES.items():
	                        #handle_uploaded_file(request.FILES[filename])
	                        #if(ctr == 0):
	                        #handle_upload_1(request.FILES[filename])
	                #        (script,div,plot2) = handle_upload_3(request.FILES[filename])
	                        #return render(request, 'home.html', {'div': div, 'script': script})
	                        #ctr = 1		
                        #handle_upload_2(request.FILES[filename])
                        #makeplot()
	                if img.is_valid():
	                        #handle_uploaded_file(request.FILES['file'])
	                        img.save()
	                        #return HttpResponseRedirect(reverse('imageupload'))
	        else:
	                img=UploadForm()
	                ctr = 1
	        images=Upload.objects.all()
	        #plot2 = "test.png"
	        return render(request,'polls/prediction.html',{'results':results})

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





def handle_upload(fn):
	patients = []
	
	patients1 = []
	patients2 = []
	genes = []
	geneNames = []
	#data = {}
	data1 = {}
	data2 = {}
	group1 = []
	group2 = []
	group_labels1 = []
	group_labels2 = []
	group1_data = []
	group2_data = []
	#patient_ids = ['3dbe99d1-e3b8-4ee2-b6a8-2e2e12c6fbe9','6F4C8D30-47FB-47DF-9EB7-4E5881E3711E','95CEF916-5545-455B-920C-773A54FC7676','67C73260-A242-4BBA-87C5-D2302556DFF7','55262FCB-1B01-4480-	B322-36570430C917','3dbe99d1-e3b8-4ee2-b6a8-2e2e12c6fbe9','6F4C8D30-47FB-47DF-9EB7-4E5881E3711E','95CEF916-5545-455B-920C-773A54FC7676','67C73260-A242-4BBA-87C5-D2302556DFF7','55262FCB-1B01-4480-	B322-36570430C917']
	#patientfilename = 'nationwidechildrens.org_clinical_patient_brca.txt'
	
	patients.append([1,2,3,4])
	patients1.append(['1','2'])
	#patients.append(['3','4'])
	patients2.append(['3','4'])
	group_labels1.append([1,1])
	group_labels2.append([2,2])
	logstring = "Creating Plots for given input files... \n\n"
	logstring = "Reading gene expression data... \n"
	line_no = 0
	patient_ids = []
	survival = []
	survival_yrs = []
	data = []
	#group2.append(group_labels2)
	group1.append([1,1])
	group2.append([2,2])
	genes = ()
	logstring = logstring + "\n\n Matching gene and protein IDs... \n"
	
	red_patch = mpatches.Patch(color='red', label='Group1')
	blue_patch = mpatches.Patch(color='blue', label='Group2')
	#lut = dict(zip(set(endData[0]), sns.hls_palette(len(set(endData[0])), l=0.5, s=0.8)))
	#col_colors = pd.DataFrame(endData[0])[0].map(lut)
	#print(col_colors)
	#colors = np.array(['#BB0000','#BB0000','#0000BB','#0000BB'])
	#df9 = pd.DataFrame(data=endData[1:,0:],index=geneNames,columns=patients)
	#df2 = pd.DataFrame(data=endData[0,0:], index='',columns=patients)
	#my_palette = dict(zip(df[.unique(), ["orange","yellow","brown"]))
	#row_colors = df2.cyl.map(my_palette)
	#fig, (ax1, ax2) = plt.subplots(1,2,sharex=True,sharey=True)
	colordict={0:'#BB0000',1:'#0000BB'}
	#logstring = logstring + str(df9)
	df2 = pd.read_csv(fn,delim_whitespace=True,header=None,index_col=0)
	#print(df2.head())
	df = df2.transpose()
	print(df)
	survival_real = df['SURVIVAL']
	df['is_train'] = np.random.uniform(0, 1, len(df)) <= .75
	train, test = df[df['is_train']==True], df[df['is_train']==False]
	features = df.columns[3:60483]
	y = pd.factorize(train['SURVIVAL'])[0]
	print(y)
	classif = RandomForestClassifier(n_jobs=2, random_state=0)
	classif.fit(train[features], y)
	print(classif.predict(df[features]))
	df3 = df
	df3['SURVIVAL'] = classif.predict(df[features])
	df4 = df3.sort_values(by=['SURVIVAL'])
	#survival_real = df['SURVIVAL']
	survival_prediction = df3['SURVIVAL']
	predictions = []
	for j in range(1,len(survival_real)):
		accur = "FALSE"
		if(int(survival_real[j]) == int(survival_prediction[j])):
			print("Foo")
			accur = "TRUE"
		predictions.append({'patient_id':j,'real_value':survival_real[j],'prediction':survival_prediction[j],'was_correct':accur})
	return(predictions)



def handle_upload_2(fn):
	patients = []
	
	patients1 = []
	patients2 = []
	genes = []
	geneNames = []
	#data = {}
	data1 = {}
	data2 = {}
	group1 = []
	group2 = []
	group_labels1 = []
	group_labels2 = []
	group1_data = []
	group2_data = []
	#patient_ids = ['3dbe99d1-e3b8-4ee2-b6a8-2e2e12c6fbe9','6F4C8D30-47FB-47DF-9EB7-4E5881E3711E','95CEF916-5545-455B-920C-773A54FC7676','67C73260-A242-4BBA-87C5-D2302556DFF7','55262FCB-1B01-4480-	B322-36570430C917','3dbe99d1-e3b8-4ee2-b6a8-2e2e12c6fbe9','6F4C8D30-47FB-47DF-9EB7-4E5881E3711E','95CEF916-5545-455B-920C-773A54FC7676','67C73260-A242-4BBA-87C5-D2302556DFF7','55262FCB-1B01-4480-	B322-36570430C917']
	#patientfilename = 'nationwidechildrens.org_clinical_patient_brca.txt'
	
	patients.append([1,2,3,4])
	patients1.append(['1','2'])
	#patients.append(['3','4'])
	patients2.append(['3','4'])
	group_labels1.append([1,1])
	group_labels2.append([2,2])
	logstring = "Creating Plots for given input files... \n\n"
	logstring = "Reading gene expression data... \n"
	line_no = 0
	patient_ids = []
	survival = []
	survival_yrs = []
	data = []
	#group2.append(group_labels2)
	group1.append([1,1])
	group2.append([2,2])
	genes = ()
	dataset = Dataset(name='hsapiens_gene_ensembl',
                  host='http://www.ensembl.org')
	conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
	genes_3 = {}
	logstring = logstring + "\n\n Matching gene and protein IDs... \n"
	
	red_patch = mpatches.Patch(color='red', label='Group1')
	blue_patch = mpatches.Patch(color='blue', label='Group2')
	#lut = dict(zip(set(endData[0]), sns.hls_palette(len(set(endData[0])), l=0.5, s=0.8)))
	#col_colors = pd.DataFrame(endData[0])[0].map(lut)
	#print(col_colors)
	#colors = np.array(['#BB0000','#BB0000','#0000BB','#0000BB'])
	#df9 = pd.DataFrame(data=endData[1:,0:],index=geneNames,columns=patients)
	#df2 = pd.DataFrame(data=endData[0,0:], index='',columns=patients)
	#my_palette = dict(zip(df[.unique(), ["orange","yellow","brown"]))
	#row_colors = df2.cyl.map(my_palette)
	#fig, (ax1, ax2) = plt.subplots(1,2,sharex=True,sharey=True)
	colordict={0:'#BB0000',1:'#0000BB'}
	#logstring = logstring + str(df9)
	df2 = pd.read_csv(fn,delim_whitespace=True,header=None,index_col=0)
	#print(df2.head())
	df = df2.transpose()
	survival_real = df['SURVIVAL']
	#df['is_train'] = np.random.uniform(0, 1, len(df)) <= .75
	#train, test = df[df['is_train']==True], df[df['is_train']==False]
	features = df.columns[3:60483]
	x = df.iloc[:,5:500]
	y = df.iloc[:,2]
	print(df)
	rm = linear_model.LinearRegression()
	rm.fit(x,y)
	#print(rm.intercept_)
	#print(rm.coef_)
	#print(rm.predict(x))
	predictions = rm.predict(x)
	real_values = df.iloc[:,2].values.tolist()
	ret = []
	for j in range(1,len(survival_real)):
		accur = "FALSE"
		if(abs(float(predictions[j]) - float(real_values[j])) < 1.0):
			print("Foo")
			accur = "TRUE"
		ret.append({'patient_id':j,'real_value':real_values[j],'prediction':predictions[j],'was_correct':accur})
	return(ret)




def handle_upload_3(fn):
	patients = []
	
	patients1 = []
	patients2 = []
	genes = []
	geneNames = []
	#data = {}
	data1 = {}
	data2 = {}
	group1 = []
	group2 = []
	group_labels1 = []
	group_labels2 = []
	group1_data = []
	group2_data = []
	#patient_ids = ['3dbe99d1-e3b8-4ee2-b6a8-2e2e12c6fbe9','6F4C8D30-47FB-47DF-9EB7-4E5881E3711E','95CEF916-5545-455B-920C-773A54FC7676','67C73260-A242-4BBA-87C5-D2302556DFF7','55262FCB-1B01-4480-	B322-36570430C917','3dbe99d1-e3b8-4ee2-b6a8-2e2e12c6fbe9','6F4C8D30-47FB-47DF-9EB7-4E5881E3711E','95CEF916-5545-455B-920C-773A54FC7676','67C73260-A242-4BBA-87C5-D2302556DFF7','55262FCB-1B01-4480-	B322-36570430C917']
	#patientfilename = 'nationwidechildrens.org_clinical_patient_brca.txt'
	
	patients.append([1,2,3,4])
	patients1.append(['1','2'])
	#patients.append(['3','4'])
	patients2.append(['3','4'])
	group_labels1.append([1,1])
	group_labels2.append([2,2])
	logstring = "Creating Plots for given input files... \n\n"
	logstring = "Reading gene expression data... \n"
	line_no = 0
	patient_ids = []
	survival = []
	survival_yrs = []
	data = []
	#group2.append(group_labels2)
	group1.append([1,1])
	group2.append([2,2])
	genes = ()
	red_patch = mpatches.Patch(color='red', label='Group1')
	blue_patch = mpatches.Patch(color='blue', label='Group2')
	#lut = dict(zip(set(endData[0]), sns.hls_palette(len(set(endData[0])), l=0.5, s=0.8)))
	#col_colors = pd.DataFrame(endData[0])[0].map(lut)
	#print(col_colors)
	#colors = np.array(['#BB0000','#BB0000','#0000BB','#0000BB'])
	#df9 = pd.DataFrame(data=endData[1:,0:],index=geneNames,columns=patients)
	#df2 = pd.DataFrame(data=endData[0,0:], index='',columns=patients)
	#my_palette = dict(zip(df[.unique(), ["orange","yellow","brown"]))
	#row_colors = df2.cyl.map(my_palette)
	#fig, (ax1, ax2) = plt.subplots(1,2,sharex=True,sharey=True)
	colordict={0:'#BB0000',1:'#0000BB'}
	#logstring = logstring + str(df9)
	df2 = pd.read_csv(fn,delim_whitespace=True,header=None,index_col=0)
	#print(df2.head())
	df = df2.transpose()
	survival_real = df['SURVIVAL']
	#df['is_train'] = np.random.uniform(0, 1, len(df)) <= .75
	#train, test = df[df['is_train']==True], df[df['is_train']==False]
	features = df.columns[3:60483]
	x = df.iloc[:,5:500]
	y = df.iloc[:,2]
	print(df)
	rm = linear_model.LinearRegression()
	rm.fit(x,y)
	#print(rm.intercept_)
	#print(rm.coef_)
	#print(rm.predict(x))
	predictions = rm.predict(x)
	real_values = df.iloc[:,2].values.tolist()
	ret = []
	for j in range(1,len(survival_real)):
		accur = "FALSE"
		if(abs(float(predictions[j]) - float(real_values[j])) < 1.0):
			print("Foo")
			accur = "TRUE"
		ret.append({'patient_id':j,'real_value':real_values[j],'prediction':predictions[j],'was_correct':accur})
	return(ret)








