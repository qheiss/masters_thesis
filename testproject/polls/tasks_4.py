import string
from io import StringIO
from multiprocessing import Pool
import time
import pandas as pd
import numpy as np
import itertools
import networkx as nx
import matplotlib.pyplot as plt
from networkx.algorithms import bipartite
import scipy.sparse as sparse
from scipy.sparse import csr_matrix
from IPython.display import Audio, display
from sklearn.cluster.bicluster import SpectralCoclustering
from collections import Counter
import collections
from sklearn import preprocessing
flatten = lambda l: [item for sublist in l for item in sublist]
from scipy import stats
import scipy.spatial.distance as ssd
from scipy.cluster import hierarchy
import seaborn as sns; sns.set(color_codes=True)
from multiprocessing import Pool
from numpy import linalg as LA
from sklearn.cluster import KMeans
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import polls.weighted_aco_lib as lib
from shutil import copyfile
import plotly
import plotly.plotly as py
import plotly.graph_objs as go
import plotly.io as pio
import plotly.offline
#from .models import Choice, Question
from datetime import datetime
from networkx.readwrite import json_graph
import json
import imp
import seaborn as sns; sns.set(color_codes=True)
import mygene
from django.contrib.auth.models import User
from django.utils.crypto import get_random_string
import os
from celery import shared_task
#import datetime
import numpy as np
import matplotlib.pyplot as plt
import mpld3

import math

import gseapy as gp

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

from ndex2.nice_cx_network import NiceCXNetwork
import ndex2.client as nc
import ndex2
#from PyEntrezId import Conversion

#from polls.script3 import algo_output_ownfile_2

@shared_task
def create_random_user_accounts(total):
    for i in range(total):
        username = 'user_{}'.format(get_random_string(10, string.ascii_letters))
        email = '{}@example.com'.format(username)
        password = get_random_string(50)
        User.objects.create_user(username=username, email=email, password=password)
    return '{} random users created with success!'.format(total)

#@shared_task(name="run_analysis")
#def run_analysis(s,L_g_min,L_g_max,fh,prot_fh,nbr_iter,nbr_ants,evap):
#	return(algo_output_ownfile_2(s,L_g_min,L_g_max,fh,prot_fh,nbr_iter,nbr_ants,evap))

def jac(x,y):
    if len(x)>0 and len(y)>0:
        return len(set(x).intersection(set(y)))/len((set(x).union(set(y))))
    else:
        return(0)
    
def jac_matrix(true,pred):
    res = np.zeros((len(true),len(true)))
    for i in range(len(true)):
        for j in range(len(true)):
            res[i,j] = jac(true[i],pred[j])
    cand1 = (res[0][0],res[1][1])
    cand2 = (res[0][1],res[1][0])
    if sum(cand1)>sum(cand2):
        return(cand1)
    else:
        return(cand2)
    
def matches(true1,true2,pred1,pred2):
    cand1 = (round(len(set(pred1).intersection(set(true2)))*100/len(true2)),
             round(len(set(pred2).intersection(set(true1)))*100/len(true1)))
    cand2 = (round(len(set(pred1).intersection(set(true1)))*100/len(true1)),
             round(len(set(pred2).intersection(set(true2)))*100/len(true2)))
    if sum(cand1)>sum(cand2):
        ans = cand1
    else:
        ans = cand2
    print(str(ans[0])+"%                  " + str(ans[1])+"%")
    


def joined_net(B,G):
    #joined net
    A_b = nx.adjacency_matrix(B).todense()
    A_b = A_b *1 ## trick to switch from boolean
    A_g = nx.adjacency_matrix(G).todense()
    n = len(A_g)
    A = A_b
    A[:n,:n] = A_g  
    return(A)





def hi(A_j,n,m):
    H = np.zeros((n+m,n+m))
    P = LA.matrix_power(A_j, 2)
    for node1 in range(n+m):
            for node2 in range(node1+1,n+m):
                if P[node1,node2] >0:
                    node_intersec = np.sum(np.multiply(A_j[node1,:],A_j[node2,:]))+A_j[node1,node2]
                    deg_node1 = np.sum(A_j[node1,:])
                    deg_node2 = np.sum(A_j[node2,:])
                    if node1<n and node2<n: #both beling to X
                        H[node1,node2] = node_intersec/(deg_node1+deg_node2+2)
                    if node1>=n and node2>=n: #both belong to Y
                        H[node1,node2] = 0.5*node_intersec/(deg_node1+deg_node2)
                    if (node2>=n) and (node1 <n): #node1 in X, node2 in Y or the other way around
                        H[node1,node2] = 4*node_intersec/(deg_node1+deg_node2)
                    H[node2,node1] = H[node1,node2]
    return(H*10)

@shared_task(name="make_empty_figure")
def make_empty_figure():
	fig = plt.figure(figsize=(10,8))
	plt.savefig("polls/static/progress.png")
	plt.close(fig)

@shared_task(name="add_loading_image")
def add_loading_image():
	copyfile("polls/static/loading.gif","polls/static/loading_1.gif")

@shared_task(name="remove_loading_image")
def remove_loading_image():
	os.unlink("polls/static/loading_1.gif")

@shared_task(name="write_to_file_1")
def write_to_file_1(text,pth):
    text_file = open(pth, "w")
    text_file.write(text)
    text_file.close()
  

@shared_task(name="write_metadata_to_file")
def write_metadata_to_file(metadata):
    text_file = open("polls/static/metadata.txt", "w")
    text_file.write("<table style=\"font-family: arial, sans-serif;border-collapse: collapse;width: 100%;\"><tr><th style=\"border: 1px solid #dddddd;text-align: left;padding: 8px;\">Patient Group</th><th style=\"border: 1px solid #dddddd;text-align: left;padding: 8px;\">Mean Survival</th><th style=\"border: 1px solid #dddddd;text-align: left;padding: 8px;\">mean age at diagnosis</th><th style=\"border: 1px solid #dddddd;text-align: left;padding: 8px;\">mean tumor size</th></tr>")
    for elem in metadata:
    	text_file.write("<tr>")
    	text_file.write("<th style=\"border: 1px solid #dddddd;text-align: left;padding: 8px;\">" + str(elem['group']) + "</th>")
    	text_file.write("<th style=\"border: 1px solid #dddddd;text-align: left;padding: 8px;\">" + str(elem['survival']) + "</th>")
    	text_file.write("<th style=\"border: 1px solid #dddddd;text-align: left;padding: 8px;\">" + str(elem['age']) + "</th>")
    	text_file.write("<th style=\"border: 1px solid #dddddd;text-align: left;padding: 8px;\">" + str(elem['size']) + "</th>")
    	text_file.write("</tr>")
    text_file.write("</table>")
    text_file.close()
@shared_task(name="empty_log_file")
def empty_log_file():
    text_file = open("/home/quirin/testproject/polls/static/output_console.txt", "w")
    text_file.write("")
    text_file.close()


@shared_task(name="metadata_to_string")
def metadata_to_string(metadata):
    ret = "<table><tr><th>Patient Group</th><th>Correlation with clinical Group</th><th>Mean Survival</th><th>percentage of metastasis</th><th>mean age at diagnosis</th></tr>"
    for elem in metadata:
    	ret = ret + "<tr>"
    	ret = ret + "<th>" + str(elem['group']) + "</th>"
    	ret = ret + "<th>" + str(elem['survival']) + "</th>"
    	ret = ret + "<th>" + str(elem['age']) + "</th>"
    	ret = ret + "<th>" + str(elem['size']) + "</th>"
    	ret = ret + "</tr>"
    ret = ret + "</table>"

@shared_task(name="list_metadata")
def list_metadata():
    fh1 = open("polls/static/metadata.txt")
    lines = fh1.read()
    lines = lines.replace('<table><tr><th>','')
    lines = lines.replace('<tr><th>','')
    lines = lines.replace('</table>','')
    lines = lines.replace('</th></tr>','\n')
    lines = lines.replace('</th><th>','\t')
    line1 = lines.split('\n')[1].split('\t')
    line2 = lines.split('\n')[2].split('\t')
    print(line1)
    ret = []
    ret.append({'group':"Group 1",'lm':line1[0],'bm':line1[1],'lymph':line1[2],'metastasis':line1[3],'path_er':line1[4],'path_pr':line1[5],'prognosis_good':line1[6]})
    ret.append({'group':"Group 2",'lm':line2[0],'bm':line2[1],'lymph':line2[2],'metastasis':line2[3],'path_er':line2[4],'path_pr':line2[5],'prognosis_good':line2[6]})
    return(ret)

@shared_task(name="list_metadata_3")
def list_metadata_3():
    fh1 = open("polls/static/metadata.txt")
    lines = fh1.read()
    lines = lines.replace('<table><tr><th>','')
    lines = lines.replace('<tr><th>','')
    lines = lines.replace('</th></table>','')
    lines = lines.replace('</table>','')
    lines = lines.replace('</th></tr>','\n')
    lines = lines.replace('</th><th>','\t')
    emptydict = {}
    print(len(lines.split('\n')))
    if(len(lines.split('\n')) < 3):
    	print("basdbasdbasdb")
    	return(emptydict,emptydict,emptydict)
    line0 = lines.split('\n')[0].split('\t')
    line1 = lines.split('\n')[1].split('\t')
    line2 = lines.split('\n')[2].split('\t')
    print(line0)
    print(line1)
    print(line2)
    ret = []
    dict3 = {}
    dict1 = {}
    dict2 = {}
    dict0 = {}
    list0 = []
    list1 = []
    list2 = []
    ctr = 0
    dict3['params'] = line0
    dict3['gr1'] = line1
    dict3['gr2'] = line2
    dict3['all'] = zip(dict3['params'],dict3['gr1'],dict3['gr2'])
    #for elem in line0:
    #	print(ctr)
    #	dict0[ctr] = line0[ctr]
    #	dict1[elem] = line1[ctr]
    #	dict2[elem] = line2[ctr]
    #	ctr = ctr + 1
    for i in range(0,len(line0)-1):
    	#print(ctr)
    	print(line0[i])
    	dict0[i] = line0[i]
    	dict1[dict0[i]] = line1[i]
    	dict2[dict0[i]] = line2[i]
    	ctr = ctr + 1
    #ret.append(dict1)
    print(dict0)
    print(dict1)
    print(dict2)
    #ret.append(dict1)
    #ret.append(dict2)
    #ret.append({'group':"Group 1",line0[0]:line1[0],'bm':line1[1],'lymph':line1[2],'metastasis':line1[3],'path_er':line1[4],'path_pr':line1[5],'prognosis_good':line1[6]})
    #ret.append({'group':"Group 2",'lm':line2[0],'bm':line2[1],'lymph':line2[2],'metastasis':line2[3],'path_er':line2[4],'path_pr':line2[5],'prognosis_good':line2[6]})
    #return(dict3['params'],dict3['gr1'],dict3['gr2'])    
    return(dict0,dict1,dict2)
@shared_task(name="list_metadata_2")
def list_metadata_2():
    fh1 = open("polls/static/metadata.txt")
    lines = fh1.read()
    lines = lines.replace('<table><tr><th>','')
    lines = lines.replace('<tr><th>','')
    lines = lines.replace('</table>','')
    lines = lines.replace('</th></tr>','\n')
    lines = lines.replace('</th><th>','\t')
    line1 = lines.split('\n')[1].split('\t')
    line2 = lines.split('\n')[2].split('\t')
    ret = []
    if(len(line1) < 5):
    	is_lungc = "true"
    	ret.append({'group':"Group 1",'gender':line1[0],'relapse':line1[1],'living':line1[2]})
    	ret.append({'group':"Group 2",'gender':line2[0],'relapse':line2[1],'living':line2[2]})
    	print(line1)
    else:
    	ret.append({'group':"Group 1",'lm':line1[0],'bm':line1[1],'lymph':line1[2],'metastasis':line1[3],'path_er':line1[4],'path_pr':line1[5],'prognosis_good':line1[6]})
    	ret.append({'group':"Group 2",'lm':line2[0],'bm':line2[1],'lymph':line2[2],'metastasis':line2[3],'path_er':line2[4],'path_pr':line2[5],'prognosis_good':line2[6]})
    	is_lungc = "false"
    return(ret,is_lungc)


@shared_task(name="list_metadata_4")
def list_metadata_4(path):
    fh1 = open(path)
    lines = fh1.read()
    lines = lines.replace('<table><tr><th>','')
    lines = lines.replace('<tr><th>','')
    lines = lines.replace('</th></table>','')
    lines = lines.replace('</table>','')
    lines = lines.replace('</th></tr>','\n')
    lines = lines.replace('</th><th>','\t')
    emptydict = {}
    print(len(lines.split('\n')))
    if(len(lines.split('\n')) < 3):
    	print("basdbasdbasdb")
    	return(emptydict,emptydict,emptydict)
    line0 = lines.split('\n')[0].split('\t')
    line1 = lines.split('\n')[1].split('\t')
    line2 = lines.split('\n')[2].split('\t')
    print(line0)
    print(line1)
    print(line2)
    ret = []
    dict3 = {}
    dict1 = {}
    dict2 = {}
    dict0 = {}
    list0 = []
    list1 = []
    list2 = []
    ctr = 0
    dict3['params'] = line0
    dict3['gr1'] = line1
    dict3['gr2'] = line2
    dict3['all'] = zip(dict3['params'],dict3['gr1'],dict3['gr2'])
    #for elem in line0:
    #	print(ctr)
    #	dict0[ctr] = line0[ctr]
    #	dict1[elem] = line1[ctr]
    #	dict2[elem] = line2[ctr]
    #	ctr = ctr + 1
    for i in range(0,len(line0)-1):
    	#print(ctr)
    	print(line0[i])
    	dict0[i] = line0[i]
    	dict1[dict0[i]] = line1[i]
    	dict2[dict0[i]] = line2[i]
    	ctr = ctr + 1
    #ret.append(dict1)
    print(dict0)
    print(dict1)
    print(dict2)
    #ret.append(dict1)
    #ret.append(dict2)
    #ret.append({'group':"Group 1",line0[0]:line1[0],'bm':line1[1],'lymph':line1[2],'metastasis':line1[3],'path_er':line1[4],'path_pr':line1[5],'prognosis_good':line1[6]})
    #ret.append({'group':"Group 2",'lm':line2[0],'bm':line2[1],'lymph':line2[2],'metastasis':line2[3],'path_er':line2[4],'path_pr':line2[5],'prognosis_good':line2[6]})
    #return(dict3['params'],dict3['gr1'],dict3['gr2'])    
    return(dict0,dict1,dict2)
@shared_task(name="list_metadata_2")
def list_metadata_2():
    fh1 = open("polls/static/metadata.txt")
    lines = fh1.read()
    lines = lines.replace('<table><tr><th>','')
    lines = lines.replace('<tr><th>','')
    lines = lines.replace('</table>','')
    lines = lines.replace('</th></tr>','\n')
    lines = lines.replace('</th><th>','\t')
    line1 = lines.split('\n')[1].split('\t')
    line2 = lines.split('\n')[2].split('\t')
    ret = []
    if(len(line1) < 5):
    	is_lungc = "true"
    	ret.append({'group':"Group 1",'gender':line1[0],'relapse':line1[1],'living':line1[2]})
    	ret.append({'group':"Group 2",'gender':line2[0],'relapse':line2[1],'living':line2[2]})
    	print(line1)
    else:
    	ret.append({'group':"Group 1",'lm':line1[0],'bm':line1[1],'lymph':line1[2],'metastasis':line1[3],'path_er':line1[4],'path_pr':line1[5],'prognosis_good':line1[6]})
    	ret.append({'group':"Group 2",'lm':line2[0],'bm':line2[1],'lymph':line2[2],'metastasis':line2[3],'path_er':line2[4],'path_pr':line2[5],'prognosis_good':line2[6]})
    	is_lungc = "false"
    return(ret,is_lungc)
    
@shared_task(name="ants_2")
def ants_2(a,b,n,m,H,GE,G,clusters,cost_limit,K,evaporation,th,L_g_min,L_g_max,eps,times,opt= None,pts = False,show_pher = True,show_plot = True, print_runs = True, save = None, show_nets = True):
    ge = GE.values
    H =H.astype(np.short)
    N = neigborhood(H,n,th)
    patients = np.arange(n,n+m)  
    #text_file = open("/home/quirin/testproject/polls/static/output_console.txt", "w")
    #text_file.write("bababa")
    #text_file.close()
    cost = H/10
    cost = np.max(cost)-cost
    scores = []
    avs = []
    count_big = 0
    max_total_score = 0
    max_round_score = -100
    av_score = 0
    st = time.time()
    t0 = np.ones((n+m,n+m))*5
    t0 = t0.astype(np.short)
    probs= prob_upd(H,t0,a,b,n,th,N)
    end = time.time()
    flag = False
    score_change = []
    print ("Running time statistics:")
    print ("###############################################################")
    print("the joint graph has "+ str(n+m) + " nodes")
    print("probability update takes "+str(round(end-st,3)))
    W = 0
    while np.abs(max_round_score-av_score)>eps and count_big<times and (W<m/3):
        av_score = 0
        W = 0
        max_round_score = 0
        scores_per_round = []

        for i in range(K):
            #for each ant
            st = time.time()
            tot_score,gene_groups,patients_groups,new_scores,wars,no_int = ant_job(GE,N,H,th,clusters,probs,a,b,cost,m,n,patients,count_big,i,cost_limit,L_g_min,L_g_max,G,ge,print_runs)
            end = time.time()
            W = W+wars
            if count_big ==0 and i ==0:
                print("one ant run takes "+str(round(end-st,3)))
            scores_per_round.append(tot_score)
            av_score = av_score + tot_score
            if tot_score > max_round_score:
                max_round_score = tot_score
                solution = (gene_groups,patients_groups)
                full_scores = new_scores
                solution_big = (no_int,patients_groups)
            if count_big ==0 and i ==K-1:
                gs = 1.5*max_round_score

                t_max = (1/evaporation)*gs
                t_min = 0
   
                t0 = np.ones((n+m,n+m))*t_max
        #after all ants have finished:
        scores.append(scores_per_round)
        
        #saving rhe best overall solution
        if max_round_score>max_total_score:
            max_total_score = max_round_score
            best_solution = solution
            max_full_scores = full_scores 
            solution_big_best = solution_big
        score_change.append(round(max_round_score,3))
        print("Iteration # "+ str(count_big+1))
        print("best round score: " + str(round(max_round_score,3)))
        print("average score: " + str(round(av_score/K,3)))
        print("foobar")
        with open("/home/quirin/testproject/polls/static/output_console.txt", "w") as text_file:
        	#print("foobar")
        	text_file.write("Iteration # "+ str(count_big+1) + " completed. \n" + "Best round score: " + str(round(max_round_score,3)) + "\n" + "Average score: " + str(round(av_score/K,3)))
        	text_file.close()
        av_score = av_score/K
        avs.append(round(av_score,2))
        #print(scores)
        print(avs)
        #pher. and prob. updates
        t = pher_upd(t0,t_min,evaporation,max_full_scores,solution_big_best,flag)
        t0 = np.copy(t)
        
        probs= prob_upd(H,t,a,b,n,th,N)
        
        #visualization options:
        
        if show_pher:
            fig = plt.figure(figsize=(18,12))
            ax = fig.add_subplot(111)
            t_max = np.max(t)   
            cax = ax.matshow(t, interpolation='nearest',cmap=plt.cm.RdPu,vmin = t_min,vmax = t_max)
            plt.colorbar(cax)
            plt.title("Pheramones")
            plt.show(block=False)
            plt.close(fig)

        count_big = count_big +1
        if show_nets:
            features(solution, GE,G)    
        if show_plot:
            fig = plt.figure(figsize=(10,8))
            plt.boxplot(np.asarray(scores).T,manage_xticks = True, patch_artist=True)
            if opt!=None:
                plt.axhline(y=opt,label = "optimal solution score", c = "r")
            #plt.ylim((0,1))
            #plt.legend()
            #this was not commented before #plt.show(block=False)
            plt.savefig("/home/quirin/testproject/polls/static/progress.png")
            plt.close(fig)
        if len(set(score_change[:3])) ==1 and len(score_change)>3:
            flag = True
    if save != None:
        fig = plt.figure(figsize=(10,8))
        plt.boxplot(np.asarray(scores).T,manage_xticks = True)
        if opt!=None:
            plt.axhline(y=opt,label = "optimal solution score", c = "r")
        #plt.legend()
        plt.savefig(save+".png")
        plt.close(fig)
        
    #after the solutution is found we make sure to cluster patients the last time with that exact solution:
    data_new = ge[solution[0][0]+solution[0][1],:]
    kmeans = KMeans(n_clusters=2, random_state=0).fit(data_new.T)
    labels = kmeans.labels_
    patients_groups =[]
    for clust in range(clusters):
        wh = np.where(labels == clust)[0]
        group_p = [patients[i] for i in wh]
        patients_groups.append(group_p)
    if np.mean(ge[best_solution[0][0],:][:,(np.asarray(patients_groups[0])-n)])<np.mean(ge[best_solution[0][1],:][:,(np.asarray(patients_groups[0])-n)]):
        patients_groups = patients_groups[::-1]
    best_solution = [best_solution[0],patients_groups]
    
    print("best total score: "+str(max_total_score))
    #print_clusters(GE,best_solution)
    #features(best_solution, GE,G)
    return(best_solution,t,max_total_score,np.asarray(scores).T)
    
def neigborhood(H,n,th):
    N_per_patient = []
    dim = len(H)
    for i in range(n,dim):
        if th<0:
            N = np.where(H[i,:]>0.001)[0]
        else:
            rad = np.mean(H[i,:]) + th*np.std(H[i,:])
            N = np.where(H[i,:]>rad)[0]
        #N = np.where(H[i,:]>0)[0]
        N_per_patient.append(N)
    return N_per_patient
    
def prob_upd(H,t,a,b,n,th,N_per_patient):
    P_per_patient = []
    dim = len(H)
    temp_t = np.power(t,a)
    temp_H = np.power(H,b)
    temp = temp_t*temp_H 
    for i in range(n,dim):
        N_temp = N_per_patient[i-n]
        P = temp[:,N_temp]
        s = np.sum(P,axis = 1)
        s[s <1.e-4] = 1
        sum_p = 1/s
        sum_p = sum_p[:,None]
        P_new = P*sum_p[:np.newaxis]
        P_per_patient.append(P_new)

    return(P_per_patient)
        



def ant_job(GE,N,H,th,clusters,probs,a,b,cost,m,n,patients,count_big,count_small,cost_limit,L_g_min,L_g_max,G,ge,print_runs):
    if print_runs:
        print(str(count_big)+"."+str(count_small) + " run")
    paths = []
    wars = 0
    #set an ant on every patient
    
    for walk in range(m):
        k = cost_limit
        path = []
        start = patients[walk]
        Nn = N[walk] #neigbohood
        path.append(start)
        go = True
        P = probs[walk]            
        while go == True:
            P_new = P[start,:]
            #if there is any node inside the radious - keep mooving
            if np.sum(P_new)> 0.5:
                #transition:
                tr = np.random.choice(Nn,1,False,p = P_new)[0]
                c = cost[patients[walk],tr]
                #if there is any cost left we keep going
                if k-c >0:
                    path.append(tr)
                    start = tr
                    k = k - c
                #if not we are done and we save only genes from the path
                else:
                    go = False
            #no node to go - we are done and we save only genes from the path
            else:
                go = False
        path = np.asarray(path)
        #we are saving only genes
        path = path[path<n]
        paths.append(path)
        if len(path) == 0:
            wars = wars+1
            print("WARNING: emply path found")
    data_new = ge[list(set(flatten(paths))),:]
    kmeans = KMeans(n_clusters=2, random_state=0).fit(data_new.T)
    labels = kmeans.labels_
    gene_groups_set =[]
    patients_groups =[]
    for clust in range(clusters):
        wh = np.where(labels == clust)[0]
        group_g = [paths[i] for i in wh]
        group_g = flatten(group_g)
        gene_groups_set.append(set(group_g))
        #save only most common genes for a group
        group_p = [patients[i] for i in wh]
        patients_groups.append(group_p)
        
    #delete intersecting genes between groups
    
    I = set.intersection(*gene_groups_set)
    no_int =[list(gene_groups_set[i].difference(I)) for i in range(clusters)]
    gene_groups = no_int
    
    # make sure that gene clusters correspond to patients clusters:
    if np.mean(ge[gene_groups[0],:][:,(np.asarray(patients_groups[0])-n)])<np.mean(ge[gene_groups[1],:][:,(np.asarray(patients_groups[0])-n)]):
        patients_groups = patients_groups[::-1]
     
    gene_groups,sizes= clean_net(gene_groups,patients_groups, clusters,L_g_min,G,GE)
    

        
    new_scores = score(G,patients_groups,gene_groups,n,m,ge,sizes,L_g_min,L_g_max)
    tot_score = new_scores[0][0]*new_scores[0][1]+new_scores[1][0]*new_scores[1][1]   
    return(tot_score,gene_groups,patients_groups,new_scores,wars,no_int)
    
def pher_upd(t,t_min,p,scores,solution,flag):
    t = t*(1-p)
    t_new = np.copy(t)
    score = scores[0][0]*scores[0][1]+scores[1][0]*scores[1][1]
    for i in range(len(solution[0])):
        group_g = solution[0][i]
        group_p = solution[1][i]
        #score = scores[i][0]*scores[i][1]
        #ge_score = new_scores[i][0]*10
        #ppi_score = new_scores[i][1]*10
        for g1 in group_g:
            for p1 in group_p:
                t_new[g1,p1] = t[g1,p1]+ score
                t_new[p1,g1] = t[p1,g1]+ score
            for g2 in group_g:
                t_new[g1,g2] = t[g1,g2]+ score


    t_new[t_new < t_min] = t_min
            
    
    return(t_new)

def h_upd(H,scores,solution,p):
    H_new = np.copy(H)
    H_new = H*(1-p)
    score = (scores[0][0]*scores[0][1]+scores[1][0]*scores[1][1])*2
    for i in range(len(solution[0])):
        group_g = solution[0][i]
        group_p = solution[1][i]
        #score = scores[i][0]*scores[i][1]
        #ge_score = new_scores[i][0]*10
        #ppi_score = new_scores[i][1]*10
        for g1 in group_g:
            for p1 in group_p:
                H_new[g1,p1] = H[g1,p1]+ score
                H_new[p1,g1] = H[p1,g1]+ score
            for g2 in group_g:
                H_new[g1,g2] = H[g1,g2]+ score


    return(H_new)
    
    
    
    
def score(G,patients_groups,gene_groups,n,m,ge,sizes,L_g_min,L_g_max):
    clusters = len(patients_groups)
    conf_matrix = np.zeros((clusters,clusters))
    conect_ppi = []
    for i in range(clusters): #over genes
        group_g = np.asarray(gene_groups[i])
        s = sizes[i]
        if len(group_g)>0:
            for j in range(clusters): #over patients
                group_p = np.asarray(patients_groups[j])
                if len(group_p)>0:
                # gene epression inside the group
                    conf_matrix[i,j] = np.mean(ge[group_g,:][:,(group_p-n)])
            #ppi score    
            con_ppi = 1
            if s<L_g_min:
                con_ppi = s/L_g_min
            elif s>L_g_max:
                con_ppi = L_g_max/s
            conect_ppi.append(con_ppi)
        else:
            conect_ppi.append(0)           
    ans = []
    for i in range(clusters):
        all_ge = np.sum(conf_matrix[i,:])
        in_group = conf_matrix[i,i]
        out_group = all_ge - in_group
        ge_con = in_group-out_group
        #scaled = scaleBetween(num,0,0.5,0,1)
        ans.append((ge_con ,conect_ppi[i]))
        
    return(ans)
    

def aco_preprocessing(path_expr, path_ppi, col,log2, gene_list = None, size = None, sample= None):
    # path_expr - path for gene expression
    # path_ppi - path for ppi
    # col - split variable name (ONLY TWO CLASSES)
    # log2 - log2 transform
    #gene_list - preselected genes (if any)
    #size -  if genes are not preselected specify size of the gene set  for standard deviation selection
    # sample = None - all patients, otherwise specify fraction of patients taken
    expr = pd.read_csv(path_expr,sep = "\t") 
    expr = expr.set_index("Unnamed: 0")
    val1,val2 = list(set(expr[col]))
    group1_true = list(expr[expr[col]==val1].index)
    group2_true = list(expr[expr[col]==val2].index)
    patients_new = group1_true+group2_true
    if sample!=None:
        idx = list(expr.index)
        new_idx = np.random.choice(idx,int(sample*len(idx)),False)
        expr = expr.loc[new_idx]
        group1_true = list(expr[expr[col]==val1].index)
        group2_true = list(expr[expr[col]==val2].index)
        patients_new = group1_true+group2_true

    expr = expr.loc[patients_new]    
    net = pd.read_csv(path_ppi,sep = "\t", header= None)
    nodes_ppi = set(net[0]).union(set(net[1]))
    genes_ge = list(set(expr.columns) - set([col]))
    new_genes = [int(x) for x in genes_ge]
    intersec_genes = set.intersection(set(new_genes), set(nodes_ppi))
    genes_for_expr = [str(x) for x in list(intersec_genes)]
    expr = expr[genes_for_expr]
    #20188 genes
    if log2:
        expr = np.log2(expr)
    z_scores = stats.zscore(expr) 
    z_scores = pd.DataFrame(z_scores,columns = expr.columns, index = expr.index)
    if gene_list !=None and size == None:# gene list is given
        new_genes = [str(gene) for gene in gene_list] 
        
    elif gene_list == None and size!= None: #std selection
        std_genes = expr[genes_for_expr].std()
        std_genes, genes_for_expr = zip(*sorted(zip(std_genes, genes_for_expr)))
        genes_for_expr = genes_for_expr[len(std_genes)-size:]
        new_genes = list(genes_for_expr)
    elif gene_list == None and size == None: #all genes
        new_genes = genes_for_expr
    else:
        print("please specify gene selection method: predifined list, standart deviation filtering or none of them")
        return()

    expr = expr[new_genes]
    z_scores = z_scores[new_genes].values
    
    labels_B = dict()
    rev_labels_B = dict()
    node = 0
    #nodes = set(deg_nodes + genes_aco)
    for g in new_genes:
       labels_B[node] = g
       rev_labels_B[g] = node
       node = node+1
    for p in patients_new:
       labels_B[node] = p
       rev_labels_B[p] = node
       node = node+1
    

    #scaler = preprocessing.MinMaxScaler(feature_range=(0, 1))
    #sim = scaler.fit_transform(expr)
    data_aco = pd.DataFrame(z_scores,columns= new_genes, index= patients_new)
    data_aco = data_aco.T
    n,m = data_aco.shape
    
    GE = pd.DataFrame(data_aco.values,index = np.arange(n), columns=np.arange(n,n+m))
    t = 2
    b = np.matrix(data_aco>t)
    b_sp = csr_matrix(b)
    B = bipartite.from_biadjacency_matrix(b_sp)
    
    
    G = nx.Graph()
    G.add_nodes_from(np.arange(n))
    for row in net.itertuples():
        node1 = str(row[1])
        node2 = str(row[2])
        if node1 in set(new_genes) and node2 in set(new_genes):    
            G.add_edge(rev_labels_B[node1],rev_labels_B[node2])
    A_new= nx.adj_matrix(G).todense()


def aco_preprocessing_ownfile(fh, fh_ppi, col,log2, gene_list = None, size = None, sample= None):
    # path_expr - path for gene expression
    # path_ppi - path for ppi
    # col - split variable name (ONLY TWO CLASSES)
    # log2 - log2 transform
    #gene_list - preselected genes (if any)
    #size -  if genes are not preselected specify size of the gene set  for standard deviation selection
    # sample = None - all patients, otherwise specify fraction of patients taken
    expr = pd.read_csv(fh,sep = "\t") 
    expr = expr.set_index("Unnamed: 0")
	#TODO: check if column 'prognosis' or 'cancer type' exists, set column based on this info
    if('cancer_type' in list(expr)):
    	col = 'cancer_type'
    else:
    	col = 'prognosis'
    val1,val2 = list(set(expr[col]))
    group1_true = list(expr[expr[col]==val1].index)
    group2_true = list(expr[expr[col]==val2].index)
    patients_new = group1_true+group2_true
    if sample!=None:
        idx = list(expr.index)
        new_idx = np.random.choice(idx,int(sample*len(idx)),False)
        expr = expr.loc[new_idx]
        group1_true = list(expr[expr[col]==val1].index)
        group2_true = list(expr[expr[col]==val2].index)
        patients_new = group1_true+group2_true

    expr = expr.loc[patients_new]    
    net = pd.read_csv(fh_ppi,sep = "\t", header= None)
    nodes_ppi = set(net[0]).union(set(net[1]))
    genes_ge = list(set(expr.columns) - set([col]))
    new_genes = [int(x) for x in genes_ge]
    intersec_genes = set.intersection(set(new_genes), set(nodes_ppi))
    genes_for_expr = [str(x) for x in list(intersec_genes)]
    expr = expr[genes_for_expr]
    #20188 genes
    if log2:
        expr = np.log2(expr)
    z_scores = stats.zscore(expr) 
    z_scores = pd.DataFrame(z_scores,columns = expr.columns, index = expr.index)
    if gene_list !=None and size == None:# gene list is given
        new_genes = [str(gene) for gene in gene_list] 
        
    elif gene_list == None and size!= None: #std selection
        std_genes = expr[genes_for_expr].std()
        std_genes, genes_for_expr = zip(*sorted(zip(std_genes, genes_for_expr)))
        genes_for_expr = genes_for_expr[len(std_genes)-size:]
        new_genes = list(genes_for_expr)
    elif gene_list == None and size == None: #all genes
        new_genes = genes_for_expr
    else:
        print("please specify gene selection method: predifined list, standart deviation filtering or none of them")
        return()

    expr = expr[new_genes]
    z_scores = z_scores[new_genes].values
    
    labels_B = dict()
    rev_labels_B = dict()
    node = 0
    #nodes = set(deg_nodes + genes_aco)
    for g in new_genes:
       labels_B[node] = g
       rev_labels_B[g] = node
       node = node+1
    for p in patients_new:
       labels_B[node] = p
       rev_labels_B[p] = node
       node = node+1
    

    #scaler = preprocessing.MinMaxScaler(feature_range=(0, 1))
    #sim = scaler.fit_transform(expr)
    data_aco = pd.DataFrame(z_scores,columns= new_genes, index= patients_new)
    data_aco = data_aco.T
    n,m = data_aco.shape
    
    GE = pd.DataFrame(data_aco.values,index = np.arange(n), columns=np.arange(n,n+m))
    t = 2
    b = np.matrix(data_aco>t)
    b_sp = csr_matrix(b)
    B = bipartite.from_biadjacency_matrix(b_sp)
    
    
    G = nx.Graph()
    G.add_nodes_from(np.arange(n))
    for row in net.itertuples():
        node1 = str(row[1])
        node2 = str(row[2])
        if node1 in set(new_genes) and node2 in set(new_genes):    
            G.add_edge(rev_labels_B[node1],rev_labels_B[node2])
    A_new= nx.adj_matrix(G).todense()

    H = HI_big(data_aco, gtg_weight = 1, gtp_weight=1 ,ptp_weight = 1)
    
    group1_true_ids= [rev_labels_B[x] for x in group1_true]
    group2_true_ids= [rev_labels_B[x] for x in group2_true]
    
    return B,G,H,n,m,GE,A_new,group1_true_ids,group2_true_ids,labels_B,rev_labels_B,val1,val2


def HI_big(data_aco, gtg_weight = 1, gtp_weight=1 ,ptp_weight = 1):
    scaler = preprocessing.MinMaxScaler(feature_range=(0, 1))#
    H_g_to_g = (data_aco.T.corr())*gtg_weight
    H_p_to_p = data_aco.corr()*ptp_weight
    H_g_to_g = scaler.fit_transform(H_g_to_g)
    H_p_to_p = scaler.fit_transform(H_p_to_p)
    H_g_to_p = scaler.fit_transform(data_aco)
    H_full_up = np.concatenate([H_g_to_g,H_g_to_p*gtp_weight], axis = 1)
    H_full_down = np.concatenate([H_g_to_p.T*gtp_weight,H_p_to_p], axis = 1)
    H_full =  np.concatenate([H_full_up,H_full_down], axis = 0)*10
#    H_full[H_full < 1] = 1
#    np.fill_diagonal(H_full, 1)
    np.fill_diagonal(H_full, 0)
    return(H_full)
    

def most_common(lst,top,L_g):
    data = Counter(lst)
    l = len(data)
    take = int(top*l)
    if take == 0:
        take = 1
    if take <L_g:
        take = L_g
    count = data.most_common(take)
    genes = [x[0] for x in count]
    return genes



def scaleBetween(unscaledNum, minAllowed, maxAllowed, min_cur, max_cur):
  return (maxAllowed - minAllowed) * (unscaledNum - min_cur) / (max_cur
         - min_cur) + minAllowed
def print_clusters(GE,solution):
    grouping_p = []
    p_num = list(GE.columns)
    for p in p_num:
        if p in solution[1][0]:
            grouping_p.append(1)
        else:
            grouping_p.append(2)
    grouping_p = pd.DataFrame(grouping_p,index = p_num)
    grouping_g = []
    g_num = list(GE.index)
    for g in g_num:
        if g in solution[0][0]:
            grouping_g.append(1)
        elif  g in solution[0][1]:
            grouping_g.append(2)
        else:
            grouping_g.append(3)
            
    grouping_g = pd.DataFrame(grouping_g,index = g_num)
    species = grouping_p[0]
    lut = {1: '#A52A2A', 2: '#7FFFD4'}
    row_colors = species.map(lut)
    species = grouping_g[0]
    lut = {1: '#A52A2A', 2: '#7FFFD4', 3:'#FAEBD7'}
    col_colors = species.map(lut)
    sns.clustermap(GE.T, row_colors=row_colors, col_colors = col_colors,figsize=(15, 10))
    
def features(solution, GE,G,pos = None):
    genes1,genes2 = solution[0]
    patients1, patients2 = solution[1]
    
    means1 = list(np.mean(GE[patients1].loc[genes1],axis = 1)-np.mean(GE[patients2].loc[genes1],axis = 1).values)
    means2 = list(np.mean(GE[patients1].loc[genes2],axis = 1)-np.mean(GE[patients2].loc[genes2],axis = 1).values)
    G_small = nx.subgraph(G,genes1+genes2)
    
    fig = plt.figure(figsize=(15,10))
    vmin = min(means1+means2)
    vmax = max(means1+means2)
    if pos == None:
        pos = nx.spring_layout(G_small)
    ec = nx.draw_networkx_edges(G_small,pos)
    nc1 = nx.draw_networkx_nodes(G_small,nodelist =genes1, pos = pos,node_color=means1, node_size=200,alpha=1.0,
                                 vmin=vmin, vmax=vmax,node_shape = "^",cmap =plt.cm.viridis)
    nc2 = nx.draw_networkx_nodes(G_small,nodelist =genes2, pos = pos,node_color=means2, node_size=200,
                                 alpha=1.0,
                                 vmin=vmin, vmax=vmax,node_shape = "o",cmap =plt.cm.viridis)
    nx.draw_networkx_labels(G_small,pos)
    plt.colorbar(nc1)
    plt.axis('off')
    
    plt.show(block=False)
    plt.close(fig)
    
def stability_plot(data, labels, name = None):
    jaccards = []
    for categ in data:       
        jk = []
        for i in range(len(categ)):
            for j in range(i+1,len(categ)):
                jk.append(jac(categ[i],categ[j]))
        jaccards.append(jk)
    fig, ax = plt.subplots(figsize=(15, 10))
    bplot1 = ax.boxplot(jaccards,
                             vert=True,  # vertical box alignment
                             patch_artist=True,  # fill with color
                             labels=labels)  # will be used to label x-ticks
    plt.ylabel("jaccard index")
    plt.xlabel("required gene module")
    plt.ylim(0,1)
    if name!=None:
        plt.savefig(name+".png")
    plt.show(block=False)
    plt.close(fig)


    
    
def clean_net(gene_groups,patients_groups, clusters,L_g,G,GE):    
    genes_components = []
    sizes = []
    for clust in range(clusters):
        group_g = gene_groups[clust]
        if clust == 0:
            not_clust = 1
        else:
            not_clust = 0
        if len(group_g)>=L_g:
            g = nx.subgraph(G,group_g)
            comp_big = max(nx.connected_component_subgraphs(g), key=len)
            dg = dict(nx.degree(comp_big))
            ones = [x for x in dg if dg[x]==1]
            nodes = list(comp_big.nodes)
            size_comp = len(nodes)
            max_out = len(nodes)- L_g
            while max_out >0:
                dif = np.mean(GE[patients_groups[clust]].loc[ones],axis = 1)-np.mean(GE[patients_groups[not_clust]].loc[ones],axis = 1)
                dif = dif.sort_values()
                ones = list(dif[dif<2].index)
                if len(ones)>0:
                    if len(ones)<=max_out:
                        outsiders = list(ones)
                    if len(ones) > max_out:
                        outsiders = list(ones)[:max_out]
     
                    nodes  = list(set(nodes) - set(outsiders))
                    g = nx.subgraph(G,nodes)
                    comp_big = g
                    dg = dict(nx.degree(comp_big))
                    ones = [x for x in dg if dg[x]==1]
                    nodes = list(comp_big.nodes)
                    size_comp = len(nodes)
                    max_out = len(nodes)- L_g
                else:
                    max_out = 0
                    
            group_g = nodes
        elif len(group_g)>0:
            g = nx.subgraph(G,group_g)
            comp_big = max(nx.connected_component_subgraphs(g), key=len)
            nodes = list(comp_big.nodes)
            size_comp = len(nodes)
        else:
            size_comp = 0
            
        genes_components.append(group_g)
        sizes.append(size_comp)
    return genes_components,sizes

def sim_data(genes1,genes2,background,patients1,patients2,dens):
    n = genes1+genes2+background
    m = patients1 +patients2
    
    
    genes = np.arange(n)
    groups_genes = list(np.ones(genes1))+list(np.ones(genes2)*2)+list(np.ones(background)*3)
    groups_p = [1 if node<patients1 else 2 for node in range(m)]
    
    to_sparce = 0.3 #to sparcify bipartite
    to_mix = 0.99 # to mix edges berween groups
    b = np.zeros((n,m))
    ge = np.random.normal(0,1,n*m).reshape(n,m)

    for patient in range(m):
        for gene in range(n):
            p_gr = groups_p[patient]
            g_gr = groups_genes[gene]
            if p_gr ==1 and g_gr == 1: #all up
                ge[gene,patient] = np.random.normal(1,0.35,1)
            elif p_gr ==2 and g_gr == 2:
                ge[gene,patient] = np.random.normal(1,0.35,1) #also up
            elif p_gr ==1 and g_gr == 2:
                ge[gene,patient] = np.random.normal(-1,0.35,1) #down
            elif p_gr ==2 and g_gr == 1:
                ge[gene,patient] = np.random.normal(-1,0.35,1) #down
    for patient in range(m):
        for gene in range(genes1+genes2):
            prob = np.random.uniform(0,1)
            if prob>0.9:
                ge[gene,patient] = np.random.normal(0,1,1)
                
    for gene in range(genes1+genes2,n):
        prob = np.random.uniform(0,1)
        if prob<0.05:


            for patient in range(m):
                if groups_p[patient] ==1: #all up
                    ge[gene,patient] = np.random.normal(0.3,0.35,1)
                else:
                    ge[gene,patient] = np.random.normal(-0.3,0.35,1)
        if prob>0.05 and prob<0.1:
            for patient in range(m):
                if groups_p[patient] ==1: #all up
                    ge[gene,patient] = np.random.normal(-0.3,0.35,1)
                else:
                    ge[gene,patient] = np.random.normal(0.3,0.35,1)
                    
                
    g1 = nx.barabasi_albert_graph(genes1,1)
    g2 = nx.barabasi_albert_graph(genes2,1)
    g3 = nx.barabasi_albert_graph(background, 1)
    G = nx.disjoint_union(g1,g2)
    G = nx.disjoint_union(G,g3)
    for _ in range( int(dens*n)):
        node1 = np.random.randint(0,genes1)
        node2 = np.random.randint(genes1,genes1+genes2)
        node3_1 = np.random.randint(genes1+genes2,n)
        node3_2 = np.random.randint(genes1+genes2,n)
        G.add_edges_from([(node1,node3_1),
                          (node2,node3_2)])
     
    d =  nx.density(G)
    count = 0 
    while d>0.002 and count<10:
        
        node3_1 = np.random.randint(genes1+genes2,n)
        node3_2 = np.random.randint(genes1+genes2,n)
        count = count+1
        if G.has_edge(node3_1,node3_2):
            G.remove_edge(node3_1,node3_2)
            d =  nx.density(G) 
        
    #A_g = nx.adj_matrix(G).todense() *1
    b_sp = csr_matrix(b) #sparse matrix for making bipartite graph
    B = bipartite.from_biadjacency_matrix(b_sp)
    
    GE = pd.DataFrame(ge,index = np.arange(n),columns = np.arange(n,n+m))
    H = HI_big(GE,1,1,1)
    return(B,GE,G,H,d,n,m)


@shared_task(name="algo_output_task")
def algo_output_task(s,L_g_min,L_g_max,expr_str,ppi_str,nbr_iter,nbr_ants,evap,epsilon,hi_sig,pher_sig):
	col = "cancer_type"
	size = 2000
	log2 = True
	with open("polls/static/output_console.txt", "w") as text_file:
   		text_file.write("Your files are being processed...")
	#B,G,H,n,m,GE,A_g,group1,group2,labels_B,rev_labels_B,val1,val2 = lib.aco_preprocessing(fh, prot_fh, col,log2 = True, gene_list = None, size = size, sample= None)
	B,G,H,n,m,GE,A_g,group1,group2,labels_B,rev_labels_B,val1,val2,group1_ids,group2_ids = lib.aco_preprocessing_strings(expr_str, ppi_str, col,log2 = True, gene_list = None, size = size, sample= None)
	print(group1_ids)	
	with open("polls/static/output_console.txt", "w") as text_file:
   		text_file.write("Starting model run...")	
	print("How many genes you want per cluster (minimum):")
	#L_g_min = int(input())
	print("How many genes you want per cluster (maximum):")
	#L_g_max = int(input())
	imp.reload(lib)
	
	# =============================================================================
	# #GENERAL PARAMETERS:
	# =============================================================================
	clusters = 2 # other options are currently unavailable
	K = int(nbr_ants) # number of ants
	eps = float(epsilon) # stopping criteria: score_max-score_av<eps
	b = float(hi_sig) #HI significance
	evaporation  = float(evap)
	a = float(pher_sig) #pheramone significance
	times =int(nbr_iter) #max amount of iterations
	#times =45 #max amount of iterations
	# =============================================================================
	# #NETWORK SIZE PARAMETERS:
	# =============================================================================
	cost_limit = 20 # will be determined authmatically soon. Aproximately the rule is that cost_limit = (geneset size)/100 but minimum  3
	#L_g_min = 10 # minimum # of genes per group
	#L_g_max = 20 # minimum # of genes per group
	
	th = 1# the coefficient to define the search radipus which is supposed to be bigger than 
	#mean(heruistic_information[patient]+th*std(heruistic_information[patient])
	#bigger th - less genes are considered (can lead to empty paths if th is too high)
	# will be also etermined authmatically soon. Right now the rule is such that th = 1 for genesets >1000
	#print(a)
	#print(type(a))
	#print(type(b))
	#print(type(n))
	#print(type(m))
	#print(type(H))
	#print(H)
	#print(type(GE))	
	#print(GE)
	#print(type(G))
	#print(G)
	#print(type(clusters))
	#print(clusters)
	
	
	with open("polls/static/output_console.txt", "w") as text_file:	
		text_file.write("Progress of the algorithm is shown below...")
	start = time.time()
	solution,t_best,sc,conv= lib.ants(a,b,n,m,H,GE,G,clusters,cost_limit,K,evaporation,th,L_g_min,L_g_max,eps,times,opt= None,pts = False,show_pher = False,show_plot = True, print_runs = False, save 	= None, show_nets = False)
	#result_99= ants_2.delay(a,b,n,m,H,GE,G,clusters,cost_limit,K,evaporation,th,L_g_min,L_g_max,eps,times,opt= None,pts = False,show_pher = False,show_plot = True, print_runs = True, save 	= None, show_nets = False)
	#print(result_99.get())
	#solution,t_best,sc,conv= result_99.get()
	#solution,t_best,sc,conv= ants_2.delay(a,b,n,m,H,GE,G,clusters,cost_limit,K,evaporation,th,L_g_min,L_g_max,eps,times,opt= None,pts = False,show_pher = False,show_plot = True, print_runs = False, save 	= None, show_nets = False)
	end = time.time()
	
	
	
	
	print("######################################################################")
	print("RESULTS ANALYSIS")
	print("total time " + str(round((end - start)/60,2))+ " minutes")
	print("jaccard indexes:")
	jacindices = lib.jac_matrix(solution[1],[group1,group2])
	print(jacindices)
	jac_1_ret = jacindices[0]
	jac_2_ret = jacindices[1]
	if lib.jac(group1,solution[1][0])>lib.jac(group1,solution[1][1]):
	    values = [val1,val2]
	else:
	    values = [val2,val1]
	# mapping to gene names (for now with API)
	mg = mygene.MyGeneInfo()
	new_genes = solution[0][0]+solution[0][1]
	new_genes_entrez = [labels_B[x] for x in new_genes]
	out = mg.querymany(new_genes_entrez, scopes='entrezgene', fields='symbol', species='human')
	mapping =dict()
	for line in out:
	    mapping[rev_labels_B[line["query"]]] = line["symbol"]
	###m plotting networks
	new_genes1 = [mapping[key] for key in mapping if key in solution[0][0] ]     
	new_genes2 = [mapping[key] for key in mapping if key in solution[0][1] ]    
	
	genes1,genes2 = solution[0]
	patients1, patients2 = solution[1]
	#print(patients1)
	#print(group1)
	means1 = [np.mean(GE[patients1].loc[gene])-np.mean(GE[patients2].loc[gene]) for gene in genes1]
	means2 = [np.mean(GE[patients1].loc[gene])-np.mean(GE[patients2].loc[gene]) for gene in genes2]
	
	G_small = nx.subgraph(G,genes1+genes2)
	G_small = nx.relabel_nodes(G_small,mapping)
	
	plt.figure(figsize=(15,15))
	cmap=plt.cm.RdYlGn
	vmin = -2
	vmax = 2
	pos = nx.spring_layout(G_small)
	ec = nx.draw_networkx_edges(G_small,pos)
	nc1 = nx.draw_networkx_nodes(G_small,nodelist =new_genes1, pos = pos,node_color=means1, node_size=600,alpha=1.0,
	                             vmin=vmin, vmax=vmax,node_shape = "^",cmap =cmap, label = values[0] )
	nc2 = nx.draw_networkx_nodes(G_small,nodelist =new_genes2, pos = pos,node_color=means2, node_size=600,
	                             alpha=1.0,
	                             vmin=vmin, vmax=vmax,node_shape = "o",cmap =cmap, label = values[1])
	nx.draw_networkx_labels(G_small,pos,font_size = 15,font_weight='bold')
	ret2 = means1 + means2
	ret3 = new_genes1 + new_genes2
	adjlist = []
	for line99 in nx.generate_edgelist(G_small,data=False):	
		lineSplit = line99.split()
		adjlist.append([lineSplit[0],lineSplit[1]])
	plt.legend(frameon  = True)
	plt.colorbar(nc1)
	plt.axis('off')
	print(list(G_small.nodes()))
	### plotting expression data
	plt.rc('font', size=30)          # controls default text sizes
	plt.rc('axes', titlesize=20)     # fontsize of the axes title
	plt.rc('axes', labelsize=20)    # fontsize of the x and y labels
	plt.rc('xtick', labelsize=15)    # fontsize of the tick labels
	plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
	plt.rc('legend', fontsize=30)
	
	
	grouping_p = []
	p_num = list(GE.columns)
	for p in p_num:
	    if p in solution[1][0]:
	        grouping_p.append(values[0])
	    else:
	        grouping_p.append(values[1])
	grouping_p = pd.DataFrame(grouping_p,index = p_num)
	grouping_g = []
	
	GE_small = GE.T[genes1+genes2]
	GE_small. rename(columns=mapping, inplace=True)
	GE_small = GE_small.T
	g_num = list(GE_small.index)
	for g in g_num:
	    if g in new_genes1 :
	        grouping_g.append(values[0])
	    elif  g in new_genes2 :
	        grouping_g.append(values[1])
	    else:
	        grouping_g.append(3)
	        
	grouping_g = pd.DataFrame(grouping_g,index = g_num)
	
	species = grouping_g[grouping_g[0]!=3][0]
	lut = {values[0]: '#4FB6D3', values[1]: '#22863E'}
	col_colors = species.map(lut)
	
	species = grouping_p[0]
	lut = {values[0]: '#4FB6D3',values[1]: '#22863E'}
	row_colors = species.map(lut)
	plt.savefig("/home/quirin/testproject/polls/static/ntw.png")
	#g = sns.clustermap(GE_small.T, row_colors=row_colors,col_colors = col_colors,figsize=(15, 10))
	#for label in values:
	#    g.ax_col_dendrogram.bar(0, 0, color=lut[label],
	#                            label=label, linewidth=0)
	#g.ax_col_dendrogram.legend(loc="upper center", ncol=2,bbox_to_anchor=(0.55, 1.5),
	#                            borderaxespad=0.)
	#ax = g.ax_heatmap
	#ax.set_xlabel("Genes")
	#ax.set_ylabel("Patients")
	
	#plotting convergence
	#fig, ax = plt.subplots(figsize=(10, 7))
	plt.clf()
	plt.boxplot(conv/2,
	                        vert=True,  # vertical box alignment
	                        patch_artist=True)  # will be used to label x-ticks
	
	plt.xlabel("iterations")
	plt.ylabel("score per subnetwork")
	#plt.show(bplot1)
	plt.savefig("/home/quirin/testproject/polls/static/conv.png")
	return(GE_small.T,row_colors,col_colors,G_small, ret2, ret3, adjlist,new_genes1,group1_ids,group2_ids,jac_1_ret,jac_2_ret)

@shared_task(name="script_output_task")
def script_output_task(T,row_colors1,col_colors1,G2,means,genes_all,adjlist,genes1):
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
	plot.add_tools(HoverTool(tooltips=None),TapTool(),BoxSelectTool())
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
	source = ColumnDataSource(data=dict(x=x, y=y,Name=[node_labels[i] for i in genes.keys()]))
	labels2 = LabelSet(x='x', y='y', text='Name',text_font_size="8pt",x_offset=-20, source=source,background_fill_color='white',background_fill_alpha=0.0,level='glyph',render_mode='canvas')
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
	with open("polls/static/output_console.txt", "w") as text_file:
        	text_file.write("Done!")
	return(script,div,plot_1)


@shared_task(name="script_output_task_2")
def script_output_task_2(T,row_colors1,col_colors1,G2,means,genes_all,adjlist,genes1,save_data,username):
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
	plot.add_tools(HoverTool(tooltips=None),TapTool(),BoxSelectTool())
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
	source = ColumnDataSource(data=dict(x=x, y=y,Name=[node_labels[i] for i in genes.keys()]))
	labels2 = LabelSet(x='x', y='y', text='Name',text_font_size="8pt",x_offset=-20, source=source,background_fill_color='white',background_fill_alpha=0.0,level='glyph',render_mode='canvas')
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
	if(save_data == "true"):
		foobar = "user_uploaded_files/" + username
		if not(os.path.isdir(foobar)):
			os.mkdir(foobar)
		filename_1 = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S_%f')[:-3]
		filename_2 = foobar + "/" + username + filename_1 + "_json.json"
		filename_3 = foobar + "/" + username + filename_1 + "_heatmap.png"
		copyfile("polls/static/test.png",filename_3)
		copyfile("polls/static/test15.json",filename_2)
	with open("polls/static/output_console.txt", "w") as text_file:
        	text_file.write("Done!")
	return(script,div,plot_1)


@shared_task(name="script_output_task_3")
def script_output_task_3(T,row_colors1,col_colors1,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr):
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
	plot.add_tools(HoverTool(tooltips=None),TapTool(),BoxSelectTool())
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
	source = ColumnDataSource(data=dict(x=x, y=y,Name=[node_labels[i] for i in genes.keys()]))
	labels2 = LabelSet(x='x', y='y', text='Name',text_font_size="8pt",x_offset=-20, source=source,background_fill_color='white',background_fill_alpha=0.0,level='glyph',render_mode='canvas')
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
	patientData = {}
	patientData_misc = {}
	clinicalLines = clinicalstr.split("\n")
	for line in clinicalLines:
		#bababa = line.replace('"','')
		bababa2 = line.split("\",\"")
		bababa2[0] = bababa2[0].replace('"','')
		if(len(bababa2) > 64):
			print(bababa2[57])
			#print(bababa2[72])
			if(bababa2[64] != "--" and bababa2[64] != "name:ch1" and bababa2[64] != "NA" and bababa2[64].replace('.','',1).isdigit()):
				patientData.update({bababa2[0]:bababa2[64]})
				print(patientData)
		if(len(bababa2) > 71):
			if(bababa2[57].replace('.','',1).isdigit() and bababa2[72].replace('.','',1).isdigit()):
				tmp97 = [float(bababa2[57]),float(bababa2[72])]
				patientData_misc.update({bababa2[0]:tmp97})
	age1 = 0
	ct1 = 0.001
	age2 = 0
	ct2 = 0.001
	size1 = 0
	ct3 = 0.001
	size2 = 0
	ct4 = 0.001
	for key in patientData_misc:
		if key in group1_ids:
			age1 = age1 + patientData_misc[key][0]
			ct1 = ct1 + 1
			size1 = size1 + patientData_misc[key][1]
			ct3 = ct3 + 1
		elif key in group2_ids:
			age2 = age2 + patientData_misc[key][0]
			ct2 = ct2 + 1
			size2 = size2 + patientData_misc[key][1]
			ct4 = ct4 + 1
	mean_age_1 = age1 / ct1
	mean_age_2 = age2 / ct2
	mean_size_1 = size1 / ct3
	mean_size_2 = size2 / ct4
	ret_metadata = []
	survival_1 = []
	survival_2 = []
	survival_mean_1 = 0
	survival_mean_2 = 0
	ctr_surv_1 = 0.001
	sum_surv_1 = 0
	ctr_surv_2 = 0.001
	sum_surv_2 = 0
	for key in patientData:
		if key in group1_ids:
			survival_1.append(patientData[key])
		elif key in group2_ids:
			survival_2.append(patientData[key])
	for elem in survival_1:
		sum_surv_1 = sum_surv_1 + float(elem)
		ctr_surv_1 = ctr_surv_1 + 1
	for elem in survival_2:
		sum_surv_2 = sum_surv_2 + float(elem)
		ctr_surv_2 = ctr_surv_2 + 1
	survival_mean_1 = sum_surv_1 / ctr_surv_1
	survival_mean_2 = sum_surv_2 / ctr_surv_2
	ret_metadata.append({'group':"Group 1",'survival':survival_mean_1,'age':mean_age_1,'size':mean_size_1})
	ret_metadata.append({'group':"Group 2",'survival':survival_mean_2,'age':mean_age_2,'size':mean_size_2})
	print(survival_mean_1)
	text_file_3 = open("polls/static/metadata.txt", "w")
	text_file_3.write("<table><tr><th>Patient Group</th><th>Mean Survival</th><th>mean age at diagnosis</th><th>mean tumor size</th></tr>")
	for elem in ret_metadata:
		text_file_3.write("<tr>")
		text_file_3.write("<th>" + str(elem['group']) + "</th>")
		text_file_3.write("<th>" + str(elem['survival']) + "</th>")
		text_file_3.write("<th>" + str(elem['age']) + "</th>")
		text_file_3.write("<th>" + str(elem['size']) + "</th>")
		text_file_3.write("</tr>")
	text_file_3.write("</table>")
	text_file_3.close()
	survival_perc_1 = {0:1}
	survival_perc_2 = {0:1}
	for i in range(1,10):	
		tmp1 = 1.0
		tmp2 = 1.0
		for k in survival_1:
			if(float(k) < float(i)):
				tmp1 = tmp1 - (1.0/len(survival_1))
		for k in survival_2:
			if(float(k) < float(i)):
				tmp2 = tmp2 - (1.0/len(survival_2))
		survival_perc_1.update({i:tmp1})
		survival_perc_2.update({i:tmp2})
	#output_notebook()	
	trace1 = go.Scatter(
	x=list(survival_perc_1.keys()),
	y=list(survival_perc_1.values()),
	mode='lines+markers',
	name="'Group 1'",
	hoverinfo='name',
	line=dict(
	shape='hv'))
	trace2 = go.Scatter(
	x=list(survival_perc_2.keys()),
	y=list(survival_perc_2.values()),
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
	#plot_div=plotly.offline.plot(fig, auto_open=False,include_plotlyjs = False, output_type='div')
	plot_div=plotly.offline.plot(fig, auto_open=False,output_type='div')
	with open("polls/static/output_plotly.html", "w") as text_file_2:
        	text_file_2.write(plot_div)
	#plotly.offline.plot(fig,auto_open=False, image_filename="tempimage.png", image='png')
	with open("polls/static/output_console.txt", "w") as text_file:
        	text_file.write("Done!")
	print(ret_metadata)
	return(script,div,plot_1,plot_div,ret_metadata)


@shared_task(name="script_output_task_4")
def script_output_task_4(T,row_colors1,col_colors1,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1_ret,jac_2_ret):
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
	plot.add_tools(HoverTool(tooltips=None),TapTool(),BoxSelectTool())
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
	source = ColumnDataSource(data=dict(x=x, y=y,Name=[node_labels[i] for i in genes.keys()]))
	labels2 = LabelSet(x='x', y='y', text='Name',text_font_size="8pt",x_offset=-20, source=source,background_fill_color='white',background_fill_alpha=0.0,level='glyph',render_mode='canvas')
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
	patientData = {}
	patientData_misc = {}
	group1_in_metadata = []
	group2_in_metadata = []
	metastasis = []
	no_metastasis = []
	group1_has_bm = []
	group2_has_bm = []
	bm = []
	no_bm = []
	group1_has_lm = []
	group2_has_lm = []
	lm = []
	no_lm = []
	group1_has_lymph = []
	group2_has_lymph = []
	lymph = []
	no_lymph = []
	group1_has_path_er = []
	group2_has_path_er = []
	path_er = []
	no_path_er = []
	group1_has_path_pr = []
	group2_has_path_pr = []
	path_pr = []
	no_path_pr = []
	group1_has_prognosis = []
	group2_has_prognosis = []
	prognosis_good = []
	prognosis_bad = []
	if(clinicalstr == "empty"):
		ret_metadata = []
		ret_metadata.append({'group':"Group 1",'lm':"NA",'bm':"NA",'lymph':"NA",'metastasis':"NA",'path_er':"NA",'path_pr':"NA",'good_prognosis':"NA"})
		ret_metadata.append({'group':"Group 2",'lm':"NA",'bm':"NA",'lymph':"NA",'metastasis':"NA",'path_er':"NA",'path_pr':"NA",'good_prognosis':"NA"})
		#ret_metadata.append({'group':"Group 1",'jac':jac_1_ret,'survival':"NA",'age':"NA",'size':"NA"})
		#ret_metadata.append({'group':"Group 2",'jac':jac_2_ret,'survival':"NA",'age':"NA",'size':"NA"})
		text_file_3 = open("polls/static/metadata.txt", "w")
		text_file_3.write("<table><tr><th>Patient Cluster</th><th>Lung metastasis</th><th>Breast metastasis</th><th>Lymph nodes</th><th>Metastasis</th><th>Path_ER</th><th>Path_PR</th><th>Good Prognosis</th></tr>")
		for elem in ret_metadata:
			text_file_3.write("<tr>")
			text_file_3.write("<th>" + str(elem['lm']) + "</th>")
			text_file_3.write("<th>" + str(elem['bm']) + "</th>")
			text_file_3.write("<th>" + str(elem['lymph']) + "</th>")
			text_file_3.write("<th>" + str(elem['metastasis']) + "</th>")
			text_file_3.write("<th>" + str(elem['path_er']) + "</th>")
			text_file_3.write("<th>" + str(elem['path_pr']) + "</th>")
			text_file_3.write("<th>" + str(elem['good_prognosis']) + "</th>")
			text_file_3.write("</tr>")
		text_file_3.write("</table>")
		text_file_3.close()	
		plot_div = ""
		with open("polls/static/output_plotly.html", "w") as text_file_2:
        		text_file_2.write("")
	
	else:
		clinicalLines = clinicalstr.split("\n")
		for line in clinicalLines:
			#bababa = line.replace('"','')
			bababa2 = line.split("\",\"")
			bababa2[0] = bababa2[0].replace('"','')
			if(len(bababa2) > 64):
				#print(bababa2[53])	
				#print(bababa2[63])
				#print(bababa2[72])
				if(bababa2[64] != "--" and bababa2[64] != "name:ch1" and bababa2[64] != "NA" and bababa2[64].replace('.','',1).isdigit()):
					patientData.update({bababa2[0]:bababa2[64]})
				#	print(patientData)
			#if(len(bababa2) > 71):
			if(len(bababa2) > 63):
				if(bababa2[63].isdigit()):
					print("patient id" + str(bababa2[0]))
					if(bababa2[0] in group1_ids):
						print("foo")
						group1_in_metadata.append(bababa2[0])
					elif(bababa2[0] in group2_ids):
						group2_in_metadata.append(bababa2[0])
					if(int(bababa2[63]) == 1):
						metastasis.append(bababa2[0])
					elif(int(bababa2[63]) == 0):
						no_metastasis.append(bababa2[0])
					tmp97 = [0.0,float(bababa2[63]),0.0]
					#patientData_misc.update({bababa2[0]:tmp97})
				if(bababa2[61].isdigit()):
					if(bababa2[0] in group1_ids):
						group1_has_lm.append(bababa2[0])
					elif(bababa2[0] in group2_ids):
						group1_has_lm.append(bababa2[0])
					if(int(bababa2[61]) == 1):
						lm.append(bababa2[0])
					elif(int(bababa2[61]) == 0):
						no_lm.append(bababa2[0])
				if(bababa2[58].isdigit()):
					if(bababa2[0] in group1_ids):
						group1_has_bm.append(bababa2[0])
					elif(bababa2[0] in group2_ids):
						group1_has_bm.append(bababa2[0])
					if(int(bababa2[58]) == 1):
						bm.append(bababa2[0])
					elif(int(bababa2[58]) == 0):
						no_bm.append(bababa2[0])
				if(len(bababa2[15]) > 1 and len(bababa2[15]) < 20):
					if(bababa2[15].split(": ")[1] == "N"):
						if(bababa2[0] in group1_ids):
							group1_has_path_er.append(bababa2[0])
						elif(bababa2[0] in group2_ids):
							group2_has_path_er.append(bababa2[0])
						no_path_er.append(bababa2[0])
					if(bababa2[15].split(": ")[1] == "P"):
						if(bababa2[0] in group1_ids):
							group1_has_path_er.append(bababa2[0])
						elif(bababa2[0] in group2_ids):
							group2_has_path_er.append(bababa2[0])
						path_er.append(bababa2[0])
				if(len(bababa2[16]) > 1 and len(bababa2[16]) < 20):
					if(bababa2[16].split(": ")[1] == "N"):
						if(bababa2[0] in group1_ids):
							group1_has_path_pr.append(bababa2[0])
						elif(bababa2[0] in group2_ids):
							group2_has_path_pr.append(bababa2[0])
						no_path_pr.append(bababa2[0])
					if(bababa2[16].split(": ")[1] == "P"):
						if(bababa2[0] in group1_ids):
							group1_has_path_pr.append(bababa2[0])
						elif(bababa2[0] in group2_ids):
							group2_has_path_pr.append(bababa2[0])
						path_pr.append(bababa2[0])
				if(len(bababa2) > 71):	
					if(bababa2[0] in group1_ids):
						group1_has_prognosis.append(bababa2[0])
					elif(bababa2[0] in group2_ids):
						group2_has_prognosis.append(bababa2[0])
					if(bababa2[71] == "Good Prognosis"):
						prognosis_good.append(bababa2[0])
					elif(bababa2[71] == "Bad Prognosis"):
						prognosis_bad.append(bababa2[0])
				if(len(bababa2) > 70):				
					if(bababa2[70].isdigit()):
						if(bababa2[0] in group1_ids):
							group1_has_lymph.append(bababa2[0])
						elif(bababa2[0] in group2_ids):
							group1_has_lymph.append(bababa2[0])
						if(int(bababa2[70]) == 1):
							lymph.append(bababa2[0])
						elif(int(bababa2[70]) == 0):
							no_lymph.append(bababa2[0])
				elif(bababa2[53].replace('.','',1).isdigit() and bababa2[63].isdigit()):
					print("patient id" + str(bababa2[0]))
					if(bababa2[0] in group1_ids):
						print("foo")
						group1_in_metadata.append(bababa2[0])
					elif(bababa2[0] in group2_ids):
						group2_in_metadata.append(bababa2[0])
					if(int(bababa2[63]) == 1):
						metastasis.append(bababa2[0])
					elif(int(bababa2[63]) == 0):
						no_metastasis.append(bababa2[0])
					#tmp97 = [float(bababa2[53]),float(bababa2[63]),0.0]
					#patientData_misc.update({bababa2[0]:tmp97})
				
		age1 = 0
		ct1 = 0.001
		age2 = 0
		ct2 = 0.001
		print(patientData_misc)
		#for key in patientData_misc:
		#	if key in group1_ids:
		#		age1 = age1 + patientData_misc[key][0]
		#		ct1 = ct1 + 1
		#		size1 = size1 + patientData_misc[key][2]
		#		ct3 = ct3 + 1
		#		ct5 = ct5 + 1
		#		metas1 = metas1 + patientData_misc[key][1]
		#	elif key in group2_ids:
		#		age2 = age2 + patientData_misc[key][0]
		#		ct2 = ct2 + 1
		#		size2 = size2 + patientData_misc[key][2]
		#		ct4 = ct4 + 1
		#		ct6 = ct6 + 1
		#		metas2 = metas2 + patientData_misc[key][1]
		ret_metadata = []
		nbr_patients = float(len(group1_ids) + len(group2_ids))
		survival_1 = []
		survival_2 = []
		survival_mean_1 = 0
		survival_mean_2 = 0
		ctr_surv_1 = 0.001
		sum_surv_1 = 0
		ctr_surv_2 = 0.001
		sum_surv_2 = 0
		for key in patientData:
			if key in group1_ids:
				survival_1.append(patientData[key])
			elif key in group2_ids:
				survival_2.append(patientData[key])
		for elem in survival_1:
			sum_surv_1 = sum_surv_1 + float(elem)
			ctr_surv_1 = ctr_surv_1 + 1
		for elem in survival_2:
			sum_surv_2 = sum_surv_2 + float(elem)
			ctr_surv_2 = ctr_surv_2 + 1
		survival_mean_1 = sum_surv_1 / ctr_surv_1
		survival_mean_2 = sum_surv_2 / ctr_surv_2
		#print(group1_in_metadata)
		#print(metastasis)

		if(float(len(group1_in_metadata) + len(group2_in_metadata)) / nbr_patients > 0.8):
			metas_jac = lib.jac(group1_in_metadata,metastasis)
			metas_jac_2 = lib.jac(group2_in_metadata,metastasis)
		else:
			metas_jac = "NA"
			metas_jac_2 = "NA"
		if(float(len(group1_has_lm) + len(group2_has_lm)) / nbr_patients > 0.8):
			lm_jac = lib.jac(group1_has_lm,lm)
			lm_jac_2 = lib.jac(group2_has_lm,lm)
		else:
			lm_jac = "NA"
			lm_jac_2 = "NA"
		if(float(len(group1_has_bm) + len(group2_has_bm)) / nbr_patients > 0.8):
			bm_jac = lib.jac(group1_has_bm,bm)
			bm_jac_2 = lib.jac(group2_has_bm,bm)
		else:
			bm_jac = "NA"
			bm_jac_2 = "NA"
		if(float(len(group1_has_lymph) + len(group2_has_lymph)) / nbr_patients > 0.8):
			lymph_jac = lib.jac(group1_has_lymph,lymph)
			lymph_jac_2 = lib.jac(group2_has_lymph,lymph)
		else:
			lymph_jac = "NA"
			lymph_jac_2 = "NA"
		if(float(len(group1_has_path_er) + len(group2_has_path_er)) / nbr_patients > 0.8):
			path_er_jac = lib.jac(group1_has_path_er,path_er)
			path_er_jac_2 = lib.jac(group2_has_path_er,path_er)
		else:
			path_er_jac = "NA"
			path_er_jac_2 = "NA"
		if(float(len(group1_has_path_er) + len(group2_has_path_er)) / nbr_patients > 0.8):
			path_pr_jac = lib.jac(group1_has_path_pr,path_pr)
			path_pr_jac_2 = lib.jac(group2_has_path_pr,path_pr)
		else:
			path_pr_jac = "NA"
			path_pr_jac_2 = "NA"
		if(float(len(group1_has_prognosis) + len(group2_has_prognosis)) / nbr_patients > 0.8):
			prognosis_jac = lib.jac(group1_has_prognosis,prognosis_good)
			prognosis_jac_2 = lib.jac(group2_has_prognosis,prognosis_good)
		else:
			prognosis_jac = "NA"
			prognosis_jac_2 = "NA"
		ret_metadata.append({'group':"Group 1",'lm':lm_jac,'bm':bm_jac,'lymph':lymph_jac,'metastasis':metas_jac,'path_er':path_er_jac,'path_pr':path_pr_jac,'good_prognosis':prognosis_jac})
		ret_metadata.append({'group':"Group 2",'lm':lm_jac_2,'bm':bm_jac_2,'lymph':lymph_jac_2,'metastasis':metas_jac_2,'path_er':path_er_jac_2,'path_pr':path_pr_jac_2,'good_prognosis':prognosis_jac_2})
		#ret_metadata.append({'group':"Group 1",'jac':jac_1_ret,'survival':survival_mean_1,'age':mean_age_1,'metastasis':metas_mean_1,'size':mean_size_1})
		#ret_metadata.append({'group':"Group 2",'jac':jac_2_ret,'survival':survival_mean_2,'age':mean_age_2,'metastasis':metas_mean_2,'size':mean_size_2})
		print(survival_mean_1)
		text_file_3 = open("polls/static/metadata.txt", "w")
		#text_file_3.write("<table><tr><th>Patient Cluster</th><th>Correlation with clinical Group</th><th>Mean Survival</th><th>mean age at diagnosis</th><th>mean tumor size</th></tr>")
		text_file_3.write("<table><tr><th>Patient Cluster</th><th>Lung metastasis</th><th>Breast metastasis</th><th>Lymph nodes</th><th>Metastasis</th><th>Path_ER</th><th>Path_PR</th><th>Good Prognosis</th></tr>")
		for elem in ret_metadata:
			text_file_3.write("<tr>")
			text_file_3.write("<th>" + str(elem['lm']) + "</th>")
			text_file_3.write("<th>" + str(elem['bm']) + "</th>")
			text_file_3.write("<th>" + str(elem['lymph']) + "</th>")
			text_file_3.write("<th>" + str(elem['metastasis']) + "</th>")
			text_file_3.write("<th>" + str(elem['path_er']) + "</th>")
			text_file_3.write("<th>" + str(elem['path_pr']) + "</th>")
			text_file_3.write("<th>" + str(elem['good_prognosis']) + "</th>")
			text_file_3.write("</tr>")
		text_file_3.write("</table>")
		text_file_3.close()
		survival_perc_1 = {0:1}
		survival_perc_2 = {0:1}
		for i in range(1,10):	
			tmp1 = 1.0
			tmp2 = 1.0
			for k in survival_1:
				if(float(k) < float(i)):
					tmp1 = tmp1 - (1.0/len(survival_1))
			for k in survival_2:
				if(float(k) < float(i)):
					tmp2 = tmp2 - (1.0/len(survival_2))
			survival_perc_1.update({i:tmp1})
			survival_perc_2.update({i:tmp2})
		#output_notebook()	
		trace1 = go.Scatter(
		x=list(survival_perc_1.keys()),
		y=list(survival_perc_1.values()),
		mode='lines+markers',
		name="'Group 1'",
		hoverinfo='name',
		line=dict(
		shape='hv'))
		trace2 = go.Scatter(
		x=list(survival_perc_2.keys()),
		y=list(survival_perc_2.values()),
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
		#plot_div=plotly.offline.plot(fig, auto_open=False,include_plotlyjs = False, output_type='div')
		plot_div=plotly.offline.plot(fig, auto_open=False,output_type='div')
		with open("polls/static/output_plotly.html", "w") as text_file_2:
        		text_file_2.write(plot_div)
		#plotly.offline.plot(fig,auto_open=False, image_filename="tempimage.png", image='png')
	with open("polls/static/output_console.txt", "w") as text_file:
        	text_file.write("Done!")
	print(ret_metadata)
	return(script,div,plot_1,plot_div,ret_metadata)


@shared_task(name="script_output_task_8")
def script_output_task_8(T,row_colors1,col_colors1,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1_ret,jac_2_ret,survival_col):
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
	plot.add_tools(HoverTool(tooltips=None),TapTool(),BoxSelectTool())
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
	source = ColumnDataSource(data=dict(x=x, y=y,Name=[node_labels[i] for i in genes.keys()]))
	labels2 = LabelSet(x='x', y='y', text='Name',text_font_size="8pt",x_offset=-20, source=source,background_fill_color='white',background_fill_alpha=0.0,level='glyph',render_mode='canvas')
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
	patientData = {}
	patientData_misc = {}
	group1_in_metadata = []
	group2_in_metadata = []
	metastasis = []
	no_metastasis = []
	group1_has_bm = []
	group2_has_bm = []
	bm = []
	no_bm = []
	group1_has_lm = []
	group2_has_lm = []
	lm = []
	no_lm = []
	group1_has_lymph = []
	group2_has_lymph = []
	lymph = []
	no_lymph = []
	group1_has_path_er = []
	group2_has_path_er = []
	path_er = []
	no_path_er = []
	group1_has_path_pr = []
	group2_has_path_pr = []
	path_pr = []
	no_path_pr = []
	group1_has_prognosis = []
	group2_has_prognosis = []
	prognosis_good = []
	prognosis_bad = []
	if(clinicalstr == "empty"):
		ret_metadata = []
		ret_metadata.append({'group':"Group 1",'lm':"NA",'bm':"NA",'lymph':"NA",'metastasis':"NA",'path_er':"NA",'path_pr':"NA",'good_prognosis':"NA"})
		ret_metadata.append({'group':"Group 2",'lm':"NA",'bm':"NA",'lymph':"NA",'metastasis':"NA",'path_er':"NA",'path_pr':"NA",'good_prognosis':"NA"})
		#ret_metadata.append({'group':"Group 1",'jac':jac_1_ret,'survival':"NA",'age':"NA",'size':"NA"})
		#ret_metadata.append({'group':"Group 2",'jac':jac_2_ret,'survival':"NA",'age':"NA",'size':"NA"})
		text_file_3 = open("polls/static/metadata.txt", "w")
		text_file_3.write("<table><tr><th>Patient Cluster</th><th>Lung metastasis</th><th>Breast metastasis</th><th>Lymph nodes</th><th>Metastasis</th><th>Path_ER</th><th>Path_PR</th><th>Good Prognosis</th></tr>")
		for elem in ret_metadata:
			text_file_3.write("<tr>")
			text_file_3.write("<th>" + str(elem['lm']) + "</th>")
			text_file_3.write("<th>" + str(elem['bm']) + "</th>")
			text_file_3.write("<th>" + str(elem['lymph']) + "</th>")
			text_file_3.write("<th>" + str(elem['metastasis']) + "</th>")
			text_file_3.write("<th>" + str(elem['path_er']) + "</th>")
			text_file_3.write("<th>" + str(elem['path_pr']) + "</th>")
			text_file_3.write("<th>" + str(elem['good_prognosis']) + "</th>")
			text_file_3.write("</tr>")
		text_file_3.write("</table>")
		text_file_3.close()	
		plot_div = ""
		with open("polls/static/output_plotly.html", "w") as text_file_2:
        		text_file_2.write("")
	
	else:
		clinicalLines = clinicalstr.split("\n")
		title_col = clinicalLines[0].split(",")
		print(title_col)
		survival_col_nbr = 64
		if(survival_col in title_col):
			survival_col_nbr = title_col.index(survival_col) + 1
			print("survival col")
		for line in clinicalLines:
			#bababa = line.replace('"','')
			bababa2 = line.split("\",\"")
			bababa2[0] = bababa2[0].replace('"','')
			if(len(bababa2) > survival_col_nbr):
				#print(bababa2[53])	
				#print(bababa2[63])
				#print(bababa2[72])
				if(bababa2[survival_col_nbr] != "--" and bababa2[survival_col_nbr] != "name:ch1" and bababa2[survival_col_nbr] != "NA" and bababa2[survival_col_nbr].replace('.','',1).isdigit()):
					patientData.update({bababa2[0]:bababa2[survival_col_nbr]})
				#	print(patientData)
			#if(len(bababa2) > 71):
			if(len(bababa2) > 63):
				if(bababa2[63].isdigit()):
					print("patient id" + str(bababa2[0]))
					if(bababa2[0] in group1_ids):
						print("foo")
						group1_in_metadata.append(bababa2[0])
					elif(bababa2[0] in group2_ids):
						group2_in_metadata.append(bababa2[0])
					if(int(bababa2[63]) == 1):
						metastasis.append(bababa2[0])
					elif(int(bababa2[63]) == 0):
						no_metastasis.append(bababa2[0])
					tmp97 = [0.0,float(bababa2[63]),0.0]
					#patientData_misc.update({bababa2[0]:tmp97})
				if(bababa2[61].isdigit()):
					if(bababa2[0] in group1_ids):
						group1_has_lm.append(bababa2[0])
					elif(bababa2[0] in group2_ids):
						group1_has_lm.append(bababa2[0])
					if(int(bababa2[61]) == 1):
						lm.append(bababa2[0])
					elif(int(bababa2[61]) == 0):
						no_lm.append(bababa2[0])
				if(bababa2[58].isdigit()):
					if(bababa2[0] in group1_ids):
						group1_has_bm.append(bababa2[0])
					elif(bababa2[0] in group2_ids):
						group1_has_bm.append(bababa2[0])
					if(int(bababa2[58]) == 1):
						bm.append(bababa2[0])
					elif(int(bababa2[58]) == 0):
						no_bm.append(bababa2[0])
				if(len(bababa2[15]) > 1 and len(bababa2[15]) < 20):
					if(bababa2[15].split(": ")[1] == "N"):
						if(bababa2[0] in group1_ids):
							group1_has_path_er.append(bababa2[0])
						elif(bababa2[0] in group2_ids):
							group2_has_path_er.append(bababa2[0])
						no_path_er.append(bababa2[0])
					if(bababa2[15].split(": ")[1] == "P"):
						if(bababa2[0] in group1_ids):
							group1_has_path_er.append(bababa2[0])
						elif(bababa2[0] in group2_ids):
							group2_has_path_er.append(bababa2[0])
						path_er.append(bababa2[0])
				if(len(bababa2[16]) > 1 and len(bababa2[16]) < 20):
					if(bababa2[16].split(": ")[1] == "N"):
						if(bababa2[0] in group1_ids):
							group1_has_path_pr.append(bababa2[0])
						elif(bababa2[0] in group2_ids):
							group2_has_path_pr.append(bababa2[0])
						no_path_pr.append(bababa2[0])
					if(bababa2[16].split(": ")[1] == "P"):
						if(bababa2[0] in group1_ids):
							group1_has_path_pr.append(bababa2[0])
						elif(bababa2[0] in group2_ids):
							group2_has_path_pr.append(bababa2[0])
						path_pr.append(bababa2[0])
				if(len(bababa2) > 71):	
					if(bababa2[0] in group1_ids):
						group1_has_prognosis.append(bababa2[0])
					elif(bababa2[0] in group2_ids):
						group2_has_prognosis.append(bababa2[0])
					if(bababa2[71] == "Good Prognosis"):
						prognosis_good.append(bababa2[0])
					elif(bababa2[71] == "Bad Prognosis"):
						prognosis_bad.append(bababa2[0])
				if(len(bababa2) > 70):				
					if(bababa2[70].isdigit()):
						if(bababa2[0] in group1_ids):
							group1_has_lymph.append(bababa2[0])
						elif(bababa2[0] in group2_ids):
							group1_has_lymph.append(bababa2[0])
						if(int(bababa2[70]) == 1):
							lymph.append(bababa2[0])
						elif(int(bababa2[70]) == 0):
							no_lymph.append(bababa2[0])
				elif(bababa2[53].replace('.','',1).isdigit() and bababa2[63].isdigit()):
					print("patient id" + str(bababa2[0]))
					if(bababa2[0] in group1_ids):
						print("foo")
						group1_in_metadata.append(bababa2[0])
					elif(bababa2[0] in group2_ids):
						group2_in_metadata.append(bababa2[0])
					if(int(bababa2[63]) == 1):
						metastasis.append(bababa2[0])
					elif(int(bababa2[63]) == 0):
						no_metastasis.append(bababa2[0])
					#tmp97 = [float(bababa2[53]),float(bababa2[63]),0.0]
					#patientData_misc.update({bababa2[0]:tmp97})
				
		age1 = 0
		ct1 = 0.001
		age2 = 0
		ct2 = 0.001
		print(patientData_misc)
		#for key in patientData_misc:
		#	if key in group1_ids:
		#		age1 = age1 + patientData_misc[key][0]
		#		ct1 = ct1 + 1
		#		size1 = size1 + patientData_misc[key][2]
		#		ct3 = ct3 + 1
		#		ct5 = ct5 + 1
		#		metas1 = metas1 + patientData_misc[key][1]
		#	elif key in group2_ids:
		#		age2 = age2 + patientData_misc[key][0]
		#		ct2 = ct2 + 1
		#		size2 = size2 + patientData_misc[key][2]
		#		ct4 = ct4 + 1
		#		ct6 = ct6 + 1
		#		metas2 = metas2 + patientData_misc[key][1]
		ret_metadata = []
		nbr_patients = float(len(group1_ids) + len(group2_ids))
		survival_1 = []
		survival_2 = []
		survival_mean_1 = 0
		survival_mean_2 = 0
		ctr_surv_1 = 0.001
		sum_surv_1 = 0
		ctr_surv_2 = 0.001
		sum_surv_2 = 0
		for key in patientData:
			if key in group1_ids:
				survival_1.append(patientData[key])
			elif key in group2_ids:
				survival_2.append(patientData[key])
		for elem in survival_1:
			sum_surv_1 = sum_surv_1 + float(elem)
			ctr_surv_1 = ctr_surv_1 + 1
		for elem in survival_2:
			sum_surv_2 = sum_surv_2 + float(elem)
			ctr_surv_2 = ctr_surv_2 + 1
		survival_mean_1 = sum_surv_1 / ctr_surv_1
		survival_mean_2 = sum_surv_2 / ctr_surv_2
		#print(group1_in_metadata)
		#print(metastasis)
		if(float(len(group1_in_metadata) + len(group2_in_metadata)) / nbr_patients > 0.8):
			metas_jac = lib.jac(group1_in_metadata,metastasis)
			metas_jac_2 = lib.jac(group2_in_metadata,metastasis)
		else:
			metas_jac = "NA"
			metas_jac_2 = "NA"
		if(float(len(group1_has_lm) + len(group2_has_lm)) / nbr_patients > 0.8):
			lm_jac = lib.jac(group1_has_lm,lm)
			lm_jac_2 = lib.jac(group2_has_lm,lm)
		else:
			lm_jac = "NA"
			lm_jac_2 = "NA"
		if(float(len(group1_has_bm) + len(group2_has_bm)) / nbr_patients > 0.8):
			bm_jac = lib.jac(group1_has_bm,bm)
			bm_jac_2 = lib.jac(group2_has_bm,bm)
		else:
			bm_jac = "NA"
			bm_jac_2 = "NA"
		if(float(len(group1_has_lymph) + len(group2_has_lymph)) / nbr_patients > 0.8):
			lymph_jac = lib.jac(group1_has_lymph,lymph)
			lymph_jac_2 = lib.jac(group2_has_lymph,lymph)
		else:
			lymph_jac = "NA"
			lymph_jac_2 = "NA"
		if(float(len(group1_has_path_er) + len(group2_has_path_er)) / nbr_patients > 0.8):
			path_er_jac = lib.jac(group1_has_path_er,path_er)
			path_er_jac_2 = lib.jac(group2_has_path_er,path_er)
		else:
			path_er_jac = "NA"
			path_er_jac_2 = "NA"
		if(float(len(group1_has_path_er) + len(group2_has_path_er)) / nbr_patients > 0.8):
			path_pr_jac = lib.jac(group1_has_path_pr,path_pr)
			path_pr_jac_2 = lib.jac(group2_has_path_pr,path_pr)
		else:
			path_pr_jac = "NA"
			path_pr_jac_2 = "NA"
		if(float(len(group1_has_prognosis) + len(group2_has_prognosis)) / nbr_patients > 0.8):
			prognosis_jac = lib.jac(group1_has_prognosis,prognosis_good)
			prognosis_jac_2 = lib.jac(group2_has_prognosis,prognosis_good)
		else:
			prognosis_jac = "NA"
			prognosis_jac_2 = "NA"
		ret_metadata.append({'group':"Group 1",'lm':lm_jac,'bm':bm_jac,'lymph':lymph_jac,'metastasis':metas_jac,'path_er':path_er_jac,'path_pr':path_pr_jac,'good_prognosis':prognosis_jac})
		ret_metadata.append({'group':"Group 2",'lm':lm_jac_2,'bm':bm_jac_2,'lymph':lymph_jac_2,'metastasis':metas_jac_2,'path_er':path_er_jac_2,'path_pr':path_pr_jac_2,'good_prognosis':prognosis_jac_2})
		#ret_metadata.append({'group':"Group 1",'jac':jac_1_ret,'survival':survival_mean_1,'age':mean_age_1,'metastasis':metas_mean_1,'size':mean_size_1})
		#ret_metadata.append({'group':"Group 2",'jac':jac_2_ret,'survival':survival_mean_2,'age':mean_age_2,'metastasis':metas_mean_2,'size':mean_size_2})
		print(survival_mean_1)
		text_file_3 = open("polls/static/metadata.txt", "w")
		#text_file_3.write("<table><tr><th>Patient Cluster</th><th>Correlation with clinical Group</th><th>Mean Survival</th><th>mean age at diagnosis</th><th>mean tumor size</th></tr>")
		text_file_3.write("<table><tr><th>Patient Cluster</th><th>Lung metastasis</th><th>Breast metastasis</th><th>Lymph nodes</th><th>Metastasis</th><th>Path_ER</th><th>Path_PR</th><th>Good Prognosis</th></tr>")
		for elem in ret_metadata:
			text_file_3.write("<tr>")
			text_file_3.write("<th>" + str(elem['lm']) + "</th>")
			text_file_3.write("<th>" + str(elem['bm']) + "</th>")
			text_file_3.write("<th>" + str(elem['lymph']) + "</th>")
			text_file_3.write("<th>" + str(elem['metastasis']) + "</th>")
			text_file_3.write("<th>" + str(elem['path_er']) + "</th>")
			text_file_3.write("<th>" + str(elem['path_pr']) + "</th>")
			text_file_3.write("<th>" + str(elem['good_prognosis']) + "</th>")
			text_file_3.write("</tr>")
		text_file_3.write("</table>")
		text_file_3.close()
		survival_perc_1 = {0:1}
		survival_perc_2 = {0:1}
		for i in range(1,10):	
			tmp1 = 1.0
			tmp2 = 1.0
			for k in survival_1:
				if(float(k) < float(i)):
					tmp1 = tmp1 - (1.0/len(survival_1))
			for k in survival_2:
				if(float(k) < float(i)):
					tmp2 = tmp2 - (1.0/len(survival_2))
			survival_perc_1.update({i:tmp1})
			survival_perc_2.update({i:tmp2})
		#output_notebook()	
		trace1 = go.Scatter(
		x=list(survival_perc_1.keys()),
		y=list(survival_perc_1.values()),
		mode='lines+markers',
		name="'Group 1'",
		hoverinfo='name',
		line=dict(
		shape='hv'))
		trace2 = go.Scatter(
		x=list(survival_perc_2.keys()),
		y=list(survival_perc_2.values()),
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
		#plot_div=plotly.offline.plot(fig, auto_open=False,include_plotlyjs = False, output_type='div')
		plot_div=plotly.offline.plot(fig, auto_open=False,output_type='div')
		with open("polls/static/output_plotly.html", "w") as text_file_2:
        		text_file_2.write(plot_div)
		#plotly.offline.plot(fig,auto_open=False, image_filename="tempimage.png", image='png')
	with open("polls/static/output_console.txt", "w") as text_file:
        	text_file.write("Done!")
	print(ret_metadata)
	return(script,div,plot_1,plot_div,ret_metadata)


@shared_task(name="script_output_task_9")
def script_output_task_9(T,row_colors1,col_colors1,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1_ret,jac_2_ret,survival_col,clinicaldf):
	#G = nx.Graph()
	#G_2 = nx.Graph()
	print("basdbrasdbawrubaerbsrbjsdb")
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
	plot.add_tools(HoverTool(tooltips=None),TapTool(),BoxSelectTool())
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
	source = ColumnDataSource(data=dict(x=x, y=y,Name=[node_labels[i] for i in genes.keys()]))
	labels2 = LabelSet(x='x', y='y', text='Name',text_font_size="8pt",x_offset=-20, source=source,background_fill_color='white',background_fill_alpha=0.0,level='glyph',render_mode='canvas')
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
	print("in method")
	patientData = {}
	patientData_misc = {}
	group1_in_metadata = []
	group2_in_metadata = []
	metastasis = []
	no_metastasis = []
	group1_has_bm = []
	group2_has_bm = []
	bm = []
	no_bm = []
	group1_has_lm = []
	group2_has_lm = []
	lm = []
	no_lm = []
	group1_has_lymph = []
	group2_has_lymph = []
	lymph = []
	no_lymph = []
	group1_has_path_er = []
	group2_has_path_er = []
	path_er = []
	no_path_er = []
	group1_has_path_pr = []
	group2_has_path_pr = []
	path_pr = []
	no_path_pr = []
	group1_has_prognosis = []
	group2_has_prognosis = []
	prognosis_good = []
	prognosis_bad = []
	with open("genelist.txt","w") as text_file_4:
		for i in G_list:
			text_file_4.write(str(i) + "\n")
	text_file_4.close()
	with open("genelist_1.txt","w") as text_file_5:
		for i in G_list:
			if(i in genes1):
				text_file_5.write(str(i) + "\n")
	text_file_5.close()
	with open("genelist_2.txt","w") as text_file_6:
		for i in G_list:
			if(i not in genes1):
				text_file_6.write(str(i) + "\n")
	text_file_6.close()
	if(clinicalstr == "empty"):
		ret_metadata = []
		ret_metadata.append({'group':"Group 1",'lm':"NA",'bm':"NA",'lymph':"NA",'metastasis':"NA",'path_er':"NA",'path_pr':"NA",'good_prognosis':"NA"})
		ret_metadata.append({'group':"Group 2",'lm':"NA",'bm':"NA",'lymph':"NA",'metastasis':"NA",'path_er':"NA",'path_pr':"NA",'good_prognosis':"NA"})
		ret_metadata.append({'group':"Group 1",'jac':jac_1_ret,'survival':"NA",'age':"NA",'size':"NA"})
		ret_metadata.append({'group':"Group 2",'jac':jac_2_ret,'survival':"NA",'age':"NA",'size':"NA"})
		text_file_3 = open("polls/static/metadata.txt", "w")
		text_file_3.write("<table><tr><th>Patient Cluster</th></tr>")
		text_file_3.write("<tr><th>Group1</th></tr>")
		text_file_3.write("</table>")
		#text_file_3.write("<table><tr><th>Patient Cluster</th><th>Lung metastasis</th><th>Breast metastasis</th><th>Lymph nodes</th><th>Metastasis</th><th>Path_ER</th><th>Path_PR</th><th>Good Prognosis</th></tr>")
		#for elem in ret_metadata:
		#	text_file_3.write("<tr>")
		#	text_file_3.write("<th>" + str(elem['lm']) + "</th>")
		#	text_file_3.write("<th>" + str(elem['bm']) + "</th>")
		#	text_file_3.write("<th>" + str(elem['lymph']) + "</th>")
		#	text_file_3.write("<th>" + str(elem['metastasis']) + "</th>")
		#	text_file_3.write("<th>" + str(elem['path_er']) + "</th>")
		#	text_file_3.write("<th>" + str(elem['path_pr']) + "</th>")
		#	text_file_3.write("<th>" + str(elem['good_prognosis']) + "</th>")
		#	text_file_3.write("</tr>")
		#text_file_3.write("</table>")
		text_file_3.close()	
		plot_div = ""
		with open("polls/static/output_plotly.html", "w") as text_file_2:
        		text_file_2.write("")
	
	else:
		print("barfbasdbhasdb")
		clinicalLines = clinicalstr.split("\n")
		title_col = clinicalLines[0].split(",")
		print(title_col)
		survival_col_nbr = 64
		clinicaldf.replace(['NaN','nan','?','--'],['NA','NA','NA','NA'], inplace=True)
		clinicaldf.replace(['NTL'],['NA'], inplace=True)
		clinicaldf.replace(['na'],['NA'], inplace=True)
		#print(clinicaldf.index)
		print(clinicaldf.columns)
		clinicaldf_col_names = list(clinicaldf.columns)
		print(clinicaldf_col_names)
		patientids_metadata = clinicaldf.iloc[:,0].values.tolist()
		if("Unnamed" in "\t".join(list(clinicaldf.columns))):
			print("basdbasdfbasfbafs")
			clinicaldf_col_names_temp = ['empty']
			clinicaldf_col_names_new = clinicaldf_col_names_temp + clinicaldf_col_names
			clinicaldf.columns = list(clinicaldf_col_names_new[:-1])
			print(clinicaldf_col_names_new)
		if("GSM" not in patientids_metadata[0]):
			patientids_metadata = list(clinicaldf.index)
		#print(clinicaldf["relapse (event=1; no event=0):ch1"])
		param_names = []
		param_values = []
		param_cols = []
		ctr = 0
		print(patientids_metadata)
		patients_0 = []
		patients_1 = []
		group1_has = []
		group2_has = []
		for column_name, column in clinicaldf.transpose().iterrows():
			column.fillna("NA",inplace=True)
			#print(column_name)
			#print(column)
			coluniq = column.unique()
			#print(coluniq)
			#print(column_name)
			#print(list(column))
			for elem in coluniq:
				#print(elem)
				if ": " in str(elem):
					elem = elem.split(": ")[1]
			for elem in column:
				#print(elem)
				if ": " in str(elem):
					elem = elem.split(": ")[1]
					#print(elem)
			for elem in column:
				if(str(elem) == "nan" or elem == np.nan or pd.isna(elem)):
					elem = "NA"
				elif(elem == float('nan')):
					elem = "NA"
				#elif(type(elem) == float):
				#	if(math.isnan(elem)):
				#		elem = "NA"
			for elem in coluniq:
				if(str(elem) == "nan" or elem == np.nan or pd.isna(elem)):
					elem = "NA"
				elif(elem == float('nan')):
					elem = "NA"
				#elif(type(elem) == float):
				#	if(math.isnan(elem)):
				#		elem = "NA"
			#for elem in coluniq:
			#	print(str(elem))
			if "NA" in coluniq:
				if(len(coluniq) == 3):
					#print(coluniq)
					patients_temp_0 = []
					patients_temp_1 = []
					#param_names.append(column_name)
					#param_cols.append(ctr)
					column = column.replace('Good Prognosis','1')
					column = column.replace('Bad Prognosis','0')
					column = column.replace('yes','1')
					column = column.replace('no','0')
					column = column.replace('P','1') 
					column = column.replace('N','0')
					column = column.replace('0A','NA')
					column = column.replace('relapse (event=1; no event=0): 0','0')
					column = column.replace('relapse (event=1; no event=0): 1','0')
					column = column.replace('relapse (event=1; no event=0): na','NA')
					column = column.replace('status: ALIVE','1')
					column = column.replace('status: DEAD','0')
					column = column.replace('status: NTL','NA')
					column = column.replace('ALIVE','1')
					column = column.replace('DEAD','0')
					column = column.replace('NTL','NA')
					#column2 = [str(w).replace(r'^P','1') for w in column]
					#column2 = [re.sub(r'^P$','1',str(w)) for w in column]
					#column2 = [re.sub(r'^N$','0',str(w)) for w in column2]
					#column2 = [re.sub(r'^yes$','1',str(w)) for w in column2]
					#column2 = [re.sub(r'^no$','0',str(w)) for w in column2]
					#column2 = [re.sub(r'^Good Prognosis$','1') for w in column2]
					#column2 = [re.sub(r'^Bad Prognosis$','0') for w in column2]
					#column2 = [str(w).replace('Good Prognosis','1') for w in column]
					#column2 = [str(w).replace('Bad Prognosis','0') for w in column2]
					#column2 = [str(w).replace('yes','1') for w in column2]
					#column2 = [str(w).replace('no','0') for w in column2]
					#column2 = [str(w).replace('P','1') for w in column2]
					#column2 = [str(w).replace('N','0') for w in column2]
					#column2 = [str(w).replace('0A','NA') for w in column2]
					#if(column.unique() 
					coluniq2 = column.unique()
					coluniq3 = [str(w) for w in coluniq2]
					print(column_name)
					#print(patients_temp_0)
					if(sorted(coluniq3) == ['0','1','NA'] or sorted(coluniq3) == ['0','1']):
						print(sorted(coluniq3))
						col_as_list = [str(i) for i in column]
						for i in range(0,len(col_as_list)-1):
							if(col_as_list[i] == '0'):
								patients_temp_0.append(patientids_metadata[i])
							elif(col_as_list[i] == '1'):
								patients_temp_1.append(patientids_metadata[i])
						patients_0.append(patients_temp_0)
						patients_1.append(patients_temp_1)										
						param_names.append(column_name)
						param_cols.append(ctr)
						all_patients = patients_temp_0 + patients_temp_1
						#print(all_patients)
						tmp95 = []
						tmp94 = []
						for i in range(0,len(all_patients)-1):
						#	print(all_patients[i])
							if(all_patients[i] in group1_ids):
								tmp95.append(all_patients[i])
							elif(all_patients[i] in group2_ids):	
								tmp94.append(all_patients[i])
						group1_has.append(tmp95)				
						group2_has.append(tmp94)
						#print(group1_has)									
					#print(column2)
					#param_values.append(column)
			else:
				if(len(coluniq) == 2):
					#print(coluniq)
					column = column.replace('Good Prognosis','1')
					column = column.replace('Bad Prognosis','0')
					column = column.replace('yes','1')
					column = column.replace('no','0')
					column = column.replace('P','1') 
					column = column.replace('N','0')
					column = column.replace('0A','NA')
					column = column.replace('relapse (event=1; no event=0): 0','0')
					column = column.replace('relapse (event=1; no event=0): 1','0')
					column = column.replace('relapse (event=1; no event=0): na','NA')
					column = column.replace('status: ALIVE','1')
					column = column.replace('status: DEAD','0')
					column = column.replace('status: NTL','NA')
					column = column.replace('ALIVE','1')
					column = column.replace('DEAD','0')
					column = column.replace('NTL','NA')
					#column2 = [str(w).replace('Good Prognosis','1') for w in column]
					#column2 = [str(w).replace('Bad Prognosis','0') for w in column2]
					#column2 = [str(w).replace('yes','1') for w in column2]
					#column2 = [str(w).replace('no','0') for w in column2]
					#column2 = [str(w).replace('P','1') for w in column2]
					#column2 = [str(w).replace('N','0') for w in column2]
					#column2 = [str(w).replace('0A','NA') for w in column2]
					#print(column.unique())
					coluniq2 = column.unique()
					coluniq3 = [str(w) for w in coluniq2]
					#print(sorted(coluniq3))
					if(sorted(coluniq3) == ['0','1','NA'] or sorted(coluniq3) == ['0','1']):
						col_as_list = [str(i) for i in column]
						for i in range(0,len(col_as_list)-1):
							if(col_as_list[i] == '0'):
								patients_temp_0.append(patientids_metadata[i])
							elif(col_as_list[i] == '1'):
								patients_temp_1.append(patientids_metadata[i])
						patients_0.append(patients_temp_0)
						patients_1.append(patients_temp_1)										
						param_names.append(column_name)
						param_cols.append(ctr)
						all_patients = patients_temp_0 + patients_temp_1
						#print(all_patients)
						tmp95 = []
						tmp94 = []
						for i in range(0,len(all_patients)-1):
							#print(all_patients[i])
							if(all_patients[i] in group1_ids):
								tmp95.append(all_patients[i])
							elif(all_patients[i] in group2_ids):	
								tmp94.append(all_patients[i])
						group1_has.append(tmp95)				
						group2_has.append(tmp94)
						#print(group1_has)									
								
					#param_values.append(column)
			ctr = ctr + 1
		print(param_names)
		print(patients_0)
		print(patients_1)
		print(param_cols)
		jaccards_1 = []
		jaccards_2 = []
		param_names_final = []
		nbr_patients = float(len(group1_ids) + len(group2_ids))
		for i in range(0,len(param_names)-1):
			if((float(len(patients_0[i])+len(patients_1[i]))/nbr_patients) > 0.8):
				param_names_final.append(param_names[i])
				jaccards_1.append(lib.jac(group1_has[i],patients_0[i]))
				jaccards_2.append(lib.jac(group2_has[i],patients_1[i]))
		print(param_names_final)
		#param_names_final.replace(":ch1","")
		#param_names_final.replace("bm event","Breast metastasis")
		#param_names_final.replace("lm event","Lung metastasis")
		#param_names_final.replace("met event","Metastasis")
		print(group1_has)
		print(jaccards_1)
		print(jaccards_2)
				#jaccards_1.append(lib.jac(patients_0[i],grou
			#col_uniq = clinicaldf.column_name.unique()
			#print(col_uniq)
		print(list(clinicaldf.columns).index(survival_col))
		print(list(clinicaldf.iloc[:,42]))
		if(survival_col in list(clinicaldf.columns)):
			survival_col_nbr = list(clinicaldf.columns).index(survival_col)
			#survival_col_nbr = title_col.index(survival_col) + 1
			print("survival col")
			print(survival_col_nbr)
			print(list(clinicaldf.iloc[:,survival_col_nbr]))
		#clinical_stringio = StringIO(clinicalstr)
		#print("finished stringio")		
		#print(clinicaldf)
		#print(clinicalLines[0])
		patient_id_list = list(clinicaldf.iloc[:,0])
		clinicaldf.iloc[:,survival_col_nbr].fillna("NA",inplace=True)
		survivalcol_list = list(clinicaldf.iloc[:,survival_col_nbr])
		for i in range(0,len(patient_id_list)):
			#bababa = line.replace('"','')
			#bababa2 = line.split("\",\"")
			#bababa2 = line.split(",")
			#print(bababa2)
			#bababa2[0] = bababa2[0].replace('"','')
			#print(bababa2[0])
			#print(bababa2[64])
			#if(len(bababa2) > survival_col_nbr):
				#print(bababa2[53])	
				#print(bababa2[63])
				#print(bababa2[72])
			if(survivalcol_list[i] != "--" and survivalcol_list[i] != "name:ch1" and survivalcol_list[i] != "NA" and survivalcol_list[i].replace('.','',1).isdigit()):
					patientData.update({patient_id_list[i]:survivalcol_list[i]})
					print(survivalcol_list[i])
				#	print(patientData)
			#if(len(bababa2) > 71):
				
		age1 = 0
		ct1 = 0.001
		age2 = 0
		ct2 = 0.001
		print(patientData_misc)
		#for key in patientData_misc:
		#	if key in group1_ids:
		#		age1 = age1 + patientData_misc[key][0]
		#		ct1 = ct1 + 1
		#		size1 = size1 + patientData_misc[key][2]
		#		ct3 = ct3 + 1
		#		ct5 = ct5 + 1
		#		metas1 = metas1 + patientData_misc[key][1]
		#	elif key in group2_ids:
		#		age2 = age2 + patientData_misc[key][0]
		#		ct2 = ct2 + 1
		#		size2 = size2 + patientData_misc[key][2]
		#		ct4 = ct4 + 1
		#		ct6 = ct6 + 1
		#		metas2 = metas2 + patientData_misc[key][1]
		ret_metadata = []
		survival_1 = []
		survival_2 = []
		survival_mean_1 = 0
		survival_mean_2 = 0
		ctr_surv_1 = 0.001
		sum_surv_1 = 0
		ctr_surv_2 = 0.001
		sum_surv_2 = 0
		for key in patientData:
			if key in group1_ids:
				survival_1.append(patientData[key])
			elif key in group2_ids:
				survival_2.append(patientData[key])
		for elem in survival_1:
			sum_surv_1 = sum_surv_1 + float(elem)
			ctr_surv_1 = ctr_surv_1 + 1
		for elem in survival_2:
			sum_surv_2 = sum_surv_2 + float(elem)
			ctr_surv_2 = ctr_surv_2 + 1
		survival_mean_1 = sum_surv_1 / ctr_surv_1
		survival_mean_2 = sum_surv_2 / ctr_surv_2
		#print(group1_in_metadata)
		#print(metastasis) 
		param_names = [elem.replace("bm event:ch1","Breast Metastasis") for elem in param_names]
		param_names = [elem.replace("lm event:ch1","Lung Metastasis") for elem in param_names]
		param_names = [elem.replace("met event:ch1","Metastasis") for elem in param_names]
		param_names = [elem.replace("relapse (event=1; no event=0):ch1","Relapse") for elem in param_names]
		text_file_3 = open("polls/static/metadata.txt", "w")
		#text_file_3.write("<table><tr><th>Patient Cluster</th><th>Correlation with clinical Group</th><th>Mean Survival</th><th>mean age at diagnosis</th><th>mean tumor size</th></tr>")
		#text_file_3.write("<table><tr><th>Patient Cluster</th><th>Lung metastasis</th><th>Breast metastasis</th><th>Lymph nodes</th><th>Metastasis</th><th>Path_ER</th><th>Path_PR</th><th>Good Prognosis</th></tr>")
		text_file_3.write("<table><tr>")
		if(len(jaccards_1) < len(param_names)):
			for i in range(len(jaccards_1),len(param_names)):
				jaccards_1.append(0.0)
				jaccards_2.append(0.0)
		for elem in param_names:
			text_file_3.write("<th>" + str(elem) + "</th>")
		text_file_3.write("</tr><tr>")
		for elem in jaccards_1:
			text_file_3.write("<th>" + str(elem) + "</th>")
		text_file_3.write("</tr><tr>")
		for elem in jaccards_2:
			text_file_3.write("<th>" + str(elem) + "</th>")
		text_file_3.write("</table>")
		text_file_3.close()
		survival_perc_1 = {0:1}
		survival_perc_2 = {0:1}
		for i in range(1,10):	
			tmp1 = 1.0
			tmp2 = 1.0
			for k in survival_1:
				if(float(k) < float(i)):
					tmp1 = tmp1 - (1.0/len(survival_1))
			for k in survival_2:
				if(float(k) < float(i)):
					tmp2 = tmp2 - (1.0/len(survival_2))
			survival_perc_1.update({i:tmp1})
			survival_perc_2.update({i:tmp2})
		#output_notebook()	
		trace1 = go.Scatter(
		x=list(survival_perc_1.keys()),
		y=list(survival_perc_1.values()),
		mode='lines+markers',
		name="'Group 1'",
		hoverinfo='name',
		line=dict(
		shape='hv'))
		trace2 = go.Scatter(
		x=list(survival_perc_2.keys()),
		y=list(survival_perc_2.values()),
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
		#plot_div=plotly.offline.plot(fig, auto_open=False,include_plotlyjs = False, output_type='div')
		plot_div=plotly.offline.plot(fig, auto_open=False,output_type='div')
		with open("polls/static/output_plotly.html", "w") as text_file_2:
        		text_file_2.write(plot_div)
		#plotly.offline.plot(fig,auto_open=False, image_filename="tempimage.png", image='png')
	with open("polls/static/output_console.txt", "w") as text_file:
        	text_file.write("Done!")
	print(ret_metadata)
	return(script,div,plot_1,plot_div,ret_metadata)

@shared_task(name="script_output_task_10")
def script_output_task_10(T,row_colors1,col_colors1,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1_ret,jac_2_ret,survival_col,clinicaldf,session_id):
	#G = nx.Graph()
	#G_2 = nx.Graph()
	print("basdbrasdbawrubaerbsrbjsdb")
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
	json_path = "polls/static/test15_" + session_id + ".json"
	with open(json_path, "w") as text_file:
		text_file.write(jsn3)		
	#nx.write_gexf(G,"test.gexf")
	output_notebook()
	plot = figure(x_range=(-1.5, 1.5), y_range=(-1.5, 1.5))
	# add tools to the plot
	plot.add_tools(HoverTool(tooltips=None),TapTool(),BoxSelectTool())
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
	source = ColumnDataSource(data=dict(x=x, y=y,Name=[node_labels[i] for i in genes.keys()]))
	labels2 = LabelSet(x='x', y='y', text='Name',text_font_size="8pt",x_offset=-20, source=source,background_fill_color='white',background_fill_alpha=0.0,level='glyph',render_mode='canvas')
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
	#plt.savefig("/home/quirin/testproject/polls/static/test.png")
	path99 = "/home/quirin/testproject/polls/static/test_" + session_id + ".png"
	plt.savefig(path99)
	script, div = components(plot)	
	plot_1=plt.gcf()
	plt.clf()
	print("in method")
	patientData = {}
	patientData_misc = {}
	group1_in_metadata = []
	group2_in_metadata = []
	metastasis = []
	no_metastasis = []
	group1_has_bm = []
	group2_has_bm = []
	bm = []
	no_bm = []
	group1_has_lm = []
	group2_has_lm = []
	lm = []
	no_lm = []
	group1_has_lymph = []
	group2_has_lymph = []
	lymph = []
	no_lymph = []
	group1_has_path_er = []
	group2_has_path_er = []
	path_er = []
	no_path_er = []
	group1_has_path_pr = []
	group2_has_path_pr = []
	path_pr = []
	no_path_pr = []
	group1_has_prognosis = []
	group2_has_prognosis = []
	prognosis_good = []
	prognosis_bad = []
	path_genelist = "genelist_" + session_id + ".txt"
	path_genelist_1 = "genelist_1_" + session_id + ".txt"
	path_genelist_2 = "genelist_2_" + session_id + ".txt"
	with open(path_genelist,"w") as text_file_4:
		for i in G_list:
			text_file_4.write(str(i) + "\n")
	text_file_4.close()
	with open(path_genelist_1,"w") as text_file_5:
		for i in G_list:
			if(i in genes1):
				text_file_5.write(str(i) + "\n")
	text_file_5.close()
	with open(path_genelist_2,"w") as text_file_6:
		for i in G_list:
			if(i not in genes1):
				text_file_6.write(str(i) + "\n")
	text_file_6.close()
	path_metadata = "polls/static/metadata_" + session_id + ".txt"
	if(clinicalstr == "empty"):
		ret_metadata = []
		ret_metadata.append({'group':"Group 1",'lm':"NA",'bm':"NA",'lymph':"NA",'metastasis':"NA",'path_er':"NA",'path_pr':"NA",'good_prognosis':"NA"})
		ret_metadata.append({'group':"Group 2",'lm':"NA",'bm':"NA",'lymph':"NA",'metastasis':"NA",'path_er':"NA",'path_pr':"NA",'good_prognosis':"NA"})
		ret_metadata.append({'group':"Group 1",'jac':jac_1_ret,'survival':"NA",'age':"NA",'size':"NA"})
		ret_metadata.append({'group':"Group 2",'jac':jac_2_ret,'survival':"NA",'age':"NA",'size':"NA"})
		text_file_3 = open(path_metadata, "w")
		text_file_3.write("<table><tr><th>Patient Cluster</th></tr>")
		text_file_3.write("<tr><th>Group1</th></tr>")
		text_file_3.write("</table>")
		#text_file_3.write("<table><tr><th>Patient Cluster</th><th>Lung metastasis</th><th>Breast metastasis</th><th>Lymph nodes</th><th>Metastasis</th><th>Path_ER</th><th>Path_PR</th><th>Good Prognosis</th></tr>")
		#for elem in ret_metadata:
		#	text_file_3.write("<tr>")
		#	text_file_3.write("<th>" + str(elem['lm']) + "</th>")
		#	text_file_3.write("<th>" + str(elem['bm']) + "</th>")
		#	text_file_3.write("<th>" + str(elem['lymph']) + "</th>")
		#	text_file_3.write("<th>" + str(elem['metastasis']) + "</th>")
		#	text_file_3.write("<th>" + str(elem['path_er']) + "</th>")
		#	text_file_3.write("<th>" + str(elem['path_pr']) + "</th>")
		#	text_file_3.write("<th>" + str(elem['good_prognosis']) + "</th>")
		#	text_file_3.write("</tr>")
		#text_file_3.write("</table>")
		text_file_3.close()	
		plot_div = ""	
		output_plot_path = "polls/static/output_plotly_" + session_id + ".html"
		with open(output_plot_path, "w") as text_file_2:
        		text_file_2.write("")
	
	else:
		print("barfbasdbhasdb")
		clinicalLines = clinicalstr.split("\n")
		title_col = clinicalLines[0].split(",")
		print(title_col)
		survival_col_nbr = 64
		clinicaldf.replace(['NaN','nan','?','--'],['NA','NA','NA','NA'], inplace=True)
		clinicaldf.replace(['NTL'],['NA'], inplace=True)
		clinicaldf.replace(['na'],['NA'], inplace=True)
		#print(clinicaldf.index)
		print(clinicaldf.columns)
		clinicaldf_col_names = list(clinicaldf.columns)
		print(clinicaldf_col_names)
		patientids_metadata = clinicaldf.iloc[:,0].values.tolist()
		if("Unnamed" in "\t".join(list(clinicaldf.columns))):
			print("basdbasdfbasfbafs")
			clinicaldf_col_names_temp = ['empty']
			clinicaldf_col_names_new = clinicaldf_col_names_temp + clinicaldf_col_names
			clinicaldf.columns = list(clinicaldf_col_names_new[:-1])
			print(clinicaldf_col_names_new)
		if("GSM" not in patientids_metadata[0]):
			patientids_metadata = list(clinicaldf.index)
		#print(clinicaldf["relapse (event=1; no event=0):ch1"])
		param_names = []
		param_values = []
		param_cols = []
		ctr = 0
		print(patientids_metadata)
		patients_0 = []
		patients_1 = []
		group1_has = []
		group2_has = []
		for column_name, column in clinicaldf.transpose().iterrows():
			column.fillna("NA",inplace=True)
			#print(column_name)
			#print(column)
			coluniq = column.unique()
			#print(coluniq)
			#print(column_name)
			#print(list(column))
			for elem in coluniq:
				#print(elem)
				if ": " in str(elem):
					elem = elem.split(": ")[1]
			for elem in column:
				#print(elem)
				if ": " in str(elem):
					elem = elem.split(": ")[1]
					#print(elem)
			for elem in column:
				if(str(elem) == "nan" or elem == np.nan or pd.isna(elem)):
					elem = "NA"
				elif(elem == float('nan')):
					elem = "NA"
				#elif(type(elem) == float):
				#	if(math.isnan(elem)):
				#		elem = "NA"
			for elem in coluniq:
				if(str(elem) == "nan" or elem == np.nan or pd.isna(elem)):
					elem = "NA"
				elif(elem == float('nan')):
					elem = "NA"
				#elif(type(elem) == float):
				#	if(math.isnan(elem)):
				#		elem = "NA"
			#for elem in coluniq:
			#	print(str(elem))
			if "NA" in coluniq:
				if(len(coluniq) == 3):
					#print(coluniq)
					patients_temp_0 = []
					patients_temp_1 = []
					#param_names.append(column_name)
					#param_cols.append(ctr)
					column = column.replace('Good Prognosis','1')
					column = column.replace('Bad Prognosis','0')
					column = column.replace('yes','1')
					column = column.replace('no','0')
					column = column.replace('P','1') 
					column = column.replace('N','0')
					column = column.replace('0A','NA')
					column = column.replace('relapse (event=1; no event=0): 0','0')
					column = column.replace('relapse (event=1; no event=0): 1','0')
					column = column.replace('relapse (event=1; no event=0): na','NA')
					column = column.replace('status: ALIVE','1')
					column = column.replace('status: DEAD','0')
					column = column.replace('status: NTL','NA')
					column = column.replace('ALIVE','1')
					column = column.replace('DEAD','0')
					column = column.replace('NTL','NA')
					#column2 = [str(w).replace(r'^P','1') for w in column]
					#column2 = [re.sub(r'^P$','1',str(w)) for w in column]
					#column2 = [re.sub(r'^N$','0',str(w)) for w in column2]
					#column2 = [re.sub(r'^yes$','1',str(w)) for w in column2]
					#column2 = [re.sub(r'^no$','0',str(w)) for w in column2]
					#column2 = [re.sub(r'^Good Prognosis$','1') for w in column2]
					#column2 = [re.sub(r'^Bad Prognosis$','0') for w in column2]
					#column2 = [str(w).replace('Good Prognosis','1') for w in column]
					#column2 = [str(w).replace('Bad Prognosis','0') for w in column2]
					#column2 = [str(w).replace('yes','1') for w in column2]
					#column2 = [str(w).replace('no','0') for w in column2]
					#column2 = [str(w).replace('P','1') for w in column2]
					#column2 = [str(w).replace('N','0') for w in column2]
					#column2 = [str(w).replace('0A','NA') for w in column2]
					#if(column.unique() 
					coluniq2 = column.unique()
					coluniq3 = [str(w) for w in coluniq2]
					print(column_name)
					#print(patients_temp_0)
					if(sorted(coluniq3) == ['0','1','NA'] or sorted(coluniq3) == ['0','1']):
						print(sorted(coluniq3))
						col_as_list = [str(i) for i in column]
						for i in range(0,len(col_as_list)-1):
							if(col_as_list[i] == '0'):
								patients_temp_0.append(patientids_metadata[i])
							elif(col_as_list[i] == '1'):
								patients_temp_1.append(patientids_metadata[i])
						patients_0.append(patients_temp_0)
						patients_1.append(patients_temp_1)										
						param_names.append(column_name)
						param_cols.append(ctr)
						all_patients = patients_temp_0 + patients_temp_1
						#print(all_patients)
						tmp95 = []
						tmp94 = []
						for i in range(0,len(all_patients)-1):
						#	print(all_patients[i])
							if(all_patients[i] in group1_ids):
								tmp95.append(all_patients[i])
							elif(all_patients[i] in group2_ids):	
								tmp94.append(all_patients[i])
						group1_has.append(tmp95)				
						group2_has.append(tmp94)
						#print(group1_has)									
					#print(column2)
					#param_values.append(column)
			else:
				if(len(coluniq) == 2):
					#print(coluniq)
					column = column.replace('Good Prognosis','1')
					column = column.replace('Bad Prognosis','0')
					column = column.replace('yes','1')
					column = column.replace('no','0')
					column = column.replace('P','1') 
					column = column.replace('N','0')
					column = column.replace('0A','NA')
					column = column.replace('relapse (event=1; no event=0): 0','0')
					column = column.replace('relapse (event=1; no event=0): 1','0')
					column = column.replace('relapse (event=1; no event=0): na','NA')
					column = column.replace('status: ALIVE','1')
					column = column.replace('status: DEAD','0')
					column = column.replace('status: NTL','NA')
					column = column.replace('ALIVE','1')
					column = column.replace('DEAD','0')
					column = column.replace('NTL','NA')
					#column2 = [str(w).replace('Good Prognosis','1') for w in column]
					#column2 = [str(w).replace('Bad Prognosis','0') for w in column2]
					#column2 = [str(w).replace('yes','1') for w in column2]
					#column2 = [str(w).replace('no','0') for w in column2]
					#column2 = [str(w).replace('P','1') for w in column2]
					#column2 = [str(w).replace('N','0') for w in column2]
					#column2 = [str(w).replace('0A','NA') for w in column2]
					#print(column.unique())
					coluniq2 = column.unique()
					coluniq3 = [str(w) for w in coluniq2]
					#print(sorted(coluniq3))
					if(sorted(coluniq3) == ['0','1','NA'] or sorted(coluniq3) == ['0','1']):
						col_as_list = [str(i) for i in column]
						for i in range(0,len(col_as_list)-1):
							if(col_as_list[i] == '0'):
								patients_temp_0.append(patientids_metadata[i])
							elif(col_as_list[i] == '1'):
								patients_temp_1.append(patientids_metadata[i])
						patients_0.append(patients_temp_0)
						patients_1.append(patients_temp_1)										
						param_names.append(column_name)
						param_cols.append(ctr)
						all_patients = patients_temp_0 + patients_temp_1
						#print(all_patients)
						tmp95 = []
						tmp94 = []
						for i in range(0,len(all_patients)-1):
							#print(all_patients[i])
							if(all_patients[i] in group1_ids):
								tmp95.append(all_patients[i])
							elif(all_patients[i] in group2_ids):	
								tmp94.append(all_patients[i])
						group1_has.append(tmp95)				
						group2_has.append(tmp94)
						#print(group1_has)									
								
					#param_values.append(column)
			ctr = ctr + 1
		print(param_names)
		print(patients_0)
		print(patients_1)
		print(param_cols)
		jaccards_1 = []
		jaccards_2 = []
		param_names_final = []
		nbr_patients = float(len(group1_ids) + len(group2_ids))
		for i in range(0,len(param_names)-1):
			if((float(len(patients_0[i])+len(patients_1[i]))/nbr_patients) > 0.8):
				param_names_final.append(param_names[i])
				jaccards_1.append(lib.jac(group1_has[i],patients_0[i]))
				jaccards_2.append(lib.jac(group2_has[i],patients_1[i]))
		print(param_names_final)
		#param_names_final.replace(":ch1","")
		#param_names_final.replace("bm event","Breast metastasis")
		#param_names_final.replace("lm event","Lung metastasis")
		#param_names_final.replace("met event","Metastasis")
		print(group1_has)
		print(jaccards_1)
		print(jaccards_2)
				#jaccards_1.append(lib.jac(patients_0[i],grou
			#col_uniq = clinicaldf.column_name.unique()
			#print(col_uniq)
		print(list(clinicaldf.columns).index(survival_col))
		print(list(clinicaldf.iloc[:,42]))
		if(survival_col in list(clinicaldf.columns)):
			survival_col_nbr = list(clinicaldf.columns).index(survival_col)
			#survival_col_nbr = title_col.index(survival_col) + 1
			print("survival col")
			print(survival_col_nbr)
			print(list(clinicaldf.iloc[:,survival_col_nbr]))
		#clinical_stringio = StringIO(clinicalstr)
		#print("finished stringio")		
		#print(clinicaldf)
		#print(clinicalLines[0])
		patient_id_list = list(clinicaldf.iloc[:,0])
		clinicaldf.iloc[:,survival_col_nbr].fillna("NA",inplace=True)
		survivalcol_list = list(clinicaldf.iloc[:,survival_col_nbr])
		for i in range(0,len(patient_id_list)):
			#bababa = line.replace('"','')
			#bababa2 = line.split("\",\"")
			#bababa2 = line.split(",")
			#print(bababa2)
			#bababa2[0] = bababa2[0].replace('"','')
			#print(bababa2[0])
			#print(bababa2[64])
			#if(len(bababa2) > survival_col_nbr):
				#print(bababa2[53])	
				#print(bababa2[63])
				#print(bababa2[72])
			if(survivalcol_list[i] != "--" and survivalcol_list[i] != "name:ch1" and survivalcol_list[i] != "NA" and survivalcol_list[i].replace('.','',1).isdigit()):
					patientData.update({patient_id_list[i]:survivalcol_list[i]})
					print(survivalcol_list[i])
				#	print(patientData)
			#if(len(bababa2) > 71):
				
		age1 = 0
		ct1 = 0.001
		age2 = 0
		ct2 = 0.001
		print(patientData_misc)
		#for key in patientData_misc:
		#	if key in group1_ids:
		#		age1 = age1 + patientData_misc[key][0]
		#		ct1 = ct1 + 1
		#		size1 = size1 + patientData_misc[key][2]
		#		ct3 = ct3 + 1
		#		ct5 = ct5 + 1
		#		metas1 = metas1 + patientData_misc[key][1]
		#	elif key in group2_ids:
		#		age2 = age2 + patientData_misc[key][0]
		#		ct2 = ct2 + 1
		#		size2 = size2 + patientData_misc[key][2]
		#		ct4 = ct4 + 1
		#		ct6 = ct6 + 1
		#		metas2 = metas2 + patientData_misc[key][1]
		ret_metadata = []
		survival_1 = []
		survival_2 = []
		survival_mean_1 = 0
		survival_mean_2 = 0
		ctr_surv_1 = 0.001
		sum_surv_1 = 0
		ctr_surv_2 = 0.001
		sum_surv_2 = 0
		for key in patientData:
			if key in group1_ids:
				survival_1.append(patientData[key])
			elif key in group2_ids:
				survival_2.append(patientData[key])
		for elem in survival_1:
			sum_surv_1 = sum_surv_1 + float(elem)
			ctr_surv_1 = ctr_surv_1 + 1
		for elem in survival_2:
			sum_surv_2 = sum_surv_2 + float(elem)
			ctr_surv_2 = ctr_surv_2 + 1
		survival_mean_1 = sum_surv_1 / ctr_surv_1
		survival_mean_2 = sum_surv_2 / ctr_surv_2
		#print(group1_in_metadata)
		#print(metastasis) 
		param_names = [elem.replace("bm event:ch1","Breast Metastasis") for elem in param_names]
		param_names = [elem.replace("lm event:ch1","Lung Metastasis") for elem in param_names]
		param_names = [elem.replace("met event:ch1","Metastasis") for elem in param_names]
		param_names = [elem.replace("relapse (event=1; no event=0):ch1","Relapse") for elem in param_names]
		text_file_3 = open(path_metadata, "w")
		#text_file_3.write("<table><tr><th>Patient Cluster</th><th>Correlation with clinical Group</th><th>Mean Survival</th><th>mean age at diagnosis</th><th>mean tumor size</th></tr>")
		#text_file_3.write("<table><tr><th>Patient Cluster</th><th>Lung metastasis</th><th>Breast metastasis</th><th>Lymph nodes</th><th>Metastasis</th><th>Path_ER</th><th>Path_PR</th><th>Good Prognosis</th></tr>")
		text_file_3.write("<table><tr>")
		if(len(jaccards_1) < len(param_names)):
			for i in range(len(jaccards_1),len(param_names)):
				jaccards_1.append(0.0)
				jaccards_2.append(0.0)
		for elem in param_names:
			text_file_3.write("<th>" + str(elem) + "</th>")
		text_file_3.write("</tr><tr>")
		for elem in jaccards_1:
			text_file_3.write("<th>" + str(elem) + "</th>")
		text_file_3.write("</tr><tr>")
		for elem in jaccards_2:
			text_file_3.write("<th>" + str(elem) + "</th>")
		text_file_3.write("</table>")
		text_file_3.close()
		survival_perc_1 = {0:1}
		survival_perc_2 = {0:1}
		for i in range(1,10):	
			tmp1 = 1.0
			tmp2 = 1.0
			for k in survival_1:
				if(float(k) < float(i)):
					tmp1 = tmp1 - (1.0/len(survival_1))
			for k in survival_2:
				if(float(k) < float(i)):
					tmp2 = tmp2 - (1.0/len(survival_2))
			survival_perc_1.update({i:tmp1})
			survival_perc_2.update({i:tmp2})
		#output_notebook()	
		trace1 = go.Scatter(
		x=list(survival_perc_1.keys()),
		y=list(survival_perc_1.values()),
		mode='lines+markers',
		name="'Group 1'",
		hoverinfo='name',
		line=dict(
		shape='hv'))
		trace2 = go.Scatter(
		x=list(survival_perc_2.keys()),
		y=list(survival_perc_2.values()),
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
		#plot_div=plotly.offline.plot(fig, auto_open=False,include_plotlyjs = False, output_type='div')
		plot_div=plotly.offline.plot(fig, auto_open=False,output_type='div')
		output_plot_path = "polls/static/output_plotly_" + session_id + ".html"
		with open(output_plot_path, "w") as text_file_2:
        		text_file_2.write(plot_div)
		#plotly.offline.plot(fig,auto_open=False, image_filename="tempimage.png", image='png')
	with open("polls/static/output_console.txt", "w") as text_file:
        	text_file.write("Done!")
	print(ret_metadata)
	return(script,div,plot_1,plot_div,ret_metadata,path99,path_metadata,output_plot_path,json_path)




@shared_task(name="script_output_task_7")
def script_output_task_7(T,row_colors1,col_colors1,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1_ret,jac_2_ret,metadata_is_lungcancer,survival_col):
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
	plot.add_tools(HoverTool(tooltips=None),TapTool(),BoxSelectTool())
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
	source = ColumnDataSource(data=dict(x=x, y=y,Name=[node_labels[i] for i in genes.keys()]))
	labels2 = LabelSet(x='x', y='y', text='Name',text_font_size="8pt",x_offset=-20, source=source,background_fill_color='white',background_fill_alpha=0.0,level='glyph',render_mode='canvas')
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
	patientData = {}
	patientData_misc = {}
	
	if(clinicalstr == "empty"):
		ret_metadata = []
		ret_metadata.append({'group':"Group 1",'lm':"NA",'bm':"NA",'lymph':"NA",'metastasis':"NA",'path_er':"NA",'path_pr':"NA",'good_prognosis':"NA"})
		ret_metadata.append({'group':"Group 2",'lm':"NA",'bm':"NA",'lymph':"NA",'metastasis':"NA",'path_er':"NA",'path_pr':"NA",'good_prognosis':"NA"})
		#ret_metadata.append({'group':"Group 1",'jac':jac_1_ret,'survival':"NA",'age':"NA",'size':"NA"})
		#ret_metadata.append({'group':"Group 2",'jac':jac_2_ret,'survival':"NA",'age':"NA",'size':"NA"})
		text_file_3 = open("polls/static/metadata.txt", "w")
		text_file_3.write("<table><tr><th>Patient Cluster</th><th>Lung metastasis</th><th>Breast metastasis</th><th>Lymph nodes</th><th>Metastasis</th><th>Path_ER</th><th>Path_PR</th><th>Good Prognosis</th></tr>")
		for elem in ret_metadata:
			text_file_3.write("<tr>")
			text_file_3.write("<th>" + str(elem['lm']) + "</th>")
			text_file_3.write("<th>" + str(elem['bm']) + "</th>")
			text_file_3.write("<th>" + str(elem['lymph']) + "</th>")
			text_file_3.write("<th>" + str(elem['metastasis']) + "</th>")
			text_file_3.write("<th>" + str(elem['path_er']) + "</th>")
			text_file_3.write("<th>" + str(elem['path_pr']) + "</th>")
			text_file_3.write("<th>" + str(elem['good_prognosis']) + "</th>")
			text_file_3.write("</tr>")
		text_file_3.write("</table>")
		text_file_3.close()	
		plot_div = ""
		with open("polls/static/output_plotly.html", "w") as text_file_2:
        		text_file_2.write("")
	
	elif(metadata_is_lungcancer == "false"):
		ret_metadata = []
		nbr_patients = float(len(group1_ids) + len(group2_ids))
		group1_in_metadata = []
		group2_in_metadata = []
		metastasis = []
		no_metastasis = []
		group1_has_bm = []
		group2_has_bm = []
		bm = []
		no_bm = []
		group1_has_lm = []
		group2_has_lm = []
		lm = []
		no_lm = []
		group1_has_lymph = []
		group2_has_lymph = []
		lymph = []
		no_lymph = []
		group1_has_path_er = []
		group2_has_path_er = []
		path_er = []
		no_path_er = []
		group1_has_path_pr = []
		group2_has_path_pr = []
		path_pr = []
		no_path_pr = []
		group1_has_prognosis = []
		group2_has_prognosis = []
		prognosis_good = []
		prognosis_bad = []
		clinicalLines = clinicalstr.split("\n")
		title_col = clinicalLines[0].split(",")
		print(title_col)
		survival_col_nbr = 64
		if(survival_col in title_col):
			survival_col_nbr = title_col.index(survival_col) + 1
			print("survival col")
		for line in clinicalLines:
			#bababa = line.replace('"','')
			bababa2 = line.split("\",\"")
			bababa2[0] = bababa2[0].replace('"','')
			if(len(bababa2) > survival_col_nbr):
				#print(bababa2[53])	
				#print(bababa2[63])
				#print(bababa2[72])
				if(bababa2[survival_col_nbr] != "--" and bababa2[survival_col_nbr] != "name:ch1" and bababa2[survival_col_nbr] != "NA" and bababa2[survival_col_nbr].replace('.','',1).isdigit()):
					patientData.update({bababa2[0]:bababa2[survival_col_nbr]})
				#	print(patientData)
			#if(len(bababa2) > 71):
			if(len(bababa2) > 63):
				if(bababa2[63].isdigit()):
					print("patient id" + str(bababa2[0]))
					if(bababa2[0] in group1_ids):
						print("foo")
						group1_in_metadata.append(bababa2[0])
					elif(bababa2[0] in group2_ids):
						group2_in_metadata.append(bababa2[0])
					if(int(bababa2[63]) == 1):
						metastasis.append(bababa2[0])
					elif(int(bababa2[63]) == 0):
						no_metastasis.append(bababa2[0])
					tmp97 = [0.0,float(bababa2[63]),0.0]
					#patientData_misc.update({bababa2[0]:tmp97})
				if(bababa2[61].isdigit()):
					if(bababa2[0] in group1_ids):
						group1_has_lm.append(bababa2[0])
					elif(bababa2[0] in group2_ids):
						group1_has_lm.append(bababa2[0])
					if(int(bababa2[61]) == 1):
						lm.append(bababa2[0])
					elif(int(bababa2[61]) == 0):
						no_lm.append(bababa2[0])
				if(bababa2[58].isdigit()):
					if(bababa2[0] in group1_ids):
						group1_has_bm.append(bababa2[0])
					elif(bababa2[0] in group2_ids):
						group1_has_bm.append(bababa2[0])
					if(int(bababa2[58]) == 1):
						bm.append(bababa2[0])
					elif(int(bababa2[58]) == 0):
						no_bm.append(bababa2[0])
				if(len(bababa2[15]) > 1 and len(bababa2[15]) < 20):
					if(bababa2[15].split(": ")[1] == "N"):
						if(bababa2[0] in group1_ids):
							group1_has_path_er.append(bababa2[0])
						elif(bababa2[0] in group2_ids):
							group2_has_path_er.append(bababa2[0])
						no_path_er.append(bababa2[0])
					if(bababa2[15].split(": ")[1] == "P"):
						if(bababa2[0] in group1_ids):
							group1_has_path_er.append(bababa2[0])
						elif(bababa2[0] in group2_ids):
							group2_has_path_er.append(bababa2[0])
						path_er.append(bababa2[0])
				if(len(bababa2[16]) > 1 and len(bababa2[16]) < 20):
					if(bababa2[16].split(": ")[1] == "N"):
						if(bababa2[0] in group1_ids):
							group1_has_path_pr.append(bababa2[0])
						elif(bababa2[0] in group2_ids):
							group2_has_path_pr.append(bababa2[0])
						no_path_pr.append(bababa2[0])
					if(bababa2[16].split(": ")[1] == "P"):
						if(bababa2[0] in group1_ids):
							group1_has_path_pr.append(bababa2[0])
						elif(bababa2[0] in group2_ids):
							group2_has_path_pr.append(bababa2[0])
						path_pr.append(bababa2[0])
				if(len(bababa2) > 71):	
					if(bababa2[0] in group1_ids):
						group1_has_prognosis.append(bababa2[0])
					elif(bababa2[0] in group2_ids):
						group2_has_prognosis.append(bababa2[0])
					if(bababa2[71] == "Good Prognosis"):
						prognosis_good.append(bababa2[0])
					elif(bababa2[71] == "Bad Prognosis"):
						prognosis_bad.append(bababa2[0])
				if(len(bababa2) > 70):				
					if(bababa2[70].isdigit()):
						if(bababa2[0] in group1_ids):
							group1_has_lymph.append(bababa2[0])
						elif(bababa2[0] in group2_ids):
							group1_has_lymph.append(bababa2[0])
						if(int(bababa2[70]) == 1):
							lymph.append(bababa2[0])
						elif(int(bababa2[70]) == 0):
							no_lymph.append(bababa2[0])
				elif(bababa2[53].replace('.','',1).isdigit() and bababa2[63].isdigit()):
					print("patient id" + str(bababa2[0]))
					if(bababa2[0] in group1_ids):
						print("foo")
						group1_in_metadata.append(bababa2[0])
					elif(bababa2[0] in group2_ids):
						group2_in_metadata.append(bababa2[0])
					if(int(bababa2[63]) == 1):
						metastasis.append(bababa2[0])
					elif(int(bababa2[63]) == 0):
						no_metastasis.append(bababa2[0])
					#tmp97 = [float(bababa2[53]),float(bababa2[63]),0.0]
					#patientData_misc.update({bababa2[0]:tmp97})
				
		age1 = 0
		ct1 = 0.001
		age2 = 0
		ct2 = 0.001
		print(patientData_misc)
		
		#print(group1_in_metadata)
		#print(metastasis)
		if(float(len(group1_in_metadata) + len(group2_in_metadata)) / nbr_patients > 0.8):
			metas_jac = lib.jac(group1_in_metadata,metastasis)
			metas_jac_2 = lib.jac(group2_in_metadata,metastasis)
		else:
			metas_jac = "NA"
			metas_jac_2 = "NA"
		if(float(len(group1_has_lm) + len(group2_has_lm)) / nbr_patients > 0.8):
			lm_jac = lib.jac(group1_has_lm,lm)
			lm_jac_2 = lib.jac(group2_has_lm,lm)
		else:
			lm_jac = "NA"
			lm_jac_2 = "NA"
		if(float(len(group1_has_bm) + len(group2_has_bm)) / nbr_patients > 0.8):
			bm_jac = lib.jac(group1_has_bm,bm)
			bm_jac_2 = lib.jac(group2_has_bm,bm)
		else:
			bm_jac = "NA"
			bm_jac_2 = "NA"
		if(float(len(group1_has_lymph) + len(group2_has_lymph)) / nbr_patients > 0.8):
			lymph_jac = lib.jac(group1_has_lymph,lymph)
			lymph_jac_2 = lib.jac(group2_has_lymph,lymph)
		else:
			lymph_jac = "NA"
			lymph_jac_2 = "NA"
		if(float(len(group1_has_path_er) + len(group2_has_path_er)) / nbr_patients > 0.8):
			path_er_jac = lib.jac(group1_has_path_er,path_er)
			path_er_jac_2 = lib.jac(group2_has_path_er,path_er)
		else:
			path_er_jac = "NA"
			path_er_jac_2 = "NA"
		if(float(len(group1_has_path_er) + len(group2_has_path_er)) / nbr_patients > 0.8):
			path_pr_jac = lib.jac(group1_has_path_pr,path_pr)
			path_pr_jac_2 = lib.jac(group2_has_path_pr,path_pr)
		else:
			path_pr_jac = "NA"
			path_pr_jac_2 = "NA"
		if(float(len(group1_has_prognosis) + len(group2_has_prognosis)) / nbr_patients > 0.8):
			prognosis_jac = lib.jac(group1_has_prognosis,prognosis_good)
			prognosis_jac_2 = lib.jac(group2_has_prognosis,prognosis_good)
		else:
			prognosis_jac = "NA"
			prognosis_jac_2 = "NA"
		ret_metadata.append({'group':"Group 1",'lm':lm_jac,'bm':bm_jac,'lymph':lymph_jac,'metastasis':metas_jac,'path_er':path_er_jac,'path_pr':path_pr_jac,'good_prognosis':prognosis_jac})
		ret_metadata.append({'group':"Group 2",'lm':lm_jac_2,'bm':bm_jac_2,'lymph':lymph_jac_2,'metastasis':metas_jac_2,'path_er':path_er_jac_2,'path_pr':path_pr_jac_2,'good_prognosis':prognosis_jac_2})
		#ret_metadata.append({'group':"Group 1",'jac':jac_1_ret,'survival':survival_mean_1,'age':mean_age_1,'metastasis':metas_mean_1,'size':mean_size_1})
		#ret_metadata.append({'group':"Group 2",'jac':jac_2_ret,'survival':survival_mean_2,'age':mean_age_2,'metastasis':metas_mean_2,'size':mean_size_2})
		#print(survival_mean_1)
		text_file_3 = open("polls/static/metadata.txt", "w")
		#text_file_3.write("<table><tr><th>Patient Cluster</th><th>Correlation with clinical Group</th><th>Mean Survival</th><th>mean age at diagnosis</th><th>mean tumor size</th></tr>")
		text_file_3.write("<table><tr><th>Patient Cluster</th><th>Lung metastasis</th><th>Breast metastasis</th><th>Lymph nodes</th><th>Metastasis</th><th>Path_ER</th><th>Path_PR</th><th>Good Prognosis</th></tr>")
		for elem in ret_metadata:				
			text_file_3.write("<tr>")
			text_file_3.write("<th>" + str(elem['lm']) + "</th>")
			text_file_3.write("<th>" + str(elem['bm']) + "</th>")
			text_file_3.write("<th>" + str(elem['lymph']) + "</th>")
			text_file_3.write("<th>" + str(elem['metastasis']) + "</th>")
			text_file_3.write("<th>" + str(elem['path_er']) + "</th>")
			text_file_3.write("<th>" + str(elem['path_pr']) + "</th>")
			text_file_3.write("<th>" + str(elem['good_prognosis']) + "</th>")
			text_file_3.write("</tr>")
		text_file_3.write("</table>")
		text_file_3.close()
	else:
		ret_metadata = []
		nbr_patients = float(len(group1_ids) + len(group2_ids))
		group1_has_gender = []
		group2_has_gender = []
		gender_male = []
		gender_female = []
		group1_has_living = []
		group2_has_living = []
		living = []
		no_living = []
		group1_has_relapse = []
		group2_has_relapse = []
		relapse = []
		no_relapse = []
		clinicalLines = clinicalstr.split("\n")
		title_col = clinicalLines[0].split(",")
		print(title_col)
		survival_col_nbr = 42
		if(survival_col in title_col):
			survival_col_nbr = title_col.index(survival_col) + 1
			print("survival col")
		for line in clinicalLines:
			#bababa = line.replace('"','')
			bababa2 = line.split("\",\"")
			bababa2[0] = bababa2[0].replace('"','')
			if(len(bababa2) > survival_col_nbr):
				#print(bababa2[53])	
				#print(bababa2[63])
				#print(bababa2[72])
				if(bababa2[survival_col_nbr] != "--" and bababa2[survival_col_nbr] != "name:ch1" and bababa2[survival_col_nbr] != "NA" and bababa2[survival_col_nbr].replace('.','',1).isdigit()):
					surv_time = float(bababa2[survival_col_nbr])/12				
					patientData.update({bababa2[0]:surv_time})
				#	print(patientData)
			#if(len(bababa2) > 71):
			if(len(bababa2) > 51):
				if(bababa2[44] == "F" or bababa2[44] == "M"):
					if(bababa2[0] in group1_ids):
						print("foo")
						group1_has_gender.append(bababa2[0])
					elif(bababa2[0] in group2_ids):
						group2_has_gender.append(bababa2[0])
					if(bababa2[44] == "F"):
						gender_female.append(bababa2[0])
					elif(bababa2[44] == "M"):
						gender_male.append(bababa2[0])
					#tmp97 = [0.0,float(bababa2[63]),0.0]
					#patientData_misc.update({bababa2[0]:tmp97})
				if(bababa2[50].isdigit()):
					if(bababa2[0] in group1_ids):
						group1_has_relapse.append(bababa2[0])
					elif(bababa2[0] in group2_ids):
						group2_has_relapse.append(bababa2[0])
					if(int(bababa2[50]) == 1):
						relapse.append(bababa2[0])
					elif(int(bababa2[50]) == 0):
						no_relapse.append(bababa2[0])
				if(bababa2[51] == "ALIVE" or bababa2[51] == "DEAD"):
					if(bababa2[0] in group1_ids):
						print("foo")
						group1_has_living.append(bababa2[0])
					elif(bababa2[0] in group2_ids):
						group2_has_living.append(bababa2[0])
					if(bababa2[51] == "ALIVE"):
						living.append(bababa2[0])
					elif(bababa2[51] == "DEAD"):
						no_living.append(bababa2[0])
					#tmp97 = [0.0,float(bababa2[63]),0.0]
					#patientData_misc.update({bababa2[0]:tmp97})
				
				
		age1 = 0
		ct1 = 0.001
		age2 = 0
		ct2 = 0.001
		print(patientData_misc)
		
		#print(group1_in_metadata)
		#print(metastasis)
		if(float(len(group1_has_gender) + len(group2_has_gender)) / nbr_patients > 0.8):
			gender_jac = lib.jac(group1_has_gender,gender_female)
			gender_jac_2 = lib.jac(group2_has_gender,gender_female)
		else:
			gender_jac = "NA"
			gender_jac_2 = "NA"
		if(float(len(group1_has_relapse) + len(group2_has_relapse)) / nbr_patients > 0.8):
			relapse_jac = lib.jac(group1_has_relapse,relapse)
			relapse_jac_2 = lib.jac(group2_has_relapse,relapse)
		else:
			relapse_jac = "NA"
			relapse_jac_2 = "NA"
		if(float(len(group1_has_living) + len(group2_has_living)) / nbr_patients > 0.8):
			living_jac = lib.jac(group1_has_living,living)
			living_jac_2 = lib.jac(group2_has_living,living)
		else:
			living_jac = "NA"
			living_jac_2 = "NA"
		ret_metadata.append({'group':"Group 1",'gender':gender_jac,'relapse':relapse_jac,'living':living_jac})
		ret_metadata.append({'group':"Group 2",'gender':gender_jac_2,'relapse':relapse_jac_2,'living':living_jac_2})
		#ret_metadata.append({'group':"Group 1",'jac':jac_1_ret,'survival':survival_mean_1,'age':mean_age_1,'metastasis':metas_mean_1,'size':mean_size_1})
		#ret_metadata.append({'group':"Group 2",'jac':jac_2_ret,'survival':survival_mean_2,'age':mean_age_2,'metastasis':metas_mean_2,'size':mean_size_2})
		#print(survival_mean_1)
		print(ret_metadata)
		text_file_3 = open("polls/static/metadata.txt", "w")
		#text_file_3.write("<table><tr><th>Patient Cluster</th><th>Correlation with clinical Group</th><th>Mean Survival</th><th>mean age at diagnosis</th><th>mean tumor size</th></tr>")
		text_file_3.write("<table><tr><th>Patient Cluster</th><th>Gender=female</th><th>Relapse</th><th>Living</th></tr>")
		for elem in ret_metadata:				
			text_file_3.write("<tr>")
			text_file_3.write("<th>" + str(elem['gender']) + "</th>")
			text_file_3.write("<th>" + str(elem['relapse']) + "</th>")
			text_file_3.write("<th>" + str(elem['living']) + "</th>")
			text_file_3.write("</tr>")
		text_file_3.write("</table>")
		text_file_3.close()
	if(clinicalstr != "empty"):		
		ret_metadata = []
		nbr_patients = float(len(group1_ids) + len(group2_ids))
		survival_1 = []
		survival_2 = []
		survival_mean_1 = 0
		survival_mean_2 = 0
		ctr_surv_1 = 0.001
		sum_surv_1 = 0
		ctr_surv_2 = 0.001
		sum_surv_2 = 0	
		for key in patientData:
			if key in group1_ids:
				survival_1.append(patientData[key])
			elif key in group2_ids:
				survival_2.append(patientData[key])
		for elem in survival_1:
			sum_surv_1 = sum_surv_1 + float(elem)
			ctr_surv_1 = ctr_surv_1 + 1
		for elem in survival_2:
			sum_surv_2 = sum_surv_2 + float(elem)
			ctr_surv_2 = ctr_surv_2 + 1
		survival_mean_1 = sum_surv_1 / ctr_surv_1
		survival_mean_2 = sum_surv_2 / ctr_surv_2
		survival_perc_1 = {0:1}
		survival_perc_2 = {0:1}
		for i in range(1,10):	
			tmp1 = 1.0
			tmp2 = 1.0
			for k in survival_1:
				if(float(k) < float(i)):
					tmp1 = tmp1 - (1.0/len(survival_1))
			for k in survival_2:
				if(float(k) < float(i)):
					tmp2 = tmp2 - (1.0/len(survival_2))
			survival_perc_1.update({i:tmp1})
			survival_perc_2.update({i:tmp2})
		trace1 = go.Scatter(
		x=list(survival_perc_1.keys()),
		y=list(survival_perc_1.values()),
		mode='lines+markers',
		name="'Group 1'",
		hoverinfo='name',
		line=dict(
		shape='hv'))
		trace2 = go.Scatter(
		x=list(survival_perc_2.keys()),
		y=list(survival_perc_2.values()),
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
		#plot_div=plotly.offline.plot(fig, auto_open=False,include_plotlyjs = False, output_type='div')
		plot_div=plotly.offline.plot(fig, auto_open=False,output_type='div')
		with open("polls/static/output_plotly.html", "w") as text_file_2:
       			text_file_2.write(plot_div)
		#plotly.offline.plot(fig,auto_open=False, image_filename="tempimage.png", image='png')
	with open("polls/static/output_console.txt", "w") as text_file:
        	text_file.write("Done!")
	print(ret_metadata)
	return(script,div,plot_1,plot_div,ret_metadata)

@shared_task(name="script_output_task_6")
def script_output_task_6(T,row_colors1,col_colors1,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1_ret,jac_2_ret,save_data,username):
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
	plot.add_tools(HoverTool(tooltips=None),TapTool(),BoxSelectTool())
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
	source = ColumnDataSource(data=dict(x=x, y=y,Name=[node_labels[i] for i in genes.keys()]))
	labels2 = LabelSet(x='x', y='y', text='Name',text_font_size="8pt",x_offset=-20, source=source,background_fill_color='white',background_fill_alpha=0.0,level='glyph',render_mode='canvas')
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
	patientData = {}
	patientData_misc = {}
	group1_in_metadata = []
	group2_in_metadata = []
	metastasis = []
	no_metastasis = []
	group1_has_bm = []
	group2_has_bm = []
	bm = []
	no_bm = []
	group1_has_lm = []
	group2_has_lm = []
	lm = []
	no_lm = []
	group1_has_lymph = []
	group2_has_lymph = []
	lymph = []
	no_lymph = []
	group1_has_path_er = []
	group2_has_path_er = []
	path_er = []
	no_path_er = []
	group1_has_path_pr = []
	group2_has_path_pr = []
	path_pr = []
	no_path_pr = []
	group1_has_prognosis = []
	group2_has_prognosis = []
	prognosis_good = []
	prognosis_bad = []
	if(clinicalstr == "empty"):
		ret_metadata = []
		ret_metadata.append({'group':"Group 1",'lm':"NA",'bm':"NA",'lymph':"NA",'metastasis':"NA",'path_er':"NA",'path_pr':"NA",'good_prognosis':"NA"})
		ret_metadata.append({'group':"Group 2",'lm':"NA",'bm':"NA",'lymph':"NA",'metastasis':"NA",'path_er':"NA",'path_pr':"NA",'good_prognosis':"NA"})
		#ret_metadata.append({'group':"Group 1",'jac':jac_1_ret,'survival':"NA",'age':"NA",'size':"NA"})
		#ret_metadata.append({'group':"Group 2",'jac':jac_2_ret,'survival':"NA",'age':"NA",'size':"NA"})
		text_file_3 = open("polls/static/metadata.txt", "w")
		text_file_3.write("<table><tr><th>Patient Cluster</th><th>Lung metastasis</th><th>Breast metastasis</th><th>Lymph nodes</th><th>Metastasis</th><th>Path_ER</th><th>Path_PR</th><th>Good Prognosis</th></tr>")
		for elem in ret_metadata:
			text_file_3.write("<tr>")
			text_file_3.write("<th>" + str(elem['lm']) + "</th>")
			text_file_3.write("<th>" + str(elem['bm']) + "</th>")
			text_file_3.write("<th>" + str(elem['lymph']) + "</th>")
			text_file_3.write("<th>" + str(elem['metastasis']) + "</th>")
			text_file_3.write("<th>" + str(elem['path_er']) + "</th>")
			text_file_3.write("<th>" + str(elem['path_pr']) + "</th>")
			text_file_3.write("<th>" + str(elem['good_prognosis']) + "</th>")
			text_file_3.write("</tr>")
		text_file_3.write("</table>")
		text_file_3.close()	
		plot_div = ""
		with open("polls/static/output_plotly.html", "w") as text_file_2:
        		text_file_2.write("")
	
	else:
		clinicalLines = clinicalstr.split("\n")
		for line in clinicalLines:
			#bababa = line.replace('"','')
			bababa2 = line.split("\",\"")
			bababa2[0] = bababa2[0].replace('"','')
			if(len(bababa2) > 64):
				#print(bababa2[53])	
				#print(bababa2[63])
				#print(bababa2[72])
				if(bababa2[64] != "--" and bababa2[64] != "name:ch1" and bababa2[64] != "NA" and bababa2[64].replace('.','',1).isdigit()):
					patientData.update({bababa2[0]:bababa2[64]})
				#	print(patientData)
			#if(len(bababa2) > 71):
			if(len(bababa2) > 63):
				if(bababa2[63].isdigit()):
					print("patient id" + str(bababa2[0]))
					if(bababa2[0] in group1_ids):
						print("foo")
						group1_in_metadata.append(bababa2[0])
					elif(bababa2[0] in group2_ids):
						group2_in_metadata.append(bababa2[0])
					if(int(bababa2[63]) == 1):
						metastasis.append(bababa2[0])
					elif(int(bababa2[63]) == 0):
						no_metastasis.append(bababa2[0])
					tmp97 = [0.0,float(bababa2[63]),0.0]
					#patientData_misc.update({bababa2[0]:tmp97})
				if(bababa2[61].isdigit()):
					if(bababa2[0] in group1_ids):
						group1_has_lm.append(bababa2[0])
					elif(bababa2[0] in group2_ids):
						group1_has_lm.append(bababa2[0])
					if(int(bababa2[61]) == 1):
						lm.append(bababa2[0])
					elif(int(bababa2[61]) == 0):
						no_lm.append(bababa2[0])
				if(bababa2[58].isdigit()):
					if(bababa2[0] in group1_ids):
						group1_has_bm.append(bababa2[0])
					elif(bababa2[0] in group2_ids):
						group1_has_bm.append(bababa2[0])
					if(int(bababa2[58]) == 1):
						bm.append(bababa2[0])
					elif(int(bababa2[58]) == 0):
						no_bm.append(bababa2[0])
				if(len(bababa2[15]) > 1 and len(bababa2[15]) < 20):
					if(bababa2[15].split(": ")[1] == "N"):
						if(bababa2[0] in group1_ids):
							group1_has_path_er.append(bababa2[0])
						elif(bababa2[0] in group2_ids):
							group2_has_path_er.append(bababa2[0])
						no_path_er.append(bababa2[0])
					if(bababa2[15].split(": ")[1] == "P"):
						if(bababa2[0] in group1_ids):
							group1_has_path_er.append(bababa2[0])
						elif(bababa2[0] in group2_ids):
							group2_has_path_er.append(bababa2[0])
						path_er.append(bababa2[0])
				if(len(bababa2[16]) > 1 and len(bababa2[16]) < 20):
					if(bababa2[16].split(": ")[1] == "N"):
						if(bababa2[0] in group1_ids):
							group1_has_path_pr.append(bababa2[0])
						elif(bababa2[0] in group2_ids):
							group2_has_path_pr.append(bababa2[0])
						no_path_pr.append(bababa2[0])
					if(bababa2[16].split(": ")[1] == "P"):
						if(bababa2[0] in group1_ids):
							group1_has_path_pr.append(bababa2[0])
						elif(bababa2[0] in group2_ids):
							group2_has_path_pr.append(bababa2[0])
						path_pr.append(bababa2[0])
				if(len(bababa2) > 71):	
					if(bababa2[0] in group1_ids):
						group1_has_prognosis.append(bababa2[0])
					elif(bababa2[0] in group2_ids):
						group2_has_prognosis.append(bababa2[0])
					if(bababa2[71] == "Good Prognosis"):
						prognosis_good.append(bababa2[0])
					elif(bababa2[71] == "Bad Prognosis"):
						prognosis_bad.append(bababa2[0])
				if(len(bababa2) > 70):				
					if(bababa2[70].isdigit()):
						if(bababa2[0] in group1_ids):
							group1_has_lymph.append(bababa2[0])
						elif(bababa2[0] in group2_ids):
							group1_has_lymph.append(bababa2[0])
						if(int(bababa2[70]) == 1):
							lymph.append(bababa2[0])
						elif(int(bababa2[70]) == 0):
							no_lymph.append(bababa2[0])
				elif(bababa2[53].replace('.','',1).isdigit() and bababa2[63].isdigit()):
					print("patient id" + str(bababa2[0]))
					if(bababa2[0] in group1_ids):
						print("foo")
						group1_in_metadata.append(bababa2[0])
					elif(bababa2[0] in group2_ids):
						group2_in_metadata.append(bababa2[0])
					if(int(bababa2[63]) == 1):
						metastasis.append(bababa2[0])
					elif(int(bababa2[63]) == 0):
						no_metastasis.append(bababa2[0])
					#tmp97 = [float(bababa2[53]),float(bababa2[63]),0.0]
					#patientData_misc.update({bababa2[0]:tmp97})
				
		age1 = 0
		ct1 = 0.001
		age2 = 0
		ct2 = 0.001
		print(patientData_misc)
		#for key in patientData_misc:
		#	if key in group1_ids:
		#		age1 = age1 + patientData_misc[key][0]
		#		ct1 = ct1 + 1
		#		size1 = size1 + patientData_misc[key][2]
		#		ct3 = ct3 + 1
		#		ct5 = ct5 + 1
		#		metas1 = metas1 + patientData_misc[key][1]
		#	elif key in group2_ids:
		#		age2 = age2 + patientData_misc[key][0]
		#		ct2 = ct2 + 1
		#		size2 = size2 + patientData_misc[key][2]
		#		ct4 = ct4 + 1
		#		ct6 = ct6 + 1
		#		metas2 = metas2 + patientData_misc[key][1]
		ret_metadata = []
		nbr_patients = float(len(group1_ids) + len(group2_ids))
		survival_1 = []
		survival_2 = []
		survival_mean_1 = 0
		survival_mean_2 = 0
		ctr_surv_1 = 0.001
		sum_surv_1 = 0
		ctr_surv_2 = 0.001
		sum_surv_2 = 0
		for key in patientData:
			if key in group1_ids:
				survival_1.append(patientData[key])
			elif key in group2_ids:
				survival_2.append(patientData[key])
		for elem in survival_1:
			sum_surv_1 = sum_surv_1 + float(elem)
			ctr_surv_1 = ctr_surv_1 + 1
		for elem in survival_2:
			sum_surv_2 = sum_surv_2 + float(elem)
			ctr_surv_2 = ctr_surv_2 + 1
		survival_mean_1 = sum_surv_1 / ctr_surv_1
		survival_mean_2 = sum_surv_2 / ctr_surv_2
		#print(group1_in_metadata)
		#print(metastasis)
		if(float(len(group1_in_metadata) + len(group2_in_metadata)) / nbr_patients > 0.8):
			metas_jac = lib.jac(group1_in_metadata,metastasis)
			metas_jac_2 = lib.jac(group2_in_metadata,metastasis)
		else:
			metas_jac = "NA"
			metas_jac_2 = "NA"
		if(float(len(group1_has_lm) + len(group2_has_lm)) / nbr_patients > 0.8):
			lm_jac = lib.jac(group1_has_lm,lm)
			lm_jac_2 = lib.jac(group2_has_lm,lm)
		else:
			lm_jac = "NA"
			lm_jac_2 = "NA"
		if(float(len(group1_has_bm) + len(group2_has_bm)) / nbr_patients > 0.8):
			bm_jac = lib.jac(group1_has_bm,bm)
			bm_jac_2 = lib.jac(group2_has_bm,bm)
		else:
			bm_jac = "NA"
			bm_jac_2 = "NA"
		if(float(len(group1_has_lymph) + len(group2_has_lymph)) / nbr_patients > 0.8):
			lymph_jac = lib.jac(group1_has_lymph,lymph)
			lymph_jac_2 = lib.jac(group2_has_lymph,lymph)
		else:
			lymph_jac = "NA"
			lymph_jac_2 = "NA"
		if(float(len(group1_has_path_er) + len(group2_has_path_er)) / nbr_patients > 0.8):
			path_er_jac = lib.jac(group1_has_path_er,path_er)
			path_er_jac_2 = lib.jac(group2_has_path_er,path_er)
		else:
			path_er_jac = "NA"
			path_er_jac_2 = "NA"
		if(float(len(group1_has_path_er) + len(group2_has_path_er)) / nbr_patients > 0.8):
			path_pr_jac = lib.jac(group1_has_path_pr,path_pr)
			path_pr_jac_2 = lib.jac(group2_has_path_pr,path_pr)
		else:
			path_pr_jac = "NA"
			path_pr_jac_2 = "NA"
		if(float(len(group1_has_prognosis) + len(group2_has_prognosis)) / nbr_patients > 0.8):
			prognosis_jac = lib.jac(group1_has_prognosis,prognosis_good)
			prognosis_jac_2 = lib.jac(group2_has_prognosis,prognosis_good)
		else:
			prognosis_jac = "NA"
			prognosis_jac_2 = "NA"
		ret_metadata.append({'group':"Group 1",'lm':lm_jac,'bm':bm_jac,'lymph':lymph_jac,'metastasis':metas_jac,'path_er':path_er_jac,'path_pr':path_pr_jac,'good_prognosis':prognosis_jac})
		ret_metadata.append({'group':"Group 2",'lm':lm_jac_2,'bm':bm_jac_2,'lymph':lymph_jac_2,'metastasis':metas_jac_2,'path_er':path_er_jac_2,'path_pr':path_pr_jac_2,'good_prognosis':prognosis_jac_2})
		#ret_metadata.append({'group':"Group 1",'jac':jac_1_ret,'survival':survival_mean_1,'age':mean_age_1,'metastasis':metas_mean_1,'size':mean_size_1})
		#ret_metadata.append({'group':"Group 2",'jac':jac_2_ret,'survival':survival_mean_2,'age':mean_age_2,'metastasis':metas_mean_2,'size':mean_size_2})
		print(survival_mean_1)
		text_file_3 = open("polls/static/metadata.txt", "w")
		#text_file_3.write("<table><tr><th>Patient Cluster</th><th>Correlation with clinical Group</th><th>Mean Survival</th><th>mean age at diagnosis</th><th>mean tumor size</th></tr>")
		text_file_3.write("<table><tr><th>Patient Cluster</th><th>Lung metastasis</th><th>Breast metastasis</th><th>Lymph nodes</th><th>Metastasis</th><th>Path_ER</th><th>Path_PR</th><th>Good Prognosis</th></tr>")
		for elem in ret_metadata:
			text_file_3.write("<tr>")
			text_file_3.write("<th>" + str(elem['lm']) + "</th>")
			text_file_3.write("<th>" + str(elem['bm']) + "</th>")
			text_file_3.write("<th>" + str(elem['lymph']) + "</th>")
			text_file_3.write("<th>" + str(elem['metastasis']) + "</th>")
			text_file_3.write("<th>" + str(elem['path_er']) + "</th>")
			text_file_3.write("<th>" + str(elem['path_pr']) + "</th>")
			text_file_3.write("<th>" + str(elem['good_prognosis']) + "</th>")
			text_file_3.write("</tr>")
		text_file_3.write("</table>")
		text_file_3.close()
		survival_perc_1 = {0:1}
		survival_perc_2 = {0:1}
		for i in range(1,10):	
			tmp1 = 1.0
			tmp2 = 1.0
			for k in survival_1:
				if(float(k) < float(i)):
					tmp1 = tmp1 - (1.0/len(survival_1))
			for k in survival_2:
				if(float(k) < float(i)):
					tmp2 = tmp2 - (1.0/len(survival_2))
			survival_perc_1.update({i:tmp1})
			survival_perc_2.update({i:tmp2})
		#output_notebook()	
		trace1 = go.Scatter(
		x=list(survival_perc_1.keys()),
		y=list(survival_perc_1.values()),
		mode='lines+markers',
		name="'Group 1'",
		hoverinfo='name',
		line=dict(
		shape='hv'))
		trace2 = go.Scatter(
		x=list(survival_perc_2.keys()),
		y=list(survival_perc_2.values()),
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
		#plot_div=plotly.offline.plot(fig, auto_open=False,include_plotlyjs = False, output_type='div')
		plot_div=plotly.offline.plot(fig, auto_open=False,output_type='div')
		with open("polls/static/output_plotly.html", "w") as text_file_2:
        		text_file_2.write(plot_div)
		#plotly.offline.plot(fig,auto_open=False, image_filename="tempimage.png", image='png')
	if(save_data == "true"):
		foobar = "user_uploaded_files/" + username
		if not(os.path.isdir(foobar)):
			os.mkdir(foobar)
		filename_1 = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S_%f')[:-3]
		filename_2 = foobar + "/" + username + filename_1 + "_json.json"
		filename_3 = foobar + "/" + username + filename_1 + "_heatmap.png"
		copyfile("polls/static/test.png",filename_3)
		copyfile("polls/static/test15.json",filename_2)
	
	with open("polls/static/output_console.txt", "w") as text_file:
        	text_file.write("Done!")
	print(ret_metadata)
	return(script,div,plot_1,plot_div,ret_metadata)

@shared_task(name="script_output_task_5")
def script_output_task_5(T,row_colors1,col_colors1,G2,means,genes_all,adjlist,genes1,group1_ids,group2_ids,clinicalstr,jac_1_ret,jac_2_ret):
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
	plot.add_tools(HoverTool(tooltips=None),TapTool(),BoxSelectTool())
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
	source = ColumnDataSource(data=dict(x=x, y=y,Name=[node_labels[i] for i in genes.keys()]))
	labels2 = LabelSet(x='x', y='y', text='Name',text_font_size="8pt",x_offset=-20, source=source,background_fill_color='white',background_fill_alpha=0.0,level='glyph',render_mode='canvas')
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
	patientData = {}
	patientData_misc = {}
	group1_in_metadata = []
	group2_in_metadata = []
	metastasis = []
	no_metastasis = []
	if(clinicalstr == "empty"):
		ret_metadata = []
		ret_metadata.append({'group':"Group 1",'jac':jac_1_ret,'survival':"NA",'age':"NA",'size':"NA"})
		ret_metadata.append({'group':"Group 2",'jac':jac_2_ret,'survival':"NA",'age':"NA",'size':"NA"})
		text_file_3 = open("polls/static/metadata.txt", "w")
		text_file_3.write("<table><tr><th>Patient Cluster</th><th>Correlation with clinical Group</th><th>Mean Survival</th><th>mean age at diagnosis</th><th>mean tumor size</th></tr>")
		for elem in ret_metadata:
			text_file_3.write("<tr>")
			text_file_3.write("<th>" + str(elem['group']) + "</th>")
			text_file_3.write("<th>" + str(elem['jac']) + "</th>")
			text_file_3.write("<th>" + str(elem['survival']) + "</th>")
			text_file_3.write("<th>" + str(elem['age']) + "</th>")
			text_file_3.write("<th>" + str(elem['size']) + "</th>")
			text_file_3.write("</tr>")
		text_file_3.write("</table>")
		text_file_3.close()	
		plot_div = ""
		with open("polls/static/output_plotly.html", "w") as text_file_2:
        		text_file_2.write("")
	
	else:
		clinicalLines = clinicalstr.split("\n")
		for line in clinicalLines:
			#bababa = line.replace('"','')
			bababa2 = line.split("\",\"")
			bababa2[0] = bababa2[0].replace('"','')
			if(len(bababa2) > 64):
				print(bababa2[53])	
				print(bababa2[63])
				#print(bababa2[72])
				if(bababa2[64] != "--" and bababa2[64] != "name:ch1" and bababa2[64] != "NA" and bababa2[64].replace('.','',1).isdigit()):
					patientData.update({bababa2[0]:bababa2[64]})
					print(patientData)
			#if(len(bababa2) > 71):
			if(len(bababa2) > 63):
				#print(bababa2[72])
				#if(bababa2[57].replace('.','',1).isdigit() and bababa2[72].replace('.','',1).isdigit() and bababa2[63].isdigit()):
				#	tmp97 = [float(bababa2[57]),float(bababa2[63]),float(bababa2[72])]
				#	patientData_misc.update({bababa2[0]:tmp97})

				if(bababa2[63].isdigit()):
					print("patient id" + str(bababa2[0]))
					if(bababa2[0] in group1_ids):
						print("foo")
						group1_in_metadata.append(bababa2[0])
					elif(bababa2[0] in group2_ids):
						group2_in_metadata.append(bababa2[0])
					if(int(bababa2[63]) == 1):
						metastasis.append(bababa2[0])
					elif(int(bababa2[63]) == 0):
						no_metastasis.append(bababa2[0])
					tmp97 = [0.0,float(bababa2[63]),0.0]
					patientData_misc.update({bababa2[0]:tmp97})
				elif(bababa2[53].replace('.','',1).isdigit() and bababa2[63].isdigit()):
					print("patient id" + str(bababa2[0]))
					if(bababa2[0] in group1_ids):
						print("foo")
						group1_in_metadata.append(bababa2[0])
					elif(bababa2[0] in group2_ids):
						group2_in_metadata.append(bababa2[0])
					if(int(bababa2[63]) == 1):
						metastasis.append(bababa2[0])
					elif(int(bababa2[63]) == 0):
						no_metastasis.append(bababa2[0])
					tmp97 = [float(bababa2[53]),float(bababa2[63]),0.0]
					patientData_misc.update({bababa2[0]:tmp97})
				
		age1 = 0
		ct1 = 0.001
		age2 = 0
		ct2 = 0.001
		size1 = 0
		ct3 = 0.001
		size2 = 0
		ct4 = 0.001
		metas1 = 0
		ct5 = 0.001
		metas2 = 0
		ct6 = 0.001
		print(patientData_misc)
		for key in patientData_misc:
			if key in group1_ids:
				age1 = age1 + patientData_misc[key][0]
				ct1 = ct1 + 1
				size1 = size1 + patientData_misc[key][2]
				ct3 = ct3 + 1
				ct5 = ct5 + 1
				metas1 = metas1 + patientData_misc[key][1]
			elif key in group2_ids:
				age2 = age2 + patientData_misc[key][0]
				ct2 = ct2 + 1
				size2 = size2 + patientData_misc[key][2]
				ct4 = ct4 + 1
				ct6 = ct6 + 1
				metas2 = metas2 + patientData_misc[key][1]
		mean_age_1 = age1 / ct1
		mean_age_2 = age2 / ct2
		mean_size_1 = size1 / ct3
		mean_size_2 = size2 / ct4
		metas_mean_1 = metas1 / ct5
		metas_mean_2 = metas2 / ct6
		ret_metadata = []
		survival_1 = []
		survival_2 = []
		survival_mean_1 = 0
		survival_mean_2 = 0
		ctr_surv_1 = 0.001
		sum_surv_1 = 0
		ctr_surv_2 = 0.001
		sum_surv_2 = 0
		for key in patientData:
			if key in group1_ids:
				survival_1.append(patientData[key])
			elif key in group2_ids:
				survival_2.append(patientData[key])
		for elem in survival_1:
			sum_surv_1 = sum_surv_1 + float(elem)
			ctr_surv_1 = ctr_surv_1 + 1
		for elem in survival_2:
			sum_surv_2 = sum_surv_2 + float(elem)
			ctr_surv_2 = ctr_surv_2 + 1
		survival_mean_1 = sum_surv_1 / ctr_surv_1
		survival_mean_2 = sum_surv_2 / ctr_surv_2
		print(group1_in_metadata)
		print(metastasis)
		metas_jac = lib.jac(group1_in_metadata,metastasis)
		metas_jac_2 = lib.jac(group2_in_metadata,metastasis)
		ret_metadata.append({'group':"Group 1",'jac':jac_1_ret,'survival':survival_mean_1,'age':mean_age_1,'metastasis':metas_jac})
		ret_metadata.append({'group':"Group 2",'jac':jac_2_ret,'survival':survival_mean_2,'age':mean_age_2,'metastasis':metas_jac_2})
		#ret_metadata.append({'group':"Group 1",'jac':jac_1_ret,'survival':survival_mean_1,'age':mean_age_1,'metastasis':metas_mean_1,'size':mean_size_1})
		#ret_metadata.append({'group':"Group 2",'jac':jac_2_ret,'survival':survival_mean_2,'age':mean_age_2,'metastasis':metas_mean_2,'size':mean_size_2})
		print(survival_mean_1)
		text_file_3 = open("polls/static/metadata.txt", "w")
		#text_file_3.write("<table><tr><th>Patient Cluster</th><th>Correlation with clinical Group</th><th>Mean Survival</th><th>mean age at diagnosis</th><th>mean tumor size</th></tr>")
		text_file_3.write("<table><tr><th>Patient Cluster</th><th>Correlation with clinical Group</th><th>Mean Survival</th><th>percent of metastasis</th><th>mean age at diagnosis</th></tr>")
		for elem in ret_metadata:
			text_file_3.write("<tr>")
			text_file_3.write("<th>" + str(elem['group']) + "</th>")
			text_file_3.write("<th>" + str(elem['jac']) + "</th>")
			text_file_3.write("<th>" + str(elem['survival']) + "</th>")
			text_file_3.write("<th>" + str(elem['metastasis']) + "</th>")
			text_file_3.write("<th>" + str(elem['age']) + "</th>")
			text_file_3.write("</tr>")
		text_file_3.write("</table>")
		text_file_3.close()
		survival_perc_1 = {0:1}
		survival_perc_2 = {0:1}
		for i in range(1,10):	
			tmp1 = 1.0
			tmp2 = 1.0
			for k in survival_1:
				if(float(k) < float(i)):
					tmp1 = tmp1 - (1.0/len(survival_1))
			for k in survival_2:
				if(float(k) < float(i)):
					tmp2 = tmp2 - (1.0/len(survival_2))
			survival_perc_1.update({i:tmp1})
			survival_perc_2.update({i:tmp2})
		#output_notebook()	
		trace1 = go.Scatter(
		x=list(survival_perc_1.keys()),
		y=list(survival_perc_1.values()),
		mode='lines+markers',
		name="'Group 1'",
		hoverinfo='name',
		line=dict(
		shape='hv'))
		trace2 = go.Scatter(
		x=list(survival_perc_2.keys()),
		y=list(survival_perc_2.values()),
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
		#plot_div=plotly.offline.plot(fig, auto_open=False,include_plotlyjs = False, output_type='div')
		plot_div=plotly.offline.plot(fig, auto_open=False,output_type='div')
		with open("polls/static/output_plotly.html", "w") as text_file_2:
        		text_file_2.write(plot_div)
		#plotly.offline.plot(fig,auto_open=False, image_filename="tempimage.png", image='png')
	with open("polls/static/output_console.txt", "w") as text_file:
        	text_file.write("Done!")
	print(ret_metadata)
	return(script,div,plot_1,plot_div,ret_metadata)


@shared_task(name="run_enrichment")
def run_enrichment(path,pval_enr):
	fh1 = open(path)
	gene_list = []
	lines = fh1.readlines()
	#lines = lines[1:]
	for line in lines:
		line.replace("\\n","")
		gene_list.append(line)
	print("baabababa")
	print(pval_enr)
	enr = gp.enrichr(gene_list=gene_list,
                 description='test_name',
                 # gene_sets='KEGG_2016',
                 # or gene_sets='KEGG_2016,KEGG_2013',
                 gene_sets=['KEGG_2016','KEGG_2013'],
                 outdir="polls/data/test/enrichr_kegg",
                 cutoff=float(pval_enr) # test dataset, use lower value of range(0,1)
                )
	print(enr.results)
	return(enr.results)

@shared_task(name="run_enrichment_2")
def run_enrichment_2(path,pval_enr,out_dir):
	fh1 = open(path)
	gene_list = []
	lines = fh1.readlines()
	#lines = lines[1:]
	for line in lines:
		line.replace("\\n","")
		gene_list.append(line)
	print("baabababa")
	print(pval_enr)
	enr = gp.enrichr(gene_list=gene_list,
                 description='test_name',
                 # gene_sets='KEGG_2016',
                 # or gene_sets='KEGG_2016,KEGG_2013',
                 gene_sets=['KEGG_2016','KEGG_2013'],
                 outdir=out_dir,
                 cutoff=float(pval_enr) # test dataset, use lower value of range(0,1)
                )
	print(enr.results)
	return(enr.results)
@shared_task(name="run_go_enrichment")
def run_go_enrichment(path,pval_enr):
	fh1 = open(path)
	gene_list = []
	lines = fh1.readlines()
	#lines = lines[1:]
	for line in lines:
		line.replace("\\n","")
		gene_list.append(line)
	print("baabababa")
	print(pval_enr)
	enr = gp.enrichr(gene_list=gene_list,
                 description='test_name',
                 # gene_sets='KEGG_2016',
                 # or gene_sets='KEGG_2016,KEGG_2013',
                 gene_sets=['GO_Biological_Process_2018','GO_Cellular_Component_2018','GO_Molecular_Function_2018'],
                 outdir="polls/data/test/enrichr_go",
                 cutoff=float(pval_enr) # test dataset, use lower value of range(0,1)
                )
	print(type(enr.results))
	return(enr.results)

@shared_task(name="run_go_enrichment_2")
def run_go_enrichment_2(path,pval_enr,out_dir):
	fh1 = open(path)
	gene_list = []
	lines = fh1.readlines()
	#lines = lines[1:]
	for line in lines:
		line.replace("\\n","")
		gene_list.append(line)
	print("baabababa")
	print(pval_enr)
	enr = gp.enrichr(gene_list=gene_list,
                 description='test_name',
                 # gene_sets='KEGG_2016',
                 # or gene_sets='KEGG_2016,KEGG_2013',
                 gene_sets=['GO_Biological_Process_2018','GO_Cellular_Component_2018','GO_Molecular_Function_2018'],
                 outdir=out_dir,
                 cutoff=float(pval_enr) # test dataset, use lower value of range(0,1)
                )
	print(type(enr.results))
	return(enr.results)


@shared_task(name="run_reac_enrichment")
def run_reac_enrichment(path,pval_enr,out_dir):
	fh1 = open(path)
	gene_list = []
	lines = fh1.readlines()
	#lines = lines[1:]
	for line in lines:
		line.replace("\\n","")
		gene_list.append(line)
	enr = gp.enrichr(gene_list=gene_list,
                 description='test_name',
                 # gene_sets='KEGG_2016',
                 # or gene_sets='KEGG_2016,KEGG_2013',
                 gene_sets=['Reactome_2013','Reactome_2016'],
                 outdir=out_dir,
                 cutoff=float(pval_enr) # test dataset, use lower value of range(0,1)
                )
	print(type(enr.results))
	return(enr.results)



@shared_task(name="read_kegg_enrichment")
def read_kegg_enrichment(path,pval_enr):
	result_file = open(path)
	ret_dict = []
	ctr = 0
	for line in result_file:
		tmp = {}
		lineSplit = line.split("\t")
		if(ctr > 0):	
			if(float(lineSplit[3]) < float(pval_enr)):
				#print(lineSplit[3])
				#for i in range(0,len(lineSplit)):
				#	tmp[i] = lineSplit[i]
				#ret_dict.append(tmp)
				for i in range(0,5):
					tmp[i] = lineSplit[i]
				tmp[5] = lineSplit[9]
				ret_dict.append(tmp)
		ctr = ctr + 1
	return(ret_dict)


@shared_task(name="read_kegg_enrichment_2")
def read_kegg_enrichment_2(path1,path2,pval_enr):
	result_file_1 = open(path1)
	result_file_2 = open(path2)
	temp_dict = {}
	temp_dict_2 = {}
	ret_dict = []
	ret_dict_2 = []
	ctr = 0
	for line in result_file_1:
		tmp = {}
		lineSplit = line.split("\t")
		if(ctr > 0):	
			if(float(lineSplit[3]) < float(pval_enr)):
				#print(lineSplit[3])
				#for i in range(0,len(lineSplit)):
				#	tmp[i] = lineSplit[i]
				#ret_dict.append(tmp)
				#temp_dict.update({lineSplit[0]:tmp})
				for i in range(0,5):
					tmp[i] = lineSplit[i]
				tmp[5] = lineSplit[9]
				#ret_dict.append(tmp)
				temp_dict.update({lineSplit[0]:tmp})
		ctr = ctr + 1

	ctr2 = 0
	for line in result_file_2:
		tmp = {}
		lineSplit = line.split("\t")
		if(ctr2 > 0):	
			if(float(lineSplit[3]) < float(pval_enr)):
				#print(lineSplit[3])
				#for i in range(0,len(lineSplit)):
				#	tmp[i] = lineSplit[i]
				#ret_dict.append(tmp)
				for i in range(0,5):
					tmp[i] = lineSplit[i]
				tmp[5] = lineSplit[9]
				temp_dict_2.update({lineSplit[0]:tmp})
		ctr = ctr + 1
	for key in temp_dict:
		if(key not in temp_dict_2):
			ret_dict.append(temp_dict[key])
	for key in temp_dict_2:
		if(key not in temp_dict):
			ret_dict_2.append(temp_dict[key])

	return(ret_dict,ret_dict_2)
		


@shared_task(name="read_ndex_file")
def read_ndex_file(fn):
	#fh1 = open(fn,encoding ='utf-8')
	#lines = fh1.readlines()
	lines = fn.split("\n")
	lines[3].replace("@","")
	lines[3].replace("uniprot:","uniprot")
	lines[3].replace("signor:","signor")
	lines[3].replace(" ","")
	lines[3].replace("\\n","")
	#lines[3].replace("
	node_line = lines[3][11:-3]
	nodelinesplit = node_line.split(", ")
	#print(nodelinesplit)
	dictlist = []
	node_dict = {}
	print(nodelinesplit[len(nodelinesplit)-1])
	for item in nodelinesplit:
		print(item)
		tmp = json.loads(item)
		dictlist.append(tmp)
		if not(any(c.islower() for c in tmp['n'])):
			node_dict[tmp['@id']] = tmp['n']
			#print(tmp['@id'])
		#print(dictlist)
	lines[4].replace("@","")
	lines[4].replace("uniprot:","uniprot")
	lines[4].replace("signor:","signor")
	lines[4].replace(" ","")
	lines[4].replace("\\n","")
	edge_line = lines[4][11:-4].rstrip()
	edgelinesplit = edge_line.split(", ")
	#print(edgelinesplit)
	edgelist = []
	dataset = Dataset(name='hsapiens_gene_ensembl',
		                  host='http://www.ensembl.org')
	conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name','entrezgene'])
	ret = []
	for item in edgelinesplit:
		print(item)
		tmp = json.loads(item)
		#dictlist.append(tmp)
		if(tmp['s'] in node_dict and tmp['t'] in node_dict):
			source = node_dict[tmp['s']]
			target = node_dict[tmp['t']]
			#print(source + "\t" + target)
			prot1_nbr = conv.index[conv['Gene name'] == source]
			prot2_nbr = conv.index[conv['Gene name'] == target]
			if(source in conv['Gene name'].tolist() and target in conv['Gene name'].tolist()):
				prot1 = conv.loc[prot1_nbr,'NCBI gene ID'].values[0]
				prot2 = conv.loc[prot2_nbr,'NCBI gene ID'].values[0]
				if(prot1 != prot2):
					baz = [str(int(prot1)), str(int(prot2))]
					ret.append("\t".join(baz))
	return("\n".join(ret))



@shared_task(name="read_ndex_file_2")
def read_ndex_file_2(fn):
	#fh1 = open(fn,encoding ='utf-8')
	#lines = fh1.readlines()
	lines99 = fn.split(", {\"edges\":[")
	lines98 = lines99[1].split(", {\"nodes\":[")
	line4=lines98[0]
	lines97 = lines98[1].split(", {\"status\"")
	line3=lines97[0]
	#print(line3)
	#print(line4)
	#return ""
	line3.replace("@","")
	line3.replace("uniprot:","uniprot")
	line3.replace("signor:","signor")
	line3.replace(" ","")
	line3.replace("\\n","")
	#lines[3].replace("
	node_line = line3[:-2]
	nodelinesplit = node_line.split(", ")
	#print(nodelinesplit)
	dictlist = []
	node_dict = {}
	#print(nodelinesplit[len(nodelinesplit)-1])
	for item in nodelinesplit:
		#print(item)
		tmp = json.loads(item)
		dictlist.append(tmp)
		if not(any(c.islower() for c in tmp['n'])):
			node_dict[tmp['@id']] = tmp['n']
			#print(tmp['@id'])
		#print(dictlist)
	line4.replace("@","")
	line4.replace("uniprot:","uniprot")
	line4.replace("signor:","signor")
	line4.replace(" ","")
	line4.replace("\\n","")
	edge_line = line4[:-2].rstrip()
	edgelinesplit = edge_line.split(", ")
	#print(edgelinesplit)
	edgelist = []
	dataset = Dataset(name='hsapiens_gene_ensembl',
		                  host='http://www.ensembl.org')
	conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name','entrezgene'])
	ret = []
	ctr = 0
	conv_genelist = conv['Gene name'].tolist()
	for item in edgelinesplit:
		ctr = ctr + 1
		if(ctr > 23000):
			print(item)
		tmp = json.loads(item)
		#dictlist.append(tmp)
		if(tmp['s'] in node_dict and tmp['t'] in node_dict):
			source = node_dict[tmp['s']]
			target = node_dict[tmp['t']]
			#print(source + "\t" + target)
			prot1_nbr = conv.index[conv['Gene name'] == source]
			prot2_nbr = conv.index[conv['Gene name'] == target]
			if(source in conv_genelist and target in conv_genelist):
				prot1 = conv.loc[prot1_nbr,'NCBI gene ID'].values[0]
				prot2 = conv.loc[prot2_nbr,'NCBI gene ID'].values[0]
				if(prot1 != prot2 and str(prot1).replace('.','',1).isdigit() and str(prot2).replace('.','',1).isdigit()):
					baz = [str(int(prot1)), str(int(prot2))]
					ret.append("\t".join(baz))
					#print("\t".join(baz))
	
	with open("polls/static/ppi_list.txt", "w") as text_file_97:
		text_file_97.write("\n".join(ret))		
	return("\n".join(ret))

@shared_task(name="read_ndex_file_3")
def read_ndex_file_3(fn):
	#fh1 = open(fn,encoding ='utf-8')
	#lines = fh1.readlines()
	lines99 = fn.split(" {\"nodes\":[")
	lines98 = lines99[1].split("{\"edges\":[")
	line3=lines98[0]
	#lines97 = lines98[1].split(" {\"networkAttributes\"")
	line4=lines98[1]
	#print(line3)
	#print(line4)
	#return ""
	line3.replace("@","")
	line3.replace("uniprot:","uniprot")
	line3.replace("signor:","signor")
	line3.replace(" ","")
	line3.replace("\\n","")
	#lines[3].replace("
	node_line = line3[:-5]
	nodelinesplit = node_line.split(", ")
	#print(nodelinesplit)
	dictlist = []
	node_dict = {}
	#print(nodelinesplit[len(nodelinesplit)-1])
	for item in nodelinesplit:
		print(item)
		tmp = json.loads(item)
		dictlist.append(tmp)
		if not(any(c.islower() for c in tmp['n'])):
			node_dict[tmp['@id']] = tmp['n']
			print(tmp['@id'])
		#print(dictlist)
	line4.replace("@","")
	line4.replace("uniprot:","uniprot")
	line4.replace("signor:","signor")
	line4.replace(" ","")
	line4.replace("\\n","")
	edge_line = line4[:-4].rstrip()
	edgelinesplit = edge_line.split(", ")
	#print(edgelinesplit)
	edgelist = []
	dataset = Dataset(name='hsapiens_gene_ensembl',
		                  host='http://www.ensembl.org')
	#print(dataset
	conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name','entrezgene'])
	#Id = Conversion('dummyemail@dummybunny.info')

	print(conv)
	ret = []
	ctr = 0
	ctr2 = 0
	conv_genelist = conv['Gene name'].tolist()
	#print(edgelinesplit)
	with open("polls/static/ppi_list_2.txt", "w") as text_file_97:	
		for i in range(1,100000):
			item = edgelinesplit[i]
			#print(edgelinesplit[i])
			ctr = ctr + 1
			if(ctr == 100000):
				print(item)
				print(ctr2)
				ctr = 0
				ctr2 = ctr2 + 100000
				#if(ctr2 == 320000):
			tmp = json.loads(item)
			#dictlist.append(tmp)
			if(tmp['s'] in node_dict and tmp['t'] in node_dict):
				source = node_dict[tmp['s']]
				target = node_dict[tmp['t']]
				#id1 = Id.convert_uniprot_to_entrez(source)
				#id2 = Id.convert_uniprot_to_entrez(target)
				#text_file_97.write("\t".join([id1,id2]))
				#print("\t".join([id1,id2]))
				#print(source + "\t" + target)
				prot1_nbr = conv.index[conv['Gene name'] == source]
				prot2_nbr = conv.index[conv['Gene name'] == target]
				if(source in conv_genelist and target in conv_genelist):
					prot1 = conv.loc[prot1_nbr,'NCBI gene ID'].values[0]
					prot2 = conv.loc[prot2_nbr,'NCBI gene ID'].values[0]
					if(prot1 != prot2 and str(prot1).replace('.','',1).isdigit() and str(prot2).replace('.','',1).isdigit()):
						baz = [str(int(prot1)), str(int(prot2))]
						ret.append("\t".join(baz))
						#print("\t".join(baz))
						text_file_97.write("\t".join(baz))	
						text_file_97.write("\n")	
		
	#for item in edgelinesplit:
	#	ctr = ctr + 1
	#	if(ctr == 10000):
	#		print(item)
	#		print(ctr2)
	#		ctr = 0
	#		ctr2 = ctr2 + 10000
	#		if(ctr2 == 320000):
	#			with open("polls/static/ppi_list_2.txt", "w") as text_file_97:
	#				text_file_97.write("\n".join(ret))		
	#	tmp = json.loads(item)
	#	#dictlist.append(tmp)
	#	if(tmp['s'] in node_dict and tmp['t'] in node_dict):
	#		source = node_dict[tmp['s']]
	#		target = node_dict[tmp['t']]
	#		#print(source + "\t" + target)
	#		prot1_nbr = conv.index[conv['Gene name'] == source]
	#		prot2_nbr = conv.index[conv['Gene name'] == target]
	#		if(source in conv_genelist and target in conv_genelist):
	#			prot1 = conv.loc[prot1_nbr,'NCBI gene ID'].values[0]
	#			prot2 = conv.loc[prot2_nbr,'NCBI gene ID'].values[0]
	#			if(prot1 != prot2 and str(prot1).replace('.','',1).isdigit() and str(prot2).replace('.','',1).isdigit()):
	#				baz = [str(int(prot1)), str(int(prot2))]
	#				ret.append("\t".join(baz))
	#				#print("\t".join(baz))
	
	#with open("polls/static/ppi_list_2.txt", "w") as text_file_97:
	#	text_file_97.write("\n".join(ret))		
	return("\n".join(ret))



@shared_task(name="read_ndex_file_4")
def read_ndex_file_4(fn):
	lines6 = ""
	if("edges" in fn.split("nodes")[1]):
		lines5 = fn.split("{\"nodes\":[")
		lines3 = lines5[1].split("{\"edges\":[")[0]
		if("cyTableColumn" in lines5[1]):
			lines4 = lines5[1].split("{\"edges\":[")[1].split("{\"cyTableColumn\":[")[0]
			lines4 = lines4[:-4]
		elif("networkAttributes" in lines5[1]):
			lines4 = lines5[1].split("{\"edges\":[")[1].split("{\"networkAttributes\":[")[0]
			lines4 = lines4[:-4]
			if("nodeAttributes" in lines5[1].split("{\"edges\":[")[1].split("{\"networkAttributes\":[")[1] and "UniprotName" in lines5[1].split("{\"edges\":[")[1].split("{\"networkAttributes\":[")[1]):
				lines6_temp = lines5[1].split("{\"edges\":[")[1].split("{\"networkAttributes\":[")[1].split("{\"nodeAttributes\":[")[1]
				lines6 = lines6_temp.split("{\"edgeAttributes\":[")[0]
		else:
			lines4 = lines5[1].split("{\"edges\":[")[1]
	elif("edges" in fn.split("nodes")[0]):
		print("badbaisrgoasrogasog")
		lines5 = fn.split("{\"nodes\":[")
		lines3 = lines5[1].split("]},")[0] + "]]]"
		lines4 = lines5[0].split("{\"edges\":[")[1][:-4]
	#lines4 = lines4 + "}"
	#lines4 = lines4[:-4]
	#lines = fn.split("\n")
	#lines[3].replace("@","")
	#lines[3].replace("uniprot:","uniprot")
	#lines[3].replace("signor:","signor")
	#lines[3].replace(" ","")
	#lines[3].replace("ncbigene:","")
	#lines[3].replace("\\n","")
	lines3.replace("@","")
	lines3.replace("uniprot:","uniprot")
	lines3.replace("signor:","signor")
	lines3.replace(" ","")
	lines3.replace("ncbigene:","")
	lines3.replace("\\n","")
	lines33 = lines3[:-3].replace("}]","")
	#print(lines33)
	#lines[3].replace("
	#node_line = lines[3][11:-3].replace("ncbigene:","")
	node_line = lines33.replace("ncbigene:","")
	#print(node_line)
	nodelinesplit = node_line.split(", ")
	#print(nodelinesplit)
	dictlist = []
	node_dict = {}
	if not(node_line.endswith("}")):
		node_line = node_line + "}"
	node_line_2 = "[" + node_line + "]"
	#print(nodelinesplit[len(nodelinesplit)-1])
	tmp2 = json.loads(node_line_2)
	#print(tmp2)
	node_dict_2 = {}
	if not(lines6 == ""):
		lines6 = "[" + lines6
		#print(lines6[:-4])
		tmp4 = json.loads(lines6[:-4])
		for item in tmp4:
			if(item['n'] == "GeneName_A"):
				node_dict_2[item['po']] = item['v']
				print(str(item['po']) + " " + str(item['v']))
	print(node_dict_2)
	dataset = Dataset(name='hsapiens_gene_ensembl',
		                  host='http://www.ensembl.org')
	conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name','entrezgene'])
	conv_genelist = conv['Gene name'].tolist()
	for item in tmp2:
		#print(item)
		#tmp = json.loads(item)
		dictlist.append(item)
		if('r' in item):
			if(any(c.islower() for c in item['r'])):
				gene_name = item['n']
				if(gene_name in conv_genelist):
					gene_nbr = conv.index[conv['Gene name'] == gene_name]
					gene_nbr1 = conv.loc[gene_nbr,'NCBI gene ID'].values[0]
					node_dict[item['@id']] = gene_nbr1
					print(item)
			else:
				node_dict[item['@id']] = item['r']
				print(item)
				#print(tmp['@id'])
		else:
			if(item['n'].isdigit()):
				print(item)
				node_dict[item['@id']] = item['n']
			elif (item['n'] in node_dict_2):
				gene_name = node_dict_2[item['n']]
				print(gene_name)
				if(gene_name in conv_genelist):
					gene_nbr = conv.index[conv['Gene name'] == gene_name]
					gene_nbr1 = conv.loc[gene_nbr,'NCBI gene ID'].values[0]
					node_dict[item['@id']] = gene_nbr1
					print(gene_nbr1)
	print(node_dict)
	#print(dictlist)
	#lines[4].replace("@","")
	#lines[4].replace("uniprot:","uniprot")
	#lines[4].replace("signor:","signor")
	#lines[4].replace(" ","")
	#lines[4].replace("\\n","")
	lines4.replace("@","")
	lines4.replace("uniprot:","uniprot")
	lines4.replace("signor:","signor")
	lines4.replace(" ","")
	#lines4 = lines4.replace("\\n","")
	lines4 = lines4.replace("]","")
	edge_line = lines4.rstrip()
	edge_line_2 = "[" + edge_line + "]"
	#print(edge_line_2)
	edgelinesplit = edge_line.split(", ")	
	#print(edgelinesplit)
	edgelist = []
	#tmp4 = json.loads(edge_line_2.replace('\r\n', '\\r\\n').rstrip())
	#tmp4 = json.loads(edge_line_2.replace('\r\n', '').rstrip())
	tmp4 = json.loads(edge_line_2)
	dataset = Dataset(name='hsapiens_gene_ensembl',
		                  host='http://www.ensembl.org')
	conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name','entrezgene'])
	ret = []
	for item in tmp4:
		print(item)
		#dictlist.append(tmp)
		if(item['s'] in node_dict and item['t'] in node_dict):
			source = node_dict[item['s']]
			target = node_dict[item['t']]
			print(source)
			print(target)
			if(source != target and not(math.isnan(float(source))) and not(math.isnan(float(target)))):
				baz = [str(int(source)),str(int(target))]
				ret.append("\t".join(baz))
	print("\n".join(ret))
	return("\n".join(ret))
	

@shared_task(name="check_input_files")
def check_input_files(ppistr,exprstr):
	errstr = ""
	ppi_stringio = StringIO(ppistr)
	ppidf = pd.read_csv(ppi_stringio,sep='\t')
	print(ppidf)
	print(ppidf.iloc[0,:])
	if(len(ppidf.columns) < 2):
		errstr = errstr + "Input file must contain two columns with interaction partners.\n"
		#errstr = errstr + "\n \n To avoid this error, go to <a href=\"infopage.html\">the infopage</a> and make sure that your input data has the specified format."
		return(errstr)
	contains_numbers = "false"
	#for i in range(len(ppidf.index)):
	for i in range(len(ppidf.index)):
		if(contains_numbers == "false"):
			print(i)
			curr_elem = str(ppidf.iloc[[i], 0].values[0])
			#curr_elem = ppidf.iloc[[i], 0].split("\t")
			#curr_elem_1 = curr_elem[0]
			#curr_elem_2 = curr_elem[len(curr_elem)-1]
			curr_elem_2 = str(ppidf.iloc[[i], 1].values[0])
			print(curr_elem)
			print(curr_elem_2)
			if(curr_elem.isdigit() and curr_elem_2.isdigit()):
				contains_numbers = "true"
	if(contains_numbers == "false"):
		errstr = errstr + "\n" + "Input file must contain lines with Entrez IDs of interaction partners.\n"
		#errstr = errstr + "\n \n "
	return(errstr)
		#for j in range(len(df.columns)):
				
@shared_task(name="convert_gene_list")
def convert_gene_list(adjlist,filename):
	dataset = Dataset(name='hsapiens_gene_ensembl',
	                  host='http://www.ensembl.org')
	conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name','entrezgene'])
	conv_genelist = conv['Gene name'].tolist()
	retstr = ""
	for elem in adjlist:
		prot_1 = elem[0]
		prot_2 = elem[1]
		gene_nbr_1 = conv.index[conv['Gene name'] == prot_1]
		gene_nbr_1_2 = conv.loc[gene_nbr_1,'NCBI gene ID'].values[0]
		gene_nbr_2 = conv.index[conv['Gene name'] == prot_2]
		gene_nbr_2_2 = conv.loc[gene_nbr_2,'NCBI gene ID'].values[0]
		retstr = retstr + str(gene_nbr_1_2).split(".")[0] + "\t" + str(gene_nbr_2_2).split(".")[0] + "\n"
	#return(retstr)
	with open(filename, "w") as text_file:
		text_file.write(retstr)

@shared_task(name="import_ndex")
def import_ndex(name):
	nice_cx_network = ndex2.create_nice_cx_from_server(server='public.ndexbio.org', uuid=name)
	tmp4 = []
	node_dict = {}
	dataset = Dataset(name='hsapiens_gene_ensembl',
		                  host='http://www.ensembl.org')
	conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name','entrezgene'])
	conv_genelist = conv['Gene name'].tolist()
	for node_id, node in nice_cx_network.get_nodes():
		current_node = {}
		current_node['id'] = node_id
		current_node['n'] = node.get('n')
		#current_node['r'] = node.get('r')
		if(name=="9c38ce6e-c564-11e8-aaa6-0ac135e8bacf"):
			curr_gene_name = nice_cx_network.get_node_attribute_value(node_id,'GeneName_A')
			if(curr_gene_name in conv_genelist):
				gene_nbr = conv.index[conv['Gene name'] == curr_gene_name]
				gene_nbr1 = conv.loc[gene_nbr,'NCBI gene ID'].values[0]
				if not(math.isnan(float(gene_nbr1))):
					node_dict[node_id] = str(int(gene_nbr1))	
				#print(gene_nbr1)
		else:
			curr_gene_name = current_node['n']				
			if(curr_gene_name in conv_genelist):
				gene_nbr = conv.index[conv['Gene name'] == curr_gene_name]
				gene_nbr1 = conv.loc[gene_nbr,'NCBI gene ID'].values[0]
				if not(math.isnan(float(gene_nbr1))):
					node_dict[node_id] = str(int(gene_nbr1))
				#print(gene_nbr1)
		tmp4.append(current_node)
	edgelist = []
	ret = ""
	for edge_id, edge in nice_cx_network.get_edges():
		source = edge.get('s')
		target = edge.get('t')
		if(source in node_dict and target in node_dict and source != target):
			curr_edge_str = str(node_dict[source]) + "\t" + str(node_dict[target]) + "\n"
			edgelist.append([node_dict[source],node_dict[target]])
			ret = ret + curr_edge_str
	print(node_dict)
	return(ret)

@shared_task(name="import_ndex_2")
def import_ndex_2(name):
	for item in tmp2:
		#print(item)
		#tmp = json.loads(item)
		dictlist.append(item)
		if('r' in item):
			if(any(c.islower() for c in item['r'])):
				gene_name = item['n']
				if(gene_name in conv_genelist):
					gene_nbr = conv.index[conv['Gene name'] == gene_name]
					gene_nbr1 = conv.loc[gene_nbr,'NCBI gene ID'].values[0]
					node_dict[item['@id']] = gene_nbr1
					print(item)
			else:
				node_dict[item['@id']] = item['r']
				print(item)
				#print(tmp['@id'])
		else:
			if(item['n'].isdigit()):
				print(item)
				node_dict[item['@id']] = item['n']
			elif (item['n'] in node_dict_2):
				gene_name = node_dict_2[item['n']]
				print(gene_name)
				if(gene_name in conv_genelist):
					gene_nbr = conv.index[conv['Gene name'] == gene_name]
					gene_nbr1 = conv.loc[gene_nbr,'NCBI gene ID'].values[0]
					node_dict[item['@id']] = gene_nbr1
					print(gene_nbr1)
	print(node_dict)
	#print(dictlist)
	#lines[4].replace("@","")
	#lines[4].replace("uniprot:","uniprot")
	#lines[4].replace("signor:","signor")
	#lines[4].replace(" ","")
	#lines[4].replace("\\n","")
	lines4.replace("@","")
	lines4.replace("uniprot:","uniprot")
	lines4.replace("signor:","signor")
	lines4.replace(" ","")
	#lines4 = lines4.replace("\\n","")
	lines4 = lines4.replace("]","")
	edge_line = lines4.rstrip()
	edge_line_2 = "[" + edge_line + "]"
	#print(edge_line_2)
	edgelinesplit = edge_line.split(", ")	
	#print(edgelinesplit)
	edgelist = []
	#tmp4 = json.loads(edge_line_2.replace('\r\n', '\\r\\n').rstrip())
	#tmp4 = json.loads(edge_line_2.replace('\r\n', '').rstrip())
	tmp4 = json.loads(edge_line_2)
	dataset = Dataset(name='hsapiens_gene_ensembl',
		                  host='http://www.ensembl.org')
	conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name','entrezgene'])
	ret = []
	for item in tmp4:
		print(item)
		#dictlist.append(tmp)
		if(item['s'] in node_dict and item['t'] in node_dict):
			source = node_dict[item['s']]
			target = node_dict[item['t']]
			print(source)
			print(target)
			if(source != target and not(math.isnan(float(source))) and not(math.isnan(float(target)))):
				baz = [str(int(source)),str(int(target))]
				ret.append("\t".join(baz))
	print("\n".join(ret))
	return("\n".join(ret))
	


@shared_task(name="show_old_data")
def show_old_data(path_json,path_heatmap):
	#G = nx.Graph()
	#G_2 = nx.Graph()
	copyfile(path_heatmap,"polls/static/test.png")
	copyfile(path_json,"polls/static/test15.json")
	


