#For Linux
from itertools import *
import pandas as pd
# import tkinter as tk
# from tkinter import *
# from tkinter import ttk
# from tkinter import scrolledtext
import subprocess as sub
import os
from argparse import ArgumentParser
import argparse

#Find Cell types  python Document
#find the most likely cell type base on the literature mining and "CELL MARKER"
#Usage: 
#    print_cool_to_screen(marker_file,pub_db_file,cluster_marker_file,top_n=40,speci_thred=2)
#Arguments:
##top n: default=20. The top genes of the conserve markers in each cluster.
##speci_thred: default=2. The threthold of the specificity, how many cell type the marker map, the bigger threshold means the lower specificity
##demo thredshold: the threshold should be 2 3 4 5,this parameter will takes alots of minutes,Because of the random selection process.
##############################################setting####################################################    
pj = lambda *paths: os.path.abspath(os.path.join(*paths))
cur_dir=os.getcwd()
need_file_dir=pj(cur_dir,"need_file")
marker_file=pj(need_file_dir,"marker_sheet.xlsx")
cluster_marker_file=pj(need_file_dir,"all_conserved_markers.xlsx")
pub_db_file=pj(need_file_dir+"/cellmarekerDB_new.xlsx")
########################################################################################################
parser = ArgumentParser(description='This Script is to find the optimal cell marker for single cell analysing')
parser.add_argument('-s', action='store', dest='strict_degree', default=2,type=int,help='the strict degree of gene filter, e.g gene A appears in two types of cell, if strict is 2,this gene is filtered,if 3 retain it otherwise')
parser.add_argument('-t', action='store', dest='top_n', default=20, type=int, help='top n gene of the concserved markers, default=20')
#######################################################################################################
def load_marker(marker_file):
    ##对marker数据进行数据预处理,将相同类型的细胞聚集在一起，包括其细胞亚类，并且做了每一个元素去除空格的操作
    ####导入所有cluster的细胞markers
    all_marker_df=pd.read_excel(marker_file)
    return all_marker_df
def load_cluster_marker(cluster_marker_file):
    all_cluster_marker_df=pd.read_excel(cluster_marker_file)
    return all_cluster_marker_df
def remove_marker_from_dict(marker_dict,marker_list):
    ##find the keys
    new_dict={}
    for key,value in marker_dict.items():
        new_value=[i for i in value if i not in marker_list]
        new_dict[key]=new_value
    marker_dict.update(new_dict)
    return marker_dict
    
            

def filter_marker(marker_dict,speci_thred):
    ##if the marker was not specificity,these marker should be remove.
    cell_types=list(marker_dict.keys())
    all_combination=permutations(cell_types, speci_thred)
#     print(all_combination)
    marker_set_list=[]
    if speci_thred == 4:
        for combination in all_combination:
            for element in combination:
                this_marker_set=set(marker_dict[element])
                marker_set_list.append(this_marker_set)
            if set.intersection(marker_set_list[0],marker_set_list[1],marker_set_list[2]):
                select_marker=list(set.intersection(marker_set_list[0],marker_set_list[1],marker_set_list[2],marker_set_list[3]))
                marker_dict=remove_marker_from_dict(marker_dict,select_marker)
        return marker_dict
    if speci_thred == 3:
        for combination in all_combination:
            for element in combination:
                this_marker_set=set(marker_dict[element])
                marker_set_list.append(this_marker_set)
            if set.intersection(marker_set_list[0],marker_set_list[1],marker_set_list[2]):
                select_marker=list(set.intersection(marker_set_list[0],marker_set_list[1],marker_set_list[2]))
                marker_dict=remove_marker_from_dict(marker_dict,select_marker)
        return marker_dict
    if speci_thred == 2:
        for combination in all_combination:
#             print(combination)
            marker_set_list=[]
            for element in combination:
                this_marker_set=set(marker_dict[element])
                marker_set_list.append(this_marker_set)
            if set.intersection(marker_set_list[0],marker_set_list[1]):
#                 print(combination)
#                 print(marker_set_list[0],marker_set_list[1])
                select_marker=list(set.intersection(marker_set_list[0],marker_set_list[1]))
                marker_dict=remove_marker_from_dict(marker_dict,select_marker)
#         print(select_marker)
        return marker_dict

            
            
def get_cluster_dict(all_cluster_marker_df):
    max_cluster=all_cluster_marker_df["cluster"].max()
    assert max_cluster==17 , "the clusters was not equal to 17, please check your cluster numbers as the Seurat results"
    gene_dict={}
    for cluster in range(max_cluster+1):
        gene_list=all_cluster_marker_df[all_cluster_marker_df["cluster"]==cluster]["gene"]
        gene_dict[cluster]=gene_list
    return gene_dict
def get_marker_dict(all_marker_df):
    marker_dict={}
    for cell_type in all_marker_df.columns:
        marker_dict[cell_type]=list(all_marker_df[cell_type])
    return marker_dict

#-----------------------------------------------owned marker repository------------------------------------------------ 
#-----------------------------------------------------//\\-------------------------------------------------------------
#----------------------------------------------------///\\\------------------------------------------------------------
#------------------------------------------------------||--------------------------------------------------------------
#------------------------------------------------------||--------------------------------------------------------------
#----------------------------------------------------\\\///------------------------------------------------------------
#-----------------------------------------------------\\//-------------------------------------------------------------
#------------------------------------------------used public dataset---------------------------------------------------
def load_pub_db(pub_db_file):
    pub_db_df=pd.read_excel(pub_db_file)
    return pub_db_df 
def get_db_dict(pub_db_df):
    pub_marker_dict={}
    need_col=["Cell Type","Gene"]
    for line in range(pub_db_df.shape[0]):
        cell_type=pub_db_df["Cell Type"][line]
        gene=pub_db_df["Gene"][line]
        pub_marker_dict.setdefault(cell_type,[]).append(gene)
    return pub_marker_dict
    


def find_cell_type(marker_file,pub_db_file,cluster_marker_file,top_n,speci_thred):
    all_marker_df=load_marker(marker_file)
    all_cluster_marker_df=load_cluster_marker(cluster_marker_file)
    gene_dict=get_cluster_dict(all_cluster_marker_df)
    marker_dict=get_marker_dict(all_marker_df)
    marker_dict=filter_marker(marker_dict,speci_thred)
    pub_db_df=load_pub_db(pub_db_file)
    pub_dict=get_db_dict(pub_db_df)
    pub_dict=filter_marker(pub_dict,speci_thred)
    result_dict={}
    for cluster,gene_list in gene_dict.items():
        if len(gene_list)<top_n:
            top_gene=set(gene_list)
        else:
            top_gene=set(gene_list[0:top_n])
        for cell_type,marker_list in marker_dict.items():
            cell_type=cell_type.strip()
            if top_gene.intersection(marker_list):
                
                match_gene=top_gene.intersection(marker_list)
                result_dict.setdefault(cluster,{}).setdefault(cell_type,match_gene)
    ##if the cluster not match let the remain cluster to match the public database
    pub_result_dict={}
    not_match=[]
    for cluster,gene_list in gene_dict.items():
        if len(gene_list)<top_n:
            top_gene=set(gene_list)
        else:
            top_gene=set(gene_list[0:top_n])
        match_num=0
        if cluster not in result_dict.keys():
            
            not_match.append(cluster)
            for cell_type,marker_list in pub_dict.items():
                cell_type=cell_type.strip()
                if top_gene.intersection(marker_list):
                    match_gene=top_gene.intersection(marker_list)
                    pub_result_dict.setdefault(cluster,{}).setdefault(cell_type,match_gene)
                    result_dict.update(pub_result_dict)
    
    return result_dict,not_match
def print_cool_to_screen(marker_file,pub_db_file,cluster_marker_file,top_n,speci_thred):
    result_dict,not_match=find_cell_type(marker_file,pub_db_file,cluster_marker_file,top_n,speci_thred)
    str_result=""
    for cluster,marker_dict in sorted(result_dict.items()):
        list_marker_sorted=[]
        list_cell_type_sorted=[]
        most_likely_cell_type=sorted(marker_dict, key=lambda key: len(marker_dict[key]),reverse=True)[0]       
        str_result+="The cluster of %s match the %d possible cell type:%s, the most likely is the %s \n" % (cluster,len(marker_dict),marker_dict,most_likely_cell_type)
    str1="==========================================WARNNING======================================="
    str2="The cluster of %s were not matched in our own database; Their result was derived from the database from the website of \"CELL MARKER\",\nThese cell types should be think over. If they were really odded,which means you just found the new cell types,you should define them by yourself,otherwise!!!" % str(not_match)
    return str_result+"\n"+str1+"\n"+str2
        ##marker 文件##待会把他换成conserve markers
arg = parser.parse_args()
top_n = arg.top_n 
strict_degree=arg.strict_degree
if __name__ == "__main__":
    result=print_cool_to_screen(marker_file,pub_db_file,cluster_marker_file,top_n,strict_degree)
    print(result)
# win=tk.Tk()
# win.title("Cell marker")
# label1=ttk.Label(win,text="Enter top n:")
# label1.grid(column=0,row=0)
# label2=ttk.Label(win,text="Strict degree")
# label2.grid(column=1,row=0)
# def ClickMe():
#     button.configure(text="BingGo"+top_n.get()+strict_degree.get())
# def insert_text():
# #     text = Text(win,width=400,height=40)
#     scr=scrolledtext.ScrolledText(win,width=400,height=40,wrap=tk.WORD)
#     scr.grid(column=0,row=3,columnspan=3)
#     result=print_cool_to_screen(marker_file,pub_db_file,cluster_marker_file,int(top_n.get()),int(strict_degree.get()))
#     scr.insert(INSERT,result)
# button=ttk.Button(win,text="Run",command=insert_text)
# button.grid(column=2,row=1)
# ##add textbox
# top_n=tk.StringVar()
# textbox=ttk.Entry(win,width=12,textvariable=top_n)
# textbox.grid(column=0,row=1)
# textbox.focus()
# #add choosed box
# strict_degree=tk.StringVar()
# numberChoosen=ttk.Combobox(win,width=12,textvariable=strict_degree,state="readonly")
# numberChoosen["values"]=(2,3,4)
# numberChoosen.grid(column=1,row=1)
# numberChoosen.current(0)
# win.mainloop()