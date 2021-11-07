#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 13:05:23 2020

@author: ofir
"""
import matplotlib._color_data as mcd
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import numpy as np
import os, sys
import matplotlib.patches as mpatches
import random

#For reproduceability
seed = 123
random.seed(seed)
np.random.seed(seed)

"""
Run as:
python3 Network_Figures.py <This is base folder> <drop fragments with size smaller then> \
                          <integer of hubness to filter> <prefix> <All_compounds_B_input txt file> \
                          <Seeds_B_input txt file> <ECs_All_input txt file> \
                          <patches_Compounds> <patches_Enzymes> <Pathways_Compounds>
                          

"""

BaseFolder = "C:/Users/removal_network/order/"
drop_fragmen_with_size = 1
filter_hubness = 100
outputprefix = "BjSa_Order_"


All_compounds_B_input = "allCompounds_BjSa.txt"
f = open(BaseFolder+All_compounds_B_input, "r")
All_compounds_B_input = f.readline()
All_compounds_B_input = All_compounds_B_input.split(";")
f.close()

Seeds_B_input = "Seeds_BjSa.txt"
f = open(BaseFolder+Seeds_B_input, "r")
Seeds_B_input = f.readline()
Seeds_B_input = Seeds_B_input.split(";")
f.close()

ECs_All_input = "EC_ALL.txt"
f = open(BaseFolder+ECs_All_input, "r")
ECs_All_input = f.readline()
ECs_All_input = ECs_All_input.split(";")
f.close()


patches_Compounds_ = "Compounds_BjSa_Order.txt"
f = open(BaseFolder+patches_Compounds_, "r")
patches_Compounds_ = f.readlines()
patches_Compounds_ = [i.strip("\n") for i in patches_Compounds_]
patches_Compounds_ = [i.split(" ") for i in patches_Compounds_]
f.close()


patches_Enzymes_ = "Enzymes_BjSa_Order_1.txt"
f = open(BaseFolder+patches_Enzymes_, "r")
patches_Enzymes_ = f.readlines()
patches_Enzymes_ = [i.strip("\n") for i in patches_Enzymes_]
patches_Enzymes_ = [i.split(" ") for i in patches_Enzymes_]
f.close()

Pathway_Compounds_ = "Pathways_BjSa_order.txt"
f = open(BaseFolder+Pathway_Compounds_, "r")
Pathway_Compounds_ = f.readlines()
Pathway_Compounds_ = [i.strip("\n").split(";") for i in Pathway_Compounds_]
for i, vali in enumerate(Pathway_Compounds_):
    for j, valj in enumerate(vali):
        Pathway_Compounds_[i][j] = valj.strip()
f.close()
print(Pathway_Compounds_)


High_alpha_pathways = ['Glucosinolate biosynthesis', 'Linoleic acid metabolism', 'Tyrosine metabolism', 'Geraniol degradation', 'Porphyrin and chlorophyll metabolism']

#Load database files
#########################
df_el_ = pd.read_csv(BaseFolder+"full_enzymes_labels_jun.txt", sep="|")
df_ecMapping_ = pd.read_csv(BaseFolder+"ec_reac_mapping_jun.txt", sep="|")
df_reactions_ = pd.read_csv(BaseFolder+"reactions_3_balanced.txt", sep="|")
df_ec_to_compoundIndex_ = pd.read_csv(BaseFolder+"compounds_lables_jun_1.txt", sep="|")



def CreateCompoundsNetwork(All_compounds_B=All_compounds_B_input,
                           Seeds_B=Seeds_B_input,
                           ECs_All=ECs_All_input,
                           df_el = df_el_,
                           df_ecMapping = df_ecMapping_,
                           df_reactions = df_reactions_,
                           df_ec_to_compoundIndex = df_ec_to_compoundIndex_,
                           drop_fragmen_with_size = 1,#size of fragments to drop
                           filter_hubness = filter_hubness,#filter by max number of node connections
                           FinalFigureName = outputprefix,
                           patches_Compounds = patches_Compounds_,
                           patches_Enzymes = patches_Enzymes_,
                           Pathway_Compounds = Pathway_Compounds_):


    df_el = df_el_
    df_ecMapping = df_ecMapping_
    df_reactions = df_reactions_
    df_ec_to_compoundIndex = df_ec_to_compoundIndex_
    
    
    All_compounds_B = All_compounds_B+Seeds_B    
    
    IndexEnzymeList = df_el[["Index", "Enzyme_Code"]].values.tolist()
    DictEL = {}
    
    for i in IndexEnzymeList:
        DictEL[i[1]] = i[0]
        DictEL[i[0]] = i[1]
        
    #EClist = [DictEL[i] for i in EClist]
    EClist_tmp = []
    Errors = []
    
    for i in ECs_All:
        try:
            EClist_tmp.append(DictEL[i])
        except:
            Errors.append(i)
            print("Cant find EC "+i)
    #Write ECs couldnt be found
    with open(BaseFolder+'unfound_ECs.txt', 'w') as filehandle:
        filehandle.writelines("%s\n" % x for x in Errors)
    
    df_ecMapping = df_ecMapping[df_ecMapping["Enzyme_index"].isin(EClist_tmp)]
    ListOfRelevantReactions = df_ecMapping["Reactions"].values.tolist()
    
    df_ecMapping["Reactions_list"] = [ i.split(":")[1].strip().split(" ") for i in ListOfRelevantReactions]
    
    #This to be used as Reaction to Enzyme_index 
    df_ecMapping_expolded = df_ecMapping[["Enzyme_index", "Reactions_list"]].explode("Reactions_list")
    
    flat_list = []
            
    for i in ListOfRelevantReactions:
        for j in i.split(":")[1].strip().split(" "):
            try:
                flat_list.append(int(j.strip()))
            except:
                print("invalid for int")
    
    #Save enzyme information of each reaction.
    
    df_reactions = df_reactions[df_reactions["Reactions"].isin(flat_list)]
    
    #Fix reaction directions - flip =
    
    tempList = []
    
    for i in df_reactions.values.tolist():
        if i[1] == "=":
            x = i
            x[3], x[2] = x[2], x[3]
            tempList.append(x)
    
    df_reactions = pd.concat([df_reactions, pd.DataFrame(tempList, columns=list(df_reactions.columns))])
    
    #Save the Enzymes linked to reaction. 
    #reactions = df_reactions["Reactions"].values.tolist()
    
    
    l = [i.split(":")[1].lstrip().rstrip().split(",") for i in df_reactions["Left"].values.tolist()]
    r = [i.split(":")[1].lstrip().rstrip().split(",") for i in df_reactions["Right"].values.tolist()]
    
      
    DictComp = {}
    
    for i in df_ec_to_compoundIndex[["Index", "Compound_Code"]].values.tolist():
        DictComp[i[1]] = i[0]
        DictComp[i[0]] = i[1]
        
      
    df = pd.DataFrame([l, r, df_reactions["Reactions"].values.tolist()]).T           
    df = df.explode(0)         
    df = df.explode(1)
    
    #map reaction to enzyme (col 4)
    ##Might be a problem here, when multiple enzymes are related to a single reaction
    
    ReactionToEnzymeDict = { i[1]:i[0] for i in df_ecMapping_expolded.values.tolist() }
    cols = ["Left", "Right", "Reaction"]
    df.columns = cols
    
    def map_enz(x):
        return ReactionToEnzymeDict[str(x)]
    
    df["Enzyme"] = df["Reaction"].apply(map_enz)
    
    df_grouped = df.groupby(["Left","Right"])["Enzyme"].apply(list)
    df_grouped = df_grouped.to_frame().reset_index()
    
    def enzInd_to_enzEC(l):
        return [DictEL[i] for i in l] 
    
    df_grouped["Enzyme"] = df_grouped["Enzyme"].apply(enzInd_to_enzEC)
    
    
    #Fix Compound IDs
    
    def CompInd_to_CompC(x):
        try:
            return DictComp[int(x)]
        except:
            return np.nan
    
    df_grouped["Left"] = df_grouped["Left"].apply(CompInd_to_CompC)
    df_grouped["Right"] = df_grouped["Right"].apply(CompInd_to_CompC)
    
    #Save only rows contain compounds from the all_compounds list
    #df_grouped = df_grouped[(df_grouped["Left"].isin(All_compounds_B+All_compounds_A)) &
    #                        (df_grouped["Right"].isin(All_compounds_B+All_compounds_A))]    
    
    df_grouped = df_grouped[(df_grouped["Left"].isin(All_compounds_B)) &
                            (df_grouped["Right"].isin(All_compounds_B))]    
    
  
  
    
    all_patch_colors = [
        'orange', 'green', 'blue' , 'red', 'aqua', 'magenta', 'lightskyblue', 'teal', 'darkviolet', 'limegreen', 'peachpuff', 'peachpuff', 'peachpuff',
        'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff',
        'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff',
        'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff',
        'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff',
        'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff',
        'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff', 'peachpuff','peachpuff'
    ]


 
    #Set node color
    node_colors_Dict = {}
    
    all_nodes = list(set(df_grouped["Left"].values.tolist()+df_grouped["Right"].values.tolist()))

    for i in all_nodes:
        node_colors_Dict[i] = "gray"

    for i, vali in enumerate(patches_Compounds):
        for j in vali:
            if j in all_nodes:
                node_colors_Dict[j] = all_patch_colors[i]
    
   
    

    
    def ColorEnzyme(l):
        #print(l)
        c = "gray"
        l = [str(i) for i in l]
        for i, vali in enumerate(patches_Enzymes):
            if [x for x in l if x in vali]:
                c = all_patch_colors[i]
                
        return c
        
    df_grouped["Color"] = df_grouped["Enzyme"].apply(ColorEnzyme)
    
    
    def setEdgeWidthByColor(x):
        if x != "gray":
            return 1.5
        else:
            return 0.75
        
    df_grouped["Edge_width"] = df_grouped["Color"].apply(setEdgeWidthByColor)
    
    df_grouped.to_csv(BaseFolder+"_"+outputprefix+"compounds_network_notFiltered_order.csv")
    
    #df_grouped["count"] = df.groupby(["Left","Right"]).count().reset_index()["Reaction"]
    
    
    #Filter df by nodes appearance
    temp_df = df_grouped["Left"].value_counts().reset_index()
    keeplist = temp_df[temp_df["Left"] <= filter_hubness]["index"].values.tolist()
    temp_df = df_grouped["Right"].value_counts().reset_index()
    keeplist += temp_df[temp_df["Right"] <= filter_hubness]["index"].values.tolist()
    keeplist = list(set(keeplist))
    
    df_grouped = df_grouped[(df_grouped["Left"].isin(keeplist)) & (df_grouped["Right"].isin(keeplist))]
    
    
    
    
    
    #Network
    G = nx.Graph()
    
    for i in df_grouped[["Left", "Right", "Color","Edge_width"]].values.tolist():
        G.add_edge(i[0], i[1], color=i[2], width=i[3])
        
    #Filter out isolates
    
    #Filter out fragments
    for component in list(nx.connected_components(G)):
        if len(component)<drop_fragmen_with_size:
            for node in component:
                G.remove_node(node)
    
    
    
    edges,colors = zip(*nx.get_edge_attributes(G,'color').items())
    edges,widths = zip(*nx.get_edge_attributes(G,'width').items())

    node_color = []
    
    for node in list(G.nodes()):
        node_color.append(node_colors_Dict[node])
          
    nodes_of_largest_component  = max(nx.connected_components(G), key = len)
    largest_component = G.subgraph(nodes_of_largest_component)
    
    pos1 = nx.spring_layout(G,k=0.075,iterations=200, seed=seed)
    
    pos2 = nx.spring_layout(G, pos=pos1,fixed=nodes_of_largest_component,
                           k=0.001,iterations=0,seed=seed)
    #Label colored nodes
    
    flatten_compounds = [item for sublist in patches_Compounds for item in sublist]
    
    nodes_labeldict = {}
    for i in G.nodes():
        if i in flatten_compounds:
            nodes_labeldict[i] = i
        else: nodes_labeldict[i] = ""
  

    fig = plt.figure(figsize=(50,50))
    nx.draw_networkx_edges(G, pos2, alpha=0.25, width=widths, edge_color=colors, arrows=False)
    nx.draw_networkx_nodes(G, pos2, node_color=node_color, node_size=30,
                                   alpha=0.5)
    nx.draw_networkx_labels(G,pos2,nodes_labeldict,font_size=4,font_color='black')
    
# =============================================================================
    PathwayDict = {}
    for i, vali in enumerate(Pathway_Compounds):
        for j in vali[1:]:
            if j in PathwayDict.keys():
                PathwayDict[j].append(vali[0])
            else:        
                PathwayDict[j] = [vali[0]]
    
    
    ax=plt.gca()

    #Create node-position and colors
    patchNodePositions = []
    PatchNodeNames = []
    handles = []
    handles2 = []

    #Legend 1
    for i, vali in enumerate(Pathway_Compounds):
        #node to position
        for j in vali:
            if j in list(nodes_of_largest_component):
                patchNodePositions.append(pos2[j])
                PatchNodeNames.append(j)
             
        if len(PatchNodeNames) > 0:
            handles.append(mpatches.Patch(color=all_patch_colors[i+len(patches_Compounds)], label=vali[0]))
        
        for cl in PatchNodeNames:
            # plot circles using the RGBA colors
            if len(set(PathwayDict[cl]).intersection(set(High_alpha_pathways))) > 0:
                circle = plt.Circle(pos2[cl], 0.01, color=all_patch_colors[i+len(patches_Compounds)], fill=True, alpha=0.8)
                ax.add_artist(circle)
                print(PathwayDict[cl])

            else:
                circle = plt.Circle(pos2[cl], 0.01, color=all_patch_colors[i+len(patches_Compounds)], fill=True, alpha=0.2)
                ax.add_artist(circle)
        
        PatchNodeNames = []



    #legend 2
    for i, vali in enumerate(patches_Compounds):  
        handles2.append(mpatches.Patch(color=all_patch_colors[i], label=vali[0]))
    

    # plot the legend


    leg1 = plt.legend(handles=handles, bbox_to_anchor=(1.05, 1), loc=('lower left'), fontsize=45, title="Pathways")
    

    
    leg1.get_title().set_fontsize('45') #legend 'Title' fontsize
    ax.add_artist(leg1)    
    leg2 = plt.legend(handles=handles2, bbox_to_anchor=(1.05, 0.15), loc=('upper left'), fontsize=45, title="Order")
    leg2.get_title().set_fontsize('45') #legend 'Title' fontsize
    ax.add_artist(leg2)    


    plt.axis('equal')    
    fig.savefig(BaseFolder+FinalFigureName+'Removal_Network_order_2.24.png', bbox_inches='tight')
    fig.savefig(BaseFolder+FinalFigureName+'Removal_Network_order_2.24.pdf', bbox_inches='tight')
    plt.close()
  
CreateCompoundsNetwork(All_compounds_B = All_compounds_B_input,
                           Seeds_B = Seeds_B_input,
                           ECs_All = ECs_All_input,
                           df_el = df_el_,
                           df_ecMapping = df_ecMapping_,
                           df_reactions = df_reactions_,
                           df_ec_to_compoundIndex = df_ec_to_compoundIndex_,
                           drop_fragmen_with_size = 1,#size of fragments to drop
                           filter_hubness = filter_hubness,#filter by max number of node connections                          
                           FinalFigureName = outputprefix,
                           patches_Compounds = patches_Compounds_,
                           patches_Enzymes = patches_Enzymes_,
                           Pathway_Compounds = Pathway_Compounds_)



