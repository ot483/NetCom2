#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys, os



BaseFolder = "C:/Users/.../"
ResultsFolder = BaseFolder
FULL_PATH = os.path.abspath(BaseFolder)+"/"


sys.path.append(FULL_PATH)
from netcom import simulation

input_1 = BaseFolder+"for_simulation_conv.txt"
f = open(input_1, "r")
input_1_list = f.readlines()
input_1_list = [i.strip("\n").strip().split(" ") for i in input_1_list]
f.close()


input_2 = BaseFolder+"env.txt"
f = open(input_2, "r")
input_2_list = f.readlines()
input_2_list = [i.strip("\n").strip().split(" ") for i in input_2_list]
f.close()

input_3 = BaseFolder+"for_simulation_org.txt"
f = open(input_3, "r")
input_3_list = f.readlines()
input_3_list = [i.strip("\n").strip().split(" ") for i in input_3_list]
f.close()




Results = []
Results2 = []
#Simulation conventional
for i in input_1_list:
    sim_df_T1, steps_df_T1 = simulation(input1 = i.copy(), input2=[input_2_list[0][1::]], resfolder=ResultsFolder, prefix="test_", basefolder=BaseFolder)
    All_compounds_A_sim = sim_df_T1.values.tolist()[0][0]
    Results.append(i[0]+" "+" ".join(All_compounds_A_sim))

#save list of strings as text file
resFile = ResultsFolder+"simmulation_compounds_conventional.txt"
with open(resFile, 'w') as f:
    for item in Results:
        f.write("%s\n" % item)


#Simulation organic
for j in input_3_list:
    sim_df_T2, steps_df_T2 = simulation(input1 = j.copy(), input2=[input_2_list[1][1::]], resfolder=ResultsFolder, prefix="test_", basefolder=BaseFolder)
    All_compounds_B_sim = sim_df_T2.values.tolist()[0][0]
    Results2.append(j[0]+" "+" ".join(All_compounds_B_sim))

#save list of strings as text file
resFile = ResultsFolder+"simmulation_compounds_organic.txt"
with open(resFile, 'w') as f:
    for item in Results2:
        f.write("%s\n" % item)











