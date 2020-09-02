#!/usr/bin/env python3
import os
import json
import numpy as np
import volestipy


current_directory = os.getcwd()
input_file_json = current_directory +  '/bigg_files/e_coli_core.json'
input_file_mat = current_directory +  '/bigg_files/e_coli_core.mat'


met_net = volestipy.read_json_file(input_file_json)

print("A")
print(met_net[0])
print("------------------------")
print("b")
print(met_net[1])
print("-------------------------")
print("Aeq")
print(met_net[2])
print("-------------------------")
print("beq")
print(met_net[3])
print("-------------------------")
print("Metabolites are the following " + str(len(met_net[4])) + ":")
print(met_net[4])
print("-------------------------")
print("Reactions are the following " + str(len(met_net[5])) + ":")
print(met_net[5])

# -----------------------------------------------------------------

met_net = volestipy.read_mat_file(input_file_mat)

print("A")
print(met_net[0])
print("------------------------")
print("b")
print(met_net[1])
print("-------------------------")
print("Aeq")
print(met_net[2])
print("-------------------------")
print("beq")
print(met_net[3])
print("-------------------------")
print("Metabolites are the following " + str(len(met_net[4])) + ":")
print(met_net[4])
print("-------------------------")
print("Reactions are the following " + str(len(met_net[5])) + ":")
print(met_net[5])

