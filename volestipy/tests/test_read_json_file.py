#!/usr/bin/env python3
import os
import json
import numpy as np
import volestipy


current_directory = os.getcwd()
input_file_json = current_directory +  '/e_coli_core.json'



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

