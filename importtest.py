# -*- coding: utf-8 -*-
"""
Created on Sun Mar 27 17:09:15 2016

@author: Billy
"""

import csv 

isotope_array=[]
counter = -1
prevz = 0
preva = 0

with open ('mathematicaexport1.csv','rb') as csvfile:
    reader = csv.DictReader(csvfile,fieldnames=['z','a','e_level','e_gamma','t_half','level_width','brj','br0','j0','jr','sigma_int'])
    for row in reader:
        if row['z'] != prevz or row['a'] != preva:
            prevz = row['z']
            preva = row['a']
            isotope_array.append([])
            counter += 1
            isotope_array[counter].append(row)
        else:
            isotope_array[counter].append(row)
            
isotope_array[0][0][None][0]

#for index1, isotope in enumerate(isotope_array): #convert everything to floats
#    for index2, state in enumerate(isotope):
#        for index3, key in enumerate(state):
#            print key
#            for index4, thing in enumerate(state[None]):
#                isotope_array[index1][index2][None][index4] = float(thing)
#            isotope_array[index1][index2][key] = float(state[key])
