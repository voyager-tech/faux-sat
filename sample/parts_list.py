# -*- coding: utf-8 -*-
"""
Created on Thu Jul  5 18:40:09 2018

@author: Darren
"""

# This file shows classifies all the data for the components in Fox Sat.
# Still has a lot of variables that need to be added!

# Utilized Scripts
import sample.parts

# %% Parts List

# Part Directory (list)
Parts_Dir = []
# Parts List (Dictionary)
Parts_List = {}

for variable in dir(parts):
    if '__' in variable:
        pass
    elif variable == 'u' or variable == 'pint':
        pass
    else:
        Parts_Dir.append(variable)

for part in Parts_Dir:
    Parts_List[part] = eval('parts.' + str(part))














