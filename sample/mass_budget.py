# -*- coding: utf-8 -*-
"""
Created on Thu Jul  5 19:21:50 2018

@author: Darren
"""

# Utilized Scripts
import parts_list

# Import dictionary of all SC parts & attatched data
Parts_List = parts_list.Parts_List
mass = 0
Parts_Dir = parts_list.Parts_Dir

for part in Parts_Dir:
    mass = mass + (Parts_List[part]['mass'] * Parts_List[part]['quantity'])

print("Mass is: " + str(mass) + ' g')
