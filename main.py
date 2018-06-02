#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 24 20:51:30 2018

@author: lunar
"""

import os
import json
import pickle
os.chdir('/home/lunar/Desktop/IgoR_crawling/CODE_Project/')
from load_data import *
from utils import *

import project_maintainment as main_pj
import project_results as res_pj

### primers 
Primers = json.load(open('/home/lunar/Desktop/IgoR_crawling/projects/Primers.json', 'r'))
Core = '/home/lunar/Desktop/IgoR_crawling/projects/'
project_names = ['project_mouse_2', 'project_mouse_3', 'Immunized', 'PD1_BALBC', 'project_mouse_memoryT']
Multi_1 = main_pj.Multi(core=Core, primers=Primers, projects=project_names)
Multi_1.load_frames()
Multi_1.load_results()