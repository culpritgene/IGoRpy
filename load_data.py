#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 24 20:52:27 2018

@author: lunar
"""



import numpy as np
import pandas as pd
from scipy import stats
from skimage import measure

import matplotlib.pyplot as plt
import seaborn as sns

import glob
import subprocess
import warnings
import json
import pickle

import re
import datetime
import itertools
from collections import OrderedDict, Mapping
from functools import reduce

## Deprecated!
from sklearn.utils import shuffle
from scipy import stats
from copy import deepcopy
####

import importlib
import importlib.util

### Import Pygor module
path_to_pygor_m = '/home/lunar/IGoR/igor_from_source/IGoR/pygor/models/'
spec = importlib.util.spec_from_file_location("genmodel", 
                                              path_to_pygor_m+'genmodel.py')
genmod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(genmod)

# make standart subprocess call function because we will need it later on many occasions
sub_run = lambda x: subprocess.run('{}'.format(x), stdout=subprocess.PIPE, stderr=subprocess.STDOUT)








