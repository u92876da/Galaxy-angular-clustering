#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 13:13:15 2021

@author: u92876da
"""

import config as cf
from field_geometry import field_geometry

class integral_constraint:
    
    def __init__(self, fieldname):
        
        self.fieldname = fieldname
        self.field_geometry = field_geometry(fieldname)