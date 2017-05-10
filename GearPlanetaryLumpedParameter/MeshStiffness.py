#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May  5 22:18:55 2017

@author: wehrle
"""
from __future__ import absolute_import, division, print_function
from numpy import cos, linspace, pi
import matplotlib.pyplot as plt
kt = 10
dpw = 10
alpha = linspace(0, pi/3, 100)
kg = 4*kt/((dpw*cos(alpha))**2)
plt.clf()
plt.plot(alpha, kg)
plt.show()

KA = 