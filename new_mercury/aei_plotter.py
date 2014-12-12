"""
Plots things

(1) Reads .aei files
(2) Plots orbital elements
(3) Makes movies of orbital elements?
(4) Plot planet(s) + stars on x-y grid (print inclination with each frame :-P)
(5) Plot planet(s) + stars on x-y/y-z grid?

"""

import numpy as np
from matplotlib import pyplot as pyplot
from useful import *

import pickle
import os
import sys


class Plotter():

	def __init__():
		pass

	def store_planets(planets = []):
		"""
		For a planet that will be frequently re-used, store the data to avoid the problem of re-reading it over and over
		"""

	def save_stored_planets(planets = []):
		"""
		Writes stored planets to pickle files
		"""


	def plot_sma(planets = [], separate_plots = False):
		""" plot sma of planets over time """
		pass

	def plot_ecc(planets = []):
		pass

	def plot_inc(planets = []):
		pass

	def plot_xy(planets = []):
		pass