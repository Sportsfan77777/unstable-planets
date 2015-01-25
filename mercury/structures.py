"""
Vector Data Structures
"""

import numpy as np

class Vector():

    def __init__(self, *args):
        self._vector = np.array(args, dtype = np.ndarray)
        self._size = len(args)
        self._dict = {}

    def pretty_print():
    	print self._vector
    	if len(self._dict):
    		print ,
    	print

    def get(self, name):
        return self._dict[name] # maybe handle exception also?

    def name_arg(self, name, arg_num):
        if 0 <= arg_num < self._size:
            self._dict[name] = self._vector[arg_num]

    def add_entry(self, entry, arg_num = -1):
    	raise Exception("Not implemented yet")

    def length(self, override = False):
        if _size == 3 or override:
            return np.sqrt(self.norm(override))
        else:
            return -1

    def norm(self, override = False):
        if _size == 3 or override:
            return np.linalg.norm(self._vector)
        else:
            return -1


class Position(Vector):

    def __init__(self, x, y, z):
        Vector.__init__(self, x, y, z)
        self.name_arg('x', 0)
        self.name_arg('y', 1)
        self.name_arg('z', 2)
        pass

class Velocity(Vector):

    def __init__(self, vx, vy, vz):
        Vector.__init__(self, vx, vy, vz)
        self.name_arg('vx', 0)
        self.name_arg('vy', 1)
        self.name_arg('vz', 2)
        pass


