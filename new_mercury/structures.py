"""
Vector Data Structures
"""

import numpy as np

class Vector():

    def __init__(self, *args):
        self._vector = np.array(args, dtype = np.ndarray)
        self._size = len(args)
        self._dict = {}

    def pretty_print(self, pretty = True):
        print self._vector
        if len(self._dict) and pretty:
            for key in sorted(self._dict):
                print "| %s = %s |" % (key, str(self._dict[key])),
        print

    def set(self, entry, arg_num = -1):
        raise Exception("Not implemented yet")

    def get(self, name_or_num):
        if isinstance(name_or_num, int):
            num = name_or_num
            return self._vector[num]
        else:
            name = name_or_num
            return self._dict[name] # maybe handle exception also?

    def name_arg(self, name, arg_num):
        if 0 <= arg_num < self._size:
            self._dict[name] = self._vector[arg_num]

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

    def __init__(self, *args):
        if len(args) == 3:
            (x,y,z) = args
            Vector.__init__(self, x, y, z)
            self.name_arg('x', 0)
            self.name_arg('y', 1)
            self.name_arg('z', 2)
        else:
            Vector.__init__(self)

    def set_x(self, x):
        self.set(x, arg_num = 0)
        self.name_arg('x', 0)

    def set_y(self, y):
        self.set(y, arg_num = 1)
        self.name_arg('y', 1)

    def set_z(self, z):
        self.set(z, arg_num = 2)
        self.name_arg('z', 2)

class Velocity(Vector):

    def __init__(self, *args):
        if len(args) == 3:
            (vx,vy,vz) = args
            Vector.__init__(self, vx, vy, vz)
            self.name_arg('vx', 0)
            self.name_arg('vy', 1)
            self.name_arg('vz', 2)
        else:
            Vector.__init__(self)

    def set_vx(self, vx):
        self.set(x, arg_num = 0)
        self.name_arg('vx', 0)

    def set_vy(self, vy):
        self.set(y, arg_num = 1)
        self.name_arg('vy', 1)

    def set_vz(self, vz):
        self.set(z, arg_num = 2)
        self.name_arg('vz', 2)

class AngularMomentum(Vector):
    """ specific angular momentum """

    def __init__(self, *args):
        if len(args) == 3:
            (hx, hy, hz) = args
            Vector.__init__(self, hx, hy, hz)
            self.name_arg('hx', 0)
            self.name_arg('hy', 1)
            self.name_arg('hz', 2)
        else:
            Vector.__init__(self)

    def set_hx(self, hx):
        self.set(x, arg_num = 0)
        self.name_arg('hx', 0)

    def set_hy(self, hy):
        self.set(y, arg_num = 1)
        self.name_arg('hy', 1)

    def set_hz(self, hz):
        self.set(hz, arg_num = 2)
        self.name_arg('hz', 2)


        

