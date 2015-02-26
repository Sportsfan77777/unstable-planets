"""
Vector Data Structures
"""

import numpy as np
import collections

class Vector():

    def __init__(self, *args):
        if len(args) == 1 and isinstance(args, collections.Iterable):
            # Unpack Array or Tuple
            self._vector = np.array(args)
        else:
            self._vector = np.array(args, dtype = np.ndarray)
        self._size = len(self._vector)
        self._dict = {}

    def pretty_print(self, pretty = True):
        """ prints vector and pretty prints vector"""
        print self._vector
        if len(self._dict) and pretty:
            for key in sorted(self._dict):
                print "| %s = %s |" % (key, str(self._dict[key])),
            print

    def size(self):
        return self._size

    def set(self, entry, arg_num = -1):
        """ set specific entry (call by number) """
        if arg_num < 0:
            raise Exception("Must specify arg_num!")
        if arg_num >= len(self._vector):
            new_vector = np.zeros(arg_num, dtype = np.ndarray)
            new_vector[:len(self._vector)] = self._vector[:]
            self._vector = new_vector
            self._size = len(self._vector)

        # Check if arg_num is a named_arg

        self._vector[arg_num] = entry

    def get(self, name_or_num):
        """ call vector entry by name or number """
        if isinstance(name_or_num, int):
            num = name_or_num
            return self._vector[num]
        else:
            name = name_or_num
            return self.get(self._dict[name]) # maybe handle exception also?

    def get_vector(self):
        return self._vector

    def trim(self, max_length):
        """ reduces vector to size of 'max_length' """
        if max_length < len(self._vector):
           self._vector = self._vector[:max_length]
           self._size = len(self._vector)

    def name_arg(self, name, arg_num):
        if 0 <= arg_num < self._size:
            self._dict[name] = arg_num

    def dot(self, vector, ret_type = None):
        """ dot product """
        if len(self._vector) == len(vector):
            if ret_type:
                return ret_type(self._vector.dot(vector))
            else:
                return Vector(self._vector.dot(vector))
        else:
            return -1

    def cross(self, vector, ret_type = None):
        """ cross product """
        if self._size == 3:
            x = self.get(1) * vector.get(2) - self.get(2) * vector.get(1)
            y = self.get(2) * vector.get(0) - self.get(0) * vector.get(2)
            z = self.get(0) * vector.get(1) - self.get(1) * vector.get(0)
            if ret_type:
                return ret_type(x, y, z)
            else:
                return Vector(x, y, z)
        else:
            return -1

    def length(self, override = False):
        if self._size == 3 or override:
            return np.sqrt(self.norm_sq(override))
        else:
            return -1

    def norm_sq(self, override = False):
        if self._size == 3 or override:
            x = self._vector[:]
            #self.pretty_print()
            sum_x = 0
            for entry in x:
                sum_x += (entry * entry)
            #return x.dot(x) # <---- doesn't work for old version of python that exo2a has
            return sum_x
        else:
            return -1

    def length_sq(self, override = False):
        return self.norm_sq(override)


class Position(Vector):
    """ specifically, a position vector """

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
        self.set(vx, arg_num = 0)
        self.name_arg('vx', 0)

    def set_vy(self, vy):
        self.set(vy, arg_num = 1)
        self.name_arg('vy', 1)

    def set_vz(self, vz):
        self.set(vz, arg_num = 2)
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
        self.set(hx, arg_num = 0)
        self.name_arg('hx', 0)

    def set_hy(self, hy):
        self.set(hy, arg_num = 1)
        self.name_arg('hy', 1)

    def set_hz(self, hz):
        self.set(hz, arg_num = 2)
        self.name_arg('hz', 2)

class OrbitalElements(Vector):
    """ 3-D orbital elements (6 of them) """

    def __init__(self, *args):
        if len(args) == 6:
            (a, e, i, argw, node, mean_anom) = args
            Vector.__init__(a, e, i, argw, node, mean_anom)
            self.name_arg('a', 0)
            self.name_arg('sma', 0)

            self.name_arg('e', 1)
            self.name_arg('ecc', 1)

            self.name_arg('i', 2)
            self.name_arg('inc', 2)

            self.name_arg('w', 3)
            self.name_arg('argw', 3)

            self.name_arg('W', 4)
            self.name_arg('node', 4)

            self.name_arg('M', 5)
            self.name_arg('mean_anom', 5)
        else:
            Vector.__init__()

    def set_a(self, a):
        self.set(a, arg_num = 0)
        self.name_arg('a', 0)
        self.name_arg('sma', 0)

    def set_sma(self, sma):
        self.set_a(sma)

    def set_e(self, e):
        self.set(e, arg_num = 1)
        self.name_arg('e', 1)
        self.name_arg('ecc', 1)

    def set_ecc(self, ecc):
        self.set_e(ecc)

    def set_i(self, i):
        self.set(i, arg_num = 2)
        self.name_arg('i', 2)
        self.name_arg('inc', 2)

    def set_inc(self, inc):
        self.set_i(inc)

    def set_w(self, w):
        self.set(w, arg_num = 3)
        self.name_arg('w', 3)
        self.name_arg('argw', 3)

    def set_argw(self, argw):
        self.set_w(argw)

    def set_W(self, W):
        self.set(W, arg_num = 4)
        self.name_arg('W', 4)
        self.name_arg('node', 4)

    def set_node(self, node):
        self.set_W(node)

    def set_M(self, M):
        self.set(M, arg_num = 5)
        self.name_arg('M', 5)
        self.name_arg('mean_anom', 5)

    def set_mean_anom(self, mean_anom):
        self.set_M(mean_anom)

    def get_q(self):
        return self.get('a') * (1 - self.get('e'))

    def get_ecc_anom(self):
        raise Exception('not implemented yet...')

    def get_true_anom(self):
        raise Exception('not implemented yet...')        
        

