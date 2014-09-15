"""
a simple leapfrog integrator
(uses cartesian coordinates, not the Hamiltonian approach)
"""

import numpy as np
import random


def advance_x(x_0, v_half, delta_t):
   """ go to the next x """
   return x_0 + v_half * delta_t
   
def advance_v(v_0, a_half, delta_t):
   """ go to the next v """
   # Note: add 0.5 to all number subscripts
   return v_0 + a_half * delta_t
   
   
def run_integration(F, delta_t, time, x_init = 0, v_init = 0):
   """ run an integration given the second derivative of x: F(x) """
   t = 0.0
   n_steps = int(round(time * 1.0 / delta_t, 0))
   print n_steps
   
   x_arr = np.zeros(n_steps + 1)
   v_arr = np.zeros(n_steps + 1) # NOTE: Off by 0.5 steps compared to x_arr
   
   x_arr[0] = x_init # x_0
   v_arr[0] = advance_v(v_init, F(x_init), delta_t / 2.0) # v_0.5
   
   # Replace this with function(s)
   print "Time: %.2f,  Position: %.2f" % (t, x_arr[0])
   print "Time: %.2f,  Velocity: %.2f" % (t + delta_t/2.0, v_arr[0])
   
   for i in range(n_steps):
       t += delta_t # useful for? printing? not input?
       
       x_now = x_arr[i]
       v_plus_half = v_arr[i]
       a_now = F(x_now)
       
       x_arr[i+1] = advance_x(x_now, v_plus_half, delta_t)
       v_arr[i+1] = advance_v(v_plus_half, a_now, delta_t)
       
       # Change to random print (or once every n_steps / something)
       print "Time: %.2f,  Position: %.4f" % (t, x_arr[i+1])
       print "Time: %.2f,  Velocity: %.3f" % (t + delta_t/2.0, v_arr[i+1])


if __name__ in ('__main__'):
   def F(x):
       return 400.0/(x**(1.5))
       
   # Incorporate args
       
   run_integration(F, 0.5, 200000, x_init = 50)
       