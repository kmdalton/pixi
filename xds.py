#!/usr/bin/env python
import numpy as np



class xparm(dict):
  def __init__(self, inFN):
    self.inFN = inFN
    lines = open(inFN, 'r').readlines()
    #self.parms = None
    while len(lines) > 5:
        p = [lines.pop(5) for i in range(9)]
        self[p[0].strip()] = parm(p)
    self.header = lines

  def align_parms(self):
    for v in self.values():
        if v.sign == '-':
            v.flip_axes()

  def write(self, outFN):
    with open(outFN, 'w') as out:
      out.writelines(self.header)
      for k in sorted(self):
        out.writelines(self[k].lines)

class parm():
  def __init__(self, lines):
    #TODO: throw an error if len(lines) is wrong
    self.lines = lines
    self.A = np.array(map(float, lines[2].split()))
    self.B = np.array(map(float, lines[3].split()))
    self.C = np.array(map(float, lines[4].split()))
    self.vert_axis_name = None
    self.sign           = None
    self.degrees        = None
    self._findvertical__()

  def _findvertical__(self):
    a  = self.A/np.linalg.norm(self.A)
    b  = self.B/np.linalg.norm(self.B)
    c  = self.C/np.linalg.norm(self.C)

    pa = 180.*np.arccos(np.dot(a, [0., 1., 0.]))/np.pi
    pb = 180.*np.arccos(np.dot(b, [0., 1., 0.]))/np.pi
    pc = 180.*np.arccos(np.dot(c, [0., 1., 0.]))/np.pi

    na = np.abs(pa-180.)
    nb = np.abs(pb-180.)
    nc = np.abs(pc-180.)

    sign = np.argmin([min(pa,pb,pc), min(na,nb,nc)])
    ax   = np.argmin([min(pa,na), min(pb,nb), min(pc,nc)])
    self.vert_axis_name = ['A', 'B', 'C'][ax]
    self.sign = ['+', '-'][sign]
    self.degrees = [pa, pb, pc][ax]

  def flip_axes(self):
    #X-->-X, Y-->-Y, Z-->Z
    self.A = self.A*[-1., -1., 1.]
    self.B = self.B*[-1., -1., 1.]
    self.C = self.C*[-1., -1., 1.]
    self._findvertical__()
    self._update_lines()

  def _update_lines(self):
    fun = lambda x: "{: 0.6f}E{:+03d}".format(x/10**(np.floor(np.log10(np.abs(x)))+1), int(np.floor(np.log10(np.abs(x))) +1))
    self.lines[2] = " "+" ".join(map(fun, self.A)) + "\n"
    self.lines[3] = " "+" ".join(map(fun, self.B)) + "\n"
    self.lines[4] = " "+" ".join(map(fun, self.C)) + "\n"

  def __str__(self):
    return """parm instance:
      A: {}
      B: {}
      C: {}
      Axis {}: {} degrees to vertical""".format(self.A, self.B, self.C, self.vert_axis_name, self.degrees)
