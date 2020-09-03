#!/usr/bin/python3

import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

for time in range(0,201):
  print('Time: '+ str(time))
  extStr = '_' + str(time ) + '.h5'

  fRho = h5py.File('Rho' + extStr, 'r')
  fMx = h5py.File('Mx' + extStr, 'r')
  fE = h5py.File('E' + extStr, 'r')

  rho = fRho['data0']
  mx = fMx['data0']
  eng = fE['data0']

  x = np.arange(rho.shape[0]) / rho.shape[0]

  plt.figure(figsize=(8,10))
  plt.subplot(3, 1, 1)
  lnRho = plt.plot(x, rho, 'k-')
  axes = plt.gca()
  plt.xlim(0, 1)
  plt.ylim(0, 1.05)
  axes.set_ylabel('$\\rho$', fontsize=20)
  axes.tick_params(axis='both', which='major', labelsize=14)

  plt.subplot(3, 1, 2)
  lnMx = plt.plot(x, mx, 'k-')
  axes = plt.gca()
  plt.xlim(0, 1)
  plt.ylim(-0.01, 0.41)
  axes.set_ylabel('$m_x$', fontsize=20)
  axes.tick_params(axis='both', which='major', labelsize=14)

  plt.subplot(3, 1, 3)
  lnE = plt.plot(x, eng, 'k-')
  axes = plt.gca()
  plt.xlim(0, 1)
  plt.ylim(0, 2.7)
  axes.set_xlabel('$x$', fontsize=20)
  axes.set_ylabel('$E$', fontsize=20)
  axes.tick_params(axis='both', which='major', labelsize=14)

  # plt.show()
  plt.savefig('plots/time'+str(time)+'.png')
  plt.close()