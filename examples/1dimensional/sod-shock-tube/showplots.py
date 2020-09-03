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

  rho = fRho['data0'][()]
  mx = fMx['data0'][()]
  eng = fE['data0'][()]

  x = np.arange(rho.shape[0]) / rho.shape[0]

  vx = mx/rho
  pressure = (eng - 0.5*mx*vx)/0.4
  temp = pressure/rho

  plt.figure(figsize=(16,9))
  plt.subplot(3, 2, 1)
  lnRho = plt.plot(x, rho, 'k-')
  axes = plt.gca()
  plt.xlim(0, 1)
  plt.ylim(0, 1.05)
  axes.set_ylabel('$\\rho$', fontsize=20)
  axes.tick_params(axis='both', which='major', labelsize=14)

  plt.subplot(3, 2, 2)
  lnMx = plt.plot(x, temp, 'k-')
  axes = plt.gca()
  plt.xlim(0, 1)
  plt.ylim(3.8, 8)
  axes.set_ylabel('$T$', fontsize=20)
  axes.tick_params(axis='both', which='major', labelsize=14)

  plt.subplot(3, 2, 3)
  lnMx = plt.plot(x, mx, 'k-')
  axes = plt.gca()
  plt.xlim(0, 1)
  plt.ylim(-0.01, 0.41)
  axes.set_ylabel('$\\rho v$', fontsize=20)
  axes.tick_params(axis='both', which='major', labelsize=14)

  plt.subplot(3, 2, 4)
  lnMx = plt.plot(x, vx, 'k-')
  axes = plt.gca()
  plt.xlim(0, 1)
  plt.ylim(-0.01, 1.2)
  axes.set_ylabel('$v$', fontsize=20)
  axes.tick_params(axis='both', which='major', labelsize=14)

  plt.subplot(3, 2, 5)
  lnE = plt.plot(x, eng, 'k-')
  axes = plt.gca()
  plt.xlim(0, 1)
  plt.ylim(0, 2.7)
  axes.set_xlabel('$x$', fontsize=20)
  axes.set_ylabel('$E$', fontsize=20)
  axes.tick_params(axis='both', which='major', labelsize=14)

  plt.subplot(3, 2, 6)
  lnMx = plt.plot(x, pressure, 'k-')
  axes = plt.gca()
  plt.xlim(0, 1)
  plt.ylim(0, 7)
  axes.set_ylabel('$p$', fontsize=20)
  axes.tick_params(axis='both', which='major', labelsize=14)

  # plt.show()
  plt.savefig('plots/time'+str(time)+'.png', dpi=120)
  plt.close()