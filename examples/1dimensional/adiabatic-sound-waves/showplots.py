#!/usr/bin/python3

import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

gamma = 1.4

for time in range(0,101):
  print('Time: '+ str(time))
  extStr = '_' + str(time ) + '.h5'

  fRho = h5py.File('Rho' + extStr, 'r')
  fMx = h5py.File('Mx' + extStr, 'r')

  rho = fRho['data'][()]
  mx = fMx['data'][()]

  x = 5*np.arange(rho.shape[0]) / rho.shape[0]

  vx = mx/rho
  pressure = rho**gamma
  temp = pressure/rho

  plt.figure(figsize=(16,9))
  plt.subplot(2, 2, 1)
  lnRho = plt.plot(x, rho, 'k-')
  axes = plt.gca()
  plt.xlim(0, 5)
  plt.ylim(0.99, 1.01)
  axes.set_ylabel('$\\rho$', fontsize=20)
  axes.tick_params(axis='both', which='major', labelsize=14)

  plt.subplot(2, 2, 2)
  lnMx = plt.plot(x, temp, 'k-')
  axes = plt.gca()
  plt.xlim(0, 5)
  plt.ylim(0.998, 1.002)
  axes.set_ylabel('$T$', fontsize=20)
  axes.tick_params(axis='both', which='major', labelsize=14)

  plt.subplot(2, 2, 3)
  lnMx = plt.plot(x, mx, 'k-')
  axes = plt.gca()
  plt.xlim(0, 5)
  plt.ylim(-0.002, 0.002)
  axes.set_ylabel('$\\rho v$', fontsize=20)
  axes.tick_params(axis='both', which='major', labelsize=14)

  plt.subplot(2, 2, 4)
  lnMx = plt.plot(x, vx, 'k-')
  axes = plt.gca()
  plt.xlim(0, 5)
  plt.ylim(-0.002, 0.002)
  axes.set_ylabel('$v$', fontsize=20)
  axes.tick_params(axis='both', which='major', labelsize=14)

  # plt.show()
  plt.savefig('plots/time'+str(time)+'.png', dpi=120)
  plt.close()