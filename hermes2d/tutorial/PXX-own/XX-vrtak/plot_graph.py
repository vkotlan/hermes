#!/usr/bin/python

# import libraries
import pylab as pl
import numpy as np

time_stop = 14
time_step = 4

markers = ['o', 'v', '^', '*', 's', 'd', 'h', 'p', '<', '>']

pl.close()

# min. potrebny posuv
pl.plot([0, 0.0135], [19e-6, 19e-6], '--r', linewidth=3)

zr = []
r = []
xT = []
T = []
index = 1
for i in range(time_step, time_stop*2+time_step, time_step*1): 
    data = np.loadtxt("chart_" + str(i) + ".dat")
    zr.append(data[:, 0])
    r.append(data[:, 1])

    marker = markers[index]
    index += 1
      
    # pl.subplot(2, 1, 1)
    # pl.plot(zr[-1], r[-1], marker=marker, markevery=3, color='k', label="$"+ ("%03d" % (i/2)) + "\,\mathrm{s}$")
    pl.plot(zr[-1], r[-1], markevery=1, label="$"+ ("%02d" % (i/2)) + "\,\mathrm{s}$")
   
    data = np.loadtxt("chart_temp_" + str(i) + ".dat")
    xT.append(data[:, 0])
    T.append(data[:, 1])

    # marker = markers[index]
    # index += 1
        
    # pl.subplot(2, 1, 2)
    # pl.plot(x[-1], y[-1], marker=marker, markevery=3, color='k', label="$"+ ("%03d" % (i/2)) + "\,\mathrm{s}$")
    # pl.plot(xT[-1], T[-1], markevery=3, label="$"+ ("%03d" % (i/2)) + "\,\mathrm{s}$")

# pl.subplot(2, 1, 1)


pl.xlabel("$z~\mathrm{(m)}$")
pl.ylabel("$u_r~\mathrm{(m)}$")
pl.grid(1)
pl.xlim([0, 1.6e-2])
pl.legend(loc = 'lower right')
ax = pl.gca()
ax.ticklabel_format(style='sci', scilimits=(0,0), axis='x') 
ax.ticklabel_format(style='sci', scilimits=(0,0), axis='y') 
"""
pl.subplot(2, 1, 2)
pl.xlabel("$z~\mathrm{(m)}$")
pl.ylabel("$T\,\mathrm{(}^{\circ}\mathrm{C}\mathrm{)}$")
pl.grid(1)
pl.legend()
"""
# pl.savefig("teplota_delka_bw.pdf")

"""
pl.close()
pl.figure(figsize=[8, 5.5])
pl.xlabel("$t\,\mathrm{(s)}$")
pl.ylabel("$T\,\mathrm{(}^{\circ}\mathrm{C}\mathrm{)}$")

index = 1
for i in range(1, len(xT[1]), 5):
    xTt = []
    Tt = []
    for j in range(1, len(xT)):
        xTt.append(j*time_step/2.0)
        Tt.append(T[j][i])

    # marker = markers[index - 8]
    index += 1
    
    # pl.plot(xt, yt, marker=marker, markevery=3, color='black', label="$"+ ("%5.3f" % x[1][i]) + "\,\mathrm{m}$")
    pl.plot(xTt, Tt, '-*', markevery=3, label="$"+ ("%5.3f" % xT[1][i]) + "\,\mathrm{m}$")

pl.legend()
"""