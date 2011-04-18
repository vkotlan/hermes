#!/usr/bin/python

# import libraries
import pylab as pl
import numpy as np

time_stop = 70
time_step = 1

markers = ['o', 'v', '^', '*', 's', 'd', 'h', 'p', '<', '>']

pl.close()

zr = []
r = []
xT = []
rt = []
T = []
index = 1
for i in range(time_step, time_stop*2+time_step, time_step*1): 
    data = np.loadtxt("chart_" + str(i) + ".dat")
    zr.append(data[:, 0])
    r.append(data[:, 1])

    # marker = markers[index]
    index += 1
      
    data = np.loadtxt("chart_temp_" + str(i) + ".dat")
    xT.append(data[:, 0])
    T.append(data[:, 1])

    # marker = markers[index]
    # index += 1


pl.close()
pl.figure(figsize=[8, 5.5])
pl.xlabel("$t\,\mathrm{(s)}$")
# pl.ylabel("$T\,\mathrm{(}^{\circ}\mathrm{C}\mathrm{)}$")
pl.ylabel("$u_r~\mathrm{(m)}$")
ax = pl.gca()
ax.ticklabel_format(style='sci', scilimits=(0,0), axis='y') 
    
# min. potrebny posuv
pl.plot([0, 70], [19e-6, 19e-6], '--r', linewidth=3)
    
index = 1
for i in range(1, len(xT[1]), 5):
    xTt = []
    Tt = []
    rt = []
    for j in range(1, len(xT)):
        xTt.append(j*time_step/2.0)
        Tt.append(T[j][i])
        rt.append(r[j][i])

    marker = markers[index - 8]
    index += 1
    
    # pl.plot(xTt, Tt, marker=marker, markevery=3, color='black', label="$"+ ("%5.3f" % xT[1][i]) + "\,\mathrm{m}$")
    # pl.plot(xTt, Tt, '-', label="$"+ ("%5.3f" % xT[1][i]) + "\,\mathrm{m}$")
    # pl.plot(xTt, rt, '-k', marker=marker, markevery=3, label="$"+ ("%5.3f" % xT[1][i]) + "\,\mathrm{m}$")
    pl.plot(xTt, rt, '-', label="$"+ ("%5.3f" % xT[1][i]) + "\,\mathrm{m}$")
    
pl.legend()

# pl.plot([14, 14], [0, 300], '--k', linewidth=1)
pl.plot([14, 14], [0, 3e-5], '--k', linewidth=1)
pl.ylim([0, 3e-5])
pl.text(5.5, 1.2e-6, "$\mathrm{ohrev}$", fontsize="20")
pl.text(15.9, 1.2e-6, "$\mathrm{chlazeni}$", fontsize="20")
# pl.text(5.5, 10, "$\mathrm{ohrev}$", fontsize="20")
# pl.text(15.9, 10, "$\mathrm{chlazeni}$", fontsize="20")