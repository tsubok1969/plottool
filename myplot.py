import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mc
import datahandle as dh

# plotting routine
def mycontour(data, x, y, xmax=None, xmin=None, ymax=None, ymin=None, cmap='jet', title=None, xlabel=None, ylabel=None):
    X, Y = np.meshgrid(x, y)
    plt.clf()
    plt.pcolormesh(X, Y, data, cmap=cmap, rasterized=True)
# plt.imshow(data[::-1], aspect='auto', extent=(min(x), max(x), min(y), max(y))    
    putlabel(title, xlabel, ylabel)
    setlimit(xmin, xmax, ymin, ymax)
    plt.colorbar()

def sliceplt(x, data, cut, xmax=None, xmin=None, ymax=None, ymin=None, xcut=True, title=None, xlabel=None, ylabel=None):
    plt.clf()
    if xcut:
        y = data[cut,:]
    else:
        y = data[:,cut]
    plt.plot(x, y)
    putlabel(title, xlabel, ylabel)
    setlimit(xmin, xmax, ymin, ymax)

def ptl_hist2d(job, time, mx, my, xbins, ybins, xmin, xmax, ymin, ymax, \
                   pui=False, cmap='jet',title=None,xlabel=None,ylabel=None, \
                   logscale=False):
    X, Y, h = make_ptl_hist2d_whole(job, time, mx, my, xbins, ybins, xmin, xmax, ymin, ymax, job.endian, pui)
    plt.clf()
    if logscale:
        norm = mc.LogNorm()
    else:
        norm = None
    plt.pcolormesh(X, Y, h, cmap=cmap, norm=norm, rasterized=True)
    plt.colorbar()
    putlabel(title, xlabel, ylabel)

def pseudoIBEX(job,time,xbins,ybins,emin,emax,ymin,ymax,cmap='jet',logscale=False,title=None):
    X, Y, ibex = makeIBEX(job,time,xbins,ybins,emin,emax,ymin,ymax,endian=job.endian)
    plt.clf()
    if logscale:
        norm = mc.LogNorm()
    else:
        norm = None
    plt.pcolormesh(X, Y, ibex, cmap=cmap, norm=norm, rasterized=True)
    plt.colorbar()
    putlabel(title, 'x', 'Energy')

def multiIBEX(job,ts,te,xbins,ybins,emin,emax,ymin,ymax,cmap='jet',logscale=False):
    tnum = te-ts+1
    for i in range(tnum):
        index = ts+i
        title = 'Time: ' + str(index*job.tp)
        pseudoIBEX(job,index,xbins,ybins,emin,emax,ymin,ymax,cmap,logscale,title)
        fname = 'p%05.f'%(index) + '.png'
        print("Time:"+str(index))
        plt.savefig(fname)

def makeanime(job,ts,te,iph=1,xmin=None,xmax=None,ymin=None,ymax=None,cmap='jet',pui=True):
    tnum = te-ts+1
    for i in range(tnum):
        index = ts+i
        plt.clf()
        df = dh.fld2d(job,index,pui=pui)
        X, Y = np.meshgrid(df.x, df.y)
        title = 'Time: ' + str(index*job.tu)
        if iph==1:
            z = df.bx
        elif iph==2:
            z = df.by
        elif iph==3:
            z = df.bz
        elif iph==4:
            z = df.bf
        elif iph==5:
            z = df.np
        elif iph==6:
            z = df.pp
        elif iph==7:
            z = df.pl

        plt.pcolormesh(X, Y, z, cmap=cmap)

        putlabel(title,'x','y')
        print('Time: ' + str(index))

        
# utility
def putlabel(title, xlabel, ylabel):
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
def setlimit(xmin, xmax, ymin, ymax):
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)

def make_ptl_hist2d_whole(job, time, mx, my, xbins, ybins, xmin, xmax, ymin, ymax, \
                            endian = False, pui=False):
# mx, my : 1 -> xp, 2 -> yp, 3 -> vx, 4 -> vy, 5 -> vz, 6 -> ke
    htmp = {}
    for i in range(job.mp):
        dp = dh.ptl2D(job, time, mpi=True, proc=i, pui=pui)
        if mx == 1:
            x = dp.xp
        elif mx == 2:
            x = dp.yp
        elif mx == 3:
            x = dp.vx
        elif mx == 4:
            x = dp.vy
        elif mx == 5:
            x = dp.vz
        elif mx == 6:
            x = dp.ke

        if my == 1:
            y = dp.xp
        elif my == 2:
            y = dp.yp
        elif my == 3:
            y = dp.vx
        elif my == 4:
            y = dp.vy
        elif my == 5:
            y = dp.vz
        elif my == 6:
            y = dp.ke

        htmp[i] = np.histogram2d(x, y, bins=[xbins, ybins], range=[[xmin,xmax],[ymin,ymax]])
        if i%8==0:
            print('Proc:' + str(i) + ' finished')
    X, Y, h = hist2d_mpi(htmp, job.mp)
    return (X, Y, h)

def make_ptl_hist2d_partial(job, time, mx, my, xbins, ybins, xmin, xmax, ymin, ymax, \
                            rng1, rng2, endian = False, pui=False):
# mx, my : 1 -> xp, 2 -> yp, 3 -> vx, 4 -> vy, 5 -> vz, 6 -> ke
# (rng1, rng2) -> y-sampling range
    htmp = {}
    impi = 0
    for i in range(job.mp):
        dp = dh.ptl2D(job, time, mpi=True, proc=i, pui=pui)
        index = np.where( (dp.yp > rng1) & (dp.yp < rng2) )
        if i%8==0:
            print('Proc:' + str(i) + ' finished')
        if dp.xp[index].size == 0:
            continue
        if mx == 1:
            x = dp.xp[index]
        elif mx == 2:
            x = dp.yp[index]
        elif mx == 3:
            x = dp.vx[index]
        elif mx == 4:
            x = dp.vy[index]
        elif mx == 5:
            x = dp.vz[index]
        elif mx == 6:
            x = dp.ke[index]

        if my == 1:
            y = dp.xp[index]
        elif my == 2:
            y = dp.yp[index]
        elif my == 3:
            y = dp.vx[index]
        elif my == 4:
            y = dp.vy[index]
        elif my == 5:
            y = dp.vz[index]
        elif my == 6:
            y = dp.ke[index]

        htmp[impi] = np.histogram2d(x, y, bins=[xbins, ybins], range=[[xmin,xmax],[ymin,ymax]])
        impi += 1
    X, Y, h = hist2d_mpi(htmp, impi)
    return (X, Y, h)

def makeIBEX(job,time,xbins,ybins,emin,emax,ymin,ymax,pui=True,endian=False):
    xmin = 0.0
    xmax = job.nx*job.dx
    htmp = {}
    impi = 0
    for i in range(job.mp):
        dp = dh.ptl2D(job,time,mpi=True,proc=i,pui=pui)
        index = np.where( (dp.yp >= ymin) & (dp.yp <= ymax) )
        if i%8==0:
            print('Proc:' + str(i) + ' complete')
        if dp.xp[index].size == 0:
            continue
        htmp[impi] = np.histogram2d(dp.xp[index],dp.ke[index],bins=[xbins,ybins],range=[[xmin,xmax], [emin,emax]])
        impi += 1
    X, Y, ibex = hist2d_mpi(htmp, impi)
    return (X, Y, ibex)

def hist2d_mpi(h, mp):
    hretn = h[0][0]
    for i in range(mp-1):
        hretn += h[i+1][0]
    X, Y = np.meshgrid(h[0][1], h[0][2])
    return (X, Y, hretn.T)
