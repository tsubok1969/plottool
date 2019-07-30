import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mc
import datahandle as dh

# plotting routine
def mycontour(data, x, y, xmax=None, xmin=None, ymax=None, ymin=None, cmin=None, cmax=None, cmap='jet', title=None, xlabel=None, ylabel=None, logscale=False):
    X, Y = np.meshgrid(x, y)
    plt.clf()
    if logscale:
        norm = mc.LogNorm()
    else:
        norm = None
    plt.pcolormesh(X, Y, data, cmap=cmap, norm=norm, rasterized=True)
# plt.imshow(data[::-1], aspect='auto', extent=(min(x), max(x), min(y), max(y))    
    putlabel(title, xlabel, ylabel)
    setlimit(xmin, xmax, ymin, ymax)
    plt.colorbar()
    plt.clim(cmin, cmax)
    plt.tight_layout()

def sliceplt(x, data, cut, xmax=None, xmin=None, ymax=None, ymin=None, xcut=True, title=None, xlabel=None, ylabel=None):
    plt.clf()
    if xcut:
        y = data[cut,:]
    else:
        y = data[:,cut]
    plt.plot(x, y)
    putlabel(title, xlabel, ylabel)
    setlimit(xmin, xmax, ymin, ymax)
    plt.tight_layout()

def ptl_hist2d(job, time, mx, my, xbins, ybins, xmin, xmax, ymin, ymax, \
                   pui=False, cmap='jet',title=None,xlabel=None,ylabel=None, cmin=None, cmax=None, \
                   logscale=False):
    X, Y, h = make_ptl_hist2d_whole(job, time, mx, my, xbins, ybins, xmin, xmax, ymin, ymax, job.endian, pui)
    plt.clf()
    if logscale:
        norm = mc.LogNorm()
    else:
        norm = None
    plt.pcolormesh(X, Y, h, cmap=cmap, norm=norm, rasterized=True)
    plt.colorbar()
    plt.clim(cmin, cmax)
    putlabel(title, xlabel, ylabel)
    plt.tight_layout()

def pseudoIBEX(job,time,xbins,ybins,emin,emax,ymin,ymax,cmap='jet',logscale=False,title=None,cmin=None,cmax=None):
    X, Y, ibex = makeIBEX(job,time,xbins,ybins,emin,emax,ymin,ymax,endian=job.endian)
    plt.clf()
    if logscale:
        norm = mc.LogNorm()
    else:
        norm = None
    plt.pcolormesh(X, Y, ibex, cmap=cmap, norm=norm, rasterized=True)
    plt.colorbar()
    plt.clim(cmin,cmax)
    putlabel(title, 'x', 'Energy')
    plt.tight_layout()

def multiIBEX(job,ts,te,xbins,ybins,emin,emax,ymin,ymax,cmap='jet',logscale=False,cmin=None,cmax=None):
    tnum = te-ts+1
    for i in range(tnum):
        index = ts+i
        title = 'Time: ' + str(index*job.tp)
        pseudoIBEX(job,index,xbins,ybins,emin,emax,ymin,ymax,cmap,logscale,title,cmin,cmax)
#        fname = '~/tmp/figure/p%05.f'%(index) + '.png'
        fname = './figure/p%05.f'%(index) + '.png'
        print("Time:"+str(index))
        plt.tight_layout()
        plt.savefig(fname)

def makeanime2D(job,ts,te,quant,xmin=None,xmax=None,ymin=None,ymax=None,cmin=None,cmax=None,cmap='jet',pui=True):
    tnum = te-ts+1
    for i in range(tnum):
        index = ts+i
        plt.clf()
        df = dh.fld2D(job,index,pui=pui)
        X, Y = np.meshgrid(df.x, df.y)
        if quant=='bx':
            z = df.bx
            title_arg = 'Bx'
        elif quant=='by':
            z = df.by
            title_arg = 'By'
        elif quant=='bz':
            z = df.bz
            title_arg = 'Bz'
        elif quant=='bf':
            z = df.bf
            title_arg = '|B|'
        elif quant=='np':
            z = df.np
            title_arg = r'$N_p$'
        elif quant=='pp':
            z = df.pp
            title_arg = r'$P_\perp$'
        elif quant=='pl':
            z = df.pl
            title_arg = r'$P_\parallel$'
        elif quant=='npui':
            z = df.npui
            title_arg = r'$N_{PUI}$'
        elif quant=='ppui':
            z = df.ppui
            title_arg = r'$P_{\perp, PUI}$'
        elif quant=='prui':
            z = df.prui
            title_arg = r'$P_{\parallel, PUI}$'

        title = title_arg + ' at Time: ' + str(index*job.tu)

        plt.pcolormesh(X, Y, z, cmap=cmap)
        plt.colorbar()
        plt.clim(cmin,cmax)

        setlimit(xmin,xmax,ymin,ymax)
        putlabel(title,'x','y')
        print('Time: ' + str(index))
        plt.tight_layout()
        fname = './figure/f%05.f'%(index) + '.png'
        plt.savefig(fname)
        
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
        if mx == 'xp':
            x = dp.xp
        elif mx == 'yp':
            x = dp.yp
        elif mx == 'vx':
            x = dp.vx
        elif mx == 'vy':
            x = dp.vy
        elif mx == 'vz':
            x = dp.vz
        elif mx == 'ke':
            x = dp.ke

        if my == 'xp':
            y = dp.xp
        elif my == 'yp':
            y = dp.yp
        elif my == 'vx':
            y = dp.vx
        elif my == 'vy':
            y = dp.vy
        elif my == 'vz':
            y = dp.vz
        elif my == 'ke':
            y = dp.ke

        htmp[i] = np.histogram2d(x, y, bins=[xbins, ybins], range=[[xmin,xmax],[ymin,ymax]])
        if i%8==0:
            print('Proc:' + str(i) + ' finished')
    X, Y, h = hist2d_mpi(htmp, job.mp)
    return (X, Y, h)

def make_ptl_hist2d_partial(job, time, mx, my, xbins, ybins, xmin, xmax, ymin, ymax, \
                            sampling_ymin=None, sampling_ymax=None, sampling_xmin=None, sampling_xmax=None, endian = False, pui=False):
# (yrng1, yrng2) -> y-sampling range
    htmp = {}
    impi = 0
    for i in range(job.mp):
        dp = dh.ptl2D(job, time, mpi=True, proc=i, pui=pui)
        if sampling_xmin==None:
            sampling_xmin = np.min(dp.xp)
        if sampling_xmax==None:
            sampling_xmax = np.max(dp.xp)
        if sampling_ymin==None:
            sampling_ymin = np.min(dp.yp)
        if sampling_ymax==None:
            sampling_ymax = np.max(dp.yp)
        index = np.where( (dp.yp >= sampling_ymin) & (dp.yp <= sampling_ymax) & (dp.xp >= sampling_xmin) & (dp.xp <= sampling_xmax) )

        if i%8==0:
            print('Proc:' + str(i) + ' finished')
        if dp.xp[index].size == 0:
            continue
        if mx == 'xp':
            x = dp.xp[index]
        elif mx == 'yp':
            x = dp.yp[index]
        elif mx == 'vx':
            x = dp.vx[index]
        elif mx == 'vy':
            x = dp.vy[index]
        elif mx == 'vz':
            x = dp.vz[index]
        elif mx == 'ke':
            x = dp.ke[index]

        if my == 'xp':
            y = dp.xp[index]
        elif my == 'yp':
            y = dp.yp[index]
        elif my == 'vx':
            y = dp.vx[index]
        elif my == 'vy':
            y = dp.vy[index]
        elif my == 'vz':
            y = dp.vz[index]
        elif my == 'ke':
            y = dp.ke[index]

        htmp[impi] = np.histogram2d(x, y, bins=[xbins, ybins], range=[[xmin,xmax],[ymin,ymax]])
        impi += 1
    X, Y, h = hist2d_mpi(htmp, impi)
    return (X, Y, h)

def make_ptl_hist1d(job, time, mx, my, xbins, ybins, xmin, xmax, ymin, ymax, \
                    sampling_xmin=None, sampling_xmax=None, endian = False, pui=False):
    htmp = {}
    impi = 0
    for i in range(job.mp):
        dp = dh.ptl1D(job, time, mpi=True, proc=i, pui=pui)
        if sampling_xmin==None:
            sampling_xmin = np.min(dp.xp)
        if sampling_xmax==None:
            sampling_xmax = np.max(dp.xp)
        index = np.where( (dp.xp >= sampling_xmin) & (dp.xp <= sampling_xmax) )

        if i%8==0:
            print('Proc:' + str(i) + ' finished')
        if dp.xp[index].size == 0:
            continue

        if mx == 'xp':
            x = dp.xp[index]
        elif mx == 'vx':
            x = dp.vx[index]
        elif mx == 'vy':
            x = dp.vy[index]
        elif mx == 'vz':
            x = dp.vz[index]
        elif mx == 'ke':
            x = dp.ke[index]

        if my == 'xp':
            y = dp.xp[index]
        elif my == 'vx':
            y = dp.vx[index]
        elif my == 'vy':
            y = dp.vy[index]
        elif my == 'vz':
            y = dp.vz[index]
        elif my == 'ke':
            y = dp.ke[index]

        htmp[impi] = np.histogram2d(x, y, bins=[xbins, ybins], range=[[xmin,xmax],[ymin,ymax]])
        impi += 1

    if impi==0:
        print("there are no data in the sampling box")
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

def distfunc(job, time, quantity, ndiv, rmin, rmax, sampling_ymin, sampling_ymax, hlog=False, pui=False):
    ## making the distribution of the particle velocity or energy
# quantity: assigned by the quantity index, e.g., 'vx', 'xp', 'ke', etc
    # ndiv: dividing number for the histogram
    # rmin, rmax: minimum and maximum ranges of the histogram
    # hlog: the range to be divided by the logarithmic scale
    numproc = 0
    dtmp = {}
    impi = 0
    x = np.arange(ndiv+1,dtype=np.float64)
    if hlog:
        i1 = np.log10(rmin)
        i2 = np.log10(rmax)
        x = 10.**(i1 + x*(i2-i1)/ndiv)
    else:
        i1 = rmin
        i2 = rmax
        x = x*(i2-i1)/ndiv + i1

    for i in range(job.mp):
        dp = dh.ptl2D(job, time, mpi=True, proc=i, pui=pui)
        index = np.where( (dp.yp > sampling_ymin) & (dp.yp < sampling_ymax) )
        if i%8 == 0:
            print('Proc:' + str(i) + ' finished')
        if dp.xp[index].size == 0:
            continue
        
        if quantity == 'vx':
            v = dp.vx[index]
        elif quantity == 'vy':
            v = dp.vy[index]
        elif quantity == 'vz':
            v = dp.vz[index]
        elif quantity == 'v':
            v = np.sqrt(dp.vx[index]**2 + dp.vy[index]**2 + dp.vz[index]**2)
        else:
            v = dp.ke[index]

        if hlog:
            v = np.log10(v)

        dtmp[impi] = np.histogram(v, bins=ndiv, range=[i1, i2])
        impi += 1
        numproc += v.size

    cum_dtmp = dtmp[0][0]
    for i in range(impi-1):
        cum_dtmp += dtmp[i+1][0]
    print(float(numproc))
    return(x[0:ndiv], cum_dtmp.astype(np.float64)/numproc)
