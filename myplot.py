import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mc
import datahandle as dh

# plotting routine
def mycontour(data, x, y, xmax=None, xmin=None, ymax=None, ymin=None, cmin=None, cmax=None, cmap='jet', \
              title=None, xlabel=None, ylabel=None, logscale=False):
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

def sliceplt(x, data, cut, xmax=None, xmin=None, ymax=None, ymin=None, xcut=True, \
             title=None, xlabel=None, ylabel=None):
    plt.clf()
    if xcut:
        y = data[cut,:]
    else:
        y = data[:,cut]
    plt.plot(x, y)
    putlabel(title, xlabel, ylabel)
    setlimit(xmin, xmax, ymin, ymax)
    plt.tight_layout()

def ptl_hist2d_plot(job, time, mx, my, xbins, ybins, xmin, xmax, ymin, ymax, dim=1, \
                    pui=False, cmap='jet',title=None,xlabel=None,ylabel=None, cmin=None, cmax=None, \
                    logscale=False):
    X, Y, h = dh.make_ptl_hist(job, time, mx, my, xbins, ybins, xmin, xmax, ymin, ymax, dim=dim, pui=pui)
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
    X, Y, ibex = dh.makeIBEX(job,time,xbins,ybins,emin,emax,ymin,ymax)
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

def mylineplot(x, y, lw=.5, cls=False, label=None, ls='solid', col=None):
    if cls:
        plt.clf()
    plt.plot(x, y, lw=lw, label=label, linestyle=ls, color=col)
        
# utility
def putlabel(title, xlabel, ylabel):
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
def setlimit(xmin, xmax, ymin, ymax):
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
