import numpy as np
import myplot as myp
import matplotlib.pyplot as plt

##### General data plot #####
def idfuncplot(job, rmin, rmax, time=None, quant='ke', ndiv=200, \
               xr=True, yr=True, hlog=False, pui=False, xlog=False, ylog=False, cls=False):
    if cls:
        plt.clf()
    sampling_xmin, sampling_xmax, sampling_ymin, sampling_ymax = None, None, None, None
    if time==None:
        time = job.get_time(job.nt)
    if xr:
        print('xmin (>0), xmax (<%f):' % (job.nx*job.dx))
        sampling_xmin, sampling_xmax = [float(i) for i in input().split()]
    if job.dim==2:
        if yr:
            print('ymin (>0), ymax (<%f):' % (job.ny*job.dx))
            sampling_ymin, sampling_ymax = [float(i) for i in input().split()]

    x, y = job.distfunc(time, quant, ndiv, rmin, rmax, \
                        sampling_xmin=sampling_xmin, sampling_xmax=sampling_xmax, \
                        sampling_ymin=sampling_ymin, sampling_ymax=sampling_ymax, \
                        hlog=hlog, pui=pui)
    plt.step(x, y)
    if xlog:
        plt.xscale('log')
    if ylog:
        plt.yscale('log')


##### 1D data plot #####
def prof1dplt(job, time=None, quant='np', cls=True):
    if time==None:
        time = job.get_time(job.ns)
    df = job.fld(time=time)
    x = df.x
    title = 'T='+str(df.time)

    y, ylabel = job.fldquant(df, quant)

    myp.mylineplot(x, y, cls=cls)
    myp.putlabel(title, 'x', ylabel)

def basic1dplt(job, time=None, cls=True, xmin=None, xmax=None, pui=False):
    if cls:
        plt.clf()
    if time==None:
        time = job.get_time(job.ns)
    df = job.fld(time=time)
    title = 'T=' + str(time)
    if pui:
        pnum = 5
    else:
        pnum = 3
    plt.subplot(pnum, 1, 1)
    myp.mylineplot(df.x, df.bf, xmin=xmin, xmax=xmax)
    plt.title(title, fontsize=18)
    plt.ylabel(r'$|B|$', fontsize=18)
    plt.subplot(pnum, 1, 2)
    myp.mylineplot(df.x, df.np, xmin=xmin, xmax=xmax)
    plt.ylabel(r'$N_p$', fontsize=18)
    plt.subplot(pnum, 1, 3)
    myp.mylineplot(df.x, df.pp, xmin=xmin, xmax=xmax)
    myp.mylineplot(df.x, df.pl)
    plt.ylabel(r'$P$', fontsize=18)
    if pui:
        plt.subplot(pnum, 1, 4)
        myp.mylineplot(df.x, df.npui, xmin=xmin, xmax=xmax)
        plt.ylabel(r'$N_{PUI}$', fontsize=18)
        plt.subplot(pnum, 1, 5)
        myp.mylineplot(df.x, df.ppui, xmin=xmin, xmax=xmax)
        myp.mylineplot(df.x, df.prui)
        plt.ylabel(r'$P_{PUI}$', fontsize=18)
    plt.xlabel('$x$')
    plt.tight_layout()

def ialignedplot(jlist, time=None, quant='np', yr=False, cls=True):
    if cls:
        plt.clf()
    if time==None:
        time = jlist[0].get_time(jlist[0].ns)
    jnum = len(jlist)
    d = []
    for i, j in enumerate(jlist):
        df = j.fld(time=time, label=quant)
        d.append(df)
    print("xmin (>%f), xmax (<%f):" % (min(d[0].x), max(d[0].x)))
    xmin, xmax = [float(i) for i in input().split()]
    if yr:
        print("ymin, ymax:")
        ymin, ymax = [float(i) for i in input().split()]
    for i, di in enumerate(d):
        plt.subplot(jnum, 1, i+1)
        plt.plot(di.x, di.dat, lw=.5)
        plt.title(di.id + ' : ' + di.label)
        plt.xlim(xmin, xmax)
        if yr:
            plt.ylim(ymin, ymax)
    plt.tight_layout()

##### 2D data plot #####
def icontourmap(job, time=None, quant='np', xr=True, yr=True, cr=False, cmap='jet', logscale=False):
    if time==None:
        time = job.get_time(job.ns)
    df = job.fld(time=time)
    x = df.x
    y = df.y
    if xr:
        print("xmin (>%f), xmax (<%f):" % (min(df.x), max(df.x)))
        xmin, xmax = [float(i) for i in input().split()]
    else:
        xmin = None
        xmax= None

    if yr:
        print("ymin (>%f), ymax (<%f):" % (min(df.y), max(df.y)))
        ymin, ymax = [float(i) for i in input().split()]
    else:
        ymin = None
        ymax = None
        
    if cr:
        print("cmin, cmax")
        cmin, cmax = [float(i) for i in input().split()]
    else:
        cmin = None
        cmax = None
        
    z, title = job.fldquant(df, quant)

    myp.mycontour(z, x, y, \
                  xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin, \
                  cmin=cmin, cmax=cmax, cmap=cmap, \
                  title=title, xlabel='x', ylabel='y', logscale=logscale)

