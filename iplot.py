import numpy as np
import myplot as myp

##### 1D data plot #####


##### 2D data plot #####
def icontourmap(df, quant='np', xr=True, yr=True, cr=False, cmap='jet', logscale=False):
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
        
    if quant=='bf':
        z = df.bf
        title = '|B|'
    elif quant=='np':
        z = df.np
        title = r'$N_p$'
    elif quant=='npui':
        z = df.npui
        title = r'$N_{PUI}$'
    elif quant=='ux':
        z = df.ux
        title = r'$U_x$'
    elif quant=='uy':
        z = df.uy
        title = r'$U_y$'
    elif quant=='uz':
        z = df.uz
        title = r'$U_z$'
    elif quant=='pp':
        z = df.pp
        title = r'$P_\perp$'
    elif quant=='pl':
        z = df.pl
        title = r'$P_\parallel$'
    elif quant=='pre':
        z = (df.pp*2.+df.pl)/3.
        title = 'P'
    elif quant=='ppui':
        z = (df.ppui*2.+df.prui)/3.
        title = r'$P_{PUI}$'

    myp.mycontour(z, x, y, \
                  xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin, \
                  cmin=cmin, cmax=cmax, cmap=cmap, \
                  title=title, xlabel='x', ylabel='y', logscale=logscale)
