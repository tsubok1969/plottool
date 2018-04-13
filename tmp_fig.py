from datahandle import *
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick

def makefig(jobid, mpi=True, pui=True, endian=True, npui=False, ymin=0, ymax=800, t1=200, t2=600, t3=1400, vmin=0, vmax=1.):
    job = job2d(jobid, mpi=mpi, pui=pui, endian=endian)

    df = []
    c  = []
    tm = [t1, t2, t3]

    for i in range(3):
        tmp = fld2D(job, tm[i], pui=pui)
        df.append(tmp)

    fig, axes = plt.subplots(figsize=(10,6), nrows=1, ncols=3)
    X, Y = np.meshgrid(df[0].x, df[0].x)

    for i in range(3):
        tmp = ax[i].pcolormesh(X, Y, df[i].np, cmap='jet', vmin=vmin, vmax=vmax, rasterized=True)
        c.append(tmp)
        ax[i].set_ylim(ymin, ymax)
        ax[i].set_xlabel('x')
        ax[i].set_title('$\Omega_p t='+str(tm[i]*job.nt))
    ax[0].set_ylabel('y')
        
    cax = fig.add_axes([0.91, 0.11, 0.01, 0.77])
    cbar = fig.colorbar(c3, cax = cax, orientation='vertical', format=ptick.FuncFormatter(fmt))
    plt.savefig(jobid+'.tiff')

def fmt(x, pos):
   a, b = '{:.1e}'.format(x).split('e')
   b = int(b)
   return r'${}\times 10^{{{}}}$'.format(a, b)
