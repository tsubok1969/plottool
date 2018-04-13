def makefig(jobid, mpi=True, pui=True, endian=True, npui=False, ymin=0, ymax=800, t1=200, t2=600, t3=1400, vmin=0, vmax=1.):
    from datahandle import *
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ptick
    job = job2d(jobid, mpi=mpi, pui=pui, endian=endian)
    df1 = fld2D(job, t1, pui=pui)
    df2 = fld2D(job, t2, pui=pui)
    df3 = fld2D(job, t3, pui=pui)
    fig, axes = plt.subplots(figsize=(10,4), nrows=1, ncols=3)
    ax1 = axes[0]
    ax2 = axes[1]
    ax3 = axes[2]
    X, Y = np.meshgrid(df1.x, df2.x)
    c1 = ax1.pcolormesh(X, Y, df1.np, cmap='jet', vmin=vmin, vmax=vmax, rasterized=True)
    c2 = ax2.pcolormesh(X, Y, df2.np, cmap='jet', vmin=vmin, vmax=vmax, rasterized=True)
    c3 = ax3.pcolormesh(X, Y, df3.np, cmap='jet', vmin=vmin, vmax=vmax, rasterized=True)
    ax1.set_ylim(ymin, ymax)
    ax2.set_ylim(ymin, ymax)
    ax3.set_ylim(ymin, ymax)
    ax1.set_ylabel('y')
    ax1.set_xlabel('x')
    ax2.set_xlabel('x')
    ax3.set_xlabel('x')
    ax1.set_title('$\Omega t='+str(t1*job.nt))
    ax2.set_title('$\Omega t='+str(t2*job.nt))
    ax3.set_title('$\Omega t='+str(t3*job.nt))
    cax = fig.add_axes([0.91, 0.11, 0.01, 0.77])
    cbar = fig.colorbar(c3, cax = cax, orientation='vertical', format=ptick.FuncFormatter(fmt))
    plt.savefig(jobid+'.tiff')

def fmt(x, pos):
   a, b = '{:.1e}'.format(x).split('e')
   b = int(b)
   return r'${}\times 10^{{{}}}$'.format(a, b)
