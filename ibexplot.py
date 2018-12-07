import matplotlib.pyplot as plt
import myplot as myp

def ibexplot(job,cmap='Oranges'):
    fig, ax = plt.subplots(figsize=(6,8), nrows=7,ncols=1)
    for i in range(7):
        X, Y, dc  = myp.makeIBEX(job,(i+1)*4,200,200,30,70,440,700)
        ax[i].pcolormesh(X, Y, dc, cmap=cmap, rasterized=True)
        ax[i].set_title('$\Omega_p t=$'+str(4.*(i+1)*job.tp))

    for i in range(6):
        ax[i].axes.xaxis.set_ticklabels([])

    ax[6].axes.set_xlabel('x', fontsize=18)
    ax[3].axes.set_ylabel('Energy', fontsize=18)

    plt.tight_layout()
