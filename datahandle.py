import numpy as np
import copy

#####################################################################
##########              INITIALIZATION PART                ##########
#####################################################################
class simjob:
    def __init__(self, jobid, dim=None, mpi=False, pui=False, endian=False, old_format=False):
        if dim==None:
            print('simulation dimension:')
            dim = int(input())
        self.mpi = mpi
        self.pui = pui
        self.form = old_format
        self.endian = endian
        self.id = jobid
        self.dim = dim

        info = self.readinfo("data.inp")

        self.npt = int(info[0]) if dim==1 else int(info[3])
        self.nx = int(info[1]) if dim==1 else int(info[4])
        self.ns = int(info[2]) if dim==1 else int(info[6])
        self.dx = float(info[3]) if dim==1 else float(info[0])
        self.dt = float(info[4]) if dim==1 else float(info[2])
        self.tu = float(info[5]) if dim==1 else float(info[7])
        self.theta = float(info[6]) if dim==1 else float(info[8])
        self.tp = float(info[7]) if dim==1 else float(info[9])
        self.nt = int(info[8]) if dim==1 else int(info[10])
        if dim==2:
            self.dy = float(info[1])
            self.ny = int(info[5])
        if mpi:
            self.mp = int(info[9]) if dim==1 else int(info[11])
        if pui:
            info = self.readinfo("datapui.inp")
            self.npi = int(info[0])
            self.mass = float(info[1])

    ##########################################****###########################
    ##########             READING DATA FROM FILE PART             ##########
    ##########################################****###########################
    def fld(self, time=None, label=None):
        dat = copy.deepcopy(self)
        dim = dat.dim
        if time==None:
            time = self.get_time(self.ns)
        dat.time = time * self.tu
        dat.x = np.linspace(0., self.dx*(self.nx-1), self.nx)
        if dim==2:
            dat.y = np.linspace(0., self.dy*(self.ny-1), self.ny)

        size = self.nx if dim==1 else self.nx * self.ny
        if label==None:
            dat.by = self.readData(self.fileinfo(time,'by'), size)
            dat.bz = self.readData(self.fileinfo(time,'bz'), size)
            dat.ex = self.readData(self.fileinfo(time,'ex'), size)
            dat.ey = self.readData(self.fileinfo(time,'ey'), size)
            dat.ez = self.readData(self.fileinfo(time,'ez'), size)
            dat.np = self.readData(self.fileinfo(time,'np'), size)
            dat.pp = self.readData(self.fileinfo(time,'pp'), size)
            dat.pl = self.readData(self.fileinfo(time,'pl'), size)
            if self.pui:
                dat.npui = self.readData(self.fileinfo(time,'npui'), size)
                dat.ppui = self.readData(self.fileinfo(time,'ppui'), size)
                dat.prui = self.readData(self.fileinfo(time,'prui'), size)

            if dim==1:
                dat.bx = np.cos(np.deg2rad(self.theta))
            if dim==2:
                dat.bx = self.readData(self.fileinfo(time,'bx'), size).reshape(self.ny, self.nx)
                dat.by = dat.by.reshape(self.ny, self.nx)
                dat.bz = dat.bz.reshape(self.ny, self.nx)
                dat.ex = dat.ex.reshape(self.ny, self.nx)
                dat.ey = dat.ey.reshape(self.ny, self.nx)
                dat.ez = dat.ez.reshape(self.ny, self.nx)
                dat.np = dat.np.reshape(self.ny, self.nx)
                dat.pp = dat.pp.reshape(self.ny, self.nx)
                dat.pl = dat.pl.reshape(self.ny, self.nx)
                if self.pui:
                    dat.npui = dat.npui.reshape(self.ny, self.nx)
                    dat.ppui = dat.ppui.reshape(self.ny, self.nx)
                    dat.prui = dat.prui.reshape(self.ny, self.nx)

            dat.bf = np.sqrt(dat.bx**2 + dat.by**2 + dat.bz**2)
            dat.an = dat.pp / dat.pl
            dat.th = np.rad2deg(np.arccos(dat.bx/dat.bf))

        else:
            if label=='bf':
                f1 = self.fileinfo(time, 'by')
                f2 = self.fileinfo(time, 'bz')
                by = self.readData(f1, size)
                bz = self.readData(f2, size)
                if dim==1:
                    bx = np.cos(np.deg2rad(self.theta))
                else:
                    f = self.fileinfo(time,'bx')
                    bx = self.readData(f, size).reshape(self.ny, self.nx)
                    by = by.reshape(self.ny, self.nx)
                    bz = bz.reshape(self.ny, self.nx)
                dat.dat = np.sqrt(bx**2+by**2+bz**2)
                dat.label = self.get_label(label)
            else:
                dat.dat = self.readData(self.fileinfo(time, label), self.nx)
                dat.label = self.get_label(label)

        return(dat)

    def ptl(self, time=None, proc=None, label=None, pui=False):
        dat = copy.deepcopy(self)
        if time==None:
            time = self.get_time(self.nt)
        dat.time = time * self.tp

        ntmp = self.npi if pui else self.npt
        if self.mpi:
            if proc==None:
                print('Process Number: 0 - %d' % (self.mp-1) )
                proc = int(input())
            ism, iem = self.mpiset(1, ntmp, self.mp, proc)
            dat.ntmp = iem - ism + 1
        else:
            dat.ntmp = ntmp

        if label==None:
            vxlabel='vxui' if pui else 'vx'
            vylabel='vyui' if pui else 'vy'
            vzlabel='vzui' if pui else 'vz'
            xplabel='xpui' if pui else 'xp'
                
            dat.vx = self.readData(self.fileinfo(time, vxlabel, proc), dat.ntmp)
            dat.vy = self.readData(self.fileinfo(time, vylabel, proc), dat.ntmp)
            dat.vz = self.readData(self.fileinfo(time, vzlabel, proc), dat.ntmp)
            dat.xp = self.readData(self.fileinfo(time, xplabel, proc), dat.ntmp)
            dat.ke = 0.5*(dat.vx**2+dat.vy**2+dat.vz**2)
            if dat.dim==2:
                yplabel = 'ypui' if pui else 'yp'
                dat.yp = self.readData(self.fileinfo(time, yplabel, proc), dat.ntmp)
            if pui:
                dat.ke = self.mass * dat.ke
        else:
            dat.dat = self.readData(self.fileinfo(time, label, proc), dat.ntmp)

        return(dat)

    ##########################################****###########################
    ##########               FILE MANIPULATION PART                ##########
    ##############################################****#######################
    def readinfo(self, infile):
        header = "../DATA/FX/" if self.endian else "../DATA/"
        jobid = self.id
        datloc = header + jobid + "/" + infile
        fd = open(datloc, "r")
        info = []
        for line in fd:
            data = line[:-1].split()
            info += [data[1]]
        fd.close()
        return info

    def fileinfo(self, time, index, proc=None):
        jobid = self.id
        header = "../DATA/FX/" if self.endian else "../DATA/"
        itime = '{0:04d}'.format(time)
        datloc = header + jobid + "/" + index + itime + "_" + '{0:03d}'.format(proc) if proc != None else \
            header + jobid + "/" + index + itime + '.dat'
        return datloc

    def readData(self, file, nx):
        f = open(file, "rb")
        dty = self.set_datatype(nx)
        dtmp = np.fromfile(f, dtype=dty, count=-1)
        data = dtmp[0]['tmp']
        f.close()
        return data

    def set_datatype(self, size):
        if self.endian:
            #        head = ('head','>i')
            #        tail = ('tail','>i')
            head = ('head','>2i')
            tail = ('tail','>2i')
            data = ('tmp','>'+str(size)+'f')
            #        dty = np.dtype([head, ('tmp','>'+str(size)+'f'), tail])
        else:
            head = ('head','<i')
            tail = ('tail','<i')
            data = ('tmp','<'+str(size)+'f')
            #        dty = np.dtype([head, ('tmp','<'+str(size)+'f'), tail])

        if self.form:
            dty = np.dtype([head, data, tail])
        else:
            dty = np.dtype([data])

        return dty

    def mpiset(self, is0, ie0, iproc, irank):
        iwork1, iwork2 = divmod(ie0-is0+1, iproc)
        ism = irank * iwork1 + is0 + min([irank, iwork2])
        iem = ism + iwork1 - 1
        if iwork2 > irank:
            iem = iem + 1
        return (ism, iem)

    def get_time(self, itmax):
        print('Time step: 0 - %d' % (itmax-1))
        time = int(input())
        return(time)

    #####################################################################
    ##########           ADVANCED DATA MAKING PART             ##########
    #####################################################################
    def ptl_minmax(self, dp):
        pminmax = np.array([ [min(dp.ke), max(dp.ke)], \
                             [min(dp.vx), max(dp.vx)], \
                             [min(dp.vy), max(dp.vy)], \
                             [min(dp.vz), max(dp.vz)] ])
        return(pminmax)

    # histogram data for particles in two dimensions
    def make_ptl_hist(self, time, mx, my, xbins, ybins, xmin, xmax, ymin, ymax, \
                      sampling_ymin=None, sampling_ymax=None, sampling_xmin=None, sampling_xmax=None, \
                      pui=False):
        htmp = {}
        impi = 0
        for i in range(self.mp):
            dp = self.ptl(time, proc=i, pui=pui)
            if self.dim==2:
                if sampling_xmin==None:
                    sampling_xmin = np.min(dp.xp)
                if sampling_xmax==None:
                    sampling_xmax = np.max(dp.xp)
                if sampling_ymin==None:
                    sampling_ymin = np.min(dp.yp)
                if sampling_ymax==None:
                    sampling_ymax = np.max(dp.yp)
                index = np.where( (dp.yp >= sampling_ymin) & (dp.yp <= sampling_ymax) & \
                                  (dp.xp >= sampling_xmin) & (dp.xp <= sampling_xmax) )
            elif self.dim==1:
                if sampling_xmin==None:
                    sampling_xmin = np.min(dp.xp)
                if sampling_xmax==None:
                    sampling_xmax = np.max(dp.xp)
                index = np.where( (dp.xp >= sampling_xmin) & (dp.xp <= sampling_xmax) )

            if i%8==0:
                print('Proc:' + str(i) + ' finished')
            if dp.xp[index].size == 0:
                continue

            x, xlabel = self.ptlquant(dp, mx, index=index)
            y, ylabel = self.ptlquant(dp, my, index=index)

            htmp[impi] = np.histogram2d(x, y, bins=[xbins, ybins], range=[[xmin,xmax],[ymin,ymax]])
            impi += 1
        X, Y, h = self.hist2d_mpi(htmp, impi)
        return (X, Y, h)

    # histogram data for particles in (x, E) plane
    # x -> assumed dimension along the heliopause
    # using data between ymin and ymax -> indicating line-of-sight integrated data across the heliopause
    def makeIBEX(self,time,xbins,ybins,emin,emax,ymin,ymax,pui=True):
        xmin = 0.0
        xmax = self.nx*self.dx
        htmp = {}
        impi = 0
        for i in range(self.mp):
            dp = self.ptl(time=time,proc=i,pui=pui)
            index = np.where( (dp.yp >= ymin) & (dp.yp <= ymax) )
            if i%8==0:
                print('Proc:' + str(i) + ' complete')
            if dp.xp[index].size == 0:
                continue
            htmp[impi] = np.histogram2d(dp.xp[index],dp.ke[index],bins=[xbins,ybins],range=[[xmin,xmax], [emin,emax]])
            impi += 1
        X, Y, ibex = self.hist2d_mpi(htmp, impi)
        return (X, Y, ibex)

    def hist2d_mpi(self, h, mp):
        hretn = h[0][0]
        for i in range(mp-1):
            hretn += h[i+1][0]
        X, Y = np.meshgrid(h[0][1], h[0][2])
        return (X, Y, hretn.T)

    ## making the distribution of the particle velocity or energy
    # quantity: assigned by the quantity index, e.g., 'vx', 'xp', 'ke', etc
    # ndiv: dividing number for the histogram
    # rmin, rmax: minimum and maximum ranges of the histogram
    # dim: specifying the simulation dimension
    # hlog: the range to be divided by the logarithmic scale
    def distfunc(self, time, quantity, ndiv, rmin, rmax, \
                 sampling_xmin=None, sampling_xmax=None, sampling_ymin=None, sampling_ymax=None, \
                 hlog=False, pui=False):
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

        for i in range(self.mp):
            dp = self.ptl(time=time, proc=i, pui=pui)
            if self.dim==2:
                if sampling_ymin == None:
                    sampling_ymin = np.min(dp.yp)
                if sampling_ymax == None:
                    sampling_ymax = np.max(dp.yp)
                if sampling_xmin == None:
                    sampling_xmin = np.min(dp.xp)
                if sampling_xmax == None:
                    sampling_xmax = np.max(dp.xp)
                index = np.where( (dp.yp > sampling_ymin) & (dp.yp < sampling_ymax) & \
                                  (dp.xp > sampling_xmin) & (dp.xp < sampling_xmax) )
            elif self.dim==1:
                if sampling_xmin == None:
                    sampling_xmin = np.min(dp.xp)
                if sampling_xmax == None:
                    sampling_xmax = np.max(dp.xp)
                index = np.where( (dp.xp > sampling_xmin) & (dp.xp < sampling_xmax)  )

            if i%8 == 0:
                print('Proc:' + str(i) + ' finished')
            if dp.xp[index].size == 0:
                continue
        
            v, label = self.ptlquant(dp, quantity, index=index)
        
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

    #### making stacked data from 1d simulation result ###
    def datastack(self, index, xskip=None, tstart=None, tend=None, tskip=None):
        if xskip == None:
            xskip = 1
        if tstart == None:
            tstart = 0
        if tend == None:
            tend = self.ns - 1
        if tskip == None:
            tskip = 1

        ds = np.array([])
        tm = np.array([])
        for i in range(tstart, tend, tskip):
            df = self.fld(time=i)

            if index == 'bf':
                dat = df.bf[0:job.nx:xskip]
            elif index == 'np':
                dat = df.np[0:job.nx:xskip]
            elif index == 'npui':
                dat = df.npui[0:job.nx:xskip]

            ds = np.hstack((ds, dat))
            tm = np.hstack((tm, df.time))

        xg = df.x[0:job.nx:xskip]

        X, Y = np.meshgrid(xg, tm)
        xsize = xg.shape[0]
        tsize = tm.shape[0]
        ds = ds.reshape(tsize,xsize)

        return (X, Y, ds)

    ##########  quantity name input -> return the corresponding data
    def fldquant(self, data, name):
        if name == 'bf':
            z = data.bf
        elif name=='bx':
            z = data.bx
        elif name=='by':
            z = data.by
        elif name=='bz':
            z = data.bz
        elif name=='ex':
            z = data.ex
        elif name=='ey':
            z = data.ey
        elif name=='ez':
            z = data.ez
        elif name=='np':
            z = data.np
        elif name=='npui':
            z = data.npui
        elif name=='ux':
            z = data.ux
        elif name=='uy':
            z = data.uy
        elif name=='uz':
            z = data.uz
        elif name=='pp':
            z = data.pp
        elif name=='pl':
            z = data.pl
        elif name=='pre':
            z = (data.pp*2.+data.pl)/3.
        elif name=='ppui':
            z = data.ppui
        elif name=='prui':
            z = data.prui
        elif name=='preui':
            z = (data.ppui*2.+data.prui)/3.

        label = self.get_label(name)
        return(z, label)

    def ptlquant(self, data, name, index=None):
        if name=='vx':
            z = data.vx
        elif name=='vy':
            z = data.vy
        elif name=='vz':
            z = data.vz
        elif name=='v':
            z = np.sqrt(data.vx**2+data.vy**2+data.vz**2)
        elif name=='ke':
            z = data.ke
        elif name=='xp':
            z = data.xp
        elif name=='yp':
            z = data.yp

        if index != None:
            z = z[index]
        label = self.get_label(name)

        return(z, label)

    def get_label(self, index):
        if index == 'bf':
            label = r'$|B|$'
        elif index=='bx':
            label = r'$B_x$'
        elif index=='by':
            label = r'$B_y$'
        elif index=='bz':
            label = r'$B_z$'
        elif index=='ex':
            label = r'$E_x$'
        elif index=='ey':
            label = r'$E_y$'
        elif index=='ez':
            label = r'$E_z$'
        elif index=='np':
            label = r'$N_p$'
        elif index=='npui':
            label = r'$N_{PUI}$'
        elif index=='ux':
            label = r'$U_x$'
        elif index=='uy':
            label = r'$U_y$'
        elif index=='uz':
            label = r'$U_z$'
        elif index=='pp':
            label = r'$P_\perp$'
        elif index=='pl':
            label = r'$P_\parallel$'
        elif index=='pre':
            label = 'P'
        elif index=='ppui':
            label = r'$P_{PUI,\perp}$'
        elif index=='prui':
            label = r'$P_{PUI, \parallel}$'
        elif index=='preui':
            label = r'$P_{PUI}$'
        elif index=='vx':
            label = r'$v_x$'
        elif index=='vy':
            label = r'$v_y$'
        elif index=='vz':
            label = r'$v_z$'
        elif index=='v':
            label = r'$v$'
        elif index=='ke':
            label = r'$KE$'
        elif index=='xp':
            label = r'$x$'
        elif index=='yp':
            label = r'$y$'

        return(label)
        

### giving parallel and perpendicular velocity as well as pitch angle from particle data
#def particle_pitch(job, time):
