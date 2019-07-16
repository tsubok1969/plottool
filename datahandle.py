import numpy as np

#class
class job1d:
    def __init__(self,jobid,mpi=False,pui=False,endian=False,old_format=False):
        header = header_endian(endian)
        info = readinfo(jobid, header, "data.inp")
        self.form = old_format
        self.endian = endian
        self.info = info
        self.id = jobid
        self.np = int(info[0])
        self.nx = int(info[1])
        self.ns = int(info[2])
        self.dx = float(info[3])
        self.dt = float(info[4])
        self.tu = float(info[5])
        self.th = float(info[6])
        self.tp = float(info[7])
        self.nt = int(info[8])
        if mpi:
            self.mp = int(info[9])
        if pui:
            info = readinfo(jobid, header, "datapui.inp")
            self.npui = int(info[0])
            self.mass = float(info[1])

class job2d:
    def __init__(self, jobid, mpi=False, pui=False, endian=False, old_format=False):
        header = header_endian(endian)
        info = readinfo(jobid, header, "data.inp")
        self.form = old_format
        self.endian = endian
        self.info = info
        self.id = jobid
        self.dx = float(info[0])
        self.dy = float(info[1])
        self.dt = float(info[2])
        self.np = int(info[3])
        self.nx = int(info[4])
        self.ny = int(info[5])
        self.ns = int(info[6])
        self.tu = float(info[7])
        self.th = float(info[8])
        self.tp = float(info[9])
        self.nt = int(info[10])
        if mpi:
            self.mp = int(self.info[11])
        if pui:
            info = readinfo(jobid, header, "datapui.inp")
            self.npui = int(info[0])
            self.mass = float(info[1])

class fld1D():
    def __init__(self, job, time, pui=False):
        self.job = job
        self.time = time * job.tu
        self.x = np.linspace(0., self.job.dx*(self.job.nx-1), self.job.nx)

        endian = job.endian
        index = ("by", "bz", "ex","ey","ez","np","pl","pp")
        file = []
        for i in range(8):
            f = fileinfo(job, time, index[i])
            file.append(f)
        self.by = read1Ddata(file[0], endian, self.job.nx, self.job.form)
        self.bz = read1Ddata(file[1], endian, self.job.nx, self.job.form)
        self.ex = read1Ddata(file[2], endian, self.job.nx, self.job.form)
        self.ey = read1Ddata(file[3], endian, self.job.nx, self.job.form)
        self.ez = read1Ddata(file[4], endian, self.job.nx, self.job.form)
        self.np = read1Ddata(file[5], endian, self.job.nx, self.job.form)
        self.pl = read1Ddata(file[6], endian, self.job.nx, self.job.form)
        self.pp = read1Ddata(file[7], endian, self.job.nx, self.job.form)

        bx = np.cos(np.deg2rad(job.th))
        self.bf = np.sqrt(bx**2 + self.by**2 + self.bz**2)
        self.an = self.pp / self.pl
        self.th = np.rad2deg(np.arccos(bx/self.bf))

        if pui:
            index = ('npui', 'ppui', 'prui')
            file = []
            for i in range(3):
                f = fileinfo(job, time, index[i])
                file.append(f)
            self.npui = read1Ddata(file[0], endian, self.job.nx, self.job.form)
            self.ppui = read1Ddata(file[1], endian, self.job.nx, self.job.form)
            self.prui = read1Ddata(file[2], endian, self.job.nx, self.job.form)
            
class fld2D():
    def __init__(self,job,time,pui=False):
        self.job = job
        self.time = time * job.tu
        endian = job.endian
        index = ('bx','by','bz','ex','ey','ez','ux','uy','uz','np','pp','pl')
        file = []
        for i in range(12):
            f = fileinfo(job, time, index[i])
            file.append(f)
        self.bx = read2Ddata(file[0], endian, self.job.nx, self.job.ny, self.job.form)
        self.by = read2Ddata(file[1], endian, self.job.nx, self.job.ny, self.job.form)
        self.bz = read2Ddata(file[2], endian, self.job.nx, self.job.ny, self.job.form)
        self.ex = read2Ddata(file[3], endian, self.job.nx, self.job.ny, self.job.form)
        self.ey = read2Ddata(file[4], endian, self.job.nx, self.job.ny, self.job.form)
        self.ez = read2Ddata(file[5], endian, self.job.nx, self.job.ny, self.job.form)
        self.ux = read2Ddata(file[6], endian, self.job.nx, self.job.ny, self.job.form)
        self.uy = read2Ddata(file[7], endian, self.job.nx, self.job.ny, self.job.form)
        self.uz = read2Ddata(file[8], endian, self.job.nx, self.job.ny, self.job.form)
        self.np = read2Ddata(file[9], endian, self.job.nx, self.job.ny, self.job.form)
        self.pp = read2Ddata(file[10], endian, self.job.nx, self.job.ny, self.job.form)
        self.pl = read2Ddata(file[11], endian, self.job.nx, self.job.ny, self.job.form)

        self.bf = np.sqrt(self.bx**2+self.by**2+self.bz**2)

        self.x = np.linspace(0., self.job.dx*(self.job.nx-1), self.job.nx)
        self.y = np.linspace(0., self.job.dy*(self.job.ny-1), self.job.ny)

        if pui:
            index = ('npui','ppui','prui')
            file = []
            for i in range(3):
                f = fileinfo(job,time,index[i])
                file.append(f)
            self.npui = read2Ddata(file[0], endian, self.job.nx, self.job.ny, self.job.form)
            self.ppui = read2Ddata(file[1], endian, self.job.nx, self.job.ny, self.job.form)
            self.prui = read2Ddata(file[2], endian, self.job.nx, self.job.ny, self.job.form)

class ptl1D():
    def __init__(self, job, time, mpi=False, proc=None, pui=False):
        self.job = job
        self.time = time * job.tp
        endian = job.endian
        if pui:
            index = ('vxui','vyui','vzui','xpui')
            ntmp = job.npui
        else:
            index = ('vx','vy','vz','xp')
            ntmp = job.np
        file = []
        for i in range(4):
            f = fileinfo(job, time, index[i], mpi, proc)
            file.append(f)
        if mpi:
            ism, iem = mpiset(1, ntmp, job.mp, proc)
            self.np = iem - ism + 1
        else:
            self.np = ntmp
        self.vx = read1Ddata(file[0], endian, self.np, self.job.form)
        self.vy = read1Ddata(file[1], endian, self.np, self.job.form)
        self.vz = read1Ddata(file[2], endian, self.np, self.job.form)
        self.xp = read1Ddata(file[3], endian, self.np, self.job.form)
        self.ke = 0.5*(self.vx**2+self.vy**2+self.vz**2)
        if pui:
            self.ke = job.mass * self.ke

class ptl2D(ptl1D):
    def __init__(self, job, time, mpi=False, proc=None, pui=False):
        super().__init__(job,time,mpi,proc,pui)
        if pui:
            index = 'ypui'
        else:
            index= 'yp'
        f = fileinfo(job, time, index, mpi, proc)
        self.yp = read1Ddata(f, job.endian, self.np, self.job.form)
            
# function
def header_endian(endian):
    if endian:
        header = "../DATA/FX/"
    else:
        header = "../DATA/"
    return header
def readinfo(jobid, header, infile):
    datloc = header + jobid + "/" + infile
    fd = open(datloc, "r")
    info = []
    for line in fd:
        data = line[:-1].split()
        info += [data[1]]
    fd.close()
    return info

def fileinfo(job, time, index, mpiproc=False, proc=None):
    jobid = job.id
    header = header_endian(job.endian)
    itime = '{0:04d}'.format(time)
    if mpiproc:
        datloc = header + jobid + "/" + index + itime + "_" + '{0:03d}'.format(proc)
    else:
        datloc = header + jobid + "/" + index + itime + '.dat'
    return datloc

def read1Ddata(file, endian, nx, old_format):
    f = open(file, "rb")
    dty = set_datatype(endian, nx, old_format)
    dtmp = np.fromfile(f, dtype=dty, count=-1)
    data = dtmp[0]['tmp']
    f.close()
    return data

def read2Ddata(file, endian, nx, ny, old_format=False):
    size = nx*ny
    data = read1Ddata(file, endian, size, old_format).reshape(ny,nx)
    return data

def set_datatype(endian, size, old_format):
    if endian:
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

    if old_format:
        dty = np.dtype([head, data, tail])
    else:
        dty = np.dtype([data])

    return dty

def mpiset(is0, ie0, iproc, irank):
    iwork1, iwork2 = divmod(ie0-is0+1, iproc)
    ism = irank * iwork1 + is0 + min([irank, iwork2])
    iem = ism + iwork1 - 1
    if iwork2 > irank:
        iem = iem + 1
    return (ism, iem)
