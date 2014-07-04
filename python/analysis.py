#! /cwi/bin/python

from csplot import *
import sys

# Update and use the corresponding SciPy functions
#from Scientific.Functions import Interpolation
import math

import numpy as num
import scipy as sci
from scipy.optimize.optimize import fmin_powell as fmin
from scipy.optimize.minpack import leastsq
from optparse import OptionParser, Option, OptionValueError

INFTY = 1e10

try:
    BASE_PATH = os.environ['CSTREAM_DEFAULT_OUT_PATH']
except KeyError:
    BASE_PATH = '/export/scratch1/luque/cstream/'

def max_in_axis(lz, m, sign, zmin=None, zmax=None):
    """ Calculates and returns the maximum of a var in the axis and
    its position.    """

    max_i = -1

    # No infinity symbol in some old versions of Python
    max_v = 0
    min_v = 1e100
    corr = 0

    dz = lz[1] - lz[0]
    
    if zmin is None:
        zmin = 0

    if zmax is None:
        zmax = INFTY
        
    for i in range(len(lz)):
        if lz[i] < zmin or lz[i] > zmax:
            continue
        
        v = m[i]
        if abs(v) > abs(max_v) and sign * v >= 0:
            max_i = i
            max_v = v

        if abs(v) < abs(min_v) and sign * v >= 0:
            min_i = i
            min_v = v

    if max_i >= 1:
        corr, max_v = quadratic_max(m[max_i - 1: max_i + 2])

    return lz[max_i] + corr * dz, max_v, min_v


def max_off_axis(lr, lz, m, sign):
    Z, R = num.meshgrid(lz, lr)
    indx = (sign * m).argmax(axis=None)
    return R.flat[indx], Z.flat[indx]

    
def find_front(var, grid, path=BASE_PATH, z0=None, ret_z0=False,
               fmin=0, fmax=INFTY, rmin=None, rmax=INFTY,
               pre_margin=400, post_margin=20, zmin=None, zmax=None,
               shift_origin=True, sign=None):


    if sign is None:
        dict_signs = {'charge': -1, 'ez': 0, 'eabs': 0}
        sign = dict_signs.get(var, 0)
        
    data = []
    lr, lz, m = load_var(var, grid, path=path)
    
    if rmin is None:
        rmin = 0.5 * (lr[1] + lr[2])
        
    if z0 is None:
        z0, max_var, _ = max_in_axis(lz, m[0, :], sign,
                                     zmin=zmin, zmax=zmax)

    fmin, fmax = abs(max_var) * fmin, abs(max_var) * fmax
    print fmin, fmax

    
    if shift_origin:
        z0_shift = z0
    else:
        z0_shift = 0

    for j in xrange(len(lz)):
        # This is a small WTF, really
        if lz[j] < z0 - pre_margin or lz[j] > z0 + post_margin:
            continue
        
        r, f = find_max_row(lr, sign * m[:, j])
        
        if abs(f) >= fmin and abs(f) < fmax and r > rmin and r < rmax:
            data.append((lz[j] - z0_shift, r, f))
        #else:
        #    print ("f out of (%f, %f) or r out of (%f, %f)" %
        #           (fmin, fmax, rmin, rmax))
            
    if ret_z0:
        return data, z0
    else:
        return data

class ThresholdNotFound(ValueError):
    pass
    
def find_envelope(grid, path=BASE_PATH, eps=0.0005, tail=350.0):
    var = 'charge'
    lr, lz, m = [num.array(a) for a in load_var(var, grid, path=path)]
    
    m = -m
    front = num.zeros((len(lz), 2), 'd')
    found = num.array([False for i in xrange(len(lz))])

    dz = lz[1] - lz[0]
    ntail = int(tail / dz)
    
    lr_inv = lr[-1::-1]
                      
    for i in xrange(len(lz)):
        try:
            a, rmax  = find_threshold(eps, m[-1::-1, i], lr_inv)
            front[i, 0:2] = lz[i], rmax
            found[i] = True
        except ThresholdNotFound:
            pass

    front = front.compress(found, axis=0)
    zlen, _ = front.shape
    if zlen > ntail:
        zlen = ntail
        
    return front[-zlen:, :]

def find_threshold(threshold, a, r):
    for ia, ir in zip(a, r):
        if ia >= threshold:
            return ia, ir

    raise ThresholdNotFound

def simple_find_front(var, grid, path=BASE_PATH, tail=350.0, eps=1e-4):
    var = 'charge'
    lr, lz, m = [num.array(a) for a in load_var(var, grid, path=path)]
    m = num.abs(m)

    maxes = m.argmax(0).squeeze()
    front = num.array([(lz[i], lr[maxes[i]])
                       for i in xrange(len(lz))
                       if maxes[i] > 0 and m[maxes[1], i] > eps])

    zlen, _ = front.shape
    if zlen > ntail:
        zlen = ntail
        
    return front[-zlen:, :]


def fit_radius(var, grid, sign, path=BASE_PATH, zmax=None):
    front, z0 = find_front(var, grid, path=path, pre_margin=400.0, sign=sign,
                           post_margin=400.0, ret_z0=True,
                           fmax=1.0, fmin=0.5, zmax=zmax)

    if not front:
        return num.NAN
    
    front = num.array(front)
    
    def errors(R):
        return num.sqrt(front[:, 1]**2 + (front[:, 0] - R)**2) - abs(R)

    try:
        x, _ = leastsq(errors, 20)
    except TypeError:
        raise
        #ipshell()
        
    return x

    
def find_max_row(lr, row):
    max_f = 0
    min_f = 1e100
    max_i = 0

    dr = lr[1] - lr[0]

    for i in xrange(len(lr)):
        f = row[i]
        if abs(f) > abs(max_f):
            max_i = i
            max_f = f

    if max_i >= 1 and max_i < len(lr) - 1:
        corr, max_f = quadratic_max(row[max_i - 1: max_i + 2])
    else:
        corr = 0
        
    return lr[max_i] + corr * dr, max_f

    

def inward_norm(nx, ny, y):
    """ gets a normal vector from the array we created making sure that
    it points always inward"""

    if num.sign(y) == num.sign(ny):
        nx, ny = -nx, -ny

    return nx, ny

def normal_to_curve(m):
    """ given an array as [(x0, y1), ...(xN, yN)] returns
    [(nx0, ny0), ... (nxN, nyN)] of components of a vector normal to the
    curve """
    print m.shape
    data = zeros(m.shape, typecode='fd')
    n = m.shape[0]
    
    for i in xrange(1, n - 1):                
        x, y = m[i, 0], m[i, 1]
        
        deltas = [(m[j, 0] - m[i, 0], m[j, 1] - m[i, 1])
                  for j in (i - 1, i + 1)]
        nx, ny = normalize(normal(*deltas))
        nx, ny = inward_norm(nx, ny, y)
        
        data[i, 0] = nx
        data[i, 1] = ny
        
    nx, ny = normalize((-(m[1, 0] - m[0, 0]),
                        m[1, 1] - m[0, 1]))

    nx, ny = inward_norm(nx, ny, m[0, 1])

    data[0, 0] = nx
    data[0, 1] = ny

    nx, ny = normalize((m[n - 2, 1] - m[n - 1, 1],
                        -(m[n - 2, 0] - m[n - 1, 0])))

    nx, ny = inward_norm(nx, ny, m[n - 1, 1])
    data[n - 1, 0] = nx
    data[n - 1, 1] = ny

    print data
    return data

def normal((x0, y0), (x1, y1)):
    print x0, y0, x1, y1
    
    if y0 == 0 and y1 == 0:
        return (0.0, 1.0)
    
    return ((-(x1**2 * y0) + (x0**2 + y0 * (y0 - y1)) * y1)
            / (x1 * y0 - x0 * y1),
            (x0**2 * x1 + x1 * y0**2 - x0 * (x1**2 + y1**2))
            / (-(x1 * y0) + x0 * y1))


def normalize((x, y)):
    a = sqrt(x**2 + y**2)
    return (x / a, y / a)


def quadratic_max(y):
    """ Finds a parabola defined by the points (dx i, y[i]) {i = 1, 2, 3}
    and returns its maximum. """
    
    try:
        f1 = (y[2] - 2 * y[1] + y[0])
        f2 = (y[2] - y[0])
    except IndexError:
        return num.NAN, num.NAN
    
    if f1 == 0:
        # The points are aligned: no parabola fit
        return 0, max(*y)
    
    return (- f2 / (2 * f1)), y[1] - f2**2 / f1 / 8.0
    

def var_interpolate(var, grid, r_staggered=False, z_staggered=False, **kwargs):
    lr, lz, var = load_var(var, grid, **kwargs)
    
    lr = num.array(lr)
    lz = num.array(lz)

    if r_staggered:
        lr = lr + 0.5 * (lr[1] - lr[0])

    if z_staggered:
        lz = lz + 0.5 * (lz[1] - lz[0])

    var = num.array(var)

    lr, var = mirror_var(lr, var)
    
    f = Interpolation.InterpolatingFunction([lr, lz], var, default=1e20)
    return f

def find_front_iter(var, grid, n, **kwargs):
    """ After finding a front, improves it by n iterations. """
    print kwargs
    
    m = find_front(var, grid, **kwargs)
    m = num.array(m)

    interp = var_interpolate(var, grid, path=kwargs.get('path', None))

    for i in xrange(n):
        norm = normal_to_curve(m)    
        m = iter_front(m, norm, interp)

    return m


def iter_front(m, norm, interp, max_step=80, ds=0.05):
    mnext = zeros(m.shape, typecode='fd')
                  
    for p in xrange(m.shape[0]):
        z0, r0 = m[p, 0], m[p, 1]
        nz, nr = norm[p, 0], norm[p, 1]

        zs = [z0 + ds * i * nz for i in xrange(-max_step, max_step)]
        rs = [r0 + ds * i * nr for i in xrange(-max_step, max_step)]

        interp_zr = [interp(r, z) for r, z in zip(rs, zs)]

        max_i, max_e = -1, 0

        for i, e in enumerate(interp_zr):
            print i, rs[i], zs[i], e
            if e > max_e or max_i < 0:
                max_i = i
                max_e = e
                
        mnext[p, 0] = zs[max_i]
        mnext[p, 1] = rs[max_i]

    return mnext

def fit_saffman(id, grid, var='eabs', width=64.0, fix_vertex=False, rmax=300):
    path = os.path.join(BASE_PATH, id)

    front = find_front(var, grid,
                       path=path, fmin=0.0001, rmin=1.1, rmax=rmax,
                       zmin=0, zmax=4000, shift_origin=False,
                       pre_margin=100, post_margin=5)
    front = num.array(front)

    #front = find_envelope(grid, path=path)
    #print front

    vert_x, vert_y, _ = front[num.abs(front).argmax(0)[0]]
    #vert_x, vert_y = front[num.abs(front).argmax(0)[0]]
    
    print "Vertex found at (%f, %f)" % (vert_x, vert_y)

    L = 2 * width
    
    def chi2((lam, C)):
        #A = 2 * width / math.pi
        #y = A * num.arccos(num.exp(((front[:, 0] - C) / A)).clip(-10000.0, 1.0))
        x = front[:, 0]
        y = (lam * L / (2 * math.pi)
             * num.arccos((2 * num.exp(2 * math.pi * x
                                      / (L * (1 - lam)) - C) - 1.0)
                          .clip(-1.0, 1.0)))
                        
        return ((y - front[:, 1]) ** 2).sum()

    if not fix_vertex:
        lam, C =  fmin(chi2, (0.5, vert_x), maxfun=100)
    else:
        def chi2_partial(lam):
            return chi2((lam, vert_x))

        lam, C =  fmin(chi2_partial, (0.5,), maxfun=100), vert_x
        

    print "\tlambda = %f\n\tx0 = %f" % (lam, C)
    return lam, C


def find_emax(rid, grid, zinterval=None, path=BASE_PATH):
    var = rid + ':eabs'
    lr, lz, m = [num.array(a) for a in load_var(var, grid, path=path)]

    z0 = lz[0]
    dz = lz[1] - z0

    if zinterval is not None:
        zmin, zmax = zinterval
        izmin, izmax = [int((z - z0) / dz) for z in zmin, zmax]

        izmin = max(0, izmin)
        izmax = min(m.shape[1] - 1, izmax)

        m = m[:, izmin:izmax]
        
    return max(m.flat)

def find_diameters(var, grid, neg_z, pos_z, path=BASE_PATH):
    lr, lz, m = load_var(var, grid, path=path)
            
    m = m * m
    diameters = num.dot(lr, m) / num.dot(ones(lr.shape), m)
    table = num.concatenate((lz.reshape(len(lz), 1),
                             diameters.reshape(len(lz), 1)), axis=1)

    table = table.compress(num.logical_and(num.greater(lz, pos_z),
                                           num.less(lz, neg_z)), axis=0)
    #ipshell()
    dfile = os.path.join(path, 'diameters.%s.tsv' % grid)
    try:
        neg_diam = num.max(table[-40:, 1])
    except ValueError:
        neg_diam = num.NAN

    try:
        pos_diam = num.max(table[:40, 1])
    except ValueError:
        pos_diam = num.NAN

    num.savetxt(os.path.join(path, 'diameters.%s.tsv' % grid), table)

    return neg_diam, pos_diam

def posneg_fronts(grid, path=BASE_PATH, var='charge', zmax=None):
    lz, m = load_axis(var, grid, path=path)
    neg_z, neg_max, neg_min = max_in_axis(lz, m, -1, zmax=zmax)
    pos_z, pos_max, pos_min = max_in_axis(lz, m, 1, zmax=zmax)
    return [neg_z, pos_z, neg_max, pos_max]


def integrate(grid, path=BASE_PATH, var='charge', sign=1.0, threshold=None,
              zrange=None):
    lr, lz, m = load_var(var, grid, path=path)
    if zrange is not None:
        z0, z1 = zrange
        if not num.isnan(z0) and not num.isnan(z1):
            lzmask = num.less(lz, z1) * num.greater(lz, z0)
            m = m.compress(lzmask, axis=1)
            lz = lz.compress(lzmask)
        
    if threshold is not None:
        if sign > 0:
            thres = num.greater
            mmax = max(num.ravel(m))
        else:
            thres = num.less
            mmax = min(num.ravel(m))
        
        m = m * thres(m, threshold * mmax)
        
    dr = lr[1] - lr[0]
    dz = lz[1] - lz[0]
    
    dV = lr.reshape((lr.size, 1)) * dr * dz
    m = m * dV
    return sum(m.flat)
    
    
# Commands
def cmd_fields(*args, **kwargs):
    try:
        path = os.path.join(BASE_PATH, args[0])
    except IndexError:
        path = '.'

    ofile = os.path.join(path, 'field.dat')
    
    fp = open(ofile, "w")
            
    pattern = kwargs.get('pattern', 'C???')
    subs = kwargs.get('subs', '')
    last = kwargs.get('last', False)

    get_grids = all_grids
    
    if last:
        get_grids = last_grid

    
    for grid in get_grids(pattern=pattern, path=path):
        print path, grid + subs
        lz, m = load_axis('ez', grid + subs, path=path)
        z, field_max, field_min = max_in_axis(lz, m, -1)
        if grid != 'C000':
            fp.write("%d\t%f\t%.18e\t%.18e\n" % (int(grid[1:]), z,
                                                 field_max, field_min))
        
            print int(grid[1:]), z, field_max, field_min
            fp.flush()

    fp.close()

def cmd_fields_last(*args, **kwargs):
    _d = kwargs.copy()
    _d['last'] = True
    
    cmd_fields(*args, **_d)

def cmd_posneg(*args, **kwargs):
    try:
        path = os.path.join(BASE_PATH, args[0])
    except IndexError:
        path = '.'

    ofile = os.path.join(path, 'posneg.dat')
    
    fp = open(ofile, "w")
    pattern = kwargs.get('pattern', 'C???')
    for grid in all_grids(pattern=pattern, path=path):
        print path, grid
        try:
            neg_z, pos_z, neg_max, pos_max = posneg_fronts(grid, path=path)
        
            if grid != 'C000':
                fp.write("%d\t%f\t%f\t%f\t%f\n"
                         % (int(grid[1:]), neg_z, pos_z, neg_max, pos_max))
        
                print int(grid[1:]), neg_z, pos_z, neg_max, pos_max
                fp.flush()
        except IOError:
            print >> sys.stderr, "IOError exception in %s" % grid

    fp.close()


def cmd_complete(*args, **kwargs):
    try:
        path = os.path.join(BASE_PATH, args[0])
    except IndexError:
        path = '.'

    ofile = os.path.join(path, kwargs.get('ofile', 'analysis.dat'))
    
    fp = open(ofile, "w")
    pattern = kwargs.get('pattern', 'C???')
    sub = kwargs.get('sub', '')
    zmax = kwargs.get('zmax', None)
    dt = kwargs.get('dt', 1.0)
    
    try:
        set_def_rescaling(kwargs['rescaling'])
    except KeyError:
        pass
    
                          
    
    fmt_str = None

    varnames = ['t',
                'Zc1', 'Zc2', 'Qmin', 'Qmax',
                'Ze1', 'Ze2', 'Emin', 'Emax',
                'V1', 'V2', 'R1', 'R2',
                'Q', 'Q1', 'Q2']
                
    fp.write("\t".join(varnames) + "\n")
    
    old_c_fronts = None
    pos_z0, neg_z0 = None, None
    
    for grid in all_grids(pattern=pattern, path=path):
        sgrid = grid + sub
        clear_load_var_memoize()
        print path, sgrid
        try:
            it = int(grid[1:]) * dt
            itms = [it]
            c_fronts = posneg_fronts(sgrid, path=path,
                                     var='charge', zmax=zmax)

            if pos_z0 is None or neg_z0 is None:
                neg_z0, pos_z0 = [c_fronts[i] for i in 0, 1]

            c_fronts[0] -= neg_z0
            c_fronts[1] -= pos_z0
                
            if old_c_fronts:
                vminus, vplus = [(c_fronts[i] - old_c_fronts[i]) / (it - old_it)
                                 for i in 0, 1]
            else:
                vminus, vplus = num.NAN, num.NAN

            old_it = it
            old_c_fronts = c_fronts
            
            neg_r, pos_r = [fit_radius('charge', sgrid, s,
                                             path=path, zmax=zmax)
                                  for s in -1, 1]
            
            #neg_diam, pos_diam = find_diameters('electrons', grid,
            #                                    c_fronts[0], c_fronts[1],
            #                                    path=path)

            itms.extend(c_fronts)

            e_fronts = posneg_fronts(sgrid, path=path, var='ez', zmax=zmax)
            e_fronts[0] -= neg_z0
            e_fronts[1] -= pos_z0

            itms.extend(e_fronts)
            itms.extend((vminus, vplus, neg_r, pos_r))
            itms.append(integrate(sgrid, path=path, var='charge'))

            # This was the old way of calculating the integrated charge:
            # neg_zrange=(neg_z0 - 2*neg_r, neg_z0 + neg_r)
            # pos_zrange=(pos_z0 - 2*pos_r, pos_z0 + pos_r)
            # 
            # itms.append(integrate(sgrid, path=path, var='charge',
            #                      sign=-1.0, threshold=0.5, zrange=neg_zrange))
            # itms.append(integrate(sgrid, path=path, var='charge',
            #                       sign=1.0, threshold=0.5, zrange=pos_zrange))

            neg_zrange=(neg_z0 - neg_r, neg_z0 + 0.5 * neg_r)
            pos_zrange=(pos_z0 - pos_r, pos_z0 + 0.5 * pos_r)
             
            itms.append(integrate(sgrid, path=path, var='charge',
                                  zrange=neg_zrange))
            itms.append(integrate(sgrid, path=path, var='charge',
                                  zrange=pos_zrange))
            
            if fmt_str == None:
                fmt_str = "\t".join("%g" for x in itms) + "\n"
                
            if int(grid[1:4]) >= 2:
                fp.write(fmt_str % tuple(itms))
                fp.flush()
        
                print "\t".join(str(x) for x in itms)
        except IOError:
            print >> sys.stderr, "IOError exception in %s" % grid

    fp.close()


    
def cmd_diameters(runid, var, **kwargs):
    try:
        path = os.path.join(BASE_PATH, runid)
    except IndexError:
        path = '.'

    ofile = os.path.join(path, 'diameters.dat')
    fp = open(ofile, "w")

    pattern = kwargs.get('pattern', 'C???')
    for grid in all_grids(pattern=pattern, path=path):
        try:

            neg_z, pos_z, neg_max, pos_max = posneg_fronts(grid, path=path)
            neg_diam, pos_diam = find_diameters(var, grid, neg_z, pos_z, path=path)
            print grid

            if grid != 'C000':
                fp.write("%d\t%f\t%f\t%f\t%f\t%f\t%f\n"
                         % (int(grid[1:]), neg_z, pos_z,
                            neg_max, pos_max, neg_diam, pos_diam))
        
                print int(grid[1:]), neg_z, pos_z, neg_max, pos_max, neg_diam,\
                      pos_diam
                fp.flush()
        except IOError:
            print >> sys.stderr, "IOError exception in %s" % grid

    fp.close()


def cmd_voltage(*args, **kwargs):
    try:
        path = os.path.join(BASE_PATH, args[0])
    except IndexError:
        path = '.'

    ofile = os.path.join(path, 'voltage.dat')
    
    fp = open(ofile, "w")
    pattern = kwargs.get('pattern', 'C???')
    for grid in all_grids(pattern=pattern, path=path):
        try:
            lz, m = load_axis('ez', grid, path=path)
            voltage = - sum(m) * (lz[1] - lz[0])
            s = "%d\t%f\n" % (int(grid[1:]), voltage)
            fp.write(s)
            fp.flush()
            print s
            
        except IOError:
            print >> sys.stderr, "IOError exception in %s" % grid


def cmd_front(grid, ofile):
    fp = open(ofile, "w")

    id, grid  = grid.split(':')
    path = os.path.join(BASE_PATH, id)

    print grid, path
    
    front = find_front('charge', grid, path=path, fmin=0.0001, rmin=2.1, rmax=300,
                       zmin=0, zmax=2000, shift_origin=False,
                       pre_margin=200, post_margin=20)
    for z, r, f in front:
        fp.write("%f\t%f\t%f\n" % (z, r, f))
    
    fp.close()

def cmd_fit_saffman(grid):
    id, grid  = grid.split(':')
    return fit_saffman(id, grid)
    
def cmd_envelope(grid, ofile):
    id, grid  = grid.split(':')
    path = os.path.join(BASE_PATH, id)

    front = find_envelope(grid, path=path)

    fp = open(ofile, "w")
    for r, z in front:
        fp.write("%f\t%f\n" % (r, z))
    fp.close()
    
    return front

def cmd_simple_front(var, grid):
    id, grid  = grid.split(':')
    path = os.path.join(BASE_PATH, id)

    front = simple_find_front(var, grid, path=path)
    for r, z in front:
        print r, z
        
    return front

def cmd_iterfront(grid, n, ofile):
    id, grid  = grid.split(':')
    path = os.path.join(BASE_PATH, id)
    
    fp = open(ofile, "w")
    n = int(n)
    
    front = find_front_iter('eabs', grid, n, path=path,
                            fmin=0.001, rmin=1, rmax=300,
                            zmin=300, zmax=600, shift_origin=False,
                            pre_margin=250, post_margin=20)
    for z, r, f in front:
        fp.write("%f\t%f\t%f\n" % (z, r, f))
    
    fp.close()


def cmd_max_off_axis(runid, var, sgn, **kwargs):
    try:
        path = os.path.join(BASE_PATH, runid)
    except IndexError:
        path = '.'

    ofile = os.path.join(path, 'max-off-%s.dat' % var)
    fp = open(ofile, "w")

    pattern = kwargs.get('pattern', 'C???')
    sgn = int(sgn)
    for grid in all_grids(pattern=pattern, path=path):
        try:
            lr, lz, m = load_var(var, grid, path=path)
        except IOError:
            print >> sys.stderr, "IOError exception in %s" % grid
            continue

        r, z = max_off_axis(lr, lz, m, sgn)
        if grid != 'C000':
            fp.write("%d\t%f\t%f\n" % (int(grid[1:]), r, z))
        
            print int(grid[1:]), r, z
            fp.flush()

    fp.close()

def cmd_shell():
    from IPython.Shell import IPShellEmbed
    ipshell = IPShellEmbed([]) 
    ipshell()

def check_eval_float(option, opt, value):
    eval_value = eval(value)
    try:
        return float(eval_value)
    except ValueError:
        raise OptionValueError(
            "option %s: invalid float value: '%s' -> %r"
            % (opt, value, eval_value))



class ExtOption(Option):
    TYPES = Option.TYPES + ("eval_float",)
    TYPE_CHECKER = Option.TYPE_CHECKER.copy()
    TYPE_CHECKER["eval_float"] = check_eval_float

    
    
def main():
    parser = OptionParser("", option_class=ExtOption)

    parser.add_option('--dt', '', action='store',
                      dest='dt',
                      type='eval_float',
                      default=1.0,
                      help="dt")

    parser.add_option('--with-units', '', action='store',
                      dest='units',
                      default=None,
                      help="Normalize with units")


    parser.add_option('--pattern', '', action='store',
                      dest='pattern', type='str',
                      default='C???',
                      help='Pattern of the analyzed grids')
    

    parser.add_option('--sub', '', action='store',
                      dest='sub', type='str',
                      default="",
                      help='Sub-grid postfix')

    parser.add_option('--ofile', '', action='store',
                      dest='sub', type='str',
                      help='Output file')


    parser.add_option('--zmax', '', action='store',
                      dest='zmax', type='eval_float',
                      default=None,
                      help='Sub-grid postfix')

    parser.add_option('--rescaling', '', action='store',
                      dest='rescaling', type='str',
                      default=None,
                      help='Rescaling to use')


    main = sys.modules[__name__]

    opts, args = parser.parse_args()
    
    try:
        cmd = sys.argv[1].replace('-', '_')
        process_method = getattr(main, 'cmd_%s' % cmd)
    except AttributeError:
        print "Dont know how to process command `%s`" % args[0]
        sys.exit(-1)
        
    process_method(*args[1:], **opts.__dict__)



if __name__ == "__main__":
    main()
