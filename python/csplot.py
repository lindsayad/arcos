import matplotlib
from matplotlib.ticker import Formatter
matplotlib.use('Agg')

#from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']))
#rc('text', usetex=True)

from pylab import *
import os, os.path
import sys
import glob
import tarfile
import math
import scaling
import re
import numpy
from subinterpol import load_fine

rcParams['contour.negative_linestyle'] = 'solid'

output_dir = '.'
#try:
#    output_dir = os.environ['CSTREAM_OUTPUT_DIR']
#except KeyError:
#    output_dir = '/export/scratch1/luque/cstream/data/'

try:
    BASE_PATH = os.environ['CSTREAM_DEFAULT_OUT_PATH']
except KeyError:
    BASE_PATH = '/export/scratch1/luque/cstream/'

LOG_EPS = 1e-20

# Some rescalings:
nitrogen1atm_cm = {'r': 2.3e-4, 'z': 2.3e-4,
                   'electrons': 1.0 / (2.3e-4)**3}

nitrogen1atm_microm = {'r': 2.3, 'z': 2.3,
                       'electrons': 4.7e2,
                       'charge': 4.7e2,
                       'eabs': 200,
                       'ez': 200}

nitrogen1atm_microm_units = {'r': r'\mu m', 'z': r'\mu m',
                             'electrons': r'\mu m^{-3}',
                             'charge': r'e/\mu m^{3}',
                             'eabs': r'kV/cm'}

nitrogen1atm_mm = {'r': 2.3e-3, 'z': 2.3e-3,
                   'electrons': 4.7e2,
                   'charge': 4.7e2,
                   'eabs': 200}

nitrogen1atm_mm_units = {'r': r'mm', 'z': r'mm',
                         'electrons': r'\mu m^{-3}',
                         'charge': r'e/\mu m^{3}',
                         'eabs': r'kV/cm'}



nitrogen5atm_microm = scaling.rescale2(1, 5, nitrogen1atm_microm)
nitrogen5atm_microm_units = nitrogen1atm_microm_units


nitrogen15mbar_mm = {'r': 0.153, 'z': 0.153,
                     'electrons': 1.0 / (0.153)**3}

nitrogen15mbar_mm_units = {'x': r'mm', 'y': r'mm',
                           'electrons':r'mm^{-3}'}

air100mbar_mm = {'r': 0.0233, 'z': 0.0233,
                 'electrons': 1.0 / (0.0233)**3,
                 'ions': 1.0 / (0.0233)**3,
                 'charge': 1.0 / (0.0233)**3
                 }

air100mbar_mm_units = {'r': r'mm', 'z': r'mm',
                       'electrons':r'mm^{-3}',
                       'ions': r'mm^{-3}',
                       'charge': r'e/mm^{3}'
                       }


nitrogen005torr_m = {'r': 0.035, 'z': 0.035,
                     'electrons': 1.0 / (0.035)**3,
                     'charge': 1.0 / (0.035)**3,
                     }

nitrogen005torr_m_units = {'r': r'm', 'z': r'm',
                           'electrons':r'm^{-3}',
                           'charge': r'e/m^{3}',
                           }


micrombar = nitrogen1atm_microm
micrombar_units = {'pr': r'\mu m \cdot bar',
                   'px': r'\mu m \cdot bar',
                   'pz': r'\mu m \cdot bar',
                   'r': r'\mu m \cdot bar',
                   'x': r'\mu m \cdot bar',
                   'z': r'\mu m \cdot bar',
                   'electrons':r'\mu m^{-3} \cdot bar^{-3}'}

# Rescaling for sprites:
def get_sprite_f(dens_decay_len):
    def _f(lz):
        return numpy.exp(-lz / dens_decay_len)
    return _f

DENS_DECAY_LENGTH_KM = 7.2

def dens_at(h_km):
    """ Density at h_km kilometers relative to the density at ground level."""
    return numpy.exp(-h_km / DENS_DECAY_LENGTH_KM)

def rescale_sprite_at(h_km):
    dens = dens_at(h_km)
    l0 = 2.3e-9 / dens
    f = get_sprite_f(DENS_DECAY_LENGTH_KM / l0)
    n0 = dens * 1.0 / (2.3e-4)**3
    return {'r': l0, 'z': -l0, 'z0': h_km,
            'eabs_f': f, 'ez_f': f, 'er_f': f,
            'electrons': n0, 'ions': n0,
            'charge': n0            }

sprite70km = rescale_sprite_at(70)
sprite70km_units = {'r': 'km', 'z': 'km', 'electrons': 'cm^{-3}',
                    'log_electrons': 'cm^{-3}'}
sprite80km = rescale_sprite_at(80)
sprite80km_units = sprite70km_units
sprite82km = rescale_sprite_at(82)
sprite82km_units = sprite70km_units
sprite90km = rescale_sprite_at(90)
sprite90km_units = sprite70km_units


PRETTY_VARS = {
    'electrons': 'n_e',
    'ions': 'n_+',
    'eabs': 'E',
    'er': 'E',
    '_er': '-E',
    '_ez': '$-E_z$',
    'ez': 'E',
    'phi': r'\phi',
    'log_ez': r'$\rm{log}\/ E_z$',
    'log_er': r'$\rm{log}\/ E_r$',
    'log_eabs': r'$\rm{log}\/ E$',
    'log_electrons': 'n_e',
    'charge': r'Charge density',
    'los_impact': r'$I(y,z)$'
    #'charge': r'\rho - \sigma'
    }

def pretty_var(var):
    return PRETTY_VARS.get(var, var)

def flex_load(fname, path=''):
    try:
        # 1. Try with `filename'
        return load(os.path.join(path, fname))
    except IOError:
        try:
            # 2. If not, try with `filename.gz'
            return load(os.path.join(path, fname + '.gz'))
        except IOError:
            # 2. If not, extract it from an all.<grid>.tgz '
            tar = tarfile.open(os.path.join(path, tar_name(fname)), 'r:gz')
            tar.extract(os.path.basename(fname))
            r = load(fname)
            os.remove(fname)
            return r

def default_loader(var, grid, path=''):
    return flex_load(input_name(var, grid), path=path)
            
def charge_loader(var, grid, path=''):
    return (default_loader('ions', grid, path=path)
            - default_loader('electrons', grid, path=path))

def d_charge_loader(var, grid, path=''):
    return (default_loader('d_ions', grid, path=path)
            - default_loader('d_electrons', grid, path=path))

def current_loader(var, grid, path=''):
    eabs = default_loader('eabs', grid, path=path)
    electrons = default_loader('electrons', grid, path=path)

    return eabs * electrons

def log_loader(var, grid, path=''):
    loader = special_vars.get(var[4:], default_loader);
    return log10(abs(loader(var[4:], grid, path=path)) + LOG_EPS)

def minus_loader(var, grid, path=''):
    real_var = var[1:]
    loader = special_vars.get(real_var, default_loader);
    return -loader(real_var, grid, path=path)
    
def townsend(efield, efield0):
    return efield * exp(-efield0 / abs(efield))

def impact_loader(var, grid, path=''):
    eabs = default_loader('eabs', grid, path=path)
    electrons = default_loader('electrons', grid, path=path)

    return townsend(eabs, 1.0) * electrons

def los_loader(var, grid, path=''):
    # Only used here:
    from losight import losight
    real_var = var[var.index('_') + 1:]
    #pre_loader = special_vars.get(real_var, default_loader)

    #lr = flex_load(input_name('r', grid), path=path)
    #lz = flex_load(input_name('z', grid), path=path)
    #m = pre_loader(real_var, grid, path=path)
    lr, lz, m = load_fine(real_var, grid, path, uptolevel=2)
    
    return lr, lz, losight(lr, lz, m)
    
special_vars = {'ccharge': charge_loader,
                'd_charge': d_charge_loader,
                'log_electrons': log_loader,
                'log_ions': log_loader, 'log_ccharge': log_loader,
                'log_charge': log_loader,
                'log_phi': log_loader,
                'log_er': log_loader,
                '_er': minus_loader,
                '_ez': minus_loader,
                '_charge': minus_loader,
                'log_ez': log_loader,
                'log_eabs': log_loader,
                'log_error': log_loader,
                'log_d_electrons': log_loader,
                'log_d_ions': log_loader,
                'current': current_loader,
                'impact': impact_loader,
                'los_impact': los_loader,
                'log_los_impact': log_loader
                }

def input_name(var, grid):
    return "%s.%s.tsv" % (var, grid)

def output_name(var, grid, ext='png', what=""):
    return "%s%s.%s.%s" % (var, what, grid, ext)
    
RE_FNAME = re.compile(r".*\.(\w\d{3}).*\.tsv")
def tar_name(fname):
    m = RE_FNAME.match(fname)
    if not m:
        raise ValueError("File-name format %s not recognized." % fname)

    gid = m.group(1)
    return "all.%s.tgz" % gid
    
def old_tar_name(fname):
    try:
        for i in range(fname.index(".C") + 2, len(fname)):
            if not fname[i].isdigit():
                return 'all.%s.tgz' % fname[fname.index(".C") + 1:i]

    except ValueError:
        indx = fname.index(".P")
        return 'all.%s.tgz' % fname[indx + 1:indx + 5]
        
    return fname

def get_subgrids(gid, path=output_dir):
    try:
        # The '.' is a dirty hack to avoid changing tar_name, which expects
        # a filename.
        tar = tarfile.open(os.path.join(path, tar_name('.' + gid + '.tsv')))
        gridnames = tar.getnames()
        tar.close()
    except (IOError, ValueError):
        gridnames = [os.path.basename(g) for g in
                     glob.glob(os.path.join(path, "z.%s*.tsv" % gid))]
        
    grids = [g for z, g, tsv in [nam.split('.') for nam in gridnames]
             if z == 'z']
    ret = [g for g in grids if g[0:len(gid)] == gid and g != gid]

    return ret
    
def all_grids(pattern="C???", path=output_dir):
    i = glob.glob(os.path.join(path, "all.%s.tgz" % pattern))
    if not i:
        i = glob.glob(os.path.join(path, "z.%s.tsv" % pattern))

    i.sort()
    for fn in i:
        yield fn.split('.')[-2]

def last_grid(pattern="C???", path=output_dir):
    i = glob.glob(os.path.join(path, "all.%s.tgz" % pattern))
    i.sort()
    try:
        yield i[-1].split('.')[-2]
    except IndexError:
        # There are no grids
        pass
    
LOAD_VAR_MEMOIZE = {}

def clear_load_var_memoize():
    global LOAD_VAR_MEMOIZE
    del LOAD_VAR_MEMOIZE
    LOAD_VAR_MEMOIZE = {}

_default_rescaling = {}

def set_def_rescaling(rescaling):
    global _default_rescaling

    if isinstance(rescaling, str):
        mainmod = sys.modules[__name__]
        rescaling = getattr(mainmod, rescaling)

    if rescaling is None:
        rescaling = {}
    
    _default_rescaling = rescaling
    

def load_var(var, grid, path='', ntheta=None, opposite=False,
             rescaling=None, zcross=None):
    """ Loads a variable read from a grid inside a matrix.
    Returns a tuple formed by lr, lz, m, where lr and lz are 1d vectors
    representing the gridpoints in r and z. """

    global _default_rescaling
    
    try:
        id, var = var.split(':')
        path = os.path.join(BASE_PATH, id)
    except ValueError:
        pass
    
    if rescaling is None:
        rescaling = _default_rescaling
        
    k = (path, var, grid, rescaling.get(var, 1.0), ntheta, zcross)
    if k in LOAD_VAR_MEMOIZE:
        return LOAD_VAR_MEMOIZE[k]
    
    loader = special_vars.get(var, default_loader);

    # Load matrix
    m = loader(var, grid, path=path)

    try:
        ltheta = flex_load(input_name('theta', grid), path=path)
    except:
        ltheta = array([0])

    # Some loaders already provide lr, lz and reshaped m
    if isinstance(m, tuple):
        lr, lz, m = m
    else:
        lr = flex_load(input_name('r', grid), path=path)
        lz = flex_load(input_name('z', grid), path=path)

    # Reshape
    m = reshape(m, (ltheta.shape[0], lr.shape[0], lz.shape[0]))


    # Remove the theta dimension
    if ntheta is None:
        ntheta = min(ltheta.shape[0] - 1, 1)

    if zcross is None:
        mr = m[ntheta, :, :]
    else:
        iz = (zcross - lz[0]) / (lz[1] - lz[0])
        mr = m[:, :, iz]
    
        lz = numpy.outer(numpy.cos(ltheta), lr)
        lr = numpy.outer(numpy.sin(ltheta), lr)


    if opposite:
        max_theta = ltheta.shape[0] - 1
        #ntheta2 = (ntheta + max_theta / 2 - 1) % max_theta
        ntheta2 = (ntheta + max_theta / 2) % max_theta + 2
        
        mr = concatenate((flipud(m[ntheta2, :, :]), mr), axis=0)
        lr = concatenate((-lr[::-1], lr))
        

    try:
        f = rescaling['%s_f' % var]
        mr = f(lz) * mr
    except KeyError:
        if not var.startswith('log_'):
            mr = rescaling.get(var, 1.0) * mr
        else:
            mr = log10(rescaling.get(var[4:], 1.0)) + mr
            
    lr = rescaling.get('r', 1.0) * lr + rescaling.get('r0', 0.0)
    lz = rescaling.get('z', 1.0) * lz + rescaling.get('z0', 0.0)
 
    if lz[0] > lz[-1]:
        mr = numpy.fliplr(mr)
        lz = numpy.flipud(lz)

    LOAD_VAR_MEMOIZE[k] = (lr, lz, mr)
    return lr, lz, mr

    
def load_axis(var, grid, path=output_dir, col=0, **kwargs):
    """ Reads the values along the minimum r of var"""
    lr, lz, m = load_var(var, grid, path=path, **kwargs)
    return lz, m[col,:]

def save_var(var, grid, lr, lz, m):
    save(input_name('r', grid), lr)
    save(input_name('z', grid), lz)
    save(input_name(var, grid), m)
    

# For potential plots, we use a lot more contours
NCONTOURS = {'phi': 100, 'charge': 20}

def add_to_dict(thedict, key):
    def doit(f):
        thedict[key] = f
        return f
    
    return doit
    
def mirror_var(lr, m):
    m = concatenate((flipud(m), m), axis=0)
    lr = concatenate((-lr[::-1], lr))

    return lr, m

def mirrored_var(lr, m):
    m = flipud(m)
    lr = -lr[::-1]

    return lr, m

prev_contours = None

def plot_var(var, grid, mode='pcolor', aspect=True, mirror=False,
             setaxis=True, finish=True, var_range=None, path='',
             vector=None, ntheta=None, rescaling={}, units={},
             fontsize=18, rlabel='r', zlabel='z', clabel=None,
             rrange=None, zrange=None,
             saffman=None, cmap=None, color='black', linewidth=None,
             opposite=False, normranges=False, subgrids=[],
             touchaxis=True, frame=False, superimpose=None, ncontours=None,
             field=None, zcross=None, fine=False, mirrored=None,
             subframe=False):
    
    if not fine:
        lr, lz, m = load_var (var, grid, path=path, ntheta=ntheta,
                              opposite=opposite, rescaling=rescaling,
                              zcross=zcross)
    else:
        def _loader(_var, _grid, path='.'):
            return load_var (_var, _grid, path=path, ntheta=ntheta,
                              opposite=opposite, rescaling=rescaling,
                              zcross=zcross)
        lr, lz, m = load_fine(var, grid, path=path, loader=_loader)
        
        
    if mirror:
        lr, m = mirror_var(lr, m)

    if mirrored == -1:
        lr, m = mirrored_var(lr, m)
        
    if zcross is None:
        [z, r] = meshgrid(lz, lr)
    else:
        z, r = lz, lr
        
    if field is not None and var == 'phi':
        m = m - field * z

    modes = {}
    @add_to_dict(modes, 'pcolor')
    def plot_pcolor():
        if var_range is not None:
            vmin, vmax = var_range
            norm = normalize(vmin=vmin, vmax=vmax)
        else:
            vmin, vmax = None, None
            norm = normalize()

        if hasattr(cmap, 'center'):
            # A dynamic colormap
            if vmin is None or vmax is None:
                vmin = min(m.flat)
                vmax = max(m.flat)

            cmap.center = -vmin / (vmax - vmin)
            
        pcolor(r, z, m, shading='flat', vmin=vmin, vmax=vmax,
               norm=norm, cmap=cmap)
        
    @add_to_dict(modes, 'contour')
    def plot_contour():
        global prev_contours
        ncont =ncontours or NCONTOURS.get(var, 1)
        
        if not prev_contours:
            prev_contours = contour(r, z, m, ncont,
                                   colors=color, linewidths=linewidth)
        else:
            contour(r, z, m, prev_contours[0], colors='black')
            
    @add_to_dict(modes, 'neg_shape')
    def plot_neg_shape():
        themin = min(m.flat)
        contour(r, z, m, [themin / 2], colors=color, linewidth=linewidth)

    @add_to_dict(modes, 'neg_filled')
    def plot_neg_filled():
        themin = min(m.flat)
        contourf(r, z, m, [themin / 2, themin], cmap=cmap)

    @add_to_dict(modes, 'pos_shape')
    def plot_pos_shape():
        themax = max(m.flat)
        contour(r, z, m, [themax / 2], colors=color, linewidth=linewidth)

    @add_to_dict(modes, 'pos_filled')
    def plot_pos_filled():
        themax = max(m.flat)
        contourf(r, z, m, [themax / 2, themax], cmap=cmap)

    @add_to_dict(modes, 'range')
    def output_range():
        themax, themin = max(m.flat), min(m.flat)
        print "%g:%g" % (themin, themax)
    
    modes[mode]()

    if frame:
        hlines([lz[0], lz[-1]], lr[0], lr[-1], "w", linewidth=1.5)
        vlines([lr[0], lr[-1]], lz[0], lz[-1], "w", linewidth=1.5)
        
               
    if saffman is not None:
        plot_saffman(color='blue', *saffman)

    if superimpose is not None:
        try:
            fname, shift, factor = superimpose
        except ValueError:
            fname, shift, factor = superimpose, 0.0, 1.0
            
        plot_xyfile(fname, shift=shift, factor=factor, color='blue')
                    
    if aspect:
        axis('scaled')
    
    if setaxis and zcross is None:
        axis([0, 0.75, 0, 1.0])
        xlim(lr[0], lr[-1])
        ylim(lz[0], lz[-1])

    if rrange:                
        if normranges:
            factor = rescaling.get('r', 1.0)
            rrange = [factor * rl for rl in rrange]

        xlim(*rrange)
        
    if zrange:
        if normranges:
            factor = rescaling.get('z', 1.0)
            zrange = [factor * zl for zl in zrange]

        ylim(*zrange)

    #@print_ret
    def get_label(x):
        theunits = units.get(x, '')
        if theunits:
            return (r'$%s$ [$\mathdefault{%s}$]'
                    % (pretty_var(x), theunits))
        else:
            return r'%s' % pretty_var(x)


    if touchaxis:
        xlabel(get_label(rlabel), fontsize=fontsize)
        ylabel(get_label(zlabel), fontsize=fontsize)
    
        
    xmin, xmax = None, None
    
    ax = gca()
    ax.xaxis.set_major_formatter(FancyFormatter())
    ax.yaxis.set_major_formatter(FancyFormatter())



    if mode == 'pcolor' and touchaxis:

        #cax = colorbar(format="%.3g", cax=ncax)

        # This was for the Saffman-Taylor paper:
        #cax = colorbar(format="%.3g", fraction=0.20, pad=0.00, aspect=12)

        format = FancyFormatter()
        if 'log' in var:
            format = r'$\mathdefault{10^{%d}}$'
        
        #format = 'log' in var and '$10^{%d}$' or "%.3g"
        cax = colorbar(format=format)
        #for t in cax.get_yticklabels():
        #    t.set_fontsize(fontsize - 2)
        #gcf().sca(ax)
        #ipshell()
        xmin, xmax = cax.get_clim()
        
        #trans = blend_xy_sep_transform(cax.transAxes, cax.transData)
        #xpos = 2.5
        #xpos = 4.0
        #xpos = -0.5
        #cax.text(xpos, 0.5 * (xmin + xmax),
        #         get_label(var), ha='right', va='center', transform=trans,
        #         rotation='vertical', fontsize=fontsize)
        if clabel is None:
            clabel = get_label(var)
            
        cax.set_label(clabel, fontsize=fontsize)

    if vector is not None:
        if ':' in vector:
            grid, vector = vector.split(':')
            
        d = {'field': False, 'current': True}
        plot_some_vect(grid, path=path, mirror=mirror, current=d[vector],
                       rescaling=rescaling)
        

    for subgrid in subgrids:
        if mirror or mirrored is not None:
            mirrs = (-1, 1)
        else:
            mirrs = (1,)

        for mirr_ in mirrs:
            plot_var(var, subgrid, mode=mode, setaxis=False, finish=False,
                     path=path, units=units, opposite=False, cmap=cmap,
                     mirror=False, frame=subframe, var_range=[xmin, xmax],
                     rescaling=rescaling, mirrored=mirr_, rrange=rrange,
                     zrange=zrange, aspect=False, touchaxis=False,
                     subframe=subframe)


    if finish:
        show()

def plot_saffman(width, C, color='white'):
    eps = 0.05
    A = width / math.pi
    x = arange(-width / 2 + eps, width / 2 - eps, 0.01)
    y = A * log(cos(x / A))+ float(C)
    
    #lam = 0.5
    #L = 2 * width
    #y = (L * ( 1 - lam) / (2 * math.pi)
    #     * log (0.5 * (1 + cos(2 * math.pi * x / (lam * L))))) + float(C)
    
    plot(x, y, linewidth=1.5, color=color)

def plot_xyfile(fname, shift=0, factor=1.0, color='white'):
    data = load(fname)
        
    plot(factor * data[:, 0], factor * data[:, 1] + shift,
         linewidth=1.5, color=color)


def plot_some_vect(grid, path='', dr=1, dz=1,
                   rescaling={}, mirror=False, current=False):
    """ Plots either the electric field (if current == False) or the
    current vector field (if current == True)"""
    
    lr, lz, er = load_var('er', grid, path=path, rescaling=rescaling) 
    void, void, ez = load_var('ez', grid, path=path, rescaling=rescaling) 

    stepr, stepz = int(dr / (lr[1] - lr[0])), int(dr / (lz[1] - lz[0])), 

    lr = lr[stepr / 2::stepr]
    lz = lz[::stepz]

    er = er[stepr / 2::stepr, 0::stepz]
    ez = ez[stepr / 2::stepr, 0::stepz]

    if mirror:
        er = concatenate((-flipud(er), er), axis=0)
        ez = concatenate((flipud(ez), ez), axis=0)
        
        lr = concatenate((-lr[::-1], lr))

    if current:
        void, void, electrons = load_var('electrons', grid, path)
        electrons = electrons[stepr / 2::stepr, 0::stepz]
        if mirror:
            electrons = concatenate((flipud(electrons), electrons), axis=0)
    else:
        electrons = 1.0

    mr, mz = electrons * er, electrons * ez

    eunit = 1.5e-3 * dr

    plot_vect(lr, lz, mr, mz, eunit=eunit)

def plot_vect(lr, lz, mr, mz, eunit=1.0, threshold=0.0125):
    EPS = 1e-15

    mr = mr + EPS
    mz = mz + EPS
    
    N = sqrt(mr**2 + mz**2)
    if threshold is not None:
        mr[N > threshold]=0
        mz[N > threshold]=0

    Nmax = max(N.flat) * 0.00075

    
    [z, r] = meshgrid(lz, lr)
    #quiver(r, z, mr, mz, eunit * Nmax, pivot='middle',
    #       width=1.0 * Nmax, headwidth=500, color="white")

    quiver(r, z, mr, mz, scale=Nmax / eunit, pivot='middle', color='w')


def draw_arrow( x, y, dx, dy, width=1, color = 'k' ):
     ax = gca()
     a = Arrow(x, y, dx, dy, width)
     a.set_edgecolor(color)
     a.set_facecolor(color)
     ax.add_patch(a)
     draw_if_interactive()
     return a

def plot_axis(var, grid, format=None, xlim=None, col=0, rescaling={}):
    if format is not None and '@@' in format:
        fmt, linewidth = format.split('@@')
        linewidth = float(linewidth)
    else:
        fmt = format
        linewidth = 1.5
        
    lz, m = load_axis (var, grid, col=col, rescaling=rescaling)
    if format:
        plot (lz, m, fmt, linewidth=linewidth)
    else:
        plot (lz, m, linewidth=1.5)
        
    
def save_fig(var, grid, ext, what="", ofile=None, **kwargs):
    if ofile is None:
        ofile = output_name(var, grid, ext, what=what)

    savefig(ofile, **kwargs)
    return ofile
    
def print_ret(f):
    def _f(*args, **kwargs):
        ret = f(*args, **kwargs)
        print "%s -> %s" % (f.__name__, repr(ret))
        return ret
    return _f



class FancyFormatter(Formatter):
    def __call__(self, x, pos=None):
        formatted = "%.3g" % x
        if 'e' in formatted:
            m, e = formatted.split('e')
            formatted = (r'$\mathdefault{%s}\cdot \mathdefault{10^{%d}}$'
                         % (m, int(e)))

        return formatted
