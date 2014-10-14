
from collections import namedtuple

import numpy as np
import pylab

from matplotlib.path import Path
from matplotlib.patches import PathPatch

Spring = namedtuple('Spring', ['p1', 'p2', 'l', 'k'])

def draw_mesh(pfile, mfile):

    particles = []
    with open(pfile, 'r') as f:
        f.readline()
        for line in f:
            # project down to 2d by chopping off 3rd coord if present
            particles.append(np.asarray(map(float, line.split()))[:2])

    springs = []
    with open(mfile, 'r') as f:
        for line in f:
            r = line.split()
            springs.append(Spring(int(r[0]), int(r[1]), float(r[2]), float(r[3])))

    x, y = zip(*particles)
    pylab.plot(x, y, 'o')

    ax = pylab.gca()
    for spring in springs:
        p1 = particles[spring.p1-1]
        p2 = particles[spring.p2-1]
        nm = np.asarray([ p2[1]-p1[1], -(p2[0]-p1[0]) ])
        cp = 0.5*(p1+p2) + 0.2*nm
        verts = [ p1, cp, p2 ]
        codes = [ Path.MOVETO, Path.CURVE3, Path.CURVE3 ]
        patch = PathPatch(Path(verts, codes), facecolor='none', lw=2)
        ax.add_patch(patch)
    pylab.show()



if __name__ == '__main__':
    draw_mesh('particle_file', 'particle_mesh_file')
