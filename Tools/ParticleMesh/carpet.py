"""\
Generate a particle+spring mesh for a wave carpet.
"""

import itertools
import numpy as np


def carpet(nx, ny, plo, phi, pfile, mfile):
    """Generate a rectangular wave carpet with *nx* by *ny* particles from
    *plo* to *phi*."""

    dim = len(plo)
    plo = np.asarray(plo)
    phi = np.asarray(phi)
    dx  = (phi - plo) / np.asarray([nx-1, ny-1, 1])[:dim]

    ids = {}                    # particle id cache

    #
    # create particle file
    #
    with open(pfile, 'w') as particles:
        write = lambda x: particles.write(' '.join(map(str, list(x))) + '\n')
        write([nx*ny])

        for i, j in itertools.product(range(nx), range(ny)):
            x = plo[0] + i * dx[0]
            y = plo[1] + j * dx[1] if dim > 1 else 0
            z = plo[2]             if dim > 2 else 0
            r = np.asarray([x, y, z])[:dim]
            write(r)

            ids[i, j] = i*ny + j + 1


    #
    # create particle mesh (springs) file
    #
    with open(mfile, 'w') as springs:

        write = lambda x: springs.write(' '.join(map(str, list(x))) + '\n')

        def write_spring(ij1, ij2, l, k):
            id1 = ids.get(ij1, 0)
            id2 = ids.get(ij2, 0)
            if id1>0 and id2>0:
                write([id1, id2, l, k])

        # connect neighbours
        for i, j in itertools.product(range(nx), range(ny)):
            write_spring((i-1, j), (i, j), 1.0, 22.0)
            write_spring((i, j-1), (i, j), 1.0, 22.0)

        # connect next-neighbours
        for i, j in itertools.product(range(nx), range(ny)):
            write_spring((i-2, j), (i, j), 1.0, 22.0)
            write_spring((i, j-2), (i, j), 1.0, 22.0)


if __name__ == '__main__':
    carpet(10, 10, [-10.0, -10.0], [10.0, 10.0], 'particle_file', 'particle_mesh_file')
