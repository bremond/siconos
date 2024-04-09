from siconos.mechanics.collision.tools import Contactor
from nonos.mechanics_run import MechanicsHdf5Runner

import nonos
import siconos.numerics as sn
import siconos.kernel as sk

import siconos

import numpy
from math import sqrt

import random
nonos.mechanics_run.set_backend('vnative')

disk_radius = 1

with MechanicsHdf5Runner() as io:

    io.add_primitive_shape('PointR', 'Disk', [0.1])
    io.add_primitive_shape('DiskR', 'Disk', [disk_radius])

    io.add_primitive_shape('Ground1', 'Line', (-30, 100, 0))
#    io.add_primitive_shape('Ground2', 'Line', (-10, 10, 0))

    io.add_Newton_impact_friction_nsl('contact', mu=0.5, e=0)

    # io.add_object('disk0', [Contactor('DiskR')],
    #               translation=[8, 8],
    #               orientation=[0], velocity=[0, 0, -5], mass=1)

    # io.add_object('disk1', [Contactor('DiskR')],
    #               translation=[3, 5],
    #               orientation=[0], velocity=[0, 0, 0], mass=1)

    io.add_object('ground1', [Contactor('Ground1')],
                  translation=[0, 0])

    # io.add_object('ground2', [Contactor('Ground2')],
    #               translation=[0, 0])

    # io.add_object('ground3', [Contactor('Ground3')],
    #               translation=[0, 0])

    N = 20
    for i in range(N):
        for j in range(N):
            io.add_object('disk{}-{}'.format(i,j), [Contactor('DiskR')],
                          translation=[2*i-N, 2*j+0.33*N],
                          orientation=[0], velocity=[0, 0, 0], mass=1)


    # io.add_object('disk2', [Contactor('DiskR')],
    #               translation=[0.5, 4+3*disk_radius],
    #               orientation=[0], velocity=[0, 0, 0], mass=1)


options = sk.solver_options_create(sn.SICONOS_FRICTION_2D_NSGS)
options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 100000
options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-12

def noforces(body):
    pass

with MechanicsHdf5Runner(mode='r+') as io:

        io.run(with_timer=False,
               gravity_scale=1,
               t0=0,
               T=20,
               h=0.005,
               theta=0.50001,
               Newton_max_iter=1000,
               set_external_forces=None,
               solver_options=options,
               numerics_verbose=False,
               output_contact_forces=False,
               output_frequency=None)
