from siconos.mechanics.collision.tools import Contactor
from nonos.mechanics_run import MechanicsHdf5Runner

import nonos
import siconos.numerics as sn
import siconos.kernel as sk

import siconos

import numpy
from math import sqrt

nonos.mechanics_run.set_backend('vnative')

disk_radius = 1
N = 10

with MechanicsHdf5Runner() as io:

    io.add_primitive_shape('DiskR1', 'Disk', [disk_radius])
    io.add_primitive_shape('DiskR2', 'Disk', [0.1*disk_radius])

    io.add_Newton_impact_friction_nsl('contact', mu=0.3, e=0)

    for i in range(N):
        for j in range(N):
            io.add_object('disk-{}-{}'.format(i,j), [Contactor('DiskR2')],
                          translation=[0.2*i, 3*disk_radius+0.2*j],
                          orientation=[0], velocity=[0, 0, 0], mass=1, inertia=0.5)

    io.add_object('fixed-disk1', [Contactor('DiskR1')], translation=[0, 0])

    io.add_object('fixed-disk2', [Contactor('DiskR1')], translation=[2, -5])                                 

options = sk.solver_options_create(sn.SICONOS_FRICTION_2D_NSGS)
options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 10
options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-3


with MechanicsHdf5Runner(mode='r+') as io:

        io.run(with_timer=False,
               gravity_scale=1,
               t0=0,
               T=10,
               h=0.005,
               theta=0.50001,
               Newton_max_iter=1000,
               set_external_forces=None,
               solver_options=options,
               numerics_verbose=False,
               output_contact_forces=False,
               output_frequency=None)
