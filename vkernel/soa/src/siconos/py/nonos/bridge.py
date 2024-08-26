import nonos as vkernel
import numpy as np
from math import sqrt
import siconos.numerics as sn

def array(l):
    return np.array(l, dtype=np.float64)

class Stored():

    @classmethod
    def setStorage(cls, data):
        cls._data = data

    @classmethod
    def data(cls):
        return cls._data

    def handle(self):
        return self._handle

class SpaceFilter(Stored):

    def __init__(self, options):
        self._options = options
        self._initialized = False
        self._ngbh = vkernel.disks.add_neighborhood(self.data())
        self._ngbh.create(options.neighborhood_radius)
        self._interman = vkernel.disks.add_interaction_manager(self.data())
        self._handle = vkernel.disks.add_space_filter(self.data())
        self._handle.set_neighborhood(self._ngbh)
        self._handle.set_diskdisk_r(vkernel.disks.add_diskdisk_r(self.data()))
        self._fdisks = {}

    def insertSegment(self, x1, y1, x2, y2):
        segment = vkernel.disks.add_segment_shape(self.data())
        segment.set_p1(array([x1, y1, 0]))
        segment.set_p2(array([x2, y2, 0]))

        mp = int(max(3, sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)) / self._options.min_radius))

        segment.set_maxpoints(mp) # fix / size of smallest disk
        segment.initialize()

        disksegment = vkernel.disks.add_disksegment_r(self.data())
        disksegment.set_segment(segment)

        self.handle().insert_disksegment_r(disksegment)

        return segment

    def insertLine(self, a, b , c):
        line = vkernel.disks.add_line_shape(self.data())
        line.set_a(a)
        line.set_b(b)
        line.set_c(c)
        line.set_maxpoints(10000)
        line.initialize()

        print("new line,  p0:", line.p0())
        print("--------, dir:", line.direction())

        diskline = vkernel.disks.add_diskline_r(self.data())
        diskline.set_line(line)

        self.handle().insert_line(diskline)

        return line

    def insertTranslatedDisk(self, radius, translation):

        disk_shape = None
        if radius in self._fdisks:
            disk_shape = self._fdisks[radius]
        else:
            disk_shape = vkernel.disks.add_disk_shape(self.data())
            disk_shape.set_radius(radius)
            self._fdisks[radius] = disk_shape

        translated_disk_shape = \
            vkernel.disks.add_translated_disk_shape(self.data())
        translated_disk_shape.set_item(disk_shape)
        translated_disk_shape.set_translation(array([translation[0],translation[1],0]))

        diskfdisk = vkernel.disks.add_diskfdisk_r(self.data())
        diskfdisk.set_translated_disk_shape(translated_disk_shape)

        self.handle().insert_diskfdisk_r(diskfdisk)

        return translated_disk_shape

    def insertNonSmoothLaw(self, nslaw, gid1, gid2):
        self._interman.insert_nonsmooth_law(nslaw.handle(), gid1, gid2)
        self._handle.set_nslaw(self._interman.get_nonsmooth_law(gid1, gid2)) # one nslaw!!

    def updateInteractions(self):
        if not self._initialized:
            self._handle.make_points()
            self._ngbh.add_point_sets(0)
            self._ngbh.set_active(0, 0, True)
            self._ngbh.set_active(0, 1, True)
            self._ngbh.set_active(0, 2, True)
            self._ngbh.set_active(1, 1, False)
            self._ngbh.set_active(2, 2, False)
            self._ngbh.set_active(1, 0, False)
            self._ngbh.set_active(2, 0, False)
            self._ngbh.set_active(1, 2, False)
            self._ngbh.set_active(2, 1, False)

            self._initialized = True

        self._ngbh.update(0)
        self._ngbh.search()
        self._handle.update_index_set0(0);
        self._ngbh.sort()

    def removeStaticBody(self, body):
        self._ngbh.search() # needed after sort

        # only segments are removable
        if type(body) == type([]):
            # box2d
            self.handle().remove_static_segment(0, body[0])
            self.handle().remove_static_segment(0, body[1])
            self.handle().remove_static_segment(0, body[2])
            self.handle().remove_static_segment(0, body[3])
        else:
            self.handle().remove_static_segment(0, body)


class NewtonImpactFrictionNSL(Stored):

    def __init__(self, e, not_used, mu, dimension):
        self._handle = vkernel.disks.add_nslaw(self.data())
        self._handle.set_e(e)
        self._handle.set_mu(mu)
        #self._handle.set_dimension(dimension)

class Osi(Stored):

    def __init__(self, theta):
        self._handle = vkernel.disks.add_osi(self.data())
        self._handle.set_theta(theta)

    def setConstraintActivationThreshold(self, cat):
        pass # unimplemented

    def setGamma(self, cat):
        pass # unimplemented


class Topology(Stored):

    def __init__(self):
        self._handle = vkernel.disks.add_topology(self.data())

    def indexSetsSize(self):
        return 2

class NonSmoothDynamicalSystem(Stored):

    def __init__(self, t0, T):
        self._t0 = t0
        self._T = T
        self._topology = Topology()
        self._mapid = {}

    def topology(self):
        return self._topology

    def insertDynamicalSystem(self, body):
        self._mapid[int(body.number())] = body
#        self.topology().handle().set_dynamical_system_id(body.handle(), body.number())

    def dynamicalSystem(self, ds_id):
        return self._mapid[ds_id]



class TimeDiscretisation(Stored):

    def __init__(self, t0, h):
        self._t0 = t0
        self._h = h
        self._handle = \
            vkernel.disks.add_time_discretization(self.data())
        self.handle().set_t0(t0)
        self.handle().set_h(h)

class Simulation(Stored):

    def __init__(self, nsds, timedisc):
        self._need_init = True
        self._nsds = nsds
        self._timedisc = timedisc
        self._timedisc.handle().set_tmax(self._nsds._T) # vkernel does not have nsds
        self._handle = vkernel.disks.add_simulation(self.data())
        self.handle().initialize()
        self.handle().one_step_integrator().set_theta(0.50001)

    def insertIntegrator(self, osi):
        pass # unimplemented

    def insertNonSmoothProblem(self, osnspb):
        pass # unimplemented

    def insertInteractionManager(self, interman):
        self._interman = interman

    def setNewtonOptions(self, nopts):
        pass # unimplemented

    def setNewtonMaxIteration(self, maxiter):
        pass # unimplemented

    def setNewtonTolerance(self, newtontol):
        pass # unimplemented

    def setSkipLastUpdateOutput(self, skipluo):
        pass # unimplemented

    def setSkipLastUpdateInput(self, skiplui):
        pass # unimplemented

    def setSkipResetLambdas(self, skipresetlbds):
        pass # unimplemented

    def setDisplayNewtonConvergence(self, dnc):
        pass # unimplemented

    def startingTime(self):
        return self.handle().time_discretization().t0()

    def nextTime(self):
        return self.startingTime() + self.handle().current_step() *\
            self.handle().time_discretization().h()

    def hasNextEvent(self):
        return self.handle().has_next_event()

    def updateInteractions(self):
        self._interman.updateInteractions()

    def computeOneStep(self):
        self._interman.updateInteractions()
        return self.handle().compute_one_step()

    def clearNSDSChangeLog(self):
        pass # unimplemented

    def nextStep(self):
        pass # unimplemented

class Body(Stored):

    __count = 0
    __disk_shapes = {}

    def __init__(self, radius, mass, position, velocity):

        Body.__count += 1
        self._ident = Body.__count

        body = vkernel.disks.add_disk(self.data())
        self._handle = body
        body.set_id(self._ident)
        body.set_q(array(position))
        body.set_velocity(array(velocity))
        body.set_mass_matrix(array([mass, mass, mass*radius*radius/2]))

        disk_shape = None
        if radius in Body.__disk_shapes:
            disk_shape = Body.__disk_shapes[radius]
        else:
            disk_shape = vkernel.disks.add_disk_shape(self.data())
            disk_shape.set_radius(radius)
            Body.__disk_shapes[radius] = disk_shape

        body.set_shape(disk_shape)
        body.set_fext(array([0,0,0])) # default

    def scalarMass(self):
        return self.handle().mass_matrix()[0]

    def setFExtPtr(self, fext):
        self.handle().set_fext(array(fext))

    def setNumber(self, num):
        self.handle().set_id(num)

    def number(self):
        return self.handle().id()

    def setQ0Ptr(self, pos):

        # FIX: next step not necessary 1
        self.handle().set_q_at_step(array(pos), 0)
        self.handle().set_q_at_step(array(pos), 1)

    def setVelocity0Ptr(self, vel):

        # FIX: next step not necessary 1
        self.handle().set_velocity_at_step(array(vel), 0)
        self.handle().set_velocity_at_step(array(vel), 1)

    def resetToInitialState(self):
        pass # compatibility

    def swapInMemory(self):
        pass # compatibility

class OSNSPB(Stored):
    def __init__(self, dim, solvopts):

        self._so = vkernel.disks.add_solver_options(self.data())

        if solvopts is not None:
            self._so.create(solvopts.solverId)
            for i,v in enumerate(solvopts.iparam):
                self._so.set_iparam(i, v)

            for i,v in enumerate(solvopts.dparam):
                self._so.set_dparam(i, v)

        else:
            # default
            self._so.create(sn.SICONOS_FRICTION_2D_NSGS)
            self._so.set_iparam(sn.SICONOS_IPARAM_MAX_ITER, 100)
            self._so.set_dparam(sn.SICONOS_DPARAM_TOL, 1e-3)
            self._so.set_dparam(sn.SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT, 10)

        self._fc2d = vkernel.disks.add_fc2d(self.data())
        self._handle = vkernel.disks.add_osnspb(self.data())
        self._handle.set_options(self._so)
        self._fc2d.create(self._so.solver_id())
        self.handle().set_problem(self._fc2d)
        self.handle().set_mu(0.5)
        self.handle().set_verbose(False)
#        self._fc2d.instance().dimension = 2

    def setMaxSize(self, maxs):
        self._maxSize = maxs

    def setMStorageType(self, stype):
        self._storageType = stype

    def setNumericsVerboseMode(self, vm):
        self._numericsVerboseMode = vm

    def setKeepLambdaAndYState(self, kly):
        self._keepLambdaAndYState = kly

    def setAssemblyType(self, at):
        pass

    def getSizeOutput(self):
        return 0 # unimplemented


class SICONOS_TS_NONLINEAR():
    pass

MoreauJeanOSI = Osi

class Unimplemented():

    def __init__(self, *args):
        assert False

MoreauJeanGOSI = Unimplemented
TimeSteppingDirectProjection = Unimplemented

FrictionContact = OSNSPB

class MechanicsIO(Stored):

    def __init__(self):
        self._handle = vkernel.disks.add_io(self.data())

    def positions(self, nsds):
        return self.handle().positions(0)

    def velocities(self, nsds):
        return self.handle().velocities(0)

    def contactPoints(self, nsds, output_contact_index_set):
        return self.handle().contact_points(0)

class SpaceFilterOptions():

    neighborhood_radius = 2.1
    min_radius = 0.5

    def toJson(self):
        return json.dumps(self, default=lambda o: o.__dict__)
