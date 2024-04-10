import nonos as vkernel
import numpy as np


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
        self._handle.set_nslaw(self._interman.get_nonsmooth_law(0,0)) # one nslaw!!
        

    def insertLine(self, a, b , c):
        line = vkernel.disks.add_line_shape(self.data())
        line.set_a(a)
        line.set_b(b)
        line.set_c(c)
        line.set_maxpoints(5000)
        line.initialize()

        print("new line,  p0:", line.p0())
        print("--------, dir:", line.direction())

        diskline = vkernel.disks.add_diskline_r(self.data())
        diskline.set_line(line)

        self.handle().insert_line(diskline)

    def insertDisk(self, radius):
        disk_shape = vkernel.disks.add_disk_shape(self.data())
        disk_shape.set_radius(radius)
        
    def insertNonSmoothLaw(self, nslaw, gid1, gid2):
        self._interman.insert_nonsmooth_law(nslaw.handle(), gid1, gid2)

    def updateInteractions(self):
        if not self._initialized:
            self._handle.make_points()
            self._ngbh.add_point_sets(0)
            self._initialized = True
            
        self._ngbh.update(0)
        self._ngbh.search()
        self._handle.update_index_set0(0);

class NewtonImpactFrictionNSL(Stored):

    def __init__(self, e, not_used, mu, dimension):
        self._handle = vkernel.disks.add_nslaw(self.data())
        self._handle.set_e(e)
        self._handle.set_mu(mu)
        #self._handle.set_dimension(dimension)

class Osi(Stored):

    def __init__(self, theta):
        self._handle = vkernel.disks.add_osi(data())
        self._handle.set_theta(theta)

    def setConstraintActivationThreshold(self, cat):
        pass # unimplemented

    def setGamma(self, cat):
        pass # unimplemented


class Topology(Stored):

    def __init_(self):
        self._handle = vkernel.disks.add_topology(self.data())

    def indexSetSize(self):
        pass
    
class NonSmoothDynamicalSystem(Stored):

    def __init__(self, t0, T):
        self._t0 = t0
        self._T = T
        self._topology = Topology()

    def topology(self):
        return self._topology

    def insertDynamicalSystem(self, body):
        pass # compatibility

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

    def computeOneStep(self):
        self._interman.updateInteractions()
        return self.handle().compute_one_step()

    def clearNSDSChangeLog(self):
        pass # unimplemented

    def nextStep(self):
        pass # unimplemented
    
class Body(Stored):

    _ident = 1
    
    def __init__(self, radius, mass, position, velocity):

        self._ident = self._ident + 1

        body = vkernel.disks.add_disk(self.data())
        self._handle = body
        body.set_id(self._ident)
        body.set_q(array(position))
        body.set_velocity(array(velocity))
        body.set_mass_matrix(array([mass, mass, mass*radius*radius/2]))
        disk_shape = vkernel.disks.add_disk_shape(self.data())
        body.set_shape(disk_shape)
        body.shape().set_radius(radius)
        body.set_fext(array([0,0,0])) # default 

    def scalarMass(self):
        return self.handle().mass_matrix()[0]

    def setFExtPtr(self, fext):
        self.handle().set_fext(array(fext))

    def setNumber(self, num):
        self.handle().set_id(num)

class OSNSPB(Stored):
    def __init__(self, dim, solvopts):

        self._so = vkernel.disks.add_solver_options(self.data())
        self._so.create(400)
        self._so.set_iparam(0, 20)
        self._fc2d = vkernel.disks.add_fc2d(self.data())
        self._handle = vkernel.disks.add_osnspb(self.data())
        self._handle.set_options(self._so)
        self._fc2d.create(400)
        self.handle().set_problem(self._fc2d)
        self.handle().set_mu(0.3)
        self.handle().set_verbose(True)
#        self._fc2d.instance().dimension = 2
        
    def setMaxSize(self, maxs):
        self._maxSize = maxs

    def setMStorageType(self, stype):
        self._storageType = stype

    def setNumericsVerboseMode(self, vm):
        self._numericsVerboseMode = vm

    def setKeepLambdaAndYState(self, kly):
        self._keepLambdaAndYState = kly

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

class SpaceFilterOptions():

    neighborhood_radius = 2.5
