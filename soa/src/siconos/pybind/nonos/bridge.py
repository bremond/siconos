import nonos as vkernel

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
        self._handle = vkernel.disks.add_interaction_manager(self.data())
        

    def insertLine(self, a, b , c):
        line = vkernel.disks.add_line_shape(self.data())
        line.set_a(a)
        line.set_b(b)
        line.set_c(c)

    def insertDisk(self, radius):
        disk = vkernel.disks.add_disk_shape(self.data())
        disk.set_radius(radius)
        
    def insertNonSmoothLaw(self, nslaw, gid1, gid2):
        self._handle.insert_nonsmooth_law(nslaw.handle(), gid1, gid2)


class NewtonImpactFrictionNSL(Stored):

    def __init__(self, e, not_used, mu, dimension):
        self._handle = vkernel.disks.add_nslaw(self.data())
        self._handle.set_e(e)
        self._handle.set_mu(mu)
        #self._handle.set_dimension(dimension)

class Osi(Stored, vkernel.disks.osi):

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
    
class Simulation(Stored):

    def __init__(self, nsds, timedisc):
        self._nsds = nsds
        self._timedisc = timedisc
        self._handle = vkernel.disks.add_simulation(self.data())

    def insertIntegrator(self, osi):
        pass # unimplemented

    def insertNonSmoothProblem(self, osnspb):
        pass # unimplemented

    def insertInteractionManager(self, interman):
        pass # unimplemented

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
        return 0. # unimplemented

    def nextTime(self):
        return 0.
    
    def hasNextEvent(self):
        return self.handle().has_next_event()

    def computeOneStep(self):
        return self.handle().compute_one_step()

    def clearNSDSChangeLog(self):
        pass # unimplemented

    def nextStep(self):
        pass # unimplemented
    
class Body(Stored):

    _ident = 0
    
    def __init__(self, position, velocity, mass, inertia):

        self._ident = self._ident + 1
        
        body = vkernel.disks.add_disk(self.data())
        self._body = body
        body.set_id(self._ident)
        body.set_q(position)
        body.set_velocity(velocity)
        body.set_mass(mass)
        body.set_inertia(inertia)


class OSNSPB(Stored):
    def __init__(self, dim, solvopts):

        self._handle = vkernel.disks.add_osnspb(self.data())

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
