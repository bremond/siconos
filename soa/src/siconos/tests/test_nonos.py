from __future__ import annotations
import numpy as np
import nonos as m

    
def test_version():
    assert m.disks.__version__ == "dev"
    
def test_make():
    data = m.disks.make_storage()
    
def test_add():
    data = m.disks.make_storage()
    disk = m.disks.add_disk(data)

def test_q():
    import numpy
    data = m.disks.make_storage()
    disk = m.disks.add_disk(data)
    disk.set_q(np.array([1,2,3]))
    q = disk.q()
    assert list(q) == [1,2,3]

def test_nslaw():
    data = m.disks.make_storage()
    nslaw = m.disks.add_nslaw(data)
    nslaw.set_e(0.7)
    assert nslaw.e() == 0.7

def test_diskdisk_r():
    data = m.disks.make_storage()
    relation = m.disks.add_diskdisk_r(data)
    relation.disk1().set_radius(1.0)
    relation.disk2().set_radius(2.0)
    assert relation.r1() == 1.0
    assert relation.r2() == 2.0

def test_time_stepping():

    
    data = m.disks.make_storage()
    simul = m.disks.add_time_stepping(data)
