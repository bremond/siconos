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
