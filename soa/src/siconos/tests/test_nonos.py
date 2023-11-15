from __future__ import annotations

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
    q = disk.q()
    q = [1,2,3]
    assert list(q) == list(disk.q())

def test_nslaw():
    data = m.disks.make_storage()
    nslaw = m.disks.add_nslaw(data)
    e = nslaw.e()
    e=0.7
    assert nslaw.e() == 0.7
