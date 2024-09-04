#!/usr/bin/env python3
"""
Description: Export a Siconos mechanics-IO HDF5 file in VTK format.
"""

from vview import VView, VExportOptions
from mechanics_hdf5 import MechanicsHdf5

if __name__=='__main__':
    # Parse command-line
    opts = VExportOptions()
    opts.parse()

    ## Options and config already loaded above
    with MechanicsHdf5(io_filename=opts.io_filename, mode='r') as io:
        vview = VView(io, opts)
        vview.initialize_vtk()
        vview.export()