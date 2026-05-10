"""Tests that libras3d is functional for i.hyper.continuum standalone mode."""
import os, sys, tempfile
import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(__file__))
from test_ras3d_common import (
    WYVERN_PATH, TANAGER_PATH,
    skip_without_ras3d, skip_without_wyvern, skip_without_tanager,
    open_cube_checked, assert_band_valid, install_ras3d_shim, make_wl_sidecar,
)

@skip_without_ras3d
def test_shim_installs():
    install_ras3d_shim()
    import grass.script as gs
    assert hasattr(gs, 'raster3d_info')

@skip_without_ras3d
@skip_without_wyvern
def test_open_geotiff():
    import ras3d
    h, r = open_cube_checked(WYVERN_PATH)
    assert r['depths'] == 23
    ras3d.close_cube(h)

@skip_without_ras3d
@skip_without_tanager
def test_open_hdf5():
    import ras3d
    h, r = open_cube_checked(TANAGER_PATH)
    assert r['depths'] == 426
    assert r['rows'] == 732
    ras3d.close_cube(h)

@skip_without_ras3d
@skip_without_wyvern
def test_extract_z_slice_ras3d(tmp_path):
    install_ras3d_shim()
    os.environ['RAS3D_OUTDIR'] = str(tmp_path)
    sys.path.insert(0, '/home/yann/dev/i.hyper.continuum')
    from i_hyper_continuum import extract_z_slice
    extract_z_slice(WYVERN_PATH, '', 0, 'cont_slice_0')
    from ras3d_grass_shim import get_band_cache
    assert 'cont_slice_0' in get_band_cache()
    assert_band_valid(get_band_cache()['cont_slice_0'], 'continuum slice 0')

@skip_without_ras3d
@skip_without_wyvern
def test_wavelength_sidecar():
    install_ras3d_shim()
    import ras3d
    h, r = open_cube_checked(WYVERN_PATH)
    sidecar, _ = make_wl_sidecar(WYVERN_PATH, r['depths'])
    ras3d.close_cube(h)
    sys.path.insert(0, '/home/yann/dev/i.hyper.continuum')
    from i_hyper_continuum import get_all_band_wavelengths
    bands = get_all_band_wavelengths(WYVERN_PATH)
    assert len(bands) == r['depths']
    os.unlink(sidecar)
