#!/usr/bin/env python
# ── ras3d standalone detection ────────────────────────────────────────────────
import os as _os
_RAS3D = False
if not _os.environ.get('GISBASE'):
    try:
        import importlib.util as _ilu
        if _ilu.find_spec('ras3d') and _ilu.find_spec('ras3d_grass_shim'):
            from ras3d_grass_shim import install as _r3_install
            _r3_install()
            _RAS3D = True
    except Exception:
        pass
# ─────────────────────────────────────────────────────────────────────────────
##############################################################################
# MODULE:    i.hyper.continuum
# AUTHOR(S): Created for hyperspectral continuum removal
# PURPOSE:   Apply continuum removal to hyperspectral imagery
# COPYRIGHT: (C) 2025 by the GRASS Development Team
# SPDX-License-Identifier: GPL-2.0-or-later
##############################################################################

# %module
# % description: Apply continuum removal to hyperspectral imagery for absorption feature analysis
# % keyword: imagery
# % keyword: hyperspectral
# % keyword: continuum removal
# % keyword: spectral analysis
# %end

# %option G_OPT_R3_INPUT
# % key: input
# % required: yes
# % description: Input hyperspectral 3D raster map (from i.hyper.import)
# % guisection: Input
# %end

# %option G_OPT_R3_OUTPUT
# % key: output
# % required: yes
# % description: Output continuum-removed 3D raster map
# % guisection: Output
# %end

# %option
# % key: min_wavelength
# % type: double
# % required: no
# % description: Minimum wavelength for processing (nanometers). Default: use minimum available wavelength
# % guisection: Wavelength Range
# %end

# %option
# % key: max_wavelength
# % type: double
# % required: no
# % description: Maximum wavelength for processing (nanometers). Default: use maximum available wavelength
# % guisection: Wavelength Range
# %end

# %option
# % key: method
# % type: string
# % required: no
# % options: convex_hull,linear,polynomial
# % answer: convex_hull
# % description: Method for continuum calculation
# % guisection: Processing
# %end

# %option
# % key: normalization
# % type: string
# % required: no
# % options: division,subtraction
# % answer: division
# % description: Normalise the spectrum against the continuum by division (hull quotient, bounded 0-1) or subtraction (residual, keeps absolute depth)
# %end

# %option
# % key: max_gap
# % type: double
# % required: no
# % answer: 0
# % description: Width in nm above which a step between band centres is a sensor gap, so the continuum is computed on each side separately. 0 detects gaps as 3x the median band spacing; a negative value disables the split
# %end

# %option
# % key: polynomial_order
# % type: integer
# % required: no
# % answer: 2
# % description: Polynomial order for polynomial method (2-5)
# % guisection: Processing
# %end

# %flag
# % key: n
# % description: Only process bands marked as valid (valid=1)
# % guisection: Processing
# %end

# %flag
# % key: k
# % description: Keep original band values outside the specified wavelength range
# % guisection: Processing
# %end

# %flag
# % key: i
# % description: Print information about bands and processing without creating output
# % guisection: Processing
# %end

import sys
import os
import ctypes
import ctypes.util
import re
import grass.script as gs
import numpy as np


def _load_g3d_lib():
    """Load libgrass_g3d via ctypes, return None on failure.

    The GIS library has to be initialised before any GRASS library call, or
    every one of them reports "System not initialized" and fails. G_gisinit
    is a macro, so the symbol to call is G__no_gisinit, and it takes the
    headers version it was built against, which is read from version.h in
    the running GISBASE rather than hard-coded. This mirrors what
    i.hyper.sleuth does for the same library.
    """
    gisbase = os.environ.get("GISBASE")
    if not gisbase:
        return None
    lib_dir = os.path.join(gisbase, "lib")
    try:
        version_header = os.path.join(gisbase, "include", "grass", "version.h")
        with open(version_header) as handle:
            match = re.search(r'#define\s+GRASS_HEADERS_VERSION\s+"([^"]+)"',
                              handle.read())
        headers_version = match.group(1).encode() if match else b""

        gis_lib = ctypes.CDLL(os.path.join(lib_dir, "libgrass_gis.so"),
                              mode=ctypes.RTLD_GLOBAL)
        gis_lib.G__no_gisinit.restype = None
        gis_lib.G__no_gisinit.argtypes = [ctypes.c_char_p]
        gis_lib.G__no_gisinit(headers_version)

        lib = ctypes.CDLL(os.path.join(lib_dir, "libgrass_g3d.so"))
        fn = lib.Rast3d_extract_z_slice
        fn.restype = ctypes.c_int
        fn.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p]
        return fn
    except (OSError, AttributeError, KeyError):
        return None


_rast3d_extract_z_slice = _load_g3d_lib()


def extract_z_slice(name3d, mapset3d, z, name2d):
    """Extract Z-slice z (0-based) from 3D raster into 2D raster name2d.

    Uses Rast3d_extract_z_slice() — tile-bulk read via RASTER3D_NO_CACHE,
    one disk read per tile versus one function call per voxel.
    Falls back to g.copy if the C binding is unavailable.
    """
    if _RAS3D:
        import ras3d as _r3, ras3d_write as _r3w
        _h = _r3.open_cube(name3d)
        _arr = _r3.get_band(_h, z)
        from ras3d_grass_shim import get_band_cache
        get_band_cache()[name2d] = _arr
        _r3w.write_raster2d(_r3w.outpath(name2d), _arr, _h)
        _r3.close_cube(_h)
        return
    if _rast3d_extract_z_slice is not None:
        ret = _rast3d_extract_z_slice(
            name3d.encode(),
            mapset3d.encode() if mapset3d else b"",
            ctypes.c_int(z),
            name2d.encode(),
        )
        if ret != 0:
            gs.fatal(f"Rast3d_extract_z_slice failed for {name3d} z={z}")
    else:
        # Fallback: band_num is 1-based, z is 0-based
        band_name = f"{name3d}#{z + 1}" if not mapset3d else f"{name3d}@{mapset3d}#{z + 1}"
        gs.run_command("g.copy", raster=f"{band_name},{name2d}", quiet=True, overwrite=True)


def get_raster3d_info(raster3d):
    """Get information about 3D raster"""
    try:
        info = gs.raster3d_info(raster3d)
        return info
    except Exception as e:
        gs.fatal(f"Cannot get info for 3D raster {raster3d}: {e}")


def parse_wavelength_from_metadata(raster3d, band_num):
    """Parse wavelength and validity from band metadata"""
    band_name = f"{raster3d}#{band_num}"
    wavelength = None
    fwhm = None
    valid = True
    unit = "nm"
    
    try:
        result = gs.read_command('r.support', map=band_name, flags='n')
        
        for line in result.split('\n'):
            line = line.strip()
            if line.startswith('wavelength='):
                wavelength = float(line.split('=')[1])
            elif line.startswith('FWHM='):
                fwhm = float(line.split('=')[1])
            elif line.startswith('valid='):
                valid = int(line.split('=')[1]) == 1
            elif line.startswith('unit='):
                unit = line.split('=')[1].strip()
    except:
        pass
    
    return wavelength, fwhm, valid, unit


def convert_wavelength_to_nm(wavelength, unit):
    """Convert wavelength to nanometers"""
    unit = unit.lower().strip()
    
    if unit in ['nm', 'nanometer', 'nanometers']:
        return wavelength
    elif unit in ['um', 'µm', 'micrometer', 'micrometers', 'micron', 'microns']:
        return wavelength * 1000.0
    elif unit in ['m', 'meter', 'meters']:
        return wavelength * 1e9
    else:
        gs.warning(f"Unknown wavelength unit '{unit}', assuming nanometers")
        return wavelength


def _load_hyper_json_bands(raster3d):
    """Read wavelength/fwhm/validity from i.hyper.import's JSON sidecar.

    i.hyper.import (HyperMetadata) stores band metadata at
    $MAPSET/grid3/<mapname>/hyper.json rather than in r3.support history,
    so that must be checked before falling back to r.support probing.
    """
    import json as _json

    name, mapset = (raster3d.split('@', 1) if '@' in raster3d
                     else (raster3d, None))
    try:
        env = gs.gisenv()
        mapset = mapset or env['MAPSET']
        path = os.path.join(env['GISDBASE'], env['LOCATION_NAME'], mapset,
                            'grid3', name, 'hyper.json')
    except Exception:
        return []

    if not os.path.isfile(path):
        return []

    with open(path) as _fj:
        data = _json.load(_fj)

    b = data.get('bands') or {}
    wavelengths = b.get('wavelength')
    if not wavelengths:
        return []
    fwhms = b.get('fwhm') or []
    valids = b.get('validity') or []

    return [{'band_num': i + 1, 'wavelength': float(wl),
             'fwhm': float(fwhms[i]) if i < len(fwhms) else 0,
             'valid': bool(valids[i]) if i < len(valids) else True,
             'unit': 'nm'}
            for i, wl in enumerate(wavelengths)]


def _read_wl_sidecar(raster3d):
    """Read the ras3d ``<cube>.wl.json`` sidecar, or return [].

    The canonical writer is ``ras3d_write`` (and *i.hyper.atcorr* in
    direct-file mode): the path is the cube path with its extension
    stripped plus ``.wl.json``, and the payload is a flat JSON list of
    band-centre wavelengths. A list of dicts carrying ``wavelength`` is
    accepted too, so a richer sidecar written later still loads.

    :param str raster3d: cube path or 3D raster name
    :return: list of ``(wavelength_nm, fwhm_or_None, valid)`` tuples
    :rtype: list[tuple]
    """
    import json as _json

    for _sfx in ('.tif', '.tiff', '.h5', '.hdf5', '.vrt', '.img', ''):
        _base = (raster3d.removesuffix(_sfx) if _sfx and raster3d.endswith(_sfx)
                 else raster3d)
        _wlp = _base + '.wl.json'
        if not os.path.isfile(_wlp):
            continue
        with open(_wlp) as _fj:
            _data = _json.load(_fj)
        _out = []
        for _entry in _data:
            if isinstance(_entry, dict):
                _wl = float(_entry['wavelength'])
                _fwhm = (float(_entry['fwhm'])
                         if _entry.get('fwhm') is not None else None)
                _valid = bool(_entry.get('valid', True))
            else:
                _wl, _fwhm, _valid = float(_entry), None, True
            # ras3d_write stores nm, but i.hyper.atcorr has written um.
            if _wl < 10.0:
                _wl *= 1000.0
            _out.append((_wl, _fwhm, _valid))
        if _out:
            return _out
    return []


def get_all_band_wavelengths(raster3d, only_valid=False):
    """Extract all band wavelengths and metadata from 3D raster"""
    json_bands = _load_hyper_json_bands(raster3d)
    if json_bands:
        if only_valid:
            json_bands = [b for b in json_bands if b['valid']]
        json_bands.sort(key=lambda x: x['wavelength'])
        if json_bands:
            return json_bands

    # A cube opened straight from a file has no GRASS metadata, so fall back
    # to the ras3d .wl.json sidecar before probing r3.info history.
    sidecar = [{'band_num': i, 'wavelength': wl,
                'fwhm': fwhm if fwhm is not None else 0,
                'valid': valid, 'unit': 'nm'}
               for i, (wl, fwhm, valid) in
               enumerate(_read_wl_sidecar(raster3d), start=1)]
    if only_valid:
        sidecar = [b for b in sidecar if b['valid']]
    if sidecar:
        sidecar.sort(key=lambda x: x['wavelength'])
        return sidecar

    info = get_raster3d_info(raster3d)
    depths = int(info['depths'])
    
    bands = []
    
    gs.verbose(f"Scanning {depths} bands for wavelength metadata...")
    
    for i in range(1, depths + 1):
        wavelength, fwhm, valid, unit = parse_wavelength_from_metadata(raster3d, i)
        
        if wavelength is not None:
            wavelength_nm = convert_wavelength_to_nm(wavelength, unit)
            
            if only_valid and not valid:
                gs.verbose(f"Band {i}: {wavelength_nm} nm - SKIPPED (invalid)")
                continue
            
            bands.append({
                'band_num': i,
                'wavelength': wavelength_nm,
                'fwhm': fwhm if fwhm else 0,
                'valid': valid,
                'unit': unit
            })
            
            gs.verbose(f"Band {i}: {wavelength_nm} nm (FWHM: {fwhm}, valid: {valid})")
    
    if not bands:
        gs.fatal("No wavelength metadata found in 3D raster bands. "
                "Please use data imported with i.hyper.import or add wavelength metadata.")
    
    # Sort bands by wavelength
    bands.sort(key=lambda x: x['wavelength'])
    
    return bands


def filter_bands_by_wavelength(bands, min_wl, max_wl):
    """Filter bands within specified wavelength range"""
    # Get all bands for output structure
    all_bands = bands.copy()
    
    # Filter processing bands
    processing_bands = [b for b in bands if min_wl <= b['wavelength'] <= max_wl]
    
    gs.message(f"Wavelength range for continuum removal: {min_wl:.1f} - {max_wl:.1f} nm")
    gs.message(f"Total bands in input: {len(all_bands)}")
    gs.message(f"Bands to process: {len(processing_bands)}")
    
    if not processing_bands:
        gs.fatal(f"No bands found in wavelength range {min_wl:.1f} - {max_wl:.1f} nm")
    
    return all_bands, processing_bands


def calculate_convex_hull_continuum(wavelengths, reflectances):
    """Calculate continuum as the upper convex hull of the spectrum.

    Monotone chain in wavelength order: a vertex is dropped as soon as a
    later point shows it is not on the upper boundary. The first and last
    band are always kept, so the continuum spans the whole spectrum and no
    value is extrapolated.

    The hull has to be able to descend. An envelope assembled from
    monotonically rising points stops at the reflectance maximum, and every
    band beyond it is then divided by a continuum extrapolated far above the
    spectrum. On a VNIR-SWIR cube that is the entire SWIR, where the
    diagnostic mineral and polymer features are, and it turns the falling
    limb into one huge spurious absorption: measured on a drone cube, a
    vegetation spectrum peaking at 900 nm gave a continuum of 5197 against a
    spectrum of 390 at 2400 nm, and a band depth of 0.702 at 1716 nm where
    the true upper hull gives 0.447.
    """
    wl = np.asarray(wavelengths, dtype=float)
    refl = np.asarray(reflectances, dtype=float)

    finite = np.isfinite(wl) & np.isfinite(refl)
    if finite.sum() < 2:
        return calculate_linear_continuum(wavelengths, reflectances)

    order = np.argsort(wl[finite])
    xs = wl[finite][order]
    ys = refl[finite][order]

    hull = []
    for x, y in zip(xs, ys):
        # Drop the last vertex while it sits on or below the line joining
        # the one before it to the new point, which means it is not on the
        # upper hull.
        while len(hull) >= 2:
            (x1, y1), (x2, y2) = hull[-2], hull[-1]
            if (x2 - x1) * (y - y1) - (y2 - y1) * (x - x1) < 0:
                break
            hull.pop()
        hull.append((float(x), float(y)))

    if len(hull) < 2:
        return calculate_linear_continuum(wavelengths, reflectances)

    return np.interp(wl, [p[0] for p in hull], [p[1] for p in hull])


def calculate_linear_continuum(wavelengths, reflectances):
    """Calculate linear continuum between endpoints"""
    continuum = np.linspace(reflectances[0], reflectances[-1], len(wavelengths))
    return continuum


def calculate_polynomial_continuum(wavelengths, reflectances, order=2):
    """Calculate polynomial continuum"""
    try:
        # Normalize wavelengths for better numerical stability
        wl_mean = np.mean(wavelengths)
        wl_std = np.std(wavelengths)
        wl_norm = (wavelengths - wl_mean) / wl_std
        
        # Fit polynomial
        coeffs = np.polyfit(wl_norm, reflectances, order)
        continuum = np.polyval(coeffs, wl_norm)
        
        # Ensure continuum doesn't go below reflectance (optional)
        # continuum = np.maximum(continuum, reflectances)
        
        return continuum
    except Exception as e:
        gs.warning(f"Polynomial fitting failed: {e}, using linear continuum")
        return calculate_linear_continuum(wavelengths, reflectances)


def spectral_segments(wavelengths, max_gap=0.0):
    """Index ranges of contiguous spectral coverage.

    A continuum drawn across a region the sensor does not sample is an
    interpolation between two unrelated shoulders, and dividing by it
    manufactures an absorption that was never measured. Measured on the C2RD
    drone cube, whose 1784.5-1984.0 nm gap swallows the 1940 nm water band, a
    hull spanning the gap returned depths of 0.64 to 0.70 on water targets
    where nothing had been observed at all.

    :param wavelengths: band centres in nm, ascending.
    :param max_gap: a step wider than this starts a new segment. 0 uses three
        times the median step, which separates real gaps from the uneven
        spacing of a VNIR-SWIR join. Negative disables the split.
    :return: list of (start, stop) index pairs, stop exclusive.
    """
    wl = np.asarray(wavelengths, dtype=float)
    if wl.size < 2 or max_gap < 0:
        return [(0, wl.size)]
    steps = np.diff(wl)
    cut = max_gap if max_gap > 0 else 3.0 * float(np.median(steps))
    breaks = [0, *(np.where(steps > cut)[0] + 1).tolist(), wl.size]
    return [(breaks[i], breaks[i + 1]) for i in range(len(breaks) - 1)]


def calculate_continuum(wavelengths, reflectances, method, polynomial_order):
    """Dispatch to the chosen continuum, on one contiguous segment."""
    if method == 'convex_hull':
        return calculate_convex_hull_continuum(wavelengths, reflectances)
    if method == 'polynomial':
        return calculate_polynomial_continuum(wavelengths, reflectances,
                                              polynomial_order)
    return calculate_linear_continuum(wavelengths, reflectances)


def remove_continuum_spectrum(wavelengths, reflectances, method,
                              polynomial_order, normalization, segments):
    """Continuum-remove one spectrum, each coverage segment on its own.

    :param normalization: ``division`` gives the hull quotient, bounded by 1
        for a continuum at or above the spectrum; ``subtraction`` gives the
        residual, which keeps the absolute size of the absorption. Both are
        in Clark and Roush (1984); division is what ENVI reports.
    :return: the continuum-removed spectrum, null where it is undefined.
    """
    wl = np.asarray(wavelengths, dtype=float)
    refl = np.asarray(reflectances, dtype=float)
    out = np.full(refl.shape, np.nan, dtype=float)

    for start, stop in segments:
        if stop - start < 2:
            continue
        seg_wl = wl[start:stop]
        seg_refl = refl[start:stop]
        if not np.isfinite(seg_refl).any():
            continue
        continuum = calculate_continuum(seg_wl, seg_refl, method,
                                        polynomial_order)
        continuum = np.asarray(continuum, dtype=float)
        if normalization == 'subtraction':
            out[start:stop] = seg_refl - continuum
        else:
            with np.errstate(divide='ignore', invalid='ignore'):
                out[start:stop] = np.where(continuum > 0,
                                           seg_refl / continuum, np.nan)
    return out


def _remove_continuum_over_region(temp_maps, all_bands, processing_band_nums,
                                  method, polynomial_order, normalization,
                                  max_gap, keep_outside, rows, cols):
    """Replace each band slice with its continuum-removed values.

    Works in row blocks: a block holds every band for a few rows, which is
    the shape the continuum needs (one spectrum per pixel) while keeping the
    memory bounded. Each output band is streamed to its own binary file and
    read back with r.in.bin, so the full output cube is never held at once.

    :return: the list of output 2D map names, in band order.
    """
    from grass.script import array as garray

    wavelengths = np.array([b['wavelength'] for b in all_bands], dtype=float)
    processed = np.array([b['band_num'] in processing_band_nums
                          for b in all_bands], dtype=bool)
    n_bands = len(all_bands)
    if not processed.any():
        gs.fatal("No bands fall inside the requested wavelength range")

    segments = spectral_segments(wavelengths[processed], max_gap)
    gs.verbose(f"{len(segments)} contiguous coverage segment(s) over "
               f"{int(processed.sum())} bands")

    # About 400 MB of float64 per block, and at least one row.
    block_rows = max(1, min(rows, int(400e6 / max(1, n_bands * cols * 8))))

    tmpdir = gs.tempdir()
    handles = [open(os.path.join(tmpdir, f"b{i}.bin"), "wb")
               for i in range(n_bands)]
    region = gs.region()

    gs.percent(0, rows, 1)
    done = 0
    try:
        for start in range(0, rows, block_rows):
            stop = min(start + block_rows, rows)
            gs.use_temp_region()
            gs.run_command(
                'g.region',
                n=region['n'] - start * region['nsres'],
                s=region['n'] - stop * region['nsres'],
                w=region['w'], e=region['e'],
                nsres=region['nsres'], ewres=region['ewres'], quiet=True)
            try:
                cube = np.stack([np.asarray(garray.array(mapname=m),
                                            dtype=float)
                                 for m in temp_maps], axis=0)
            finally:
                gs.del_temp_region()

            nrow, ncol = cube.shape[1], cube.shape[2]
            flat = cube.reshape(n_bands, nrow * ncol)
            out = np.full_like(flat, np.nan)
            if keep_outside:
                out[~processed] = flat[~processed]
            sub_wl = wavelengths[processed]
            for pixel in range(flat.shape[1]):
                spectrum = flat[processed, pixel]
                if not np.isfinite(spectrum).any():
                    continue
                out[processed, pixel] = remove_continuum_spectrum(
                    sub_wl, spectrum, method, polynomial_order,
                    normalization, segments)

            out = out.reshape(n_bands, nrow, ncol)
            for i in range(n_bands):
                handles[i].write(np.ascontiguousarray(out[i],
                                                      dtype=np.float64).tobytes())
            done = stop
            gs.percent(done, rows, 1)
    finally:
        for handle in handles:
            handle.close()

    # A NaN in a DCELL raster is GRASS's own null, so the bands with no
    # measurement come back as null rather than as a plausible number.
    out_maps = []
    for i, band in enumerate(all_bands):
        name = f"tmp_cr_{band['band_num']}_{os.getpid()}"
        gs.run_command('r.in.bin', flags='f', input=os.path.join(tmpdir,
                       f"b{i}.bin"), output=name, bytes=8,
                       north=region['n'], south=region['s'],
                       east=region['e'], west=region['w'],
                       rows=rows, cols=cols, anull='nan',
                       overwrite=True, quiet=True)
        out_maps.append(name)

    gs.run_command('g.remove', type='raster', name=','.join(temp_maps),
                   flags='f', quiet=True)
    return out_maps


def apply_continuum_removal(input_raster, processing_bands, all_bands, method,
                            polynomial_order, output_raster, keep_outside,
                            normalization='division', max_gap=0.0):
    """Apply continuum removal to hyperspectral data"""
    
    gs.message(f"Applying continuum removal using {method} method...")
    
    # Get region info for processing
    region = gs.region()
    rows = int(region['rows'])
    cols = int(region['cols'])
    
    gs.message(f"Processing {rows} x {cols} pixels across {len(processing_bands)} bands...")
    
    # Split name3d into name and mapset components for Rast3d_extract_z_slice
    if '@' in input_raster:
        _name3d, _mapset3d = input_raster.split('@', 1)
    else:
        _name3d, _mapset3d = input_raster, ""

    # Create temporary files for each band
    temp_maps = []
    processing_band_nums = set(b['band_num'] for b in processing_bands)

    # Extract wavelengths for processing bands
    proc_wavelengths = np.array([b['wavelength'] for b in processing_bands])

    gs.percent(0, len(all_bands), 1)

    for idx, band in enumerate(all_bands):
        band_num = band['band_num']
        temp_output = f"tmp_continuum_{band_num}_{os.getpid()}"

        if band_num in processing_band_nums:
            gs.message(f"Processing band {band_num} ({band['wavelength']:.1f} nm)...")

        # extract_z_slice uses Rast3d_extract_z_slice (tile-bulk read, RASTER3D_NO_CACHE)
        # band_num is 1-based; z index is 0-based
        extract_z_slice(_name3d, _mapset3d, band_num - 1, temp_output)

        temp_maps.append(temp_output)
        gs.percent(idx + 1, len(all_bands), 1)

    # --- the continuum removal itself ------------------------------------
    # The band slices above are the input; this replaces each of them with
    # its continuum-removed values, in row blocks so that a large region
    # never needs the whole cube in memory. Bands outside the processing
    # range are left as they are when keep_outside is set, and nulled
    # otherwise, so the output never implies a measurement that was not made.
    temp_maps = _remove_continuum_over_region(
        temp_maps, all_bands, processing_band_nums, method, polynomial_order,
        normalization, max_gap, keep_outside, rows, cols)

    # Create 3D raster from temporary 2D rasters
    gs.message("Creating output 3D raster...")

    _band_wavelengths_nm = [b['wavelength'] for b in all_bands]

    if _RAS3D:
        import numpy as _np, ras3d as _r3, ras3d_write as _r3w
        from ras3d_grass_shim import get_band_cache
        _cube = _np.stack([get_band_cache()[m] for m in temp_maps], axis=0)
        _h = _r3.open_cube(input_raster)
        _r3w.write_raster3d(_r3w.outpath(output_raster), _cube, _h, wavelengths=_band_wavelengths_nm)
        _r3.close_cube(_h)
    else:
        # r.to.rast3 stacks the per-band 2D slices back into one 3D raster,
        # taking the depth order from the order of the list. r3.cross.rast is
        # the opposite operation, cutting a 2D section out of a 3D raster, so
        # it never wrote this output and the failure was swallowed by a bare
        # except that then reported success.
        gs.run_command('r.to.rast3', input=','.join(temp_maps),
                       output=output_raster, overwrite=True)
        if not gs.find_file(output_raster, element='grid3')['name']:
            gs.fatal(f"Failed to write 3D raster map <{output_raster}>")
    
    # Clean up temporary maps
    gs.message("Cleaning up temporary files...")
    for temp_map in temp_maps:
        try:
            gs.run_command('g.remove', type='raster', name=temp_map, 
                          flags='f', quiet=True)
        except:
            pass
    
    gs.message(f"Continuum removal complete: {output_raster}")


def print_band_info(all_bands, processing_bands, method):
    """Print information about bands and processing"""
    gs.message("=" * 70)
    gs.message("Continuum Removal Information:")
    gs.message("=" * 70)
    gs.message(f"Method: {method}")
    gs.message(f"Total bands in input: {len(all_bands)}")
    gs.message(f"Bands to process: {len(processing_bands)}")
    gs.message(f"{'Band':>6} | {'Wavelength (nm)':>16} | {'FWHM':>8} | {'Process':>10}")
    gs.message("-" * 70)
    
    processing_nums = set(b['band_num'] for b in processing_bands)
    
    for band in all_bands:
        fwhm_str = f"{band['fwhm']:.1f}" if band['fwhm'] else "N/A"
        process_str = "YES" if band['band_num'] in processing_nums else "NO"
        gs.message(f"{band['band_num']:>6} | {band['wavelength']:>16.2f} | "
                  f"{fwhm_str:>8} | {process_str:>10}")
    
    gs.message("=" * 70)
    if processing_bands:
        gs.message(f"Processing range: {processing_bands[0]['wavelength']:.1f} - "
                  f"{processing_bands[-1]['wavelength']:.1f} nm")
    gs.message("=" * 70)


def main(options, flags):
    """Main function"""
    input_raster = options['input']
    output_raster = options['output']
    min_wl = options.get('min_wavelength', None)
    max_wl = options.get('max_wavelength', None)
    method = options['method']
    polynomial_order = int(options['polynomial_order'])
    normalization = options['normalization']
    max_gap = float(options['max_gap'])
    only_valid = flags['n']
    keep_outside = flags['k']
    info_only = flags['i']
    
    # Validate polynomial order
    if method == 'polynomial' and (polynomial_order < 2 or polynomial_order > 5):
        gs.fatal("Polynomial order must be between 2 and 5")
    
    gs.message(f"Processing continuum removal for: {input_raster}")
    
    # Get all bands with wavelength metadata
    all_bands = get_all_band_wavelengths(input_raster, only_valid=only_valid)
    
    if not all_bands:
        gs.fatal("No valid bands found in input")
    
    # Determine wavelength range. The parser hands back an empty string for
    # an option the user left out, never None, so both of these have to test
    # for emptiness: float("") raises and the documented default of "use the
    # minimum available wavelength" was unreachable.
    if not min_wl:
        min_wl = all_bands[0]['wavelength']
        gs.message(f"Using minimum available wavelength: {min_wl:.1f} nm")
    else:
        min_wl = float(min_wl)
    
    if not max_wl:
        max_wl = all_bands[-1]['wavelength']
        gs.message(f"Using maximum available wavelength: {max_wl:.1f} nm")
    else:
        max_wl = float(max_wl)
    
    if min_wl >= max_wl:
        gs.fatal("Minimum wavelength must be less than maximum wavelength")
    
    # Filter bands by wavelength range
    all_bands_out, processing_bands = filter_bands_by_wavelength(
        all_bands, min_wl, max_wl
    )
    
    # Print information
    print_band_info(all_bands_out, processing_bands, method)
    
    if info_only:
        gs.message("Info mode: No output created.")
        return 0
    
    # Apply continuum removal
    apply_continuum_removal(
        input_raster, processing_bands, all_bands_out,
        method, polynomial_order, output_raster, keep_outside,
        normalization, max_gap
    )
    
    # Copy metadata to output
    gs.message("Setting output metadata...")
    try:
        gs.run_command('r3.support', map=output_raster,
                      title="Continuum Removed Hyperspectral Data",
                      history=f"Continuum removal applied using {method} method")
    except:
        pass
    
    gs.message(f"Continuum removal complete: {output_raster}")
    gs.message(f"Method: {method}")
    gs.message(f"Wavelength range: {min_wl:.1f} - {max_wl:.1f} nm")
    
    return 0


if __name__ == "__main__":
    options, flags = gs.parser()
    sys.exit(main(options, flags))
