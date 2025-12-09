## DESCRIPTION

*i.hyper.continuum* applies continuum removal to hyperspectral imagery
imported as 3D raster maps (`raster_3d`) by
[i.hyper.import](i.hyper.import.html).

The module performs continuum removal, a spectral processing technique
that normalizes reflectance spectra to emphasize absorption features.
Continuum removal divides the original spectrum by a continuum line (or
envelope) that connects the spectral peaks, effectively removing the
overall spectral shape and highlighting absorption bands.

*i.hyper.continuum* is part of the **i.hyper** module family designed for
hyperspectral data import, processing, and analysis in GRASS. It is
essential for mineral identification, vegetation analysis, and any
application requiring detailed analysis of spectral absorption features.

The module allows users to specify a wavelength range for processing,
producing a spectral subset of the original hyperspectral image. By
default, it processes the entire spectral range of the input data. Bands
outside the specified range can optionally be kept unchanged in the
output.

Three methods are available for calculating the continuum:

- **convex_hull** -- Uses convex hull algorithm to find the upper
  envelope of the spectrum (default, recommended for most applications)
- **linear** -- Simple linear continuum connecting the first and last
  points in the wavelength range
- **polynomial** -- Polynomial fitting to the spectral envelope (order
  2-5)

The convex hull method is most widely used in remote sensing as it
effectively captures the spectral continuum without making assumptions
about the shape. It identifies the upper envelope by constructing the
convex hull of the spectral points and extracting the upper boundary.

Continuum-removed reflectance values typically range from 0 to 1, where
1 represents points on the continuum (spectral peaks) and values less
than 1 represent absorption features. Deeper absorptions produce values
closer to 0.

## NOTES

The module expects input data to be a 3D raster map created by
*i.hyper.import* or any 3D raster with wavelength metadata stored in
band-level metadata following the *i.hyper* standard format:
**wavelength**, **FWHM**, **valid**, and **unit**.

When wavelength metadata is not found, the module will fail with an
error message. It is essential to use properly imported hyperspectral
data with wavelength information.

The wavelength range parameters (*min_wavelength* and *max_wavelength*)
define both the spectral region for continuum calculation and the extent
of the output 3D raster. If not specified, the module uses the minimum
and maximum wavelengths available in the input data.

The **-n** flag restricts processing to only bands marked as valid
(`valid=1`) in the metadata. This is useful for excluding bands with
quality issues, such as those affected by atmospheric absorption (water
vapor bands around 1400 nm and 1900 nm) or sensor artifacts.

The **-k** flag preserves original reflectance values for bands outside
the specified wavelength range. Without this flag, all bands in the
output maintain the structure of the input, but only bands within the
range undergo continuum removal. This is useful when you want to analyze
a specific absorption feature while maintaining the full spectral cube
structure.

The **-i** flag provides information mode that displays band information
and processing parameters without creating the output raster. This helps
verify which bands will be processed before running the full analysis.

For the polynomial method, the order parameter controls the degree of
the polynomial fit. Lower orders (2-3) produce smoother continua, while
higher orders (4-5) can follow more complex spectral shapes but may
overfit noise. The convex hull method is generally preferred as it
adapts naturally to the spectral shape without requiring parameter
tuning.

Continuum removal is particularly effective for:

- **Mineral identification** : Absorption features at specific
  wavelengths indicate different minerals (e.g., clay minerals at
  2200 nm, iron oxides at 900 nm)
- **Vegetation analysis** : Chlorophyll absorption around 680 nm, water
  content features at 1450 nm and 1950 nm
- **Comparison across scenes** : Normalizing for illumination and
  atmospheric effects
- **Feature depth analysis** : Quantifying absorption strength for
  compositional mapping

The output 3D raster maintains the same spatial resolution and extent as
the input, with continuum-removed reflectance values replacing the
original reflectance in the processed wavelength range.

## EXAMPLES

::: code

    # Apply continuum removal to entire spectral range using convex hull
    i.hyper.continuum input=prisma \
                      output=prisma_cr

    # The output contains continuum-removed reflectance across all wavelengths
:::

::: code

    # Focus on clay mineral absorption feature (2100-2400 nm)
    i.hyper.continuum input=enmap \
                      output=enmap_clay_cr \
                      min_wavelength=2100 \
                      max_wavelength=2400

    # Extract specific band for visualization
    # Band depth at 2200 nm indicates clay mineral presence
    r3.cross.rast input=enmap_clay_cr output=clay_2200nm@2200
:::

::: code

    # Apply linear continuum removal for faster processing
    i.hyper.continuum input=tanager \
                      output=tanager_cr_linear \
                      method=linear \
                      min_wavelength=400 \
                      max_wavelength=2500
:::

::: code

    # Polynomial continuum with 3rd order fit
    i.hyper.continuum input=prisma \
                      output=prisma_cr_poly3 \
                      method=polynomial \
                      polynomial_order=3 \
                      min_wavelength=400 \
                      max_wavelength=2500
:::

::: code

    # Process only iron oxide absorption region (800-1000 nm)
    # Keep other wavelengths unchanged
    i.hyper.continuum input=enmap \
                      output=enmap_fe_cr \
                      min_wavelength=800 \
                      max_wavelength=1000 \
                      -k

    # The -k flag preserves original values outside 800-1000 nm
:::

::: code

    # Preview processing parameters without creating output
    i.hyper.continuum input=prisma \
                      output=test \
                      min_wavelength=2000 \
                      max_wavelength=2400 \
                      -i

    # Info mode displays which bands will be processed
:::

::: code

    # Use only valid bands for processing
    i.hyper.continuum input=enmap \
                      output=enmap_cr_valid \
                      min_wavelength=400 \
                      max_wavelength=2500 \
                      -n

    # Excludes bands with valid=0 (atmospheric absorption, sensor issues)
:::

::: code

    # Complete workflow: mineral mapping with continuum removal
    # Step 1: Import hyperspectral data
    i.hyper.import input=/data/PRISMA.he5 \
                   product=prisma \
                   output=prisma

    # Step 2: Apply continuum removal to Al-OH region (clay minerals)
    i.hyper.continuum input=prisma \
                      output=prisma_aloh_cr \
                      min_wavelength=2000 \
                      max_wavelength=2400

    # Step 3: Extract band at 2200 nm (main clay absorption)
    r3.cross.rast input=prisma_aloh_cr output=clay_depth@2200

    # Step 4: Calculate absorption depth (1 - continuum_removed_reflectance)
    r.mapcalc "absorption_depth = 1.0 - clay_depth"

    # Step 5: Classify clay mineral abundance
    r.recode input=absorption_depth output=clay_class << EOF
    0.0:0.1:1:low
    0.1:0.2:2:moderate
    0.2:1.0:3:high
    EOF

    # Step 6: Visualize results
    r.colors map=absorption_depth color=viridis
    d.rast map=absorption_depth
:::

::: code

    # Vegetation water content analysis
    # Step 1: Apply continuum removal to water absorption region
    i.hyper.continuum input=enmap \
                      output=enmap_water_cr \
                      min_wavelength=1400 \
                      max_wavelength=2000

    # Step 2: Extract bands at water absorption centers
    r3.cross.rast input=enmap_water_cr output=water_1450@1450
    r3.cross.rast input=enmap_water_cr output=water_1950@1950

    # Step 3: Calculate water absorption index
    r.mapcalc "water_index = (1.0 - water_1450) + (1.0 - water_1950)"

    # Higher values indicate stronger water absorption (higher water content)
:::

::: code

    # Compare different continuum methods
    # Convex hull (default)
    i.hyper.continuum input=prisma output=prisma_cr_hull \
                      min_wavelength=2000 max_wavelength=2400

    # Linear
    i.hyper.continuum input=prisma output=prisma_cr_linear \
                      method=linear \
                      min_wavelength=2000 max_wavelength=2400

    # Polynomial order 2
    i.hyper.continuum input=prisma output=prisma_cr_poly2 \
                      method=polynomial polynomial_order=2 \
                      min_wavelength=2000 max_wavelength=2400

    # Extract same band from each for comparison
    r3.cross.rast input=prisma_cr_hull output=hull_2200@2200
    r3.cross.rast input=prisma_cr_linear output=linear_2200@2200
    r3.cross.rast input=prisma_cr_poly2 output=poly2_2200@2200
:::

## SEE ALSO

[i.hyper.import](i.hyper.import.html),
[i.hyper.rgb](i.hyper.rgb.html),
[i.hyper.albedo](i.hyper.albedo.html),
[i.hyper.composite](i.hyper.composite.html),
[i.hyper.preproc](i.hyper.preproc.html),
[i.hyper.explore](i.hyper.explore.html),
[i.hyper.export](i.hyper.export.html),
[r3.cross.rast](https://grass.osgeo.org/grass-stable/manuals/r3.cross.rast.html),
[r.mapcalc](https://grass.osgeo.org/grass-stable/manuals/r.mapcalc.html),
[r3.support](https://grass.osgeo.org/grass-stable/manuals/r3.support.html)

## REFERENCES

- Clark, R. N., & Roush, T. L. (1984). Reflectance spectroscopy:
  Quantitative analysis techniques for remote sensing applications.
  *Journal of Geophysical Research: Solid Earth*, 89(B7), 6329-6340.
- Kokaly, R. F., & Clark, R. N. (1999). Spectroscopic determination of
  leaf biochemistry using band-depth analysis of absorption features and
  stepwise multiple linear regression. *Remote Sensing of Environment*,
  67(3), 267-287.
- van der Meer, F. (2004). Analysis of spectral absorption features in
  hyperspectral imagery. *International Journal of Applied Earth
  Observation and Geoinformation*, 5(1), 55-68.
- Kruse, F. A., et al. (1993). The spectral image processing system
  (SIPS)—interactive visualization and analysis of imaging spectrometer
  data. *Remote Sensing of Environment*, 44(2-3), 145-163.

## AUTHORS

Created for the i.hyper module family

Based on work by Alen Mangafić and Tomaž Žagar, Geodetic Institute of
Slovenia
