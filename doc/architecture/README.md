---
package: rspgen
summary: gtrspgen and its library, which build an OGIP response matrix (RSP) for a LAT PHA1 spectrum from the IRFs, the spacecraft pointing history and the spectrum's DSS cuts, for a steady point source or a GRB.
layer: analysis
depends_on: [evtbin, irfs, dataSubselector, tip, astro, st_app, st_stream, facilities, CLHEP]
used_by: [sane]
last_reviewed_commit: 792b3f9
---

# rspgen Architecture

**rspgen** builds the response file that goes with a PHA1 spectrum made by
`gtbin` (evtbin), so that the pair can be fitted in XSPEC. You give it the
spectrum and an FT2 file. It reads the source position and the extraction
radius from the spectrum's DSS keywords and computes, for each true-energy
bin, the effective area, energy redistribution and PSF containment. For a
steady source (`PS`) these are averaged over the livetime distribution of
the source's off-axis angle. For a burst (`GRB`) they are evaluated at the
single pointing at the burst time.

The code was written at the HEASARC (Y. Ikebe, J. Peachey) in 2004–2009.
The last functional change was the `irfs=CALDB` support added in 2013.
`src/mainpage.h` still documents `ra`, `dec` and `psfradius` parameters
that no longer exist (they now come from the DSS keywords) and links to a
design page that is no longer online. Use `pfiles/gtrspgen.par` and
`src/RspGenApp.cxx` as the source of truth. User help is in
`src/fermitools-fhelp/fhelp_files/gtrspgen.txt`.

## Purpose and scope

**Does:**
- Read the apparent-energy channels from the spectrum's `EBOUNDS`, the
  GTI from its `GTI` extension, and the cone selection and IRF selection
  from the DSS keywords of its `SPECTRUM` HDU.
- Bin the FT2 livetime into (cos θ, φ) of the source in the instrument
  frame (`PS`), or interpolate the pointing at one time (`GRB`).
- Compute `R(E_true, channel) = ⟨Aeff · ∫_channel edisp · PSF(<r)⟩` in
  cm², summed over the IRF components (for example front + back), and
  write it as a compressed OGIP `RSP_MATRIX` with an `EBOUNDS`
  extension.
- Write `RESPFILE = <outfile>` into the input spectrum's `SPECTRUM`
  header. If the spectrum file is read-only, this is skipped silently.

**Does not:**
- Bin events or compute the spectrum's `EXPOSURE`. `gtbin` does that (see
  `src/evtbin/doc/architecture/`).
- Handle extended sources, non-circular regions, more than one region,
  or a source that is not at the region centre.
- Use any DSS cut other than the sky cone (and the event class and type
  through `irfs=CALDB`). Energy and zenith-angle cuts are not modelled.
  Zenith or ROI cuts are expected to be reflected in the GTI.
- Handle GBM data. Only the LAT `BinConfig` makes sense here.
- Produce an ARF and RMF separately.

## Position in the system

| Direction | Package | Why |
|-----------|---------|-----|
| Depends on | `irfs` (`irfInterface`, `irfLoader`) | `IrfsFactory::create`, `IAeff::value`, `IEdisp::integral`, `IPsf::angularIntegral`, and `Loader::go`/`respIds` for name lookup. |
| Depends on | `evtbin` | `Binner`/`OrderedBinner`/`LinearBinner`, `Hist1D`/`Hist2D` (differential exposure), `Gti` (GTI read and `getFraction`), and `BinConfig` (true-energy parameters: prompting and binner creation). |
| Depends on | `dataSubselector` | `Cuts(spec_file, "SPECTRUM", …)` reads the DSS keywords. `SkyConeCut` supplies the source position and radius, and `Cuts::CALDB_implied_irfs()` resolves `irfs=CALDB`. |
| Depends on | `tip` | All FITS I/O, plus `LinearInterp` for FT2 pointing interpolation. |
| Depends on | `astro` | `SkyDir` angles. CLHEP `Hep3Vector` is used in `calcPhi`. |
| Depends on | `st_app`, `st_stream` | `RspGenApp` is an `st_app::StApp` compiled **into the library**. It gets `st_app` transitively through evtbin's public link, because rspgen only links it for the `gtrspgen` executable. |
| Depends on | `facilities` | `commonUtilities::getDataPath("rspgen")` (default template) and `Util::basename`. It comes in transitively and is not linked directly. |
| Link only | `st_facilities` | Linked `PUBLIC`. The only use is an unused `#include "st_facilities/Env.h"` in `RspGenApp.cxx`. |
| Used by | `sane` | `gt_apps.rspgen = GtApp('gtrspgen', 'rspgen')` wraps the tool for Python. `test_sane` links the library but includes no rspgen header. |

No other package includes `rspgen/` headers. `Likelihood/src/FileUtils.cxx`
writes DRM tables "to look like the output of gtrspgen", using the same
column layout but none of the code. Changing the RSP layout here does not
change that code, and the reverse is also true.

## Key abstractions

```mermaid
classDiagram
    class IResponse {
      <<abstract>>
      +IResponse(resp_type, spec_file, true_en_binner)
      +writeOutput(creator, file, template)
      +compute(E_true, response&)*
      #m_true_en_binner : Binner*
      #m_app_en_binner : OrderedBinner (from EBOUNDS)
      #m_irfs : vector~Irfs*~
      #calcPhi(x, z, dir)
    }
    IResponse <|-- PointResponse
    IResponse <|-- GrbResponse
    PointResponse *-- Hist2D : Δt(θ,φ)
    PointResponse *-- IWindow
    GrbResponse *-- IWindow
    IWindow <|-- CircularWindow
    class CosineBinner {
      evtbin::Binner, uniform in cos θ, value in degrees
    }
    evtbin_Binner <|-- CosineBinner
    Hist2D o-- CosineBinner
    class SpaceCraftCalculator {
      +lookUpResponse(name, match)$
      +displayIrfNames(os)$
      (θ-only duplicate of PointResponse; tests only)
    }
    class RspGenApp {
      st_app::StApp, private header src/RspGenApp.h
    }
    RspGenApp ..> IResponse : creates
    IResponse ..> SpaceCraftCalculator : name lookup
```

| Type | Header | Role |
|------|--------|------|
| `IResponse` | `rspgen/IResponse.h` | Base class. Its constructor reads `EBOUNDS` (keV → MeV) into an `OrderedBinner` of apparent-energy channels, clones the true-energy binner, and creates one `Irfs` per component name. `writeOutput` loops over true-energy bins, calls `compute(midpoint)`, compresses each row and writes `MATRIX` and `EBOUNDS`. |
| `PointResponse` | `rspgen/PointResponse.h` | Steady source. Its constructor builds `Hist2D` Δt(θ, φ) from FT2 rows prorated by the spectrum GTI. `compute` sums over the (θ, φ) bin centres. |
| `GrbResponse` | `rspgen/GrbResponse.h` | Burst. A single (θ, φ) interpolated from FT2 at `time`. A second constructor takes θ, φ, binners, one `Irfs` and a window directly. Only the tests use it. |
| `IWindow`, `CircularWindow` | `rspgen/IWindow.h`, `CircularWindow.h` | Extraction region. `integrate` returns `IPsf::angularIntegral(E, θ, φ, radius)`, which is the PSF containment fraction within the radius. |
| `rspgen::CosineBinner` | `rspgen/CosineBinner.h` | An `evtbin::Binner` whose `computeIndex(θ°)` bins **cos θ** uniformly. `getInterval` returns θ in degrees. Not to be confused with `healpix::CosineBinner`. |
| `SpaceCraftCalculator` | `rspgen/SpaceCraftCalculator.h` | Its static `lookUpResponse`/`displayIrfNames` are used by `IResponse`. The constructor and `psf()` are a θ-only, **linear-θ** copy of `PointResponse` that only `test10` uses. The TODOs in `PointResponse.cxx` say the two were meant to be merged. |
| `RspGenApp` | `src/RspGenApp.h` (private) | The `gtrspgen` application. `src/gtrspgen/gtrspgen.cxx` is only the `StAppFactory` line. |

## Data flow and control flow

```
gtrspgen → RspGenApp::run
  irfLoader::Loader::go(); evtbin::BinConfig::load()
  prompt: respalg, specfile, scfile, sctable, outfile,
          GRB: time | PS: thetacut, dcostheta, phinumbins,
          irfs, resptpl, then BinConfig::energyParPrompt (ebinalg, efield, emin, emax, denergy, enumbins, ebinfile)
  writeResponse
    true_en_binner = BinConfig::create(specfile)->createEnergyBinner(pars)    // MeV; ebinfile ENERGYBINS keV→MeV
    resptpl == DEFAULT → $FERMI_INST_DIR/refdata/fermi/rspgen/LatResponseTemplate
    editTable(specfile, SPECTRUM)["RESPFILE"] = outfile        // TipException swallowed
    Cuts(specfile, "SPECTRUM", check_columns=false, skipTimeRangeCuts=true)
      exactly one SkyConeCut → (ra, dec, radius) else throw
    irfs == CALDB → cuts.CALDB_implied_irfs()
    PS  → PointResponse(ra, dec, thetacut, dcostheta, radius, phinumbins, irfs, specfile, scfile, sctable, binner)
    GRB → GrbResponse(ra, dec, time, radius, irfs, specfile, scfile, sctable, binner)
    response->writeOutput("gtrspgen", outfile, template)
```

### Integration over livetime and angles (`PS`)

The `PointResponse` constructor does the following:

```
gti = evtbin::Gti(specfile)                     // the spectrum's GTI extension
for each FT2 row while GTIs remain:
    f  = gti.getFraction(START, STOP, pos)      // fraction of [START,STOP) inside the GTI
    Δt = f × LIVETIME
    pointing = LinearInterp on START at (START+STOP)/2   // SCZ, SCX
    θ = angle(source, SCZ) in degrees;  φ = calcPhi(SCX, SCZ, source) in [0, 360)
    hist[θ-bin][φ-bin] += Δt                    // θ outside [0, thetacut) is not binned …
    T += Δt                                     // … but still counted in T
```

`compute(E)` then evaluates, at each θ and φ bin centre and for each IRF
component `k`:

```
R[j] = Σ_{θ,φ} Σ_k  hist[θ][φ]/T · Aeff_k(E,θ,φ) · ∫_{E'_j,lo}^{E'_j,hi} edisp_k(E'|E,θ,φ) dE' · PSF_k(<r | E,θ,φ)
```

`T` is the total GTI-prorated livetime, which is the same quantity
`gtbin` writes as the PHA `EXPOSURE`. So `R` is a livetime-averaged
effective area (cm²), and `counts = flux × EXPOSURE × R` holds. Livetime
when the source is beyond `thetacut` counts as zero response rather than
being dropped. With the default `thetacut=60` and `dcostheta=0.05` there
are 10 θ bins, from cos θ = 0.5 to 1.

`GRB` uses the same formula with a single (θ, φ) and weight 1, so `R` is
the instantaneous effective area at `time`.

### Energy conventions

| Quantity | Units inside | Units on disk |
|----------|--------------|---------------|
| Apparent channels (input `EBOUNDS` `E_MIN`/`E_MAX`) | MeV (×0.001 on read) | keV |
| True-energy binning parameters `emin`/`emax`/`denergy` | MeV | n/a |
| `ebinfile` (`ENERGYBINS`, made by `gtbindef`) | MeV (evtbin's `BinConfig` ×0.001) | keV |
| Output `MATRIX` `ENERG_LO`/`ENERG_HI`, output `EBOUNDS` | MeV | keV (×1000 on write) |
| IRF calls (`Aeff`, `Edisp`, `Psf`) | MeV, θ and φ in degrees | n/a |

This matches evtbin: MeV in code, keV in every EBOUNDS and ENERGYBINS
HDU.

### Output compression

Each `MATRIX` row stores only the channels whose response is greater than
`LO_THRES = 0.0`. Consecutive positive channels form groups: `N_GRP`,
`F_CHAN` (1-based) and `N_CHAN` hold the groups, and `MATRIX` holds only
the stored values. `F_CHAN` and `N_CHAN` are padded with zeros to at
least 3 entries. A row where every channel is zero has `N_GRP = 0` and an
empty `MATRIX`. Compression is controlled by a global
`bool rspgen::compressed = true` in `IResponse.cxx`.

## Public interfaces

### Application

| Tool | Par file | Notes |
|------|----------|-------|
| `gtrspgen` | `pfiles/gtrspgen.par` | `respalg` ∈ {`GRB`, `PS`}. `irfs` defaults to **`P6_V3_DIFFUSE`** (see gotchas). `CALDB` means "use the IRFs implied by the DSS keywords". Hidden: `sctable=SC_DATA`, `resptpl=DEFAULT`, `phinumbins=1` (prompted only when `respalg=PS`). Energy parameters are the same as in `gtbin.par`. |

For the startup, error contract and standard hidden parameters, see
`src/st_app/doc/architecture/`. For par files and `PFILES`, see
`src/ape/doc/architecture/parameter-files.md`.

### C++ library

Headers are in `rspgen/`. Nothing outside the package includes them. The
class most likely to be reused is `IResponse` with a custom `compute`.

### Python

None in this package. `sane`'s `gt_apps.rspgen` wraps the executable.

### File formats

- **Input spectrum:** PHA1 from `gtbin`, with `SPECTRUM` (DSS keywords,
  `RESPFILE`), `EBOUNDS` (`CHANNEL` 1-based, `DETCHANS`) and `GTI`. The
  primary or first extension with `TELESCOP`/`INSTRUME` selects the
  `BinConfig`.
- **Input FT2:** `START`, `STOP`, `LIVETIME`, `RA_SCZ`, `DEC_SCZ`,
  `RA_SCX` and `DEC_SCX`, sorted by `START`.
- **Output RSP:** created from `data/LatResponseTemplate`. It has a
  primary HDU, `MATRIX` (`HDUCLAS2=RSP_MATRIX`, `ENERG_LO`/`ENERG_HI` E,
  `N_GRP` I, and variable-length `F_CHAN` PI, `N_CHAN` PI and `MATRIX` PE)
  and `EBOUNDS`. `TELESCOP`, `INSTRUME`, `DATE-OBS` and `DATE-END` are
  copied from the input `EBOUNDS`. `DATE`, `CREATOR = "gtrspgen"` (no
  version) and `FILENAME` are set. `DETCHANS` is the number of input
  channels. **DSS keywords are not copied** (`NDSKEYS = 0`).

### Environment

`FERMI_INST_DIR` is used, through `getDataPath`, for the default
template. `CALDB` is used through `irfLoader` (see "Runtime environment"
in `docs/architecture/README.md`). The package reads no environment
variable itself.

## Build and test

- **CMake target:** `rspgen` is a **STATIC** library of 7 sources (about
  1,000 lines). It contains `RspGenApp`. Executables: `gtrspgen` and
  `test_rspgen`. The `install(TARGETS …)` line has `EXPORT fermiTargets`
  commented out, unlike most packages. `include(UseSWIG)` is present but
  nothing uses SWIG.
- **Install:**
  - headers go to `include/rspgen`
  - `pfiles/` goes to `syspfiles`
  - `data/` (template and test inputs, excluding `outref`) goes to
    `${FERMI_INSTALL_REFDATADIR}/rspgen`
  - `data/outref/` goes to `${FERMI_INSTALL_DATADIR}/rspgen/outref`
- **`test_rspgen`:**
  - It is an `st_app` application registered under the name
    **`gtrspgen`**, so it uses `gtrspgen.par` (there is no
    `test_rspgen.par`). `PFILES` must find that file. It is **not
    registered with CTest**.
  - It finds inputs with `getDataPath("rspgen")`, which means **the
    install tree**: `$FERMI_INST_DIR/refdata/fermi/rspgen`. The inputs
    are `PHA1.pha`, `ft2tiny.fits` (table `Ext1`) and
    `rspgen_energy_bins.fits`.
  - It uses the `testIrfs::Front`/`testIrfs::Back` and `TEST`/`TESTF`/`TESTB`
    responses from `irfs/testResponse`, so no CALDB IRFs are needed.
  - Tests 1–3 compute responses "by hand" with `tip`. Tests 4–8 use
    `GrbResponse`, `PointResponse` and `RspGenApp`. Test 9 covers DSS
    parsing, test 10 `SpaceCraftCalculator`, and test 11 checks that
    front + back equals the total to within 10 ε. `testPhiCalc` checks
    `calcPhi` in all octants. **`test12` is skipped** ("Pass 7 irfs were
    removed"), so `phinumbins > 1` with real IRFs is untested.
  - Every check prints `Unexpected: …` and sets a failure flag, and the
    test throws at the end if any check failed.
  - Outputs (`test_response1..8.rsp` and a copied, modified `PHA1.pha`)
    go to the **current directory**.
  - Test 6 makes its `PHA1.pha` copy read-only and expects `RESPFILE` to
    stay `NONE`. **That check fails when the test runs as root**, because
    root ignores mode 0444 (inferred).
- **outref consumers:** nothing in this repo reads `data/outref/`.
  `recipe/tests/ST-unit-test` (`check_test_rspgen_output`) compares
  `PHA1.pha` against
  `$CONDA_PREFIX/share/fermitools/data/rspgen/outref/PHA1.pha` with
  `ftdiff`, but checks only that `test_response1..8.rsp` **exist**. The
  `.rsp` files in `outref/` are not compared, and `test_response_p7.rsp`
  is left over from the skipped test12.
- **`sane`:** `python/tests/rspgen_tests.py` runs `gtrspgen` with
  `irfs='DC1AF'` on `Crab.pha`. Both are long obsolete.
- **Legacy:** `SConscript`, `rspgenLib.py` and `azure-pipelines.yml`.
- **Branches:** `origin/addition_creator_fits` adds a versioned `CREATOR`
  through `GenericUtils::creator_banner`. `origin/f26-cpp17-compliance`
  has small C++17 fixes. `origin/Issue17` makes the library type
  optional. None of them is merged.

## Invariants and gotchas

**Inputs and side effects**
- **gtrspgen modifies its input spectrum** by setting `RESPFILE`. If the
  file is not writable, the failure is swallowed with no message (there
  is a `TODO: add warning`).
- The spectrum must have **exactly one** cone (`SkyConeCut`) in its DSS
  keywords. The cone centre is assumed to be the true source position,
  and its radius is the PSF integration radius. Running `gtselect`
  without a cone, or twice with two cones, makes gtrspgen throw.
- `EBOUNDS` channels must be 1-based. A channel ≤ 0 throws. Rows with
  `CHANNEL > DETCHANS` are skipped. If `EBOUNDS` has fewer rows than
  `DETCHANS`, the missing channels become zero-width `(0, 0)` intervals
  and get zero response.
- The `clobber` parameter is ignored. `tip::createFile` always
  overwrites.

**Defaults that no longer work**
- The `irfs` default `P6_V3_DIFFUSE` is a Pass 6 name. Nothing in the
  current `irfs` package registers it, so the default will probably fail
  with "cannot find response function(s)" and a list of valid names.
  Use `irfs=CALDB`, as the fhelp examples do (inferred from a grep; not
  checked at runtime).
- The `phinumbins` prompt text says "Number of logarithmically uniform
  energy bins". It is the number of φ bins.

**Numerics**
- The response for each true-energy bin is evaluated **at the bin's
  midpoint only**. Nothing is averaged within the bin, so use fine
  true-energy binning where Aeff varies quickly (low energies).
- `PS` evaluates the IRFs at **(θ, φ) bin centres**. With the default
  `phinumbins=1`, φ is fixed at **180°**: the response is not averaged
  over φ, it uses the value at φ = 180°.
- `rspgen::CosineBinner`: θ = 0 exactly (cos θ = 1) falls outside
  `[cos thetacut, 1)` and is dropped. Bin 0 is the **largest** θ. The
  last bin can be narrower when `(1 − cos thetacut)/dcostheta` is not an
  integer. `dcostheta` is a cos θ step, not degrees.
- Livetime outside `thetacut` lowers the response (it stays in `T`).
  That is correct when it is paired with gtbin's `EXPOSURE`, and wrong
  if someone renormalizes `T`.
- `SpaceCraftCalculator` bins θ **linearly**, unlike `PointResponse`.
  Its `psf()`, and `PointResponse::psf()`, index the histogram with
  `computeIndex(θ)` without checking for −1, which reads out of bounds
  when θ is outside the range (from reading the code; only tests call
  them).

**FT2 coverage**
- Pointing is interpolated with `tip::LinearInterp` on `START` at each
  row's midpoint (`PS`) or at `time` (`GRB`). This throws if the time is
  before the first `START` or at or after the **last** `START`. With
  `PS`, if the last FT2 row overlaps the spectrum's GTI, the run fails.
  The FT2 file must extend at least one row past the end of the GTI
  (from reading the code and `tip/src/LinearInterp.cxx`; not checked at
  runtime).
- If no livetime overlaps the GTI, `PointResponse` throws "cannot
  continue with 0. total exposure".

**Ownership and globals**
- `IResponse` owns its binners and `Irfs`, and deletes them in its
  destructor. The test-only `GrbResponse(θ, φ, …)` constructor clones
  everything it is passed.
- `rspgen::compressed` is a mutable global with external linkage. Like
  evtbin, rspgen is static, so each shared library that links it gets
  its own copy.
- `setVersion("$Name:  $")` is a leftover CVS keyword, so the tool's
  reported version is meaningless.

## Extension points

- **New region shape:** subclass `IWindow` and implement
  `integrate(psf, E, θ, φ)` (for example an annulus, as the difference
  of two `angularIntegral` calls). Then select it in `RspGenApp` from a
  different DSS cut type. Today only `SkyConeCut` is recognized.
- **New response mode** (for example time-resolved, or for a moving
  source): subclass `IResponse`, implement `compute(E_true, response)`
  by filling `response` with `m_app_en_binner->getNumBins()` values in
  cm², add a `respalg` value in `RspGenApp::prompt`/`writeResponse`, and
  add it to the par file.
- **Averaging within true-energy bins:** change `IResponse::writeOutput`,
  which currently calls `compute(true_en.midpoint(), …)`.
- **Merging `SpaceCraftCalculator` and `PointResponse`** is an
  outstanding TODO. If you do it, keep `PointResponse`'s cos θ and φ
  binning and its normalization to the total livetime.

Before merging, run `test_rspgen` from a scratch directory with an
installed tree. Also run the `sane` rspgen test if `sane` is in use.
Changes to evtbin's `Gti`, `Binner` or `BinConfig`, and to
`dataSubselector::Cuts`, can change rspgen output.

## Further reading

- `src/evtbin/doc/architecture/`: PHA1 production, `EXPOSURE`, `Gti` and
  binners.
- `src/dataSubselector/` (`Cuts`, `SkyConeCut`, `CALDB_implied_irfs`):
  how DSS keywords become cuts.
- `src/irfs/`: `IAeff`, `IEdisp::integral`, `IPsf::angularIntegral` and
  `testResponse`.
- `src/fermitools-fhelp/fhelp_files/gtrspgen.txt`: user documentation.
- OGIP CAL/GEN/92-002 (response matrix format) and OGIP/92-007 (PHA).
- `src/mainpage.h` and `doc/release.notes`: history only, partly stale.
