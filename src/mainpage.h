/**
    \mainpage rspgen package

    \author  Yasushi Ikebe ikebe@milkyway.gsfc.nasa.gov
             James Peachey peachey@lheamail.gsfc.nasa.gov

    \section intro Introduction
    This package consists of a class library and an application. The
    library contains abstractions which facilitate combining information
    from 1) the spacecraft pointing file, 2) a spectrum, and 3)
    instrument response functions to produce an integrated, total instrument
    response in a format acceptable to Xspec. See
    <a href="http://glast.gsfc.nasa.gov/ssc/dev/binned_analysis/design_RspGen.html"> rspgen_design </a>
    for details of how the response is computed in each case. Currently, only
    responses for a stationary point source and gammaray bursts, for circular
    regions are implemented.

    \section parameters Application Parameters

    \subsection key Key To Parameter Descriptions
\verbatim
Automatic parameters:
par_name [ = value ] type

Hidden parameters:
(par_name = value ) type

Where "par_name" is the name of the parameter, "value" is the
default value, and "type" is the type of the parameter. The
type is enclosed in square brackets.

Examples:
infile [file]
    Describes an automatic (queried) file-type parameter with
    no default value.

(plot = yes) [bool]
    Describes a hidden bool-type parameter named plot, whose
    default value is yes (true).
\endverbatim

    \subsection general General Parameters
\verbatim
respalg [string]
    The type of response computation to perform, currently either
    GRB or PS. If this is set to GRB, the response will be computed
    for a single specific time supplied by the user, presumably the
    time of a gammaray burst. If this is set to PS, the GTI in the
    input spectrum file will be used with the spacecraft data file
    to pre-compute exposure as a function of spacecraft inclination
    angle with respect to the true point source direction. This
    differential exposure is then used to compute the effective total
    response to the given point source position.

specfile [file]
    Name of input spectrum file, in PHA1 format. Typically this file
    will have been created by the evtbin tool.

outfile [file]
    Name of output response file, in the OGIP "RSP" format.

scfile [file]
    Name of input spacecraft data file, FT2 format or equivalent.

ra [double]
    The true RA of the point source, in degrees.

dec [double]
    The true DEC of the point source, in degrees.

psfradius [double]
    The upper limit of integration for the point spread function. This
    is effectively the radius in degrees of the selected region, which
    presently is required to be a circle centered on the true point source
    location.

(resptype = DC1::Front [string])
    The specific set of functions used to compute the response. For
    details of what functions are available see the irfs package
    documentation. Note that the DC1 functions are pretty choppy, so
    for a better behaved response, try e.g. testIrfs::Front.

(resptpl = DEFAULT) [file]
    The full path to the template used to create the output file.
    This should not normally be changed. If DEFAULT is given, the
    rspgen application will use the standard template for LAT
    RSP files.
\endverbatim

    \subsection grb_pars GRB Response Parameters
\verbatim
time [double]
    The time of the burst.
\endverbatim

    \subsection ps_pars Point Source Response Parameters
\verbatim
thetacut [double]
    The maximum angle cutoff for the differential exposure
    computation. The differential exposure will be computed
    by binning the spacecraft pointing information into a
    histogram as a function of inclination angle covering the
    range [0, thetacut].

thetabinsize [double]
    The size of bins for the differential exposure
    computation.
\endverbatim

    \subsection energybins Energy Binning Parameters
\verbatim
energybinalg = LOG [string]
    Indicates how the energy bins will be specified. Legal values
    are FILE (bins will be read from a bin definition file), LIN
    (linearly uniform bins), and LOG (logarithmically uniform bins).
    This is only used if energy binning is required by the output
    type selected by the algorithm parameter.

(energyfield = ENERGY) [string]
    This is the name of the field containing the energy values for
    energy binning. The default value is consistent with the FT1
    format.

emin [double]
    The lowest energy of the first interval for linearly or
    logarithmically uniform bins. Only used if energybinalg
    is LIN or LOG.

emax [double]
    The highest energy of the last interval for linearly or
    logarithmically uniform bins. Only used if energybinalg
    is LIN or LOG.

enumbins [integer]
    The number of bins for logarithmically uniform bins. Only
    used if energybinalg is LOG.

deltaenergy [double]
    The width of linearly uniform bins. Only used if energybinalg
    is LIN.

energybinfile [file]
    The name of the energy bin definition file. Only used if
    energybinalg is FILE.
\endverbatim

*/
