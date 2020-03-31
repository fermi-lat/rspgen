/** \file test_rspgen.cxx
    \brief Test application for response generator code.
    \author Yasushi Ikebe, James Peachey
*/
// C++ standard inclusions.
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <utility>
#include <vector>

#ifndef WIN32
#include <sys/types.h>
#include <sys/stat.h>
#endif

// Main appliation class RspGenApp.
#include "RspGenApp.h"

// Sky directions use astro::SkyDir.
#include "astro/SkyDir.h"

// CLHEP for operations not supported directly by astro.
#include "CLHEP/Vector/ThreeVector.h"

// Include irfs.
#include "irfLoader/Loader.h"

// Handle cuts.
#include "dataSubselector/Cuts.h"
#include "dataSubselector/SkyConeCut.h"

// From evtbin, generic histograms and associated binners are used.
#include "evtbin/Hist1D.h"
#include "evtbin/Binner.h"
#include "evtbin/LinearBinner.h"
#include "evtbin/OrderedBinner.h"

// Use instrument response function abstractions from irfInterface.
#include "irfInterface/Irfs.h"
#include "irfInterface/IrfsFactory.h"

// Tests include circular region (window) abstractions.
#include "rspgen/CircularWindow.h"

// Response abstraction for burst case.
#include "rspgen/GrbResponse.h"

// Response abstraction for burst case.
#include "rspgen/PointResponse.h"

// Spacecraft position calculator.
#include "rspgen/SpaceCraftCalculator.h"

// Standard application-related code.
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

// Use st_facilities to find data files.
#include "st_facilities/Env.h"
#include "facilities/commonUtilities.h"

#include "st_stream/StreamFormatter.h"
#include "st_stream/Stream.h"
#include "st_stream/st_stream.h"

// Table access through tip.
#include "tip/Header.h"
#include "tip/IFileSvc.h"
#include "tip/LinearInterp.h"
#include "tip/Table.h"

// Obvious unit conversions.
static double s_keV_per_MeV = 1000.;
static double s_MeV_per_keV = .001;

static std::string s_cvs_id = "$Name:  $";

/** \class RspGenTestApp
    \brief The main test application itself.
*/
class RspGenTestApp : public st_app::StApp {
  public:
    /// \brief Construct a test application.
    RspGenTestApp();

    virtual ~RspGenTestApp() throw() {}

    /// \brief Run all tests, reporting failures.
    void run();

    /** \brief Simple test in which apparent energy bins are read from ebounds, but true energy bins are completely
        fabricated. Also a single random direction is chosen for sc pointing exactly coincident with the photon direction.
        This amounts to the case of a GRB. The psf is integrated over a circular region with an arbitrary radius.
        This test produces test_response1.rsp.
    */
    void test1();

    /** \brief Simple test much like test 1, except the true energy bins are read from a bin file.
        This test produces test_response2.rsp.
        \param ps_ra The true photon RA in degrees
        \param ps_dec The true photon DEC in degrees
        \param radius The radius of psf integration in degrees.
        \param file_name The name of the output response file.
    */
    void test2(double ps_ra, double ps_dec, double radius, const std::string & file_name);

    /** \brief Like test2, but consistent with a steady point source observed over a period of time. Thus the
        sc data file is read, and sc positions are binned on theta wrt point source direction.
        This test produces test_response3.rsp.
    */
    void test3();

    /** \brief Test utilizing GrbResponse class to compute response for a fictitious grb. The sc data are
        interpolated for the time of the burst to give a sc pointing. The burst was assumed to be close to this
        value so that the angle theta will be non-0 but small.
        This test produces test_response4.rsp.
    */
    void test4();

    /** \brief Repeat test 2, using photon direction identical to that in test 4. Apart from keywords, the spectrum
        in this test should be the same as that written by test 4.
        This test produces test_response5.rsp.
    */
    void test5();

    /** \brief Test RspGenApp class for GrbResponse.
    */
    void test6();

    /** \brief Test utilizing PointResponse class to compute response for a given point source direction. The sc data are
        binned to give a differential exposure. This is then integrated to produce the total response in the sc file.
        This test produces test_response7.rsp.
    */
    void test7();

    /** \brief Test RspGenApp class for PointResponse.
    */
    void test8();

    /** \brief Test extracting information from cuts which are contained in DSS keywords.
    */
    void test9();

    /** \brief Test SpaceCraftCalculator class.
    */
    void test10();

    /** \brief Test addition of several responses.
    */
    void test11();

    /** \brief Test pass 7 irfs and phi binning.
    */
    void test12();

    /// \brief Test calculation of phi.
    void testPhiCalc();

  private:
    /** \brief Compare three arrays to ensure the third is the sum of the first two.
        \param descriptor A string describing the arrays being compared, for diagnostic purposes.
        \param vec1 First array.
        \param vec2 Second array.
        \param total Array which is supposed to be the sum of vec1 and vec2.
    */
    void compare(const std::string & descriptor, const std::vector<double> & vec1, const std::vector<double> & vec2,
      const std::vector<double> & total);

    /** \brief Compare two values to see if they are equal.
        \param descriptor A string describing the comparison being performed.
        \param value The (computed) value.
        \param expected_value The expected value.
        \param abs_tol Tolerance.
    */
    void compare(const std::string & descriptor, double value, double expected_value, double abs_tol);

    /** \brief Return a standard energy binner used throughout the tests.
    */
    evtbin::Binner * createStdBinner();

    /** \brief Return a standard test Irfs object.
    */
    irfInterface::Irfs * createIrfs() const;

    /** \brief Return the full path to the given data file.
    */
    std::string findFile(const std::string & file_root) const;

    void copyFile(const std::string & in_file, const std::string & out_file) const;

    st_stream::StreamFormatter m_os;
    std::string m_data_dir;
    bool m_failed;
};

RspGenTestApp::RspGenTestApp(): m_os("RspGenTestApp", "RspGenTestApp", 2), m_data_dir(), m_failed(false) {
  // Get the directory in which to find the input data files.
  m_data_dir = facilities::commonUtilities::getDataPath("rspgen");
  setName("test_rspgen");
  setVersion(s_cvs_id);
}

void RspGenTestApp::run() {
  // Load all irfs.
  irfLoader::Loader::go();

  // Load instrument-specific binner configurations.
  evtbin::BinConfig::load();

  // Run all tests in order.
  try {
    test1();
  } catch (const std::exception & x) {
    m_failed = true;
    std::cerr << "While running test1, RspGenTestApp caught " << typeid(x).name() << ": what == " << x.what() << std::endl;
  }

  try {
    // Crab pulsar (8.3633225E+01, 2.2014458E+01). Arbitrary radius of psf integration == 1.5
    test2(8.3633225E+01, 2.2014458E+01, 1.5, "test_response2.rsp");
  } catch (const std::exception & x) {
    m_failed = true;
    std::cerr << "While running test2, RspGenTestApp caught " << typeid(x).name() << ": what == " << x.what() << std::endl;
  }

  try {
    test3();
  } catch (const std::exception & x) {
    m_failed = true;
    std::cerr << "While running test3, RspGenTestApp caught " << typeid(x).name() << ": what == " << x.what() << std::endl;
  }

  try {
    test4();
  } catch (const std::exception & x) {
    m_failed = true;
    std::cerr << "While running test4, RspGenTestApp caught " << typeid(x).name() << ": what == " << x.what() << std::endl;
  }

  try {
    test5();
  } catch (const std::exception & x) {
    m_failed = true;
    std::cerr << "While running test5, RspGenTestApp caught " << typeid(x).name() << ": what == " << x.what() << std::endl;
  }

  try {
    test6();
  } catch (const std::exception & x) {
    m_failed = true;
    std::cerr << "While running test6, RspGenTestApp caught " << typeid(x).name() << ": what == " << x.what() << std::endl;
  }

  try {
    test7();
  } catch (const std::exception & x) {
    m_failed = true;
    std::cerr << "While running test7, RspGenTestApp caught " << typeid(x).name() << ": what == " << x.what() << std::endl;
  }

  try {
    test8();
  } catch (const std::exception & x) {
    m_failed = true;
    std::cerr << "While running test8, RspGenTestApp caught " << typeid(x).name() << ": what == " << x.what() << std::endl;
  }

  try {
    test9();
  } catch (const std::exception & x) {
    m_failed = true;
    std::cerr << "While running test9, RspGenTestApp caught " << typeid(x).name() << ": what == " << x.what() << std::endl;
  }

  try {
    test10();
  } catch (const std::exception & x) {
    m_failed = true;
    std::cerr << "While running test10, RspGenTestApp caught " << typeid(x).name() << ": what == " << x.what() << std::endl;
  }

  try {
    test11();
  } catch (const std::exception & x) {
    m_failed = true;
    std::cerr << "While running test11, RspGenTestApp caught " << typeid(x).name() << ": what == " << x.what() << std::endl;
  }

  try {
    test12();
  } catch (const std::exception & x) {
    m_failed = true;
    std::cerr << "While running test12, RspGenTestApp caught " << typeid(x).name() << ": what == " << x.what() << std::endl;
  }

  try {
    testPhiCalc();
  } catch (const std::exception & x) {
    m_failed = true;
    std::cerr << "While running testPhiCalc, RspGenTestApp caught " << typeid(x).name() << ": what == " << x.what() << std::endl;
  }

  if (m_failed) throw std::runtime_error("test_rspgen failed");
}

// Test just getting values for an arbitrary simple case.
void RspGenTestApp::test1() {
  double ra_ps = 8.3633225E+01; // RA of point source
  double dec_ps = 2.2014458E+01; // DEC of point source
  double ra_scz = ra_ps; // RA of spacecraft Z axis
  double dec_scz = dec_ps; // DEC of spacecraft Z axis

  // Open input pha file.
  std::auto_ptr<const tip::Table> ebounds_ext(tip::IFileSvc::instance().readTable(findFile("PHA1.pha"), "EBOUNDS"));

  // Read the detchans keyword. This determines the number of channels used for apparent energy.
  int detchans = 0;
  ebounds_ext->getHeader()["DETCHANS"].get(detchans);

  // Make sure there are at least detchans channels in the input ebounds. This is just a basic sanity check.
  std::size_t num_rec = ebounds_ext->getNumRecords();
  if (num_rec < std::size_t(detchans)) throw std::runtime_error("test1: Channel number mismatch");
  std::vector<double> min_app_en(detchans);
  std::vector<double> max_app_en(detchans);

  std::size_t index = 0;
  for (tip::Table::ConstIterator ebounds_itor = ebounds_ext->begin(); ebounds_itor != ebounds_ext->end(); ++ebounds_itor) {
    (*ebounds_itor)["CHANNEL"].get(index);
    if (index > std::size_t(detchans)) continue; // Skip any rows with channel numbers > the number of channels.
    // Warning: This assumes first channel is 1, but that is not necessarily true. Check TLMIN/TLMAX keywords
    --index; // Arrays start with 0, channels with 1.
    (*ebounds_itor)["E_MIN"].get(min_app_en[index]);
    (*ebounds_itor)["E_MAX"].get(max_app_en[index]);
    min_app_en[index] *= s_MeV_per_keV;
    max_app_en[index] *= s_MeV_per_keV;
  }
  // Add a check here to make sure all the channels were set. We can't compute response if any are missing.

  // Make a couple directions.
  astro::SkyDir ps_dir(ra_ps, dec_ps);
  astro::SkyDir scz_dir(ra_scz, dec_scz);

  // Obtain test response functor.
  irfInterface::Irfs * irfs = createIrfs();

  // First get psf.
  irfInterface::IPsf * psf = irfs->psf();

  // Compute angle, which should be 0.
  double theta = ps_dir.difference(scz_dir) * 180. / M_PI;

  double phi = 0.;
  double radius = 0.1;

  // Create output response file from template:
  tip::IFileSvc::instance().createFile("test_response1.rsp", findFile("LatResponseTemplate"));

  // Open the response file:
  tip::Table * resp_table = tip::IFileSvc::instance().editTable("test_response1.rsp", "MATRIX");

  resp_table->getHeader()["DETCHANS"].set(detchans);

  // Next get aeff.
  irfInterface::IAeff * aeff = irfs->aeff();

  // And redistribution.
  irfInterface::IEdisp * edisp = irfs->edisp();

  tip::Table::Iterator out_itor = resp_table->begin();

  // Arbitrarily set the number of channels used for "true" energy dimension.
  int num_true_chan = 100;

  // Arrays used to write f_chan and n_chan columns, respectively.
  std::vector<int> f_chan(1, 1);
  std::vector<int> n_chan(1, detchans);

  // Create fake true energy bins, logarithmically spanning the GLAST spectrum.
  double offset = 2.;
  double factor = (6. - offset)/num_true_chan;
  for (long true_en_idx = 0; true_en_idx < num_true_chan; ++true_en_idx, ++out_itor) {
    double true_en = pow(10., (true_en_idx + .5) * factor + offset);
    double aeff_val = aeff->value(true_en, theta, phi);
    double int_psf_val = psf->angularIntegral(true_en, theta, phi, radius);

    // Populate response vector for each true energy value.
    std::vector<double> response;
    for (index = 0; index < std::size_t(detchans); ++index) {
      response.push_back(aeff_val * edisp->integral(min_app_en[index], max_app_en[index], true_en, theta, phi) * int_psf_val);
    }

    // Write response to file, converting to keV on the fly.
    (*out_itor)["ENERG_LO"].set(s_keV_per_MeV * pow(10., true_en_idx * factor + offset));
    (*out_itor)["ENERG_HI"].set(s_keV_per_MeV * pow(10., (true_en_idx + 1) * factor + offset));
    (*out_itor)["N_GRP"].set(1);
    (*out_itor)["F_CHAN"].set(f_chan);
    (*out_itor)["N_CHAN"].set(n_chan);
    (*out_itor)["MATRIX"].set(response);
  }

  // Write the ebounds extension.
  tip::Table * out_ebounds = tip::IFileSvc::instance().editTable("test_response1.rsp", "EBOUNDS");

  // Set detchans explicitly.
  out_ebounds->getHeader()["DETCHANS"].set(detchans);

  out_ebounds->setNumRecords(num_rec);

  // Just copy the input ebounds extension.
  tip::Table::ConstIterator in_itor = ebounds_ext->begin();
  for (out_itor = out_ebounds->begin(); out_itor != out_ebounds->end(); ++in_itor, ++out_itor) {
    *out_itor = *in_itor;
  }

  delete out_ebounds;
  delete resp_table;
  delete irfs;
}

// In this example. "true" energy comes from bin definition file.
void RspGenTestApp::test2(double ra_ps, double dec_ps, double radius, const std::string & file_name) {
  using namespace evtbin;

  // Make sure output will reveal small differences.
  std::cout.precision(24);

  // Create interval container for user defined bin intervals.
  OrderedBinner::IntervalCont_t intervals;

  // Open the file for true energy bin definition.
  std::auto_ptr<const tip::Table> table(tip::IFileSvc::instance().readTable(findFile("rspgen_energy_bins.fits"), "ENERGYBINS"));

  // Iterate over the file, saving the relevant values into the interval array, converting to MeV on the fly.
  for (tip::Table::ConstIterator itor = table->begin(); itor != table->end(); ++itor) {
    intervals.push_back(Binner::Interval(s_MeV_per_keV*(*itor)["E_MIN"].get(), s_MeV_per_keV*(*itor)["E_MAX"].get()));
  }

  // Create binner from these intervals.
  OrderedBinner binner(intervals);

  long num_bins = binner.getNumBins();

  double ra_scz = ra_ps; // RA of spacecraft Z axis
  double dec_scz = dec_ps; // DEC of spacecraft Z axis

  // Open input ebounds extension, used for apparent energy bins.
  std::auto_ptr<const tip::Table> ebounds_ext(tip::IFileSvc::instance().readTable(findFile("PHA1.pha"), "EBOUNDS"));

  // Read the detchans keyword. This determines the number of channels used for apparent energy.
  int detchans = 0;
  ebounds_ext->getHeader()["DETCHANS"].get(detchans);

  // Make sure there are at least detchans channels in the input ebounds. This is just a basic sanity check.
  std::size_t num_rec = ebounds_ext->getNumRecords();
  if (num_rec < std::size_t(detchans)) throw std::runtime_error("test2: Channel number mismatch");
  std::vector<double> min_app_en(detchans);
  std::vector<double> max_app_en(detchans);

  std::size_t index = 0;
  for (tip::Table::ConstIterator ebounds_itor = ebounds_ext->begin(); ebounds_itor != ebounds_ext->end(); ++ebounds_itor) {
    (*ebounds_itor)["CHANNEL"].get(index);
    if (index > std::size_t(detchans)) continue; // Skip any rows with channel numbers > the number of channels.
    // Warning: This assumes first channel is 1, but that is not necessarily true. Check TLMIN/TLMAX keywords
    --index; // Arrays start with 0, channels with 1.
    (*ebounds_itor)["E_MIN"].get(min_app_en[index]);
    (*ebounds_itor)["E_MAX"].get(max_app_en[index]);
    min_app_en[index] *= s_MeV_per_keV;
    max_app_en[index] *= s_MeV_per_keV;
  }
  // Add a check here to make sure all the channels were set. We can't compute response if any are missing.

  astro::SkyDir ps_dir(ra_ps, dec_ps);
  astro::SkyDir scz_dir(ra_scz, dec_scz);

  // Obtain test response functor.
  irfInterface::Irfs * irfs = createIrfs();

  // First get psf.
  irfInterface::IPsf * psf = irfs->psf();

  // Compute angle, which should be 0.
  double theta = ps_dir.difference(scz_dir) * 180. / M_PI;

  double phi = 0.;

  // Create output response file from template:
  tip::IFileSvc::instance().createFile(file_name, findFile("LatResponseTemplate"));

  // Open the response file:
  tip::Table * resp_table = tip::IFileSvc::instance().editTable(file_name, "MATRIX");

  resp_table->getHeader()["DETCHANS"].set(detchans);

  // Next get aeff.
  irfInterface::IAeff * aeff = irfs->aeff();

  // And redistribution.
  irfInterface::IEdisp * edisp = irfs->edisp();

  tip::Table::Iterator out_itor = resp_table->begin();

  // Number of true energy channels is not arbitrary, but is taken from the true energy binning info obtained above.
  int num_true_chan = num_bins;

  // Arrays used to write f_chan and n_chan columns, respectively.
  std::vector<int> f_chan(1, 1);
  std::vector<int> n_chan(1, detchans);
  for (long true_en_idx = 0; true_en_idx < num_true_chan; ++true_en_idx, ++out_itor) {
    double true_en = binner.getInterval(true_en_idx).midpoint();
    double aeff_val = aeff->value(true_en, theta, phi);
    double int_psf_val = psf->angularIntegral(true_en, theta, phi, radius);

    // Populate response vector.
    std::vector<double> response;
    for (index = 0; index < std::size_t(detchans); ++index) {
      response.push_back(aeff_val * edisp->integral(min_app_en[index], max_app_en[index], true_en, theta, phi) * int_psf_val);
    }

    // Write response to file, using keV.
    (*out_itor)["ENERG_LO"].set(s_keV_per_MeV * binner.getInterval(true_en_idx).begin());
    (*out_itor)["ENERG_HI"].set(s_keV_per_MeV * binner.getInterval(true_en_idx).end());
    (*out_itor)["N_GRP"].set(1);
    (*out_itor)["F_CHAN"].set(f_chan);
    (*out_itor)["N_CHAN"].set(n_chan);
    (*out_itor)["MATRIX"].set(response);
  }

  // Copy input ebounds extension to the output.
  tip::Table * out_ebounds = tip::IFileSvc::instance().editTable(file_name, "EBOUNDS");

  // Set detchans explicitly.
  out_ebounds->getHeader()["DETCHANS"].set(detchans);

  out_ebounds->setNumRecords(num_rec);

  tip::Table::ConstIterator in_itor = ebounds_ext->begin();
  for (out_itor = out_ebounds->begin(); out_itor != out_ebounds->end(); ++in_itor, ++out_itor) {
    *out_itor = *in_itor;
  }

  delete out_ebounds;
  delete resp_table;
  delete irfs;
}

// In this example. "true" energy comes from bin definition file.
// Also, the FT2 file is read and a histogram in time as a function of angle is created.
void RspGenTestApp::test3() {
  using namespace evtbin;

  double ra_ps = 8.3633225E+01; // RA of point source
  double dec_ps = 2.2014458E+01; // DEC of point source
  astro::SkyDir ps_pos(ra_ps, dec_ps);

  // Open the sc data file.
  std::auto_ptr<const tip::Table> sc_data(tip::IFileSvc::instance().readTable(findFile("ft2tiny.fits"), "Ext1"));

  // Set up a histogram to hold the binned differential exposure (theta vs. DeltaT).
  Hist1D diff_exp(LinearBinner(0., 60., 5.));

  double total_exposure = 0.;

  // Read SC Z positions, bin them into a histogram:
  for (tip::Table::ConstIterator itor = sc_data->begin(); itor != sc_data->end(); ++itor) {
    // Get size of interval.
    double delta_t = (*itor)["LIVETIME"].get();

    // Get SC coordinates.
    double ra_scz = (*itor)["RA_SCZ"].get();
    double dec_scz = (*itor)["DEC_SCZ"].get();
    astro::SkyDir scz_pos(ra_scz, dec_scz);

    // Compute inclination angle from the source position to the sc.
    double theta = ps_pos.difference(scz_pos) * 180. / M_PI;

    // Bin this angle into the histogram.
    diff_exp.fillBin(theta, delta_t);

    total_exposure += delta_t;
  }

  // Confirm that something was accumulated.
  if (0. == total_exposure) throw std::runtime_error("test3 cannot continue with 0. total exposure");

  // Make sure small differences will be printed correctly.
  std::cout.precision(24);

#ifdef foo
  std::cout << "***************** histogram values ***************" << std::endl;
  for (int ii = 0; ii < 12; ++ii)
    std::cout << diff_exp[ii]/(2.5 + ii * 5.) << std::endl;
  std::cout << "***************** histogram values ***************" << std::endl;
#endif

  // Create interval container for user defined bin intervals.
  OrderedBinner::IntervalCont_t intervals;

  // Open the data file.
  std::auto_ptr<const tip::Table> table(tip::IFileSvc::instance().readTable(findFile("rspgen_energy_bins.fits"), "ENERGYBINS"));

  // Iterate over the file, saving the relevant values into the interval array, converting to MeV on the fly.
  for (tip::Table::ConstIterator itor = table->begin(); itor != table->end(); ++itor) {
    intervals.push_back(Binner::Interval(s_MeV_per_keV * (*itor)["E_MIN"].get(), s_MeV_per_keV * (*itor)["E_MAX"].get()));
  }

  // Create binner from these intervals.
  OrderedBinner binner(intervals);

  long num_bins = binner.getNumBins();

  // Open input pha file.
  std::auto_ptr<const tip::Table> ebounds_ext(tip::IFileSvc::instance().readTable(findFile("PHA1.pha"), "EBOUNDS"));

  // Read the detchans keyword. This determines the number of channels used for apparent energy.
  int detchans = 0;
  ebounds_ext->getHeader()["DETCHANS"].get(detchans);

  // Make sure there are at least detchans channels in the input ebounds. This is just a basic sanity check.
  std::size_t num_rec = ebounds_ext->getNumRecords();
  if (num_rec < std::size_t(detchans)) throw std::runtime_error("test3: Channel number mismatch");
  std::vector<double> min_app_en(detchans);
  std::vector<double> max_app_en(detchans);

  std::size_t index = 0;
  for (tip::Table::ConstIterator ebounds_itor = ebounds_ext->begin(); ebounds_itor != ebounds_ext->end(); ++ebounds_itor) {
    (*ebounds_itor)["CHANNEL"].get(index);
    if (index > std::size_t(detchans)) continue; // Skip any rows with channel numbers > the number of channels.
    // Warning: This assumes first channel is 1, but that is not necessarily true. Check TLMIN/TLMAX keywords
    --index; // Arrays start with 0, channels with 1.
    (*ebounds_itor)["E_MIN"].get(min_app_en[index]);
    (*ebounds_itor)["E_MAX"].get(max_app_en[index]);
    min_app_en[index] *= s_MeV_per_keV;
    max_app_en[index] *= s_MeV_per_keV;
  }
  // Add a check here to make sure all the channels were set. We can't compute response if any are missing.

  // Obtain test response functor.
  irfInterface::Irfs * irfs = createIrfs();

  // First get psf.
  irfInterface::IPsf * psf = irfs->psf();

  // Theta varies over the histogram now.
  // double theta = ps_dir.difference(scz_dir) * 180. / M_PI;

  double phi = 0.;
  double radius = 0.1;

  // Create output response file from template:
  tip::IFileSvc::instance().createFile("test_response3.rsp", findFile("LatResponseTemplate"));

  // Open the response file:
  tip::Table * resp_table = tip::IFileSvc::instance().editTable("test_response3.rsp", "MATRIX");

  resp_table->getHeader()["DETCHANS"].set(detchans);

  // Next get aeff.
  irfInterface::IAeff * aeff = irfs->aeff();

  // And redistribution.
  irfInterface::IEdisp * edisp = irfs->edisp();

  tip::Table::Iterator out_itor = resp_table->begin();

  // Arbitrarily set the number of channels used for "true" energy dimension.
  int num_true_chan = num_bins;

  // Arrays used to write f_chan and n_chan columns, respectively.
  std::vector<int> f_chan(1, 1);
  std::vector<int> n_chan(1, detchans);

  // Iterate over the bins of the histogram.
  const Binner * theta_bins = diff_exp.getBinners()[0];
  long num_theta_bins = theta_bins->getNumBins();

  for (long true_en_idx = 0; true_en_idx < num_true_chan; ++true_en_idx, ++out_itor) {
    double true_en = binner.getInterval(true_en_idx).midpoint();

    // Populate response vector.
    std::vector<double> response(detchans, 0.);

    // Integrate over binned angle-dependent differential exposure.
    for (long theta_bin = 0; theta_bin < num_theta_bins; ++theta_bin) {
      double theta = theta_bins->getInterval(theta_bin).midpoint();
      double aeff_val = aeff->value(true_en, theta, phi);
      double int_psf_val = psf->angularIntegral(true_en, theta, phi, radius);

      for (index = 0; index < std::size_t(detchans); ++index) {
        response[index] += diff_exp[theta_bin] / total_exposure * aeff_val * edisp->integral(min_app_en[index], max_app_en[index], true_en, theta, phi) * int_psf_val;
      }

      // Write response to file, using keV.
      (*out_itor)["ENERG_LO"].set(s_keV_per_MeV * binner.getInterval(true_en_idx).begin());
      (*out_itor)["ENERG_HI"].set(s_keV_per_MeV * binner.getInterval(true_en_idx).end());
      (*out_itor)["N_GRP"].set(1);
      (*out_itor)["F_CHAN"].set(f_chan);
      (*out_itor)["N_CHAN"].set(n_chan);
      (*out_itor)["MATRIX"].set(response);
    }
  }

  // Copy ebounds extension from input to output.
  tip::Table * out_ebounds = tip::IFileSvc::instance().editTable("test_response3.rsp", "EBOUNDS");

  // Set detchans explicitly.
  out_ebounds->getHeader()["DETCHANS"].set(detchans);

  out_ebounds->setNumRecords(num_rec);

  tip::Table::ConstIterator in_itor = ebounds_ext->begin();
  for (out_itor = out_ebounds->begin(); out_itor != out_ebounds->end(); ++in_itor, ++out_itor) {
    *out_itor = *in_itor;
  }

  delete out_ebounds;
  delete resp_table;
  delete irfs;
}

// Test GrbResponse class constructors and compute method.
void RspGenTestApp::test4() {
  using namespace rspgen;

  try {
    // This test assumes there was a burst at (114., -30.), at t = 105. seconds past the first entry in the spacecraft file.

    // Get spacecraft data.
    std::auto_ptr<const tip::Table> sc_table(tip::IFileSvc::instance().readTable(findFile("ft2tiny.fits"), "Ext1"));

    // Get object for interpolating values from the table.
    tip::LinearInterp sc_record(sc_table->begin(), sc_table->end());

    // Interpolate values for a burst.
    sc_record.interpolate("START", 2.167440000000000E+06 + 105.);

    // Compute inclination angle from burst RA and DEC and spacecraft pointing.
    double theta = astro::SkyDir(114., -30.).difference(astro::SkyDir(sc_record.get("RA_SCZ"), sc_record.get("DEC_SCZ")))*180./M_PI;

    double phi = 0.;

    // Confirm that this value is correct.
    float correct_val = 1.18774705e02;
    if (correct_val != float(sc_record.get("RA_SCZ"))) {
      m_failed = true;
      std::cerr << "Unexpected: in test4, interpolated RA_SCZ was " << sc_record.get("RA_SCZ") << ", not " << correct_val
        << std::endl;
    }


    // Obtain test response functor.
    std::auto_ptr<irfInterface::Irfs> irfs(createIrfs());


    // Create window object for psf integration, a circle of radius 1.4 degrees.
    CircularWindow window(1.4);


    // Get input ebounds extension.
    std::auto_ptr<const tip::Table> in_ebounds(tip::IFileSvc::instance().readTable(findFile("PHA1.pha"), "EBOUNDS"));

    // Get number of channels currently in use.
    int detchans = 0;
    in_ebounds->getHeader()["DETCHANS"].get(detchans);

    // Get apparent energy binner from ebounds extension.
    int index = 0;
    evtbin::OrderedBinner::IntervalCont_t app_intervals(detchans);
    for (tip::Table::ConstIterator itor = in_ebounds->begin(); itor != in_ebounds->end() && index < detchans; ++itor, ++index)
      app_intervals[index] = evtbin::Binner::Interval(s_MeV_per_keV*(*itor)["E_MIN"].get(), s_MeV_per_keV*(*itor)["E_MAX"].get());

    // Create apparent energy binner.
    evtbin::OrderedBinner app_en_binner(app_intervals);


    // Process bin definition file to produce true energy bins.
    std::auto_ptr<const tip::Table> true_en(tip::IFileSvc::instance().readTable(findFile("rspgen_energy_bins.fits"), "ENERGYBINS"));

    // Create object to hold the intervals in the bin definition file.
    evtbin::OrderedBinner::IntervalCont_t true_intervals(true_en->getNumRecords());

    // Read true energy bin definitions into interval object.
    index = 0;
    for (tip::Table::ConstIterator itor = true_en->begin(); itor != true_en->end(); ++itor, ++index)
      true_intervals[index] = evtbin::Binner::Interval(s_MeV_per_keV*(*itor)["E_MIN"].get(), s_MeV_per_keV*(*itor)["E_MAX"].get());

    // Create a binner for true energy.
    evtbin::OrderedBinner true_en_binner(true_intervals);

    // Create response object for burst, using the first constructor.
    GrbResponse resp(theta, phi, &true_en_binner, &app_en_binner, irfs.get(), &window);

    // Sanity check: just compute one response at 137. MeV.
    std::vector<double> resp_slice(detchans, 0.);
    resp.compute(137., resp_slice);

    // Use the second constructor, and test that its results are the same.
    // RA, DEC, t, radius, response type, pha file, FT2 file, energy bin def file
    GrbResponse resp2(114., -30., 2.167440000000000E+06 + 105., 1.4, "testIrfs::Front", findFile("PHA1.pha"),
      findFile("ft2tiny.fits"), "Ext1", &true_en_binner);

    // Sanity check: just compute one response at 137. MeV.
    std::vector<double> resp_slice2(detchans, 0.);
    resp2.compute(137., resp_slice2);

    // Make sure the two response objects agree with their computations.
    if (resp_slice != resp_slice2) {
      m_failed = true;
      std::cerr << "Unexpected: GrbResponse objects constructed differently from same data give different results" << std::endl;
    }

    // Look for at least one non-zero value in the response computation.
    std::vector<double>::iterator resp_itor;
    for (resp_itor = resp_slice.begin(); resp_itor != resp_slice.end(); ++resp_itor)
      if (0. != *resp_itor) break;

    if (resp_slice.end() == resp_itor) {
      m_failed = true;
      std::cerr << "Unexpected: GrbResponse computed an all zero response slice at 137. MeV" << std::endl;
    }

    // Write output rsp file.
    resp2.writeOutput("test_rspgen", "test_response4.rsp", findFile("LatResponseTemplate"));

  } catch (const std::exception & x) {
    m_failed = true;
    std::cerr << "Unexpected: test4 caught " << typeid(x).name() << ": " << x.what() << std::endl;
  }

}

// Verify that test2 produces same result as GrbResponse class in test4.
void RspGenTestApp::test5() {
  using namespace rspgen;

  try {
    // Repeat test2 using test4 parameters.
    test2(114., -30., 1.4, "test_response5.rsp");
  } catch (const std::exception & x) {
    m_failed = true;
    std::cerr << "Unexpected: test5 caught " << typeid(x).name() << ": " << x.what() << std::endl;
  }

}

// Test application class for GrbResponse case.
void RspGenTestApp::test6() {
  using namespace rspgen;

  try {
    RspGenApp app;

    st_app::AppParGroup & pars(getParGroup("gtrspgen"));

    // Make copy of input spectrum file so that the original spectrum will not be changed by
    // this test.
    std::string orig_spec = findFile("PHA1.pha");
    copyFile(orig_spec, "PHA1.pha");

    // Write-protect copy of spectrum to make sure response still gets written even if
    // keywords in spectrum cannot be updated.
#ifndef WIN32
    // Set mode 0444.
    chmod("PHA1.pha", S_IRUSR | S_IRGRP | S_IROTH);
#endif

    // Set parameters "by hand"
    pars["respalg"] = "GRB";
    pars["specfile"] = "PHA1.pha";
    pars["scfile"] = findFile("ft2tiny.fits");
    pars["sctable"] = "Ext1";
    pars["outfile"] = "test_response6.rsp";
    pars["time"] = 2.167440000000000E+06 + 105.;
    pars["irfs"] = "testIrfs::Front";
    pars["resptpl"] = "DEFAULT";
    pars["ebinalg"] = "FILE";
    pars["ebinfile"] = findFile("rspgen_energy_bins.fits");

    // And writing the output.
    app.writeResponse(pars);

#ifndef WIN32
    // Make sure RESPFILE keyword was *not* written to spectrum file, which should be write protected.
    std::auto_ptr<const tip::Table> spec(tip::IFileSvc::instance().readTable("PHA1.pha", "SPECTRUM"));
    std::string resp_file;
    spec->getHeader()["RESPFILE"].get(resp_file);
    if (resp_file != "NONE") {
      m_failed = true;
      std::cerr << "Unexpected: after writing response, test6 found that RESPFILE keyword in PHA1.pha was " << resp_file <<
        ", not \"NONE\", as expected." << std::endl;
    }

    // Set mode 0664.
    chmod("PHA1.pha", S_IRUSR | S_IRGRP | S_IROTH | S_IWUSR | S_IWGRP);
#endif

  } catch (const std::exception & x) {
    m_failed = true;
    std::cerr << "Unexpected: test6 caught " << typeid(x).name() << ": " << x.what() << std::endl;
  }

}

// Test PointResponse class constructors and compute method. Parameters are such that
// this test's output should be identical to test3's output.
void RspGenTestApp::test7() {
  using namespace rspgen;

  std::auto_ptr<PointResponse> response_no_binner(0);
  try {
    // Construct a steady point source response, but do not supply a binner. This should succeed.
    double ra_ps = 8.3633225E+01; // RA of point source
    double dec_ps = 2.2014458E+01; // DEC of point source
    response_no_binner.reset(new PointResponse(ra_ps, dec_ps, 60., .05, 3., 1, "testIrfs::Back", findFile("PHA1.pha"),
      findFile("ft2tiny.fits"), "Ext1", 0));

    // Confirm that psf method works even without the binner.
    response_no_binner->psf(300., 15., 0.);
  } catch (const std::exception & x) {
    m_failed = true;
    std::cerr << "Unexpected: test7 failed to create or use PointResponse object with null binner:" << x.what() << std::endl;
  }

  if (0 != response_no_binner.get()) {
    try {
      // Write output, which should fail because there is no true energy binner.
      response_no_binner->writeOutput("test_rspgen", "test_response7.rsp", findFile("LatResponseTemplate"));
      m_failed = true;
      std::cerr << "Unexpected: test7 did not catch exception when writeOutput called but no true energy binner exists." <<
        std::endl;
    } catch (const std::exception &) {
    }
  }

  try {
    // Create a binner for true energy.
    std::auto_ptr<evtbin::Binner> true_en_binner(createStdBinner());

    // Construct a steady point source response for the given RA, DEC, thetabins.
    // RA, DEC, theta_cut, theta_bin_size, radius, number of phi bins, response type, pha file, FT2 file, energy bin def file.
    double ra_ps = 8.3633225E+01; // RA of point source
    double dec_ps = 2.2014458E+01; // DEC of point source
    PointResponse resp(ra_ps, dec_ps, 60., .05, 3., 80, "testIrfs::Back", findFile("PHA1.pha"), findFile("ft2tiny.fits"),
      "Ext1", true_en_binner.get());

    // Confirm binner has expected numbers of bins.
    std::pair<long, long> dims(resp.getSpatialNumBins());
    if (dims.first != 10) {
      m_failed = true;
      std::cerr << "Unexpected: PointResponse::getSpatialNumBins reports " << dims.first <<
        " theta bins, not 10 as expected." << std::endl;
    }
    if (dims.second != 80) {
      m_failed = true;
      std::cerr << "Unexpected: PointResponse::getSpatialNumBins reports " << dims.second <<
        " phi bins, not 80 as expected." << std::endl;
    }

    // Sanity check: just compute one response at 137. MeV.
    std::vector<double> resp_slice;
    resp.compute(137., resp_slice);

    // Look for at least one non-zero value in the response computation.
    std::vector<double>::iterator resp_itor;
    for (resp_itor = resp_slice.begin(); resp_itor != resp_slice.end(); ++resp_itor)
      if (0. != *resp_itor) break;

    if (resp_slice.end() == resp_itor) {
      m_failed = true;
      std::cerr << "Unexpected: PointResponse computed an all zero response slice at 137. MeV" << std::endl;
    }

    // Write output rsp file.
    resp.writeOutput("test_rspgen", "test_response7.rsp", findFile("LatResponseTemplate"));

  } catch (const std::exception & x) {
    m_failed = true;
    std::cerr << "Unexpected: test7 caught " << typeid(x).name() << ": " << x.what() << std::endl;
  }

}

// Test application class for PointResponse case.
void RspGenTestApp::test8() {
  using namespace rspgen;

  try {
    RspGenApp app;

    st_app::AppParGroup & pars(getParGroup("gtrspgen"));

    // Make copy of input spectrum file so that the original spectrum will not be changed by
    // this test.
    std::string orig_spec = findFile("PHA1.pha");
    copyFile(orig_spec, "PHA1.pha");

    // Set parameters "by hand"
    pars["respalg"] = "PS";
    pars["specfile"] = "PHA1.pha";
    pars["scfile"] = findFile("ft2tiny.fits");
    pars["sctable"] = "Ext1";
    pars["outfile"] = "test_response8.rsp";
    pars["thetacut"] = 60.;
    pars["dcostheta"] = .05;
    pars["phinumbins"] = 80;
    pars["irfs"] = "testIrfs::Front";
    pars["resptpl"] = "DEFAULT";
    pars["ebinalg"] = "FILE";
    pars["ebinfile"] = findFile("rspgen_energy_bins.fits");

    // And writing the output.
    app.writeResponse(pars);

    // Make sure RESPFILE keyword was written to spectrum file.
    std::auto_ptr<const tip::Table> spec(tip::IFileSvc::instance().readTable("PHA1.pha", "SPECTRUM"));
    std::string resp_file;
    spec->getHeader()["RESPFILE"].get(resp_file);
    if (resp_file != "test_response8.rsp") {
      m_failed = true;
      std::cerr << "Unexpected: after writing response, test8 found that RESPFILE keyword in PHA1.pha was " << resp_file <<
        ", not \"test_response8.rsp\", as expected." << std::endl;
    }

  } catch (const std::exception & x) {
    m_failed = true;
    std::cerr << "Unexpected: test8 caught " << typeid(x).name() << ": " << x.what() << std::endl;
  }

}

// Test extracting information from cuts which are contained in DSS keywords.
void RspGenTestApp::test9() {
  using namespace dataSubselector;
  using namespace rspgen;

  std::string orig_spec = findFile("PHA1.pha");

  // Read cuts from spectrum.
  Cuts cuts(orig_spec, "SPECTRUM", false, true);

  // Confirm cuts contain a single sky cone centered on the Crab.
  std::vector<Cuts *>::size_type num_cone = 0;
  double ra = 0.;
  double dec = 0.;
  double radius = 0.;
  bool cuts_failed = false;

  // Iterate over all cuts.
  for (std::vector<Cuts *>::size_type ii = 0; ii != cuts.size(); ++ii) {
    const SkyConeCut * sky_cut = 0;
    if (0 != (sky_cut = dynamic_cast<const SkyConeCut *>(&cuts[ii]))) {
      ++num_cone;
      ra = sky_cut->ra();
      dec = sky_cut->dec();
      radius = sky_cut->radius();
    }
  }

  // Confirm single sky cone.
  if (1 != num_cone) {
    cuts_failed = true;
    std::cerr << "Unexpected: DSS keywords contained " << num_cone << " sky cones, not 1 as expected." << std::endl;
  }

  double ra_expected = 83.6332;
  double dec_expected = 22.0145;
  double radius_expected = 12.;
  // Confirm sky cone parameters.
  if (ra != ra_expected) {
    cuts_failed = true;
    std::cerr << "Unexpected: DSS keywords contained ra of " << ra << " not " << ra_expected << ", as expected." << std::endl;
  }
  if (dec != dec_expected) {
    cuts_failed = true;
    std::cerr << "Unexpected: DSS keywords contained dec of " << dec << " not " << dec_expected << ", as expected." << std::endl;
  }
  if (radius != radius_expected) {
    cuts_failed = true;
    std::cerr << "Unexpected: DSS keywords contained radius of " << radius << " not " << radius_expected << ", as expected." <<
      std::endl;
  }

  // Flag global error.
  if (cuts_failed) m_failed = true;
}

void RspGenTestApp::test10() {
  using namespace rspgen;
  double src_ra = 83.6332;
  double src_dec = 22.0145;
  double theta_cut = 20.;
  double theta_bin_size = 0.5;
  double psf_radius = 30.;
  std::string resp_type = "TESTF";
  std::string gti_file = findFile("PHA1.pha");
  std::string sc_file = findFile("ft2tiny.fits");

  try {
    SpaceCraftCalculator calc(astro::SkyDir(src_ra, src_dec), theta_cut, theta_bin_size, psf_radius,
      resp_type, gti_file, sc_file, "Ext1");

    double true_energy = 100.;
    double theta = 15.;
    double phi = 0.;
    calc.psf(true_energy, theta, phi);
  } catch (const std::exception & x) {
    m_failed = true;
    std::cerr << "Unexpected: test10 failed: " << x.what() << std::endl;
  }
}

void RspGenTestApp::test11() {
  using namespace rspgen;
  m_os.setMethod("test11");

  SpaceCraftCalculator::displayIrfNames(m_os.info());

  // Create a binner for true energy.
  std::auto_ptr<evtbin::Binner> true_en_binner(createStdBinner());

  // Construct a steady point source response for the given RA, DEC, thetabins.
  // RA, DEC, theta_cut, theta_bin_size, radius, response type, pha file, FT2 file, energy bin def file.
  double ra_ps = 8.3633225E+01; // RA of point source
  double dec_ps = 2.2014458E+01; // DEC of point source
  PointResponse back(ra_ps, dec_ps, 60., .05, 3., 80, "TESTB", findFile("PHA1.pha"), findFile("ft2tiny.fits"),
    "Ext1", true_en_binner.get());
  PointResponse front(ra_ps, dec_ps, 60., .05, 3., 80, "TESTF", findFile("PHA1.pha"), findFile("ft2tiny.fits"),
    "Ext1", true_en_binner.get());
  PointResponse total(ra_ps, dec_ps, 60., .05, 3., 80, "TEST", findFile("PHA1.pha"), findFile("ft2tiny.fits"),
    "Ext1", true_en_binner.get());

  std::vector<double> back_resp;
  std::vector<double> front_resp;
  std::vector<double> total_resp;

  double energy = 800.;

  back.compute(energy, back_resp);
  front.compute(energy, front_resp);
  total.compute(energy, total_resp);

  compare("PointResponse", back_resp, front_resp, total_resp);

  GrbResponse grb_back(114., -30., 2.167440000000000E+06 + 105., 1.4, "TESTB", findFile("PHA1.pha"),
    findFile("ft2tiny.fits"), "Ext1", true_en_binner.get());
  GrbResponse grb_front(114., -30., 2.167440000000000E+06 + 105., 1.4, "testIrfs::Front", findFile("PHA1.pha"),
    findFile("ft2tiny.fits"), "Ext1", true_en_binner.get());
  GrbResponse grb_total(114., -30., 2.167440000000000E+06 + 105., 1.4, "TEST", findFile("PHA1.pha"),
    findFile("ft2tiny.fits"), "Ext1", true_en_binner.get());

  grb_back.compute(energy, back_resp);
  grb_front.compute(energy, front_resp);
  grb_total.compute(energy, total_resp);

  compare("GrbResponse", back_resp, front_resp, total_resp);
}

// Test application class for PointResponse case using Pass 7 irfs and phi binning.
void RspGenTestApp::test12() {
#if 1
    std::cerr << "Skipping test12 because Pass 7 irfs were removed." << std::endl;
#else
  using namespace rspgen;

  try {
    RspGenApp app;

    st_app::AppParGroup & pars(getParGroup("gtrspgen"));

    // Make copy of input spectrum file so that the original spectrum will not be changed by
    // this test.
    std::string orig_spec = findFile("PHA1.pha");
    copyFile(orig_spec, "PHA1.pha");

    // Set parameters "by hand"
    pars["respalg"] = "PS";
    pars["specfile"] = "PHA1.pha";
    pars["scfile"] = findFile("ft2tiny.fits");
    pars["sctable"] = "Ext1";
    pars["outfile"] = "test_response_p7.rsp";
    pars["thetacut"] = 60.;
    pars["dcostheta"] = .05;
    pars["phinumbins"] = 20;
    pars["irfs"] = "P7_V1_DIFFUSE";
    pars["resptpl"] = "DEFAULT";
    pars["ebinalg"] = "FILE";
    pars["ebinfile"] = findFile("rspgen_energy_bins.fits");

    // And writing the output.
    app.writeResponse(pars);

    // Make sure RESPFILE keyword was written to spectrum file.
    std::auto_ptr<const tip::Table> spec(tip::IFileSvc::instance().readTable("PHA1.pha", "SPECTRUM"));
    std::string resp_file;
    spec->getHeader()["RESPFILE"].get(resp_file);
    if (resp_file != "test_response_p7.rsp") {
      m_failed = true;
      std::cerr << "Unexpected: after writing response, test12 found that RESPFILE keyword in PHA1.pha was " << resp_file <<
        ", not \"test_response_p7.rsp\", as expected." << std::endl;
    }

  } catch (const std::exception & x) {
    m_failed = true;
    std::cerr << "Unexpected: test12 caught " << typeid(x).name() << ": " << x.what() << std::endl;
  }
#endif
}


namespace {
  double s_calcPhi(const astro::SkyDir & x_ref, const astro::SkyDir & z_ref, const astro::SkyDir & dir) {
    typedef CLHEP::Hep3Vector vec_t;
    //static const double pi = std::acos(-1);
    static const double pi = M_PI;
    const vec_t & x_hat = x_ref.dir();
    const vec_t y_hat = z_ref.dir().cross(x_hat);
    double phi = std::atan2(dir.dir().dot(y_hat), dir.dir().dot(x_hat));
    while (phi < 0.) phi += 2 * pi;
    return phi * 180. / pi;
  }
}

void RspGenTestApp::testPhiCalc() {
  typedef CLHEP::Hep3Vector vec_t;
  //const double pi = std::acos(-1);
  const double pi = M_PI;

  // Point spacecraft X along the X axis to make life easy.
  astro::SkyDir sc_x(vec_t(1., 0., 0.));
  astro::SkyDir sc_z(vec_t(0., 0., 1.));

  // Set source angles.
  double src_theta = pi / 4.;
  double src_phi = pi / 3.;

  // Unit vector components for source with the above angles.
  double abs_x = std::cos(src_phi) * std::sin(src_theta);
  double abs_y = std::sin(src_phi) * std::sin(src_theta);
  double abs_z = std::cos(src_theta);

  // Compute phi from product of vector: SpacecraftZ cross SourceX dot SpacecraftY.
  // Test vectors in each octant.
  double expected_phi = src_phi; // The right answer.

  double phi = s_calcPhi(sc_x, sc_z, astro::SkyDir(vec_t(abs_x, abs_y, abs_z)));
  compare("in testPhiCalc, phi in octant 1", phi * pi / 180., expected_phi, 1.e-12);

  phi = s_calcPhi(sc_x, sc_z, astro::SkyDir(vec_t(-abs_x, abs_y, abs_z)));
  compare("in testPhiCalc, phi in octant 2", phi * pi / 180., pi - expected_phi, 1.e-12);

  phi = s_calcPhi(sc_x, sc_z, astro::SkyDir(vec_t(-abs_x, -abs_y, abs_z)));
  compare("in testPhiCalc, phi in octant 3", phi * pi / 180., pi + expected_phi, 1.e-12);

  phi = s_calcPhi(sc_x, sc_z, astro::SkyDir(vec_t(abs_x, -abs_y, abs_z)));
  compare("in testPhiCalc, phi in octant 4", phi * pi / 180., 2 * pi - expected_phi, 1.e-12);

  phi = s_calcPhi(sc_x, sc_z, astro::SkyDir(vec_t(abs_x, abs_y, -abs_z)));
  compare("in testPhiCalc, phi in octant 5", phi * pi / 180., expected_phi, 1.e-12);

  phi = s_calcPhi(sc_x, sc_z, astro::SkyDir(vec_t(-abs_x, abs_y, -abs_z)));
  compare("in testPhiCalc, phi in octant 6", phi * pi / 180., pi - expected_phi, 1.e-12);

  phi = s_calcPhi(sc_x, sc_z, astro::SkyDir(vec_t(-abs_x, -abs_y, -abs_z)));
  compare("in testPhiCalc, phi in octant 7", phi * pi / 180., pi + expected_phi, 1.e-12);

  phi = s_calcPhi(sc_x, sc_z, astro::SkyDir(vec_t(abs_x, -abs_y, -abs_z)));
  compare("in testPhiCalc, phi in octant 8", phi * pi / 180., 2 * pi - expected_phi, 1.e-12);
}

void RspGenTestApp::compare(const std::string & descriptor, const std::vector<double> & vec1, const std::vector<double> & vec2,
  const std::vector<double> & total) {
  if (vec1.size() != vec2.size() || vec1.size() != total.size()) {
    m_failed = true;
    std::cerr << "Unexpected: in test11, the three " << descriptor << " responses did not have same size." << std::endl;
  } else {
    for (std::vector<double>::size_type index = 0; index != total.size(); ++index) {
      double expected_total = vec1[index] + vec2[index];
      if (expected_total == total[index]) continue;
      else if (0. == expected_total && std::fabs(total[index]) < std::numeric_limits<double>::epsilon()) continue;
      else if (0. == total[index] && std::fabs(expected_total) < std::numeric_limits<double>::epsilon()) continue;
      else if (std::fabs((total[index] - expected_total)/ expected_total) < 10. * std::numeric_limits<double>::epsilon())
        continue;
      m_failed = true;
      std::cerr.precision(std::numeric_limits<double>::digits10);
      std::cerr << "Unexpected: in test11, " << descriptor << " gave total[" << index << "] = " << total[index] << ", not " <<
        expected_total << ", as expected." << std::endl;
    }
  }
}

void RspGenTestApp::compare(const std::string & descriptor, double value, double expected_value, double abs_tol) {
  int comparison = 0;
  if ((expected_value - value) > abs_tol) comparison = 1;
  else if ((value - expected_value) > abs_tol) comparison = -1;
  if (0 != comparison) {
    m_failed = true;
    std::cerr.precision(15);
    std::cerr << "Unexpected: " << descriptor << " is " << value << ", not " << expected_value << ", as expected" << std::endl;
  }
}

evtbin::Binner * RspGenTestApp::createStdBinner() {
  // Process bin definition file to produce true energy bins.
  std::auto_ptr<const tip::Table> true_en(tip::IFileSvc::instance().readTable(findFile("rspgen_energy_bins.fits"), "ENERGYBINS"));

  // Create object to hold the intervals in the bin definition file.
  evtbin::OrderedBinner::IntervalCont_t true_intervals(true_en->getNumRecords());

  // Read true energy bin definitions into interval object.
  int index = 0;
  for (tip::Table::ConstIterator itor = true_en->begin(); itor != true_en->end(); ++itor, ++index)
    true_intervals[index] = evtbin::Binner::Interval(s_MeV_per_keV*(*itor)["E_MIN"].get(), s_MeV_per_keV*(*itor)["E_MAX"].get());

  // Create a binner for true energy.
  return new evtbin::OrderedBinner(true_intervals);
}

irfInterface::Irfs * RspGenTestApp::createIrfs() const {
  return irfInterface::IrfsFactory::instance()->create("testIrfs::Front");
}

std::string RspGenTestApp::findFile(const std::string & file_root) const {
  return facilities::commonUtilities::joinPath(m_data_dir, file_root);
}

void RspGenTestApp::copyFile(const std::string & in_file, const std::string & out_file) const {
  std::ifstream in_stream(in_file.c_str(), std::ios::in | std::ios::binary);
  std::ofstream out_stream(out_file.c_str(), std::ios::out | std::ios::binary);
  char c;
  in_stream.get(c);
  while (in_stream) {
    out_stream.put(c);
    in_stream.get(c);
  }
}

// Factory object to create this test executable.
st_app::StAppFactory<RspGenTestApp> g_factory("gtrspgen");
