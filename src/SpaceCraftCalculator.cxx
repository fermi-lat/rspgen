/** \file SpaceCraftCalculator.cxx
    \brief Interface for Point-specific response calculations.
    \author James Peachey, HEASARC
*/
#include "astro/SkyDir.h"

#include "evtbin/Gti.h"
#include "evtbin/Hist1D.h"
#include "evtbin/LinearBinner.h"

#include "irfInterface/IrfsFactory.h"
#include "irfLoader/Loader.h"

#include "rspgen/CircularWindow.h"
#include "rspgen/IWindow.h"
#include "rspgen/SpaceCraftCalculator.h"

#include "st_stream/Stream.h"

#include "tip/Header.h"
#include "tip/IFileSvc.h"
#include "tip/LinearInterp.h"
#include "tip/Table.h"

#include "CLHEP/Vector/ThreeVector.h"

#include <algorithm>
#include <cctype>
#include <map>
#include <memory>
#include <stdexcept>

namespace {
  class MatchCaseInsensitive {
    public:
      MatchCaseInsensitive(const std::string & match_string): m_string_to_match(match_string) {
        for (std::string::iterator itor = m_string_to_match.begin(); itor != m_string_to_match.end(); ++itor)
          *itor = toupper(*itor);
      }

      bool operator ()(const std::string & input) const {
        std::string::const_iterator itor1 = input.begin();
        std::string::const_iterator itor2 = m_string_to_match.begin();
        for (; itor1 != input.end() && itor2 != m_string_to_match.end() && toupper(*itor1) == *itor2; ++itor1, ++itor2);
        return itor1 == input.end() && itor2 == m_string_to_match.end();
      }
    private:
      std::string m_string_to_match;
  };

  template <typename T>
  class MatchKey{
    public:
      MatchKey(const MatchCaseInsensitive & matcher): m_matcher(matcher) {}
      bool operator ()(const T & input) const {
        return m_matcher(input.first);
      }
    private:
      const MatchCaseInsensitive & m_matcher;
  };
}

namespace rspgen {

  void SpaceCraftCalculator::lookUpResponse(const std::string & resp, std::vector<std::string> & match) {
    // Start off with no matching responses.
    match.clear();

    // Get look-up container of irf containers.
    typedef std::map<std::string, std::vector<std::string> > irf_lookup_type;
    const irf_lookup_type & resp_id(irfLoader::Loader::respIds());

    // Create object for matching to the input name.
    MatchCaseInsensitive matcher(resp);

    MatchKey<irf_lookup_type::value_type> match_key(matcher);

    irf_lookup_type::const_iterator found_cont = std::find_if(resp_id.begin(), resp_id.end(), match_key);

    // Look for upper-cased response descriptor.
    if (resp_id.end() != found_cont) {
      // Found a container of responses matching the descriptor, so assign them to the output.
      match = found_cont->second;
    } else {
      // Didn't find container, so see if the given name matches one of the irf names directly.
      irf_name_cont_type irf_names;
      irfInterface::IrfsFactory::instance()->getIrfsNames(irf_names);
      irf_name_cont_type::iterator found_string = std::find_if(irf_names.begin(), irf_names.end(), matcher);
      if (irf_names.end() != found_string) {
        match.push_back(*found_string);
      }
    }
  }

  void SpaceCraftCalculator::displayIrfNames(st_stream::OStream & os) {
    // Get look-up container of irf containers.
    typedef std::map<std::string, std::vector<std::string> > irf_lookup_type;
    const irf_lookup_type & resp_id(irfLoader::Loader::respIds());

    os.prefix() << "Valid instrument response function names are:" << std::endl;

    for (irf_lookup_type::const_iterator cont_itor = resp_id.begin(); cont_itor != resp_id.end(); ++cont_itor) {
      os.prefix() << "\t" << cont_itor->first;
      std::string connector = " = ";
      for (irf_name_cont_type::const_iterator itor = cont_itor->second.begin(); itor != cont_itor->second.end(); ++itor) {
        os << connector << *itor;
        connector = " + ";
      }
      os << std::endl;
    }

    // Get container of irf atomic names.
    irf_name_cont_type irf_simple_name;
    irfInterface::IrfsFactory::instance()->getIrfsNames(irf_simple_name);

    for (irf_name_cont_type::const_iterator itor = irf_simple_name.begin(); itor != irf_simple_name.end(); ++itor) {
      os.prefix() << "\t" << *itor << std::endl;
    }
  }

  SpaceCraftCalculator::SpaceCraftCalculator(const astro::SkyDir & src_dir, double theta_cut, double theta_bin_size,
    double psf_radius, const std::string & resp_type, const std::string & gti_file, const std::string & sc_file,
    const std::string & sc_table): m_diff_exp(0), m_irfs(), m_window(0), m_total_exposure(0.) {
    using evtbin::Gti;
    // TODO If this class is ever actually used to merge common functions from GrbReponse and PointResponse,
    // phi needs to be handled properly, including Hist2D etc. See PointResponse.

    // Get spacecraft data.
    std::auto_ptr<const tip::Table> table(tip::IFileSvc::instance().readTable(sc_file, sc_table));

    // Get GTI information.
    Gti gti(gti_file);

    // Set up a histogram to hold the binned differential exposure (theta vs. DeltaT).
    std::auto_ptr<evtbin::Hist1D> diff_exp(new evtbin::Hist1D(evtbin::LinearBinner(0., theta_cut, theta_bin_size)));

    // Start with first GTI in the GTI table.
    Gti::ConstIterator gti_pos = gti.begin();

    // Read SC positions, bin them into a histogram:
    for (tip::Table::ConstIterator itor = table->begin(); itor != table->end(); ++itor) {
      double start = (*itor)["START"].get();
      double stop = (*itor)["STOP"].get();

      double fract = gti.getFraction(start, stop, gti_pos);
      // Save some time by not computing further if fraction is 0.
      if (0. == fract) continue;

      // If we fell off the edge of the last GTI, no point in continuing this loop.
      if (gti.end() == gti_pos) break;

      // Get size of interval, multiply by the fraction of the time which overlapped the GTI.
      double delta_t = fract * (*itor)["LIVETIME"].get();
    
      // Get object for interpolating values from the table.
      tip::LinearInterp sc_record(itor, table->end());

      // Get interpolated SC coordinates.
      sc_record.interpolate("START", (start + stop) / 2.);

      // Spacecraft Z.
      double ra_scz = sc_record.get("RA_SCZ");
      double dec_scz = sc_record.get("DEC_SCZ");
      astro::SkyDir scz_pos(ra_scz, dec_scz);

      // Compute inclination angle from the source position to the sc.
      double theta = src_dir.difference(scz_pos) * 180. / M_PI;

      // Spacecraft X.
      double ra_scx = sc_record.get("RA_SCX");
      double dec_scx = sc_record.get("DEC_SCX");
      astro::SkyDir scx_pos(ra_scx, dec_scx);

      // Bin this angle into the histogram.
      diff_exp->fillBin(theta, delta_t);

      m_total_exposure += delta_t;
    }

    // Check that something was actually accumulated.
    if (0. == m_total_exposure)
      throw std::runtime_error("SpaceCraftCalculator cannot continue with 0. total exposure.");

    // Get container of irf names matching the given response name.
    irf_name_cont_type irf_name;
    lookUpResponse(resp_type, irf_name);
    if (irf_name.empty()) throw std::runtime_error("SpaceCraftCalculator cannot find response matching " + resp_type);

    // Populate the irfs container with matching irfs.
    m_irfs.resize(irf_name.size(), 0);
    for (irf_name_cont_type::size_type index = 0; index != irf_name.size(); ++index) {
      m_irfs[index] = irfInterface::IrfsFactory::instance()->create(irf_name[index]);
    }

    // Create window object for circular psf integration with the given inclination angle and psf radius.
    m_window = new CircularWindow(psf_radius);

    // Everything succeeded, so release the pointers from their auto_ptrs.
    m_diff_exp = diff_exp.release();
  }

  SpaceCraftCalculator::~SpaceCraftCalculator() {
    delete m_window;
    for (irf_cont_type::reverse_iterator itor = m_irfs.rbegin(); itor != m_irfs.rend(); ++itor) delete *itor;
    delete m_diff_exp;
  }

  double SpaceCraftCalculator::psf(double true_energy, double theta, double phi) const {
    double psf_val = 0.;

    for (irf_cont_type::const_iterator itor = m_irfs.begin(); itor != m_irfs.end(); ++itor) {
      irfInterface::Irfs * irfs = *itor;

      // Get the psf for this theta bin.
      psf_val += m_window->integrate(irfs->psf(), true_energy, theta, phi);
    }

    // Find the theta bin index corresponding to this theta.
    const evtbin::Binner * theta_bins = m_diff_exp->getBinners()[0];
    long bin_index = theta_bins->computeIndex(theta);

    // Return the psf for this theta bin, weighted by the fractional differential exposure.
    return psf_val * (*m_diff_exp)[bin_index] / m_total_exposure;
  }

  double SpaceCraftCalculator::calcPhi(const astro::SkyDir & x_ref, const astro::SkyDir & z_ref, const astro::SkyDir & dir) const {
    typedef CLHEP::Hep3Vector vec_t;
    static const double pi = std::acos(-1);
    const vec_t & x_hat = x_ref.dir();
    const vec_t y_hat = z_ref.dir().cross(x_hat);
    double phi = std::atan2(dir.dir().dot(y_hat), dir.dir().dot(x_hat));
    while (phi < 0.) phi += 2 * pi;
    return phi * 180./pi;
  }

}
