/** \file SpaceCraftCalculator.cxx
    \brief Interface for Point-specific response calculations.
    \author James Peachey, HEASARC
*/
#include "astro/SkyDir.h"

#include "evtbin/Gti.h"
#include "evtbin/Hist1D.h"
#include "evtbin/LinearBinner.h"

#include "irfInterface/IrfsFactory.h"

#include "rspgen/CircularWindow.h"
#include "rspgen/IWindow.h"
#include "rspgen/SpaceCraftCalculator.h"

#include "tip/Header.h"
#include "tip/IFileSvc.h"
#include "tip/LinearInterp.h"
#include "tip/Table.h"

#include <map>
#include <memory>
#include <stdexcept>

namespace rspgen {

  std::string SpaceCraftCalculator::lookUpResponse(const std::string & resp) {
    typedef std::map<std::string, std::string> Dict_t;
    static Dict_t s_resp_dict;

    /* For dictionary lookup, use all caps. */
    std::string uc_resp = resp;
    for (std::string::iterator itor = uc_resp.begin(); itor != uc_resp.end(); ++itor) *itor = toupper(*itor);

    if (s_resp_dict.empty()) {
      /* Use all caps for dictionary entries. */
      s_resp_dict["DC1F"] = "DC1::Front";
      s_resp_dict["DC1B"] = "DC1::Back";
      s_resp_dict["G25F"] = "Glast25::Front";
      s_resp_dict["G25B"] = "Glast25::Back";
      s_resp_dict["TESTF"] = "testIrfs::Front";
      s_resp_dict["TESTB"] = "testIrfs::Back";
    }
    Dict_t::const_iterator match = s_resp_dict.find(uc_resp);
    if (s_resp_dict.end() == match) return resp;
    return match->second;
  }

  SpaceCraftCalculator::SpaceCraftCalculator(const astro::SkyDir & src_dir, double theta_cut, double theta_bin_size,
    double psf_radius, const std::string & resp_type, const std::string & gti_file, const std::string & sc_file,
    const std::string & sc_table): m_diff_exp(0), m_irfs(0), m_window(0), m_total_exposure(0.) {
    using evtbin::Gti;

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

      // TODO: Read SCX, bin in phi as well as theta.
      double ra_scz = sc_record.get("RA_SCZ");
      double dec_scz = sc_record.get("DEC_SCZ");

      astro::SkyDir scz_pos(ra_scz, dec_scz);

      // Compute inclination angle from the source position to the sc.
      double theta = src_dir.difference(scz_pos) * 180. / M_PI;

      // Bin this angle into the histogram.
      diff_exp->fillBin(theta, delta_t);

      m_total_exposure += delta_t;
    }

    // Check that something was actually accumulated.
    if (0. == m_total_exposure)
      throw std::runtime_error("SpaceCraftCalculator cannot continue with 0. total exposure.");

    // Get irfs object.
    std::auto_ptr<irfInterface::Irfs> irfs(irfInterface::IrfsFactory::instance()->create(lookUpResponse(resp_type)));

    // Create window object for circular psf integration with the given inclination angle and psf radius.
    m_window = new CircularWindow(psf_radius);

    // Everything succeeded, so release the pointers from their auto_ptrs.
    m_diff_exp = diff_exp.release();
    m_irfs = irfs.release();
  }

  SpaceCraftCalculator::~SpaceCraftCalculator() {
    delete m_window;
    delete m_irfs;
    delete m_diff_exp;
  }

  double SpaceCraftCalculator::psf(double true_energy, double theta) const {
    // TODO: Integrate over phi.
    double phi = 0.;

    // Get the psf for this theta bin.
    double psf_val = m_window->integrate(m_irfs->psf(), true_energy, theta, phi);

    // Find the theta bin index corresponding to this theta.
    const evtbin::Binner * theta_bins = m_diff_exp->getBinners()[0];
    long bin_index = theta_bins->computeIndex(theta);

    // Return the psf for this theta bin, weighted by the fractional differential exposure.
    return psf_val * (*m_diff_exp)[bin_index] / m_total_exposure;
  }

}
