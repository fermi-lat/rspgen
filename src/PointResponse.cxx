/** \file PointResponse.cxx
    \brief Interface for Point-specific response calculations.
    \author James Peachey, HEASARC
*/
#include <iostream>
#include <memory>
#include <stdexcept>

#include "astro/SkyDir.h"

#include "evtbin/LinearBinner.h"

#include "rspgen/CircularWindow.h"
#include "rspgen/Gti.h"
#include "rspgen/PointResponse.h"

#include "tip/IFileSvc.h"
#include "tip/LinearInterp.h"
#include "tip/Table.h"
#include "tip/TipException.h"

namespace rspgen {

  PointResponse::PointResponse(double ps_ra, double ps_dec, double theta_cut, double theta_bin_size, double psf_radius,
    const std::string resp_type, const std::string & spec_file, const std::string & sc_file, const evtbin::Binner * true_en_binner):
    IResponse(resp_type, spec_file, true_en_binner), m_window(0), m_diff_exp(0), m_total_exposure(0.) {
    // Process spacecraft data.
    std::auto_ptr<const tip::Table> sc_table(tip::IFileSvc::instance().readTable(sc_file, "Ext1"));

    // Put point source direction into standard form.
    astro::SkyDir ps_pos(ps_ra, ps_dec);

    // Set up a histogram to hold the binned differential exposure (theta vs. DeltaT).
    std::auto_ptr<evtbin::Hist1D> diff_exp(new evtbin::Hist1D(evtbin::LinearBinner(0., theta_cut, theta_bin_size)));

    // Get GTI information.
    Gti gti(spec_file);

    // Start with first GTI in the GTI table.
    Gti::ConstIterator gti_pos = gti.begin();

    // Read SC Z positions, bin them into a histogram:
    for (tip::Table::ConstIterator itor = sc_table->begin(); itor != sc_table->end(); ++itor) {
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
      tip::LinearInterp sc_record(itor, sc_table->end());

      // Get interpolated SC coordinates.
      sc_record.interpolate("START", (start + stop) / 2.);

      double ra_scz = sc_record.get("RA_SCZ");
      double dec_scz = sc_record.get("DEC_SCZ");

      astro::SkyDir scz_pos(ra_scz, dec_scz);

      // Compute inclination angle from the source position to the sc.
      double theta = ps_pos.difference(scz_pos) * 180. / M_PI;

      // Bin this angle into the histogram.
      diff_exp->fillBin(theta, delta_t);

      m_total_exposure += delta_t;
    }

    // Check that something was actually accumulated.
    if (0. == m_total_exposure) throw std::runtime_error("PointReseponse constructor: cannot continue with 0. total exposure.");

    // Create window object for circular psf integration with the given inclination angle and psf radius.
    m_window = new CircularWindow(psf_radius);

    // At this point, object was correctly constructed, so release the pointer from the auto_ptr.
    m_diff_exp = diff_exp.release();
  }

  PointResponse::~PointResponse() throw() {
    delete m_window;
    delete m_diff_exp;
  }

  void PointResponse::compute(double true_energy, std::vector<double> & response) {
    double phi = 0.; // Dummy for now.

    // Iterate over the bins of the histogram.
    const evtbin::Binner * theta_bins = m_diff_exp->getBinners()[0];
    long num_theta_bins = theta_bins->getNumBins();

    // Reset the response vector to be all zeroes, and enough of them.
    response.assign(m_app_en_binner->getNumBins(), 0.);

    // Integrate over binned angle-dependent differential exposure.
    for (long bin_index = 0; bin_index < num_theta_bins; ++bin_index) {
      // Get the angle from the center of the bin.
      double theta = theta_bins->getInterval(bin_index).midpoint();

      // Compute effective area, which is a function of true_energy and sc pointing direction only.
      double aeff_val = m_irfs->aeff()->value(true_energy, theta, phi);

      // Use the window object to integrate psf over the region.
      double int_psf_val = m_window->integrate(m_irfs->psf(), true_energy, theta, phi);
  
      // For each apparent energy bin, compute integral of the redistribution coefficient.
      for (long index = 0; index < m_app_en_binner->getNumBins(); ++index) {
        // Get limits of integration over apparent energy bins
        evtbin::Binner::Interval limits = m_app_en_binner->getInterval(index);
    
        response[index] += (*m_diff_exp)[bin_index] / m_total_exposure * aeff_val *
          m_irfs->edisp()->integral(limits.begin(), limits.end(), true_energy, theta, phi) * int_psf_val;
      }
    }
  }

}
