/** \file PointResponse.cxx
    \brief Interface for Point-specific response calculations.
    \author James Peachey, HEASARC
*/
#include <iostream>
#include <memory>
#include <stdexcept>

#include "astro/SkyDir.h"

#include "evtbin/Gti.h"
#include "evtbin/LinearBinner.h"

#include "rspgen/CircularWindow.h"
#include "rspgen/CosineBinner.h"
#include "rspgen/PointResponse.h"

#include "tip/IFileSvc.h"
#include "tip/LinearInterp.h"
#include "tip/Table.h"
#include "tip/TipException.h"

namespace rspgen {

  PointResponse::PointResponse(double ps_ra, double ps_dec, double theta_cut, double theta_bin_size, double psf_radius,
    long phi_num_bins, const std::string & resp_type, const std::string & spec_file, const std::string & sc_file,
    const std::string & sc_table, const evtbin::Binner * true_en_binner): IResponse(resp_type, spec_file, true_en_binner),
    m_window(0), m_diff_exp(0), m_total_exposure(0.) {
    // TODO: This constructor is now duplicated in class SpaceCraftCalculator. PointResponse should have a member
    // of SpaceCraftCalculator and not duplicate this code here.
    using evtbin::Gti;

    // Process spacecraft data.
    std::auto_ptr<const tip::Table> table(tip::IFileSvc::instance().readTable(sc_file, sc_table));

    // Put point source direction into standard form.
    astro::SkyDir ps_pos(ps_ra, ps_dec);

    double phi_cut = 360.;
    double phi_bin_size = phi_cut / phi_num_bins;

    // Set up a histogram to hold the binned differential exposure (theta vs. DeltaT).
    std::auto_ptr<evtbin::Hist2D> diff_exp(new evtbin::Hist2D(
      CosineBinner(0., theta_cut, theta_bin_size), evtbin::LinearBinner(0., phi_cut, phi_bin_size)
    ));

    // Get GTI information.
    Gti gti(spec_file);

    // Start with first GTI in the GTI table.
    Gti::ConstIterator gti_pos = gti.begin();

    // Read SC Z positions, bin them into a histogram:
    for (tip::Table::ConstIterator itor = table->begin(); itor != table->end() && gti_pos != gti.end(); ++itor) {
      double start = (*itor)["START"].get();
      double stop = (*itor)["STOP"].get();

      double fract = gti.getFraction(start, stop, gti_pos);
      // Save some time by not computing further if fraction is 0.
      if (0. == fract) continue;

      // Get size of interval, multiply by the fraction of the time which overlapped the GTI.
      double delta_t = fract * (*itor)["LIVETIME"].get();

      // Get object for interpolating values from the table.
      tip::LinearInterp sc_record(itor, table->end());

      // Get interpolated SC coordinates.
      sc_record.interpolate("START", (start + stop) / 2.);

      double ra_scz = sc_record.get("RA_SCZ");
      double dec_scz = sc_record.get("DEC_SCZ");

      astro::SkyDir scz_pos(ra_scz, dec_scz);

      // Compute inclination angle from the source position relative to the sc z.
      double theta = ps_pos.difference(scz_pos) * 180. / M_PI;

      double ra_scx = sc_record.get("RA_SCX");
      double dec_scx = sc_record.get("DEC_SCX");

      // Compute azimuthal angle from the source position relative to the sc z and x.
      astro::SkyDir scx_pos(ra_scx, dec_scx);
      double phi = calcPhi(scx_pos, scz_pos, ps_pos);

      // Bin this angle into the histogram.
      diff_exp->fillBin(theta, phi, delta_t);

      m_total_exposure += delta_t;
    }

    // Check that something was actually accumulated.
    if (0. == m_total_exposure) throw std::runtime_error("PointResponse constructor: cannot continue with 0. total exposure.");

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
    // TODO: This method uses functionality which ought to be in class SpaceCraftCalculator. PointResponse should use
    // a member of SpaceCraftCalculator and thus simplify the code here.
    // Iterate over the bins of the histogram.
    const evtbin::Binner * theta_bins = m_diff_exp->getBinners()[0];
    long num_theta_bins = theta_bins->getNumBins();
    const evtbin::Binner * phi_bins = m_diff_exp->getBinners()[1];
    long num_phi_bins = phi_bins->getNumBins();

    // Reset the response vector to be all zeroes, and enough of them.
    response.clear();
    response.assign(m_app_en_binner->getNumBins(), 0.);

    // Integrate over binned angle-dependent differential exposure.
    for (long theta_index = 0; theta_index < num_theta_bins; ++theta_index) {
      // Get the angle from the center of the bin.
      double theta = theta_bins->getInterval(theta_index).midpoint();
      for (long phi_index = 0; phi_index < num_phi_bins; ++phi_index) {
        double phi = phi_bins->getInterval(phi_index).midpoint();
  
        for (irf_cont_type::iterator itor = m_irfs.begin(); itor != m_irfs.end(); ++itor) {
          irfInterface::Irfs * irfs = *itor;
  
          // Compute effective area, which is a function of true_energy and sc pointing direction only.
          double aeff_val = irfs->aeff()->value(true_energy, theta, phi);
          // Use the window object to integrate psf over the region.
          double int_psf_val = m_window->integrate(irfs->psf(), true_energy, theta, phi);
  
          // For each apparent energy bin, compute integral of the redistribution coefficient.
          for (long index = 0; index < m_app_en_binner->getNumBins(); ++index) {
            // Get limits of integration over apparent energy bins.
            evtbin::Binner::Interval limits = m_app_en_binner->getInterval(index);
  
            response[index] += (*m_diff_exp)[theta_index][phi_index] / m_total_exposure * aeff_val *
              irfs->edisp()->integral(limits.begin(), limits.end(), true_energy, theta, phi) * int_psf_val;
          }
        }
      }
    }
  }

  double PointResponse::psf(double true_energy, double theta, double phi) const {
    // TODO: This method is now duplicated in class SpaceCraftCalculator. PointResponse should use a member
    // of SpaceCraftCalculator and not duplicate this code here.

    double psf_val = 0.;

    for (irf_cont_type::const_iterator itor = m_irfs.begin(); itor != m_irfs.end(); ++itor) {
      irfInterface::Irfs * irfs = *itor;

      // Get the psf for this theta bin.
      psf_val += m_window->integrate(irfs->psf(), true_energy, theta, phi);
    }

    // Find the theta bin index corresponding to this theta.
    const evtbin::Binner * theta_bins = m_diff_exp->getBinners()[0];
    long theta_index = theta_bins->computeIndex(theta);

    // Find the phi bin index corresponding to this phi.
    const evtbin::Binner * phi_bins = m_diff_exp->getBinners()[1];
    long phi_index = phi_bins->computeIndex(phi);

    // Return the psf for this theta bin, weighted by the fractional differential exposure.
    return psf_val * (*m_diff_exp)[theta_index][phi_index] / m_total_exposure;
  }

  std::pair<long, long> PointResponse::getSpatialNumBins() const {
    const evtbin::Binner * theta_bins = m_diff_exp->getBinners()[0];
    long num_theta_bins = theta_bins->getNumBins();
    const evtbin::Binner * phi_bins = m_diff_exp->getBinners()[1];
    long num_phi_bins = phi_bins->getNumBins();
    return std::make_pair(num_theta_bins, num_phi_bins);
  }

}
