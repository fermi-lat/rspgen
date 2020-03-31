/** \file GrbResponse.cxx
    \brief Interface for Grb-specific response calculations.
    \author James Peachey, HEASARC
*/
#include <memory>
#include <stdexcept>

#include "astro/SkyDir.h"

#include "rspgen/CircularWindow.h"
#include "rspgen/GrbResponse.h"

#include "tip/IFileSvc.h"
#include "tip/LinearInterp.h"
#include "tip/Table.h"

namespace rspgen {

  GrbResponse::GrbResponse(double theta, double phi, const evtbin::Binner * true_en_binner, const evtbin::Binner * app_en_binner,
    irfInterface::Irfs * irfs, const IWindow * window): m_theta(theta), m_phi(phi), m_window(0) {
    // Check inputs.
    if (0 == true_en_binner || 0 == app_en_binner || 0 == irfs || 0 == window)
      throw std::logic_error("GrbResponse constructor was passed a null pointer");

    // Clone inputs to become members.
    m_true_en_binner = true_en_binner->clone();
    m_app_en_binner = app_en_binner->clone();
    m_irfs.push_back(irfs->clone());
    m_window = window->clone();
  }

  GrbResponse::GrbResponse(double grb_ra, double grb_dec, double grb_time, double psf_radius, const std::string & resp_type,
    const std::string & spec_file, const std::string & sc_file, const std::string & sc_table,
    const evtbin::Binner * true_en_binner): IResponse(resp_type, spec_file, true_en_binner), m_theta(0.), m_phi(0.), m_window(0) {
    // Process spacecraft data.
    std::auto_ptr<const tip::Table> table(tip::IFileSvc::instance().readTable(sc_file, sc_table));

    // Get object for interpolating values from the table.
    tip::LinearInterp sc_record(table->begin(), table->end());

    // Interpolate values of spacecraft parameters for the burst time.
    sc_record.interpolate("START", grb_time);

    // Compute angles from burst RA and DEC and spacecraft X and Z axes.
    astro::SkyDir z_ref(sc_record.get("RA_SCZ"), sc_record.get("DEC_SCZ"));
    astro::SkyDir x_ref(sc_record.get("RA_SCX"), sc_record.get("DEC_SCX"));
    astro::SkyDir grb_pos(grb_ra, grb_dec);
    m_theta = z_ref.difference(grb_pos)*180./M_PI;
    m_phi = calcPhi(x_ref, z_ref, grb_pos);

    // Create window object for circular psf integration with the given inclination angle and psf radius.
    m_window = new CircularWindow(psf_radius);
  }

  GrbResponse::~GrbResponse() throw() {
    delete m_window;
  }

  void GrbResponse::compute(double true_energy, std::vector<double> & response) {
    // Make sure the response vector has enough room, but all 0s.
    response.clear();
    response.resize(m_app_en_binner->getNumBins(), 0.);

    for (irf_cont_type::iterator itor = m_irfs.begin(); itor != m_irfs.end(); ++itor) {
      irfInterface::Irfs * irfs = *itor;

      // Compute effective area, which is a function of true_energy and sc pointing direction only.
      double aeff_val = irfs->aeff()->value(true_energy, m_theta, m_phi);

      // Use the window object to integrate psf over the region.
      double int_psf_val = m_window->integrate(irfs->psf(), true_energy, m_theta, m_phi);

      // For each apparent energy bin, compute integral of the redistribution coefficient.
      for (long index = 0; index < m_app_en_binner->getNumBins(); ++index) {
        // Get limits of integration over apparent energy bins
        evtbin::Binner::Interval limits = m_app_en_binner->getInterval(index);

        // Compute the response for the current app. energy bin.
        response[index] += aeff_val * irfs->edisp()->integral(limits.begin(), limits.end(), true_energy, m_theta, m_phi) *
          int_psf_val;
      }
    }

  }

}
