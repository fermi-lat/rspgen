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

  GrbResponse::GrbResponse(double theta, const evtbin::Binner * true_en_binner, const evtbin::Binner * app_en_binner,
    latResponse::Irfs * irfs, const IWindow * window): m_theta(theta), m_window(0) {
    // Check inputs.
    if (0 == true_en_binner || 0 == app_en_binner || 0 == irfs || 0 == window)
      throw std::logic_error("GrbResponse constructor was passed a null pointer");

    // Clone inputs to become members.
    m_true_en_binner = true_en_binner->clone();
    m_app_en_binner = app_en_binner->clone();
    m_irfs= irfs->clone();
    m_window = window->clone();
  }

  GrbResponse::GrbResponse(double grb_ra, double grb_dec, double grb_time, double psf_radius, const std::string resp_type,
    const std::string & spec_file, const std::string & sc_file, const evtbin::Binner * true_en_binner):
    IResponse(resp_type, spec_file, true_en_binner), m_theta(0.), m_window(0) {
    // Process spacecraft data.
    std::auto_ptr<const tip::Table> sc_table(tip::IFileSvc::instance().readTable(sc_file, "Ext1"));

    // Get object for interpolating values from the table.
    tip::LinearInterp sc_record(sc_table->begin(), sc_table->end());

    // Interpolate values of spacecraft parameters for the burst time.
    sc_record.interpolate("START", grb_time);

    // Compute inclination angle from burst RA and DEC and spacecraft Z RA and DEC.
    m_theta = astro::SkyDir(grb_ra, grb_dec).difference(astro::SkyDir(sc_record.get("RA_SCZ"), sc_record.get("DEC_SCZ")))*180./M_PI;


    // Create window object for circular psf integration with the given inclination angle and psf radius.
    m_window = new CircularWindow(psf_radius);
  }

  GrbResponse::~GrbResponse() throw() {
    delete m_window;
  }

  void GrbResponse::compute(double true_energy, std::vector<double> & response) {
    double phi = 0.; // Dummy for now.

    // Compute effective area, which is a function of true_energy and sc pointing direction only.
    double aeff_val = m_irfs->aeff()->value(true_energy, m_theta, phi);

    // Use the window object to integrate psf over the region.
    double int_psf_val = m_window->integrate(m_irfs->psf(), true_energy, m_theta, phi);
  
    // Make sure the response vector has enough room.
    response.resize(m_app_en_binner->getNumBins());

    // For each apparent energy bin, compute integral of the redistribution coefficient.
    for (long index = 0; index < m_app_en_binner->getNumBins(); ++index) {
      // Get limits of integration over apparent energy bins
      evtbin::Binner::Interval limits = m_app_en_binner->getInterval(index);

      // Compute the response for the current app. energy bin.
      response[index] = aeff_val * m_irfs->edisp()->integral(limits.begin(), limits.end(), true_energy, m_theta, phi) * int_psf_val;
    }
  
  }

}
