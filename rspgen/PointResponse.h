/** \file PointResponse.h
    \brief Interface for Point-specific response calculations.
    \author James Peachey, HEASARC
*/
#ifndef rspgen_PointResponse_h
#define rspgen_PointResponse_h

#include <string>
#include <utility>
#include <vector>

#include "evtbin/Binner.h"
#include "evtbin/Hist2D.h"

#include "rspgen/IResponse.h"
#include "rspgen/IWindow.h"

#include "tip/Header.h"

namespace rspgen {

  /** \class PointResponse
      \brief Interface for Point-specific response calculations.
  */

  class PointResponse : public IResponse {
    public:
      /** \brief Create PointResponse object for a given burst and spacecraft coordinates
          \param ps_ra The RA of the point source, in degrees.
          \param ps_dec The DEC of the .point source, in degrees.
          \param theta_cut The cutoff for spacecraft theta integration, in degrees.
          \param theta_bin_size The size of bins for spacecraft theta integration, in degrees.
          \param psf_radius The radius of the psf integration, in degrees.
          \param phi_num_bins The number of bins used in binning the azimuthal angle phi.
          \param resp_type Identifies response function type.
          \param spec_file The name of the spectrum file.
          \param sc_file The name of the file containing spacecraft data.
          \param sc_table The name of the table containing spacecraft data.
          \param true_en_binner Binner object used for true energy bin definitions.
      */
      PointResponse(double ps_ra, double ps_dec, double theta_cut, double theta_bin_size, double psf_radius,
        long phi_num_bins, const std::string & resp_type, const std::string & spec_file, const std::string & sc_file,
        const std::string & sc_table, const evtbin::Binner * true_en_binner);

      virtual ~PointResponse() throw();

      /** \brief Compute the response for the given value of true energy.
          \param true_energy The energy for which to compute response.
          \param response The response (vector of apparent energy bins).
      */
      virtual void compute(double true_energy, std::vector<double> & response);

      /** \brief Compute the psf for the given energy and inclination angle,
          weighted by the fraction of the livetime during which the spacecraft axis was
          oriented at this angle.
          \param true_energy The energy for which to compute the psf.
          \param theta Spacecraft inclination angle in degrees.
          \param phi Spacecraft azimuthal angle in degrees.
      */
      virtual double psf(double true_energy, double theta, double phi) const;

      /// \brief Return the number of bins used in binning in theta, phi.
      std::pair<long, long> getSpatialNumBins() const;

    private:
      IWindow * m_window;
      evtbin::Hist2D * m_diff_exp;
      double m_total_exposure;
  };

}

#endif
