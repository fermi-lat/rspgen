/** \file PointResponse.h
    \brief Interface for Point-specific response calculations.
    \author James Peachey, HEASARC
*/
#ifndef rspgen_PointResponse_h
#define rspgen_PointResponse_h

#include <string>
#include <vector>

#include "evtbin/Binner.h"
#include "evtbin/Hist1D.h"

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
          \param resp_type Identifies response function type.
          \param spec_file The name of the spectrum file.
          \param sc_file The name of the file containing spacecraft data.
          \param sc_table The name of the table containing spacecraft data.
          \param true_en_binner Binner object used for true energy bin definitions.
      */
      PointResponse(double ps_ra, double ps_dec, double theta_cut, double theta_bin_size, double psf_radius,
        const std::string resp_type, const std::string & spec_file, const std::string & sc_file,
        const std::string & sc_table, const evtbin::Binner * true_en_binner);

      virtual ~PointResponse() throw();

      /** \brief Compute the response for the given value of true energy.
          \param true_energy The energy for which to compute response.
          \param response The response (vector of apparent energy bins).
      */
      virtual void compute(double true_energy, std::vector<double> & response);

    private:
      IWindow * m_window;
      evtbin::Hist1D * m_diff_exp;
      double m_total_exposure;
  };

}

#endif
