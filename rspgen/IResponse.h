/** \file IResponse.h
    \brief Generic interface for response calculators.
    \author James Peachey, HEASARC
*/
#ifndef rspgen_IResponse_h
#define rspgen_IResponse_h

#include <string>
#include <vector>

#include "evtbin/Binner.h"

#include "irfInterface/Irfs.h"

#include "rspgen/IResponse.h"

#include "tip/Header.h"

namespace rspgen {

  /** \class IResponse
      \brief Generic interface for response calculators.
  */
  class IResponse {
    public:
      /** \brief Create GrbResponse object for a given burst and spacecraft coordinates
          \param resp_type Identifies response function type.
          \param spec_file The name of the spectrum file.
          \param true_en_binner Binner object used for true energy bin definitions.
      */
      IResponse(const std::string resp_type, const std::string & spec_file, const evtbin::Binner * true_en_binner);

      virtual ~IResponse() throw();

      /** \brief Compute responses and write them to eh output file.
          \param creator String to write for the CREATOR keyword.
          \param file_name Name of the output file.
          \param fits_template Name of the template file used to create the output file.
      */
      virtual void writeOutput(const std::string & creator, const std::string & file_name, const std::string & fits_template);

      /** \brief Compute the response for the given value of true energy.
          \param true_energy The energy for which to compute response.
          \param response The response (vector of apparent energy bins).
      */
      virtual void compute(double true_energy, std::vector<double> & response) = 0;

    protected:
      static const double s_keV_per_MeV;
      static const double s_MeV_per_keV;

      /** \brief Look up response nicknames in dictionary. If not found, just return original string.
          \param resp The response nickname.
      */
      static std::string lookUpResponse(const std::string & resp);

      IResponse();

      tip::Header::KeyValCont_t m_kwds;
      evtbin::Binner * m_true_en_binner;
      evtbin::Binner * m_app_en_binner;
      irfInterface::Irfs * m_irfs;
  };

}

#endif
