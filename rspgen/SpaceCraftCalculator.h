/** \file SpaceCraftCalculator.h
    \brief Declaration of SpaceCraftCalculator class.
    \author James Peachey, HEASARC
*/
#ifndef rspgen_SpaceCraftCalculator_h
#define rspgen_SpaceCraftCalculator_h

namespace astro {
  class SkyDir;
}

namespace evtbin {
  class Hist1D;
}

namespace irfInterface {
  class Irfs;
}

namespace st_stream {
  class OStream;
}

#include <string>
#include <vector>

namespace rspgen {

  class IWindow;

  /** \class SpaceCraftCalculator
      \brief 
  */
  class SpaceCraftCalculator {
    public:
      typedef std::vector<irfInterface::Irfs *> irf_cont_type;
      typedef std::vector<std::string> irf_name_cont_type;

      static void lookUpResponse(const std::string & resp, irf_name_cont_type & match);

      static void displayIrfNames(st_stream::OStream & os);

      SpaceCraftCalculator(const astro::SkyDir & src_dir, double theta_cut, double theta_bin_size, double psf_radius,
        const std::string & resp_type, const std::string & gti_file, const std::string & sc_file, const std::string & sc_table);

      virtual ~SpaceCraftCalculator();

      virtual double psf(double true_energy, double theta, double phi) const;

      double calcPhi(const astro::SkyDir & x_ref, const astro::SkyDir & z_ref, const astro::SkyDir & dir) const;

    private:
      evtbin::Hist1D * m_diff_exp;
      irf_cont_type m_irfs;
      IWindow * m_window;
      double m_total_exposure;
  };

}

#endif
