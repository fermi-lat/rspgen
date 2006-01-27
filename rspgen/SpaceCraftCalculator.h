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

#include <string>

namespace rspgen {

  class IWindow;

  /** \class SpaceCraftCalculator
      \brief 
  */
  class SpaceCraftCalculator {
    public:
      static std::string lookUpResponse(const std::string & resp);

      SpaceCraftCalculator(const astro::SkyDir & src_dir, double theta_cut, double theta_bin_size, double psf_radius,
        const std::string & resp_type, const std::string & gti_file, const std::string & sc_file, const std::string & sc_table);

      virtual ~SpaceCraftCalculator();

      virtual double psf(double true_energy, double theta) const;

    private:
      evtbin::Hist1D * m_diff_exp;
      irfInterface::Irfs * m_irfs;
      IWindow * m_window;
      double m_total_exposure;
  };

}

#endif
