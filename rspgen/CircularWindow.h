/** \file CircularWindow.h
    \brief Declaration of circular region windowing class.
    \author James Peachey, HEASARC
*/
#ifndef rspgen_CircularWindow_h
#define rspgen_CircularWindow_h

#include "rspgen/IWindow.h"

namespace rspgen {

  /** \class CircularWindow
      \brief Circular region windowing class.
  */
  class CircularWindow : public IWindow {
    public:
      /** \brief Create a circular window (region) of the given radius.
          \param radius The radius of the circular region in degrees.
      */
      CircularWindow(double radius);

      /// \brief Make a copy of this object.
      virtual IWindow * clone() const;

      /** \brief Integrate the psf over this window for the given values.
          \param psf The psf object obtained from the latResponse package.
          \param true_energy The true energy for which to evaluate the integral.
          \param theta The true inclination angle of the incident photon.
          \param phi The true asimuthal angle of the incident photon.
      */
      virtual double integrate(irfInterface::IPsf * psf, double true_energy, double theta, double phi = 0.) const;

    private:
      double m_radius;
  };

}

#endif
