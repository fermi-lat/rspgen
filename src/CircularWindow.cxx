/** \file CircularWindow.cxx
    \brief Implementation of circular window class.
    \author James Peachey, HEASARC
*/
#include <stdexcept>

#include "rspgen/CircularWindow.h"

namespace rspgen {

  CircularWindow::CircularWindow(double radius): m_radius(radius) {
    if (radius <= 0) throw std::logic_error("CircularWindow() called for non-positive radius");
  }

  IWindow * CircularWindow::clone() const { return new CircularWindow(*this); }

  double CircularWindow::integrate(latResponse::IPsf * psf, double true_energy, double theta, double phi) const {
    // This case is very simple. The region is a circle of the given radius. The psf
    // object from latResponse is immediately capable of this integration.
    return psf->angularIntegral(true_energy, theta, phi, m_radius);
  }

}
