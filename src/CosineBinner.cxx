/** \file CosineBinner.cxx
    \brief Implementation CosineBinner class.
*/

#include <cmath>

#include "rspgen/CosineBinner.h"

namespace {
  double s_rad_per_deg = M_PI / 180.;
  double s_deg_per_rad = 180. / M_PI;
}

namespace rspgen {

  CosineBinner::CosineBinner(double interval_begin, double interval_end, double bin_size, const std::string & name):
    Binner(name),
    m_interval_begin(std::cos(interval_end * s_rad_per_deg)),
    m_interval_end(std::cos(interval_begin * s_rad_per_deg)),
    m_bin_size(bin_size),
    m_num_bins(long(std::ceil((m_interval_end - m_interval_begin)/bin_size))) {}

  long CosineBinner::computeIndex(double value) const {
    value = std::cos(value * s_rad_per_deg);
    if (value < m_interval_begin || value >= m_interval_end) return -1;
    return long((value - m_interval_begin) / m_bin_size);
  }

  long CosineBinner::getNumBins() const { return m_num_bins; }

  evtbin::Binner::Interval CosineBinner::getInterval(long index) const {
    // Check bounds, and handle endpoints explicitly to avoid any round-off:
    if (index < 0 || index >= m_num_bins)
      return Binner::Interval(0., 0.);
    else if (index == 0)
      return Binner::Interval(
        std::acos(m_interval_begin + /* (1 + index) * */ m_bin_size) * s_deg_per_rad,
        std::acos(m_interval_begin) * s_deg_per_rad
      );
    else if (index == m_num_bins - 1)
      return Binner::Interval(
        std::acos(m_interval_end) * s_deg_per_rad,
        std::acos(m_interval_begin + index * m_bin_size) * s_deg_per_rad
      );

    return Binner::Interval(
      std::acos(m_interval_begin + (1 + index) * m_bin_size) * s_deg_per_rad,
      std::acos(m_interval_begin + index * m_bin_size) * s_deg_per_rad
    );
  }

  evtbin::Binner * CosineBinner::clone() const { return new CosineBinner(*this); }

}
