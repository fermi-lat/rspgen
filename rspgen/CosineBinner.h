/** \file CosineBinner.h
    \brief Declaration of CosineBinner class.
*/
#ifndef rspgen_CosineBinner_h
#define rspgen_CosineBinner_h

#include <string>

#include "evtbin/Binner.h"

namespace rspgen {
  /** \class CosineBinner
      \brief Declaration of a binner which is linearly uniform in cosine space.
  */
  class CosineBinner : public evtbin::Binner {
    public:
      /** \brief Construct the binner object.
          \param interval_begin Left boundary of the binning interval.
          \param interval_end Right boundary of the binning interval.
          \param bin_size The size of the bins.
          \param name Optional name of the quantity being binned.
      */
      CosineBinner(double interval_begin, double interval_end, double bin_size, const std::string & name = std::string());

      /** \brief Return the bin number for the given value.
          \param value The value being binned.
      */
      virtual long computeIndex(double value) const;

      /** \brief Return the number of bins currently defined.
      */
      virtual long getNumBins() const;

      /** \brief Return the interval spanned by the given bin.
          \param index The index indicating the bin number.
      */
      virtual evtbin::Binner::Interval getInterval(long index) const;

      /** \brief Create copy of this object.
      */
      virtual evtbin::Binner * clone() const;

    private:
      double m_interval_begin;
      double m_interval_end;
      double m_bin_size;
      long m_num_bins;
  };

}

#endif
