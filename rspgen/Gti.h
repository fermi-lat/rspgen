/** \file Gti.h
*/
#ifndef rspgen_Gti_h
#define rspgen_Gti_h

#include <string>
#include <utility>
#include <vector>

namespace rspgen {

  class Gti {
    public:
      typedef std::pair<double, double> Interval_t;
      typedef std::vector<Interval_t> IntervalCont_t;
      typedef IntervalCont_t::const_iterator ConstIterator;

      Gti(const std::string & file_name, const std::string & ext_name = std::string("GTI"));

      /** \brief Compute fraction of time interval bounded by tstart and tstop which falls
          within one or more of the Good time intervals.
          \param tstart Interval start time.
          \param tstop Interval start time.
          \param gti_pos Iterator pointing to the first good time interval which might contain
          part of the given interval (a hint on which GTI to use for increased efficiency).
      */
      double getFraction(double tstart, double tstop, ConstIterator & gti_pos);

      ConstIterator begin() const;
      ConstIterator end() const;

    private:
      IntervalCont_t m_intervals;
  };

}

#endif
