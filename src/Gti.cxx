#include <iostream>
#include <memory>

#include "rspgen/Gti.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

namespace rspgen {
  Gti::Gti(const std::string & file_name, const std::string & ext_name): m_intervals() {
    // Open GTI extension.
    std::auto_ptr<const tip::Table> gti_table(tip::IFileSvc::instance().readTable(file_name, ext_name));

    // Resize the interval container.
    m_intervals.resize(gti_table->getNumRecords());

    // Start with first interval.
    IntervalCont_t::iterator int_itor = m_intervals.begin();

    // Fill container with intervals from the extension.
    for (tip::Table::ConstIterator itor = gti_table->begin(); itor != gti_table->end(); ++itor, ++int_itor) {
      double start = (*itor)["START"].get();
      double stop = (*itor)["STOP"].get();
      *int_itor = Interval_t(start, stop);
    }
  }

  double Gti::getFraction(double tstart, double tstop, ConstIterator & gti_pos) {
    double fraction = 0.;

    for (; gti_pos != m_intervals.end(); ++gti_pos) { 
      // Check if this interval ends before GTI starts and return 0. fraction if it does.
      if (tstop <= gti_pos->first) break;

      // Check if this interval is completely contained in the GTI and return 1 if it does.
      if (tstart >= gti_pos->first && tstop <= gti_pos->second) {
        if (tstop == gti_pos->second) ++gti_pos;
        fraction = 1.;
        break;
      }

      // Check if there is some overlap and add that overlap.
      if (tstart < gti_pos->second) {
        double start = tstart > gti_pos->first ? tstart : gti_pos->first;
        double stop = tstop < gti_pos->second ? tstop : gti_pos->second;
        fraction += (stop - start) / (tstop - tstart);
      }

      // Check if this GTI still has some part which might overlap some future interval.
      // If it does, break to avoid incrementing the GTI iterator.
      if (tstop < gti_pos->second) break;
    }

    return fraction;
  }

  Gti::ConstIterator Gti::begin() const { return m_intervals.begin(); }

  Gti::ConstIterator Gti::end() const { return m_intervals.end(); }

}
