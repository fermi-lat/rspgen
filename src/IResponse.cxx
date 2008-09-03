/** \file IResponse.cxx
    \brief Implementation of generic response calculator.
    \author James Peachey, HEASARC
*/
#include <map>
#include <memory>
#include <stdexcept>

#include "evtbin/OrderedBinner.h"

#include "irfInterface/IrfsFactory.h"
#include "irfLoader/Loader.h"

#include "rspgen/IResponse.h"
#include "rspgen/SpaceCraftCalculator.h"

#include "facilities/Util.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

namespace rspgen {

  const double IResponse::s_keV_per_MeV = 1000.;
  const double IResponse::s_MeV_per_keV = .001;

  // Cutoff value for compress Response Matrix.
  const double IResponse::lower_threshold = 0.0;
  // Set this to false to disable writing of compressed files.
  bool compressed = true;

  IResponse::IResponse(const std::string & resp_type, const std::string & spec_file, const evtbin::Binner * true_en_binner):
    m_os("IResponse", "IResponse(const std::string &...)", 2), m_kwds(), m_true_en_binner(0), m_app_en_binner(0), m_irfs() {
    // Process input ebounds extension to get apparent energy bin definitions..
    std::auto_ptr<const tip::Table> in_ebounds(tip::IFileSvc::instance().readTable(spec_file, "EBOUNDS"));

    // Get ebounds header for keyword access.
    const tip::Header & header = in_ebounds->getHeader();

    // List of standard keywords to read from ebounds if they are present.
    const char * keywords [] = { "TELESCOP", "INSTRUME", "DATE-OBS", "DATE-END", 0 };

    // Read them.
    header.get(keywords, m_kwds);

    // Add date keyword.
    m_kwds.push_back(tip::Header::KeyValPair_t("DATE", header.formatTime(time(0))));

    // Get number of channels currently in use.
    int detchans = 0;
    header["DETCHANS"].get(detchans);

    // Read apparent energy intervals from ebounds extension.
    tip::Index_t index = 0;
    evtbin::OrderedBinner::IntervalCont_t app_intervals(detchans);
    for (tip::Table::ConstIterator itor = in_ebounds->begin(); itor != in_ebounds->end(); ++itor) {
      (*itor)["CHANNEL"].get(index);
      if (index > detchans) continue; // Skip any rows with channel numbers > the number of channels.
      if (0 >= index) throw std::logic_error("Response constructor encountered a non-positive channel number");
      --index; // Arrays start at 0, not 1. This may not always be true, hence the check above.
      app_intervals[index] = evtbin::Binner::Interval(s_MeV_per_keV*(*itor)["E_MIN"].get(), s_MeV_per_keV*(*itor)["E_MAX"].get());
    }



    // To prevent memory leaks, first allocate memory into temporary auto_ptrs, then copy (release)
    // the pointers into the *real* member pointers.
    std::auto_ptr<evtbin::Binner> true_en_auto_ptr(0);

    // Clone true energy binner, if one was supplied.
    if (0 != true_en_binner) {
      true_en_auto_ptr.reset(true_en_binner->clone());
    }

    // Create apparent energy binner.
    std::auto_ptr<evtbin::Binner> app_en_binner(new evtbin::OrderedBinner(app_intervals));

    // Get container of irf names matching the given response name.
    irf_name_cont_type irf_name;
    lookUpResponse(resp_type, irf_name);
    if (irf_name.empty()) {
      SpaceCraftCalculator::displayIrfNames(m_os.info());
      throw std::runtime_error("IResponse cannot find response function(s) matching " + resp_type);
    }

    // Populate the irfs container with matching irfs.
    m_irfs.resize(irf_name.size(), 0);
    for (irf_name_cont_type::size_type index = 0; index != irf_name.size(); ++index) {
      m_irfs[index] = irfInterface::IrfsFactory::instance()->create(irf_name[index]);
    }

    // Everything succeeded, so release the pointers from their auto_ptrs.
    m_true_en_binner = true_en_auto_ptr.release();
    m_app_en_binner = app_en_binner.release();
  }

  IResponse::~IResponse() throw() {
    for (irf_cont_type::reverse_iterator itor = m_irfs.rbegin(); itor != m_irfs.rend(); ++itor) delete *itor;
    delete m_app_en_binner;
    delete m_true_en_binner;
  }

  void IResponse::writeOutput(const std::string & creator, const std::string & file_name, const std::string & fits_template) {
    if (0 == m_true_en_binner)
      throw std::logic_error("Cannot write response matrix without true energy bins defined");

    // Create output file.
    tip::IFileSvc::instance().createFile(file_name, fits_template);

    // Update keywords in the output file, using tip's file service keyword update mechanism.
    m_kwds.push_back(tip::Header::KeyValPair_t("CREATOR", creator));

    // Update FILENAME keyword in output file
    m_kwds.push_back(tip::Header::KeyValPair_t("FILENAME", facilities::Util::basename(file_name)));

    // Update output header with these keywords.
    tip::IFileSvc::instance().updateKeywords(file_name, m_kwds);

    // Open response table for writing.
    std::auto_ptr<tip::Table> resp_table(tip::IFileSvc::instance().editTable(file_name, "MATRIX"));

    // Get dimensions of matrix.
    long true_num_elem = m_true_en_binner->getNumBins();
    long app_num_elem = m_app_en_binner->getNumBins();

    // Get header of response table.
    tip::Header & header = resp_table->getHeader();

    // Explicitly set DETCHANS.
    header["DETCHANS"].set(app_num_elem);

    // Set LO_THRES keyword based on value we are using
    header["LO_THRES"].set(lower_threshold);

    // Resize the table to hold number of records == the number of true energy bins.
    resp_table->setNumRecords(true_num_elem);

    // Create a vector to hold one row of the matrix (index on apparent energy).
//    std::vector<double> response(app_num_elem);
    std::vector<double> response;
    // And a scratch one for the values that will get written.
    std::vector<double> scratch_response;

    // Compute values for other mandatory columns. (Must be vectors for vector valued columns in tip).
    // These are set assuming non-compressed rsp files and then fixed for compress ones later.
    int n_grp = 1;
    std::vector<int> f_chan(1, 1);
    std::vector<int> n_chan(1, app_num_elem);
    
    // Some counter variables for generating xspec compressed files.
    long counter;
    bool grpcnt = false;

    // Loop over true energy bins, computing response for each and writing it to the output response table.
    long true_idx = 0;
    tip::Table::Iterator true_itor = resp_table->begin();
    for (; true_itor != resp_table->end(); ++true_idx, ++true_itor) {
      // Get details of this true energy bin.
      evtbin::Binner::Interval true_en = m_true_en_binner->getInterval(true_idx);

      // Compute the response for this true energy.
      compute(true_en.midpoint(), response);

      if (compressed) {
	// Zero counters and clear vectors before each new row.
	n_grp = 0;
	grpcnt = false;
	counter = 0;
	f_chan.clear();
	n_chan.clear();
	scratch_response.clear();
	// Populate vecotrs for compressed rsp file.
	for(std::vector<double>::iterator r_itor = response.begin(); r_itor != response.end(); ++r_itor) {
	  counter++;
	  if (*r_itor <= lower_threshold) {
	    if (grpcnt) {
	      grpcnt = false;
	    }
	  } else {
	    scratch_response.push_back(*r_itor);
	    if (!grpcnt) {
	      grpcnt = true;
	      f_chan.push_back(counter);
	      n_grp++;
	      n_chan.resize(n_grp,0);
	    }
	    // Line above this ensures we are never writing past the end of n_chan.
	    n_chan[n_grp-1]++;
	  }
	}
	// Make sure that if f_chan and n_chan are smaller than three elements 
	// they get resized and padded with 0.
	if (f_chan.size() < 3) f_chan.resize(3,0);
	if (n_chan.size() < 3) n_chan.resize(3,0);
	response.swap(scratch_response);
      }

      (*true_itor)["ENERG_LO"].set(s_keV_per_MeV * true_en.begin());
      (*true_itor)["ENERG_HI"].set(s_keV_per_MeV * true_en.end());
      (*true_itor)["N_GRP"].set(n_grp);
      (*true_itor)["F_CHAN"].set(f_chan);
      (*true_itor)["N_CHAN"].set(n_chan);
      (*true_itor)["MATRIX"].set(response);
    }

    // Open ebounds table for writing.
    std::auto_ptr<tip::Table> ebounds(tip::IFileSvc::instance().editTable(file_name, "EBOUNDS"));

    // Set detchans explicitly.
    ebounds->getHeader()["DETCHANS"].set(app_num_elem);

    // Resize the table to hold number of records == the number of apparent energy bins.
    ebounds->setNumRecords(app_num_elem);

    // Loop over apparent energy bins, writing each to the ebounds extension.
    long app_idx = 0;
    tip::Table::Iterator app_itor = ebounds->begin();
    for (; app_itor != ebounds->end(); ++app_idx, ++app_itor) {
      // Get details of this apparent energy bin.
      evtbin::Binner::Interval app_en = m_app_en_binner->getInterval(app_idx);

      // Write bin parameters to output extension.
      (*app_itor)["CHANNEL"].set(app_idx + 1);
      (*app_itor)["E_MIN"].set(s_keV_per_MeV * app_en.begin());
      (*app_itor)["E_MAX"].set(s_keV_per_MeV * app_en.end());
    }

  }

  void IResponse::lookUpResponse(const std::string & resp, irf_name_cont_type & match) {
    SpaceCraftCalculator::lookUpResponse(resp, match);
  }

  IResponse::IResponse(): m_os("IResponse", "IResponse()", 2), m_true_en_binner(0), m_app_en_binner(0), m_irfs() {}

}
