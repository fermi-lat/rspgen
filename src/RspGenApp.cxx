/** \file RspGenApp.cxx
    \brief Implementation of main rspgen application class.
    \author James Peachey, HEASARC
*/
#include <cstdlib>
#include <memory>

#include "RspGenApp.h"
#include "rspgen/GrbResponse.h"

namespace rspgen {

  void RspGenApp::run() {
    st_app::AppParGroup & pars(getParGroup("rspgen"));
    prompt(pars);
    writeResponse(pars);
  }

  void RspGenApp::prompt(st_app::AppParGroup & pars) {
    pars.Prompt("specfile");
    pars.Prompt("scfile");
    pars.Prompt("outfile");
    pars.Prompt("ra");
    pars.Prompt("dec");
    pars.Prompt("time");
    pars.Prompt("psfradius");
    pars.Prompt("resptype");
    pars.Prompt("resptpl");

    // Prompt for (true) energy binning parameters.
    m_bin_config.energyParPrompt(pars);

    pars.Save();
  }

  void RspGenApp::writeResponse(const st_app::AppParGroup & pars) const {
    // Create energy binner from related parameters.
    std::auto_ptr<evtbin::Binner> true_en_binner(m_bin_config.createEnergyBinner(pars));

    // Create response object.
    GrbResponse response(pars["ra"], pars["dec"], pars["time"], pars["psfradius"], pars["resptype"], pars["specfile"],
      pars["scfile"], true_en_binner.get());

    // Get name of template for output file.
    std::string resp_tpl = pars["resptpl"];

    // If it was not defined, get default template file name.
    if (resp_tpl.empty()) resp_tpl = getDataDir() + "LatResponseTemplate";

    // Write the output response file.
    response.writeOutput("rspgen", pars["outfile"], resp_tpl);

  }

  std::string RspGenApp::getDataDir() const {
    std::string retval;
    const char * value = getenv("RSPGENROOT");
    if (0 != value) retval = std::string(value) + "/data/";
    return retval;
  }

}
