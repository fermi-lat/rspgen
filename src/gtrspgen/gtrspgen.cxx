/** \file gtrspgen.cxx
    \brief Main gtrspgen application.
    \author James Peachey, HEASARC/GSSC
*/
#include "st_app/StAppFactory.h"
#include "RspGenApp.h"
st_app::StAppFactory<rspgen::RspGenApp> g_factory("gtrspgen");
