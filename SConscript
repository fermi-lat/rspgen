# -*- python -*-
# $Id: SConscript,v 1.11 2009/06/24 18:30:44 glastrm Exp $
# Authors: James Peachey <peachey@lheamail.gsfc.nasa.gov>
# Version: rspgen-03-00-02

Import('baseEnv')
Import('listFiles')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('rspgenLib', depsOnly = 1)
rspgenLib = libEnv.StaticLibrary('rspgen', listFiles(['src/*.cxx']))

progEnv.Tool('rspgenLib')
gtrspgenBin = progEnv.Program('gtrspgen', listFiles(['src/gtrspgen/*.cxx']))
test_rspgenBin = progEnv.Program('test_rspgen', listFiles(['src/test/*.cxx']))

progEnv.Tool('registerObjects', package = 'rspgen', libraries = [rspgenLib], binaries = [gtrspgenBin], testApps = [test_rspgenBin], includes = listFiles(['rspgen/*.h']),
             pfiles = listFiles(['pfiles/*.par']), data = listFiles(['data/*'], recursive = True))
