# -*- python -*-
# $Id: SConscript,v 1.15 2009/07/24 15:28:51 peachey Exp $
# Authors: James Peachey <peachey@lheamail.gsfc.nasa.gov>
# Version: rspgen-03-02-00

Import('baseEnv')
Import('listFiles')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('rspgenLib', depsOnly = 1)
rspgenLib = libEnv.StaticLibrary('rspgen', listFiles(['src/*.cxx']))

progEnv.Tool('rspgenLib')
gtrspgenBin = progEnv.Program('gtrspgen', listFiles(['src/gtrspgen/*.cxx']))
test_rspgenBin = progEnv.Program('test_rspgen', listFiles(['src/test/*.cxx']))

progEnv.Tool('registerTargets', package = 'rspgen', staticLibraryCxts = [[rspgenLib, libEnv]],
             binaryCxts = [[gtrspgenBin, progEnv]], testAppCxts = [[test_rspgenBin, progEnv]],
             includes = listFiles(['rspgen/*.h']),
             pfiles = listFiles(['pfiles/*.par']), data = listFiles(['data/*'], recursive = True))
