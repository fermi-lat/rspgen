# -*- python -*-
# $Id: SConscript,v 1.16 2009/10/13 17:47:33 peachey Exp $
# Authors: James Peachey <peachey@lheamail.gsfc.nasa.gov>
# Version: rspgen-03-02-00

Import('baseEnv')
Import('listFiles')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

rspgenLib = libEnv.StaticLibrary('rspgen', listFiles(['src/*.cxx']))

progEnv.Tool('rspgenLib')
gtrspgenBin = progEnv.Program('gtrspgen', listFiles(['src/gtrspgen/*.cxx']))
test_rspgenBin = progEnv.Program('test_rspgen', listFiles(['src/test/*.cxx']))

progEnv.Tool('registerTargets', package = 'rspgen', staticLibraryCxts = [[rspgenLib, libEnv]],
             binaryCxts = [[gtrspgenBin, progEnv]], testAppCxts = [[test_rspgenBin, progEnv]],
             includes = listFiles(['rspgen/*.h']),
             pfiles = listFiles(['pfiles/*.par']), data = listFiles(['data/*'], recursive = True))
