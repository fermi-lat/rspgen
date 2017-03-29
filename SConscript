# -*- python -*-
# $Id: SConscript,v 1.20 2013/10/09 22:26:17 asercion Exp $
# Authors: James Peachey <peachey@lheamail.gsfc.nasa.gov>, Joe Asercion <joseph.a.asercion@nasa.gov>
# Version: rspgen-03-02-03

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
