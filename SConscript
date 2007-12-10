import glob,os,platform

Import('baseEnv')
Import('listFiles')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

rspgenLib = libEnv.StaticLibrary('rspgen', listFiles(['src/*.cxx']))

progEnv.Tool('rspgenLib')
gtrspgenBin = progEnv.Program('gtrspgen', listFiles(['src/gtrspgen/*.cxx']))

progEnv.Tool('registerObjects', package = 'rspgen', libraries = [rspgenLib], binaries = [gtrspgenBin], includes = listFiles(['rspgen/*.h']), pfiles = listFiles(['pfiles/*.par']))