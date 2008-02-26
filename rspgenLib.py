#$Id$
def generate(env, **kw):
    if not kw.get('depsOnly',0):
        env.Tool('addLibrary', library = ['rspgen'])
    env.Tool('astroLib')
    env.Tool('dataSubselectorLib')
    env.Tool('evtbinLib')
    env.Tool('irfsLib')
    env.Tool('st_appLib')
    env.Tool('st_facilitiesLib')
    env.Tool('st_streamLib')
    env.Tool('tipLib')

def exists(env):
    return 1
