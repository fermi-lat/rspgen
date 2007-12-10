def generate(env, **kw):
    env.Tool('addLibrary', library = ['rspgen'], package = 'rspgen')
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
