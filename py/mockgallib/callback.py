import mockgallib._mockgallib as c

def standby(f):
    """Call standby for nonzero MPI nodes
    Args:
      f is a callable that takes a list of floats
    """
    c._callback_standby(f)

def sync(x):
    """Call at the beginning of function f to sync
    Arg:
        x is the argument for f, a list of floats.
    """
    c._callback_sync(x)

def release():
    """Call this function when f will no longer called"""
    c._callback_release()
