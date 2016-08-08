import mockgallib._mockgallib as c

def minimise(f, function_callback, x0, step_size):
    """Minimise function f

    Args:
        f (function)
        function_callback: callback function (may be None)
        x0 (array): starting point of parameters
        step_size (array): initial step size of parameter change
    """
    return c._minimise(f, function_callback, x0, step_size)
