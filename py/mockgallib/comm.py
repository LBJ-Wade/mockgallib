import mockgallib._mockgallib as c
import atexit

rank= -1
n_nodes= 0

def init():
    c._comm_init()
    atexit.register(c._comm_finalise)

    global rank
    rank = c._comm_this_rank()

    global n_nodes
    n_nodes = c._comm_n_nodes()

def bcast_str(src_str):
    """Broadcast string to other MPI nodes from rank=0
    Args:
        src_str (str): sourse sring at rank=0
    Return value:
        str:           src_str in all MPI nodes
    """
    return c._comm_bcast_str(src_str)

def bcast_int(src_int):
    """Broadcast an integer to other MPI nodes from rank=0
    Args:
        src_int (int): sourse sring at rank=0
    Return value:
        int:           src_int in all MPI nodes
    """
    return c._comm_bcast_int(src_int)
