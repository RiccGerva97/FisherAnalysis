import numpy as np
from dataclasses import dataclass

def to_uint64(x: int) -> int:
    """Clip an integer so that it occupies 64 bits"""
    return x & 0xffffffffffffffff


def to_uint32(x: int) -> int:
    """Clip an integer so that it occupies 32 bits"""
    return x & 0xffffffff

@dataclass
class PCG:
    """PCG Uniform Pseudo-random Number Generator"""
    state: int = 0
    inc: int = 0
    def __init__(self, init_state=42, init_seq=54):
        # 64-bit
        self.state = 0
        # 64-bit
        self.inc = (init_seq << 1) | 1
        self.random()
        # 64-bit
        self.state += init_state
        self.random()
    def random(self):
        """Return a new random number and advance PCG's internal state"""
        # 64-bit
        oldstate = self.state
        # 64-bit
        self.state = to_uint64((oldstate * 6364136223846793005 + self.inc))
        # 32-bit
        xorshifted = to_uint32((((oldstate >> 18) ^ oldstate) >> 27))
        # 32-bit
        rot = oldstate >> 59
        # 32-bit
        return to_uint32((xorshifted >> rot) | (xorshifted << ((-rot) & 31)))
    def random_float(self):
        """Return a new random number uniformly distributed over [0, 1]"""
        return self.random() / 0xffffffff

def randind(L_i:int, dil:int):
    """
    Gives a list of non-repeated array indeces of `dil` elements.

    Arguments:
    - `L_i`: `int`, length of the inizial array to be diluited
    - `dil`: `int`, final lenght of the array

    Returns:
    - `randind`: uniform distributed random indices
    """
    r = PCG()
    indlist = []
    if L_i/2<=dil:
        randind = np.arange(L_i)
        for i in range(L_i-dil):
            a = int(L_i*r.random_float())
            while a in indlist:
                a = int(L_i*r.random_float())
            indlist.append(a)
        randind = np.delete(randind, indlist)
    else:
        for i in range(dil):
            a = int(L_i*r.random_float())
            while a in indlist:
                a = int(L_i*r.random_float())
            indlist.append(a)
        randind = indlist
    
    return randind