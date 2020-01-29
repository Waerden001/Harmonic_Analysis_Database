
from collections import OrderedDict
SYMMETRIC_SPACE_KEY_DICT= {
    1: "space",
    2: "operator",
    3: "symmetry",
    4: "isotropy",
    5: "invariant measure",
    6: "discrete subgroup",
    7: "eigenfunctions",
    8: "eigenvalues",
    9: "fourier transform",
    10: "spectrum decomposition",
    11: "convolution",
    12: "differentiation",
    13: "fundamental solution to heat equation",
    14: "eigenfunctions on quotient",
    15: "eigenvalues on quotient",
    16: "fourier transform on quotient",
    17: "spectrum decomposition on quotient",
    18: "poisson summation",
    19: "partition function",
    20: "circle problem",
    21: "fundamental domain",
    22: "test function",
    23: "average",
}






class Symmetric_Space():
    def __init__(self, record):
        """Create a dictioary of information of the space

        """
        self._record = record

    @classmethod
    def from_input(cls):
        return cls({SYMMETRIC_SPACE_KEY_DICT[key]: input(SYMMETRIC_SPACE_KEY_DICT[key]+": ") for key in SYMMETRIC_SPACE_KEY_DICT})

    def save_to_database():
        pass

    def find_empty():
        pass

