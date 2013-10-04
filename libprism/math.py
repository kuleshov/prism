# Copyright (c) 2012-2013 Illumina, Inc.
#
# This software is provided under the terms and conditions of the
# Illumina Open Source Software License 1.
#
# You should have received a copy of the Illumina Open Source
# Software License 1 along with this program. If not, see
# <https://github.com/downloads/sequencing/licenses/>.
#


from numpy import asarray, sum, log, exp, rollaxis

def logsumexp(a, axis=None):
    """Computes the log of the sum of exponentials of input elements.

    This function was originally taken from a version of SciPy.
    """
    a = asarray(a)
    if axis is None:
        a = a.ravel()
    else:
        a = rollaxis(a, axis)
    a_max = a.max(axis=0)
    s = sum(exp(a - a_max), axis=0)
    # keyboard()
    if a.ndim == 1 and s == 0:
        return float("-inf")
    else:
        out = log(s)
    out += a_max
    return out

def H(n):
    """Computes the n-th harmonic number."""

    # Euler-Mascheroni constant
    gamma = 0.577215665

    return gamma + log(n) + 0.5/n - 1./(12*n**2) + 1./(120*n**4)