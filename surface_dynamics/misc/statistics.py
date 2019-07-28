# *************************************************************************
# Copyright (C) 2015-2016 Charles Fougeron <charlesfougeron@gmail.com>
#                         Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *************************************************************************

def mean_and_std_dev(l):
    r"""
    Return the mean and standard deviation of the floatting point numbers in
    the list ``l``.

    The implementation is very naive and should not be used for large list
    (>1000) of numbers.

    .. NOTE::

        mean and std are implemented in Sage but are quite buggy!
    """
    from math import sqrt
    m = sum(l) / len(l)
    if len(l) == 1:
        d = 0
    else:
        d = sum((x-m)**2 for x in l) / (len(l)-1)
    return m,sqrt(d)



