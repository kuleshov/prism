#!/usr/bin/env python
#
# Copyright (c) 2012-2013 Illumina, Inc.
#
# This software is provided under the terms and conditions of the
# Illumina Open Source Software License 1.
#
# You should have received a copy of the Illumina Open Source
# Software License 1 along with this program. If not, see
# <https://github.com/downloads/sequencing/licenses/>.
#


class Positions:
    """Represents the set of positions over which the HMM will be run.

    A Positions object is loaded at the start of a Prism run; it stores
    the set of genomic positions over which the HMM will be run and provides
    convenient helper functions for tasks such as retrieving the next 
    position, or retrieving the index of a position.
    
    """
    def __init__(self, pos_file_path, start=0, end=2**40):
        """Initializes and loads from a .pos file. """

        positions = list()
        index_by_position = dict()

        with open(pos_file_path) as pos_file:
            idx = 0
            for line in pos_file:
                i = line.find(',')
                pos = int(line[0:i])

                if not (start <= pos <= end): continue

                positions.append(pos)
                index_by_position[pos] = idx
                idx += 1

        self.positions = sorted(positions)
        self.positions_set = set(positions)
        self.index_by_position = index_by_position

    def __contains__(self, pos):
        """Indicates whether a genomic position is among those considered."""
        return pos in self.positions_set

    def __iter__(self):
        """Iterates over all genomic positions."""
        return self.positions.__iter__()

    def __reversed__(self):
        """Iterates over all genomic positions in reverse order."""
        return reversed(self.positions)

    def __len__(self):
        """Indicates number of positions."""
        return len(self.positions)

    def get_position_index(self, pos):
        """Returns the index of a position in the sorted list of positions."""
        return self.index_by_position[pos]

    def get_position_by_index(self, i):
        """Returns the i-th genomic position."""
        return self.positions[i]

    def get_previous_position(self, pos):
        """ Returns preceding genomic position. """
        if self.is_start_position(pos):
            exit("Error: Cannot get previous position")

        j = self.index_by_position[pos]
        return self.positions[j-1]

    def get_next_position(self, pos):
        """Returns next genomic position"""
        if self.is_end_position(pos):
            exit("Error: Cannot get next position")

        j = self.index_by_position[pos]
        return self.positions[j+1]

    def is_start_position(self, pos):
        """Indicates whether this is the first genomic position."""
        return (pos == self.positions[0])

    def is_end_position(self, pos):
        """Indicates whether this is the last genomic position."""
        return (pos == self.positions[-1])

    # TODO: Replace is_... by these:
    @property
    def first_position(self):
        """Returns first genomic position."""
        return self.positions[0]

    @property
    def last_position(self):
        """Returns last genomic position."""
        return self.positions[-1]        