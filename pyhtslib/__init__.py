#!/usr/bin/env python3
"""Module pyhtslib with common types and routines"""


class GenomeInterval:
    """Zero-based genome interval."""

    def __init__(self, seq, begin_pos, end_pos):
        self.seq = seq
        self.begin_pos = begin_pos
        self.end_pos = end_pos

    def __str__(self):
        return '{}:{:,}-{:,}'.format(self.seq, self.begin_pos + 1,
                                     self.end_pos)
