#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
:copyright:
    Lion Krischer (lion.krischer@gmail.com), 2018
:license:
    MIT
"""
from __future__ import absolute_import, division, print_function

import collections

import pandas as pd


def to_dataframe(filename, definition):
    """
    Helper function converting the tabular files to pandas dataframes.
    """
    with open(filename, "rt") as fh:
        lines = fh.readlines()

    items = []
    for line in lines:
        item = collections.OrderedDict()

        for name, indices, type in definition:
            item[name] = type(line[indices[0]:indices[1]])
            if type == str:
                item[name] = item[name].strip()

        items.append(item)

    return pd.DataFrame(items)
