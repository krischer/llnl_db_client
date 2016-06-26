#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Implementation of the client.

:copyright:
    Lion Krischer (krischer@geophysik.uni-muenchen.de), 2016
:license:
    GNU Lesser General Public License, Version 3
    (http://www.gnu.org/copyleft/lgpl.html)
"""
from __future__ import absolute_import, division, print_function

import os
import re
import warnings

import numpy as np
import pandas as pd

import obspy
from obspy.clients.base import BaseClient

from . import util

DTYPE = {
    # Big-endian integers
    's4': '>i',
    's2': '>h',
    # Little-endian integers
    'i4': '<i',
    'i2': '<h',
    # ASCII integers
    'c0': ('S12', np.int),
    'c#': ('S12', np.int),
    # Big-endian floating point
    't4': '>f',
    't8': '>d',
    # Little-endian floating point
    'f4': '<f',
    'f8': '<d',
    # ASCII floating point
    'a0': ('S15', np.float32),
    'a#': ('S15', np.float32),
    'b0': ('S24', np.float64),
    'b#': ('S24', np.float64),
}


class LLNLDBClient(BaseClient):
    def __init__(self, wfdisc, debug=False):
        """
        :param wfdisc: Path to the wfdisc file. Must be an existing filename.
        :type wfdisc: str
        """
        BaseClient.__init__(self, debug=debug)

        wfdisc = os.path.normpath(os.path.abspath(wfdisc))

        assert os.path.exists(wfdisc), "'%s' does not exist" % wfdisc
        assert os.path.isfile(wfdisc), "'%s' is not a file" % wfdisc

        self._basedir = os.path.dirname(wfdisc)
        self._db_name = re.sub(r"\.wfdisc$", "",
                               os.path.basename(wfdisc).strip())

        self._dataframes = {}

        self._assemble_filenames()
        self._parse_wf_disc_file()
        self._parse_site_file()

    def __str__(self):
        ret_str = "LLNL Database '%s' (%s)\n" % (self._db_name, self._basedir)
        ret_str += "\t%i waveform files\n" % len(self._dataframes["wfdisc"])
        ret_str += "\t%i stations\n" % len(self._dataframes["site"])
        return ret_str

    def _assemble_filenames(self):
        """
        Create a simple dictionary containing all the various filenames and
        make sure they all exist.
        """
        filetypes = [
            "arrival", "assoc", "etype", "event", "evids", "explosion",
            "instrument", "origin", "remark", "search_link", "sensor", "site",
            "sitechan", "springer_explosion", "wfdisc", "wftag"]

        files = {}

        for ftype in filetypes:
            filename = os.path.join(
                self._basedir, self._db_name + os.extsep + ftype)

            assert os.path.exists(filename), \
                "File '%s' does not exist." % filename
            files[ftype] = filename

        self._files = files

    def _parse_wf_disc_file(self):
        """
        Partial copy of obspy.io.css.core._read_css
        """
        with open(self._files["wfdisc"], "rt") as fh:
            lines = fh.readlines()

        items = []

        for line in lines:
            item = {}

            item["station"] = line[0:6].strip()
            item["channel"] = line[7:15].strip()
            item["starttime"] = obspy.UTCDateTime(float(line[16:33]))

            item["unknown_a"] = int(line[35:44].strip())
            item["unknown_b"] = int(line[46:51].strip())
            item["unknown_c"] = int(line[53:62].strip())
            item["unknown_d"] = float(line[63:80].strip())

            item["npts"] = int(line[82:87])
            item["sampling_rate"] = float(line[88:99])
            item["calib"] = float(line[100:116])
            item["calper"] = float(line[117:133])

            dirname = line[148:212].strip()
            filename = line[213:245].strip()
            item["filename"] = os.path.join(self._basedir, dirname, filename)

            if not os.path.exists(item["filename"]):
                warnings.warn(
                    "File '%s' does not exists. "
                    "Will not be accessible via the client.", UserWarning)
                continue

            item["offset"] = int(line[246:256])
            item["dtype"] = DTYPE[line[143:145]]
            item["read_fmt"] = np.dtype(item["dtype"])
            items.append(item)

        self._dataframes["wfdisc"] = pd.DataFrame(items)

    def plot_stations(self):
        _t = self._dataframes["site"]

        net = obspy.core.inventory.Network(code="LL")
        inv = obspy.core.inventory.Inventory(networks=[net], source="")
        for _, sta in _t.iterrows():
            net.stations.append(
                obspy.core.inventory.Station(
                    code=sta["code"],
                    latitude=sta["latitude"],
                    longitude=sta["longitude"],
                    elevation=sta["elevation_in_km"] * 1000.0))
        inv.plot(projection="local")

    def _parse_site_file(self):
        """
        Parse the site file.
        """
        definition = [
            ("code", (0, 7), str),
            ("unknown_a", (8, 16), int),
            ("unknown_b", (17, 26), int),
            ("latitude", (27, 34), float),
            ("longitude", (35, 46), float),
            ("elevation_in_km", (47, 54), float),
            ("description", (55, 105), str),
            ("unknown_c", (106, 110), str),
            ("code_2", (111, 120), str),
            ("unknown_d", (121, 130), float),
            ("unknown_e", (131, 137), float),
            ("date", (138, 148), str)
        ]

        self._dataframes["site"] = \
            util.to_dataframe(self._files["site"], definition)
