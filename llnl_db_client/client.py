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

import collections
import io
import os
import re
import warnings

import numpy as np
import pandas as pd

import obspy
from obspy.core.compatibility import from_buffer

from . import util


class LLNLDBClient(object):
    def __init__(self, wfdisc):
        """
        :param wfdisc: Path to the wfdisc file. Must be an existing filename.
        :type wfdisc: str
        """
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
        self._parse_events()
        self._parse_sensor_information()

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

            item["id"] = int(line[35:44].strip())
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

    def _parse_events(self):
        # Step 1: Get list of events.
        with open(self._files["evids"], "rt") as fh:
            event_ids = set(int(_i.strip()) for _i in fh.readlines())

        # Step 2: Get event types.
        with open(self._files["etype"], "rt") as fh:
            event_type = {
                int(_j[0]): _j[1] for _j in
                (_i.strip().split() for _i in fh.readlines())}

        assert set(event_type.keys()) == event_ids

        # Step 3: Parse the event file.
        definition = [
            ("event_id", (0, 8), int),
            ("location", (9, 25), str),
            ("unknown_a", (26, 33), int),
            ("unknown_b", (34, 55), str),
            ("unknown_c", (56, 58), int),
            ("date", (59, 70), str)
        ]

        event_info = util.to_dataframe(self._files["event"], definition)
        assert set(event_info["event_id"]) == event_ids

        # Step 4: Parse the tags but only keep the event associations for now.
        definition = [
            ("type", (0, 4), str),
            ("id", (6, 18), int),
            ("waveform_id", (19, 26), int),
            ("date", (27, 36), str)
        ]
        tags = util.to_dataframe(self._files["wftag"], definition)
        tags = tags[tags.type == "evid"]

        # Make sure all events have associations.
        assert set(tags.id) == event_ids

        # Assemble all the information in some fashion for now.
        self._events = {}
        for ev_id in event_ids:
            ev = event_info[event_info.event_id == ev_id].iloc[0]

            self._events[ev_id] = {
                "location": ev.location,
                "unknown_a": ev.unknown_a,
                "unknown_b": ev.unknown_b,
                "unknown_c": ev.unknown_c,
                "date": ev.date,
                "waveform_ids": list(tags[tags.id == ev_id].waveform_id)
            }

    def get_waveforms_for_event(self, event_id):
        wf_ids = self._events[event_id]["waveform_ids"]
        _t = self._dataframes["wfdisc"]

        st = obspy.Stream()

        for wf in wf_ids:
            wf = _t[_t.id == wf].iloc[0]

            with io.open(wf.filename, "rb") as fh:
                data = fh.read(4 * wf.npts)

            data = from_buffer(data, dtype=np.float32)
            # Data is big-endian - we just want to work with little endian.
            data.byteswap(True)

            tr = obspy.Trace(data=data)
            tr.stats.station = wf.station
            tr.stats.sampling_rate = wf.sampling_rate
            tr.stats.starttime = wf.starttime
            tr.stats.channel = wf.channel.upper()
            tr.stats.calib = wf.calib

            st.append(tr)

        return st

    def list_events(self):
        return list(self._events.keys())

    def _parse_sensor_information(self):
        # Step 1: Parse the sensor file - we just need this to get the
        # instrument id for each channel.
        definition = [
            ("station", (0, 6), str),
            ("channel", (7, 14), str),
            ("starttime", (15, 33), lambda x: obspy.UTCDateTime(float(x))),
            ("endtime", (34, 52), lambda x: obspy.UTCDateTime(float(x))),
            ("instrument_id", (53, 61), int),
            ("unknown_c", (62, 70), int),
            ("unknown_d", (71, 82), int),
            ("unknown_e", (83, 101), float),
            ("unknown_f", (102, 113), float),
            ("unknown_g", (114, 119), float),
            ("unknown_h", (120, 121), str),
            ("date", (122, 132), str)
        ]
        sensor = util.to_dataframe(self._files["sensor"], definition)

        # Step 2: Parse the instrument information to get the filename and the
        # type of response.
        definition = [
            ("id", (0, 8), int),
            ("description", (9, 58), str),
            ("unknown_a", (59, 61), str),
            ("unknown_b", (62, 66), str),
            ("unknown_c", (67, 68), str),
            ("unknown_d", (69, 70), str),
            ("unknown_e", (71, 86), float),
            ("unknown_f", (87, 105), float),
            ("unknown_g", (106, 116), float),
            ("unknown_h", (117, 180), str),
            ("filename", (181, 213), str),
            ("response_type", (214, 221), str),
            ("data", (222, 232), str)
        ]
        instrument = util.to_dataframe(self._files["instrument"], definition)

        # Step 3: Find the actual filenames of all responses and match.
        responses_dir = os.path.join(self._basedir, "responses")
        all_files = {}
        for root, _, files in os.walk(responses_dir):
            for filename in files:
                assert filename not in all_files
                all_files[filename] = \
                    os.path.abspath(os.path.join(root, filename))

        Sensor = collections.namedtuple(
            "Sensor", ["starttime", "endtime", "id", "response_type",
                       "filename"])

        # Step 4: Merge it all.
        sensors = collections.defaultdict(list)
        for _, s in sensor.iterrows():
            inst = instrument[instrument.id == s.instrument_id].iloc[0]
            sensors[(s.station, s.channel.upper())].append(Sensor(
                s.starttime, s.endtime, s.instrument_id, inst.response_type,
                all_files[inst.filename]))

        self.sensors = sensors

    def remove_response(self, tr, output="VEL", water_level=60, pre_filt=None):
        """
        Remove the response of the passed waveform trace or stream.

        The passed arguments are documented in more detail in the main ObsPy
        package. No tapering or detrending is performed - the user is
        responsible for it.

        :param tr: Stream or Trace object whose instrument response(s) will be
            removed in-place.
        :type tr: :class:`obspy.core.trace.Trace` or
            :class:`obspy.core.stream.Stream`
        :param output: The output units. One of ``"DISP"``, ``"VEL"``,
            ``"ACC"``
        :param water_level: The water level.
        :param pre_filt: Frequency domain pre-filter.
        """
        assert output in ["DISP", "VEL", "ACC"]

        if isinstance(tr, obspy.Stream):
            for i in tr:
                self.remove_response(
                    i, output=output, water_level=water_level,
                    pre_filt=pre_filt)
            return

        # Step 1: Find the corresponding response.
        assert (tr.stats.station, tr.stats.channel) in self.sensors, \
            "No response for station '%s' and channel '%s' found." % (
                tr.stats.station, tr.stats.channel)

        sensor = self.sensors[(tr.stats.station, tr.stats.channel)]

        time = tr.stats.starttime + \
            (tr.stats.endtime - tr.stats.starttime) / 2.0

        for epoch in sensor:
            if epoch.starttime <= time <= epoch.endtime:
                break
        else:
            raise ValueError("Found some responses for the channel but not "
                             "for the correct time span.")

        # Step 2: Actually correct the data.
        #
        # Case 1: sacpz files.
        if epoch.response_type == "sacpzf":
            with io.open(epoch.filename, "rt") as fh:
                cur_state = None
                poles = []
                zeros = []
                for line in fh:
                    line = line.strip()
                    if line.startswith("*"):
                        continue

                    l = line.split()
                    if l[0].upper() == "ZEROS":
                        cur_state = "ZEROS"
                        num_zeros = int(l[1])
                    elif l[0].upper() == "POLES":
                        cur_state = "POLES"
                        num_poles = int(l[1])
                    elif l[0].upper() == "CONSTANT":
                        constant = float(l[1])
                    elif cur_state:
                        v = float(l[0]) + float(l[1]) * 1j
                        if cur_state == "ZEROS":
                            zeros.append(v)
                        elif cur_state == "POLES":
                            poles.append(v)
                        else:
                            raise NotImplementedError
                    else:
                        raise NotImplementedError

            assert len(poles) == num_poles, "sacpz parsing error"
            assert len(zeros) == num_zeros, "sacpz parsing error"
            assert constant, "sacpz parsing error"

            paz = {"poles": poles, "zeros": zeros, "gain": 1.0,
                   "sensitivity": constant}

            # Assume they correct to velocity.
            if output == "DISP":
                paz["zeros"].append(0 + 0j)
            elif output == "ACC":
                paz["zeros"] = paz["zeros"][:-1]

            tr.simulate(paz_remove=paz, water_level=water_level,
                        zero_mean=False, taper=False,
                        pre_filt=pre_filt)

        # Case 2: RESP files.
        elif epoch.response_type == "evresp":
            tr.simulate(seedresp={"units": output,
                                  "filename": epoch.filename},
                        water_level=water_level,
                        zero_mean=False, taper=False,
                        pre_filt=pre_filt)

        # Case 3: Funky response/paz files.
        elif epoch.response_type == "paz":
            with io.open(epoch.filename, "rt") as fh:
                lines = [_i.strip() for _i in fh.readlines()]
            lines = [_i for _i in lines if not _i.startswith("#")]

            paz_sets = []
            cur_set = None
            cur_status = None

            for line in lines:
                if line.startswith("theoretical"):
                    if cur_set:
                        paz_sets.append(cur_set)
                    cur_set = {"poles": [], "zeros": []}
                    paz_sets.append(cur_set)
                    cur_status = None
                    continue

                line = line.split()

                if len(line) == 2:
                    cur_set["sensitivity"] = float(line[0])
                    continue
                elif len(line) == 1:
                    if cur_status is None:
                        cur_status = "poles"
                    elif cur_status == "poles":
                        cur_status = "zeros"
                    else:
                        raise NotImplementedError
                    continue
                elif len(line) == 4:
                    v = float(line[0]) + float(line[1]) * 1j
                    if cur_status == "poles":
                        cur_set["poles"].append(v)
                    elif cur_status == "zeros":
                        cur_set["zeros"].append(v)
                    else:
                        raise NotImplementedError
                else:
                    raise NotImplementedError

            paz = paz_sets[0]
            for p in paz_sets[1:]:
                paz["sensitivity"] *= p["sensitivity"]
            paz["gain"] = 1.0

            # Assume they correct to velocity.
            if output == "DISP":
                paz["zeros"].append(0 + 0j)
            elif output == "ACC":
                paz["zeros"] = paz["zeros"][:-1]

            tr.simulate(paz_remove=paz, water_level=water_level,
                        zero_mean=False, taper=False,
                        pre_filt=pre_filt)
        else:
            raise NotImplementedError(
                "Unknown response type '%s' for file '%s'." % (
                    epoch.response_type, epoch.filename))
