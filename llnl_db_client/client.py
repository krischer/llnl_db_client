#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Implementation of the client.

:copyright:
    Lion Krischer (lion.krischer@gmail.com), 2018
:license:
    MIT
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
import obspy.core.event
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
        self._parse_sitechan_file()
        self._parse_events()
        self._parse_sensor_information()
        self._parse_arrival_file()
        self._parse_assoc_file()

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
        definition = [
            ("station", (0, 6), str),
            ("channel", (7, 15), str),
            ("starttime", (16, 33), lambda x: obspy.UTCDateTime(float(x))),
            ("id", (35, 44), int),
            ("unknown_b", (46, 51), int),
            ("unknown_c", (53, 62), int),
            ("unknown_d", (63, 80), float),
            ("npts", (82, 87), int),
            ("sampling_rate", (88, 99), float),
            ("calib", (100, 116), float),
            ("calper", (117, 133), float),
            ("_dirname", (148, 212), str),
            ("_filename", (213, 245), str)]

        self._dataframes["wfdisc"] = \
            util.to_dataframe(self._files["wfdisc"], definition)

        files = []
        for _, row in self._dataframes["wfdisc"].iterrows():
            _f = os.path.join(self._basedir, row["_dirname"], row["_filename"])
            files.append(_f)
        exists = [_i for _i, _j in enumerate(files) if not os.path.exists(_j)]

        self._dataframes["wfdisc"]["filename"] = files

        for _i in exists:
            msg = ("File '%s' does not exists. "
                   "Will not be accessible via the client.") % \
                   self._dataframes["wfdisc"]["filename"][_i]
            warnings.warn(msg)

        # Get rid of all non-existant files.
        self._dataframes["wfdisc"].drop(exists, inplace=True)

    def get_inventory(self):
        _t = self._dataframes["site"]
        _sc = self._dataframes["sitechan"]

        net = obspy.core.inventory.Network(code="LL")
        inv = obspy.core.inventory.Inventory(networks=[net], source="")
        for _, sta in _t.iterrows():
            latitude = sta["latitude"]
            longitude = sta["longitude"]
            elevation = sta["elevation_in_km"] * 1000.0

            # Get all channels for that station.
            channels = []
            chas = _sc[_sc["station"] == sta["code"]]
            for c in chas.iterrows():
                c = c[1]

                azimuth = c.azimuth
                # The azimuth of the files for vertical channels is always
                # set to -1. StationXML uses 0 for that thus we patch it
                # here. But the azimuth for vertical channels is not really
                # defined in any case.
                if azimuth == -1.0 and c.dip == 0:
                    azimuth = 0
                # Make sure its always in range.
                azimuth %= 360

                # The dip is apparently defined as angle against the
                # vertical axis. SEED defines it as angle against the ground.
                dip = c.dip
                dip -= 90.0

                channels.append(obspy.core.inventory.channel.Channel(
                    code=c.channel,
                    location_code="",  # XXX: Cannot find the location.
                    latitude=latitude,
                    longitude=longitude,
                    elevation=elevation,
                    # Its in kilometer!
                    depth=c.emplacement_depth * 1000,
                    azimuth=azimuth,
                    dip=dip))

            net.stations.append(
                obspy.core.inventory.Station(
                    code=sta["code"],
                    latitude=latitude,
                    longitude=longitude,
                    elevation=elevation,
                    channels=channels))
        return inv

    def plot_stations(self):
        self.get_inventory().plot(projection="local")

    def _parse_assoc_file(self):
        definition = [
            ("arrival id", (0, 8), int),
            ("origin id", (9, 17), int),
            ("station", (18, 24), str),
            ("associated phase", (25, 33), str),
            ("phase confidence", (34, 38), float),
            ("station to event distance", (39, 47), float),
            ("station to event azimuth", (48, 55), float),
            ("event to station azimuth", (56, 63), float),
            ("time residual", (64, 72), float),
            ("time = defining, non-defining", (73, 74), str),
            ("azimuth residual", (75, 82), float),
            ("azimuth = defining, non-defining", (83, 84), str),
            ("slowness residual", (85, 92), float),
            ("slowness = defining, non-defining", (93, 94), str),
            ("incidence angle residual", (95, 102), float),
            ("location weight", (103, 109), float),
            ("velocity model", (110, 125), str),
            ("comment id", (126, 134), int),
            ("(epoch) time of last record modification", (135, 152), str)]
        self._dataframes["assoc"] = \
            util.to_dataframe(self._files["assoc"], definition)

    def _parse_arrival_file(self):
        definition = [
            ("station", (0, 6), str),
            ("epoch time", (7, 24), float),
            ("arrival id", (25, 33), int),
            ("julian date", (34, 42), int),
            ("stassoc id", (43, 51), int),
            ("channel operation id", (52, 60), int),
            ("channel", (61, 69), str),
            ("reported phase", (70, 78), str),
            ("signal type", (79, 80), str),
            ("delta time", (81, 87), float),
            ("observed azimuth", (88, 95), float),
            ("delta azimuth", (96, 103), float),
            ("observed slowness (s/deg)", (104, 111), float),
            ("delta slowness", (112, 119), float),
            ("emergence angle", (120, 127), float),
            ("rectilinearity", (128, 135), float),
            ("amplitude, infloat corrected, nm", (136, 146), float),
            ("period", (147, 154), float),
            ("log(amp/per)", (155, 162), float),
            ("clipped flag", (163, 164), str),
            ("first motion", (165, 167), str),
            ("signal to noise ratio", (168, 178), float),
            ("signal onset quality", (179, 180), str),
            ("source/originator", (181, 196), str),
            ("comment id", (197, 205), int),
            ("(epoch) time of last record modification", (206, 223), str)]
        self._dataframes["arrival"] = \
            util.to_dataframe(self._files["arrival"], definition)

    def _parse_sitechan_file(self):
        """
        Parse the site chan file.
        """
        definition = [
            ("station", (0, 6), str),
            ("channel", (7, 15), str),
            ("julday_start", (16, 24), int),
            ("channel_operation_id", (25, 33), int),
            ("julday_end", (34, 42), int),
            ("channel_type", (43, 47), str),
            ("emplacement_depth", (48, 57), float),
            ("azimuth", (58, 64), float),
            ("dip", (65, 71), float),
            ("channel_description", (72, 122), str),
            ("last_modified", (123, 140), str)
        ]
        self._dataframes["sitechan"] = \
            util.to_dataframe(self._files["sitechan"], definition)

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
                "waveform_ids": list(tags[tags.id == ev_id].waveform_id),
                "origins": collections.OrderedDict()
            }

        # Also parse all the origins.
        definition = [
            ("latitude", (0, 9), float),
            ("longitude", (10, 18), float),
            ("depth_in_km", (19, 29), float),
            ("origin_time", (30, 47), float),
            ("origin_id", (48, 56), int),
            ("event_id", (58, 66), int),
            ("body_wave_magnitude", (128, 136), float),
            ("surface_wave_magnitude", (145, 152), float),
            ("local_magnitude", (162, 169), float),
            ("mlid", (170, 178), int),
            ("agency", (193, 207), str)
        ]

        origins = util.to_dataframe(self._files["origin"], definition)
        missed_events = set()
        for i, row in origins.iterrows():
            # The database is just really incomplete :-(
            if row.event_id not in event_ids:
                missed_events.add(row.event_id)
                continue

            self._events[row.event_id]["origins"][row.agency] = {
                "origin_id": row.origin_id,
                "body_wave_magnitude": row.body_wave_magnitude,
                "surface_wave_magnitude": row.surface_wave_magnitude,
                "local_magnitude": row.local_magnitude,
                "latitude": row.latitude,
                "longitude": row.longitude,
                "depth_in_m": row.depth_in_km * 1000.0,
                "origin_time": obspy.UTCDateTime(row.origin_time)
            }
        if missed_events:
            print("Skipped origins of %i events as the corresponding events "
                  "are not part of the database." % len(missed_events))

        # Finally make sure all events have origins.
        x_o = [
            ev for ev, value in self._events.items() if not value["origins"]]
        assert not x_o

    def get_unique_channels(self):
        return set((_i[1].station.upper(),
                    _i[1].channel.upper()) for _i in
                   self._dataframes["wfdisc"].iterrows())

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
            tr.stats.network = "LL"
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
            sensors[(s.station.upper(), s.channel.upper())].append(Sensor(
                s.starttime, s.endtime, s.instrument_id, inst.response_type,
                all_files[inst.filename]))

        self.sensors = sensors

    def get_obspy_event(self, event):
        """
        Convert a given event to an ObsPy event object.
        """
        event_id = event
        event = self._events[event]

        ev_obj = obspy.core.event.Event()

        ev_obj.event_descriptions.append(obspy.core.event.EventDescription(
            str(event_id), type="earthquake name"))

        assoc = self._dataframes["assoc"]
        arr = self._dataframes["arrival"]

        picks = {}

        for key, origin in event["origins"].items():
            org = obspy.core.event.Origin(
                longitude=origin["longitude"],
                latitude=origin["latitude"],
                depth=origin["depth_in_m"],
                time=origin["origin_time"],
                creation_info=obspy.core.event.CreationInfo(agency_id=key)
            )
            ev_obj.origins.append(org)

            # Map to QuakeML shortcuts.
            mag_dict = {
                "body_wave_magnitude": "Mb",
                "surface_wave_magnitude": "MS",
                "local_magnitude": "ML"}

            for mag in ["body_wave_magnitude", "surface_wave_magnitude",
                        "local_magnitude"]:
                # Not set.
                if origin[mag] < -100:
                    continue
                ev_obj.magnitudes.append(obspy.core.event.Magnitude(
                    mag=origin[mag],
                    magnitude_type=mag_dict[mag],
                    origin_id=str(org.resource_id.resource_id)
                ))

            # Multi-step process - we first find all arrivals associations
            # for the given origin.
            _a = assoc[assoc["origin id"] == origin["origin_id"]]
            if _a.empty:
                continue
            # Now find the arrival for the given association.
            for _, _assoc in _a.iterrows():
                _arr_id = int(_assoc["arrival id"])
                _arr = arr[arr["arrival id"] == _arr_id]
                if _arr.empty:
                    continue
                if len(_arr) != 1:
                    raise

                # If the pick already exists - use it.
                if _arr_id in picks:
                    p = picks[_arr_id]
                # Otherwise create it.
                else:
                    pick_time = obspy.UTCDateTime(float(_arr["epoch time"]))
                    wf_id = obspy.core.event.base.WaveformStreamID(
                        network_code="LL",
                        station_code=_arr["station"].tolist()[0],
                        location_code="",
                        channel_code=_arr["channel"].tolist()[0])
                    phase = _arr["reported phase"].tolist()[0]

                    p = obspy.core.event.Pick(
                        time=pick_time,
                        waveform_id=wf_id,
                        phase_hint=phase)
                    picks[_arr_id] = p
                    ev_obj.picks.append(p)

                # We now use this pick and with the assoc create a new
                # arrival - the terminoloy of the database differs a bit
                # with QuakeML - don't be confused ^^
                a = obspy.core.event.Arrival(
                    pick_id=p.resource_id,
                    phase=p.phase_hint)
                org.arrivals.append(a)

        # Find one of the preferred origin ids, otherwise just give a random
        # one.
        # for agency in ['ISC', 'PDE-M']:   # original code
        # full list of agencies from the LLNL database
        for agency in ['DOE_Springer', 'CNSS', 'DEWEY_USMINES', 'LLNL_GT',
                       'PDE-M', 'UNR_2000_2001', 'UNR_2002present',
                       'UNR_SGB_reloc',  'Weston_BlMesa']:
            for origin in ev_obj.origins:
                if origin.creation_info.agency_id.lower() == agency.lower():
                    ev_obj.preferred_origin_id = str(
                        origin.resource_id.resource_id)
                    break
            else:
                continue
            break
        else:
            ev_obj.preferred_origin_id = str(ev_obj.origins[
                                                 0].resource_id.resource_id)

        return ev_obj

    def get_catalog(self):
        cat = obspy.core.event.Catalog()
        for event in self.list_events():
            cat.append(self.get_obspy_event(event))
        return cat

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

                    _l = line.split()
                    if _l[0].upper() == "ZEROS":
                        cur_state = "ZEROS"
                        num_zeros = int(_l[1])
                    elif _l[0].upper() == "POLES":
                        cur_state = "POLES"
                        num_poles = int(_l[1])
                    elif _l[0].upper() == "CONSTANT":
                        constant = float(_l[1])
                    elif cur_state:
                        v = float(_l[0]) + float(_l[1]) * 1j
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

            # XXX: This is a wild guess - it guesses they correct to
            # displacment in micrometer - the only justification for this is
            # that I looked at the seismograms corrected with RESP files and
            # this puts the amplitudes in the same ballpark.
            constant *= 1E6

            paz = {"poles": poles, "zeros": zeros, "gain": constant,
                   "sensitivity": 1.0}

            # Assume they correct to displacement. Adding a zero effectively
            # differentiates.
            if output == "DISP":
                pass
            elif output == "VEL":
                paz["zeros"].append(0 + 0j)
            elif output == "ACC":
                paz["zeros"].append(0 + 0j)
                paz["zeros"].append(0 + 0j)
            else:
                raise NotImplementedError

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
                    cur_set["gain"] = float(line[0])
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
                paz["gain"] *= p["gain"]
            paz["sensitivity"] = 1.0

            # XXX: This is a wild guess - it guesses they correct to
            # displacment in 1E5 meters for whatever reason - the only
            # justification for this is that I looked at the seismograms
            # corrected with RESP files and this puts the amplitudes in the
            # same ballpark.
            paz["gain"] *= 1E-5

            # Assume they correct to displacement. Adding a zero effectively
            # differentiates.
            if output == "DISP":
                pass
            elif output == "VEL":
                paz["zeros"].append(0 + 0j)
            elif output == "ACC":
                paz["zeros"].append(0 + 0j)
                paz["zeros"].append(0 + 0j)
            else:
                raise NotImplementedError

            tr.simulate(paz_remove=paz, water_level=water_level,
                        zero_mean=False, taper=False,
                        pre_filt=pre_filt)
        else:
            raise NotImplementedError(
                "Unknown response type '%s' for file '%s'." % (
                    epoch.response_type, epoch.filename))
