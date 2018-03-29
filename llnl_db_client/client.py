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

import obspy
import obspy.core.event
from obspy.core.compatibility import from_buffer

from . import tables, util


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

        self._parse_tables()
        self._check_wfdisc_files()
        self._parse_events()
        self._parse_sensor_information()

    def __str__(self):
        ret_str = "LLNL Database '%s' (%s)\n" % (self._db_name, self._basedir)
        ret_str += "\t%i waveform files\n" % len(self._dataframes["wfdisc"])
        ret_str += "\t%i stations\n" % len(self._dataframes["site"])
        return ret_str

    def _parse_tables(self):
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

        # Now read all the defined files.
        for name, definition in tables.TABLES.items():
            self._dataframes[name] = \
                util.to_dataframe(self._files[name], definition)

    def _check_wfdisc_files(self):
        """
        Check if the files in the wfdisc table exist.
        """
        files = []
        for _, row in self._dataframes["wfdisc"].iterrows():
            _f = os.path.join(self._basedir, row["_dirname"], row["_filename"])
            files.append(_f)
        exists = [_i for _i, _j in enumerate(files) if not os.path.exists(_j)]

        self._dataframes["wfdisc"]["filename"] = files

        for _i in exists:  # pragma: no cover
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

    def plot_stations(self):  # pragma: no cover
        self.get_inventory().plot(projection="local")

    def plot_events(self):  # pragma: no cover
        self.get_catalog().plot(projection="local")

    def plot(self):  # pragma: no cover
        fig = self.get_inventory().plot(projection="local", show=False)
        self.get_catalog().plot(fig=fig)

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

        event_info = self._dataframes["event"]
        assert set(event_info["event_id"]) == event_ids

        tags = self._dataframes["wftag"]
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

        origins = self._dataframes["origin"]
        missed_events = set()
        for i, row in origins.iterrows():
            # Safety measure.
            if row.event_id not in event_ids:  # pragma: no cover
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
        if missed_events:  # pragma: no cover
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
        sensor = self._dataframes["sensor"]

        # Step 2: Parse the instrument information to get the filename and the
        # type of response.
        instrument = self._dataframes["instrument"]

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
                if _arr.empty:  # pragma: no cover
                    continue
                if len(_arr) != 1:  # pragma: no cover
                    raise NotImplementedError

                # If the pick already exists - use it. Should not really happen
                # but who knows.
                if _arr_id in picks:  # pragma: no cover
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
        else:  # pragma: no cover
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
                    if line.startswith("*"):  # pragma: no cover
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
                        else:  # pragma: no cover
                            raise NotImplementedError
                    else:  # pragma: no cover
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
            else:  # pragma: no cover
                raise NotImplementedError

            tr.simulate(paz_remove=paz, water_level=water_level,
                        zero_mean=False, taper=False,
                        pre_filt=pre_filt)

        # Case 2: RESP files.
        elif epoch.response_type == "evresp":
            # Overwrite network temporarily so evalresp works as expected.
            original_network = tr.stats.network
            network_code = os.path.basename(epoch.filename).split(".")[1]
            tr.stats.network = network_code
            try:
                tr.simulate(seedresp={"units": output,
                                      "filename": epoch.filename},
                            water_level=water_level,
                            zero_mean=False, taper=False,
                            pre_filt=pre_filt)
            finally:
                tr.stats.network = original_network

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
                    else:  # pragma: no cover
                        raise NotImplementedError
                    continue
                elif len(line) == 4:
                    v = float(line[0]) + float(line[1]) * 1j
                    if cur_status == "poles":
                        cur_set["poles"].append(v)
                    elif cur_status == "zeros":
                        cur_set["zeros"].append(v)
                    else:  # pragma: no cover
                        raise NotImplementedError
                else:  # pragma: no cover
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
            else:  # pragma: no cover
                raise NotImplementedError

            tr.simulate(paz_remove=paz, water_level=water_level,
                        zero_mean=False, taper=False,
                        pre_filt=pre_filt)
        else:  # pragma: no cover
            raise NotImplementedError(
                "Unknown response type '%s' for file '%s'." % (
                    epoch.response_type, epoch.filename))
