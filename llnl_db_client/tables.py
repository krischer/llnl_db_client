#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Tables for the CSS database files.

XXX: This is neither complete nor necessarily correct. Any subsequent work
should first check the indices against the CSS definitions. But it works for
our one chosen purpose so far.

:copyright:
    Lion Krischer (lion.krischer@gmail.com), 2018
:license:
    MIT
"""
import obspy

TABLES = {
    "wfdisc": [
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
        ("_filename", (213, 245), str)],

    "assoc": [
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
        ("(epoch) time of last record modification", (135, 152), str)],

    "arrival": [
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
        ("(epoch) time of last record modification", (206, 223), str)],

    "sitechan": [
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
        ("last_modified", (123, 140), str)],

    "site": [
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
        ("date", (138, 148), str)],

    "event": [
        ("event_id", (0, 8), int),
        ("location", (9, 25), str),
        ("unknown_a", (26, 33), int),
        ("unknown_b", (34, 55), str),
        ("unknown_c", (56, 58), int),
        ("date", (59, 70), str)],

    "wftag": [
        ("type", (0, 4), str),
        ("id", (6, 18), int),
        ("waveform_id", (19, 26), int),
        ("date", (27, 36), str)],

    "origin": [
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
        ("agency", (193, 207), str)],

    "sensor": [
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
        ("date", (122, 132), str)],

    "instrument": [
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
        ("data", (222, 232), str)]
}
