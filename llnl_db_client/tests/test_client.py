#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Some very basic test. Requires the database to be available in a location
specified by the `LLNL_DB_PATH` environment variable.

This is largely a regression test suite.

:copyright:
    Lion Krischer (lion.krischer@gmail.com), 2018
:license:
    MIT
"""
import os

import numpy as np
import obspy
import pytest

from llnl_db_client import LLNLDBClient


@pytest.fixture(scope="session")
def db():
    if "LLNL_DB_PATH" not in os.environ:
        raise Exception("The `LLNL_DB_PATH` environemnt variable must be set.")

    return LLNLDBClient(os.environ["LLNL_DB_PATH"])


def test_str(db):
    assert "LLNL Database" in str(db)
    assert "11176 waveform files" in str(db)
    assert "147 stations" in str(db)


def test_get_inventory(db):
    inv = db.get_inventory()
    assert inv.get_contents()["networks"] == ["LL"]
    assert len(inv.get_contents()["stations"]) == 147
    assert len(inv.get_contents()["channels"]) == 1280


def test_misc_functions(db):
    assert len(db.get_unique_channels()) == 708
    assert len(db.list_events()) == 153


def test_get_waveforms(db):
    st = db.get_waveforms_for_event(522227)
    assert len(st) == 12
    assert [tr.id for tr in st] == [
        'LL.ELK..R', 'LL.ELK..T', 'LL.ELK..V', 'LL.KNB..R', 'LL.KNB..T',
        'LL.KNB..V', 'LL.LAC..R', 'LL.LAC..T', 'LL.LAC..V', 'LL.MNV..R',
        'LL.MNV..T', 'LL.MNV..V']
    assert st[0].stats.sampling_rate == 40.0
    assert st[0].stats.npts == 21002
    assert st[0].stats.starttime == \
        obspy.UTCDateTime(1982, 8, 5, 13, 57, 55, 86979)
    np.testing.assert_allclose(
        st[0].data[:10],
        np.array([-103., -106., -109., -112., -115., -118., -122., -124.,
                  -126., -128.], dtype=np.float32))


def test_get_event(db):
    ev = db.get_obspy_event(1778048)
    assert len(ev.origins, 4)
    assert len(ev.picks, 9)
    assert len(ev.magnitudes, 4)
    assert ev.preferred_origin().latitude == 36.7293
    assert ev.preferred_origin().longitude == -116.297
    assert ev.preferred_origin().time == \
        obspy.UTCDateTime(1992, 7, 4, 5, 57, 30, 830000)


def test_remove_response_does_something(db):
    st = db.get_waveforms_for_event(522227)
    np.testing.assert_allclose(st[0].data.ptp(), 167744.0)
    db.remove_response(st[0])
    np.testing.assert_allclose(st[0].data.ptp(), 0.02299499884053656)
