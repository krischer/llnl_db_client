#!/usr/bin/env python
# -*- coding: utf-8 -*-
u"""
LLNL Database Client for ObsPy.

:copyright:
    Lion Krischer (lion.krischer@gmail.com), 2018
:license:
    MIT
"""
from setuptools import find_packages, setup

import inspect
import os


DOCSTRING = __doc__.strip().split("\n")


def get_package_data():
    """
    Returns a list of all files needed for the installation relative to the
    "llnl_db_client" subfolder.
    """
    filenames = []
    # The root dir.
    root_dir = os.path.join(os.path.dirname(os.path.abspath(
        inspect.getfile(inspect.currentframe()))), "llnl_db_client")
    # Recursively include all files in these folders:
    folders = [os.path.join(root_dir, "tests", "data")]
    for folder in folders:
        for directory, _, files in os.walk(folder):
            for filename in files:
                # Exclude hidden files.
                if filename.startswith("."):
                    continue
                filenames.append(os.path.relpath(
                    os.path.join(directory, filename),
                    root_dir))
    return filenames


setup_config = dict(
    name="llnl_db_client",
    version="0.1.0",
    description=DOCSTRING[0],
    long_description="\n".join(DOCSTRING[2:]),
    author=u"Lion Krischer",
    author_email="lion.krischer@gmail.com",
    url="https://github.com/krischer/llnl_db_client",
    packages=find_packages(),
    package_data={
        "llnl_db_client": get_package_data()},
    license='MIT',
    platforms="OS Independent",
    install_requires=["obspy >= 1.1.0"],
    classifiers=[
        # complete classifier list:
        # http://pypi.python.org/pypi?%3Aaction=list_classifiers
        "Development Status :: 4 - Alpha",
        "Intended Audience :: Science/Research",
        "Operating System :: Unix",
        "Operating System :: POSIX",
        "Operating System :: MacOS",
        "Operating System :: Windows",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: Implementation :: CPython",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Physics"
        ],
)

if __name__ == "__main__":
    setup(**setup_config)
