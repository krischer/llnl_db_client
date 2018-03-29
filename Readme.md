# LLNL Database Client

### Installation

Assuming you use conda, run (in the correct `conda` environment):

```bash
$ conda install -c conda-forge obspy pandas basemap pytest
```

Assuming `ObsPy` and everything is installed, grab this repository with git, and install it with `pip`. Make sure to use the `-e` flag for an editable installation.

```bash
$ git clone https://github.com/krischer/llnl_db_client.git
$ cd llnl_db_client
$ pip install -v -e .
```

If you ever want to update the installation, just run

```bash
$ git pull
$ pip install -v -e .
```

within the repository.


### Running the Tests

This is currently only a regression test suite but it at least checks some things. You need the database to run the tests. `cd` to the `llnl_db_client` directory and run

```bash
$ LLNL_DB_PATH=/path/to/UCRL-MI-222502/westernus.wfdisc py.test
```

### Potential Pitfalls

* Essentially untested!
* This is not a proper CSS reader and should at one point probably replaced by
  something more solid, e.g. `obspy.io.css` or something else.
* The reponses are a bit tricky, except the RESP files which are easy. The
  others follow this logic (**this might very well be wrong!!!**):
	* SACPZ files:
		* The client currently assumes they correct to displacement in 10E-6 meters! The only justification I have is that this puts the amplitudes in the same ballpark as those corrected for with evalresp.
	* Funky PAZ files:
		* The client currently assumes they correct to displacement in 10E5 meters! The only justification I have is that this puts the amplitudes in the same ballpark as those corrected for with evalresp.

### Usage

First initialize a client:

```python
from llnl_db_client import LLNLDBClient

c = LLNLDBClient("./UCRL-MI-222502/westernus.wfdisc")
```

Print some basic information:


```python
>>> print(c)
LLNL Database 'westernus' (/Users/lion/temp/LLNL_Alaska/LLNL/UCRL-MI-222502)
	11176 waveform files
	147 stations
```

Get a list of all events with:

```python
>>> c.list_events()
[592897, 2021892, ...]
```

Get an ObsPy `Stream` object with all waveforms for an event:

```python
>>> st = c.get_waveforms_for_event(592897)
>>> print(st)
44 Trace(s) in Stream:

.ELK..BBE | 1988-10-13T13:58:01.143010Z - 1988-10-13T14:07:01.118010Z | 40.0 Hz, 21600 samples
...
(42 other traces)
...
.NEL..SHZ | 1988-10-13T13:59:31.660000Z - 1988-10-13T14:05:02.640000Z | 50.0 Hz, 16550 samples

[Use "print(Stream.__str__(extended=True))" to print all Traces]
```

Instrument correction.

```python
# Has to be done manually!
st.detrend("linear")
st.taper(max_percentage=0.05, type="hann")

# Water level, output, and pre_filt are the same as in Trace.remove_response()
c.remove_response(st, water_level=10, output="VEL", pre_filt=None)
```


Plot a map of all available stations:

```python
c.plot_stations()
```

![](doc/images/stations.png)
