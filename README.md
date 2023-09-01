# pyudv

pyudv is a python package allowing to read and process UDV data. The main modules are:

- `pyudv.read_mfprof`, which helps you read binary files from Met-Flow UDVs
- `pyudv.probes`, which helps you plot the various probes, calculate the intersection points, etc ..
- `pyudv.velocity`, which helps you reconstruct 2D velocity fields from the various probes measurements
- `pyudv.amplitude`, which helps you infer concentration from amplitude measurements

> [!WARNING]  
> Although tested, this is still in development, use with caution. Feedbacks welcome.

### Installation

#### User only

Using `pip3 install --upgrade https://github.com/cgadal-pythonpackages/pyudv/tarball/master`

#### If code or development

- clone the repository, e.g. `git clone https://github.com/cgadal-pythonpackages/pyudv`
- `cd pyudv && pip3 install -e ./` will install in editable mode.