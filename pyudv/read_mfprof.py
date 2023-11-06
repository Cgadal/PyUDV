"""
This module allows you to read `.mfprof` files from Met-Flow UDVs, and extract the contained data.

"""

import netCDF4
import numpy as np

Types = {
    "Frequency": np.int32,
    "StartChannel": np.float64,
    "ChannelDistance": np.float64,
    "ChannelWidth": np.float64,
    "MaximumDepth": np.float64,
    "SoundSpeed": np.int32,
    "Angle": int,
    "GainStart": np.int32,
    "GainEnd": np.int32,
    "Voltage": np.int32,
    "Iterations": np.int32,
    "NoiseLevel": np.int32,
    "CyclesPerPulse": np.int32,
    "TriggerMode": np.int32,
    "TriggerModeName": str,
    "ProfileLength": np.int32,
    "ProfilesPerBlock": np.int32,
    "Blocks": np.int32,
    "AmplitudeStored": int,
    "DoNotStoreDoppler": int,
    "RawDataMin": int,
    "RawDataMax": int,
    "RawDataRange": np.int32,
    "AmplDataMin": int,
    "AmplDataMax": int,
    "VelocityInterpretingMode": int,
    "UserSampleTime": int,
    "SampleTime": np.int32,
    "UseMultiplexer": int,
    "FlowMapping": int,
    "FirstValidChannel": np.int32,
    "LastValidChannel": np.int32,
    "FlowRateType": int,
    "PeriodEnhOffset": np.int32,
    "PeriodEnhPeriod": np.int32,
    "PeriodEnhNCycles": np.int32,
    "Comment": str,
    "MeasurementProtocol": str,
    "NumberOfCycles": np.int32,
    "CycleDelay": np.int32,
    "Version": np.int32,
    "Table": str,
    "CalcOffsetPulsesCount": np.int32,
    "MultiplexerUsingType": str,
}

Base_units = {
    "Frequency": "Hz",
    "StartChannel": "mm",
    "ChannelDistance": "mm",
    "ChannelWidth": "mm",
    "MaximumDepth": "mm",
    "SoundSpeed": "m/s",
    "Angle": "int",
    "GainStart": "int",
    "GainEnd": "int",
    "Voltage": "V",
    "Iterations": "int",
    "NoiseLevel": "int",
    "CyclesPerPulse": "int",
    "TriggerMode": "index",
    "TriggerModeName": "text",
    "ProfileLength": "nChannels",
    "ProfilesPerBlock": "nSamples",
    "Blocks": "int",
    "AmplitudeStored": "bool",
    "DoNotStoreDoppler": "-",
    "RawDataMin": "int",
    "RawDataMax": "int",
    "RawDataRange": "int",
    "AmplDataMin": "int",
    "AmplDataMax": "int",
    "VelocityInterpretingMode": "int",
    "UserSampleTime": "bool",
    # 'SampleTime': 'int',
    "SampleTime": "ms",
    "UseMultiplexer": "bool",
    "FlowMapping": "bool",
    "FirstValidChannel": "int",
    "LastValidChannel": "int",
    "FlowRateType": "bool",
    "PeriodEnhOffset": "int",
    "PeriodEnhPeriod": "int",
    "PeriodEnhNCycles": "int",
    "Comment": "text",
    "MeasurementProtocol": "text",
    "NumberOfCycles": "int",
    "CycleDelay": "int",
    "Version": "int",
    "Table": "int",
    # switched from np.int32 to str due to problem loading files with multiplexing
    "CalcOffsetPulsesCount": "int",
    "MultiplexerUsingType": "int",
}

Absolute_gains = {
    0.5: {
        3: 2.17,
        4: 4.41,
        5: 8.82,
        6: 16.67,
        7: 33.33,
        8: 60.00,
        9: 150.00
    },
    1: {
        3: 2.17,
        4: 4.41,
        5: 8.82,
        6: 16.67,
        7: 33.33,
        8: 60.00,
        9: 150.00
    },
    2: {
        3: 0.91,
        4: 1.76,
        5: 3.41,
        6: 6.67,
        7: 15.00,
        8: 25.00,
        9: 60.00
    },
    4: {
        3: 0.91,
        4: 1.76,
        5: 3.41,
        6: 6.67,
        7: 15.00,
        8: 25.00,
        9: 60.00
    },
    8: {
        3: 0.65,
        4: 1.36,
        5: 2.80,
        6: 5.26,
        7: 11.11,
        8: 23.08,
        9: 42.86
    },
}


def read_mfprof(fileName: str,
                SI_units: bool = True,
                convert_time: bool = True):
    """Read .mfprof binary files of the Met-Flow UDV. This is mostly a direct python translation of the matlab script given by Met-Flow.

    Parameters
    ----------
    fileName : str
        Path to the .mfprof file.
    SI_units : bool
        If `True`, convert units to the international system (the default is True).
    convert_time : bool
        If True, convert the time vector into seconds (the default is True).

    Returns
    -------
    Data : dict
        Dictionnary with the data stored in the mprof files. Available keys are: 'transducer', 'profileTime', 'DopplerData', 'AmplitudeData', 'DistanceAlongBeam'
    Parameters : dict
        Dictionnary with the parameters used in the UVP software when sampling the data. See UVP documentation for detail. Keys are: 'Frequency', 'StartChannel', 'ChannelDistance', 'ChannelWidth', 'MaximumDepth', 'SoundSpeed', 'Angle', 'GainStart', 'GainEnd', 'Voltage', 'Iterations', 'NoiseLevel', 'CyclesPerPulse', 'TriggerMode', 'TriggerModeName', 'ProfileLength', 'ProfilesPerBlock', 'Blocks', 'AmplitudeStored', 'DoNotStoreDoppler', 'RawDataMin', 'RawDataMax', 'RawDataRange', 'AmplDataMin', 'AmplDataMax', 'VelocityInterpretingMode', 'UserSampleTime', 'SampleTime', 'UseMultiplexer', 'FlowMapping', 'FirstValidChannel', 'LastValidChannel', 'FlowRateType', 'PeriodEnhOffset', 'PeriodEnhPeriod', 'PeriodEnhNCycles', 'Comment', 'MeasurementProtocol', 'NumberOfCycles', 'CycleDelay', 'Version', 'Table'
    Infos: dict
        Dictionnary with some informations stored during sampling of the data. Keys are 'Signum', 'measParamsOffset', 'nProfiles', 'reserved1', 'flags', 'recordSize', 'nChannels', 'reserved2', 'startTime'.
    Units: dict
        Dictionnary with the units of the variables stored in the other dictionnaries. Keys are the name of the variables.

    """
    Infos = {}
    Parameters = {}
    Data = {}
    with open(fileName, mode="rb") as file:  # b is important -> binary
        ############################################################################
        # ######################## File HEADER #####################################
        ############################################################################
        # Met-Flow 64-symbol file type signature
        Infos["Signum"] = np.fromfile(file, np.uint8, 64)  # char 64
        # 64-bit number giving position of beginning of Measurement Parameters
        Infos["measParamsOffset"] = np.fromfile(file, np.int64, 1)
        # Number of profiles saved in the file
        Infos["nProfiles"] = np.fromfile(file, np.int32, 1)
        # Always 0 parameter
        Infos["reserved1"] = np.fromfile(file, np.int32, 1)
        # Flag. If this is set (=1) then Profiles also contain Amplitude Infos
        Infos["flags"] = np.fromfile(file, np.int32, 1)
        # Size of Infos in a single Profile
        Infos["recordSize"] = np.fromfile(file, np.int32, 1)
        # Number of Infos channels in a single Profile
        Infos["nChannels"] = np.fromfile(file, np.int32, 1)
        # Always 0 parameter
        Infos["reserved2"] = np.fromfile(file, np.int32, 1)
        # Time of measurement of the start(first profile) in the format FILETIME (Win32) structure.
        Infos["startTime"] = np.fromfile(file, np.int64, 1)

        ############################################################################
        # ######### PROFILES - Read contents of PROFILES BLOCK #####################
        ############################################################################

        # Profile information (in 0 to 255, 8-bit)values follow immediately after
        # Profile header. Profiles consists of nProfiles Doppler Infos and if flags = 1
        # then profile also contain Amplitude Infos. Each block consists of:
        # Profile Header, Doppler Infos and Amplitude Infos. This information must be
        # read nProfiles times.

        # First, we allocate memory space for faster file reading
        Data["transducer"] = np.zeros((Infos["nProfiles"][0], ))
        Data["profileTime"] = np.zeros((Infos["nProfiles"][0], ))
        Data["DopplerData"] = np.zeros(
            (Infos["nChannels"][0], Infos["nProfiles"][0]))
        Data["AmplitudeData"] = np.zeros(
            (Infos["nChannels"][0], Infos["nProfiles"][0]))

        # READ CONTENTS FROM BLOCK nProfiles TIMES
        for i in range(Infos["nProfiles"][0]):
            # PROFILE HEADER BLOCK
            # Always 0 not used
            _ = np.fromfile(file, np.int32, 1)  # status
            # Data['transducer'] number -reads which Data['transducer'] Infos comes from
            Data["transducer"][i] = np.fromfile(file, np.int32, 1)
            # Time of measurement of the profile in 100 ns increments
            Data["profileTime"][i] = np.fromfile(file, np.int64, 1)
            # DOPPLER Infos BLOCK
            # Doppler Infos block contains nChannels Int16 Doppler values from UVP
            Data["DopplerData"][:, i] = np.fromfile(file, np.int16,
                                                    Infos["nChannels"][0])
            # AMPLITUDE Infos BLOCK
            # Amplitude Infos block contains nChannels Int16 Amplitude values from UVP
            if Infos["flags"]:  # Read Infos if it exists
                Data["AmplitudeData"][:,
                                      i] = np.fromfile(file, np.int16,
                                                       Infos["nChannels"][0])

        ############################################################################
        # ##### UVP PARAMETER  - Read the UVP Parameters into a Matlab struct ######
        ############################################################################

        # UVP Measurement Parameters are stored at measParamsOffset position. Measurement
        # Parameters are saved in text form. Values of parameters are saved in a form:
        # parameter name=parameter value. Each parameter is given in one line. In case of
        # text values including more lines, at the end of line (with the exception of the last
        # line) a symbol '\' is added. E.g. comment including two lines:
        # "An example of forwards-backwards oscillating 'piston' flow
        # Piston action starts only after Profile 110" will be saved as:
        # Comment=An example of forwards-backwards oscillating 'piston' flow.\
        # Piston action starts only after Profile 110.
        # Measurement Parameters consist of two blocks:
        # 1. UVP Parameters, and 2. and Multiplexer Parameters if present

        # READ CONTENTS FROM UVP PARAMETER BLOCK
        # Find position of Measurement Parameters start location
        file.seek(Infos["measParamsOffset"][0], 0)
        _ = np.fromfile(file, np.uint8, 15)  # measParamsposition
        Parameters = {}

        # READ CONTENTS FROM UVP PARAMETER BLOCK LINE BY LINE
        par = ""
        while par != "Table":
            line = file.readline().rstrip()  # Read new line
            if (len(line) > 0) & (b"=" in line):  # skip empty lines
                # print(line)
                tp = line.decode("utf-8").split("=")[0]
                par = line.decode("utf-8").split("=")[0]
                if len(tp) > 1:
                    Parameters[par] = line.decode("utf-8").split("=")[1]
                else:
                    Parameters[par] = []

    # converting time in seconds - dt = 100 ns
    if convert_time:
        Data["profileTime"] = Data["profileTime"] * 100 * 1e-9
    # Formatting parameters
    Units = Base_units.copy()
    for key in Parameters.keys():
        # Assigning good type to parameters
        Parameters[key] = Types[key](Parameters[key])
        if SI_units & (Base_units[key] == "mm"):  # Converting mm to m
            Parameters[key] = 1e-3 * Parameters[key]
            Units[key] = "m"
        if SI_units & (Base_units[key] == "ms"):  # Converting ms to s
            Parameters[key] = 1e-3 * Parameters[key]
            Units[key] = "s"

    Data["DistanceAlongBeam"] = (np.arange(Data["AmplitudeData"].shape[0]) *
                                 Parameters["ChannelDistance"] +
                                 Parameters["StartChannel"])
    return Data, Parameters, Infos, Units


# def filetime_to_dt(ft):
#     us = int(ft) // 10
#     return dt.datetime(1601, 1, 1) + dt.timedelta(microseconds=us)


def velocity_from_UVPdata(raw_data,
                          SoundSpeed,
                          MaximumDepth,
                          Angle: float,
                          Frequency: float,
                          Nbytes: int = 8):
    """Calculate velocity from UDV raw data. Units must be checked outside of this function - no conversion is done here.

    Parameters
    ----------
    raw_data : float, numpy array
        Raw data coming out of the UDV.
    SoundSpeed : float
        Sound velocity.
    MaximumDepth : float
        Maximum depth measurement, as defined by the Pulse Repetition Frequency.
    Angle : float
        Angle between the probe and the flow direction.
    Frequency : float
        Frequency of the probe.
    Nbytes : int
        Number of bytes over which the raw data are coded (the default is 8).

    Returns
    -------
    float, numpy array
        Velocity field

    """
    PRF = SoundSpeed / (2*MaximumDepth)  # Pulse Repetiion Frequency
    N_DU = 2**Nbytes  # Number of possible velocities coded on 'Nbytes'-bit
    Doppler_coeff = PRF / N_DU  # Doppler coefficient
    #
    Fd = raw_data * Doppler_coeff  # Doppler frequencies
    Angle_correction = 1 / np.sin(np.pi * Angle / 180)  # angle correction
    #
    # return radial velocity
    return Fd * SoundSpeed * Angle_correction / (2*Frequency)


def velocity_from_mfprof_reading(Data: dict,
                                 Parameters: dict,
                                 Nbytes: int = 8):
    """Calculate velocity from `Data` and `Parameters` dictionnaries as output of :func:`read_mfprof <pyudv.read_mfprof.read_mfprof>`,
    using the function :func:`velocity_from_UVPdata <pyudv.read_mfprof.velocity_from_UVPdata>`.

    Parameters
    ----------
    Data : dict
        Data dictionnary coming from :func:`read_mfprof <pyudv.read_mfprof.read_mfprof>`.
    Parameters : dict
        Parameters dictionnary coming from :func:`read_mfprof <pyudv.read_mfprof.read_mfprof>`.
    Nbytes : int
        Number of bytes over which the raw data are coded (the default is 8).

    Returns
    -------
    float, numpy array
        Velocity field

    """
    raw_data = Data["DopplerData"]
    SoundSpeed, MaximumDepth, Angle, Frequency = (
        Parameters["SoundSpeed"],
        Parameters["MaximumDepth"],
        Parameters["Angle"],
        Parameters["Frequency"],
    )
    return velocity_from_UVPdata(raw_data,
                                 SoundSpeed,
                                 MaximumDepth,
                                 Angle,
                                 Frequency,
                                 Nbytes=Nbytes)


def amplitude_from_UVPdata(
    raw_data,
    z,
    GainStart,
    GainEnd,
    zend: float,
    zstart: float = 0.37 * 1e-3,
    Nbytes: int = 14,
    deltaV: float = 5,
):
    """Correct the raw amplitude data as outut of the UVP-DUO.

    Parameters
    ----------
    raw_data : scalar, array
        raw echo signal before demodulation outputted by the UVP.
    z : scalar, array
        position corresponding to `raw_data`.
    GainStart : scalar
        start absolute gain.
    GainEnd : scalar
        end absolute gain.
    zend : float
        Maximum measurable distance by the UVP.
    zstart : float
        minimum measurable distance by the UVP (the default is 0.00037)).
    Nbytes : int
        Number of bytes over which the raw data are coded (the default is 14).
    deltaV : int
        Volt range corresponding to the number of bytes (the default is 5).

    Returns
    -------
    scalar, array
        Unamplified/corrected echo signal.

    """
    return (raw_data * (1/GainStart) * (GainStart / GainEnd)**((z-zstart) /
                                                               (zend-zstart)) *
            deltaV / 2**Nbytes)


def amplitude_from_mfprof_reading(Data: dict,
                                  Parameters: dict,
                                  Nbytes: int = 14,
                                  deltaV: float = 5):
    """Calculate the unamplified/corrected echo signal from `Data` and `Parameters` dictionnaries as output of :func:`read_mfprof <pyudv.read_mfprof.read_mfprof>`,
    using the function :func:`amplitude_from_UVPdata <pyudv.read_mfprof.amplitude_from_UVPdata>`.

    Parameters
    ----------
    Data : dict
        Data dictionnary coming from :func:`read_mfprof <pyudv.read_mfprof.read_mfprof>`.
    Parameters : dict
        Parameters dictionnary coming from :func:`read_mfprof <pyudv.read_mfprof.read_mfprof>`.
    Nbytes : int
        Number of bytes over which the raw data are coded (the default is 14).
    deltaV : int
        Volt range corresponding to the number of bytes (the default is 5).

    Returns
    -------
    scalar, array
        Unamplified/corrected echo signal.

    """
    raw_data, z = Data["AmplitudeData"], Data["DistanceAlongBeam"][:, None]
    zend = Parameters["MaximumDepth"]
    F0 = Parameters["Frequency"] / 1e6  # frequency in Mhz
    GainStart, GainEnd = (
        Absolute_gains[F0][Parameters["GainStart"]],
        Absolute_gains[F0][Parameters["GainEnd"]],
    )
    return amplitude_from_UVPdata(raw_data,
                                  z,
                                  GainStart,
                                  GainEnd,
                                  zend,
                                  Nbytes=Nbytes,
                                  deltaV=deltaV)


# #### Writing functions


def write_dictionnary(dico: dict, file: str):
    """Write parameter dictionnary to a .txt file.
    Each line of this file is a dictionnary key, followed by the corresponding entry, directly converted using :func:`str <str>`.

    Parameters
    ----------
    dico : dict
        input dictionnary
    file : str
        output txt file.

    Returns
    -------

        Nothing.

    """
    with open(file, "w") as f:
        for key in sorted(dico.keys()):
            line = key + ": " + str(dico[key]) + "\n"
            f.write(line)


def _create_variable(
    netcdf_group,
    name,
    data,
    dimensions=None,
    std=None,
    unit=None,
    comments=None,
    type="float64",
):
    if dimensions is not None:
        var = netcdf_group.createVariable(name, type, (dimensions))
    else:
        var = netcdf_group.createVariable(name, type)
    var[:] = data
    if std is not None:
        var.std = std
    if unit is not None:
        var.unit = unit
    if comments is not None:
        var.comments = comments


def mfprof_to_netcdf(input_mfprof: str,
                     output_netcdf: str | object,
                     add_attr: dict = None,
                     cut_zeros: bool = True):
    """
    Convert a .mfprof file to a netCDF file.

    Parameters
    ----------
    input_mfprof : str
        input .mfprof file path
    output_netcdf : str | object
        output netCDF file path, or netCDF4 object (Dataset, Group) in which the data will be stored.
    add_attr : dict, optional
        dictionnary containing attributes to be added to the netCDF file, by default None. For example, `add_attr = {name: 'John Smith', team: 'the best'}`
    cut_zeros : bool, optional
        If True, cut the measurements data and remove all the trailing 0 added when the measurements are stopped bedore the attributed time, by default True

    Examples
    --------
    >>> mfprof_to_netcdf('input.mfprof', 'output.nc', add_attr = {name: 'John Smith', team: 'the best'}, cut_zeros=True)
    """
    if isinstance(output_netcdf, str):
        newfile = netCDF4.Dataset(output_netcdf, "w", format="NETCDF4")
    else:
        newfile = output_netcdf
    #
    Data, Parameters, Info, Units = read_mfprof(input_mfprof)
    amplitude_data = amplitude_from_mfprof_reading(Data, Parameters)
    velocity_data = velocity_from_mfprof_reading(Data, Parameters)
    #
    tmax = (np.argwhere(Data["transducer"] == 0).squeeze()[0] if
            ((Data["transducer"] == 0).any() & cut_zeros) else -1)
    # ### fill Data
    newfile.createDimension("time", Data["profileTime"][:tmax].size)
    newfile.createDimension("z", Data["DistanceAlongBeam"].size)
    #
    for key, var in Data.items():
        if key in ["DistanceAlongBeam"]:
            dims = "z"
        elif key in ["profileTime", "transducer"]:
            dims = "time"
        else:
            dims = ("z", "time")
        #
        _create_variable(
            newfile,
            key,
            var[..., :tmax] if key != "DistanceAlongBeam" else var,
            dimensions=dims,
            unit=Units[key] if key in Units.keys() else None,
        )
    _create_variable(
        newfile,
        "amplitude_data",
        amplitude_data[..., :tmax],
        dimensions=("z", "time"),
        unit="V",
    )
    _create_variable(
        newfile,
        "velocity_data",
        velocity_data[..., :tmax],
        dimensions=("z", "time"),
        unit="m/s",
    )
    # ### fill Parameters and Info
    for key, var in Parameters.items():
        setattr(newfile, key, var)
    for key, var in Info.items():
        setattr(newfile, key, var)
    # ### fill additional attributes
    if add_attr is not None:
        for key, var in add_attr.items():
            setattr(newfile, key, var)
    #
    if cut_zeros:
        newfile.cut_zeros = "True. File has been shortened while converted to netcdf as there was trailing unused zeros."
        newfile.NumberOfCyclesShort = Data["profileTime"][:tmax].size
    #
    if isinstance(output_netcdf, str):
        newfile.close()
