"""
This code was written by Matthew Perry
https://github.com/mr-perry/sharad_pds_tools
"""


from os.path import isfile, getsize, basename
from os import access, R_OK, remove
from struct import iter_unpack
from numpy import fromfile, unique, zeros, log2, ceil
from datetime import datetime
from bitstring import BitArray, BitStream


#
# Random functions for read SHARAD PDS Data
#


def parse_sharad_pds_filename(inputFile):
    """
    Should work with either an LBL, *_A.DAT, or *_S.DAT
    :param inputFile:
    :return:
    """
    #
    # Determine to Bits per Sample
    #
    SSInstrMode = {'SS01': {'Mode': 'SS01', 'Presum': 32, 'BitsPerSample': 8, 'recLen': 3786},
                   'SS02': {'Mode': 'SS02', 'Presum': 28, 'BitsPerSample': 6, 'recLen': 2886},
                   'SS03': {'Mode': 'SS03', 'Presum': 16, 'BitsPerSample': 4, 'recLen': 1986},
                   'SS04': {'Mode': 'SS04', 'Presum': 8, 'BitsPerSample': 8, 'recLen': 3786},
                   'SS05': {'Mode': 'SS05', 'Presum': 4, 'BitsPerSample': 6, 'recLen': 2886},
                   'SS06': {'Mode': 'SS06', 'Presum': 2, 'BitsPerSample': 4, 'recLen': 1986},
                   'SS07': {'Mode': 'SS07', 'Presum': 1, 'BitsPerSample': 8, 'recLen': 3786},
                   'SS08': {'Mode': 'SS08', 'Presum': 32, 'BitsPerSample': 6, 'recLen': 2886},
                   'SS09': {'Mode': 'SS09', 'Presum': 28, 'BitsPerSample': 4, 'recLen': 1986},
                   'SS10': {'Mode': 'SS10', 'Presum': 16, 'BitsPerSample': 8, 'recLen': 3786},
                   'SS11': {'Mode': 'SS11', 'Presum': 8, 'BitsPerSample': 6, 'recLen': 2886},
                   'SS12': {'Mode': 'SS12', 'Presum': 4, 'BitsPerSample': 4, 'recLen': 1986},
                   'SS13': {'Mode': 'SS13', 'Presum': 2, 'BitsPerSample': 8, 'recLen': 3786},
                   'SS14': {'Mode': 'SS14', 'Presum': 1, 'BitsPerSample': 6, 'recLen': 2886},
                   'SS15': {'Mode': 'SS15', 'Presum': 32, 'BitsPerSample': 4, 'recLen': 1986},
                   'SS16': {'Mode': 'SS16', 'Presum': 28, 'BitsPerSample': 8, 'recLen': 3786},
                   'SS17': {'Mode': 'SS17', 'Presum': 16, 'BitsPerSample': 6, 'recLen': 2886},
                   'SS18': {'Mode': 'SS18', 'Presum': 8, 'BitsPerSample': 4, 'recLen': 1986},
                   'SS19': {'Mode': 'SS19', 'Presum': 4, 'BitsPerSample': 8, 'recLen': 3786},
                   'SS20': {'Mode': 'SS20', 'Presum': 2, 'BitsPerSample': 6, 'recLen': 2886},
                   'SS21': {'Mode': 'SS21', 'Presum': 1, 'BitsPerSample': 4, 'recLen': 1986},
                   }
    #
    # Receive only
    #
    ROInstrMode = {'RO01': {'Presum': 32, 'BitsPerSample': 8, 'recLen': 3786},
                   'RO02': {'Presum': 28, 'BitsPerSample': 6, 'recLen': 2886},
                   'RO03': {'Presum': 16, 'BitsPerSample': 4, 'recLen': 1986},
                   'RO04': {'Presum': 8, 'BitsPerSample': 8, 'recLen': 3786},
                   'RO05': {'Presum': 4, 'BitsPerSample': 6, 'recLen': 2886},
                   'RO06': {'Presum': 2, 'BitsPerSample': 4, 'recLen': 1986},
                   'RO07': {'Presum': 1, 'BitsPerSample': 8, 'recLen': 3786},
                   'RO08': {'Presum': 32, 'BitsPerSample': 6, 'recLen': 2886},
                   'RO09': {'Presum': 28, 'BitsPerSample': 4, 'recLen': 1986},
                   'RO10': {'Presum': 16, 'BitsPerSample': 8, 'recLen': 3786},
                   'RO11': {'Presum': 8, 'BitsPerSample': 6, 'recLen': 2886},
                   'RO12': {'Presum': 4, 'BitsPerSample': 4, 'recLen': 1986},
                   'RO13': {'Presum': 2, 'BitsPerSample': 8, 'recLen': 3786},
                   'RO14': {'Presum': 1, 'BitsPerSample': 6, 'recLen': 2886},
                   'RO15': {'Presum': 32, 'BitsPerSample': 4, 'recLen': 1986},
                   'RO16': {'Presum': 28, 'BitsPerSample': 8, 'recLen': 3786},
                   'RO17': {'Presum': 16, 'BitsPerSample': 6, 'recLen': 2886},
                   'RO18': {'Presum': 8, 'BitsPerSample': 4, 'recLen': 1986},
                   'RO19': {'Presum': 4, 'BitsPerSample': 8, 'recLen': 3786},
                   'RO20': {'Presum': 2, 'BitsPerSample': 6, 'recLen': 2886},
                   'RO21': {'Presum': 1, 'BitsPerSample': 4, 'recLen': 1986}
                   }
    #
    # Get the file basename
    #
    bname = basename(inputFile).split('_')
    prodType = bname[0].upper()
    TransID = bname[1]
    OSTLine = bname[2]
    OperMode = bname[3].upper()
    #
    # Deal with incorrectly formatted filenames
    #
    OperMode = OperMode[:2] + OperMode[2:].zfill(2) if len(OperMode) == 3 else OperMode
    #
    # Get the correct PRF value
    #
    PRF = sharad_pds_detPRF(bname[4])
    #
    # Extract information
    #
    if OperMode[0:2] == 'SS':
        OperMode = SSInstrMode[OperMode]
    elif OperMode[0:2] == 'RO':
        OperMode = ROInstrMode[OperMode]
    return prodType, TransID, OSTLine, OperMode, PRF


def sharad_pds_detPRF(val):
    PRF = {'335': 335.12,
           '350': 350.14,
           '387': 387.60,
           '670': 670.24,
           '700': 700.28,
           '775': 775.19,
           }
    return PRF[str(val)]


def load6bit(f, n, c):
    fmt = ['int:6'] * n
    a = zeros(int(n * c), int)
    with open(f, 'rb') as _f:
        for _i in range(c):
            s = BitStream(_f.read(2700))
            a[_i * n:(_i * n) + n] = s.unpack(fmt)
    return a


def decompressScienceData(sciData, compressionSelection, presum, BitsPerSample=None, SDI_BIT_FIELD=None):
    """
        This function decompresses the data based off page 8 of the
        SHALLOW RADAR EXPERIMENT DATA RECORD SOFTWARE INTERFACE SPECIFICATION.
        If the compression type is 'Static'
           U = C* ( 2**S / N )
             C is the Compressed Data
             S = L - R + 8
                L is base 2 log of N rounded UP to the nearest integer
                R is the bit resolution of the compressed values
             N is the number of pre-summed echoes
        If the compression type is 'dynamic'
          NOT WORKING YET
          U = C*2**S/N
            C is the compressed data
            N is the number of pre-summed echoes
            S = SDI for SDI <= 5
            S = SDI-6 for 5 < SDI <= 16
            S = SDI-16 for SDI > 16
              where SDI is the SDI_BIT_FIELD parameter from ANCILLARY DATA
    """

    if len(compressionSelection) == 1:
        if compressionSelection == 0:
            #
            # Static scaling
            #
            L = ceil(log2(presum))
            R = float(BitsPerSample)
            S = L - R + 8
            N = float(presum)
        elif compressionSelection == 1:
            SDI = unique(SDI_BIT_FIELD)
            if len(SDI) == 1:
                N = presum
                if SDI <= 5:
                    S = SDI
                elif 5 < SDI <= 16:
                    S = SDI - 6
                elif SDI > 16:
                    S = SDI - 16
        else:
            ValueError('Compression selection {} not understood'.format(compressionSelection))
        decom = (2 ** S) / N
        sciData = sciData * decom
        return sciData


def applyInstrumentResponse(manualGainControl, sciData):
        # CalVal = (VoltSat*2/(2**8-1)*10**((ANC.ManualGainControl(ind,1)-(Gain-ATx))/20)*10**3)
        VoltSat = 0.5  # Saturation Voltage (V)
        Gain = 88  # Receiver Maximum Gain (dB)
        ATx = 8.4  # TFE Rx Path Attenuation (dB)
        CalVal = (VoltSat * 2 / (2 ** 8 - 1) * 10 ** (
                (manualGainControl - (Gain - ATx)) / 20) * 10 ** 3)
        sciData = sciData * CalVal
        return sciData


#
# Functions for making empty data dictionaries
#


def make_sharad_edr_aux_dict(nCols):
    auxiliary_data = {'SCET_BLOCK_WHOLE': zeros(nCols, int),
                      'SCET_BLOCK_FRAC': zeros(nCols, float),
                      'EPHEMERIS_TIME': zeros(nCols, float),
                      'GEOMETRY_EPOCH': [],
                      'SOLAR_LONGITUDE': zeros(nCols, float),
                      'ORBIT_NUMBER': zeros(nCols, float),
                      'X_MARS_SC_POSITION_VECTOR': zeros(nCols, float),
                      'Y_MARS_SC_POSITION_VECTOR': zeros(nCols, float),
                      'Z_MARS_SC_POSITION_VECTOR': zeros(nCols, float),
                      'SPACECRAFT_ALTITUDE': zeros(nCols, float),
                      'SUB_SC_EAST_LONGITUDE': zeros(nCols, float),
                      'SUB_SC_PLANETOCENTRIC_LATITUDE': zeros(nCols, float),
                      'SUB_SC_PLANETOGRAPHIC_LATITUDE': zeros(nCols, float),
                      'X_MARS_SC_VELOCITY_VECTOR': zeros(nCols, float),
                      'Y_MARS_SC_VELOCITY_VECTOR': zeros(nCols, float),
                      'Z_MARS_SC_VELOCITY_VECTOR': zeros(nCols, float),
                      'MARS_SC_RADIAL_VELOCITY': zeros(nCols, float),
                      'MARS_SC_TANGENTIAL_VELOCITY': zeros(nCols, float),
                      'LOCAL_TRUE_SOLAR_TIME': zeros(nCols, float),
                      'SOLAR_ZENITH_ANGLE': zeros(nCols, float),
                      'SC_PITCH_ANGLE': zeros(nCols, float),
                      'SC_YAW_ANGLE': zeros(nCols, float),
                      'SC_ROLL_ANGLE': zeros(nCols, float),
                      'MRO_SAMX_INNER_GIMBAL_ANGLE': zeros(nCols, float),
                      'MRO_SAMX_OUTER_GIMBAL_ANGLE': zeros(nCols, float),
                      'MRO_SAPX_INNER_GIMBAL_ANGLE': zeros(nCols, float),
                      'MRO_SAPX_OUTER_GIMBAL_ANGLE': zeros(nCols, float),
                      'MRO_HGA_INNER_GIMBAL_ANGLE': zeros(nCols, float),
                      'MRO_HGA_OUTER_GIMBAL_ANGLE': zeros(nCols, float),
                      'DES_TEMP': zeros(nCols, float),
                      'DES_5V': zeros(nCols, float),
                      'DES_12V': zeros(nCols, float),
                      'DES_2V5': zeros(nCols, float),
                      'RX_TEMP': zeros(nCols, float),
                      'TX_TEMP': zeros(nCols, float),
                      'TX_LEV': zeros(nCols, float),
                      'TX_CURR': zeros(nCols, float),
                      'CORRUPTED_DATA_FLAG': zeros(nCols, int),
                      }
    return auxiliary_data


def make_sharad_edr_dict(nCols):
    sharad_edr_dict = {'SCET_BLOCK_WHOLE': zeros(nCols, int),
                       'SCET_BLOCK_FRAC': zeros(nCols, float),
                       'TLM_COUNTER': zeros(nCols, int),
                       'FMT_LENGTH': zeros(nCols, int),
                       'BIT_RESOLUTION': zeros(nCols, int),
                       'SCET_OST_WHOLE': zeros(nCols, int),
                       'SCET_OST_FRAC': zeros(nCols, float),
                       'OST_LINE_NUMBER': zeros(nCols, float),
                       'OST_LINE': {'PULSE_REPETITION_INTERVAL': zeros(nCols, int),
                                    'PHASE_COMPENSATION_TYPE': zeros(nCols, int),
                                    'DATA_TAKE_LENGTH': zeros(nCols, int),
                                    'OPERATIVE_MODE': zeros(nCols, int),
                                    'MANUAL_GAIN_CONTROL': zeros(nCols, int),
                                    'COMPRESSION_SELECTION': zeros(nCols, int),
                                    'CLOSED_LOOP_TRACKING': zeros(nCols, int),
                                    'TRACKING_DATA_STORAGE': zeros(nCols, int),
                                    'TRACKING_PRE_SUMMING': zeros(nCols, int),
                                    'TRACKING_LOGIC_SELECTION': zeros(nCols, int),
                                    'THRESHOLD_LOGIC_SELECTION': zeros(nCols, int),
                                    'SAMPLE_NUMBER': zeros(nCols, int),
                                    'ALPHA_BETA': zeros(nCols, int),
                                    'REFERENCE_BIT': zeros(nCols, int),
                                    'THRESHOLD': zeros(nCols, int),
                                    'THRESHOLD_INCREMENT': zeros(nCols, int),
                                    'INITIAL_ECHO_VALUE': zeros(nCols, int),
                                    'EXPECTED_ECHO_SHIFT': zeros(nCols, int),
                                    'WINDOW_LEFT_SHIFT': zeros(nCols, int),
                                    'WINDOW_RIGHT_SHIFT': zeros(nCols, int),
                                    },
                       'DATA_BLOCK_ID': zeros(nCols, int),
                       'SCIENCE_DATA_SOURCE_COUNTER': zeros(nCols, int),
                       'PACK_SEGMENTATION_AND_FPGA_STATUS': {'SCIENTIFIC_DATA_TYPE': zeros(nCols, int),
                                                             'SEGMENTATION_FLAG': zeros(nCols, int),
                                                             'DMA_ERROR': zeros(nCols, int),
                                                             'TC_OVERRUN': zeros(nCols, int),
                                                             'FIFO_FULL': zeros(nCols, int),
                                                             'TEST': zeros(nCols, int),
                                                             },
                       'DATA_BLOCK_FIRST_PRI': zeros(nCols, int),
                       'TIME_DATA_BLOCK_WHOLE': zeros(nCols, int),
                       'TIME_DATA_BLOCK_FRAC': zeros(nCols, float),
                       'SDI_BIT_FIELD': zeros(nCols, int),
                       'TIME_N': zeros(nCols, float),
                       'RADIUS_N': zeros(nCols, float),
                       'TANGENTIAL_VELOCITY_N': zeros(nCols, float),
                       'RADIAL_VELOCITY_N': zeros(nCols, float),
                       'TLP': zeros(nCols, float),
                       'TIME_WPF': zeros(nCols, float),
                       'DELTA_TIME': zeros(nCols, float),
                       'TLP_INTERPOLATE': zeros(nCols, float),
                       'RADIUS_INTERPOLATE': zeros(nCols, float),
                       'TANGENTIAL_VELOCITY_INTERPOLATE': zeros(nCols, float),
                       'RADIAL_VELOCITY_INTERPOLATE': zeros(nCols, float),
                       'END_TLP': zeros(nCols, float),
                       'S_COEFFS': zeros([8, nCols], float),
                       'C_COEFFS': zeros([7, nCols], float),
                       'SLOPE': zeros(nCols, float),
                       'TOPOGRAPHY': zeros(nCols, float),
                       'PHASE_COMPENSATION_STEP': zeros(nCols, float),
                       'RECEIVE_WINDOW_OPENING_TIME': zeros(nCols, float),
                       'RECEIVE_WINDOW_POSITION': zeros(nCols, int),
                       'SCIENCE_DATA': zeros([3600, nCols]),
                       }
    return sharad_edr_dict


def make_sharad_rdr_dict(nCols):
    sharad_rdr_dict = {'SCET_BLOCK_WHOLE': zeros(nCols, int),
                       'SCET_BLOCK_FRAC': zeros(nCols, float),
                       'TLM_COUNTER': zeros(nCols, int),
                       'FMT_LENGTH': zeros(nCols, int),
                       'SCET_OST_WHOLE': zeros(nCols, int),
                       'SCET_OST_FRAC': zeros(nCols, float),
                       'OST_LINE_NUMBER': zeros(nCols, int),
                       'PULSE_REPETITION_INTERVAL': zeros(nCols, int),
                       'PHASE_COMPENSATION_TYPE': zeros(nCols, int),
                       'DATA_TAKE_LENGTH': zeros(nCols, int),
                       'OPERATIVE_MODE': zeros(nCols, int),
                       'MANUAL_GAIN_CONTROL': zeros(nCols, int),
                       'COMPRESSION_SELECTION': zeros(nCols, bool),
                       'CLOSED_LOOP_TRACKING': zeros(nCols, bool),
                       'TRACKING_DATA_STORAGE': zeros(nCols, bool),
                       'TRACKING_PRE_SUMMING': zeros(nCols, int),
                       'TRACKING_LOGIC_SELECTION': zeros(nCols, int),
                       'THRESHOLD_LOGIC_SELECTION': zeros(nCols, int),
                       'SAMPLE_NUMBER': zeros(nCols, int),
                       'ALPHA_BETA': zeros(nCols, int),
                       'REFERENCE_BIT': zeros(nCols, int),
                       'THRESHOLD': zeros(nCols, int),
                       'THRESHOLD_INCREMENT': zeros(nCols, int),
                       'INITIAL_ECHO_VALUE': zeros(nCols, int),
                       'EXPECTED_ECHO_SHIFT': zeros(nCols, int),
                       'WINDOW_LEFT_SHIFT': zeros(nCols, int),
                       'WINDOW_RIGHT_SHIFT': zeros(nCols, int),
                       'DATA_BLOCK_ID': zeros(nCols, int),
                       'SCIENCE_DATA_SOURCE_COUNTER': zeros(nCols, int),
                       'SCIENTIFIC_DATA_TYPE': zeros(nCols, int),
                       'SEGMENTATION_FLAG': zeros(nCols, int),
                       'DMA_ERROR': zeros(nCols, int),
                       'TC_OVERRUN': zeros(nCols, int),
                       'FIFO_FULL': zeros(nCols, int),
                       'TEST': zeros(nCols, int),
                       'DATA_BLOCK_FIRST_PRI': zeros(nCols, int),
                       'TIME_DATA_BLOCK_WHOLE': zeros(nCols, int),
                       'TIME_DATA_BLOCK_FRAC': zeros(nCols, float),
                       'SDI_BIT_FIELD': zeros(nCols, int),
                       'TIME_N': zeros(nCols, float),
                       'RADIUS_N': zeros(nCols, float),
                       'TANGENTIAL_VELOCITY_N': zeros(nCols, float),
                       'RADIAL_VELOCITY_N': zeros(nCols, float),
                       'TLP': zeros(nCols, float),
                       'TIME_WPF': zeros(nCols, float),
                       'DELTA_TIME': zeros(nCols, float),
                       'TLP_INTERPOLATE': zeros(nCols, float),
                       'RADIUS_INTERPOLATE': zeros(nCols, float),
                       'TANGENTIAL_VELOCITY_INTERPOLATE': zeros(nCols, float),
                       'RADIAL_VELOCITY_INTERPOLATE': zeros(nCols, float),
                       'END_TLP': zeros(nCols, float),
                       'S_COEFFS': zeros([nCols, 8], float),
                       'C_COEFFS': zeros([nCols, 7], float),
                       'SLOPE': zeros(nCols, float),
                       'TOPOGRAPHY': zeros(nCols, float),
                       'PHASE_COMPENSATION_STEP': zeros(nCols, float),
                       'RECEIVE_WINDOW_OPENING_TIME': zeros(nCols, float),
                       'ANTENNA_RELATIVE_GAIN': zeros(nCols, float),
                       'ECHO_SAMPLES_REAL': zeros([667, nCols], float),
                       'ECHO_SAMPLES_IMAGINARY': zeros([667, nCols], float),
                       'N_PRE': zeros(nCols, int),
                       'BLOCK_NR': zeros(nCols, int),
                       'BLOCK_ROWS': zeros(nCols, int),
                       'DOPPLER_BW': zeros(nCols, float),
                       'DOPPLER_CENTROID': zeros(nCols, float),
                       'AZ_TIME_SPACING': zeros(nCols, float),
                       'AZ_RES': zeros(nCols, float),
                       'T_INT': zeros(nCols, float),
                       'AVG_TAN_VELOCITY': zeros(nCols, float),
                       'RANGE_SHIFT': zeros(nCols, int),
                       'EPHEMERIS_TIME': zeros(nCols, float),
                       'GEOMETRY_EPOCH': [],
                       'SOLAR_LONGITUDE': zeros(nCols, float),
                       'ORBIT_NUMBER': zeros(nCols, int),
                       'MARS_SC_POSITION_VECTOR': zeros([nCols, 3], float),
                       'SPACECRAFT_ALTITUDE': zeros(nCols, float),
                       'SUB_SC_EAST_LONGITUDE': zeros(nCols, float),
                       'SUB_SC_PLANETOCENTRIC_LATITUDE': zeros(nCols, float),
                       'SUB_SC_PLANETOGRAPHIC_LATITUDE': zeros(nCols, float),
                       'MARS_SC_VELOCITY_VECTOR': zeros([nCols, 3], float),
                       'MARS_SC_RADIAL_VELOCITY': zeros(nCols, float),
                       'MARS_SC_TANGENTIAL_VELOCITY': zeros(nCols, float),
                       'LOCAL_TRUE_SOLAR_TIME': zeros(nCols, float),
                       'SOLAR_ZENITH_ANGLE': zeros(nCols, float),
                       'SC_PITCH_ANGLE': zeros(nCols, float),
                       'SC_YAW_ANGLE': zeros(nCols, float),
                       'SC_ROLL_ANGLE': zeros(nCols, float),
                       'MRO_SAMX_INNER_GIMBAL_ANGLE': zeros(nCols, float),
                       'MRO_SAMX_OUTER_GIMBAL_ANGLE': zeros(nCols, float),
                       'MRO_SAPX_INNER_GIMBAL_ANGLE': zeros(nCols, float),
                       'MRO_SAPX_OUTER_GIMBAL_ANGLE': zeros(nCols, float),
                       'MRO_HGA_INNER_GIMBAL_ANGLE': zeros(nCols, float),
                       'MRO_HGA_OUTER_GIMBAL_ANGLE': zeros(nCols, float),
                       'DES_TEMP': zeros(nCols, float),
                       'DES_5V': zeros(nCols, float),
                       'DES_12V': zeros(nCols, float),
                       'DES_2V5': zeros(nCols, float),
                       'RX_TEMP': zeros(nCols, float),
                       'TX_TEMP': zeros(nCols, float),
                       'TX_LEV': zeros(nCols, float),
                       'TX_CURR': zeros(nCols, float),
                       'QUALITY_CODE': zeros(nCols, int)}
    return sharad_rdr_dict


#
# Functions for actually reading the data
#


def read_sharad_edr_aux(inputFile):
    """
    DOC STRING

    :param inputFile:
    :return:
    """
    if isfile(inputFile) and access(inputFile, R_OK):
        #
        # Before proceeding, parse the filename for specific observation details
        #
        prodType, TransID, OSTLine, OperMode, PRF = parse_sharad_pds_filename(inputFile)
        #
        # The number of record bytes for auxiliary files is data at 267
        #
        recordBytes = int(267)
        #
        # Use the file size divided by the recordBytes to calculate the number of data columns
        #
        nCols = int(getsize(inputFile)/recordBytes)
        #
        # The byte format of a column of data as described by the INDEX file
        #
        byteFormat = ">IHd23cdi23d8fh"
        #
        # Date string format
        #
        dformat = "%Y-%m-%dT%H:%M:%S.%f"
        #
        # Make empty dictionary for the data
        #
        auxiliary_data = make_sharad_edr_aux_dict(nCols)
        #
        # Read in the data
        #
        with open(inputFile, 'rb') as f:
            for frame in range(int(nCols)):
                f.seek(int(frame * recordBytes))
                column = f.read(recordBytes)
                a = iter_unpack(byteFormat, column)
                for item in a:
                    auxiliary_data['SCET_BLOCK_WHOLE'][frame] = item[0]
                    auxiliary_data['SCET_BLOCK_FRAC'][frame] = item[1] * 2E-16
                    auxiliary_data['EPHEMERIS_TIME'][frame] = item[2]
                    auxiliary_data['GEOMETRY_EPOCH'].append(datetime.strptime(
                        "".join([x.decode('UTF-8') for x in item[3:26]]),
                        dformat))
                    auxiliary_data['SOLAR_LONGITUDE'][frame] = item[26]
                    auxiliary_data['ORBIT_NUMBER'][frame] = item[27]
                    auxiliary_data['X_MARS_SC_POSITION_VECTOR'][frame] = item[28]
                    auxiliary_data['Y_MARS_SC_POSITION_VECTOR'][frame] = item[29]
                    auxiliary_data['Z_MARS_SC_POSITION_VECTOR'][frame] = item[30]
                    auxiliary_data['SPACECRAFT_ALTITUDE'][frame] = item[31]
                    auxiliary_data['SUB_SC_EAST_LONGITUDE'][frame] = item[32]
                    auxiliary_data['SUB_SC_PLANETOCENTRIC_LATITUDE'][frame] = item[33]
                    auxiliary_data['SUB_SC_PLANETOGRAPHIC_LATITUDE'][frame] = item[34]
                    auxiliary_data['X_MARS_SC_VELOCITY_VECTOR'][frame] = item[35]
                    auxiliary_data['Y_MARS_SC_VELOCITY_VECTOR'][frame] = item[36]
                    auxiliary_data['Z_MARS_SC_VELOCITY_VECTOR'][frame] = item[37]
                    auxiliary_data['MARS_SC_RADIAL_VELOCITY'][frame] = item[38]
                    auxiliary_data['MARS_SC_TANGENTIAL_VELOCITY'][frame] = item[39]
                    auxiliary_data['LOCAL_TRUE_SOLAR_TIME'][frame] = item[40]
                    auxiliary_data['SOLAR_ZENITH_ANGLE'][frame] = item[41]
                    auxiliary_data['SC_PITCH_ANGLE'][frame] = item[42]
                    auxiliary_data['SC_YAW_ANGLE'][frame] = item[43]
                    auxiliary_data['SC_ROLL_ANGLE'][frame] = item[44]
                    auxiliary_data['MRO_SAMX_INNER_GIMBAL_ANGLE'][frame] = item[45]
                    auxiliary_data['MRO_SAMX_OUTER_GIMBAL_ANGLE'][frame] = item[46]
                    auxiliary_data['MRO_SAPX_INNER_GIMBAL_ANGLE'][frame] = item[47]
                    auxiliary_data['MRO_SAPX_OUTER_GIMBAL_ANGLE'][frame] = item[48]
                    auxiliary_data['MRO_HGA_INNER_GIMBAL_ANGLE'][frame] = item[49]
                    auxiliary_data['MRO_HGA_OUTER_GIMBAL_ANGLE'][frame] = item[50]
                    auxiliary_data['DES_TEMP'][frame] = item[51]
                    auxiliary_data['DES_5V'][frame] = item[52]
                    auxiliary_data['DES_12V'][frame] = item[53]
                    auxiliary_data['DES_2V5'][frame] = item[54]
                    auxiliary_data['RX_TEMP'][frame] = item[55]
                    auxiliary_data['TX_TEMP'][frame] = item[56]
                    auxiliary_data['TX_LEV'][frame] = item[57]
                    auxiliary_data['TX_CURR'][frame] = item[58]
                    auxiliary_data['CORRUPTED_DATA_FLAG'][frame] = item[59]
        auxiliary_data['ProdType'] = prodType
        auxiliary_data['TransID'] = TransID
        auxiliary_data['OSTLine'] = OSTLine
        auxiliary_data['OperMode'] = OperMode
        auxiliary_data['PRF'] = PRF
        return auxiliary_data
    else:
        if not isfile(inputFile):
            raise FileNotFoundError('File not found.\n{}'.format(inputFile))
        if not access(inputFile, R_OK):
            raise PermissionError('File is not readable.\n{}'.format(inputFile))


def read_sharad_edr_sci(inputFile):
    """
    DOC STRING

    :param inputFile:
    :param readAnc:
    :param readSci:
    :return:
    """
    if isfile(inputFile) and access(inputFile, R_OK):
        #
        # Before proceeding, parse the filename for specific observation details
        #
        prodType, TransID, OSTLine, OperMode, PRF = parse_sharad_pds_filename(inputFile)
        #
        # Byte format for the science ancillary data
        #
        #byteFormat = ">IHIHHIHBB16sB3sH2sB3sIHHffffffffffff8f7fffffI"
        byteFormat = ">IHI2HIH2B16sB3sH2sB3sI2H31fI"
        #
        # Use the file size divided by the recordBytes to calculate the number of data columns
        #
        nCols = int(getsize(inputFile)/OperMode['recLen'])
        #
        # Make the empty dictionary
        #
        sharad_edr_dict = make_sharad_edr_dict(nCols)
        #
        # Temporary file for storing the science data
        #
        tmpSci = "sci.tmp"
        #
        # Open temporary file to store science data
        #
        _tmpSci = open(tmpSci, 'wb')
        #
        # Open science file and start reading
        #
        with open(inputFile, 'rb') as f:
            for frame in range(int(nCols)):
                f.seek(int(frame * OperMode['recLen']))
                column = f.read(OperMode['recLen'])
                #
                # Separate out the ancillary data
                #
                anc = column[:186]
                #
                # Write science data to temporary file
                #
                _tmpSci.write(column[186:])
                #
                # Unpack the ancillary data bytes
                #
                a = iter_unpack(byteFormat, anc)
                for item in a:
                    sharad_edr_dict['SCET_BLOCK_WHOLE'][frame] = item[0]
                    sharad_edr_dict['SCET_BLOCK_FRAC'][frame] = item[1] * 2E-16
                    sharad_edr_dict['TLM_COUNTER'][frame] = item[2]
                    sharad_edr_dict['FMT_LENGTH'][frame] = item[3]
                    if item[3] == 1972:
                        sharad_edr_dict['BIT_RESOLUTION'][frame] = 4
                    elif item[3] == 2872:
                        sharad_edr_dict['BIT_RESOLUTION'][frame] = 6
                    elif item[3] == 3772:
                        sharad_edr_dict['BIT_RESOLUTION'][frame] = 8
                    sharad_edr_dict['SCET_OST_WHOLE'][frame] = item[5]
                    sharad_edr_dict['SCET_OST_FRAC'][frame] = item[6] * 2E-16
                    sharad_edr_dict['OST_LINE_NUMBER'][frame] = item[8]
                    #
                    # DEAL WITH OST BITSTRING
                    #
                    ostBitString = BitArray(item[9])
                    pri = ostBitString[0:4].uint
                    if pri == 1:
                        sharad_edr_dict['OST_LINE']['PULSE_REPETITION_INTERVAL'][frame] = 1428
                    elif pri == 2:
                        sharad_edr_dict['OST_LINE']['PULSE_REPETITION_INTERVAL'][frame] = 1492
                    elif pri == 3:
                        sharad_edr_dict['OST_LINE']['PULSE_REPETITION_INTERVAL'][frame] = 1290
                    elif pri == 4:
                        sharad_edr_dict['OST_LINE']['PULSE_REPETITION_INTERVAL'][frame] = 2856
                    elif pri == 5:
                        sharad_edr_dict['OST_LINE']['PULSE_REPETITION_INTERVAL'][frame] = 2984
                    elif pri == 6:
                        sharad_edr_dict['OST_LINE']['PULSE_REPETITION_INTERVAL'][frame] = 2580
                    sharad_edr_dict['OST_LINE']['PHASE_COMPENSATION_TYPE'][frame] = ostBitString[4:8].uint
                    sharad_edr_dict['OST_LINE']['DATA_TAKE_LENGTH'][frame] = ostBitString[10:32].uint
                    sharad_edr_dict['OST_LINE']['OPERATIVE_MODE'][frame] = ostBitString[32:40].uint
                    sharad_edr_dict['OST_LINE']['MANUAL_GAIN_CONTROL'][frame] = ostBitString[40:48].uint
                    sharad_edr_dict['OST_LINE']['COMPRESSION_SELECTION'][frame] = ostBitString[48]
                    sharad_edr_dict['OST_LINE']['CLOSED_LOOP_TRACKING'][frame] = ostBitString[49]
                    sharad_edr_dict['OST_LINE']['TRACKING_DATA_STORAGE'][frame] = ostBitString[50]
                    sharad_edr_dict['OST_LINE']['TRACKING_PRE_SUMMING'][frame] = ostBitString[51:54].uint
                    sharad_edr_dict['OST_LINE']['TRACKING_LOGIC_SELECTION'][frame] = ostBitString[54]
                    sharad_edr_dict['OST_LINE']['TRACKING_LOGIC_SELECTION'][frame] = ostBitString[54]
                    sharad_edr_dict['OST_LINE']['THRESHOLD_LOGIC_SELECTION'][frame] = ostBitString[55]
                    sharad_edr_dict['OST_LINE']['SAMPLE_NUMBER'][frame] = ostBitString[56:60].uint
                    sharad_edr_dict['OST_LINE']['ALPHA_BETA'][frame] = ostBitString[61:63].uint
                    sharad_edr_dict['OST_LINE']['REFERENCE_BIT'][frame] = ostBitString[63]
                    sharad_edr_dict['OST_LINE']['THRESHOLD'][frame] = ostBitString[64:72].uint
                    sharad_edr_dict['OST_LINE']['THRESHOLD_INCREMENT'][frame] = ostBitString[72:80].uint
                    sharad_edr_dict['OST_LINE']['INITIAL_ECHO_VALUE'][frame] = ostBitString[84:87].uint
                    sharad_edr_dict['OST_LINE']['EXPECTED_ECHO_SHIFT'][frame] = ostBitString[87:90].uint
                    sharad_edr_dict['OST_LINE']['WINDOW_LEFT_SHIFT'][frame] = ostBitString[90:93].uint
                    sharad_edr_dict['OST_LINE']['WINDOW_RIGHT_SHIFT'][frame] = ostBitString[93:96].uint
                    #
                    # Done with OST Line; item[10] is a SPARE
                    #
                    sharad_edr_dict['DATA_BLOCK_ID'][frame] = BitArray(item[11]).uint
                    sharad_edr_dict['SCIENCE_DATA_SOURCE_COUNTER'][frame] = item[12]
                    #
                    # Deal with PACK_SEGMENTATION_AND_FPGA_STATUS bitstring
                    #
                    psfs = BitArray(item[13])
                    sharad_edr_dict['PACK_SEGMENTATION_AND_FPGA_STATUS']['SCIENTIFIC_DATA_TYPE'][frame] = psfs[0]
                    sharad_edr_dict['PACK_SEGMENTATION_AND_FPGA_STATUS']['SEGMENTATION_FLAG'][frame] = psfs[1:3].uint
                    sharad_edr_dict['PACK_SEGMENTATION_AND_FPGA_STATUS']['DMA_ERROR'][frame] = psfs[12]
                    sharad_edr_dict['PACK_SEGMENTATION_AND_FPGA_STATUS']['TC_OVERRUN'][frame] = psfs[13]
                    sharad_edr_dict['PACK_SEGMENTATION_AND_FPGA_STATUS']['FIFO_FULL'][frame] = psfs[14]
                    sharad_edr_dict['PACK_SEGMENTATION_AND_FPGA_STATUS']['TEST'][frame] = psfs[15]
                    #
                    # End PACK_SEGMENTATION_AND_FPGA_STATUS bitstring
                    # item[14] is SPARE
                    sharad_edr_dict['DATA_BLOCK_FIRST_PRI'][frame] = BitArray(item[15]).uint
                    sharad_edr_dict['TIME_DATA_BLOCK_WHOLE'][frame] = item[16]
                    sharad_edr_dict['TIME_DATA_BLOCK_FRAC'][frame] = item[17] * 2E-16
                    sharad_edr_dict['SDI_BIT_FIELD'][frame] = item[18]
                    sharad_edr_dict['TIME_N'][frame] = item[19]
                    sharad_edr_dict['RADIUS_N'][frame] = item[20]
                    sharad_edr_dict['TANGENTIAL_VELOCITY_N'][frame] = item[21]
                    sharad_edr_dict['RADIAL_VELOCITY_N'][frame] = item[22]
                    sharad_edr_dict['TLP'][frame] = item[23]
                    sharad_edr_dict['TIME_WPF'][frame] = item[24]
                    sharad_edr_dict['DELTA_TIME'][frame] = item[25]
                    sharad_edr_dict['TLP_INTERPOLATE'][frame] = item[26]
                    sharad_edr_dict['RADIUS_INTERPOLATE'][frame] = item[27]
                    sharad_edr_dict['TANGENTIAL_VELOCITY_INTERPOLATE'][frame] = item[28]
                    sharad_edr_dict['RADIAL_VELOCITY_INTERPOLATE'][frame] = item[29]
                    sharad_edr_dict['END_TLP'][frame] = item[30]
                    sharad_edr_dict['S_COEFFS'][:, frame] = item[31:39]
                    sharad_edr_dict['C_COEFFS'][:, frame] = item[39:46]
                    sharad_edr_dict['SLOPE'][frame] = item[46]
                    sharad_edr_dict['TOPOGRAPHY'][frame] = item[47]
                    sharad_edr_dict['PHASE_COMPENSATION_STEP'][frame] = item[48]
                    sharad_edr_dict['RECEIVE_WINDOW_OPENING_TIME'][frame] = item[49]
                    sharad_edr_dict['RECEIVE_WINDOW_POSITION'][frame] = item[50]
        _tmpSci.close()
        #
        # Now, deal with the actual science data
        #
        if sharad_edr_dict['BIT_RESOLUTION'][-1] == 8:
            #sharad_edr_dict['SCIENCE_DATA'] = fromfile(tmpSci, dtype='int8').reshape([3600, nCols])
            sharad_edr_dict['SCIENCE_DATA'] = fromfile(tmpSci, dtype='int8').reshape([nCols, 3600]).T
        elif sharad_edr_dict['BIT_RESOLUTION'][-1] == 4:
            #sharad_edr_dict['SCIENCE_DATA'] = fromfile(tmpSci, dtype='int4').reshape([3600, nCols])
            sharad_edr_dict['SCIENCE_DATA'] = fromfile(tmpSci, dtype='int4').reshape([nCols, 3600]).T
        elif sharad_edr_dict['BIT_RESOLUTION'][-1] == 6:
            #sharad_edr_dict['SCIENCE_DATA'] = load6bit(tmpSci, 3600, nCols)
            sharad_edr_dict['SCIENCE_DATA'] = load6bit(tmpSci, 3600, nCols)
        remove('sci.tmp')
        #
        # Decompress the Science Data
        #
        compressionSelection = unique(sharad_edr_dict['OST_LINE']['COMPRESSION_SELECTION'])
        sharad_edr_dict['SCIENCE_DATA'] = decompressScienceData(sharad_edr_dict['SCIENCE_DATA'],
                                                                compressionSelection,
                                                                OperMode['Presum'],
                                                                BitsPerSample=OperMode['BitsPerSample'],
                                                                SDI_BIT_FIELD=sharad_edr_dict['SDI_BIT_FIELD'])
        sharad_edr_dict['ProdType'] = prodType
        sharad_edr_dict['TransID'] = TransID
        sharad_edr_dict['OSTLine'] = OSTLine
        sharad_edr_dict['OperMode'] = OperMode
        sharad_edr_dict['PRF'] = PRF
        sharad_edr_dict['Sample_Rate'] = 0.0375e-6
        sharad_edr_dict['Latency'] = 11.12e-6
        #
        # Apply Instrument Response; should be moved elsewhere...
        #
        sharad_edr_dict['SCIENCE_DATA'] = applyInstrumentResponse(sharad_edr_dict['OST_LINE']['MANUAL_GAIN_CONTROL'],
                                                                  sharad_edr_dict['SCIENCE_DATA'])
        return sharad_edr_dict
    else:
        if not isfile(inputFile):
            raise FileNotFoundError('File not found.\n{}'.format(inputFile))
        if not access(inputFile, R_OK):
            raise PermissionError('File is not readable.\n{}'.format(inputFile))


def read_sharad_rdr(inputFile):
    if isfile(inputFile) and access(inputFile, R_OK):
        #
        # Before proceeding, parse the filename for specific observation details
        #
        prodType, TransID, OSTLine, OperMode, PRF = parse_sharad_pds_filename(inputFile)
        #
        # Record Length, this is data
        #
        recLen = 5822
        #
        # Byte format for the science ancillary data
        #
        byteFormat = "<IHIHIH3BI3B2?12BIH6B2I2H12f8f7f5f667f667f3H6fhd23cdi3d4d3d13d8fB"
        #
        # Use the file size divided by the recordBytes to calculate the number of data columns
        #
        nCols = int(getsize(inputFile) / recLen)
        #
        # Make the empty dictionary
        #
        sharad_rdr_dict = make_sharad_rdr_dict(nCols)
        #
        #
        with open(inputFile, 'rb') as f:
            for frame in range(int(nCols)):
                f.seek(frame * recLen)
                column = f.read(recLen)
                a = iter_unpack(byteFormat, column)
                for item in a:
                    sharad_rdr_dict['SCET_BLOCK_WHOLE'][frame] = item[0]
                    sharad_rdr_dict['SCET_BLOCK_FRAC'][frame] = item[1] * 1e-16
                    sharad_rdr_dict['TLM_COUNTER'][frame] = item[2]
                    sharad_rdr_dict['FMT_LENGTH'][frame] = item[3]
                    sharad_rdr_dict['SCET_OST_WHOLE'][frame] = item[4]
                    sharad_rdr_dict['SCET_OST_FRAC'][frame] = item[5] * 1e-16
                    sharad_rdr_dict['OST_LINE_NUMBER'][frame] = item[6]
                    sharad_rdr_dict['PULSE_REPETITION_INTERVAL'][frame] = item[7]
                    sharad_rdr_dict['PHASE_COMPENSATION_TYPE'][frame] = item[8]
                    sharad_rdr_dict['DATA_TAKE_LENGTH'][frame] = item[9]
                    sharad_rdr_dict['OPERATIVE_MODE'][frame] = item[10]
                    sharad_rdr_dict['MANUAL_GAIN_CONTROL'][frame] = item[11]
                    sharad_rdr_dict['COMPRESSION_SELECTION'][frame] = item[12]
                    sharad_rdr_dict['CLOSED_LOOP_TRACKING'][frame] = item[13]
                    sharad_rdr_dict['TRACKING_DATA_STORAGE'][frame] = item[14]
                    sharad_rdr_dict['TRACKING_PRE_SUMMING'][frame] = item[15]
                    sharad_rdr_dict['TRACKING_LOGIC_SELECTION'][frame] = item[16]
                    sharad_rdr_dict['THRESHOLD_LOGIC_SELECTION'][frame] = item[17]
                    sharad_rdr_dict['SAMPLE_NUMBER'][frame] = item[18]
                    sharad_rdr_dict['ALPHA_BETA'][frame] = item[19]
                    sharad_rdr_dict['REFERENCE_BIT'][frame] = item[20]
                    sharad_rdr_dict['THRESHOLD'][frame] = item[21]
                    sharad_rdr_dict['THRESHOLD_INCREMENT'][frame] = item[22]
                    sharad_rdr_dict['INITIAL_ECHO_VALUE'][frame] = item[23]
                    sharad_rdr_dict['EXPECTED_ECHO_SHIFT'][frame] = item[24]
                    sharad_rdr_dict['WINDOW_LEFT_SHIFT'][frame] = item[25]
                    sharad_rdr_dict['WINDOW_RIGHT_SHIFT'][frame] = item[26]
                    sharad_rdr_dict['DATA_BLOCK_ID'][frame] = item[27]
                    sharad_rdr_dict['SCIENCE_DATA_SOURCE_COUNTER'][frame] = item[28]
                    sharad_rdr_dict['SCIENTIFIC_DATA_TYPE'][frame] = item[29]
                    sharad_rdr_dict['SEGMENTATION_FLAG'][frame] = item[30]
                    sharad_rdr_dict['DMA_ERROR'][frame] = item[31]
                    sharad_rdr_dict['TC_OVERRUN'][frame] = item[32]
                    sharad_rdr_dict['FIFO_FULL'][frame] = item[33]
                    sharad_rdr_dict['TEST'][frame] = item[34]
                    sharad_rdr_dict['DATA_BLOCK_FIRST_PRI'][frame] = item[35]
                    sharad_rdr_dict['TIME_DATA_BLOCK_WHOLE'][frame] = item[36]
                    sharad_rdr_dict['TIME_DATA_BLOCK_FRAC'][frame] = item[37]
                    sharad_rdr_dict['SDI_BIT_FIELD'][frame] = item[38]
                    sharad_rdr_dict['TIME_N'][frame] = item[39]
                    sharad_rdr_dict['RADIUS_N'][frame] = item[40]
                    sharad_rdr_dict['TANGENTIAL_VELOCITY_N'][frame] = item[41]
                    sharad_rdr_dict['RADIAL_VELOCITY_N'][frame] = item[42]
                    sharad_rdr_dict['TLP'][frame] = item[43]
                    sharad_rdr_dict['TIME_WPF'][frame] = item[44]
                    sharad_rdr_dict['DELTA_TIME'][frame] = item[45]
                    sharad_rdr_dict['TLP_INTERPOLATE'][frame] = item[46]
                    sharad_rdr_dict['RADIUS_INTERPOLATE'][frame] = item[47]
                    sharad_rdr_dict['TANGENTIAL_VELOCITY_INTERPOLATE'][frame] = item[48]
                    sharad_rdr_dict['RADIAL_VELOCITY_INTERPOLATE'][frame] = item[49]
                    sharad_rdr_dict['END_TLP'][frame] = item[50]
                    sharad_rdr_dict['S_COEFFS'][frame, :] = item[51:59]
                    sharad_rdr_dict['C_COEFFS'][frame, :] = item[59:66]
                    sharad_rdr_dict['SLOPE'][frame] = item[66]
                    sharad_rdr_dict['TOPOGRAPHY'][frame] = item[67]
                    sharad_rdr_dict['PHASE_COMPENSATION_STEP'][frame] = item[68]
                    sharad_rdr_dict['RECEIVE_WINDOW_OPENING_TIME'][frame] = item[69]
                    sharad_rdr_dict['ANTENNA_RELATIVE_GAIN'][frame] = item[70]
                    sharad_rdr_dict['ECHO_SAMPLES_REAL'][:, frame] = item[71:738]
                    sharad_rdr_dict['ECHO_SAMPLES_IMAGINARY'][:, frame] = item[738:1405]
                    sharad_rdr_dict['N_PRE'][frame] = item[1405]
                    sharad_rdr_dict['BLOCK_NR'][frame] = item[1406]
                    sharad_rdr_dict['BLOCK_ROWS'][frame] = item[1407]
                    sharad_rdr_dict['DOPPLER_BW'][frame] = item[1408]
                    sharad_rdr_dict['DOPPLER_CENTROID'][frame] = item[1409]
                    sharad_rdr_dict['AZ_TIME_SPACING'][frame] = item[1410]
                    sharad_rdr_dict['AZ_RES'][frame] = item[1411]
                    sharad_rdr_dict['T_INT'][frame] = item[1412]
                    sharad_rdr_dict['AVG_TAN_VELOCITY'][frame] = item[1413]
                    sharad_rdr_dict['RANGE_SHIFT'][frame] = item[1414]
                    sharad_rdr_dict['EPHEMERIS_TIME'][frame] = item[1415]
                    sharad_rdr_dict['GEOMETRY_EPOCH'].append(''.join([x.decode('utf-8') for x in item[1416:1439]]))
                    sharad_rdr_dict['SOLAR_LONGITUDE'][frame] = item[1439]
                    sharad_rdr_dict['ORBIT_NUMBER'][frame] = item[1440]
                    sharad_rdr_dict['MARS_SC_POSITION_VECTOR'][frame, :] = item[1441:1444]
                    sharad_rdr_dict['SPACECRAFT_ALTITUDE'][frame] = item[1444]
                    sharad_rdr_dict['SUB_SC_EAST_LONGITUDE'][frame] = item[1445]
                    sharad_rdr_dict['SUB_SC_PLANETOCENTRIC_LATITUDE'][frame] = item[1446]
                    sharad_rdr_dict['SUB_SC_PLANETOGRAPHIC_LATITUDE'][frame] = item[1447]
                    sharad_rdr_dict['MARS_SC_VELOCITY_VECTOR'][frame, :] = item[1448:1451]
                    sharad_rdr_dict['MARS_SC_RADIAL_VELOCITY'][frame] = item[1451]
                    sharad_rdr_dict['MARS_SC_TANGENTIAL_VELOCITY'][frame] = item[1452]
                    sharad_rdr_dict['LOCAL_TRUE_SOLAR_TIME'][frame] = item[1453]
                    sharad_rdr_dict['SOLAR_ZENITH_ANGLE'][frame] = item[1454]
                    sharad_rdr_dict['SC_PITCH_ANGLE'][frame] = item[1455]
                    sharad_rdr_dict['SC_YAW_ANGLE'][frame] = item[1456]
                    sharad_rdr_dict['SC_ROLL_ANGLE'][frame] = item[1457]
                    sharad_rdr_dict['MRO_SAMX_INNER_GIMBAL_ANGLE'][frame] = item[1458]
                    sharad_rdr_dict['MRO_SAMX_OUTER_GIMBAL_ANGLE'][frame] = item[1459]
                    sharad_rdr_dict['MRO_SAPX_INNER_GIMBAL_ANGLE'][frame] = item[1460]
                    sharad_rdr_dict['MRO_SAPX_OUTER_GIMBAL_ANGLE'][frame] = item[1461]
                    sharad_rdr_dict['MRO_HGA_INNER_GIMBAL_ANGLE'][frame] = item[1462]
                    sharad_rdr_dict['MRO_HGA_OUTER_GIMBAL_ANGLE'][frame] = item[1463]
                    sharad_rdr_dict['DES_TEMP'][frame] = item[1464]
                    sharad_rdr_dict['DES_5V'][frame] = item[1465]
                    sharad_rdr_dict['DES_12V'][frame] = item[1466]
                    sharad_rdr_dict['DES_2V5'][frame] = item[1467]
                    sharad_rdr_dict['RX_TEMP'][frame] = item[1468]
                    sharad_rdr_dict['TX_TEMP'][frame] = item[1469]
                    sharad_rdr_dict['TX_LEV'][frame] = item[1470]
                    sharad_rdr_dict['TX_CURR'][frame] = item[1471]
                    sharad_rdr_dict['QUALITY_CODE'][frame] = item[1472]
        sharad_rdr_dict['ProdType'] = prodType
        sharad_rdr_dict['TransID'] = TransID
        sharad_rdr_dict['OSTLine'] = OSTLine
        sharad_rdr_dict['OperMode'] = OperMode
        sharad_rdr_dict['PRF'] = PRF
        return sharad_rdr_dict
    else:
        if not isfile(inputFile):
            raise FileNotFoundError('File not found.\n{}'.format(inputFile))
        if not access(inputFile, R_OK):
            raise PermissionError('File is not readable.\n{}'.format(inputFile))


def read_sharad_us_rdr(inputFile):
    """
    Reads in SHARAD US Science Data.
    Works with both CO-SHARPS .IMG files and PDS SHARAD US RDR

    :param str inputFile: The full path to SHARAD US Science Files
    :return: 2D Numpy Array [3600, n] of the SHARAD Science Data
    """
    if isfile(inputFile) and access(inputFile, R_OK):
        return fromfile(inputFile, "<f").reshape((3600, -1))
    else:
        if not isfile(inputFile):
            raise FileNotFoundError('CO-SHARPS Science file not found.\n{}'.format(inputFile))
        if not access(inputFile, R_OK):
            raise PermissionError('CO-SHARPS Science file is not readable.\n{}'.format(inputFile))


def read_sharad_us_geom(inputFile):
    dateFMT = "%Y-%m-%dT%H:%M:%S.%f"
    if isfile(inputFile) and access(inputFile, R_OK):
        geomDict = {'RADARGRAM_COLUMN': [], 'TIME': [], 'LATITUDE': [], 'LONGITUDE': [],
                    'MARS_RADIUS': [], 'SPACECRAFT_RADIUS': [], 'RADIAL_VELOCITY': [],
                    'TANGENTIAL_VELOCITY': [], 'SZA': [], 'PHASE/1.0E16': [],}
        with open(inputFile, 'r') as f:
            content = f.readlines()
            content = [line.strip() for line in content]
        for _l, line in enumerate(content):
            csvs = line.split(',')
            geomDict['RADARGRAM_COLUMN'].append(int(csvs[0]))
            try:
                geomDict['TIME'].append(datetime.strptime(csvs[1], dateFMT))
            except Exception:
                print('[WARNING] GEOM File contains TIME Error at line {}\n{}'.format(_l, csvs[1]))
                dt = geomDict['TIME'][-1] - geomDict['TIME'][-2]
                geomDict['TIME'].append(geomDict['TIME'][-1] + dt)
            geomDict['LATITUDE'].append(float(csvs[2]))
            geomDict['LONGITUDE'].append(float(csvs[3]))
            geomDict['MARS_RADIUS'].append(float(csvs[4]))
            geomDict['SPACECRAFT_RADIUS'].append(float(csvs[5]))
            geomDict['RADIAL_VELOCITY'].append(float(csvs[6]))
            geomDict['TANGENTIAL_VELOCITY'].append(float(csvs[7]))
            geomDict['SZA'].append(float(csvs[8]))
            geomDict['PHASE/1.0E16'].append(float(csvs[9]))
        return geomDict
    else:
        if not isfile(inputFile):
            raise FileNotFoundError('File not found.\n{}'.format(inputFile))
        if not access(inputFile, R_OK):
            raise PermissionError('File is not readable.\n{}'.format(inputFile))


#
# Wrapper Function
#


def read_sharad_pds(inputFile):
    if isfile(inputFile) and access(inputFile, R_OK):
        basenameParts = basename(inputFile).split('_')
        if basenameParts[0].upper() == 'E':
            #
            # SHARAD EDR Files
            #
            if basenameParts[-1].upper() == 'A.DAT':
                return read_sharad_edr_aux(inputFile)
            elif basenameParts[-1].upper() == 'S.DAT':
                return read_sharad_edr_sci(inputFile)
            else:
                raise RuntimeError('File type not understood.\n{}'.format(inputFile))
        elif basenameParts[0].upper() == 'R':
            #
            # SHARAD RDR File
            #
            if basenameParts[-1].upper() == 'A.DAT':
                return read_sharad_rdr(inputFile)
            else:
                raise RuntimeError('File type not understood.\n{}'.format(inputFile))
        elif basenameParts[0].upper() == 'S':
            #
            # SHARAD US RDR Files
            #
            if basenameParts[-1].upper() == 'RGRAM.IMG':
                #
                # SHARAD US RDR File
                #
                return read_sharad_us_rdr(inputFile)
            elif basenameParts[-1].upper() == 'GEOM.TAB':
                #
                # SHARAD US RDR GEOM File
                #
                return read_sharad_us_geom(inputFile)
            else:
                raise RuntimeError('File type not understood.\n{}'.format(inputFile))
        else:
            raise RuntimeError('File type not understood.\n{}'.format(inputFile))
    else:
        if not isfile(inputFile):
            raise FileNotFoundError('File not found.\n{}'.format(inputFile))
        if not access(inputFile, R_OK):
            raise PermissionError('File is not readable.\n{}'.format(inputFile))


#
# Example on how to read various data types
#
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from numpy import abs, log10, percentile
    #
    # Reading in Test EDR Data both Auxilary and Science Data. The RAW SCIENCE DATA is contained in
    # testSci['SCIENCE_DATA']
    #
    testAux = read_sharad_pds('e_0323403_001_ss19_700_a_a.dat')
    testSci = read_sharad_pds('e_0323403_001_ss19_700_a_s.dat')
    del testAux, testSci
    #
    # Read in Test RDR
    #
    testRDR = read_sharad_pds('r_1002101_001_ss05_700_a.dat')
    toPlot = testRDR['ECHO_SAMPLES_REAL'] + 1j*testRDR['ECHO_SAMPLES_IMAGINARY']
    toPlot = 20*log10(abs(toPlot))
    vmin = percentile(toPlot, 50)
    vmax = percentile(toPlot, 99)
    fig, ax = plt.subplots(1, 1, figsize=[12, 6])
    c00 = ax.imshow(toPlot, aspect='auto', cmap='gray', vmin=vmin, vmax=vmax)
    ax.set_title('RDR Example -- RX WIN Datum', fontsize=16)
    ax.set_xlabel('Radargram Frame', fontsize=16)
    ax.set_ylabel('Frame Sample', fontsize=16)
    cbar00 = plt.colorbar(c00, ax=ax)
    cbar00.set_label('Power (dB)', fontsize=16)
    plt.tight_layout()
    plt.show()
    #
    # Test read US RDR and GEOM File
    #
    testUSRDR = read_sharad_pds('s_03662702_rgram.img')
    testGEOM = read_sharad_pds('s_03662702_geom.tab')
    toPlot = 10*log10(abs(testUSRDR))
    vmin = percentile(toPlot, 50)
    vmax = percentile(toPlot, 99.99)
    fig, ax = plt.subplots(1, 1, figsize=[12, 6])
    c00 = ax.imshow(toPlot, aspect='auto', cmap='gray', vmin=vmin, vmax=vmax)
    ax.set_title('US RDR Example -- FPB Datum', fontsize=16)
    ax.set_xlabel('Radargram Frame', fontsize=16)
    ax.set_ylabel('Frame Sample', fontsize=16)
    cbar00 = plt.colorbar(c00, ax=ax)
    cbar00.set_label('Power (dB)', fontsize=16)
    plt.tight_layout()
    plt.show()
