from h5py import File
import numpy
from doy_to_ymd import datetime_to_timestamp, timestamp_to_datetime

def save_to_hdf(output_directory,
                rfilevel,
                time_datetime,
                ffreq,
                flag,
                distribution_observation_numebeam_azimuth, bins_distribution_observation_numebeam_azimuth,
                distribution_observation_numebeam_elevation, bins_distribution_observation_numebeam_elevation,
                distribution_observation_numebeam_frequency, bins_distribution_observation_numebeam_frequency,
                distribution_observation_anabeam_UT, bins_distribution_observation_anabeam_UT,
                distribution_observation_anabeam_LT, bins_distribution_observation_anabeam_LT,
                distribution_flag_data_universaltime_all_ana_beam, distribution_observations_data_universaltime_all_ana_beam,
                distribution_flag_data_localtime_all_ana_beam,distribution_observations_data_localtime_all_ana_beam,
                distribution_flag_data_azimuth_all_ana_beam,distribution_observations_data_azimuth_all_ana_beam,
                distribution_flag_data_elevation_all_ana_beam,distribution_observations_data_elevation_all_ana_beam,
                distribution_flag_data_elevation_all_nume_beam,distribution_observations_data_elevation_all_nume_beam,
                distribution_flag_data_azimuth_all_nume_beam,distribution_observations_data_azimuth_all_nume_beam,
                distribution_flag_data_elaz_all_nume_beam,distribution_observations_data_elaz_all_nume_beam
                ):
        """
        Saves the data to disk as a HDF5 file.
        """
        output_file = File(
            output_directory+'rfistat_distributions_rfilevel'+str(rfilevel)+'.hdf5', 'w'
        )
        #output_file.attrs.create('Description', 'RFI flag of NenuFAR observations')
        #output_file.attrs.create('rfilevel', rfilevel)
        #output_file.attrs.create('number of 6 minutes observations', len())
        output_file.create_dataset('distribution_observation_numebeam_azimuth',
                                   data = distribution_observation_numebeam_azimuth)
        output_file.create_dataset('bins_distribution_observation_numebeam_azimuth',
                                   data = bins_distribution_observation_numebeam_azimuth)
        output_file.create_dataset('distribution_observation_numebeam_elevation',
                                   data = distribution_observation_numebeam_elevation)
        output_file.create_dataset('bins_distribution_observation_numebeam_elevation',
                                   data = bins_distribution_observation_numebeam_elevation)
        output_file.create_dataset('distribution_observation_numebeam_frequency',
                                   data = distribution_observation_numebeam_frequency)
        output_file.create_dataset('bins_distribution_observation_numebeam_frequency',
                                   data = bins_distribution_observation_numebeam_frequency)
        output_file.create_dataset('distribution_flag_data_universaltime_all_ana_beam',
                                   data=distribution_flag_data_universaltime_all_ana_beam)
        output_file.create_dataset('distribution_flag_data_localtime_all_ana_beam',
                                   data=distribution_flag_data_localtime_all_ana_beam)
        output_file.create_dataset('distribution_flag_data_azimuth_all_ana_beam',
                                   data=distribution_flag_data_azimuth_all_ana_beam)
        output_file.create_dataset('distribution_flag_data_elevation_all_ana_beam',
                                   data=distribution_flag_data_elevation_all_ana_beam)
        output_file.create_dataset('distribution_flag_data_elevation_all_nume_beam',
                                   data=distribution_flag_data_elevation_all_nume_beam)
        output_file.create_dataset('distribution_flag_data_azimuth_all_nume_beam',
                                   data=distribution_flag_data_azimuth_all_nume_beam)
        output_file.create_dataset('distribution_flag_data_elaz_all_nume_beam',
                                   data=distribution_flag_data_elaz_all_nume_beam)
        
        output_file.create_dataset('distribution_observations_data_universaltime_all_ana_beam',
                                   data=distribution_observations_data_universaltime_all_ana_beam)
        output_file.create_dataset('distribution_observations_data_localtime_all_ana_beam',
                                   data=distribution_observations_data_localtime_all_ana_beam)
        output_file.create_dataset('distribution_observations_data_azimuth_all_ana_beam',
                                   data=distribution_observations_data_azimuth_all_ana_beam)
        output_file.create_dataset('distribution_observations_data_elevation_all_ana_beam',
                                   data=distribution_observations_data_elevation_all_ana_beam)
        output_file.create_dataset('distribution_observations_data_elevation_all_nume_beam',
                                   data=distribution_observations_data_elevation_all_nume_beam)
        output_file.create_dataset('distribution_observations_data_azimuth_all_nume_beam',
                                   data=distribution_observations_data_azimuth_all_nume_beam)
        output_file.create_dataset('distribution_observations_data_elaz_all_nume_beam',
                                   data=distribution_observations_data_elaz_all_nume_beam)

        output_file.create_dataset('distribution_observation_anabeam_UT',
                                   data=distribution_observation_anabeam_UT)
        output_file.create_dataset('bins_distribution_observation_anabeam_UT',
                                   data=bins_distribution_observation_anabeam_UT)
        output_file.create_dataset('distribution_observation_anabeam_LT',
                                   data=distribution_observation_anabeam_LT)
        output_file.create_dataset('bins_distribution_observation_anabeam_LT',
                                   data=bins_distribution_observation_anabeam_LT)
        time_timestamp = datetime_to_timestamp(time_datetime)
        output_file.create_dataset('Time', data = time_timestamp)
        output_file.create_dataset('Frequency', data=ffreq)
        output_file['Frequency'].attrs.create('units', 'MHz')
        output_file.create_dataset('flag', data = flag)

        output_file.close()


def read_hdf5_file(input_directory, rfilevel):
        
        with File(input_directory+'rfistat_distributions_rfilevel'+str(rfilevel)+'.hdf5', 'r') as file_hdf5:
                ffreq = numpy.array(file_hdf5['Frequency'])
                flag = numpy.array(file_hdf5['flag'])
                time_timestamp = numpy.array(file_hdf5['Time'])
                time_datetime = timestamp_to_datetime(time_timestamp)
                distribution_observation_numebeam_azimuth = numpy.array(file_hdf5['distribution_observation_numebeam_azimuth'])
                bins_distribution_observation_numebeam_azimuth = numpy.array(file_hdf5['bins_distribution_observation_numebeam_azimuth'])
                distribution_observation_numebeam_elevation = numpy.array(file_hdf5['distribution_observation_numebeam_elevation'])
                bins_distribution_observation_numebeam_elevation = numpy.array(file_hdf5['bins_distribution_observation_numebeam_elevation'])
                distribution_observation_numebeam_frequency = numpy.array(file_hdf5['distribution_observation_numebeam_frequency'])
                bins_distribution_observation_numebeam_frequency = numpy.array(file_hdf5['bins_distribution_observation_numebeam_frequency'])
                distribution_observation_anabeam_UT = numpy.array(file_hdf5['distribution_observation_anabeam_UT'])
                bins_distribution_observation_anabeam_UT = numpy.array(file_hdf5['bins_distribution_observation_anabeam_UT'])
                distribution_observation_anabeam_LT = numpy.array(file_hdf5['distribution_observation_anabeam_LT'])
                bins_distribution_observation_anabeam_LT = numpy.array(file_hdf5['bins_distribution_observation_anabeam_LT'])
                distribution_flag_data_universaltime_all_ana_beam = numpy.array(file_hdf5['distribution_flag_data_universaltime_all_ana_beam'])
                distribution_flag_data_localtime_all_ana_beam = numpy.array(file_hdf5['distribution_flag_data_localtime_all_ana_beam'])
                distribution_flag_data_azimuth_all_ana_beam = numpy.array(file_hdf5['distribution_flag_data_azimuth_all_ana_beam'])
                distribution_flag_data_elevation_all_ana_beam = numpy.array(file_hdf5['distribution_flag_data_elevation_all_ana_beam'])
                distribution_flag_data_elevation_all_nume_beam = numpy.array(file_hdf5['distribution_flag_data_elevation_all_nume_beam'])
                distribution_flag_data_azimuth_all_nume_beam = numpy.array(file_hdf5['distribution_flag_data_azimuth_all_nume_beam'])
                distribution_flag_data_elaz_all_nume_beam = numpy.array(file_hdf5['distribution_flag_data_elaz_all_nume_beam'])
                distribution_observations_data_universaltime_all_ana_beam = numpy.array(file_hdf5['distribution_observations_data_universaltime_all_ana_beam'])
                distribution_observations_data_localtime_all_ana_beam = numpy.array(file_hdf5['distribution_observations_data_localtime_all_ana_beam'])
                distribution_observations_data_azimuth_all_ana_beam = numpy.array(file_hdf5['distribution_observations_data_azimuth_all_ana_beam'])
                distribution_observations_data_elevation_all_ana_beam = numpy.array(file_hdf5['distribution_observations_data_elevation_all_ana_beam'])
                distribution_observations_data_elevation_all_nume_beam = numpy.array(file_hdf5['distribution_observations_data_elevation_all_nume_beam'])
                distribution_observations_data_azimuth_all_nume_beam = numpy.array(file_hdf5['distribution_observations_data_azimuth_all_nume_beam'])
                distribution_observations_data_elaz_all_nume_beam = numpy.array(file_hdf5['distribution_observations_data_elaz_all_nume_beam'])

        return(time_datetime, ffreq, flag,
               distribution_observation_numebeam_azimuth, bins_distribution_observation_numebeam_azimuth,
               distribution_observation_numebeam_elevation, bins_distribution_observation_numebeam_elevation,
               distribution_observation_numebeam_frequency, bins_distribution_observation_numebeam_frequency,
               distribution_observation_anabeam_UT, bins_distribution_observation_anabeam_UT,
               distribution_observation_anabeam_LT, bins_distribution_observation_anabeam_LT,
                distribution_flag_data_universaltime_all_ana_beam, distribution_observations_data_universaltime_all_ana_beam,
                distribution_flag_data_localtime_all_ana_beam,distribution_observations_data_localtime_all_ana_beam,
                distribution_flag_data_azimuth_all_ana_beam,distribution_observations_data_azimuth_all_ana_beam,
                distribution_flag_data_elevation_all_ana_beam,distribution_observations_data_elevation_all_ana_beam,
                distribution_flag_data_elevation_all_nume_beam,distribution_observations_data_elevation_all_nume_beam,
                distribution_flag_data_azimuth_all_nume_beam,distribution_observations_data_azimuth_all_nume_beam,
                distribution_flag_data_elaz_all_nume_beam,distribution_observations_data_elaz_all_nume_beam
         )