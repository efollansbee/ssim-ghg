import rdata
import numpy as np
from netCDF4 import Dataset
import os, pickle, yaml, time
from collections import defaultdict

class Timer(object):
    indent_levels = [] # data to be shared across all instances to keep track of indentation level
    def __init__(self, msg, **kwargs):
        self.msg = msg
        self.print_msg = kwargs['print'] if 'print' in kwargs else True
        self.addendum = {}

    def __enter__(self):
        self.t1 = time.time()
        self.indent_levels.append(1)
        return self.addendum

    def __exit__(self, ex_type, ex_value, ex_traceback):
        t2 = time.time()
        dt = t2-self.t1
        num_indents = len(self.indent_levels)
        if self.print_msg:
            if 'prefix' in self.addendum:
                print_line = "  "*num_indents + "%s %s %s"%(self.addendum['prefix'], self.msg, self.format_dt(dt))
            else:
                print_line = "  "*num_indents + "%s %s"%(self.msg, self.format_dt(dt))
            if 'postfix' in self.addendum:
                print_line = '%s %s'%(print_line, self.addendum['postfix'])
            print(print_line)
        _ = self.indent_levels.pop(0)

    def format_dt(self, dt):
        if dt > 60.0:
            minutes = dt // 60
            seconds = dt % 60.0
            return "%im %.2fs"%(minutes, seconds)
        else:
            return "%.2fs"%dt

class Paths(object):

    def __init__(self):
        super(Paths, self).__init__()
        # setting to be stored in ../site_settings.yml
        with open('../site_settings.yml', 'r') as fid:
            conf = yaml.safe_load(fid)
        self.data_root = conf['global_paths']['input_folder']
        self.output_root = conf['global_paths']['output_folder']
        try:
            self.figure_font = conf['plotting']['font']
        except KeyError:
            self.figure_font = 'Sans Serif'

        self.jacobi_rda = os.path.join(self.data_root, 'jacobians/trunc_full_jacob_032624_with_dimnames_unit_pulse_4x5_mask_hour_timestamping.rda')
        self.jacobi_nc = os.path.join(self.data_root, 'jacobians/trunc_full_jacob_032624_with_dimnames_unit_pulse_4x5_mask.nc')
        self.obs_rda = os.path.join(self.data_root, 'obs/obs_catalog_041724_unit_pulse_hour_timestamp_witherrors.rda')
        self.obs_nc = os.path.join(self.data_root, 'obs/obs_catalog_041724_unit_pulse_hour_timestamp_witherrors.nc')
        # self.obs_rda = os.path.join(self.data_root, 'obs/obs_catalog_042424_unit_pulse_hour_timestamp_witherrors_withdates.rda')
        # self.obs_nc = os.path.join(self.data_root, 'obs/obs_catalog_042424_unit_pulse_hour_timestamp_witherrors_withdates.nc')
        self.obs_by_dataset = os.path.join(self.data_root, 'obs/obs_by_dataset.pickle')
        self.bg_rda = os.path.join(self.data_root, 'jacobians/jacob_bgd_060524.rda')
        self.bg_nc = os.path.join(self.data_root, 'jacobians/jacob_bgd_060524.nc')

        self.region_aggregates = {
            'North America': ['North American Boreal', 'North American Temperate'],
            'South America': ['South American Tropical', 'South American Temperate'],
            'Africa': ['Northern Africa', 'Southern Africa'],
            'Eurasia': ['Eurasian Boreal', 'Eurasian Temperate', 'Tropical Asia'],
            'Pacific Ocean': ['North Pacific Temperate', 'West Pacific Tropical', 'East Pacific Tropical', 'South Pacific Temperate'],
            'Atlantic Ocean': ['North Atlantic Temperate', 'Atlantic Tropical', 'South Atlantic Temperate'],
            'Indian Ocean': ['Indian Tropical', 'South Indian Temperate'],
            'Global Land': ['North American Boreal', 'North American Temperate', 'South American Tropical', 'South American Temperate',
                'Northern Africa', 'Southern Africa', 'Eurasian Boreal', 'Eurasian Temperate', 'Tropical Asia', 'Australia', 'Europe'],
            'Global Ocean': ['North Pacific Temperate', 'West Pacific Tropical', 'East Pacific Tropical', 'South Pacific Temperate',
                'Northern Ocean', 'North Atlantic Temperate', 'Atlantic Tropical', 'South Atlantic Temperate', 'Southern Ocean',
                'Indian Tropical', 'South Indian Temperate'],
            }
        self.region_aggregates['Globe'] = self.region_aggregates['Global Land'] + self.region_aggregates['Global Ocean']

class Convert_RDA(Paths):

    def __init__(self, *args, **kwargs):
        super(Convert_RDA, self).__init__()
        self.error_scaling = { # IS data errors are in mol/mol, need to multiply by 10^6 to get ppm
            'IS': 1.0E6,
            }

    def convert_bg(self):
        with Timer("Read background from RDA in "):
            rda_data = rdata.read_rda(self.bg_rda)
        bg = rda_data['jacob_bgd']
        nobs, ncomp = bg.shape
        with Timer('Wrote background to netcdf in '):
            with Dataset(self.bg_nc, 'w') as fid:
                fid.createDimension('n_obs', nobs)
                fid.createDimension('n_comp', ncomp)
                v = fid.createVariable('BG', bg.dtype, ('n_obs','n_comp'), zlib=True, shuffle=True, complevel=5)
                v[:] = bg

    def convert_jacobian(self):
        with Timer('Read Jacobian from RDA in '):
            rda_data = rdata.read_rda(self.jacobi_rda)
        J = rda_data['jacob'].values
        nobs, nstate = J.shape
        with Timer('Wrote Jacobian to netcdf in '):
            with Dataset(self.jacobi_nc, 'w') as fid:
                fid.createDimension('n_obs', nobs)
                fid.createDimension('n_state', nstate)
                v = fid.createVariable('Jacobian', J.dtype, ('n_obs','n_state'), zlib=True, shuffle=True, complevel=5, chunksizes=(900,nstate))
                v[:] = J
                setattr(fid, 'original_file', self.jacobi_rda)

    def convert_obs(self):
        # first read the rda file
        with Timer('Read observations from RDA in '):
            rda_data = rdata.read_rda(self.obs_rda)
        obs_data = rda_data['obs_catalog']
        obs_by_dataset = defaultdict(list)
        comp_dict = {'zlib': True, 'complevel': 5, 'shuffle': True}
        with Timer("Wrote the observations to netcdf in "):
            with Dataset(self.obs_nc, 'w') as fid:
                fid.createDimension('nobs', len(obs_data))
                # write the coordinates first
                v = fid.createVariable('latitude', np.float32, ('nobs',), **comp_dict)
                v[:] = obs_data['LAT']
                v.units = 'degrees north'
                v = fid.createVariable('longitude', np.float32, ('nobs',), **comp_dict)
                v[:] = obs_data['LON']
                v.units = 'degrees east'
                v = fid.createVariable('time', np.int32, ('nobs',), **comp_dict)
                v[:] = obs_data['TIME']
                v.units = 'Minutes since 2000-01-01 00:00z'

                # convert the dataframe into array
                obs_data = obs_data.iloc[:].values

                # scale the errors to be in same units (hack to fix Andrew's file)
                data_errors = obs_data[:,5]
                for data_type, err_scale in self.error_scaling.items():
                    data_idx = obs_data[:,0] == data_type
                    data_errors[data_idx] = err_scale * data_errors[data_idx]
                obs_data[:,5] = data_errors

                # store the errors
                v = fid.createVariable('mip_mdm', np.float32, ('nobs',), **comp_dict)
                v[:] = obs_data[:,5]
                v.units = 'micromol/mol'
                v.description = 'Error on obs used in the OCO2 MIP'

                # now store the obs type
                possible_data_types, type_counts = np.unique(obs_data[:,0], return_counts=True)
                for data_type, count in zip(possible_data_types, type_counts):
                    with Timer('Storing %s data in '%data_type):
                        dim_name = 'n_%s'%data_type
                        fid.createDimension(dim_name, count)

                        data_idx = np.where(obs_data[:,0] == data_type)[0]
                        var_name = '%s_idx'%data_type
                        v = fid.createVariable(var_name, np.int32, (dim_name,), **comp_dict)
                        v[:] = data_idx
                        v.comment = 'Indices of the original array for the data type'

                        var_name = '%s_id'%(data_type)
                        if data_type == 'IS':
                            data_dtype = str
                            write_arr = obs_data[:,1][data_idx]
                            v = fid.createVariable(var_name, data_dtype, (dim_name,))
                            # a good place to segregate obs by dataset
                            for obspack_id, i in zip(write_arr, data_idx):
                                dataset = obspack_id.split('~')[1]
                                obs_by_dataset[dataset].append(i)
                        else:
                            data_dtype = np.int64
                            write_arr = np.array([int(s) for s in obs_data[:,1][data_idx]], np.int64)
                            v = fid.createVariable(var_name, data_dtype, (dim_name,), **comp_dict)
                        v[:] = write_arr

        with open(self.obs_by_dataset, 'wb') as fid:
            pickle.dump(obs_by_dataset, fid, pickle.HIGHEST_PROTOCOL)
