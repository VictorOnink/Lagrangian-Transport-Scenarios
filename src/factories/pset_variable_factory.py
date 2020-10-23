import numpy as np
import settings as settings


class PsetVariableFactory:
    @staticmethod
    def initialize_variable_dict_from_varlist(var_list: list, start_files: dict) -> dict:
        beach = True if "beach" in var_list else False
        age = True if "age" in var_list else False
        prox = True if "prox" in var_list else False
        size = True if "size" in var_list else False
        return PsetVariableFactory.initialize_variable_dict(start_files=start_files,
                                                            beach=beach, age=age, prox=prox, size=size)


    @staticmethod
    def initialize_variable_dict(start_files: dict,
                                 beach: bool = True, age: bool = True,
                                 prox: bool = False, size: bool = False) -> dict:
        var_dict = {}
        for variable in ['lon', 'lat', 'weight']:
            var_dict[variable] = np.load(start_files[variable])
        # if input_scenario == 'Jambeck':
        #     var_dict = _create_Jambeck(input_dir)
        # elif input_scenario == 'Lebreton':
        #     var_dict = _create_Lebreton(input_dir)
        # else:
        #     raise ValueError("invalid input scenario")

        if beach:
            var_dict['beach'] = np.zeros(len(var_dict['lon']), dtype=np.int32)
        if age:
            var_dict['age'] = np.zeros(len(var_dict['lon']), dtype=np.int32)
        if prox:
            var_dict['prox'] = np.zeros(len(var_dict['lon']), dtype=np.int32)
        if size:
            var_dict['size'] = np.ones(len(var_dict['lon']), dtype=np.float32)
        return var_dict


# def _create_Jambeck(input_dir: str) -> dict:
#     var_dict = {}
#     var_dict['lon'] = np.load(
#         input_dir + 'Jambeck2010/Jam' + str(2010) + 'Lons' + str(settings.RUN) + '.npy')
#     var_dict['lat'] = np.load(
#         input_dir + 'Jambeck2010/Jam' + str(2010) + 'Lats' + str(settings.RUN) + '.npy')
#     var_dict['weights'] = np.load(
#         input_dir + 'Jambeck2010/Jam' + str(2010) + 'Weight' + str(settings.RUN) + '.npy')
#     return var_dict
#
#
# def _create_Lebreton(input_dir: str) -> dict:
#     var_dict = {}
#     var_dict['lon'] = np.load(
#         input_dir + 'Lebreton2010/Leb' + str(2010) + 'Lons' + str(settings.RUN) + '.npy')
#     var_dict['lat'] = np.load(
#         input_dir + 'Lebreton2010/Leb' + str(2010) + 'Lats' + str(settings.RUN) + '.npy')
#     var_dict['weights'] = np.load(
#         input_dir + 'Lebreton2010/Leb' + str(2010) + 'Weight' + str(settings.RUN) + '.npy')
#     return var_dict



