import logging
import os
from datetime import date
from UnitTesting.create_dict_string import create_dict_string

# [first_time_print] takes in a module [mod], a value dictionary [value_dict], a path [path], and a boolean [write].
# It prints to the console the properly formatted trusted_values_dict entry based on [mod] and [value_dict].
# Additionally, if [write] is [True], it appends this output to the file [path]/trusted_values_dict.py

# Called by run_test

# Uses self.module_name, self.trusted_values_dict_name, self.calculated_dict, self.path


def first_time_print(self, write=True):
    dict_name = self.trusted_values_dict_name

    output_string = r"""
# Generated on: """ + str(date.today()) + r"""
trusted_values_dict['""" + dict_name + r"""'] = """ + \
                  str(create_dict_string(self.calculated_dict))

    error = r"""
Module: """ + self.module_name + r"""
Please copy the following code between the ##### and paste it into your trusted_values_dict.py file for this module:

#####

""" + output_string + """

#####
"""
    logging.error(error)

    # If [write] is [True], write to [trusted_values_dict]
    if write:
        logging.debug(' Writing trusted_values_dict entry to trusted_values_dict.py...')
        with open(os.path.join(self.path, 'trusted_values_dict.py'), 'a') as file:
            file.write(output_string)
        logging.debug(' ...Success: entry written to trusted_values_dict.py\n')
