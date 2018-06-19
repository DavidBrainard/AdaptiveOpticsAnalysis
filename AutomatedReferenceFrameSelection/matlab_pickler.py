
import numpy
import os,pickle


def pickle_mat(pickle_path, this_dict):
	

    this_dict['desinusoid_matrix']=numpy.reshape(this_dict['desinusoid_matrix'], (int(this_dict['n_columns_desinusoided']), int(this_dict['n_columns_raw_sequence']) ) )
    this_dict['desinusoid_matrix']=numpy.transpose(this_dict['desinusoid_matrix'])

    pickle_file = open(pickle_path, 'wb')
    
    pickle.dump(this_dict,pickle_file)

    pickle_file.close()
