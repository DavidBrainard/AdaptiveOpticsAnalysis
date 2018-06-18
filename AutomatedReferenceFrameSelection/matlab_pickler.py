
import numpy
import os,pickle


def pickle_mat(pickle_path, this_dict):

    wut=numpy.reshape(this_dict['desinusoid_matrix'], (this_dict['n_columns_raw_sequence'], this_dict['n_columns_desinusoided'] ) )
    print len(wut)
    #pickle_file = open(pickle_path, 'wb')
    
    #pickle.dump(this_dict,pickle_file)

    #pickle_file.close()
