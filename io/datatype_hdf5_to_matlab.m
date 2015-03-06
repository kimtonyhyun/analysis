function datatype_matlab = datatype_hdf5_to_matlab(datatype_hdf5)

datatype_matlab = '';
switch (datatype_hdf5)
    case 'H5T_STD_U16LE'
        datatype_matlab = 'uint16';
end