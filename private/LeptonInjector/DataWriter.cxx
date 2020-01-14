#include "DataWriter.h"


namespace LeptonInjector{

DataWriter::DataWriter(){
}

// close all the hdf5 things that are open! 
DataWriter::~DataWriter(){
    //avoid the awkwatd situation where the deconstructot is called before you even do anything
    if(this->opened){
        herr_t status;
        status = H5Tclose( particleTable );
        status = H5Tclose( rangedPropertiesTable );
        status = H5Tclose( volumePropertiesTable );

        status = H5Gclose( group_handle ); 
        status = H5Dclose( initials );
        status = H5Dclose( final_1 );
        status = H5Dclose( final_2 );
        status = H5Dclose( properties );

        status = H5Fclose( fileHandle );
    }
    // deconstructor
    // should kill all the hdf5 stuff
}

// open the file, establish the datastructure, and create the datatypes that will be written
void DataWriter::OpenFile( std::string filename ){
    this->opened = true;
    // open a new file. If one already exists by the same name - overwrite it (H5F_ACC_TRUNC)
    // leave the defaults for property lists
    // the last two can be adjusted to allow for parallel file access as I understand, which may allow for multi-processing
    //      the different injectors, and therefore much quicker injection 
    fileHandle = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
    makeTables();

    event_count = 0;

}



void DataWriter::AddInjector( std::string injector_name , bool ranged){
    herr_t status;

    event_count = 0;

    if(name_iterator>0){
        status = H5Dclose(initials);
        status = H5Dclose(final_1);
        status = H5Dclose(final_2);
        status = H5Dclose(properties);
        status = H5Gclose( group_handle );
    }
    
    group_handle = H5Gcreate( fileHandle, (injector_name + std::to_string(name_iterator)).c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    name_iterator++;

    // prepare the event properties writing dataspace! 
    const hsize_t ndims = 1;
    hsize_t dims[ndims]={0}; // current dimensionality of the file space is zero
    hsize_t max_dims[ndims] = {H5S_UNLIMITED}; // but can be expanded infinitely 
    hsize_t file_space = H5Screate_simple( ndims, dims, max_dims);

    hid_t plist = H5Pcreate( H5P_DATASET_CREATE);
    H5Pset_layout( plist, H5D_CHUNKED);
    // 32768
    // any time you have an infinite sized dimension, you need to chunk it. 1000 is a nice round number 
    const int nice_round_number = 32768;
    hsize_t chunk_dims[ndims] = {nice_round_number};
    H5Pset_chunk(plist, ndims, chunk_dims);

    const char* init_name = "initial";
    const char* final_1_name = "final_1";
    const char* final_2_name = "final_2";
    initials  = H5Dcreate(group_handle,init_name, particleTable, file_space, H5P_DEFAULT, plist, H5P_DEFAULT );
    final_1 = H5Dcreate(group_handle,final_1_name, particleTable, file_space, H5P_DEFAULT, plist, H5P_DEFAULT); 
    final_2 = H5Dcreate(group_handle,final_2_name, particleTable, file_space, H5P_DEFAULT, plist, H5P_DEFAULT); 
    H5Sclose(file_space);


    const hsize_t ndims2 =1; 
    hsize_t dims2[ndims2] = {0};
    hsize_t max_dims2[ndims2] = {H5S_UNLIMITED};
    hsize_t file_space2 = H5Screate_simple( ndims2, dims2, max_dims2);

    hid_t plist2 = H5Pcreate( H5P_DATASET_CREATE);
    H5Pset_layout( plist2, H5D_CHUNKED);
    hsize_t chunk_dims2[ndims2] = {nice_round_number};
    H5Pset_chunk(plist2, ndims2, chunk_dims2);

    const char* props = "properties";
    if (ranged){
        write_ranged = true;
        properties = H5Dcreate(group_handle, props , rangedPropertiesTable, file_space2, H5P_DEFAULT, plist2, H5P_DEFAULT );
    }else{
        write_ranged = false;
        properties = H5Dcreate(group_handle, props , volumePropertiesTable, file_space2, H5P_DEFAULT, plist2, H5P_DEFAULT );
    }
    std::cout << "tried making second dataset" << std::endl;

    H5Sclose(file_space2);
}

void DataWriter::WriteEvent( BasicEventProperties& props, h5Particle& part1, h5Particle& part2, h5Particle& part3 ){
    // std::cout << " writing an event" << std::endl;
    hid_t memspace, file_space;
    const hsize_t n_dims = 1;
    hsize_t dims[n_dims] = {1};
    memspace = H5Screate_simple(n_dims, dims, NULL);

    //Extend dataset
    dims[0] = event_count+1;
    H5Dset_extent(initials, dims);
    H5Dset_extent(final_1, dims);
    H5Dset_extent(final_2, dims);
    H5Dset_extent(properties, dims);

    //Write waveforms
    file_space = H5Dget_space(initials);
    hsize_t start[n_dims] = {event_count};
    hsize_t count[n_dims] = {1};
    H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, count, NULL);


    // h5Particle temp_data[3] = {part1, part2, part3};
    H5Dwrite(initials, particleTable, memspace, file_space, H5P_DEFAULT, &part1);
    file_space = H5Dget_space( final_1 );
    H5Sselect_hyperslab( file_space, H5S_SELECT_SET, start, NULL, count, NULL);
    H5Dwrite(final_1, particleTable, memspace, file_space, H5P_DEFAULT, &part2);
    file_space = H5Dget_space( final_2 );
    H5Sselect_hyperslab( file_space, H5S_SELECT_SET, start, NULL, count, NULL);
    H5Dwrite(final_2, particleTable, memspace, file_space, H5P_DEFAULT, &part3);

    file_space = H5Dget_space( properties );
    H5Sselect_hyperslab( file_space, H5S_SELECT_SET, start, NULL, count, NULL);
    if (write_ranged){
        H5Dwrite(properties, rangedPropertiesTable, memspace, file_space, H5P_DEFAULT, &props);
    }else{
        H5Dwrite(properties, volumePropertiesTable, memspace, file_space, H5P_DEFAULT, &props);
    }

    // delete(&temp_data);

    // std::cout << "cleaning up" << std::endl;
    H5Sclose(file_space);
    H5Sclose(memspace);
    event_count++; 
}

void DataWriter::makeTables(){
    // Build the event properties table
    size_t dataSize = 0;
    dataSize += 8; // double - totalEnergy 
    dataSize += 8; // double - zenith 
    dataSize += 8; // double - azimuth 
    dataSize += 8; // double - finalStateX 
    dataSize += 8; // double - finalStateY 
    dataSize += 8; // int32_t - finalType1 // yeah, this should be 4 bytes. BUT hdf5 is /really/ dumb about data sizes 
    dataSize += 8; // int32_t - finalType2
    dataSize += 8; // int32_t - initialType
    dataSize += 8; // double - impactParam OR radius
    dataSize += 8; // double - totalColumnDepth or z
    
    herr_t status; // hdf5 error type. 
    hid_t basicPropertiesTable = H5Tcreate(H5T_COMPOUND, dataSize);
    size_t offset = 0;
    hid_t real_long = H5Tcopy( H5T_NATIVE_LONG );
    status = H5Tset_size( real_long, 4);
    hid_t real_bool = H5Tcopy( H5T_NATIVE_HBOOL);
    status = H5Tset_size( real_bool, 1);

    status = H5Tinsert(basicPropertiesTable, "totalEnergy", HOFFSET(BasicEventProperties, totalEnergy), H5T_NATIVE_DOUBLE);
    status = H5Tinsert(basicPropertiesTable, "zenith", HOFFSET(BasicEventProperties, zenith), H5T_NATIVE_DOUBLE); 
    status = H5Tinsert(basicPropertiesTable, "azimuth", HOFFSET(BasicEventProperties, azimuth) , H5T_NATIVE_DOUBLE); 
    status = H5Tinsert(basicPropertiesTable, "finalStateX", HOFFSET(BasicEventProperties, finalStateX) , H5T_NATIVE_DOUBLE);
    status = H5Tinsert(basicPropertiesTable, "finalStateY", HOFFSET(BasicEventProperties, finalStateY) , H5T_NATIVE_DOUBLE); 
    status = H5Tinsert(basicPropertiesTable, "finalType1", HOFFSET(BasicEventProperties, finalType1) , real_long); 
    status = H5Tinsert(basicPropertiesTable, "finalType2", HOFFSET(BasicEventProperties, finalType2) , real_long);
    status = H5Tinsert(basicPropertiesTable, "initialType", HOFFSET(BasicEventProperties, initialType) , real_long); 

    // we want tables for volume and ranged, so let's copy that basic one and make the (slightly) different ones below
    rangedPropertiesTable = H5Tcopy( basicPropertiesTable );
    status = H5Tinsert(rangedPropertiesTable, "impactParameter", HOFFSET(RangedEventProperties, impactParameter) , H5T_NATIVE_DOUBLE); 
    status = H5Tinsert(rangedPropertiesTable, "totalColumnDepth", HOFFSET(RangedEventProperties, totalColumnDepth) , H5T_NATIVE_DOUBLE); 

    volumePropertiesTable = H5Tcopy( basicPropertiesTable );
    status = H5Tinsert(volumePropertiesTable, "radius", HOFFSET(VolumeEventProperties, radius) , H5T_NATIVE_DOUBLE); 
    status = H5Tinsert(volumePropertiesTable, "z", HOFFSET(VolumeEventProperties, z) , H5T_NATIVE_DOUBLE); 

    H5Tclose( basicPropertiesTable );

    hsize_t point3d_dim[1] = {3};
    hid_t   point3d     = H5Tarray_create( H5T_NATIVE_DOUBLE, 1, point3d_dim);
    hsize_t dir3d_dim[1] = {2};
    hid_t   direction   = H5Tarray_create( H5T_NATIVE_DOUBLE, 1, dir3d_dim);
    
    dataSize = 8 + 8 + 3*8 +2*8 +8 ; // native bool, int32, position, direction, energy = 64 bytes
    //  hdf5 has asinine datasizes. Everything is 8 bytes!!! 
    offset = 0;
    particleTable = H5Tcreate( H5T_COMPOUND, dataSize);
    status = H5Tinsert(particleTable, "initial",        HOFFSET(h5Particle, initial), real_bool);
    status = H5Tinsert(particleTable, "ParticleType",   HOFFSET(h5Particle, ptype),  real_long); 
    status = H5Tinsert(particleTable, "Position",       HOFFSET(h5Particle, pos) , point3d);     
    status = H5Tinsert(particleTable, "Direction",      HOFFSET(h5Particle, dir) , direction); 
    status = H5Tinsert(particleTable, "Energy",         HOFFSET(h5Particle, energy) , H5T_NATIVE_DOUBLE); 

    /*
    status = H5Tinsert(particleTable, "initial", offset, H5T_NATIVE_HBOOL); offset +=8;
    status = H5Tinsert(particleTable, "ParticleType", offset,  H5T_NATIVE_LONG); offset += 8;
    status = H5Tinsert(particleTable, "Position", offset , point3d); offset += 24;
    status = H5Tinsert(particleTable, "Direction", offset , direction); offset += 16;
    status = H5Tinsert(particleTable, "Energy", offset , H5T_NATIVE_DOUBLE); offset += 8;
    */

}

} // end namespace LeptonInjector
