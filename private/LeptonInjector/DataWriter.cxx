#include "DataWriter.h"

#include<hdf5/serial/hdf5.h>

namespace LeptonInjector{

DataWriter::DataWriter(){
}

DataWriter::~DataWriter(){
    // deconstructor
    // should kill all the hdf5 stuff
}

// open the file, establish the datastructure, and create the datatypes that will be written
hid_t DataWriter::OpenFile( std::string filename ){
    // open a new file. If one already exists by the same name - overwrite it (H5F_ACC_TRUNC)
    // leave the defaults for property lists
    // the last two can be adjusted to allow for parallel file access as I understand 
    fileHandle = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
    makeTables();

    event_count = 0;

}



void DataWriter::AddInjector( std::string injector_name , bool ranged){
    if(events!=NULL){
        H5Dclose(events);
    }
    herr_t status = H5Gclose( group_handle );

    group_handle = H5Gcreate( fileHandle, injector_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // prepare the event properties writing dataspace! 
    const hsize_t ndims = 3;
    hsize_t dims[ndims]={0,0,0};
    hsize_t max_dims[ndims] = {H5S_UNLIMITED,H5S_UNLIMITED,H5S_UNLIMITED};
    hsize_t file_space = H5Screate_simple( ndims, dims, max_dims);

    hid_t plist = H5Pcreate( H5P_DATASET_CREATE);
    H5Pset_layout( plist, H5D_CHUNKED);
    hsize_t chunk_dims[ndims] = {32768,32768,32768};
    H5Pset_chunk(plist, ndims, chunk_dims);

    const char* name = "events";
    events = H5Dcreate(group_handle,name, particleTable, file_space, H5P_DEFAULT, plist, H5P_DEFAULT );
    

    delete(&ndims);
    delete(&dims);
    delete(&max_dims);

    const hsize_t ndims =1; 
    hsize_t dims[ndims] = {0};
    hsize_t max_dims[ndims] = {H5S_UNLIMITED};
    file_space = H5Screate_simple( ndims, dims, max_dims);

    plist = H5Pcreate( H5P_DATASET_CREATE);
    H5Pset_layout( plist, H5D_CHUNKED);
    chunk_dims[ndims] = {32768};
    H5Pset_chunk(plist, ndims, chunk_dims);

    const char* props = "properties";
    if (ranged){
        events = H5Dcreate(group_handle, props , rangedPropertiesTable, file_space, H5P_DEFAULT, plist, H5P_DEFAULT );
    }else{
        events = H5Dcreate(group_handle, props , volumePropertiesTable, file_space, H5P_DEFAULT, plist, H5P_DEFAULT );
    }
}

void DataWriter::WriteEvent( BasicEventProperties& props, h5Particle& part1, h5Particle& part2, h5Particle& part3 ){
    hid_t memspace, file_space;
    const hsize_t n_dims = 3;
    hsize_t dims[n_dims] = {1,1,1};
    memspace = H5Screate_simple(n_dims, dims, NULL);

    //Extend dataset
    dims[0] = event_count+1;
    dims[1] = event_count+1;
    dims[2] = event_count+1;    
    H5Dset_extent(event_count, dims);

    //Write waveforms
    file_space = H5Dget_space(events);
    hsize_t start[3] = {event_count,event_count,event_count};
    hsize_t count[3] = {1,1,1};
    H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, count, NULL);

    std::array<h5Particle,3>  temp_data= {part1, part2, part3};
    H5Dwrite(events, particleTable, memspace, file_space, H5P_DEFAULT, &temp_data);
    delete(&temp_data);

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
    status = H5Tinsert(basicPropertiesTable, "totalEnergy", HOFFSET(BasicEventProperties, totalEnergy), H5T_NATIVE_DOUBLE);
    status = H5Tinsert(basicPropertiesTable, "zenith", HOFFSET(BasicEventProperties, zenith), H5T_NATIVE_DOUBLE);
    status = H5Tinsert(basicPropertiesTable, "azimuth", HOFFSET(BasicEventProperties, azimuth), H5T_NATIVE_DOUBLE);
    status = H5Tinsert(basicPropertiesTable, "finalStateX", HOFFSET(BasicEventProperties, finalStateX), H5T_NATIVE_DOUBLE);
    status = H5Tinsert(basicPropertiesTable, "finalStateY", HOFFSET(BasicEventProperties, finalStateY), H5T_NATIVE_DOUBLE);
    status = H5Tinsert(basicPropertiesTable, "finalType1", HOFFSET(BasicEventProperties, finalType1), H5T_NATIVE_LONG);
    status = H5Tinsert(basicPropertiesTable, "finalType2", HOFFSET(BasicEventProperties, finalType2), H5T_NATIVE_LONG);
    status = H5Tinsert(basicPropertiesTable, "initialType", HOFFSET(BasicEventProperties, initialType), H5T_NATIVE_LONG);

    // we want tables for volume and ranged, so let's copy that basic one and make the (slightly) different ones below
    rangedPropertiesTable = H5Tcopy( basicPropertiesTable );
    status = H5Tinsert(rangedPropertiesTable, "impactParameter", HOFFSET(RangedEventProperties, impactParameter), H5T_NATIVE_DOUBLE);
    status = H5Tinsert(rangedPropertiesTable, "totalColumnDepth", HOFFSET(RangedEventProperties, totalColumnDepth), H5T_NATIVE_DOUBLE);

    volumePropertiesTable = H5Tcopy( basicPropertiesTable );
    status = H5Tinsert(volumePropertiesTable, "radius", HOFFSET(VolumeEventProperties, radius), H5T_NATIVE_DOUBLE);
    status = H5Tinsert(volumePropertiesTable, "z", HOFFSET(VolumeEventProperties, z), H5T_NATIVE_DOUBLE);

    hsize_t point3d_dim[1] = {3};
    hid_t   point3d     = H5Tarray_create( H5T_NATIVE_DOUBLE, 1, point3d_dim);
    hsize_t dir3d_dim[1] = {2};
    hid_t   direction   = H5Tarray_create( H5T_NATIVE_DOUBLE, 1, dir3d_dim);
    

    dataSize = 8 + 8 + 3*8 +2*8 +8 ; // native bool, int32, position, direction, energy = 64 bytes
    //  hdf5 has asinine datasizes. Everything is 8 bytes!!! 
    particleTable = H5Tcreate( H5T_COMPOUND, dataSize);
    status = H5Tinsert(particleTable, "initial", HOFFSET(h5Particle, initial), H5T_NATIVE_HBOOL);
    status = H5Tinsert(particleTable, "ParticleType", HOFFSET(h5Particle, ptype), H5T_NATIVE_INT32);
    status = H5Tinsert(particleTable, "Position", HOFFSET(h5Particle, pos), point3d);
    status = H5Tinsert(particleTable, "Direction", HOFFSET(h5Particle, dir), direction);
    status = H5Tinsert(particleTable, "Energy", HOFFSET(h5Particle, energy), H5T_NATIVE_DOUBLE);

    

}

} // end namespace LeptonInjector