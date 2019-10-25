#include <LeptonInjector/DataWriter.h>

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

}

void DataWriter::AddInjector( std::string injector_name ){
    hid_t injector_handle = H5Gcreate( fileHandle, injector_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

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