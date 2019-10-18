#include <LeptonInjector/DataWriter.h>

DataWriter::DataWriter(){

}

DataWriter::~DataWriter(){
    // deconstructor
}

// open the file, establish the datastructure, and create the datatypes that will be written
DataWriter::Open( std::string filename ){
    // open a new file. If one already exists by the same name - overwrite it (H5F_ACC_TRUNC)
    // leave the defaults for property lists
    // the last two can be adjusted to allow for parallel file access as I understand 
    fileHandle = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

    // Build the event properties table
    size_t dataSize = 0;
    herr_t status; // hdf5 error type. 

    dataSize += 8; // double - totalEnergy 
    dataSize += 8; // double - zenith 
    dataSize += 8; // double - azimuth 
    dataSize += 8; // double - finalStateX 
    dataSize += 8; // double - finalStateY 
    dataSize += 8 // int32_t - finalType1
    dataSize += 8 // int32_t - finalType2
    dataSize += 8 // int32_t - initialType
    dataSize += 8 // double - impactParam OR radius
    dataSize += 8 // double - totalColumnDepth or z
    
    basicPropertiesTable = H5T_CREATE(H5T_COMPOUND, dataSize);
    status = H5Tinsert(eventPropertiesTable, "totalEnergy", HOFFSET(BasicEventProperties, totalEnergy), H5T_NATIVE_DOUBLE);
    status = H5Tinsert(eventPropertiesTable, "zenith", HOFFSET(BasicEventProperties, zenith), H5T_NATIVE_DOUBLE);
    status = H5Tinsert(eventPropertiesTable, "azimuth", HOFFSET(BasicEventProperties, azimuth), H5T_NATIVE_DOUBLE);
    status = H5Tinsert(eventPropertiesTable, "finalStateX", HOFFSET(BasicEventProperties, finalStateX), H5T_NATIVE_DOUBLE);
    status = H5Tinsert(eventPropertiesTable, "finalStateY", HOFFSET(BasicEventProperties, finalStateY), H5T_NATIVE_DOUBLE);
    status = H5Tinsert(eventPropertiesTable, "finalType1", HOFFSET(BasicEventProperties, finalType1), H5T_NATIVE_LONG);
    status = H5Tinsert(eventPropertiesTable, "finalType2", HOFFSET(BasicEventProperties, finalType2), H5T_NATIVE_LONG);
    status = H5Tinsert(eventPropertiesTable, "initialType", HOFFSET(BasicEventProperties, initialType), H5T_NATIVE_LONG);

    // we want tables for volume and ranged, so let's copy that basic one and make the (slightly) different ones below
    rangedPropertiesTable = H5Tcopy( basicPropertiesTable );
    status = H5Tinsert(rangedPropertiesTable, "impactParameter", HOFFSET(RangedEventProperties, impactParameter), H5T_NATIVE_DOUBLE);
    status = H5Tinsert(rangedPropertiesTable, "totalColumnDepth", HOFFSET(RangedEventProperties, totalColumnDepth), H5T_NATIVE_DOUBLE);

    volumePropertiesTable = H5Tcopy( basicPropertiesTable );
    status = H5Tinsert(volumePropertiesTable, "radius", HOFFSET(RangedEventProperties, radius), H5T_NATIVE_DOUBLE);
    status = H5Tinsert(volumePropertiesTable, "z", HOFFSET(RangedEventProperties, z), H5T_NATIVE_DOUBLE);

}