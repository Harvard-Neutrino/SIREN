#ifndef LI_H5WRITE
#define LI_H5WRITE

#include <hdf5.h>

#include <string>
#include <iostream>
#include <EventProps.h>
#include <array>

// See https://portal.hdfgroup.org/display/HDF5/HDF5
// for hdf5 documentation! 

namespace LeptonInjector{

    class DataWriter{
        public:
            DataWriter();
            ~DataWriter();

            void OpenFile(std::string filename);

            // close the active group and its sub-datasets
            // add a new group for the next injector 
            void AddInjector(std::string injector_name , bool ranged);
            void WriteEvent( BasicEventProperties& props, h5Particle& part1, h5Particle& part2, h5Particle& part3 );

        private:
            bool opened = false; 
            uint32_t event_count;

            // utility function for constructing the datatypes 
            void makeTables();

            // a 2D dataset for particles
            // [event][particle]
            hid_t initials;
            hid_t final_1;
            hid_t final_2;

            hid_t properties;

            hid_t group_handle;

            // handle for the file itself 
            hid_t fileHandle;

            // handle for the particle datatype 
            hid_t particleTable;
            hid_t rangedPropertiesTable;
            hid_t volumePropertiesTable;

            bool write_ranged;


    };

}// end namespace LeptonInjector

#endif 
