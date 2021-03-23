#ifndef LI_DataWriter_H
#define LI_DataWriter_H

#include <array>
#include <string> // strings
#include <vector> //MAKE_ENUM_VECTOR
#include <fstream> // std::ostream 
#include <iostream> // std::cout 

#include <hdf5.h> // all the data writing
#include <boost/assign/list_inserter.hpp>
#include <boost/preprocessor/seq/transform.hpp>

#include "LeptonInjector/EventProps.h"
#include "LeptonInjector/BasicInjectionConfiguration.h"

// See https://portal.hdfgroup.org/display/HDF5/HDF5
// for hdf5 documentation! 

namespace LeptonInjector{

    class DataWriter{
        public:
            DataWriter();
            ~DataWriter();

            void OpenFile(std::string filename);
            void OpenLICFile( std::string filename);
            void SetOverwrite( bool overwrite_ ){this->overwrite = overwrite_;}


            // close the active group and its sub-datasets
            // add a new group for the next injector 
            void AddInjector(std::string injector_name , bool ranged);
            void WriteEvent( BasicEventProperties& props, h5Particle& part1, h5Particle& part2, h5Particle& part3 );
            void WriteConfig( BasicInjectionConfiguration& config, bool ranged );

        private:
            bool opened = false; 
            bool overwrite = false;
            uint8_t name_iterator = 0;
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

            std::ofstream lic_file_output;

            void writeRangedConfig( BasicInjectionConfiguration& config );
            void writeVolumeConfig( BasicInjectionConfiguration& config );
    };

}// end namespace LeptonInjector

#endif // LI_DataWriter_H

