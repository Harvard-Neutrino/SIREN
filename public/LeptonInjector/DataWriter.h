#ifndef LI_H5WRITE
#define LI_H5WRITE

#include "hdf5/serial/hdf5.h"
#include <string>
#include <LeptonInjector/EventProps.h>

// See https://portal.hdfgroup.org/display/HDF5/HDF5
// for hdf5 documentation! 

namespace LeptonInjector{

    class DataWriter{
        public:
            DataWriter();
            ~DataWriter();

            hid_t OpenFile(std::string filename);
            void AddInjector(std::string injector_name );
            
            void WriteEventProperties(BasicEventProperties& props);

        private:
            void makeTables();


            hid_t fileHandle;
            
            hid_t particleTable;
            hid_t rangedPropertiesTable;
            hid_t volumePropertiesTable;

    };

}// end namespace LeptonInjector

#endif 