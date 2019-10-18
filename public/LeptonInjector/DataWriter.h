#ifndef LI_H5WRITE
#define LI_H5WRITE

#include <hdf5.h>

// See https://portal.hdfgroup.org/display/HDF5/HDF5
// for hdf5 documentation! 

namespace LeptonInjector{

    class DataWriter{
        public:
            DataWriter();
            ~DataWriter();

            void OpenFile(std::string filename);

            void WriteEventProperties(&BasicEventProperties props);

        private:
            hid_t fileHandle;
            hid_t basicPropertiesTable;
            hid_t rangedPropertiesTable;
            hid_t volumePropertiesTable;

    };

}// end namespace LeptonInjector

#endif 