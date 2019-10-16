#include <LeptonInjector/LeptonInjector.h>

#include <icetray/I3Logging.h>
#include <iostream>
#include <string>

#include <fstream>
#include <tableio/I3TableRowDescription.h>

#include <hdf5.h> //used for writing hdf5 files

#include <boost/detail/endian.hpp>

// implements icetray module for writing LeptonInjector configuration to an hdf5 file
// Benjamin Smithers <bsmithers@mail.icecube.wisc.edu>

// TODO:
//  + verify that the blobs are being written correctly
//  + combine the volume and ranged injectors into one: maybe use a wrapper function? 
//  + write these hdf5 files so that you don't have to index[0] to get the /only/ thing written... (optional)


namespace LeptonInjector{
    // need a little interface to trim the fat when writing the hdf5 file
    // otherwise it throws in other stuff and messes up... everything.    
    struct simple_holder{
        uint32_t events;
        double energyMinimum;
        double energyMaximum;
        double powerlawIndex;
        double azimuthMinimum;
        double azimuthMaximum;
        double zenithMinimum;
        double zenithMaximum;
        int32_t finalType1;
        int32_t finalType2;
//        std::string crossSectionBlob;
 //       std::string totalCrossSectionBlob;
    };

//template<class InputIterator> string (InputIterator begin, InputIterator end);
    //helper function to fill in the injectors 
    void construct_holder( const BasicInjectionConfiguration& config, simple_holder* target){
        target->events       = config.events;
        target->energyMinimum    = config.energyMinimum;
        target->energyMaximum    = config.energyMaximum;
        target->powerlawIndex= config.powerlawIndex;
        target->azimuthMinimum   = config.azimuthMinimum;
        target->azimuthMaximum   = config.azimuthMaximum;
        target->zenithMinimum    = config.zenithMinimum;
        target->zenithMaximum    = config.zenithMaximum;
        target->finalType1   = config.finalType1;
        target->finalType2   = config.finalType2;

        //target->crossSectionBlob = config.crossSectionBlob;
        //target->crossSectionBlob.push_back('\0');

        //target->totalCrossSectionBlob = config.totalCrossSectionBlob;
        //target->crossSectionBlob.push_back('\0');

        //log_trace("filling the hdf5 interface object with the cross sectinos");
        //make sure the strings are big enough for what's about to hit them
   //     target->crossSectionBlob.reserve( config.crossSectionBlob.size() + 1); //add one byte for the null terminator
     //   target->totalCrossSectionBlob.reserve(config.totalCrossSectionBlob.size() + 1);
        
        // because these iterate over the vector when filling the string. Yuck! 
    //    target->crossSectionBlob = std::string( config.crossSectionBlob.begin(), config.crossSectionBlob.end() );
    //    target->crossSectionBlob.push_back('\0');
    //    target->totalCrossSectionBlob = std::string( config.totalCrossSectionBlob.begin(), config.totalCrossSectionBlob.end() );
    //    target->totalCrossSectionBlob.push_back('\0');
       //   target->crossSectionBlob = "test1";
     //     target->totalCrossSectionBlob = "test2323";
    }

    // derived structs for ranged/volume
    struct simple_ranged_holder : public simple_holder{
        double injectionRadius;
        double endcapLength;
    };
    struct simple_volume_holder : public simple_holder{
        double cylinderRadius;
        double cylinderHeight;
    };

    
	// writes a table onto the file for a volume injector configuration
	void writeConfiguration(hid_t file, const VolumeInjectionConfiguration& config, unsigned int configurator){
        simple_volume_holder temp;
        construct_holder( config, &temp);
        temp.cylinderRadius = config.cylinderRadius;
        temp.cylinderHeight = config.cylinderHeight;

        uint64_t dataSize=0;
		dataSize+=4; //events
		dataSize+=8; //energyMinimum
		dataSize+=8; //energyMaximum
		dataSize+=8; //powerlawIndex
		dataSize+=8; //azimuthMinimum
		dataSize+=8; //azimuthMaximum
		dataSize+=8; //zenithMinimum
		dataSize+=8; //zenithMaximum
		dataSize+=4; //finalType1
		dataSize+=4; //finalType2
		//dataSize+=temp.crossSectionBlob.length(); 
        //dataSize+=temp.totalCrossSectionBlob.length(); //crossSection
		dataSize+=8; // cylinder radius
		dataSize+=8; // cylinder height   

        // create group name group 
        log_trace("creating the group");
        std::string table_name = "/Volume_Injector" + std::to_string(configurator);

        hid_t vartype = H5Tcopy( H5T_C_S1 );
        herr_t status;
        status = H5Tset_size( vartype, H5T_VARIABLE);

        hsize_t data_table = H5Tcreate( H5T_COMPOUND, dataSize );
        status = H5Tinsert( data_table, "events",        HOFFSET( simple_volume_holder, events),          H5T_NATIVE_UINT32);
        status = H5Tinsert( data_table, "energyMinimum", HOFFSET( simple_volume_holder, energyMinimum),   H5T_NATIVE_DOUBLE);
        status = H5Tinsert( data_table, "energyMaximum", HOFFSET( simple_volume_holder, energyMaximum),   H5T_NATIVE_DOUBLE);
        status = H5Tinsert( data_table, "powerlawIndex", HOFFSET( simple_volume_holder, powerlawIndex),   H5T_NATIVE_DOUBLE);
        status = H5Tinsert( data_table, "azimuthMinimum",HOFFSET( simple_volume_holder, azimuthMinimum),  H5T_NATIVE_DOUBLE);
        status = H5Tinsert( data_table, "azimuthMaximum",HOFFSET( simple_volume_holder, azimuthMaximum),  H5T_NATIVE_DOUBLE);
        status = H5Tinsert( data_table, "zenithMinimum", HOFFSET( simple_volume_holder, zenithMinimum),   H5T_NATIVE_DOUBLE);
        status = H5Tinsert( data_table, "zenithMaximum", HOFFSET( simple_volume_holder, zenithMaximum),   H5T_NATIVE_DOUBLE);
        status = H5Tinsert( data_table, "finalType1",    HOFFSET( simple_volume_holder, finalType1),      H5T_NATIVE_INT32); //final type is an enum relating names to ints, so here we have it as an int
        status = H5Tinsert( data_table, "finalType2",    HOFFSET( simple_volume_holder, finalType2),      H5T_NATIVE_INT32);
        //status = H5Tinsert( data_table, "crossSectionBlob",  HOFFSET( simple_volume_holder, crossSectionBlob), vartype );
        //status = H5Tinsert( data_table, "totalCrossSectionBlob",HOFFSET(simple_volume_holder, totalCrossSectionBlob), vartype );
        status = H5Tinsert( data_table,"cylinderRadius", HOFFSET( simple_volume_holder, cylinderRadius), H5T_NATIVE_DOUBLE);
        status = H5Tinsert( data_table,"cylinderHeight", HOFFSET( simple_volume_holder, cylinderHeight), H5T_NATIVE_DOUBLE);

        
        hid_t plist     = H5Pcreate(H5P_DATASET_CREATE);
        hsize_t dims[1] = {1};

        const long long unsigned int max_int = 1;
        hid_t dataspace = H5Screate_simple(1 , dims, &max_int);
        hid_t dataset   = H5Dcreate( file, table_name.c_str(), data_table, dataspace, H5P_DEFAULT, plist, H5P_DEFAULT);
        hid_t file_space= H5Dget_space(dataset);

        H5Dwrite( dataset, data_table, dataspace, file_space, H5P_DEFAULT, &temp);

        H5Sclose( file_space );
        H5Sclose( dataspace  );
        

	}
    
    // writes a datatable to the file for a ranged injection configuration
	void writeConfiguration(hid_t file, const RangedInjectionConfiguration& config, unsigned int configurator){
		//compute data size
        // 1  uint8_t version
        // 2  uint32_t events
        // 3  double energyMinimum
        // 4  double energyMaximum
        // 5  double powerlawIndex
        // 6  double azimuthMinimum
        // 7  double azimuthMaximum
        // 8  double zenithMinimum
        // 9  double zenithMaximum
        // 10 int32_t finalType1.GetType()
        // 11 int32_t finalType2.GetType()
        // 12 vector<char> crossSectionBlob
        // 13 vector<char> totalCrossSectionBlob
        // 14 double injectionRadius
        // 15 double endcapLength
            
        simple_ranged_holder temp;
        construct_holder( config, &temp);
        temp.injectionRadius = config.injectionRadius;
        temp.endcapLength = config.endcapLength;


 		uint64_t dataSize=0;
		dataSize+=4+4; //events
		dataSize+=8; //energyMinimum
		dataSize+=8; //energyMaximum
		dataSize+=8; //powerlawIndex
		dataSize+=8; //azimuthMinimum
		dataSize+=8; //azimuthMaximum
		dataSize+=8; //zenithMinimum
		dataSize+=8; //zenithMaximum
		dataSize+=4+4; //finalType1
		dataSize+=4+4; //finalType2
//		dataSize+=temp.crossSectionBlob.size(); //crossSection
//		dataSize+=temp.totalCrossSectionBlob.size(); //crossSection
		dataSize+=8; //injection radius
		dataSize+=8; // endcap length     

        // create group name group 
        log_trace("creating the group");
        std::string table_name = "/Ranged_Injector" + std::to_string(configurator);
        
        hid_t tottype = H5Tcopy( H5T_C_S1 );
        hid_t diftype = H5Tcopy( H5T_C_S1 );
        herr_t status = 0;
        status = H5Tset_size( tottype , H5T_VARIABLE); //temp.totalCrossSectionBlob.length() );
        if(status<0){throw("Bad status!!");}
        status = H5Tset_size( diftype , H5T_VARIABLE); //temp.crossSectionBlob.length() );
        if(status<0){throw("Bad status!!");}

        hid_t true_double = H5Tcopy( H5T_NATIVE_DOUBLE);
        status = H5Tset_size( true_double, 8 );
        if(status<0){throw("Bad Status!!!");}

        hsize_t data_table = H5Tcreate( H5T_COMPOUND, dataSize );
        uint64_t current_offset = 0;
        log_trace("creating the table type");
        
        /*
        status = H5Tinsert( data_table, "events",        current_offset,   H5T_NATIVE_UINT32); current_offset+=4; 
        status = H5Tinsert( data_table, "energyMinimum", current_offset,   true_double); current_offset+=8;
        status = H5Tinsert( data_table, "energyMaximum", current_offset,   true_double); current_offset+=8;
        status = H5Tinsert( data_table, "powerlawIndex", current_offset,   true_double); current_offset+=8;
        status = H5Tinsert( data_table, "azimuthMinimum",current_offset,   true_double); current_offset+=8;
        status = H5Tinsert( data_table, "azimuthMaximum",current_offset,   true_double); current_offset+=8;
        status = H5Tinsert( data_table, "zenithMinimum", current_offset,   true_double); current_offset+=8;
        status = H5Tinsert( data_table, "zenithMaximum", current_offset,   true_double); current_offset+=8;
        status = H5Tinsert( data_table, "finalType1",    current_offset,   H5T_NATIVE_INT32);  current_offset+=4; //final type is an enum relating names to ints, so here we have it as an int
        status = H5Tinsert( data_table, "finalType2",    current_offset,   H5T_NATIVE_INT32);  current_offset+=4;
        status = H5Tinsert( data_table, "crossSectionBlob",  current_offset, diftype );        current_offset+=temp.crossSectionBlob.length();
        status = H5Tinsert( data_table, "totalCrossSectoinBlob",current_offset, tottype );     current_offset+=temp.totalCrossSectionBlob.length();
        status = H5Tinsert( data_table,"injectionRadius",current_offset, true_double);   current_offset+=8;
        status = H5Tinsert( data_table,"endcapLength",   current_offset,    true_double);current_offset+=8;
        */
        
        status = H5Tinsert( data_table, "events",        HOFFSET( simple_ranged_holder, events),          H5T_NATIVE_UINT32);
        status = H5Tinsert( data_table, "energyMinimum", HOFFSET( simple_ranged_holder, energyMinimum),   H5T_NATIVE_DOUBLE);
        status = H5Tinsert( data_table, "energyMaximum", HOFFSET( simple_ranged_holder, energyMaximum),   H5T_NATIVE_DOUBLE);
        status = H5Tinsert( data_table, "powerlawIndex", HOFFSET( simple_ranged_holder, powerlawIndex),   H5T_NATIVE_DOUBLE);
        status = H5Tinsert( data_table, "azimuthMinimum",HOFFSET( simple_ranged_holder, azimuthMinimum),  H5T_NATIVE_DOUBLE);
        status = H5Tinsert( data_table, "azimuthMaximum",HOFFSET( simple_ranged_holder, azimuthMaximum),  H5T_NATIVE_DOUBLE);
        status = H5Tinsert( data_table, "zenithMinimum", HOFFSET( simple_ranged_holder, zenithMinimum),   H5T_NATIVE_DOUBLE);
        status = H5Tinsert( data_table, "zenithMaximum", HOFFSET( simple_ranged_holder, zenithMaximum),   H5T_NATIVE_DOUBLE);
        status = H5Tinsert( data_table, "finalType1",    HOFFSET( simple_ranged_holder, finalType1),      H5T_NATIVE_INT32); //final type is an enum relating names to ints, so here we have it as an int
        status = H5Tinsert( data_table, "finalType2",    HOFFSET( simple_ranged_holder, finalType2),      H5T_NATIVE_INT32);
      //  status = H5Tinsert( data_table, "crossSectionBlob",  HOFFSET( simple_ranged_holder, crossSectionBlob), diftype );
      //  status = H5Tinsert( data_table, "totalCrossSectoinBlob",HOFFSET(simple_ranged_holder, totalCrossSectionBlob), tottype );
        status = H5Tinsert( data_table,"injectionRadius",HOFFSET( simple_ranged_holder, injectionRadius), H5T_NATIVE_DOUBLE);
        status = H5Tinsert( data_table,"endcapLength",   HOFFSET( simple_ranged_holder, endcapLength),    H5T_NATIVE_DOUBLE);
        

        hid_t plist     = H5Pcreate(H5P_DATASET_CREATE);
        hsize_t dims[3] = {1,1,1};
        const long long unsigned int max_int = 1;

        hid_t dataspace = H5Screate_simple(3 , dims, &max_int);
        hid_t dataset   = H5Dcreate( file, table_name.c_str(), data_table, dataspace, H5P_DEFAULT, plist, H5P_DEFAULT);
//        hid_t file_space= H5Dget_space(dataset);
        hid_t dataset2  = H5Dcreate( file, table_name.c_str(), diftype   , dataspace, H5P_DEFAULT, plist, H5P_DEFAULT);
        hid_t dataset3  = H5Dcreate( file, table_name.c_str(), diftype   , dataspace, H5P_DEFAULT, plist, H5P_DEFAULT);

        log_trace("writing injector to hdf5 file");
//        status = H5Dwrite( dataset, data_table, dataspace, file_space, H5P_DEFAULT, &temp);
        status = H5Dwrite( dataset, data_table, H5S_ALL, H5S_ALL, H5P_DEFAULT, &temp);
        status = H5Dwrite( dataset2, diftype,   H5S_ALL, H5S_ALL, H5P_DEFAULT, &(config.crossSectionBlob));
        status = H5Dwrite( dataset3, diftype,   H5S_ALL, H5S_ALL, H5P_DEFAULT, &(config.totalCrossSectionBlob));

        if (status<0){throw("bad status");}

//        status = H5Sclose( file_space );
 //       status = H5Sclose( dataspace  );

       	}

	
	class InjectionConfigH5Write : public I3Module{
	public:
		InjectionConfigH5Write(const I3Context& ctx):
		I3Module(ctx),wroteParticleTypes(false){
			Register(I3Frame::Stream('S'),&InjectionConfigH5Write::S);
			AddParameter("OutputPath","");
			AddOutBox("OutBox");
		}
        unsigned int configurator = 0;
        hid_t file;
		
        // called once as the module is instantiated. Creates a handle to the hdf5 file, and opens the hdf5 file. 
		void Configure(){
            log_trace("configuring injection config h5 write");
			std::string outputPath;
			GetParameter("OutputPath",outputPath);		    
            file = H5Fcreate(outputPath.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

            log_trace("we opened up an hdf5 file");
            //output.open(outputPath.c_str());
		}
		
        // called for every S frame. (once for each injector made)
		void S(boost::shared_ptr<I3Frame> frame){
            // in a multi LI situation, there will be multiple frames. Makes things compicated... 
            log_trace("Okay, so now we're processing an S frame");
			if(frame->Has("LeptonInjectorProperties")){
				if(!wroteParticleTypes){
					MAKE_ENUM_VECTOR(type,I3Particle,I3Particle::ParticleType,I3PARTICLE_H_I3Particle_ParticleType);
					wroteParticleTypes=true;
				}
				boost::shared_ptr<const I3FrameObject> config
				  =frame->Get<boost::shared_ptr<const I3FrameObject> >("LeptonInjectorProperties");
				boost::shared_ptr<const RangedInjectionConfiguration> rconfig;
				boost::shared_ptr<const VolumeInjectionConfiguration> vconfig;
				if((rconfig=boost::dynamic_pointer_cast<const RangedInjectionConfiguration>(config)))
					writeConfiguration(file, *rconfig, configurator);
				else if((vconfig=boost::dynamic_pointer_cast<const VolumeInjectionConfiguration>(config)))
                    writeConfiguration(file, *vconfig, configurator);
				else
					log_fatal("Unable to write LeptonInjectorProperties with unrecognized type");
			}
            configurator++; // remember to increase this so that each configuration object in the hdf5 file has a unique name
			PushFrame(frame);
		}
		
	private:
		std::ofstream output;
		bool wroteParticleTypes;
	};
	
	I3_MODULE(InjectionConfigH5Write);
	
} //namespace LeptonInjector

/*
I3_SERIALIZABLE(LeptonInjector::RangedInjectionConfiguration);
I3_SERIALIZABLE(LeptonInjector::BasicInjectionConfiguration);
I3_SERIALIZABLE(LeptonInjector::VolumeInjectionConfiguration);
I3_SERIALIZABLE(LeptonInjector::BasicEventProperties);
I3_SERIALIZABLE(LeptonInjector::RangedEventProperties);
I3_SERIALIZABLE(LeptonInjector::VolumeEventProperties);
*/

