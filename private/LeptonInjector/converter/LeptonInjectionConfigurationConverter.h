#ifndef LEPTONINJECTIONCONFIGURATIONCONVERTER_H
#define LEPTONINJECTIONCONFIGURATIONCONVERTER_H

#include <tableio/I3Converter.h>
#include <LeptonInjector/LeptonInjector.h>

namespace LeptonInjector{
	
class EventPropertiesConverter : public I3ConverterImplementation<BasicEventProperties>{
public:
	I3TableRowDescriptionPtr CreateDescription(const BasicEventProperties& p);
	size_t FillRows(const BasicEventProperties& p, I3TableRowPtr rows);
};

} //namespace LeptonInjector

#endif //LEPTONINJECTIONCONFIGURATIONCONVERTER_H
