#ifndef LI_CONTROLLER
#define LI_CONTROLLER

#include <LeptonInjector.h>

namespace LeptonInjector{

class Controller{
    private:
        void Generate();
        LeptonInjectorBase& ActiveGenerator;
    public:
        Controller();

        void Execute(); 
        std::deque<LeptonInjectorBase*> generators;

};

} // end namespace LeptonInjector

#endif
