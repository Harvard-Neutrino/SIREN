#include <dlfcn.h>
#include <iostream>

// Do not link this executable to SIREN or Python: doing so would supply the
// symbols whose missing library dependency this test is intended to detect.
int main(int argc, char **argv) {
    if (argc != 2) {
        std::cerr << "Usage: StandaloneLoad_TEST /path/to/libSIREN\n";
        return 2;
    }
    void *library = dlopen(argv[1], RTLD_NOW | RTLD_LOCAL);
    if (!library) {
        std::cerr << dlerror() << '\n';
        return 1;
    }
    std::cout << "Loaded standalone SIREN successfully\n";
    // Keep the library loaded until process exit, as a native plugin host can.
    return 0;
}
