// serialization.cxx
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>

// Register the archives
CEREAL_REGISTER_ARCHIVE(cereal::BinaryOutputArchive);
CEREAL_REGISTER_ARCHIVE(cereal::BinaryInputArchive);

// Ensure dynamic initialization is not optimized away
CEREAL_REGISTER_DYNAMIC_INIT(siren);
CEREAL_FORCE_DYNAMIC_INIT(siren);
