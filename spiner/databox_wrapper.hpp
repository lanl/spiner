#ifndef SINGE_UTIL_DATABOXWRAPPER_HPP
#define SINGE_UTIL_DATABOXWRAPPER_HPP

#include "singe/util/databox.hpp"
#include "singe/util/portability_macros.hpp"

#include <cmath>

namespace singe {
namespace util {

// ================================================================================================

// a wrapper to (a) adapt the interface and (b) handle extrapolation off the table
// -- note that we assume this is a 1D databox
template <typename DataType_t, typename LowerExtrapolation, typename UpperExtrapolation>
class DataBoxWrapper
{
private:
    DataBox<DataType_t> databox_;
    // Avoid re-calculating temperature bounds on every call
    DataType_t Tlo_;
    DataType_t Thi_;

public:
    SINGE_PORTABLE_FUNCTION DataBoxWrapper(DataBox<DataType_t> databox)
        : databox_{databox} // MUST be initialized before Tlo_ and Thi_
        , Tlo_{get_temperature(0)}
        , Thi_{get_temperature(size() - 1)}
    {}
    SINGE_PORTABLE_FUNCTION constexpr DataType_t get_logT(int const index) const
    {
        return databox_.range(0).x(index);
    }
    SINGE_PORTABLE_FUNCTION constexpr DataType_t get_temperature(int const index) const
    {
        return std::exp(get_logT(index));
    }
    SINGE_PORTABLE_FUNCTION constexpr DataType_t get_value(int const index) const
    {
        return databox_(index);
    }
    SINGE_PORTABLE_FUNCTION constexpr DataType_t get_logV(int const index) const
    {
        return std::log(get_value(index));
    }
    SINGE_PORTABLE_FUNCTION constexpr int size() const
    {
        return databox_.dim(1);
    }
    SINGE_PORTABLE_FUNCTION constexpr DataType_t operator()(DataType_t const temperature) const
    {
        if (temperature < Tlo_) {
            return LowerExtrapolation::extrapolate(temperature, *this);
        } else if (temperature > Thi_) {
            return UpperExtrapolation::extrapolate(temperature, *this);
        } else {
            return databox_.interpToReal(std::log(temperature));
        }
    }
    SINGE_PORTABLE_FUNCTION constexpr void finalize()
    {
        databox_.finalize();
    }
};

// ================================================================================================

} // end namespace util
} // end namespace singe

#endif // ifndef SINGE_UTIL_DATABOXWRAPPER_HPP
