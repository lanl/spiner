#ifndef _SPINER_DATABOXWRAPPER_HPP_
#define _SPINER_DATABOXWRAPPER_HPP_

#include "databox.hpp"

#include "ports-of-call//portability.hpp"

#include <cmath>

namespace Spiner {

// ================================================================================================

// The original use-case did three things:
// * Add transformations of the data.  This was originally hard-coded, but we'll make it general.
// * Add extrapolations off the edge of the table.
// * Adapt the interface.  We should instead follow the Databox interface.
//   * DataBox() = default;
//   * DataBox(T * data, Args... args) noexcept;
//   * DataBox(AllocationTarget t, Args... args) noexcept;
//   * DataBox(Args... args) noexcept;
//   * DataBox(PortableMDArray<T> A) noexcept;
//   * DataBox(PortableMDArray<T>& A) noexcept;
//   * DataBox(const DataBox& src) noexcept;
//   * DataBox(const DataBox& b, const int dim, const int indx, const int nvar) noexcept
//   #ifdef SPINER_USE_HDF
//   * DataBox(const std::string& filename);
//   * DataBox(hid_t loc, const std::string& groupname);
//   #endif // SPINER_USE_HDF
//   * void setArray(PortableMDArray<T>& A);
//   * void resize(AllocationTarget t, Args... args);
//   * void resize(Args... args);
//   * T& operator()(Args... args);
//   * T& operator()(Args... args) const;
//   * DataBox slice(const int dim, const int indx, const int nvar) const;
//   * DataBox slice(const int indx) const;
//   * DataBox slice(const int ix2, const int ix1) const;
//   * void reshape(Args... args);
//   * T interpToReal(const T x) const noexcept;
//   * T interpToReal(const T x2, const T x1) const noexcept;
//   * T interpToReal(const T x3, const T x2, const T x1) const noexcept;
//   * T interpToReal(const T x4, const T x3, const T x2, const T x1) const noexcept;
//   * void interpFromDB(const DataBox& db, const T x);
//   * void interpFromDB(const DataBox& db, const T x2, const T x1);
//   * DataBox interpToDB(Args... args);
//   * void setIndexType(int i, IndexType t);
//   * void setRange(int i, Grid_t g);
//   * void setRange(int i, Args&&... args);
//   * void copyShape(const DataBox& db, const int ndims = 0);
//   * void copyMetadata(const DataBox& src);
//   #ifdef SPINER_USE_HDF
//   * herr_t saveHDF() const;
//   * herr_t saveHDF(const std::string& filename) const;
//   * herr_t saveHDF(hid_t loc, const std::string& groupname) const;
//   * herr_t loadHDF();
//   * herr_t loadHDF(const std::string& filename);
//   * herr_t loadHDF(hid_t loc, const std::string& groupname);
//   #endif // SPINER_USE_HDF
//   * IndexType& indexType(const int i);
//   * Grid_t& range(const int i);
//   * DataBox& operator=(const DataBox& other);
//   * void copy(const DataBox& src);
//   * DataStatus dataStatus() const;
//   * bool isReference();
//   * bool ownsAllocatedMemory();
//   * bool operator==(const DataBox& other) const;
//   * bool operator!=(const DataBox& other) const;
//   * void makeShallow();
//   * void reset();
//   * T* data() const;
//   * T min() const;
//   * T max() const;
//   * int rank() const;
//   * int size() const;
//   * int sizeBytes() const;
//   * int dim() const;
//   * Grid_t range(int i) const;
//   * IndexType indexType(const int i) const;
//   * std::size_t serializedSizeInBytes() const;
//   * std::size_t serialize(char* dst) const;
//   * std::size_t setPointer(T* src);
//   * std::size_t setPointer(char* src);
//   * std::size_t deSerialize(char* src);
//   * DataBox getOnDevice() const;
//   * void finalize();


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
    PORTABLE_FUNCTION DataBoxWrapper(DataBox<DataType_t> databox)
        : databox_{databox} // MUST be initialized before Tlo_ and Thi_
        , Tlo_{get_temperature(0)}
        , Thi_{get_temperature(size() - 1)}
    {}
    PORTABLE_FUNCTION constexpr DataType_t get_logT(int const index) const
    {
        return databox_.range(0).x(index);
    }
    PORTABLE_FUNCTION constexpr DataType_t get_temperature(int const index) const
    {
        return std::exp(get_logT(index));
    }
    PORTABLE_FUNCTION constexpr DataType_t get_value(int const index) const
    {
        return databox_(index);
    }
    PORTABLE_FUNCTION constexpr DataType_t get_logV(int const index) const
    {
        return std::log(get_value(index));
    }
    PORTABLE_FUNCTION constexpr int size() const
    {
        return databox_.dim(1);
    }
    PORTABLE_FUNCTION constexpr DataType_t operator()(DataType_t const temperature) const
    {
        if (temperature < Tlo_) {
            return LowerExtrapolation::extrapolate(temperature, *this);
        } else if (temperature > Thi_) {
            return UpperExtrapolation::extrapolate(temperature, *this);
        } else {
            return databox_.interpToReal(std::log(temperature));
        }
    }
    PORTABLE_FUNCTION constexpr void finalize()
    {
        databox_.finalize();
    }
};

// ================================================================================================

} // end namespace Spiner

#endif // ifndef _SPINER_DATABOXWRAPPER_HPP_
