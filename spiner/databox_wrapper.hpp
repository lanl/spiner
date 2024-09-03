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

// TODO:
// * In theory the extra features could be split:
//   * One wrapper that adds transformations
//   * One wrapper that adds extrapolation
// * If you want both, you can have Extrapolation<Transformation<DataBox<T>>>
// * But would it also work if you did Transformation<Extrapolation<DataBox<T>>>?
// * This idea seems overly complicated for little benefit.

// TODO: Name?
// DataBoxPlusPlus
// ExtendedDataBox
// DataBoxButCooler
// FancyPantsDataBox
// DataBoxNowWithMoreFeatures
// DataBoxBetterStrongerFaster

// TODO: Template defaults?
// * The obvious/naive defaults would be to match an unwrapped DataBox
//   * All extrapolations default to "generate an error"
//   * All transformations are the identity transformation
// * But given that the entire purpose of this wrapper is to add features beyond that of the basic
//   DataBox, does it make sense for it to have defaults that make this a DataBox but with extra
//   layers of complexity?
// * It's also not clear what order the template arguments should be in to manage defaults.
template <
    typename DataType_t,
    // TODO: Right now all axes will use the same lower and upper extrapolations, which may or may
    //       not be the desired behavior.
    typename LowerExtrapolation, typename UpperExtrapolation,
    typename TransformX, typename TransformY>
class DataBoxWrapper
{
private:
    // TODO: Deal with the other templates.
    using DB_t = DataBox<DataType_t>;
    DB_t databox_;
    // Avoid re-calculating independent variable bounds on every call
    // TODO: My original implementation in Singe did this, but after thinking about it I think we
    //       should throw this out.
    //       * When interpolating, we already have to transform the independent variable(s), so no
    //         extra work.
    //       * When extrapolating, we don't use the transformed independent variable(s), so we do
    //         an extra transformation.
    //       * Finding the _exact_ bounds, that's going to be defined by the underlying DataBox,
    //         which is in transformed space.  If we do bounds checking in the untransformed space,
    //         we have the possibility of numerical error creating a small space between the
    //         untransformed bound and the transformed bound, which could create errors.
    DataType_t [DB_t::MAXRANK] lower_bounds_;
    DataType_t [DB_t::MAXRANK] upper_bounds_;

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
