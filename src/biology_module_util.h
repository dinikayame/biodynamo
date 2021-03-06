#ifndef BIOLOGY_MODULE_UTIL_H_
#define BIOLOGY_MODULE_UTIL_H_

#include <TError.h>
#include <limits>
#include <mutex>
#include "variant.h"

namespace bdm {

/// BmEvent is used inside biology modules to determine if a biology module
/// should be copied if a new simulation object is created.
/// Possible events are cell division, neurite branching, ...\n
/// BmEvents invariant: the number of bits set to 1 must be 1.
///
///     // This is how a new event can be defined:
///     // Declare new event in header file
///     extern const BmEvent gCellDivision;
///     // Define it in source file
///     const BmEvent gCellDivision =
///     UniqueBmEventFactory::Get()->NewUniqueBmEvent();
using BmEvent = uint64_t;

/// Biology module event representing the union of all events.\n
/// Used to create a biology module  which is copied for every event.
/// @see `BaseBiologyModule`
const BmEvent gAllBmEvents = std::numeric_limits<uint64_t>::max();

/// This class generates unique ids for biology module events satisfying the
/// BmEvent invariant. Thread safe.
class UniqueBmEventFactory {
 public:
  UniqueBmEventFactory(const UniqueBmEventFactory&) = delete;

  static UniqueBmEventFactory* Get() {
    static UniqueBmEventFactory kInstance;
    return &kInstance;
  }

  BmEvent NewUniqueBmEvent() {
    std::lock_guard<std::recursive_mutex> lock(mutex_);
    constexpr uint64_t kOne = 1;
    if (counter_ == 64) {
      Fatal("UniqueBmEventFactory",
            "BioDynaMo only supports 64 unique BmEvents."
            " You requested a 65th one.");
    }
    return kOne << counter_++;
  }

 private:
  UniqueBmEventFactory() {}
  std::recursive_mutex mutex_;
  uint64_t counter_ = 0;
};

/// BaseBiologyModule encapsulates logic to decide for which BmEvents
/// a biology module should be copied.
struct BaseBiologyModule {
  /// Default ctor sets copy_mask_ to 0; meaning that IsCopied will always
  /// return false
  BaseBiologyModule() : copy_mask_(0) {}
  explicit BaseBiologyModule(BmEvent event) : copy_mask_(event) {}
  BaseBiologyModule(std::initializer_list<BmEvent> events) {
    copy_mask_ = 0;
    for (BmEvent event : events) {
      copy_mask_ |= event;
    }
  }

  BaseBiologyModule(const BaseBiologyModule& other)
      : copy_mask_(other.copy_mask_) {}

  /// Function returns whether the biology module should be copied for the
  /// given event.
  bool IsCopied(BmEvent event) const { return (event & copy_mask_) != 0; }

 private:
  BmEvent copy_mask_;
  ClassDefNV(BaseBiologyModule, 1);
};

/// \brief Used for simulation objects where biology modules are not used.
/// Variant implementation does not allow `Variant<>`
/// -> `Variant<NullBiologyModule>`
struct NullBiologyModule : public BaseBiologyModule {
  template <typename T>
  void Run(T* t) {}

  ClassDefNV(NullBiologyModule, 1);
};

/// \brief Visitor to execute the `Run` method of a biology module
/// @tparam TSimulationObject type of simulation object that owns the biology
///         module
template <typename TSimulationObject>
struct RunVisitor {
  /// @param so pointer to the simulation object on which the biology module
  ///        should be executed
  explicit RunVisitor(TSimulationObject* const so) : kSimulationObject(so) {}

  template <typename T>
  void operator()(T& t) const {
    t.Run(kSimulationObject);
  }

 private:
  TSimulationObject* const kSimulationObject;
};

/// \brief Visitor to copy biology modules from one structure to another
/// @tparam TVector type of the destination container where the biology modules
///         are stored
template <typename TVector>
struct CopyVisitor {
  /// @param event that lead to the copy operation - e.g. cell division:
  ///        biology modules should be copied from mother to daughter cell
  /// @param vector biology module vector in which the copied module will be
  ///        inserted
  CopyVisitor(BmEvent event, TVector* vector)
      : kEvent(event), vector_(vector) {}

  template <typename T>
  void operator()(const T& from) const {
    if (from.IsCopied(kEvent)) {
      T copy(from);  // NOLINT
      vector_->emplace_back(std::move(copy));
    }
  }

  const BmEvent kEvent;
  TVector* vector_;
};

}  // namespace bdm

#endif  // BIOLOGY_MODULE_UTIL_H_
