#include "biology_module_util.h"
#include "cell.h"
#include "resource_manager.h"
#include "scheduler.h"
#include "simulation_object_util.h"

namespace bdm {

// 1. Define growth behaviour
struct GrowthModule {
  template <typename T>
  void Run(T* cell) {
    if (cell->GetDiameter() <= 40) {
      cell->ChangeVolume(300);
    } else {
      Divide(*cell, ResourceManager<SoaCell>::Get()->Get<SoaCell>());
    }
  }

  bool IsCopied(Event event) const { return true; }
};

// 2. Define biology modules that should be used in this simulation
//    Simulation objects automatically pick up this definition
BDM_DEFINE_BIOLOGY_MODULES(GrowthModule);

void Simulate(size_t cells_per_dim = 128) {
  // 3. Get cell container
  auto cells = ResourceManager<SoaCell>::Get()->Get<SoaCell>();
  cells->reserve(cells_per_dim * cells_per_dim * cells_per_dim);

  // 4. Define initial model - in this case 3D grid of cells
  double space = 20;
  for (size_t i = 0; i < cells_per_dim; i++) {
    for (size_t j = 0; j < cells_per_dim; j++) {
      for (size_t k = 0; k < cells_per_dim; k++) {
        Cell cell({i * space, j * space, k * space});
        cell.SetDiameter(30);
        cell.SetAdherence(0.4);
        cell.SetMass(1.0);
        cell.UpdateVolume();
        cell.AddBiologyModule(GrowthModule());
        cells->push_back(cell);
      }
    }
  }

  // 5. Run simulation for one timestep
  Scheduler scheduler;
  scheduler.Simulate<SoaCell>(1);
}

}  // namespace bdm

int main() { bdm::Simulate(); }
