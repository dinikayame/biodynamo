#ifndef RETINA_PROJ2_H_
#define RETINA_PROJ2_H_

#include <vector>
#include <math.h>
#include "biodynamo.h"
#include "math_util.h"
#include "substance_initializers.h"

namespace bdm {

/* Extend into cell class to define the cell types
*/
  BDM_SIM_OBJECT(MyCell, Cell) {
    BDM_SIM_OBJECT_HEADER(MyCellExt, 1, cellType);

  public:
    MyCellExt() {}
    MyCellExt(const std::array<double, 3>& position) : Base(position) {}

    void SetCellType(int t) { cellType[kIdx] = t; }

    int GetCellType() { return cellType[kIdx]; }
    // This function is used by ParaView for coloring the cells by their type
    int* GetCellTypePtr() { return cellType.data(); }

  private:
    vec<int> cellType;
  };

  /*
  External substance will lead cell migration based on the concentration
  levels. Substance behavior is stated in a biology module
  */
 enum Substances { kSubstance };

  /*
    Created own function to create cells
    so all cells will start from a fixed bound and then migrate
    based on the substance concentration
  */

 template <typename Function, typename TResourceManager = ResourceManager<>>
 static void MyCellCreator(double min, double max, int num_cells, Function cell_builder) {
  auto rm = TResourceManager::Get();
  // Determine simulation object type which is returned by the cell_builder
  using FunctionReturnType = decltype(cell_builder({0, 0, 0}));

  auto container = rm->template Get<FunctionReturnType>();
  container->reserve(num_cells);

  // so cells will be created at random only on the x and y axis
  // z axis is used to move cells to final resting position
  for (int i = 0; i < num_cells; i++) {
    double x = gTRandom.Uniform(min, max);
    double y = gTRandom.Uniform(min, max);
    //stop cells from moving in the z axis when generated
    double z = 0;
    auto new_simulation_object = cell_builder({x, y, z});
    container->push_back(new_simulation_object);
  }
  container->Commit();
}

/*
  BiologyModules for each cell type - specifies the behaviours of each cell type
*/

/*
  Ganglion layer cells -> Assume only 1 type for the simplicity of model
  Ganglion cells will link to amacrine and bipolar cells
*/
    struct ganglionCell : public BaseBiologyModule {
      ganglionCell() : BaseBiologyModule(gAllBmEvents){}

      template <typename T>
      void Run(T* cell){
        //initialisation of diffusion grid
        if (!init_) {
          dg_= GetDiffusionGrid(kSubstance);
          init_ = true;
        }
        //get the position of cell and grad
        auto& position = cell->GetPosition();
        dg_->GetGradient(position, &gradient_);
        double concentration = dg_->GetConcentration(position);
        //use movement to make the cells move slower
        movement[0] = gradient_[0]*0.5;
        movement[1] = gradient_[1]*0.5;
        movement[2] = gradient_[2]*0.5;

        if(concentration < 0.0000025){
          cell->UpdatePosition(movement);
          cell->SetPosition(cell->GetMassLocation());
        }
      }
      private:
        bool init_ = false;
        DiffusionGrid* dg_ = nullptr;
        std::array<double, 3> gradient_;
        std::array<double, 3> movement;
    };

    /*
    Amacrine cells that are in the INNER LIMITING LAYER
    link bipolar and ganglion cells
    cell takes input from ganglion cell to bipolar cell
    */
    struct amacrineCell : public BaseBiologyModule {
      amacrineCell() : BaseBiologyModule(gAllBmEvents){}

      template <typename T>
      void Run(T* cell){
        //initialisation of diffusion grid
        if (!init_) {
          dg_= GetDiffusionGrid(kSubstance);
          init_ = true;
        }
        //get the position of cell and grad
        auto& position = cell->GetPosition();
        dg_->GetGradient(position, &gradient_);
        double concentration = dg_->GetConcentration(position);
        //use movement to make the cells move slower
        movement[0] = gradient_[0]*0.5;
        movement[1] = gradient_[1]*0.5;
        movement[2] = gradient_[2]*0.5;

          if(concentration < 0.0000002){
          cell->UpdatePosition(movement);
          cell->SetPosition(cell->GetMassLocation());
        }
      }
    private:
      bool init_ = false;
      DiffusionGrid* dg_ = nullptr;
      std::array<double, 3> gradient_;
      std::array<double, 3> movement;
    };

    /*
    Bipolar cells
    exists between photoreceptors and ganglion cells
    so it will either synpase with:
    1. photoreceptors -> rods/ cones (i.e. Parasol cells)
    2. horizontal cells
    */
    struct bipolarCell : public BaseBiologyModule {
      bipolarCell() : BaseBiologyModule(gAllBmEvents){}

      template <typename T>
      void Run(T* cell){
        //initialisation of diffusion grid
        if (!init_) {
          dg_= GetDiffusionGrid(kSubstance);
          init_ = true;
        }
        //get the position of cell and grad
        auto& position = cell->GetPosition();
        dg_->GetGradient(position, &gradient_);
        double concentration = dg_->GetConcentration(position);
        //use movement to make the cells move slower
        movement[0] = gradient_[0]*0.5;
        movement[1] = gradient_[1]*0.5;
        movement[2] = gradient_[2]*0.5;

        if (concentration < 0.000000045) {
          cell->UpdatePosition(movement);
          cell->SetPosition(cell->GetMassLocation());
        }
      }
    private:
      bool init_ = false;
      DiffusionGrid* dg_ = nullptr;
      std::array<double, 3> gradient_;
      std::array<double, 3> movement;
    };

    /*horizontal cells
    part 2:
    cells work laterally
    connected to outputs from rods and cones
    */
    struct horizontalCell : public BaseBiologyModule {
      horizontalCell() : BaseBiologyModule(gAllBmEvents){}

      template <typename T>
      void Run(T* cell){
        //initialisation of diffusion grid
        if (!init_) {
          dg_= GetDiffusionGrid(kSubstance);
          init_ = true;
        }
        //get the position of cell and grad
        auto& position = cell->GetPosition();
        dg_->GetGradient(position, &gradient_);
        double concentration = dg_->GetConcentration(position);
        //use movement to make the cells move slower
        movement[0] = gradient_[0]*0.5;
        movement[1] = gradient_[1]*0.5;
        movement[2] = gradient_[2]*0.5;

        if (concentration < 0.000000037) {
          cell->UpdatePosition(movement);
          cell->SetPosition(cell->GetMassLocation());
        }
      }
    private:
      bool init_ = false;
      DiffusionGrid* dg_ = nullptr;
      std::array<double, 3> gradient_;
      std::array<double, 3> movement;
    };

    /*cones
    part 2:
    cells are long but rounder and wider than rods
    connect to horizontal cells + bipolar cells
    */
    struct coneCell : public BaseBiologyModule {
      coneCell() : BaseBiologyModule(gAllBmEvents){}

      template <typename T>
      void Run(T* cell){
        //initialisation of diffusion grid
        if (!init_) {
          dg_= GetDiffusionGrid(kSubstance);
          init_ = true;
        }
        //get the position of cell and grad
        auto& position = cell->GetPosition();
        dg_->GetGradient(position, &gradient_);
        double concentration = dg_->GetConcentration(position);
        //use movement to make the cells move slower
        movement[0] = gradient_[0]*0.5;
        movement[1] = gradient_[1]*0.5;
        movement[2] = gradient_[2]*0.5;

        //cones do not migrate far off from basal
        if (concentration < 0.000000002) {
          cell->UpdatePosition(movement);
          cell->SetPosition(cell->GetMassLocation());
        }
      }
    private:
      bool init_ = false;
      DiffusionGrid* dg_ = nullptr;
      std::array<double, 3> gradient_;
      std::array<double, 3> movement;
    };

    /*rods
    part 2:
    cells are longish -> outer seg + inner seg + nucleus
    connect to horizontal cells + bipolar cells
    */
    struct rodCell : public BaseBiologyModule {
      rodCell() : BaseBiologyModule(gAllBmEvents){}

      template <typename T>
      void Run(T* cell){
        //initialisation of diffusion grid
        if (!init_) {
          dg_= GetDiffusionGrid(kSubstance);
          init_ = true;
        }
        //get the position of cell and grad
        auto& position = cell->GetPosition();
        dg_->GetGradient(position, &gradient_);
        double concentration = dg_->GetConcentration(position);
        //use movement to make the cells move slower
        movement[0] = gradient_[0]*0.5;
        movement[1] = gradient_[1]*0.5;
        movement[2] = gradient_[2]*0.5;

        //rods do not migrate far off from basal
        if (concentration < 0.000000002) {
          cell->UpdatePosition(movement);
          cell->SetPosition(cell->GetMassLocation());
        }
      }
    private:
      bool init_ = false;
      DiffusionGrid* dg_ = nullptr;
      std::array<double, 3> gradient_;
      std::array<double, 3> movement;
    };

// Define compile time parameter
template <typename Backend>
struct CompileTimeParam : public DefaultCompileTimeParam<Backend> {
  using BiologyModules = Variant<ganglionCell, amacrineCell, bipolarCell, horizontalCell, coneCell, rodCell>;
  using AtomicTypes = VariadicTypedef <MyCell>;
};

inline int Simulate(int argc, const char** argv) {
  InitializeBioDynamo(argc, argv);

  /*
    Params for Paraview
  */
  Param::bound_space_ = true;
  Param::min_bound_ = 0;
  //max bound is 250um
  Param::max_bound_ = 250;
  Param::run_mechanical_interactions_ = true;

  /* this is to ensure that the model is reproducible for a
  specific seed
  */
  int randSeed = rand() % 1000;    gTRandom.SetSeed(randSeed);
  cout << "modelling seed: " << randSeed <<endl;

  /* Define cells to simulate:
  the cell diameters, cell type and name will be defined here
  */

  auto construct_ganglion = [](const std::array<double, 3>& position) {
        MyCell cell(position);
        //estimate gangliong cell diameter to be 11um
        cell.SetDiameter(11);
        cell.AddBiologyModule(ganglionCell());
        cell.SetCellType(6);
        return cell;
      };
    cout << "Ganglion cells created" << endl;
    MyCellCreator(Param::min_bound_, Param::max_bound_, 400, construct_ganglion);

    auto construct_amacrine = [](const std::array<double, 3>& position){
        MyCell cell(position);
        //assume average dendritic field to be small/medium field as we are looking
        //at central ret so around 100um
        //170418 --> ignore the dendtritic fields, estimate cell diameter
        //approx 8.75
        cell.SetDiameter(9);
        cell.AddBiologyModule(amacrineCell());
        cell.SetCellType(5);
        return cell;
      };
    cout << "Amacrine cells created" << endl;
    MyCellCreator(Param::min_bound_, Param::max_bound_, 400, construct_amacrine);

      auto construct_bipolar = [](const std::array<double, 3>& position){
        MyCell cell(position);
        //assume avg diamter using midget cone bipolar fmB and imB of 12um in 4.5mm region

        cell.SetDiameter(9);
        cell.AddBiologyModule(bipolarCell());
        cell.SetCellType(4);
        return cell;
      };
      cout << "Bipolar cells created" << endl;
      MyCellCreator(Param::min_bound_, Param::max_bound_, 400, construct_bipolar);

      auto construct_horizontal = [](const std::array<double, 3>& position){
          MyCell cell(position);
          //assume avg diamter using H1 of 25um dentritic field in 2.5mm region
          //170418 --> changed to estimated cell diameter 7.5
          cell.SetDiameter(8);
          cell.AddBiologyModule(horizontalCell());
          cell.SetCellType(3);
          return cell;
      };
      cout << "Horizontal cells created" << endl;
      MyCellCreator(Param::min_bound_, Param::max_bound_, 200, construct_horizontal);

      auto construct_cone = [](const std::array<double, 3>& position){
        MyCell cell(position);
        //approximately 2um in foveal area
        cell.SetDiameter(2);
        cell.AddBiologyModule(coneCell());
        cell.SetCellType(1);
        return cell;
      };
      cout << "Cone cells created" << endl;
      MyCellCreator(Param::min_bound_, Param::max_bound_, 250, construct_cone);

      auto construct_rod = [](const std::array<double, 3>& position){
        MyCell cell(position);
        //approximately 2umn in diameter
        cell.SetDiameter(2);
        cell.AddBiologyModule(rodCell());
        cell.SetCellType(2);
        return cell;
      };
      cout << "Rod cells created" << endl;
      MyCellCreator(Param::min_bound_, Param::max_bound_, 250, construct_rod);

    //defining substances in simulation
    //diffusion coefficient of 0.5, a decay constant 0f 0.1 and a resolution of 1
    ModelInitializer::DefineSubstance(kSubstance, "kSubstance", 0.5, 0.1, 4);
    //initialise substance: enum of substance, name, function type used
    //mean value of 200 along the z-axis, and a variance of 100
    ModelInitializer::InitializeSubstance(kSubstance, "kSubstance",
                                            GaussianBand(200, 100, Axis::kZAxis));


  //link to paraview to show visualization
    Param::live_visualization_ = true;
    Param::export_visualization_ = true;
    Param::visualization_export_interval_ = 2;
    Param::visualize_sim_objects_["MyCell"] = std::set<std::string>{"cellType"};

  // Run simulation for one timestep
  Scheduler<> scheduler;
  int maxStep = 1800;
  for (int i = 0; i < maxStep; i++){
    scheduler.Simulate(1);
  }


  std::cout << "Simulation completed successfully!" << std::endl;
  return 0;
}

} // namespace bdm

#endif // RETINA_PROJ2_H_
