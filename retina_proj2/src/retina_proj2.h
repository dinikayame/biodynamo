#ifndef RETINA_PROJ2_H_
#define RETINA_PROJ2_H_

#include <vector>
#include <math.h>
#include "biodynamo.h"
#include "math_util.h"
#include "substance_initializers.h"

namespace bdm {
  /*
  part 2:
  TODO add substance to this model
  use substance as a way to control the way the cells migrate

  part 3:
  TODO test against real data for ret thickness
  correlate diameter of each cell type to what is in the program
  thickness can be tested from paraview using the axis thickness
  or set the lowest bound cell and highest bound cell and take difference
  for the thickness 
  */

/* extend into cell class to define the cell types
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

  /* List the extracellular substances

  Part 2: substance will lead cell migration
  use various concentration to lead the cell to meet

  substance behavior is stated in a biology module

  */
 enum Substances { kSubstance };

  /*
    part 3: create own function to create cells
    so all cells will start from a fixed bound and then migrate
    based on the substance concentration

    260318: revert back to this since the migration sort of works
  */

 template <typename Function, typename TResourceManager = ResourceManager<>>
 static void MyCellCreator(double min, double max, int num_cells, Function cell_builder) {
   //int num_cells = 100;
  auto rm = TResourceManager::Get();
  // Determine simulation object type which is returned by the cell_builder
  using FunctionReturnType = decltype(cell_builder({0, 0, 0}));

  auto container = rm->template Get<FunctionReturnType>();
  container->reserve(num_cells);
  //double min = 0;
  //double max_x = 520;
  //double max_y = 520;

  // so cells will be created at random only on the x and y axis
  // in this project, ignore z axis
  for (int i = 0; i < num_cells; i++) {
    double x = gTRandom.Uniform(min, max);
    double y = gTRandom.Uniform(min, max);
    //stop cells from moving in the z axis
    double z = 0;
    auto new_simulation_object = cell_builder({x, y, z});
    container->push_back(new_simulation_object);
  }
  container->Commit();
}

/*
EDITS:
210318
simplifying model to check for memory allocation error

220318
simplifying model further to just 1 cell as the gradient thing is not correct
cells only migrating along x axis and not y

1643
try adding another cell --> BUT CELL NOT MOVING
Q: should try adding substance for each cell type?

changing the concentration on the constructor changes the way the cells behave
amacrine moves but not as expected...
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

        /*
          260318: TODO add stopping criteria for cells to stop migrating
        */
        if(concentration < 0.000000002){
          cell->UpdatePosition(gradient_);
          cell->SetPosition(cell->GetMassLocation());
        }
      }
      private:
        bool init_ = false;
        DiffusionGrid* dg_ = nullptr;
        std::array<double, 3> gradient_;
    };

    /*
    Amacrine cells that are in the INNER LIMITING LAYER
    link bipolar and ganglion cells
    cell takes input from ganglion cell to bipolar cell
    So for ganglion cells --> link to the MIDGET & BISTRATIFIED CELLS
    Cells work laterally
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
          /*
            260318: TODO add stopping criteria for cells to stop migrating
          */
          if(concentration < 0.00000003){
          cell->UpdatePosition(gradient_);
          cell->SetPosition(cell->GetMassLocation());
        }
      }
    private:
      bool init_ = false;
      DiffusionGrid* dg_ = nullptr;
      std::array<double, 3> gradient_;
    };

    /*
    260318: re-add and edited
    bipolar cells
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

        if (concentration < 0.000000045) {
          cell->UpdatePosition(gradient_);
          cell->SetPosition(cell->GetMassLocation());
        }
      }
    private:
      bool init_ = false;
      DiffusionGrid* dg_ = nullptr;
      std::array<double, 3> gradient_;
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

        if (concentration < 0.00000006) {
          cell->UpdatePosition(gradient_);
          cell->SetPosition(cell->GetMassLocation());
        }
      }
    private:
      bool init_ = false;
      DiffusionGrid* dg_ = nullptr;
      std::array<double, 3> gradient_;
    };

    /*cones
    part 2:
    cells are long but rounder and wider than rods
    connect to horizontal cells + bipolar cells + Parasol cells? <check this>
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

        if (concentration < 0.000000078) {
          cell->UpdatePosition(gradient_);
          cell->SetPosition(cell->GetMassLocation());
        }
      }
    private:
      bool init_ = false;
      DiffusionGrid* dg_ = nullptr;
      std::array<double, 3> gradient_;
    };

    /*rods
    part 2:
    cells are longish -> outer seg + inner seg + nucleus
    connect to horizontal cells + bipolar cells + Parasol cells? <check this>
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

        if (concentration < 0.000000078) {
          cell->UpdatePosition(gradient_);
          cell->SetPosition(cell->GetMassLocation());
        }
      }
    private:
      bool init_ = false;
      DiffusionGrid* dg_ = nullptr;
      std::array<double, 3> gradient_;
    };

// Define compile time parameter
template <typename Backend>
struct CompileTimeParam : public DefaultCompileTimeParam<Backend> {
  using BiologyModules = Variant<ganglionCell, amacrineCell, bipolarCell, horizontalCell, coneCell, rodCell>;
  using AtomicTypes = VariadicTypedef <MyCell>;
};

inline int Simulate(int argc, const char** argv) {
  InitializeBioDynamo(argc, argv);

  /* this is to ensure that the model is reproducible for a
  specific seed
  part 2: TODO check what is a good value range for a seed to occur
  */

  Param::bound_space_ = true;
  Param::min_bound_ = 0;
  Param::max_bound_ = 520;
  Param::run_mechanical_interactions_ = true;

    int randSeed = rand() % 1000;
    gTRandom.SetSeed(randSeed);
    cout << "modelling seed: " << randSeed <<endl;

  /* Define initial model - in this example: single cell at origin
  <s>Assume if we start from outer to inner retina --> start with RPE</s>

  Part 2: should instead work from ganglion cells downwards
  assume how the light rays fall onto the retina
  */

  auto construct_ganglion = [](const std::array<double, 3>& position) {
        MyCell cell(position);
        cell.SetDiameter(20);
        cell.AddBiologyModule(ganglionCell());
        cell.SetCellType(6);
        return cell;
      };
      //CellCreator(Param::min_bound_, Param::max_bound_, 50, construct_bistratified);
    cout << "Ganglion cells created" << endl;
    MyCellCreator(Param::min_bound_, Param::max_bound_, 200, construct_ganglion);

    /*amacrine cells
    part 2: This is in the inner limiting layer where it links to the
    bipolar and ganglion cells
    so it takes the input from the ganglion cell o the bipolar cell
    So for ganglion cells --> link to the MIDGET & BISTRATIFIED CELLS
    Cells work laterally
    */
    auto construct_amacrine = [](const std::array<double, 3>& position){
        MyCell cell(position);
        cell.SetDiameter(15);
        cell.AddBiologyModule(amacrineCell());
        cell.SetCellType(5);
        return cell;
      };
    cout << "Amacrine cells created" << endl;
    MyCellCreator(Param::min_bound_, Param::max_bound_, 200, construct_amacrine);

    /*bipolar cells
    part 2: exists between photoreceptors and ganglion cells
    so it will either synpase with:
    1. photoreceptors -> rods/ cones (i.e. Parasol cells)
    2. horizontal cells
    */
      auto construct_bipolar = [](const std::array<double, 3>& position){
        MyCell cell(position);
        cell.SetDiameter(12);
        cell.AddBiologyModule(bipolarCell());
        cell.SetCellType(4);
        return cell;
      };
      cout << "Bipolar cells created" << endl;
      MyCellCreator(Param::min_bound_, Param::max_bound_, 200, construct_bipolar);

      /*horizontal cells
      part 2:
      cells work laterally
      connected to outputs from rods and cones
      */
      auto construct_horizontal = [](const std::array<double, 3>& position){
          MyCell cell(position);
          cell.SetDiameter(10);
          cell.AddBiologyModule(horizontalCell());
          cell.SetCellType(3);
          return cell;
      };
      cout << "Horizontal cells created" << endl;
      MyCellCreator(Param::min_bound_, Param::max_bound_, 200, construct_horizontal);

    /*cones
    part 2:
    cells are long but rounder and wider than rods
    connect to horizontal cells + bipolar cells + Parasol cells? <check this>
    */
      auto construct_cone = [](const std::array<double, 3>& position){
        MyCell cell(position);
        cell.SetDiameter(6);
        cell.AddBiologyModule(coneCell());
        cell.SetCellType(2);
        return cell;
      };
      cout << "Cone cells created" << endl;
      MyCellCreator(Param::min_bound_, Param::max_bound_, 100, construct_cone);

    /*rods
      part 2:
      cells are longish -> outer seg + inner seg + nucleus
      connect to horizontal cells + bipolar cells + Parasol cells? <check this>
    */
      auto construct_rod = [](const std::array<double, 3>& position){
        MyCell cell(position);
        cell.SetDiameter(5);
        cell.AddBiologyModule(rodCell());
        cell.SetCellType(1);
        return cell;
      };
      cout << "Rod cells created" << endl;
      MyCellCreator(Param::min_bound_, Param::max_bound_, 100, construct_rod);

    //defining substances in simulation
    //diffusion coefficient of 0.5, a decay constant 0f 0.1 and a resolution of 1
    ModelInitializer::DefineSubstance(kSubstance, "kSubstance", 0.5, 0.1, 4);
    //initialise substance: enum of substance, name, function type used
    //mean value of 120 along the x-axis, and a variance of 5
    //along which axis (0 = x, 1 = y, 2 = z). See the documentation of `GaussianBand` for
    //information about its arguments
    // ModelInitializer::InitializeSubstance(kSubstance, "Substance",
    //                                         GaussianBand(120, 5, Axis::kXAxis));
    ModelInitializer::InitializeSubstance(kSubstance, "kSubstance",
                                            GaussianBand(200, 100, Axis::kZAxis));


  //link to paraview to show visualization
    Param::live_visualization_ = true;
    Param::export_visualization_ = true;
    Param::visualization_export_interval_ = 2;
    //Param::visualize_sim_objects_["MyCell"] = std::set<std::string>{"diameter_"};
    Param::visualize_sim_objects_["MyCell"] = std::set<std::string>{"cellType"};

  // Run simulation for one timestep
  Scheduler<> scheduler;
  int maxStep = 1000;
  for (int i = 0; i < maxStep; i++){
    scheduler.Simulate(1);
  }


  std::cout << "Simulation completed successfully!" << std::endl;
  return 0;
}

} // namespace bdm

#endif // RETINA_PROJ2_H_
