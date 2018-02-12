#ifndef PARAM_H_
#define PARAM_H_

#include <cinttypes>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>
#include "cpptoml/cpptoml.h"

namespace bdm {

struct Param {

  // simulation values ---------------------------------------------------------

  ///protein_amount specifies the amount of different proteins in the simulation
  static const int protein_amount = 3;

  ///Lambda function which contains functions by which GeneCalculation module calculate concentration of proteins on each step
  static std::function<std::vector<double>(double, std::vector<std::array<double, protein_amount>>)> functions;

  /// Variable which is tracking current step of simulation
  static size_t step_global_;

  ///Variable which specifies method using for solving differential equation {"Euler", "RK4"}.
  static std::string dE_solve_method;

  /// Backup file name for full simulation backups\n
  /// Default value: `""` (no backups will be made)\n
  /// TOML config file:
  ///
  ///     [simulation]
  ///     backup_file = <path>/<filename>.root
  /// Command line argument: `-b, --backup`
  static std::string backup_file_;

  /// File name to restore simulation from\n
  /// Default value: `""` (no restore will be made)\n
  /// TOML config file:
  ///
  ///     [simulation]
  ///     restore_file = <path>/<filename>.root
  /// Command line argument: `-r, --restore`
  static std::string restore_file_;

  /// Specifies the interval (in seconds) in which backups will be performed.\n
  /// Default Value: `1800` (every half an hour)\n
  /// TOML config file:
  ///
  ///     [simulation]
  ///     backup_interval = 1800  # backup every half an hour
  static uint32_t backup_interval_;

  /// Time between two simulation steps, in hours.
  /// Default value: `0.01`\n
  /// TOML config file:
  ///
  ///     [simulation]
  ///     time_step = 0.0125
  static double simulation_time_step_;

  /// Maximum jump that a point mass can do in one time step. Useful to
  /// stabilize the simulation\n
  /// Default value: `3.0`\n
  /// TOML config file:
  ///
  ///     [simulation]
  ///     max_displacement = 3.0
  static double simulation_max_displacement_;

  /// Calculate mechanical interactions between simulation objects.\n
  /// Default value: `true`\n
  /// TOML config file:
  ///
  ///     [simulation]
  ///     run_mechanical_interactions = true
  static bool run_mechanical_interactions_;

  /// Enforce an artificial cubic bounds around the simulation space.
  /// Simulation objects cannot move outside this cube. Dimensions of this cube
  /// are determined by parameter `lbound` and `rbound`.\n
  /// Default value: `false` (simulation space is "infinite")\n
  /// TOML config file:
  ///
  ///     [simulation]
  ///     bound_space = false
  static bool bound_space_;

  /// Minimum allowed value for x-, y- and z-position if simulation space is
  /// bound (@see `bound_space_`).\n
  /// Default value: `0`\n
  /// TOML config file:
  ///
  ///     [simulation]
  ///     min_bound = 0
  static double min_bound_;

  /// Maximum allowed value for x-, y- and z-position if simulation space is
  /// bound (@see `bound_space_`).\n
  /// Default value: `100`\n
  /// TOML config file:
  ///
  ///     [simulation]
  ///     max_bound = 100
  static double max_bound_;

  // visualization values ------------------------------------------------------

  /// Use ParaView Catalyst for live visualization.\n
  /// Defaut value: `false`\n
  /// TOML config file:
  ///
  ///     [visualization]
  ///     live = false
  static bool live_visualization_;

  /// Write data to file for post-simulation visualization
  /// Defaut value: `false`\n
  /// TOML config file:
  ///
  ///     [visualization]
  ///     export = false
  static bool export_visualization_;

  /// If `export_visualization_` is set to true, this parameter specifies
  /// how often it should be exported. 1 = every timestep, 10: every 10
  /// time steps.\n
  /// Defaut value: `1`\n
  /// TOML config file:
  ///
  ///     [visualization]
  ///     export_interval = 1
  static uint32_t visualization_export_interval_;

  /// Specifies which simulation objects should be visualized. \n
  /// Every simulation object defines the minimum set of data members which
  /// are required to visualize it. (e.g. Cell: `position_` and `diameter_`).\n
  /// With this parameter it is also possible to extend the number of data
  /// members that are sent to the visualization engine.
  /// Default value: empty (no simulation object will be visualized)\n
  /// TOML config file:
  ///
  ///     [visualization]
  ///     # turn on live or export
  ///     export = true
  ///
  ///       [[visualize_sim_object]]
  ///       name = "Cell"
  ///       # the following entry is optional
  ///       additional_data_members = [ "density_" ]
  ///
  ///       # The former block can be repeated for further simulation objects
  ///       [[visualize_sim_object]]
  ///       name = "Neurite"
  static std::unordered_map<std::string, std::set<std::string>>
      visualize_sim_objects_;

  struct VisualizeDiffusion {
    std::string name_;
    bool concentration_ = true;
    bool gradient_ = false;
  };

  /// Spefifies if for which substances extracellular diffusion should be
  /// visualized.\n
  /// Default value: empty (no diffusion will be visualized)\n
  /// TOML config file:
  ///
  ///     [visualization]
  ///     # turn on live or export
  ///     export = true
  ///
  ///       [[visualize_diffusion]]
  ///       # Name of the substance
  ///       name = "Na"
  ///       # the following two entries are optional
  ///       #   default value for concentration is true
  ///       concentration = true
  ///       #   default value for gradient is false
  ///       gradient = false
  ///
  ///       # The former block can be repeated for further substances
  ///       [[visualize_diffusion]]
  ///       name = "K"
  ///       # default values: concentration = true and gradient = false
  static std::vector<VisualizeDiffusion> visualize_diffusion_;

  // development values --------------------------------------------------------
  /// Statistics of profiling data; keeps track of the execution time of each
  /// operation at every timestep.\n
  /// Default Value: `false`\n
  /// TOML config file:
  ///
  ///     [development]
  ///     statistics = false
  static bool statistics_;

  /// Use the python script (simple_pipeline.py) to do Live Visualization with
  /// ParaView. If false, we use the C++ pipeline
  /// Defautl value: `false`\n
  /// TOML config file:
  ///     [development]
  ///     python_catalyst_pipeline_ = false
  static bool python_catalyst_pipeline_;

  /// Resets the static variables to its default values
  static void Reset();

 private:
  friend void InitializeBioDynamo(int, const char**);

  /// Assign values from config file to static variables
  static void AssignFromConfig(const std::shared_ptr<cpptoml::table>&);
};

}  // namespace bdm

#endif  // PARAM_H_
