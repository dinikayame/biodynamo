#ifndef UNIT_SEPARATE_BINARY_CELL_DIVISION_TEST_H_
#define UNIT_SEPARATE_BINARY_CELL_DIVISION_TEST_H_

#include "biology_module_op.h"
#include "cell.h"
#include "compile_time_param.h"
#include "gtest/gtest.h"
#include "transactional_vector.h"
#include "debug.h"

namespace bdm {

inline void EXPECT_ARR_NEAR(const std::array<double, 3>& actual, const std::array<double, 3>& expected) {
  for (size_t i = 0; i < actual.size(); i++) {
    EXPECT_NEAR(expected[i], actual[i], 1e-9);
  }
}

struct GrowDivide : public BaseBiologyModule {
  double growth_rate_;
  double threshold_;

  GrowDivide() : growth_rate_(5000), threshold_(30.05) {}
  explicit GrowDivide(double growth_rate) : growth_rate_(growth_rate) {}

  template <typename T>
  void Run(T* cell) {
    if (cell->GetDiameter() <= threshold_) {
      cell->ChangeVolume(growth_rate_);
    } else {
      Divide(*cell);
    }
  }
};

template <typename Backend>
struct CompileTimeParam : public DefaultCompileTimeParam<Backend> {
  // use predefined biology module GrowDivide
  using BiologyModules = Variant<GrowDivide>;
};

template <typename T>
inline void RunTest(T* cells) {
  size_t cells_per_dim = 2;
  auto construct = [](const std::array<double, 3>& position) {
    Cell cell(position);
    cell.SetDiameter(30);
    cell.SetAdherence(0.4);
    cell.SetMass(1.0);
    cell.AddBiologyModule(GrowDivide());
    return cell;
  };

  for (size_t x = 0; x < cells_per_dim; x++) {
    double x_pos = x * 20.0;
    for (size_t y = 0; y < cells_per_dim; y++) {
      double y_pos = y * 20.0;
      for (size_t z = 0; z < cells_per_dim; z++) {
        auto new_simulation_object = construct({x_pos, y_pos, z * 20.0});
        cells->push_back(new_simulation_object);
      }
    }
  }

  BiologyModuleOp op;

  // At step 2 a division should take place; In step 3 these new cells are
  // created
  for (int time = 0; time < 10; time++) {
    op(cells, 0);
  }

  EXPECT_ARR_NEAR((*cells)[0].GetPosition(), {-6.22335511811398944815, -1.66604831760015636988, 1.81231512675687889136});
  EXPECT_ARR_NEAR((*cells)[1].GetPosition(), {5.6398633086479312837, 0.956959214342484321136, 23.6096212621024079681});
  EXPECT_ARR_NEAR((*cells)[2].GetPosition(), {5.83721343145111593032, 22.1469771658114389368, -3.92925001658284855921});
  EXPECT_ARR_NEAR((*cells)[3].GetPosition(), {-2.65380535403253503546, 22.6212948300807141777, 14.2130203890232333919});
  EXPECT_ARR_NEAR((*cells)[4].GetPosition(), {17.0724094351280903936, 2.53835632812800193747, 6.15967749034800338137});
  EXPECT_ARR_NEAR((*cells)[5].GetPosition(), {16.5005638568587400528, 2.75317893380783740298, 25.8993128689914371421});
  EXPECT_ARR_NEAR((*cells)[6].GetPosition(), {22.5013890209255080777, 25.2198353273051303347, -3.39034181797633227262});
  EXPECT_ARR_NEAR((*cells)[7].GetPosition(), {24.119543938806053518, 23.7508536867231256906, 24.4832735992416345994});
}

inline void RunTestSoa() {
  auto cells = Cell::NewEmptySoa();
  RunTest(&cells);
}

}  // namespace bdm

#endif  // UNIT_SEPARATE_BINARY_CELL_DIVISION_TEST_H_
