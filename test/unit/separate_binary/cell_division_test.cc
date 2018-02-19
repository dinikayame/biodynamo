#include "unit/separate_binary/cell_division_test.h"

namespace bdm {

// Disable test until the trello card (Fix inconsistency in cell state due to 
// direct updates in Biology Modules) is fixed
// TEST(CellDivisionTest, ComputeSoa) { RunTestSoa(); }

}  // namespace bdm

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

