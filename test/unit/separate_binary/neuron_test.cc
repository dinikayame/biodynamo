#include "neuroscience/neuron.h"
#include "gtest/gtest.h"

#include "compile_time_param.h"
#include "neuroscience/compile_time_param.h"

#include "unit/test_util.h"
// FIXME move to neuroscience directory

namespace bdm {

namespace neuroscience {
BDM_SIM_OBJECT(SpecializedNeurite, Neurite) {
  BDM_SIM_OBJECT_HEADER(SpecializedNeuriteExt, 1, foo_);
public:
  SpecializedNeuriteExt() {}
private:
  vec<int> foo_;
};
}

template <typename TBackend>
struct CompileTimeParam
    : public DefaultCompileTimeParam<TBackend>,
      public neuroscience::DefaultCompileTimeParam<TBackend> {
  using TNeuron = neuroscience::SpecializedNeuron;
  using TNeurite = neuroscience::SpecializedNeurite;
  using AtomicTypes = VariadicTypedef<neuroscience::SpecializedNeuron, neuroscience::SpecializedNeurite>;
};

namespace neuroscience {

TEST(NeuronTest, Scalar) {
  Neurite neurite;
  Neuron neuron;
  typename Neuron::template Self<Scalar> neuron1;
  SpecializedNeuron sneuron;
}

TEST(NeuronTest, Soa) {
  SoaNeuron neuron;
  SoaSpecializedNeuron sneuron;
  typename SpecializedNeuron::template Self<Soa> soan;
  typename CompileTimeParam<>::TNeuron soan1;
}

TEST(NeuronTest, ExtendNewNeurite) {
  auto* rm = Rm();
  Rm()->Clear();
  const double kEpsilon = abs_error<double>::value;

  // create neuron
  std::array<double, 3> origin = {0, 0, 0};
  auto neuron = rm->New<SpecializedNeuron>(origin);
  neuron.SetDiameter(20);

  // new neurite
  auto neurite = neuron.ExtendNewNeurite({0, 0, 1}).Get();
  neurite.SetDiameter(2);

  // verify
  EXPECT_ARR_NEAR(neurite.GetPosition(), {0, 0, 10.5});
  EXPECT_ARR_NEAR(neurite.GetMassLocation(), {0, 0, 11});
  EXPECT_ARR_NEAR(neurite.GetXAxis(), {0, 0, 1});
  EXPECT_ARR_NEAR(neurite.GetYAxis(), {0, 1, 0});
  EXPECT_ARR_NEAR(neurite.GetZAxis(), {-1, 0, 0});
  EXPECT_ARR_NEAR(neurite.GetSpringAxis(), {0, 0, 1});
  EXPECT_NEAR(3.1415926535897931, neurite.GetVolume(), kEpsilon);
  EXPECT_NEAR(2, neurite.GetDiameter(), kEpsilon);
  EXPECT_NEAR(0, neurite.GetBranchOrder(), kEpsilon);
  EXPECT_NEAR(1, neurite.GetActualLength(), kEpsilon);
  EXPECT_NEAR(0, neurite.GetTension(), kEpsilon);
  EXPECT_NEAR(10, neurite.GetSpringConstant(), kEpsilon);
  EXPECT_NEAR(1, neurite.GetRestingLength(), kEpsilon);
  EXPECT_TRUE(neurite.GetDaughterLeft().IsNullPtr());
  EXPECT_TRUE(neurite.GetDaughterRight().IsNullPtr());
  EXPECT_TRUE(neurite.GetMother().IsNeuron());

  EXPECT_EQ(1, rm->Get<SpecializedNeuron>()->size());
  EXPECT_EQ(1, rm->Get<SpecializedNeurite>()->size());
}

TEST(NeuronTest, ExtendNeuriteAndElongate) {
  auto* rm = Rm();
  Rm()->Clear();
  const double kEpsilon = abs_error<double>::value;
  std::array<double, 3> origin = {0, 0, 0};

  auto neuron = rm->New<SpecializedNeuron>(origin);
  neuron.SetDiameter(20);

  auto neurite_segment = neuron.ExtendNewNeurite({0, 0, 1}).Get();
  neurite_segment.SetDiameter(2);

  // will create a new neurite segment at iteration 139
  for (int i = 0; i < 200; ++i) {
    neurite_segment.ElongateTerminalEnd(10, {0, 0, 1});
    neurite_segment.RunDiscretization();
  }

  // verify
  //   distal segment
  EXPECT_ARR_NEAR(neurite_segment.GetMassLocation(), {0, 0, 31});
  EXPECT_ARR_NEAR(neurite_segment.GetPosition(), {0, 0, 27.25});
  EXPECT_ARR_NEAR(neurite_segment.GetXAxis(), {0, 0, 1});
  EXPECT_ARR_NEAR(neurite_segment.GetYAxis(), {0, 1, 0});
  EXPECT_ARR_NEAR(neurite_segment.GetZAxis(), {-1, 0, 0});
  EXPECT_ARR_NEAR(neurite_segment.GetSpringAxis(), {0, 0, 7.5});
  EXPECT_NEAR(23.561944901923749, neurite_segment.GetVolume(), kEpsilon);
  EXPECT_NEAR(2, neurite_segment.GetDiameter(), kEpsilon);
  EXPECT_NEAR(0, neurite_segment.GetBranchOrder(), kEpsilon);
  EXPECT_NEAR(7.5, neurite_segment.GetActualLength(), kEpsilon);
  EXPECT_NEAR(0, neurite_segment.GetTension(), kEpsilon);
  EXPECT_NEAR(10, neurite_segment.GetSpringConstant(), kEpsilon);
  EXPECT_NEAR(7.5, neurite_segment.GetRestingLength(), kEpsilon);
  EXPECT_TRUE(neurite_segment.GetDaughterLeft().IsNullPtr());
  EXPECT_TRUE(neurite_segment.GetDaughterRight().IsNullPtr());
  EXPECT_TRUE(neurite_segment.GetMother().IsNeurite());

  //   proximal segment
  auto proximal_segment = neurite_segment.GetMother().GetNeuriteSoPtr().Get();
  EXPECT_ARR_NEAR(proximal_segment.GetMassLocation(), {0, 0, 23.5});
  EXPECT_ARR_NEAR(proximal_segment.GetPosition(), {0, 0, 16.75});
  EXPECT_ARR_NEAR(proximal_segment.GetXAxis(), {0, 0, 1});
  EXPECT_ARR_NEAR(proximal_segment.GetYAxis(), {0, 1, 0});
  EXPECT_ARR_NEAR(proximal_segment.GetZAxis(), {-1, 0, 0});
  EXPECT_ARR_NEAR(proximal_segment.GetSpringAxis(), {0, 0, 13.5});
  EXPECT_NEAR(42.411500823462518, proximal_segment.GetVolume(), kEpsilon);
  EXPECT_NEAR(2, proximal_segment.GetDiameter(), kEpsilon);
  EXPECT_NEAR(0, proximal_segment.GetBranchOrder(), kEpsilon);
  EXPECT_NEAR(13.5, proximal_segment.GetActualLength(), kEpsilon);
  EXPECT_NEAR(0, proximal_segment.GetTension(), kEpsilon);
  EXPECT_NEAR(10, proximal_segment.GetSpringConstant(), kEpsilon);
  EXPECT_NEAR(13.5, proximal_segment.GetRestingLength(), kEpsilon);
  EXPECT_FALSE(proximal_segment.GetDaughterLeft().IsNullPtr());
  EXPECT_TRUE(proximal_segment.GetDaughterRight().IsNullPtr());
  EXPECT_TRUE(proximal_segment.GetMother().IsNeuron());

  EXPECT_EQ(1, rm->Get<SpecializedNeuron>()->size());
  EXPECT_EQ(2, rm->Get<SpecializedNeurite>()->size());
}

TEST(NeuriteTest, PartialRetraction) {
  auto* rm = Rm();
  Rm()->Clear();
  const double kEpsilon = abs_error<double>::value;
  std::array<double, 3> origin = {0, 0, 0};
  auto commit = [](auto* container, uint16_t type_idx){
    container->Commit();
  };

  auto neuron = rm->New<SpecializedNeuron>(origin);
  neuron.SetDiameter(20);

  auto neurite_segment = neuron.ExtendNewNeurite({0, 0, 1}).Get();
  neurite_segment.SetDiameter(2);

  // will create a new neurite segment at iteration 139
  for (int i = 0; i < 200; ++i) {
    neurite_segment.ElongateTerminalEnd(10, {0, 0, 1});
    neurite_segment.RunDiscretization();
  }

  // will remove the proximal segment
  for (int i = 0; i < 140; ++i) {
    neurite_segment.RetractTerminalEnd(10);
    neurite_segment.RunDiscretization();
    rm->ApplyOnAllTypes(commit);
  }

  // verify
  EXPECT_ARR_NEAR(neurite_segment.GetMassLocation(), {0, 0, 17});
  EXPECT_ARR_NEAR(neurite_segment.GetPosition(), {0, 0, 13.5});
  EXPECT_ARR_NEAR(neurite_segment.GetXAxis(), {0, 0, 1});
  EXPECT_ARR_NEAR(neurite_segment.GetYAxis(), {0, 1, 0});
  EXPECT_ARR_NEAR(neurite_segment.GetZAxis(), {-1, 0, 0});
  EXPECT_ARR_NEAR(neurite_segment.GetSpringAxis(), {0, 0, 7});
  EXPECT_NEAR(21.991148575129266, neurite_segment.GetVolume(), kEpsilon);
  EXPECT_NEAR(2, neurite_segment.GetDiameter(), kEpsilon);
  EXPECT_NEAR(0, neurite_segment.GetBranchOrder(), kEpsilon);
  EXPECT_NEAR(7, neurite_segment.GetActualLength(), kEpsilon);
  EXPECT_NEAR(0, neurite_segment.GetTension(), kEpsilon);
  EXPECT_NEAR(10, neurite_segment.GetSpringConstant(), kEpsilon);
  EXPECT_NEAR(7, neurite_segment.GetRestingLength(), kEpsilon);
  EXPECT_TRUE(neurite_segment.GetDaughterLeft().IsNullPtr());
  EXPECT_TRUE(neurite_segment.GetDaughterRight().IsNullPtr());
  EXPECT_TRUE(neurite_segment.GetMother().IsNeuron());

  EXPECT_EQ(1, rm->Get<SpecializedNeuron>()->size());
  EXPECT_EQ(1, rm->Get<SpecializedNeurite>()->size());
}

TEST(NeuriteTest, TotalRetraction) {
  auto* rm = Rm();
  Rm()->Clear();
  const double kEpsilon = abs_error<double>::value;
  std::array<double, 3> origin = {0, 0, 0};
  auto commit = [](auto* container, uint16_t type_idx){
    container->Commit();
  };

  auto neuron = rm->New<SpecializedNeuron>(origin);
  neuron.SetDiameter(20);

  auto neurite_segment = neuron.ExtendNewNeurite({0, 0, 1}).Get();
  neurite_segment.SetDiameter(2);

  // will create a new neurite segment at iteration 139
  for (int i = 0; i < 200; ++i) {
    neurite_segment.ElongateTerminalEnd(10, {0, 0, 1});
    neurite_segment.RunDiscretization();
  }

  // will remove the entire neurite
  // neurite_segment will be removed in iteration 209
  for (int i = 0; i < 210; ++i) {
    neurite_segment.RetractTerminalEnd(10);
    neurite_segment.RunDiscretization();
    rm->ApplyOnAllTypes(commit);
  }

  // verify
  EXPECT_EQ(1, rm->Get<SpecializedNeuron>()->size());
  EXPECT_EQ(0, rm->Get<SpecializedNeurite>()->size());
  EXPECT_EQ(0, neuron.GetDaughters().size());
}

TEST(NeuriteTest, Branch) {
  auto* rm = Rm();
  Rm()->Clear();
  const double kEpsilon = abs_error<double>::value;
  std::array<double, 3> origin = {0, 0, 0};
  auto commit = [](auto* container, uint16_t type_idx){
    container->Commit();
  };

  auto neuron = rm->New<SpecializedNeuron>(origin);
  neuron.SetDiameter(20);

  auto neurite_segment = neuron.ExtendNewNeurite({0, 0, 1}).Get();
  neurite_segment.SetDiameter(2);

  // will create a new neurite segment at iteration 139
  for (int i = 0; i < 200; ++i) {
    neurite_segment.ElongateTerminalEnd(10, {0, 0.5, 0.5});
    neurite_segment.RunDiscretization();
  }

  auto branch = neurite_segment.Branch({0, 1, 0}).Get();

  // verify
  //  neurite segment
  EXPECT_ARR_NEAR(neurite_segment.GetMassLocation(), {0, 14.142135623730928, 25.142135623730923});
  EXPECT_ARR_NEAR(neurite_segment.GetPosition(), {0,  12.881717786265909, 23.856717786265904});
  EXPECT_ARR_NEAR(neurite_segment.GetXAxis(), {0, 0.70012926565611089, 0.7140161142242063});
  EXPECT_ARR_NEAR(neurite_segment.GetYAxis(), {0, 0.71401611422420619, -0.70012926565611078});
  EXPECT_ARR_NEAR(neurite_segment.GetZAxis(), {-1, 0, 0});
  EXPECT_ARR_NEAR(neurite_segment.GetSpringAxis(), {0, 2.5208356749300371, 2.5708356749300378});
  EXPECT_NEAR(11.311395231915842, neurite_segment.GetVolume(), kEpsilon);
  EXPECT_NEAR(2, neurite_segment.GetDiameter(), kEpsilon);
  EXPECT_NEAR(0, neurite_segment.GetBranchOrder(), kEpsilon);
  EXPECT_NEAR(3.6005289288510043, neurite_segment.GetActualLength(), kEpsilon);
  EXPECT_NEAR(0, neurite_segment.GetTension(), kEpsilon);
  EXPECT_NEAR(10, neurite_segment.GetSpringConstant(), kEpsilon);
  EXPECT_NEAR(3.6005289288510043, neurite_segment.GetRestingLength(), kEpsilon);
  EXPECT_TRUE(neurite_segment.GetDaughterLeft().IsNullPtr());
  EXPECT_TRUE(neurite_segment.GetDaughterRight().IsNullPtr());
  EXPECT_TRUE(neurite_segment.GetMother().IsNeurite());

  //  proximal segment
  auto proximal_segment = neurite_segment.GetMother().GetNeuriteSoPtr().Get();
  EXPECT_ARR_NEAR(proximal_segment.GetMassLocation(), {0, 11.621299948800891, 22.571299948800885});
  EXPECT_ARR_NEAR(proximal_segment.GetPosition(), {0, 10.360882111335872, 21.285882111335866});
  EXPECT_ARR_NEAR(proximal_segment.GetXAxis(), {0, 0.70012926565611089, 0.7140161142242063});
  EXPECT_ARR_NEAR(proximal_segment.GetYAxis(), {0, 0.71401611422420619, -0.70012926565611078});
  EXPECT_ARR_NEAR(proximal_segment.GetZAxis(), {-1, 0, 0});
  EXPECT_ARR_NEAR(proximal_segment.GetSpringAxis(), {0, 2.5208356749300371, 2.5708356749300378});
  EXPECT_NEAR(11.311395231915842, proximal_segment.GetVolume(), kEpsilon);
  EXPECT_NEAR(2, proximal_segment.GetDiameter(), kEpsilon);
  EXPECT_NEAR(0, proximal_segment.GetBranchOrder(), kEpsilon);
  EXPECT_NEAR(3.6005289288510043, proximal_segment.GetActualLength(), kEpsilon);
  EXPECT_NEAR(0, proximal_segment.GetTension(), kEpsilon);
  EXPECT_NEAR(10, proximal_segment.GetSpringConstant(), kEpsilon);
  EXPECT_NEAR(3.6005289288510043, proximal_segment.GetRestingLength(), kEpsilon);
  EXPECT_FALSE(proximal_segment.GetDaughterLeft().IsNullPtr());
  EXPECT_FALSE(proximal_segment.GetDaughterRight().IsNullPtr());
  EXPECT_TRUE(proximal_segment.GetMother().IsNeurite());

  //  new branch
  EXPECT_ARR_NEAR(branch.GetMassLocation(), {0, 12.621299948800891, 22.571299948800885});
  EXPECT_ARR_NEAR(branch.GetPosition(), {0, 12.121299948800891, 22.571299948800885});
  EXPECT_ARR_NEAR(branch.GetXAxis(), {0, 1, 0});
  EXPECT_ARR_NEAR(branch.GetYAxis(), {0, 0, -1});
  EXPECT_ARR_NEAR(branch.GetZAxis(), {-1, 0, 0});
  EXPECT_ARR_NEAR(branch.GetSpringAxis(), {0, 1, 0});
  EXPECT_NEAR(3.1415926535897931, branch.GetVolume(), kEpsilon);
  EXPECT_NEAR(2, branch.GetDiameter(), kEpsilon);
  EXPECT_NEAR(1, branch.GetBranchOrder(), kEpsilon);
  EXPECT_NEAR(1, branch.GetActualLength(), kEpsilon);
  EXPECT_NEAR(0, branch.GetTension(), kEpsilon);
  EXPECT_NEAR(10, branch.GetSpringConstant(), kEpsilon);
  EXPECT_NEAR(1, branch.GetRestingLength(), kEpsilon);
  EXPECT_TRUE(branch.GetDaughterLeft().IsNullPtr());
  EXPECT_TRUE(branch.GetDaughterRight().IsNullPtr());
  EXPECT_TRUE(branch.GetMother().IsNeurite());

  rm->ApplyOnAllTypes(commit);
  EXPECT_EQ(1, rm->Get<SpecializedNeuron>()->size());
  EXPECT_EQ(4, rm->Get<SpecializedNeurite>()->size());
}

}  // namespace neuroscience
}  // namespace bdm


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
