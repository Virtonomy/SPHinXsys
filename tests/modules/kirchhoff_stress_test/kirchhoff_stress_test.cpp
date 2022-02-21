#include <gtest/gtest.h>
#include "structural_simulation_class.h"

void MockStent_2()
{
	Real scale_stl = 0.001;
	Real end_time = 0.08;
	/* MATERIAL PARAMETERS */
	Real rho_0 = 1000.0;
	Real poisson = 0.35;
	Real physical_viscosity = 200;
	/** STL IMPORT PARAMETERS */
	string relative_input_path = "./input/"; //path definition for linux
	vector<string> imported_stl_list = { "mock_stent.stl", "vessel_cylinder.stl" };
	vector<Vec3d> translation_list = { Vec3d(0), Vec3d(0) };
	vector<Real> resolution_list = { 3, 1 };
	LinearElasticSolid material = LinearElasticSolid(rho_0, 1e4, poisson);
	LinearElasticSolid material_stiff = LinearElasticSolid(rho_0, 1e6, poisson);
	StdVec<shared_ptr<LinearElasticSolid>> material_model_list = { make_shared<LinearElasticSolid>(material_stiff), make_shared<LinearElasticSolid>(material) };
	StructuralSimulationInput input
	{
		relative_input_path,
		imported_stl_list,
		scale_stl,
		translation_list,
		resolution_list,
		material_model_list,
		StdVec<Real>(3, physical_viscosity),
		{}
	};
	pair<array<int, 2>, array<Real, 2>> contact_pair_1( {0, 1}, {0.001, end_time} );
	input.time_dep_contacting_body_pairs_list_ = { contact_pair_1 };

	input.non_zero_gravity_ = { GravityPair(0, Vec3d(0, 75, 0)) };

	// spring damper constraint
	Real spring_stiffness = 100;
	Real damping_ratio = 0.02;
	input.spring_damper_tuple_ = { SpringDamperTuple(1, Vec3d(spring_stiffness), damping_ratio) };

	/** SIMULATION MODEL */
	StructuralSimulation sim(input);

	try
	{
		sim.runSimulation(end_time);
	}
	catch(const std::exception& e)
	{
		std::cout << e.what() << std::endl;
	}

	StdLargeVec<Vecd> pos_n = sim.get_solid_body_list_()[1]->getElasticSolidParticles()->pos_n_;

	Real max_y = 0.0;
	for (auto vec: pos_n)
		if (vec[1] > max_y) max_y = vec[1];
	std::cout << " ============ max_y: " << max_y << std::endl;

	// EXPECT_EQ(sim.get_contacting_body_pairs_list_().size(), 0);
	// EXPECT_EQ(sim.get_time_dep_contacting_body_pairs_list_().size(), 1);
	// EXPECT_EQ(sim.get_contact_list_().size(), 2);
	// EXPECT_EQ(sim.get_contact_density_list_().size(), 2);
	// EXPECT_EQ(sim.get_contact_force_list_().size(), 2);
	
	// /** START SIMULATION */
	// try
	// {
	// 	sim.setKirchhofCorrection(1, correction_scalar);
	// 	sim.TestRunSimulation(end_time);
		
	// 	std::vector<Real> max_strain_over_time = sim.getVonMisesStrainMax();
	// 	std::cout << "max_strain_over_time.size(): " << max_strain_over_time.size() << std::endl;
	// 	std::cout << "max_strain_over_time.back(): " << max_strain_over_time.back() << std::endl;

	// 	return {max_strain_over_time.size(), max_strain_over_time.back()};
	// }
	// catch(const std::exception& e)
	// {
	// 	std::vector<Real> max_strain_over_time = sim.getVonMisesStrainMax();
	// 	std::cout << "max_strain_over_time.size(): " << max_strain_over_time.size() << std::endl;
	// 	std::cout << "max_strain_over_time.back(): " << max_strain_over_time.back() << std::endl;

	// 	return {max_strain_over_time.size(), max_strain_over_time.back()};
	// }

	
}

TEST(StructuralSimulation, kirchhoff_stress_test)
{
	// StdVec<Real> correction_values = {0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5};
	// StdVec<Real> full_corr_list = {};
	// StdVec<Real> steps = {};
	// StdVec<Real> strains = {};

	// for (int i = 0; i<3; ++i)
	// 	for (Real correction: correction_values)
	// 	{
	// 		std::array<Real, 2> res = MockStent_2(correction);
	// 		full_corr_list.push_back(correction);
	// 		steps.push_back(res[0]);
	// 		strains.push_back(res[1]);
	// 	}

	// matplotlibcpp::figure();
	// //matplotlibcpp::figure_size(1200, 780);
	// matplotlibcpp::scatter(full_corr_list, steps);
	// std::string plot_name_1 = "Steps.png";
	// matplotlibcpp::save(plot_name_1);

	// matplotlibcpp::figure();
	// //matplotlibcpp::figure_size(1200, 780);
	// matplotlibcpp::scatter(full_corr_list, strains);
	// std::string plot_name_2 = "Strains.png";
	// matplotlibcpp::save(plot_name_2);

	MockStent_2();
}

int main(int argc, char* argv[])
{	
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
