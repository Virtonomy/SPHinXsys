#include <gtest/gtest.h>
#include "solid_structural_simulation_class.h"

Real tolerance = 1e-6;

class TestStructuralSimulation : StructuralSimulation
{
public:
	TestStructuralSimulation(StructuralSimulationInput& input) : StructuralSimulation(input){};
	void TestRunSimulation(Real end_time){ RunSimulation(end_time); };
	// input members
	string Get_relative_input_path_(){ return relative_input_path_; };
	vector<string> Get_imported_stl_list_(){ return imported_stl_list_; };
	Real Get_scale_stl_(){ return scale_stl_; };
	vector<Vec3d> Get_translation_list_(){ return translation_list_; };
	Real Get_system_resolution_(){ return system_resolution_; };
	vector<Real> Get_resolution_list_(){ return resolution_list_; };
	vector<LinearElasticSolid> Get_material_model_list_(){ return material_model_list_; };
	Real Get_physical_viscosity_(){ return physical_viscosity_; };
	// internal members
	SPHSystem Get_system_(){ return system_; };
	Real Get_scale_system_boundaries_(){ return scale_system_boundaries_; };
	In_Output Get_in_output_(){ return in_output_; };
	// other
	vector<TriangleMeshShape> Get_body_mesh_list_(){ return body_mesh_list_; };
	vector<shared_ptr<SolidBodyForSimulation>> Get_solid_body_list_(){ return solid_body_list_; };
	vector<array<int, 2>> Get_contacting_body_pairs_list_(){ return contacting_body_pairs_list_; };
	vector<pair<array<int, 2>, array<Real, 2>>> Get_time_dep_contacting_body_pairs_list_(){ return time_dep_contacting_body_pairs_list_; };
	vector<shared_ptr<SolidContactBodyRelation>> Get_contact_list_(){ return contact_list_; };
	vector<shared_ptr<solid_dynamics::ContactDensitySummation>> Get_contact_density_list_(){ return contact_density_list_; };
	vector<shared_ptr<solid_dynamics::ContactForce>> Get_contact_force_list_(){ return contact_force_list_; };
	// for InitializeATimeStep
	vector<shared_ptr<TimeStepInitialization>> Get_initialize_gravity_(){ return initialize_gravity_; };
	vector<GravityPair> Get_non_zero_gravity_(){ return non_zero_gravity_; };
	// for AccelerationForBodyPartInBoundingBox
	vector<shared_ptr<solid_dynamics::AccelerationForBodyPartInBoundingBox>> Get_acceleration_bounding_box_(){ return acceleration_bounding_box_; };
	vector<AccelTuple> Get_acceleration_bounding_box_tuple_(){ return acceleration_bounding_box_tuple_; };
	// for SpringDamperConstraintParticleWise
	vector<shared_ptr<solid_dynamics::SpringDamperConstraintParticleWise>> Get_spring_damper_constraint_(){ return spring_damper_constraint_; };
	vector<SpringDamperTuple> Get_spring_damper_tuple_(){ return spring_damper_tuple_; };
	// for ConstrainSolidBodyRegion
	vector<shared_ptr<solid_dynamics::ConstrainSolidBodyRegion>> Get_fixed_constraint_(){ return fixed_constraint_; };
	vector<int> Get_body_indeces_fixed_constraint_(){ return body_indeces_fixed_constraint_; };
	// for PositionSolidBody
	vector<shared_ptr<solid_dynamics::PositionSolidBody>> Get_position_solid_body_(){ return position_solid_body_; };
	vector<PositionSolidBodyTuple> Get_position_solid_body_tuple_(){ return position_solid_body_tuple_; };
	// for PositionScaleSolidBody
	vector<shared_ptr<solid_dynamics::PositionScaleSolidBody>> Get_position_scale_solid_body_(){ return position_scale_solid_body_; };
	vector<PositionScaleSolidBodyTuple> Get_position_scale_solid_body_tuple_(){ return position_scale_solid_body_tuple_; };
};

TEST(StructuralSimulation, TimeDependentContact)
{
	Real scale_stl = 0.001;
	Real res = 1.5;
	Real end_time = 0.1;
	/* MATERIAL PARAMETERS */
	Real rho_0 = 1000.0;
	Real poisson = 0.35;
	Real Youngs_modulus = 1e5;
	Real physical_viscosity = 200;
	/** STL IMPORT PARAMETERS */
	string relative_input_path = "./input/"; //path definition for linux
	string plate_stl = "plate.stl";
	string plate2_stl = "plate2.stl";
	vector<string> imported_stl_list = { plate_stl, plate2_stl };
	vector<Vec3d> translation_list = { Vec3d(0), Vec3d(0) };
	vector<Real> resolution_list = { res, res };
	LinearElasticSolid material = LinearElasticSolid(rho_0, Youngs_modulus, poisson);
	vector<LinearElasticSolid> material_model_list = { material, material };
	StructuralSimulationInput input
	{
		relative_input_path,
		imported_stl_list,
		scale_stl,
		translation_list,
		resolution_list,
		material_model_list,
		physical_viscosity,
		{}
	};
	pair<array<int, 2>, array<Real, 2>> contact_pair( {0, 1}, {end_time * 0.5, end_time} );
	input.time_dep_contacting_body_pairs_list_ = { contact_pair };
	// position the plate
	input.position_solid_body_tuple_ = { PositionSolidBodyTuple(1, 0.0, end_time * 0.5, Vec3d(0, 0, 10.0) * scale_stl),
											PositionSolidBodyTuple(1, end_time * 0.5, end_time, Vec3d(0, 0, 0.0) * scale_stl) };
	
	// spring damper constraint
	Real spring_stiffness = 500;
	Real damping_ratio = 0.02;
	input.spring_damper_tuple_ = { SpringDamperTuple(0, Vec3d(spring_stiffness), damping_ratio)};

	/** SIMULATION MODEL */
	TestStructuralSimulation sim(input);
	/** START SIMULATION */
	sim.TestRunSimulation(end_time);

	// check that the two plates don't overlap, this will check if the contact works properly
	StdLargeVec<Vecd>& pos_n_1 = sim.Get_solid_body_list_()[0].get()->GetElasticSolidParticles()->pos_n_; // plate 1 particles
	StdLargeVec<Vecd>& pos_n_2 = sim.Get_solid_body_list_()[1].get()->GetElasticSolidParticles()->pos_n_; // plate 2 particles

	// plate 1 max z-coordinate
	Real max_z_1 = -1000.0;
	for (size_t index = 0; index < pos_n_1.size(); index++)
	{
		if (max_z_1 < pos_n_1[index][2]) max_z_1 = pos_n_1[index][2];
	}
	// plate 2 min z-coordinate
	Real min_z_2 = 1000.0;
	for (size_t index = 0; index < pos_n_1.size(); index++)
	{
		if (min_z_2 > pos_n_2[index][2]) min_z_2 = pos_n_2[index][2];
	}
	// check that min_z_2, minimum z of plate 1 is larger than max_z_1, maximum z of plate 2
	EXPECT_GT(min_z_2, max_z_1);
}

TEST(StructuralSimulation, TimeDependentContactStiff)
{
	Real scale_stl = 0.001;
	Real res = 1.5;
	Real end_time = 0.04;
	/* MATERIAL PARAMETERS */
	Real rho_0 = 1000.0;
	Real poisson = 0.35;
	Real physical_viscosity = 200;
	/** STL IMPORT PARAMETERS */
	string relative_input_path = "./input/"; //path definition for linux
	string plate_stl = "plate.stl";
	string plate2_stl = "plate2.stl";
	vector<string> imported_stl_list = { plate_stl, plate2_stl };
	vector<Vec3d> translation_list = { Vec3d(0), Vec3d(0, 0, 8) };
	vector<Real> resolution_list = { res, res };
	LinearElasticSolid material = LinearElasticSolid(rho_0, 1e4, poisson);
	LinearElasticSolid material_stiff = LinearElasticSolid(rho_0, 1e8, poisson);
	vector<LinearElasticSolid> material_model_list = { material, material_stiff };
	StructuralSimulationInput input
	{
		relative_input_path,
		imported_stl_list,
		scale_stl,
		translation_list,
		resolution_list,
		material_model_list,
		physical_viscosity,
		{}
	};
	pair<array<int, 2>, array<Real, 2>> contact_pair( {0, 1}, {end_time * 0.05, end_time} );
	input.time_dep_contacting_body_pairs_list_ = { contact_pair };
	// position the plate
	// input.position_solid_body_tuple_ = { PositionSolidBodyTuple(1, 0.0, end_time * 0.5, Vec3d(0, 0, 10.0) * scale_stl),
	// 										PositionSolidBodyTuple(1, end_time * 0.5, end_time, Vec3d(0, 0, 0.0) * scale_stl) };
	
	input.non_zero_gravity_ = { GravityPair(1, Vec3d(0, 0, -100)) };

	// spring damper constraint
	Real spring_stiffness = 500;
	Real damping_ratio = 0.02;
	input.spring_damper_tuple_ = { SpringDamperTuple(0, Vec3d(spring_stiffness), damping_ratio)};

	/** SIMULATION MODEL */
	TestStructuralSimulation sim(input);
	/** START SIMULATION */
	sim.TestRunSimulation(end_time);

	// check that the two plates don't overlap, this will check if the contact works properly
	StdLargeVec<Vecd>& pos_n_1 = sim.Get_solid_body_list_()[0].get()->GetElasticSolidParticles()->pos_n_; // plate 1 particles
	StdLargeVec<Vecd>& pos_n_2 = sim.Get_solid_body_list_()[1].get()->GetElasticSolidParticles()->pos_n_; // plate 2 particles

	// plate 1 max z-coordinate
	Real max_z_1 = -1000.0;
	for (size_t index = 0; index < pos_n_1.size(); index++)
	{
		if (max_z_1 < pos_n_1[index][2]) max_z_1 = pos_n_1[index][2];
	}
	// plate 2 min z-coordinate
	Real min_z_2 = 1000.0;
	for (size_t index = 0; index < pos_n_1.size(); index++)
	{
		if (min_z_2 > pos_n_2[index][2]) min_z_2 = pos_n_2[index][2];
	}
	// check that min_z_2, minimum z of plate 1 is larger than max_z_1, maximum z of plate 2
	EXPECT_GT(min_z_2, max_z_1);
}

TEST(StructuralSimulation, MockStent)
{
	Real scale_stl = 0.001;
	Real end_time = 0.08;
	/* MATERIAL PARAMETERS */
	Real rho_0 = 1000.0;
	Real poisson = 0.35;
	Real physical_viscosity = 200;
	/** STL IMPORT PARAMETERS */
	string relative_input_path = "./input/"; //path definition for linux
	vector<string> imported_stl_list = { "mock_stent.stl", "plate_stent.stl", "vessel_cylinder.stl" };
	vector<Vec3d> translation_list = { Vec3d(0, 5, 0), Vec3d(0, 28, 0), Vec3d(0) };
	vector<Real> resolution_list = { 3, 2, 1 };
	LinearElasticSolid material = LinearElasticSolid(rho_0, 1e4, poisson);
	LinearElasticSolid material_stiff = LinearElasticSolid(rho_0, 1e6, poisson);
	vector<LinearElasticSolid> material_model_list = { material_stiff, material, material };
	StructuralSimulationInput input
	{
		relative_input_path,
		imported_stl_list,
		scale_stl,
		translation_list,
		resolution_list,
		material_model_list,
		physical_viscosity,
		{}
	};
	pair<array<int, 2>, array<Real, 2>> contact_pair_1( {0, 1}, {0.0, end_time * 0.4} );
	pair<array<int, 2>, array<Real, 2>> contact_pair_2( {0, 2}, {end_time * 0.5, end_time} );
	input.time_dep_contacting_body_pairs_list_ = { contact_pair_1, contact_pair_2 };

	input.non_zero_gravity_ = { GravityPair(0, Vec3d(0, 100, 0)) };
	Vecd translation_vector = Vec3d(0, -7.5, 0) * scale_stl;
	input.translation_solid_body_tuple_ = { TranslateSolidBodyTuple(1, 0.0, end_time * 0.5, translation_vector) };

	// spring damper constraint
	Real spring_stiffness = 100;
	Real damping_ratio = 0.02;
	input.spring_damper_tuple_ = { SpringDamperTuple(2, Vec3d(spring_stiffness), damping_ratio) };

	/** SIMULATION MODEL */
	TestStructuralSimulation sim(input);

	EXPECT_EQ(sim.Get_contacting_body_pairs_list_().size(), 0);
	EXPECT_EQ(sim.Get_time_dep_contacting_body_pairs_list_().size(), 2);
	EXPECT_EQ(sim.Get_contact_list_().size(), 4);
	EXPECT_EQ(sim.Get_contact_density_list_().size(), 4);
	EXPECT_EQ(sim.Get_contact_force_list_().size(), 4);
	
	/** START SIMULATION */
	sim.TestRunSimulation(end_time);

	// check the plate position, that the translation is correct
	StdLargeVec<Vecd>& pos_0 = sim.Get_solid_body_list_()[1].get()->GetElasticSolidParticles()->pos_0_;
	StdLargeVec<Vecd>& pos_n = sim.Get_solid_body_list_()[1].get()->GetElasticSolidParticles()->pos_n_;
	for (size_t index = 0; index < pos_0.size(); index++)
	{
		Vec3d end_pos = pos_0[index] + translation_vector;
		EXPECT_NEAR(pos_n[index][2], end_pos[2], end_pos.norm() * 1e-2);
	}

	// check that the stent is inside the vessel:
	// 1. the max y pos is larger for the vessel than the stent
	// 1. the min y pos is smaller for the vessel than the stent
	StdLargeVec<Vecd>& pos_n_stent = sim.Get_solid_body_list_()[0].get()->GetElasticSolidParticles()->pos_n_; // stent particles
	StdLargeVec<Vecd>& pos_n_vessel = sim.Get_solid_body_list_()[2].get()->GetElasticSolidParticles()->pos_n_; // vessel particles

	Real max_y_stent = -1000.0;
	Real min_y_stent = 1000.0;
	for (size_t index = 0; index < pos_n_stent.size(); index++)
	{
		if (max_y_stent < pos_n_stent[index][1]) max_y_stent = pos_n_stent[index][1];
		if (min_y_stent > pos_n_stent[index][1]) min_y_stent = pos_n_stent[index][1];
	}
	Real max_y_vessel = -1000.0;
	Real min_y_vessel = 1000.0;
	for (size_t index = 0; index < pos_n_vessel.size(); index++)
	{
		if (max_y_vessel < pos_n_vessel[index][1]) max_y_vessel = pos_n_vessel[index][1];
		if (min_y_vessel > pos_n_vessel[index][1]) min_y_vessel = pos_n_vessel[index][1];
	}
	EXPECT_GT(max_y_vessel, max_y_stent);
	EXPECT_GT(min_y_stent, min_y_vessel);
}

int main(int argc, char* argv[])
{	
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
