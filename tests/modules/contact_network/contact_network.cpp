#include <gtest/gtest.h>
#include "test_structural_simulation_class.h"

// these are tests to make sure all the contacts are properly executed for all the bodies
// there was a bug that some contacts were not executed, if the contact pairs were defined
// change: now the contacts are defined body-wise, meaning all target bodies are included for the contact body at once, instead of pair-wise

TEST(ContactNetwork, 2BodiesWith2Bodies)
{	
	// two bodies collide with other two bodies
	// both have a contact with the other two
	
	/** INPUT PARAMETERS */
	Real scale_stl = 0.001;

	Real rho_0 = 1000.0;
	Real poisson = 0.35;
	Real Youngs_modulus = 5e4;
	Real physical_viscosity = Youngs_modulus/100;

	Real end_time = 0.5;

	/** STL IMPORT PARAMETERS */
	string relative_input_path = "./input/"; //path definition for linux
	string cube = "cube.stl";
	string cube2 = "cube2.stl";
	string ball = "ball.stl";
	string ball2 = "ball2.stl";

	vector<string> imported_stl_list = {
		cube, cube2, ball, ball2
	};
	vector<Vec3d> translation_list(4, Vec3d(0));
	Real res = 6.0;
	vector<Real> resolution_list(4, res);

	LinearElasticSolid material_tah = LinearElasticSolid(rho_0, Youngs_modulus, poisson);
	NeoHookeanSolid material_tissue = NeoHookeanSolid(rho_0, Youngs_modulus, poisson);
	vector<LinearElasticSolid> material_model_list = {
		material_tah, material_tah, material_tissue, material_tissue
	};

	StdVec<IndexVector> contacting_bodies_list = {
		{2, 3},
		{2, 3},
		{0, 1},
		{0, 1}
	};

	/** INPUT DECLERATION */
	StructuralSimulationInput input
	{
		relative_input_path,
		imported_stl_list,
		scale_stl,
		translation_list,
		resolution_list,
		material_model_list,
		StdVec<Real>(4, physical_viscosity),
		contacting_bodies_list
	};

	/** CALCULATION OF TAH CENTER POINT*/
	string name_left = relative_input_path;
	string name_right = relative_input_path;
	TriangleMeshShape cube_left_mesh(name_left.append(cube), translation_list[0] * scale_stl, scale_stl);
	TriangleMeshShape cube_right_mesh(name_right.append(cube2), translation_list[1] * scale_stl, scale_stl);
	BoundingBox bb_left = cube_left_mesh.findBounds();
	BoundingBox bb_right = cube_right_mesh.findBounds();
	Vec3d translate(0.0, 30.0, 0.0);
	Vec3d pos_final_left = (bb_left.first + bb_left.second) * 0.5 + translate * scale_stl;
	Vec3d pos_final_right = (bb_right.first + bb_right.second) * 0.5 + translate * scale_stl;

	input.position_solid_body_tuple_ = {
		PositionSolidBodyTuple(0, 0.0, end_time, pos_final_left),
		PositionSolidBodyTuple(1, 0.0, end_time, pos_final_right)
	};
	
	//spring damper constraint
	Real spring_stiffness = 100;
	Real damping_ratio = 0.0;
	input.spring_damper_tuple_ = {
		SpringDamperTuple(2, Vec3d(spring_stiffness), damping_ratio),
		SpringDamperTuple(3, Vec3d(spring_stiffness), damping_ratio)
	};
	
	//=================================================================================================//
	TestStructuralSimulation sim (input);
	sim.TestRunSimulation(end_time);
	//=================================================================================================//

	StdLargeVec<Vecd>& ball_pos_n = sim.get_solid_body_list_()[2].get()->getElasticSolidParticles()->pos_n_;
	StdLargeVec<Vecd>& ball2_pos_n = sim.get_solid_body_list_()[3].get()->getElasticSolidParticles()->pos_n_;

	for (Vec3d pos: ball_pos_n)
	{
		// all the particles of the ball should outside of the cube
		// 0.05 is the y coordinate of the closest side of the cube at the end position, residual overlap is possible
		EXPECT_GT(pos[1], 0.048);
	}
	for (Vec3d pos: ball2_pos_n)
	{
		// all the particles of the ball should outside of the cube
		// 0.05 is the y coordinate of the closest side of the cube at the end position, residual overlap is possible
		EXPECT_GT(pos[1], 0.048);
	}
}

TEST(ContactNetwork, 4BodiesAllWithAll)
{	
	// two bodies collide with other two bodies
	// both have a contact with the other two
	
	/** INPUT PARAMETERS */
	Real scale_stl = 0.001;

	Real rho_0 = 1000.0;
	Real poisson = 0.35;
	Real Youngs_modulus = 5e4;
	Real physical_viscosity = Youngs_modulus/100;

	Real end_time = 0.5;

	/** STL IMPORT PARAMETERS */
	string relative_input_path = "./input/"; //path definition for linux
	string cube = "cube_b.stl"; // tilted, moving brick
	string cube2 = "cube2_b.stl"; // vertical brick
	string ball = "ball_b.stl";
	string ball2 = "ball2_b.stl";

	vector<string> imported_stl_list = {
		cube, cube2, ball, ball2
	};
	vector<Vec3d> translation_list(4, Vec3d(0));
	Real res = 5.0;
	vector<Real> resolution_list(4, res);

	LinearElasticSolid material_tah = LinearElasticSolid(rho_0, Youngs_modulus, poisson);
	NeoHookeanSolid material_tissue = NeoHookeanSolid(rho_0, Youngs_modulus, poisson);
	vector<LinearElasticSolid> material_model_list = {
		material_tah, material_tah, material_tissue, material_tissue
	};

	StdVec<IndexVector> contacting_bodies_list = {
		{1, 2, 3},
		{0, 2, 3},
		{0, 1, 3},
		{0, 1, 2}
	};

	/** INPUT DECLERATION */
	StructuralSimulationInput input
	{
		relative_input_path,
		imported_stl_list,
		scale_stl,
		translation_list,
		resolution_list,
		material_model_list,
		StdVec<Real>(4, physical_viscosity),
		contacting_bodies_list
	};

	/** CALCULATION OF TAH CENTER POINT*/
	string name_left = relative_input_path;
	TriangleMeshShape cube_left_mesh(name_left.append(cube), translation_list[0] * scale_stl, scale_stl);
	BoundingBox bb_left = cube_left_mesh.findBounds();
	Vec3d translate(0.0, -30.0, 0.0);
	Vec3d pos_final_left = (bb_left.first + bb_left.second) * 0.5 + translate * scale_stl;
	input.position_solid_body_tuple_ = {
		PositionSolidBodyTuple(0, 0.0, end_time, pos_final_left)
	};

	//spring damper constraint
	Real spring_stiffness = 100;
	Real damping_ratio = 0.0;
	input.spring_damper_tuple_ = {
		SpringDamperTuple(1, Vec3d(spring_stiffness), damping_ratio),
		SpringDamperTuple(2, Vec3d(spring_stiffness), damping_ratio),
		SpringDamperTuple(3, Vec3d(spring_stiffness), damping_ratio)
	};
	
	//=================================================================================================//
	TestStructuralSimulation sim (input);
	sim.TestRunSimulation(end_time);
	//=================================================================================================//

	StdLargeVec<Vecd>& brick_tilted_pos_n = sim.get_solid_body_list_()[0].get()->getElasticSolidParticles()->pos_n_;
	StdLargeVec<Vecd>& brick_vertical_pos_n = sim.get_solid_body_list_()[1].get()->getElasticSolidParticles()->pos_n_;
	StdLargeVec<Vecd>& ball_pos_n = sim.get_solid_body_list_()[2].get()->getElasticSolidParticles()->pos_n_;
	StdLargeVec<Vecd>& ball2_pos_n = sim.get_solid_body_list_()[3].get()->getElasticSolidParticles()->pos_n_;

	// check that the balls stay on the left side of the tilted brick because of the contact
	// check that all the ball particles are "under" the brick
	// z < -2y + 0.047
	for (Vec3d pos: ball_pos_n) EXPECT_LT(pos[2], -0.5 * pos[1] + 0.047);
	for (Vec3d pos: ball2_pos_n) EXPECT_LT(pos[2], -0.5 * pos[1] + 0.047);

	// chech that the balls are on the right side of the vertical brick
	// find the right most particle y coordinate
	Real y_limit = -1e3;
	for (Vec3d pos: brick_vertical_pos_n)
	{
		if (pos[1] > y_limit) y_limit = pos[1];
	}
	for (Vec3d pos: ball_pos_n) EXPECT_GT(pos[1], y_limit);
	for (Vec3d pos: ball2_pos_n) EXPECT_GT(pos[1], y_limit);
}


int main(int argc, char* argv[])
{	
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
