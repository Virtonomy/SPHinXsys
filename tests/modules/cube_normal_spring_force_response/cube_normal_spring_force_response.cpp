#include <gtest/gtest.h>
#include "test_structural_simulation_class.h"

void NormalSpringForceTest(Real end_time, Real scale_stl, bool particle_relaxation, Real particles_in_thickness, Real gravity_force, 
							Real spring_surface_pressure, Real spring_damper_ratio, Vecd source_point, Real tolerance)
{
	/** INPUT PARAMETERS */
	Real rho_0 = 1000.0;
	Real poisson = 0.45;
	Real Youngs_modulus = 5e4;
	Real physical_viscosity = 200;

	Real resolution_cube = 100 / particles_in_thickness;

	/** STL IMPORT PARAMETERS */
	std::string relative_input_path = "./input/"; //path definition for linux
	std::vector<std::string> imported_stl_list = { "cube.stl" };
	std::vector<Vec3d> translation_list = { Vec3d(0) };
	std::vector<Real> resolution_list = { resolution_cube};
	LinearElasticSolid material = LinearElasticSolid(rho_0, Youngs_modulus, poisson);
	std::vector<LinearElasticSolid> material_model_list = { material };

	/** INPUT DECLERATION */
	StructuralSimulationInput input
	{
		relative_input_path,
		imported_stl_list,
		scale_stl,
		translation_list,
		resolution_list,
		material_model_list,
		{physical_viscosity},
		{}
	};
	input.particle_relaxation_list_ = { particle_relaxation };
	input.non_zero_gravity_ = vector<GravityPair>{ GravityPair(0, Vec3d(0.0, 0.0, gravity_force)) };
	TriangleMeshShape cube("./input/cube.stl", Vec3d(0), scale_stl);
	input.surface_spring_tuple_ = StdVec<SurfaceSpringTuple>{ SurfaceSpringTuple(0, &cube, false, source_point, spring_surface_pressure, spring_damper_ratio ) };

	//=================================================================================================//
	TestStructuralSimulation sim(input);

	// check the particle relaxation, if the initial position is set up correctly	
	StdLargeVec<Vecd>& pos_0 = sim.get_solid_body_list_()[0].get()->getElasticSolidParticles()->pos_0_;
	StdLargeVec<Vecd>& pos_n = sim.get_solid_body_list_()[0].get()->getElasticSolidParticles()->pos_n_;
	for (size_t index = 0; index < pos_0.size(); index++)
	{	
		EXPECT_EQ(pos_0[index], pos_n[index]);
	}

	sim.TestRunSimulation(end_time);
	//=================================================================================================//

	// check which particles have spring force
	// only one layer of particles on one side should have
	int particles_ref = (particles_in_thickness * particles_in_thickness) - 4;
	StdLargeVec<bool>& apply_spring_force_to_particle_ = sim.get_surface_spring_()[0].get()->GetApplySpringForceToParticle();
	int particles_with_spring_force = 0;
	for (size_t i = 0; i < apply_spring_force_to_particle_.size(); i++)
	{
		if (apply_spring_force_to_particle_[i]) particles_with_spring_force++;
	}
	// if (!particle_relaxation)
	// {
	// 	EXPECT_EQ(particles_ref, particles_with_spring_force);
	// }

	//check if the resulting displacement is correct
	Real end_displ = 0.05;
	//check mass
	StdLargeVec<Real>& mass = sim.get_solid_body_list_()[0].get()->getElasticSolidParticles()->mass_;
	StdLargeVec<Vecd>& normal = sim.get_solid_body_list_()[0].get()->getElasticSolidParticles()->n_0_;
	for (size_t index = 0; index < pos_0.size(); index++)
	{
		if(pos_0[index][2] < 0.01)
		{
			// Real end_pos = pos_0[index][2] + end_displ;
			// EXPECT_NEAR(pos_n[index][2], end_pos, end_displ * tolerance); // 1% tolerance without particle relaxation
			
			EXPECT_DOUBLE_EQ(mass[index], (1.0 / std::pow(particles_in_thickness, 3.0)));
			// if (pos_0[index][0] < 0.01 && pos_0[index][0] > 0.09 &&
			// pos_0[index][1] < 0.01 && pos_0[index][1] > 0.09) //printing out bottom particles without edge points
			if (pos_0[index][0] < 0.01) //printing out some edge points
			{
				std::cout <<"bottom point: " << index << std::endl;
				std::cout <<"normal vector: " << normal[index] << std::endl;
			}
		}
		
		//check for sideway movement
		Real final_x_position = pos_0[index][0];
		Real final_y_position = pos_0[index][1];
		EXPECT_NEAR(final_x_position, pos_n[index][0], tolerance);
		EXPECT_NEAR(final_y_position, pos_n[index][1], tolerance);
	}

	// TODO: fix sideway normal vectors for particle relaxation and improve test with particle relaxation

}

// TEST(NormalSpringForceResponse, NoParticleRelaxation)
// {
// 	Real end_time = 3.0;
// 	Real scale_stl = 0.001;
// 	bool particle_relaxation = false;
// 	Real particles_in_thickness = 10;

// 	// SI setup - designed
// 	Real gravity_force = 8.0;
// 	Real spring_surface_pressure = 2e4;
// 	Real spring_damper_ratio = 0.05;
// 	Vecd source_point = Vec3d(50, 50, -50) * scale_stl;

// 	Real tolerance = 1e-2; // 1% tolerance without particle relaxation
// 	NormalSpringForceTest(end_time, scale_stl, particle_relaxation, particles_in_thickness, gravity_force, spring_surface_pressure, spring_damper_ratio, source_point, tolerance);
// }

// TEST(NormalSpringForceResponse, WithParticleRelaxation)
// {
// 	Real end_time = 3.0;
// 	Real scale_stl = 0.001;
// 	bool particle_relaxation = true;
// 	Real particles_in_thickness = 10;

// 	// SI setup - designed
// 	Real gravity_force = 8.0;
// 	Real spring_surface_pressure = 2e4;
// 	Real spring_damper_ratio = 0.05;
// 	Vecd source_point = Vec3d(50, 50, -50) * scale_stl;

// 	Real tolerance = 1e-1; // 10% tolerance with particle relaxation
// 	NormalSpringForceTest(end_time, scale_stl, particle_relaxation, particles_in_thickness, gravity_force, spring_surface_pressure, spring_damper_ratio, source_point, tolerance);
// }

TEST(NormalSpringForceResponse, DifferentSourcePoint)
{
	Real end_time = 3.0;
	Real scale_stl = 0.001;
	bool particle_relaxation = false;
	Real particles_in_thickness = 10;

	// SI setup - designed
	Real gravity_force = 8.0;
	Real spring_surface_pressure = 2e4;
	Real spring_damper_ratio = 0.05;
	Vecd source_point = Vec3d(50, 50, -25) * scale_stl;

	Real tolerance = 1e-1; // 1% tolerance without particle relaxation
	NormalSpringForceTest(end_time, scale_stl, particle_relaxation, particles_in_thickness, gravity_force, spring_surface_pressure, spring_damper_ratio, source_point, tolerance);
}

int main(int argc, char* argv[])
{	
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
