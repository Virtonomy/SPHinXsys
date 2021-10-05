#include "structural_simulation_class.h"
#include <gtest/gtest.h>

Real tolerance = 1e-6;

class TestStructuralSimulation : StructuralSimulation
{
public:
	TestStructuralSimulation(StructuralSimulationInput& input) : StructuralSimulation(input){};
	void TestRunSimulation(Real end_time){ runSimulation(end_time); };
	// input members
	string get_relative_input_path_(){ return relative_input_path_; };
	vector<string> get_imported_stl_list_(){ return imported_stl_list_; };
	Real get_scale_stl_(){ return scale_stl_; };
	vector<Vec3d> get_translation_list_(){ return translation_list_; };
	Real get_system_resolution_(){ return system_resolution_; };
	vector<Real> get_resolution_list_(){ return resolution_list_; };
	vector<LinearElasticSolid> get_material_model_list_(){ return material_model_list_; };
	vector<Real> get_physical_viscosity_(){ return physical_viscosity_; };
	// internal members
	SPHSystem get_system_(){ return system_; };
	Real get_scale_system_boundaries_(){ return scale_system_boundaries_; };
	In_Output get_in_output_(){ return in_output_; };
	// other
	vector<TriangleMeshShape> get_body_mesh_list_(){ return body_mesh_list_; };
	vector<shared_ptr<SolidBodyForSimulation>> get_solid_body_list_(){ return solid_body_list_; };
	StdVec<IndexVector> get_contacting_body_pairs_list_(){ return contacting_body_pairs_list_; };
	vector<pair<array<int, 2>, array<Real, 2>>> get_time_dep_contacting_body_pairs_list_(){ return time_dep_contacting_body_pairs_list_; };
	vector<shared_ptr<SolidBodyRelationContact>> get_contact_list_(){ return contact_list_; };
	vector<shared_ptr<solid_dynamics::ContactDensitySummation>> get_contact_density_list_(){ return contact_density_list_; };
	vector<shared_ptr<solid_dynamics::ContactForce>> get_contact_force_list_(){ return contact_force_list_; };
	// for initializeATimeStep
	vector<shared_ptr<TimeStepInitialization>> get_initialize_gravity_(){ return initialize_gravity_; };
	vector<GravityPair> get_non_zero_gravity_(){ return non_zero_gravity_; };
	// for AccelerationForBodyPartInBoundingBox
	vector<shared_ptr<solid_dynamics::AccelerationForBodyPartInBoundingBox>> get_acceleration_bounding_box_(){ return acceleration_bounding_box_; };
	vector<AccelTuple> get_acceleration_bounding_box_tuple_(){ return acceleration_bounding_box_tuple_; };
	// for SurfacePressureFromSource
	vector<shared_ptr<solid_dynamics::SurfacePressureFromSource>> get_surface_pressure_(){ return surface_pressure_; };
	vector<PressureTuple> get_surface_pressure_tuple_(){ return surface_pressure_tuple_; };
	// for SpringDamperConstraintParticleWise
	vector<shared_ptr<solid_dynamics::SpringDamperConstraintParticleWise>> get_spring_damper_constraint_(){ return spring_damper_constraint_; };
	vector<SpringDamperTuple> get_spring_damper_tuple_(){ return spring_damper_tuple_; };
	// for SurfaceSpring
	vector<shared_ptr<solid_dynamics::SpringNormalOnSurfaceParticles>> get_surface_spring_(){ return surface_spring_; };
	vector<SurfaceSpringTuple> get_surface_spring_tuple_(){ return surface_spring_tuple_; };
	// for ConstrainSolidBody
	vector<shared_ptr<solid_dynamics::ConstrainSolidBodyRegion>> get_fixed_constraint_(){ return fixed_constraint_body_; };
	vector<int> get_body_indeces_fixed_constraint_(){ return body_indeces_fixed_constraint_; };
	// for PositionSolidBody
	vector<shared_ptr<solid_dynamics::PositionSolidBody>> get_position_solid_body_(){ return position_solid_body_; };
	vector<PositionSolidBodyTuple> get_position_solid_body_tuple_(){ return position_solid_body_tuple_; };
	// for PositionScaleSolidBody
	vector<shared_ptr<solid_dynamics::PositionScaleSolidBody>> get_position_scale_solid_body_(){ return position_scale_solid_body_; };
	vector<PositionScaleSolidBodyTuple> get_position_scale_solid_body_tuple_(){ return position_scale_solid_body_tuple_; };
    // for TranslateSolidBody
	vector<shared_ptr<solid_dynamics::TranslateSolidBody>> get_translation_solid_body_(){ return translation_solid_body_; };
	vector<TranslateSolidBodyTuple> get_translation_solid_body_tuple_(){ return translation_solid_body_tuple_; };

	// get data
	vector<Real> getVonMisesStressMax(){ return von_mises_stress_max_; };

	StdLargeVec<StdLargeVec<Real>> getVonMisesStressOfParticlesOverTime(){
		// get the time series of von Mises stress of all particles
		return von_mises_stress_particles_;
	};

	vector<Real> getVonMisesStrainMax(){ return von_mises_strain_max_; };

	StdLargeVec<StdLargeVec<Real>> getVonMisesStrainOfParticlesOverTime(){
		// get the time series of von Mises strain of all particles
		return von_mises_strain_particles_;
	};

	StdLargeVec<Real> getVonMisesStressAverageOfParticlesOverTime(StdLargeVec<int> particle_IDs){
		// get the time series of average stress of the input particles
		// input: particle IDs to be averaged
		// output: time series of average stress
		
		StdLargeVec<StdLargeVec<Real>> von_mises_stress_particles = getVonMisesStressOfParticlesOverTime();

		StdLargeVec<Real> avg_stress_over_time = {};
		// loop over time steps
		// von_mises_stress_particles.size() = all time steps
		for(size_t i = 0; i < von_mises_stress_particles.size(); i++)
		{
			Real current_avg = 0.0;
			// loop over particles
			// get the cumulative stress of all particles in a time step
			for(size_t j = 0; j < particle_IDs.size(); j++)
			{	
				// get the current ID of a particle
				size_t current_ID = particle_IDs[j];
				// time step "i" and particle ID "current_ID"
				current_avg += von_mises_stress_particles[i][current_ID];
			}
			// average out the particle stress
			current_avg /= particle_IDs.size();
			// record the average into the avg_stress_over_time
			avg_stress_over_time.push_back(current_avg);
		}

		return avg_stress_over_time;
	};

	StdLargeVec<Real> getVonMisesStrainAverageOfParticlesOverTime(StdLargeVec<int> particle_IDs){
		// get the time series of average strain of the input particles
		// input: particle IDs to be averaged
		// output: time series of average stress
		
		StdLargeVec<StdLargeVec<Real>> von_mises_strain_particles = getVonMisesStrainOfParticlesOverTime();

		StdLargeVec<Real> avg_strain_over_time = {};
		// loop over time steps
		// von_mises_strain_particles.size() = all time steps
		for(size_t i = 0; i < von_mises_strain_particles.size(); i++)
		{
			Real current_avg = 0.0;
			// loop over particles
			// get the cumulative stress of all particles in a time step
			for(size_t j = 0; j < particle_IDs.size(); j++)
			{	
				// get the current ID of a particle
				size_t current_ID = particle_IDs[j];
				// time step "i" and particle ID "current_ID"
				current_avg += von_mises_strain_particles[i][current_ID];
			}
			// average out the particle stress
			current_avg /= particle_IDs.size();
			// record the average into the avg_stress_over_time
			avg_strain_over_time.push_back(current_avg);
		}

		return avg_strain_over_time;
	};

	Real getVonMisesStressCloseToPoint(int body_index, Vec3d point)
	{
		// find particle id closest to point
		StdLargeVec<Vecd>& pos_0 = solid_body_list_[body_index].get()->getElasticSolidParticles()->pos_0_;
		size_t particle_id = 0;
		Real min_distance = (pos_0[0] - point).norm();
		for (size_t i = 0; i < pos_0.size(); i++)
		{
			Real distance = (pos_0[i] - point).norm();
			if (min_distance > distance)
			{
				min_distance = distance;
				particle_id = i;
			}
		}
		Real stress = solid_body_list_[body_index].get()->getElasticSolidParticles()->von_Mises_stress(particle_id);
		return stress;
	};
	Real getMaxDisplacement(int body_index)
	{
		StdLargeVec<Vecd>& pos_0 = solid_body_list_[body_index].get()->getElasticSolidParticles()->pos_0_;
		StdLargeVec<Vecd>& pos_n = solid_body_list_[body_index].get()->getElasticSolidParticles()->pos_n_;
		Real displ_max = 0;
		for (size_t i = 0; i < pos_0.size(); i++)
		{
			Real displ = (pos_n[i] - pos_0[i]).norm();
			if (displ > displ_max) displ_max = displ;
		}
		return displ_max;
	}
};