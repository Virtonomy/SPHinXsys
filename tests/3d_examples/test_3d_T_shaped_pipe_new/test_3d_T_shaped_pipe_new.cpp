/**
 * @file 	T_shaped_pipe.cpp
 * @brief 	This is the benchmark test of multi-inlet and multi-outlet.
 * @details We consider a flow with one inlet and two outlets in a T-shaped pipe in 2D.
 * @author 	Bence Rochlitz
 */

#include "sphinxsys.h"

using namespace SPH;

Real scale = 0.001;
Real End_Time = 1.0;			/**< End time. */

Real resolution_ref = 7.5 * scale; // particle spacing
Real wall_resolution = 3.0 * scale;
Real emitter_length = resolution_ref * 4;		  /**< Reference size of the emitter. */
Real buffer_length = resolution_ref * 20; /**< Reference size of the emitter buffer to impose inflow condition. */

Real vessel_radius = 100.0 * 0.5 * scale;
Real fluid_radius = 80.0 * 0.5 * scale;
Real vessel_length_half = 200 * scale;

Real vertical_cylinder_length_half = vessel_length_half * 0.5;
Real cannula_outer_radius = 100.0 * 0.5 * scale;
Real cannula_inner_radius = 80.0 * 0.5 * scale;

Real cannula_angle = 0.0;

int resolution_SimTK = 10;

// for material properties of the fluid
Real rho0_f = 1.0;
Real U_f = 1.0;
Real L_f = 0.5; // reference length
Real c_f = 10.0 * U_f;
Real Re = 100.0;
Real mu_f = rho0_f * U_f * L_f / Re;

// material wall
Real youngs_modulus = 1e5;
Real rho_wall = 1e3;
Real possion = 0.35;

/** Domain bounds of the system. */
BoundingBox system_domain_bounds(
	Vecd(-vessel_radius, -vessel_length_half, -vessel_radius),
	Vecd(vessel_radius, vessel_length_half, vertical_cylinder_length_half * 2) // max z coordinate with 0° rotation
);

//	import the fluid body
class WaterBlock : public FluidBody
{
public:
	WaterBlock(SPHSystem &system, const std::string &body_name)
		: FluidBody(system, body_name, makeShared<SPHAdaptation>(1.3, 1.0, 0.75, resolution_ref/wall_resolution))
	{
		// it's the nagative of the WallBoundary 
		ComplexShape _temp_shape;
		// vessel inner fluid
		_temp_shape.add<TriangleMeshShapeCylinder>(SimTK::UnitVec3(0, 1.0, 0), fluid_radius, vessel_length_half, resolution_SimTK, Vec3d(0));
		// angled cannula part
		_temp_shape.add<TriangleMeshShapeCylinder>(SimTK::UnitVec3(0, 0, 1.0), cannula_inner_radius, vertical_cylinder_length_half, resolution_SimTK, Vec3d(0, 0, vertical_cylinder_length_half));

		// actual LevelSetShape
		body_shape_.add<LevelSetShape>(this, _temp_shape, true, false);
	}
};

//	import the static solid wall boundary
class WallBoundary : public SolidBody
{
public:
	WallBoundary(SPHSystem &system, const std::string &body_name)
		: SolidBody(system, body_name, makeShared<SPHAdaptation>(1.15, resolution_ref/wall_resolution, 0.75, 1.0))
	{
		ComplexShape _temp_shape;
		// vessel outer
		_temp_shape.add<TriangleMeshShapeCylinder>(SimTK::UnitVec3(0, 1.0, 0), vessel_radius, vessel_length_half, resolution_SimTK, Vec3d(0));
		// angled cannula part
		_temp_shape.add<TriangleMeshShapeCylinder>(SimTK::UnitVec3(0, 0, 1.0), cannula_outer_radius, vertical_cylinder_length_half, resolution_SimTK, Vec3d(0, 0, vertical_cylinder_length_half));
		// vessel inner fluid
		_temp_shape.substract<TriangleMeshShapeCylinder>(SimTK::UnitVec3(0, 1.0, 0), fluid_radius, vessel_length_half, resolution_SimTK, Vec3d(0));
		// angled cannula part
		_temp_shape.substract<TriangleMeshShapeCylinder>(SimTK::UnitVec3(0, 0, 1.0), cannula_inner_radius, vertical_cylinder_length_half, resolution_SimTK, Vec3d(0, 0, vertical_cylinder_length_half));

		// actual LevelSetShape
		body_shape_.add<LevelSetShape>(this, _temp_shape, true, false);
	}
};
//----------------------------------------------------------------------
//	Define emitter buffer inflow boundary condition
//----------------------------------------------------------------------
class EmitterBufferInflowCondition : public fluid_dynamics::InflowBoundaryCondition
{
	Real u_ave_, u_ref_, t_ref_;

public:
	EmitterBufferInflowCondition(FluidBody &body, BodyPartByCell &body_part)
		: InflowBoundaryCondition(body, body_part),
		  u_ave_(0), u_ref_(U_f), t_ref_(4.0) {}

	Vecd getTargetVelocity(Vecd &position, Vecd &velocity) override
	{
		Real u = velocity[0];
		Real v = velocity[1];
		Real w = velocity[1];

		if (position[1] > vertical_cylinder_length_half)
		{
			u = 0.0;
			v = -u_ave_;
			w = 0.0;
		}
		return Vecd(u, v, w);
	}

	void setupDynamics(Real dt = 0.0) override
	{
		Real run_time = GlobalStaticVariables::physical_time_;
		u_ave_ = run_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / t_ref_)) : u_ref_;
	}
};
//-----------------------------------------------------------------------------------------------------------
//	Main program starts here.
//-----------------------------------------------------------------------------------------------------------
int main(int ac, char *av[])
{
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem with global controls.
	//----------------------------------------------------------------------
	SPHSystem system(system_domain_bounds, resolution_ref);
	/** Tag for computation from restart files. 0: not from restart files. */
	system.restart_step_ = 0;
	//handle command line arguments
	system.handleCommandlineOptions(ac, av);
	In_Output in_output(system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.cd
	//----------------------------------------------------------------------
	WaterBlock water_block(system, "WaterBody");
	// water_block.setBodyDomainBounds(fluid_body_domain_bounds);
	FluidParticles fluid_particles(water_block, makeShared<WeaklyCompressibleFluid>(rho0_f, c_f, mu_f));

	WallBoundary wall_boundary(system, "Wall");
	SolidParticles wall_particles(wall_boundary);
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	BodyRelationInner water_block_inner(water_block);
	ComplexBodyRelation water_block_complex_relation(water_block_inner, {&wall_boundary});
	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simulation.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	/** Initialize particle acceleration. */
	TimeStepInitialization initialize_a_fluid_step(water_block);
	/** Emmiter. */
	TriangleMeshShapeCylinder emitter_shape(SimTK::UnitVec3(0, 1.0, 0), fluid_radius, emitter_length * 0.5, resolution_SimTK, Vec3d(0, 0, vessel_length_half*2 - emitter_length * 0.5));
	BodyRegionByParticle emitter(water_block, "Emitter", emitter_shape);
	fluid_dynamics::EmitterInflowInjecting emitter_inflow_injecting(water_block, emitter, 300, 1, false);
	/** Emitter condition. */
	TriangleMeshShapeCylinder emitter_buffer_shape(SimTK::UnitVec3(0, 1.0, 0), fluid_radius, buffer_length * 0.5, resolution_SimTK, Vec3d(0, 0, vessel_length_half*2 - buffer_length * 0.5));
	BodyRegionByCell emitter_buffer(water_block, "EmitterBuffer", emitter_buffer_shape);
	EmitterBufferInflowCondition emitter_buffer_inflow_condition(water_block, emitter_buffer);
	/** time-space method to detect surface particles. */
	fluid_dynamics::SpatialTemporalFreeSurfaceIdentificationComplex
		inlet_outlet_surface_particle_indicator(water_block_complex_relation);
	/** Evaluation of density by freestream approach. */
	fluid_dynamics::DensitySummationFreeStreamComplex update_density_by_summation(water_block_complex_relation);
	/** We can output a method-specific particle data for debug */
	fluid_particles.addAVariableToWrite<Real>("Pressure");
	fluid_particles.addAVariableToWrite<int>("SurfaceIndicator");
	/** Time step size without considering sound wave speed. */
	fluid_dynamics::AdvectionTimeStepSize get_fluid_advection_time_step_size(water_block, U_f);
	/** Time step size with considering sound wave speed. */
	fluid_dynamics::AcousticTimeStepSize get_fluid_time_step_size(water_block);
	/** Pressure relaxation. */
	fluid_dynamics::PressureRelaxationWithWall pressure_relaxation(water_block_complex_relation);
	/** Density relaxation. */
	fluid_dynamics::DensityRelaxationRiemannWithWall density_relaxation(water_block_complex_relation);
	/** Computing viscous acceleration. */
	fluid_dynamics::ViscousAccelerationWithWall viscous_acceleration(water_block_complex_relation);
	/** Impose transport velocity. */
	fluid_dynamics::TransportVelocityCorrectionComplex transport_velocity_correction(water_block_complex_relation);
	/** recycle real fluid particle to buffer particles at outlet. */
	OpenBoundaryConditionInAxisDirection tansfer_to_buffer_particles_lower_bound(water_block, yAxis, negativeDirection);
	OpenBoundaryConditionInAxisDirection tansfer_to_buffer_particles_upper_bound(water_block, yAxis, positiveDirection);
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations and observations of the simulation.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp write_body_states(in_output, system.real_bodies_);
	RestartIO restart_io(in_output, system.real_bodies_);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	wall_particles.initializeNormalDirectionFromBodyShape();
	inlet_outlet_surface_particle_indicator.parallel_exec();
	//----------------------------------------------------------------------
	//	Load restart file if necessary.
	//----------------------------------------------------------------------
	/** If the starting time is not zero, please setup the restart time step ro read in restart states. */
	if (system.restart_step_ != 0)
	{
		GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(system.restart_step_);
		water_block.updateCellLinkedList();
		water_block_complex_relation.updateConfiguration();
	}
	//----------------------------------------------------------------------
	//	Setup computing and initial conditions.
	//----------------------------------------------------------------------
	size_t number_of_iterations = system.restart_step_;
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval * 10;

	Real D_Time = End_Time / 100.0; /**< Time stamps for output of body states. */
	Real dt = 0.0;					/**< Default acoustic time step sizes. */
	//----------------------------------------------------------------------
	//	Statistics for CPU time
	//----------------------------------------------------------------------
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	//----------------------------------------------------------------------
	//	First output before the main loop.
	//----------------------------------------------------------------------
	write_body_states.writeToFile();
	//----------------------------------------------------------------------------------------------------
	//	Main loop starts here.
	//----------------------------------------------------------------------------------------------------
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < D_Time)
		{
			initialize_a_fluid_step.parallel_exec();
			Real Dt = get_fluid_advection_time_step_size.parallel_exec();
			update_density_by_summation.parallel_exec();
			viscous_acceleration.parallel_exec();
			transport_velocity_correction.parallel_exec(Dt);

			/** Dynamics including pressure relaxation. */
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				dt = SMIN(get_fluid_time_step_size.parallel_exec(), Dt - relaxation_time);
				pressure_relaxation.parallel_exec(dt);
				emitter_buffer_inflow_condition.parallel_exec();
				density_relaxation.parallel_exec(dt);

				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}

			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
						  << GlobalStaticVariables::physical_time_
						  << "	Dt = " << Dt << "	dt = " << dt << "\n";

				if (number_of_iterations % restart_output_interval == 0 && number_of_iterations != system.restart_step_)
					restart_io.writeToFile(Real(number_of_iterations));
			}
			number_of_iterations++;

			/** inflow injecting*/
			emitter_inflow_injecting.exec();
			tansfer_to_buffer_particles_lower_bound.particle_type_transfer.parallel_exec();
			tansfer_to_buffer_particles_upper_bound.particle_type_transfer.parallel_exec();

			/** Update cell linked list and configuration. */
			water_block.updateCellLinkedList();
			water_block_complex_relation.updateConfiguration();
			inlet_outlet_surface_particle_indicator.parallel_exec();
		}

		tick_count t2 = tick_count::now();
		write_body_states.writeToFile();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds()
			  << " seconds." << std::endl;

	return 0;
}
