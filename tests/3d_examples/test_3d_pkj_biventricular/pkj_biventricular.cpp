/**
 * @file 	pkj_lv_electrocontraction.cpp
 * @brief 	This is the case studying the electromechanics on a biventricular heart model in 3D.
 * @author 	Chi Zhang and Xiangyu Hu
 * @version 0.2.1
 * 			Chi Zhang
 * 			Unit :
 *			time t = ms = 12.9 [-]
 * 			length l = mm
 * 			mass m = g
 *			density rho = g * (mm)^(-3)
 *			Pressure pa = g * (mm)^(-1) * (ms)^(-2)
 *			diffusion d = (mm)^(2) * (ms)^(-2)
 *@version 0.3
 *			Here, the coupling with Purkinje network will be condcuted.
 */
/**  SPHinXsys Library. */
#include "pkj_biventricular.h"
/** Namespace cite here. */
using namespace SPH;
/** 
 * The main program. 
 */
int main(int ac, char *av[])
{
	/** 
	 * Build up context -- a SPHSystem. 
	 */
	SPHSystem system(system_domain_bounds, dp_0);
	/** Set the starting time. */
	GlobalStaticVariables::physical_time_ = 0.0;
	/** Tag for run particle relaxation for the initial body fitted distribution. */
	system.run_particle_relaxation_ = true;
	/** Tag for reload initially repaxed particles. */
	system.reload_particles_ = false;
	/** Tag for computation from restart files. 0: not from restart files. */
	system.restart_step_ = 0;
//handle command line arguments
#ifdef BOOST_AVAILABLE
	system.handleCommandlineOptions(ac, av);
#endif
	/** in- and output environment. */
	InOutput in_output(system);
	//----------------------------------------------------------------------
	//	SPH Particle relaxation section
	//----------------------------------------------------------------------
	/** check whether run particle relaxation for body fitted particle distribution. */
	if (system.run_particle_relaxation_)
	{
		SolidBody heart_model(system, makeShared<Heart>("HeartModel"));
		heart_model.defineBodyLevelSetShape()->writeLevelSet(heart_model);
		heart_model.defineParticlesAndMaterial<DiffusionReactionParticles<ElasticSolidParticles>, FiberDirectionDiffusion>();
		heart_model.generateParticles<ParticleGeneratorLattice>();
		/** topology */
		BodyRelationInner heart_model_inner(heart_model);
		/** Random reset the relax solid particle position. */
		RandomizeParticlePosition random_particles(heart_model);
		/** A  Physics relaxation step. */
		relax_dynamics::RelaxationStepInner relaxation_step_inner(heart_model_inner);
		/** Time step for diffusion. */
		GetDiffusionTimeStepSize<SolidBody, ElasticSolidParticles, LocallyOrthotropicMuscle> get_time_step_size(heart_model);
		/** Diffusion process for diffusion body. */
		DiffusionRelaxation diffusion_relaxation(heart_model_inner);
		/** Compute the fiber and sheet after diffusion. */
		ComputeFiberandSheetDirections compute_fiber_sheet(heart_model);
		/** Write the body state to Vtp file. */
		BodyStatesRecordingToVtp write_heart_model_state_to_vtp(in_output, {heart_model});
		/** Write the particle reload files. */
		ReloadParticleIO write_particle_reload_files(in_output, {&heart_model}, {"HeartModel"});
		/** Write material property to xml file. */
		ReloadMaterialParameterIO write_material_property(in_output, heart_model.base_material_, "FiberDirection");
		//----------------------------------------------------------------------
		//	Physics relaxation starts here.
		//----------------------------------------------------------------------
		random_particles.parallel_exec(0.25);
		relaxation_step_inner.surface_bounding_.parallel_exec();
		write_heart_model_state_to_vtp.writeToFile(0.0);
		//----------------------------------------------------------------------
		//From here the time stepping begins.
		//----------------------------------------------------------------------
		int ite = 0;
		int relax_step = 1000;
		int diffusion_step = 100;
		while (ite < relax_step)
		{
			relaxation_step_inner.parallel_exec();
			ite++;
			if (ite % 100 == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite << "\n";
				write_heart_model_state_to_vtp.writeToFile(ite);
			}
		}

		BodySurface surface_part(heart_model);
		/** constraint boundary condition for diffusion. */
		DiffusionBCs impose_diffusion_bc(heart_model, surface_part);
		impose_diffusion_bc.parallel_exec();

		write_heart_model_state_to_vtp.writeToFile(ite);

		Real dt = get_time_step_size.parallel_exec();
		while (ite <= diffusion_step + relax_step)
		{
			diffusion_relaxation.parallel_exec(dt);
			impose_diffusion_bc.parallel_exec();
			if (ite % 10 == 0)
			{
				std::cout << "Diffusion steps N=" << ite - relax_step << "	dt: " << dt << "\n";
				write_heart_model_state_to_vtp.writeToFile(ite);
			}
			ite++;
		}
		compute_fiber_sheet.exec();
		ite++;
		write_heart_model_state_to_vtp.writeToFile(ite);
		compute_fiber_sheet.parallel_exec();
		write_material_property.writeToFile(0);
		write_particle_reload_files.writeToFile(0);

		return 0;
	}
	//----------------------------------------------------------------------
	//	SPH simultion section
	//----------------------------------------------------------------------
	/** create a SPH body, material and particles */
	SolidBody physiology_heart(system, makeShared<Heart>("PhysiologyHeart"));
	AlievPanfilowModel muscle_reaction_model(k_a, c_m, k, a, b, mu_1, mu_2, epsilon);
	physiology_heart.defineParticlesAndMaterial<
		ElectroPhysiologyParticles, LocalMonoFieldElectroPhysiology>(muscle_reaction_model, diffusion_coff, bias_coff, fiber_direction);
	(!system.run_particle_relaxation_ && system.reload_particles_)
		? physiology_heart.generateParticles<ParticleGeneratorReload>(in_output, "HeartModel")
		: physiology_heart.generateParticles<ParticleGeneratorLattice>();

	/** create a SPH body, material and particles */
	SolidBody mechanics_heart(system,  makeShared<Heart>("MechanicalHeart"));
	mechanics_heart.defineParticlesAndMaterial<
		ElasticSolidParticles, ActiveMuscle<LocallyOrthotropicMuscle>>(rho0_s, bulk_modulus, fiber_direction, sheet_direction, a0, b0);
	(!system.run_particle_relaxation_ && system.reload_particles_)
		? mechanics_heart.generateParticles<ParticleGeneratorReload>(in_output, "HeartModel")
		: mechanics_heart.generateParticles<ParticleGeneratorLattice>();

	/** check whether reload material properties. */
	if (!system.run_particle_relaxation_ && system.reload_particles_)
	{
		ReloadMaterialParameterIO read_physiology_heart_fiber(in_output, physiology_heart.base_material_, "FiberDirection");
		ReloadMaterialParameterIO read_mechanics_heart_fiber(in_output, mechanics_heart.base_material_, "FiberDirection");
		read_mechanics_heart_fiber.readFromFile();
		read_physiology_heart_fiber.readFromFile();
	} 

	/** Creat a Purkinje network for fast diffusion, material and particles - Right Ventricle */
	TreeBody pkj_body_RV(system, makeShared<Heart>("Purkinje"));
	AlievPanfilowModel pkj_reaction_model_RV(k_a, c_m, k, a, b, mu_1, mu_2, epsilon);
	pkj_body_RV.defineParticlesAndMaterial<
		ElectroPhysiologyReducedParticles, MonoFieldElectroPhysiology>(
			pkj_reaction_model_RV, diffusion_coff * acceleration_factor, bias_coff, fiber_direction);
	pkj_body_RV.generateParticles<NetworkGeneratorWithExtraCheck>(starting_point_RV, second_point_RV, 50, 1.0, 1.79, 10);
	TreeTerminates pkj_leaves_RV(pkj_body_RV);

	/** Creat a Purkinje network for fast diffusion, material and particles - Left Ventricle */
	TreeBody pkj_body_LV(system, makeShared<Heart>("Purkinje"));
	AlievPanfilowModel pkj_reaction_model_LV(k_a, c_m, k, a, b, mu_1, mu_2, epsilon);
	pkj_body_LV.defineParticlesAndMaterial<
		ElectroPhysiologyReducedParticles, MonoFieldElectroPhysiology>(
			pkj_reaction_model_LV, diffusion_coff * acceleration_factor, bias_coff, fiber_direction);
	pkj_body_LV.generateParticles<NetworkGeneratorWithExtraCheck>(starting_point_LV, second_point_LV, 50, 1.0, 1.79, 10);
	TreeTerminates pkj_leaves_LV(pkj_body_LV);
	//----------------------------------------------------------------------
	//	SPH Observation section
	//----------------------------------------------------------------------
	ObserverBody voltage_observer(system, "VoltageObserver");
	voltage_observer.generateParticles<HeartObserverParticleGenerator>();
	ObserverBody myocardium_observer(system, "MyocardiumObserver");
	myocardium_observer.generateParticles<HeartObserverParticleGenerator>();

	/** topology */
	BodyRelationInner physiology_heart_inner(physiology_heart);
	BodyRelationInner mechanics_heart_inner(mechanics_heart);
	BodyRelationContact physiology_heart_contact(physiology_heart, {&mechanics_heart});
	BodyRelationContact mechanics_heart_contact(mechanics_heart, {&physiology_heart});
	BodyRelationContact voltage_observer_contact(voltage_observer, {&physiology_heart});
	BodyRelationContact myocardium_observer_contact(myocardium_observer, {&mechanics_heart});
	ComplexBodyRelation physiology_heart_complex_RV(physiology_heart, {&pkj_leaves_RV});
	TreeBodyRelationInner pkj_inner_RV(pkj_body_RV);
	ComplexBodyRelation physiology_heart_complex_LV(physiology_heart, {&pkj_leaves_LV});
	TreeBodyRelationInner pkj_inner_LV(pkj_body_LV);

	/** Corrected configuration. */
	solid_dynamics::CorrectConfiguration correct_configuration_excitation(physiology_heart_inner);
	/** Time step size calculation. */
	electro_physiology::GetElectroPhysiologyTimeStepSize get_myocardium_physiology_time_step(physiology_heart);
	/** Diffusion process for diffusion body. */
	electro_physiology::ElectroPhysiologyDiffusionRelaxationComplex myocardium_diffusion_relaxation_RV(physiology_heart_complex_RV);
	electro_physiology::ElectroPhysiologyDiffusionRelaxationComplex myocardium_diffusion_relaxation_LV(physiology_heart_complex_LV);
	/** Solvers for ODE system */
	electro_physiology::ElectroPhysiologyReactionRelaxationForward myocardium_reaction_relaxation_forward(physiology_heart);
	electro_physiology::ElectroPhysiologyReactionRelaxationBackward myocardium_reaction_relaxation_backward(physiology_heart);
	/** Physiology for PKJ*/
	/** Time step size calculation. */
	electro_physiology::GetElectroPhysiologyTimeStepSize get_pkj_physiology_time_step(pkj_body_RV);
	electro_physiology::ElectroPhysiologyDiffusionRelaxationInner pkj_diffusion_relaxation_RV(pkj_inner_RV);
	electro_physiology::ElectroPhysiologyDiffusionRelaxationInner pkj_diffusion_relaxation_LV(pkj_inner_LV);
	/** Solvers for ODE system */
	electro_physiology::ElectroPhysiologyReactionRelaxationForward pkj_reaction_relaxation_forward_RV(pkj_body_RV);
	electro_physiology::ElectroPhysiologyReactionRelaxationBackward pkj_reaction_relaxation_backward_RV(pkj_body_RV);
	electro_physiology::ElectroPhysiologyReactionRelaxationForward pkj_reaction_relaxation_forward_LV(pkj_body_LV);
	electro_physiology::ElectroPhysiologyReactionRelaxationBackward pkj_reaction_relaxation_backward_LV(pkj_body_LV);
	/**IO for observer.*/
	BodyStatesRecordingToVtp write_states(in_output, system.real_bodies_);
	ObservedQuantityRecording<Real> write_voltage("Voltage", in_output, voltage_observer_contact);
	ObservedQuantityRecording<Vecd> write_displacement("Position", in_output, myocardium_observer_contact);
	/**Apply the Iron stimulus.*/
	ApplyStimulusCurrentToMmyocardium apply_stimulus_myocardium(physiology_heart);
	ApplyStimulusCurrentToPKJ apply_stimulus_pkj_RV(pkj_body_RV);
	ApplyStimulusCurrentToPKJ apply_stimulus_pkj_LV(pkj_body_LV);
	/** Active mechanics. */
	solid_dynamics::CorrectConfiguration correct_configuration_contraction(mechanics_heart_inner);
	/** Observer Dynamics */
	observer_dynamics::CorrectInterpolationKernelWeights
		correct_kernel_weights_for_interpolation(mechanics_heart_contact);
	/** Interpolate the active contract stress from electrophysiology body. */
	observer_dynamics::InterpolatingAQuantity<Real>
		active_stress_interpolation(mechanics_heart_contact, "ActiveContractionStress", "ActiveContractionStress");
	/** Interpolate the particle position in physiology_heart  from mechanics_heart. */
	observer_dynamics::InterpolatingAQuantity<Vecd>
		interpolation_particle_position(physiology_heart_contact, "Position", "Position");
	/** Time step size calculation. */
	solid_dynamics::AcousticTimeStepSize get_mechanics_time_step(mechanics_heart);
	/** active and passive stress relaxation. */
	solid_dynamics::StressRelaxationFirstHalf stress_relaxation_first_half(mechanics_heart_inner);
	solid_dynamics::StressRelaxationSecondHalf stress_relaxation_second_half(mechanics_heart_inner);
	/** Constrain region of the inserted body. */
	MuscleBaseShapeParameters muscle_base_parameters;
	BodyRegionByParticle muscle_base(mechanics_heart,  makeShared<TriangleMeshShapeBrick>(muscle_base_parameters, "Holder"));
	solid_dynamics::ConstrainSolidBodyRegion constrain_holder(mechanics_heart, muscle_base);
	/** 
	 * Pre-simultion. 
	 */
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	correct_configuration_excitation.parallel_exec();
	correct_configuration_contraction.parallel_exec();
	correct_kernel_weights_for_interpolation.parallel_exec();
	/**Output global basic parameters. */
	write_states.writeToFile(0);
	write_voltage.writeToFile(0);
	write_displacement.writeToFile(0);
	write_states.writeToFile(0);
	/**
	 * main loop. 
	 */
	int screen_output_interval = 10;
	int ite = 0;
	int reaction_step = 2;
	Real End_Time = 80;
	Real Ouput_T = End_Time / 200.0;
	Real Observer_time = 0.01 * Ouput_T;
	Real dt_myocardium = 0.0;
	Real dt_pkj = 0.0;
	Real dt_muscle = 0.0;
	/** Statistics for computing time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	cout << "Main Loop Starts Here : "
		 << "\n";
	/** Main loop starts here. */
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		while (integration_time < Ouput_T)
		{
			Real relaxation_time = 0.0;
			while (relaxation_time < Observer_time)
			{
				if (ite % screen_output_interval == 0)
				{
					cout << fixed << setprecision(9) << "N=" << ite << "	Time = "
						 << GlobalStaticVariables::physical_time_
						 << "	dt_pkj = " << dt_pkj
						 << "	dt_myocardium = " << dt_myocardium
						 << "	dt_muscle = " << dt_muscle << "\n";
				}
				/** Apply stimulus excitation. */
				// if( 0 <= GlobalStaticVariables::physical_time_
				// 	&&  GlobalStaticVariables::physical_time_ <= 0.5)
				// {
				// 	apply_stimulus_myocardium.parallel_exec(dt_myocardium);
				// }

				Real dt_pkj_sum = 0.0;
				while (dt_pkj_sum < dt_myocardium)
				{
					/**
					 * When network generates particles, the final particle spacing, which is after particle projected in to 
					 * complex geometry, may small than the reference one, therefore, a smaller time step size is required. 
					 */
					dt_pkj = 0.5 * get_pkj_physiology_time_step.parallel_exec();
					if (dt_myocardium - dt_pkj_sum < dt_pkj)
						dt_pkj = dt_myocardium - dt_pkj_sum;

					if (0 <= GlobalStaticVariables::physical_time_ && GlobalStaticVariables::physical_time_ <= 0.5)
					{
						apply_stimulus_pkj_RV.parallel_exec(dt_pkj);
						apply_stimulus_pkj_LV.parallel_exec(dt_pkj);
					}
					/**Strang splitting method. */
					int ite_pkj_forward = 0;
					while (ite_pkj_forward < reaction_step)
					{
						pkj_reaction_relaxation_forward_RV.parallel_exec(0.5 * dt_pkj / Real(reaction_step));
						pkj_reaction_relaxation_forward_LV.parallel_exec(0.5 * dt_pkj / Real(reaction_step));
						ite_pkj_forward++;
					}
					/** 2nd Runge-Kutta scheme for diffusion. */
					pkj_diffusion_relaxation_RV.parallel_exec(dt_pkj);
					pkj_diffusion_relaxation_LV.parallel_exec(dt_pkj);
					//backward reaction
					int ite_pkj_backward = 0;
					while (ite_pkj_backward < reaction_step)
					{
						pkj_reaction_relaxation_backward_RV.parallel_exec(0.5 * dt_pkj / Real(reaction_step));
						pkj_reaction_relaxation_backward_LV.parallel_exec(0.5 * dt_pkj / Real(reaction_step));
						ite_pkj_backward++;
					}

					dt_pkj_sum += dt_pkj;
				}

				/**Strang splitting method. */
				int ite_forward = 0;
				while (ite_forward < reaction_step)
				{
					myocardium_reaction_relaxation_forward.parallel_exec(0.5 * dt_myocardium / Real(reaction_step));
					ite_forward++;
				}
				/** 2nd Runge-Kutta scheme for diffusion. */
				myocardium_diffusion_relaxation_RV.parallel_exec(dt_myocardium);
				myocardium_diffusion_relaxation_LV.parallel_exec(dt_myocardium);

				//backward reaction
				int ite_backward = 0;
				while (ite_backward < reaction_step)
				{
					myocardium_reaction_relaxation_backward.parallel_exec(0.5 * dt_myocardium / Real(reaction_step));
					ite_backward++;
				}

				active_stress_interpolation.parallel_exec();
				Real dt_muscle_sum = 0.0;
				while (dt_muscle_sum < dt_myocardium)
				{
					dt_muscle = get_mechanics_time_step.parallel_exec();
					if (dt_myocardium - dt_muscle_sum < dt_muscle)
						dt_muscle = dt_myocardium - dt_muscle_sum;
					stress_relaxation_first_half.parallel_exec(dt_muscle);
					constrain_holder.parallel_exec(dt_muscle);
					stress_relaxation_second_half.parallel_exec(dt_muscle);
					dt_muscle_sum += dt_muscle;
				}

				ite++;
				dt_myocardium = get_myocardium_physiology_time_step.parallel_exec();

				relaxation_time += dt_myocardium;
				integration_time += dt_myocardium;
				GlobalStaticVariables::physical_time_ += dt_myocardium;
			}
			write_voltage.writeToFile(ite);
			write_displacement.writeToFile(ite);
		}
		tick_count t2 = tick_count::now();
		interpolation_particle_position.parallel_exec();
		write_states.writeToFile();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	cout << "Total wall time for computation: " << tt.seconds() << " seconds." << endl;

	return 0;
}