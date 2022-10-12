/**
 * @file 	excitation-contraction.cpp
 * @brief 	This is the case studying the electromechanics on a biventricular heart model in 3D.
 * @author 	Chi Zhang and Xiangyu Hu
 * 			Unit :
 *			time t = ms = 12.9 [-]
 * 			length l = mm
 * 			mass m = g
 *			density rho = g * (mm)^(-3)
 *			Pressure pa = g * (mm)^(-1) * (ms)^(-2)
 *			diffusion d = (mm)^(2) * (ms)^(-2)
 */
#include "sphinxsys.h" // SPHinXsys Library.
using namespace SPH; // Namespace cite here.
/** Geometry parameter. */
/** Set the file path to the stl file. */
std::string full_path_to_stl_file = "./input/heart-new.stl";
Real length_scale = 1.0;
Real time_scale = 1.0 / 12.9;
Real stress_scale = 1.0e-6;
/** Parameters and physical properties. */
Vec3d domain_lower_bound(-55.0 * length_scale, -75.0 * length_scale, -35.0 * length_scale);
Vec3d domain_upper_bound(35.0 * length_scale, 5.0 * length_scale, 35.0 * length_scale);
Real dp_0 = (domain_upper_bound[0] - domain_lower_bound[0]) / 45.0; /**< Initial particle spacing. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(domain_lower_bound, domain_upper_bound);

/** Material properties. */
Real rho0_s = 1.06e-3;
/** Active stress factor */
Real k_a = 1.0e2 * stress_scale;
// Real a0[4] = {496.0 * stress_scale, 15196.0 * stress_scale, 3283.0 * stress_scale, 662.0 * stress_scale};
// Real b0[4] = {7.209, 20.417, 11.176, 9.466};
Real a0[4] = {59.0 * stress_scale, 18472.0 * stress_scale, 2841.0 * stress_scale, 216.0 * stress_scale};
Real b0[4] = {8.023, 16.026, 11.12, 11.436};
/** reference stress to achieve weakly compressible condition */
Real poisson = 0.4995;
Real bulk_modulus = 5.0e6 * stress_scale; // 2.0 * a0[0] * (1.0 + poisson) / (3.0 * (1.0 - 2.0 * poisson));
/** Electrophysiology parameters. */
StdVec<std::string> species_name_list{"Phi"};
Real diffusion_coff = 1.0;
Real bias_coff = 0.0;
/** Electrophysiology parameters. */
Real c_m = 1.0;
Real k = 8.0;
Real a = 0.01;
Real b = 0.15;
Real mu_1 = 0.2;
Real mu_2 = 0.3;
Real epsilon = 0.002;
/** Fibers and sheet. */
Vec3d fiber_direction(1.0, 0.0, 0.0);
Vec3d sheet_direction(0.0, 1.0, 0.0);

/**
 * Define heart geometry
 */
class Heart : public ComplexShape
{
public:
	explicit Heart(const std::string &shape_name) : ComplexShape(shape_name)
	{
		Vecd translation(-53.5 * length_scale, -70.0 * length_scale, -32.5 * length_scale);
		add<TriangleMeshShapeSTL>(full_path_to_stl_file, translation, length_scale);
	}
};
//----------------------------------------------------------------------
//	Setup diffusion material properties.
//----------------------------------------------------------------------
class FiberDirectionDiffusion : public DiffusionReaction<LocallyOrthotropicMuscle>
{
public:
	FiberDirectionDiffusion() : DiffusionReaction<LocallyOrthotropicMuscle>(
									species_name_list, rho0_s, bulk_modulus, fiber_direction, sheet_direction, a0, b0)
	{
		initializeAnDiffusion<IsotropicDiffusion>("Phi", "Phi", diffusion_coff);
	};
};
/** Set diffusion relaxation method. */
class DiffusionRelaxation
	: public RelaxationOfAllDiffusionSpeciesRK2<
		  RelaxationOfAllDiffussionSpeciesInner<SolidBody, ElasticSolidParticles, LocallyOrthotropicMuscle>>
{
public:
	explicit DiffusionRelaxation(BodyRelationInner &body_inner_relation)
		: RelaxationOfAllDiffusionSpeciesRK2(body_inner_relation){};
	virtual ~DiffusionRelaxation(){};
};
/** Imposing diffusion boundary condition */
class DiffusionBCs
	: public ConstrainDiffusionBodyRegion<SolidBody, ElasticSolidParticles, BodySurface, LocallyOrthotropicMuscle>
{
protected:
	size_t phi_;
	virtual void Update(size_t index_i, Real dt = 0.0) override
	{
		Vecd dist_2_face = body_->body_shape_->findNormalDirection(pos_n_[index_i]);
		Vecd face_norm = dist_2_face / (dist_2_face.norm() + 1.0e-15);

		Vecd center_norm = pos_n_[index_i] / (pos_n_[index_i].norm() + 1.0e-15);

		Real angle = dot(face_norm, center_norm);
		if (angle >= 0.0)
		{
			species_n_[phi_][index_i] = 1.0;
		}
		else
		{
			if (pos_n_[index_i][1] < -body_->sph_adaptation_->ReferenceSpacing())
				species_n_[phi_][index_i] = 0.0;
		}
	};

public:
	DiffusionBCs(SolidBody &body, BodySurface &body_part)
		: ConstrainDiffusionBodyRegion<SolidBody, ElasticSolidParticles, BodySurface, LocallyOrthotropicMuscle>(body, body_part)
	{
		phi_ = material_->SpeciesIndexMap()["Phi"];
	};
	virtual ~DiffusionBCs(){};
};
/** Compute Fiber and Sheet direction after diffusion */
class ComputeFiberandSheetDirections
	: public DiffusionBasedMapping<SolidBody, ElasticSolidParticles, LocallyOrthotropicMuscle>
{
protected:
	size_t phi_;
	Real beta_epi_, beta_endo_;
	/** We define the centerline vector, which is parallel to the ventricular centerline and pointing  apex-to-base.*/
	Vecd center_line_;
	virtual void Update(size_t index_i, Real dt = 0.0) override
	{
		/**
		 * Ref: original doi.org/10.1016/j.euromechsol.2013.10.009
		 * 		Present  doi.org/10.1016/j.cma.2016.05.031
		 */
		/** Probe the face norm from Levelset field. */
		Vecd dist_2_face = body_->body_shape_->findNormalDirection(pos_n_[index_i]);
		Vecd face_norm = dist_2_face / (dist_2_face.norm() + 1.0e-15);
		Vecd center_norm = pos_n_[index_i] / (pos_n_[index_i].norm() + 1.0e-15);
		if (dot(face_norm, center_norm) <= 0.0)
		{
			face_norm = -face_norm;
		}
		/** Compute the centerline's projection on the plane orthogonal to face norm. */
		Vecd circumferential_direction = SimTK::cross(center_line_, face_norm);
		Vecd cd_norm = circumferential_direction / (circumferential_direction.norm() + 1.0e-15);
		/** The rotation angle is given by beta = (beta_epi - beta_endo) phi + beta_endo */
		Real beta = (beta_epi_ - beta_endo_) * species_n_[phi_][index_i] + beta_endo_;
		/** Compute the rotation matrix through Rodrigues rotation formulation. */
		Vecd f_0 = cos(beta) * cd_norm + sin(beta) * SimTK::cross(face_norm, cd_norm) +
				   dot(face_norm, cd_norm) * (1.0 - cos(beta)) * face_norm;

		if (pos_n_[index_i][1] < -body_->sph_adaptation_->ReferenceSpacing())
		{
			material_->local_f0_[index_i] = f_0 / (f_0.norm() + 1.0e-15);
			material_->local_s0_[index_i] = face_norm;
		}
		else
		{
			material_->local_f0_[index_i] = Vecd(0);
			material_->local_s0_[index_i] = Vecd(0);
		}
	};

public:
	explicit ComputeFiberandSheetDirections(SolidBody &body)
		: DiffusionBasedMapping<SolidBody, ElasticSolidParticles, LocallyOrthotropicMuscle>(body)
	{
		phi_ = material_->SpeciesIndexMap()["Phi"];
		center_line_ = Vecd(0.0, 1.0, 0.0);
		beta_epi_ = -(70.0 / 180.0) * M_PI;
		beta_endo_ = (80.0 / 180.0) * M_PI;
	};
	virtual ~ComputeFiberandSheetDirections(){};
};
//	define shape parameters which will be used for the constrained body part.
class MuscleBaseShapeParameters : public TriangleMeshShapeBrick::ShapeParameters
{
public:
	MuscleBaseShapeParameters() : TriangleMeshShapeBrick::ShapeParameters()
	{
		Real l = domain_upper_bound[0] - domain_lower_bound[0];
		Real w = domain_upper_bound[2] - domain_lower_bound[2];
		halfsize_ = Vec3d(0.5 * l, 1.0 * dp_0, 0.5 * w);
		resolution_ = 20;
		translation_ = Vec3d(-10.0 * length_scale, -1.0 * dp_0, 0.0);
	}
};
//	application dependent initial condition
class ApplyStimulusCurrentSI
	: public electro_physiology::ElectroPhysiologyInitialCondition
{
protected:
	size_t voltage_;

	void Update(size_t index_i, Real dt) override
	{
		// if (-30.0 * length_scale <= pos_n_[index_i][0] && pos_n_[index_i][0] <= -15.0 * length_scale)
		if (-11.0 * length_scale <= pos_n_[index_i][0] && pos_n_[index_i][0] <= -4.0 * length_scale)
		{
			// if (-2.0 * length_scale <= pos_n_[index_i][1] && pos_n_[index_i][1] <= 0.0)
			if (-68.7042 * length_scale <= pos_n_[index_i][1] && pos_n_[index_i][1] <= -66.7042)
			{
				// if (-3.0 * length_scale <= pos_n_[index_i][2] && pos_n_[index_i][2] <= 3.0 * length_scale)
				if (-1.9806 * length_scale <= pos_n_[index_i][2] && pos_n_[index_i][2] <= 4.0194 * length_scale)
				{
					species_n_[voltage_][index_i] = 1.2;
				}
			}
		}
	};

public:
	explicit ApplyStimulusCurrentSI(SolidBody &muscle)
		: electro_physiology::ElectroPhysiologyInitialCondition(muscle)
	{
		voltage_ = material_->SpeciesIndexMap()["Voltage"];
	};
};
/**
 * application dependent initial condition
 */
class ApplyStimulusCurrentSII
	: public electro_physiology::ElectroPhysiologyInitialCondition
{
protected:
	size_t voltage_;

	void Update(size_t index_i, Real dt) override
	{
		if (0.0 <= pos_n_[index_i][0] && pos_n_[index_i][0] <= 6.0 * length_scale)
		{
			if (-6.0 * length_scale <= pos_n_[index_i][1])
			{
				if (12.0 * length_scale <= pos_n_[index_i][2])
				{
					species_n_[voltage_][index_i] = 0.95;
				}
			}
		}
	};

public:
	explicit ApplyStimulusCurrentSII(SolidBody &muscle)
		: electro_physiology::ElectroPhysiologyInitialCondition(muscle)
	{
		voltage_ = material_->SpeciesIndexMap()["Voltage"];
	};
};
/**
 * define observer particle generator.
 */
class HeartObserverParticleGenerator : public ObserverParticleGenerator
{
public:
	explicit HeartObserverParticleGenerator(SPHBody &sph_body) : ObserverParticleGenerator(sph_body)
	{
		/** position and volume. */
		positions_.push_back(Vecd(-45.0 * length_scale, -30.0 * length_scale, 0.0));
		positions_.push_back(Vecd(0.0, -30.0 * length_scale, 26.0 * length_scale));
		positions_.push_back(Vecd(-30.0 * length_scale, -50.0 * length_scale, 0.0));
		positions_.push_back(Vecd(0.0, -50.0 * length_scale, 20.0 * length_scale));
		positions_.push_back(Vecd(0.0, -70.0 * length_scale, 0.0));
	}
};
/**
 * The main program.
 */
int main(int ac, char *av[])
{
	//----------------------------------------------------------------------
	//	SPHSystem section
	//----------------------------------------------------------------------
	SPHSystem system(system_domain_bounds, dp_0);
	/** Set the starting time. */
	GlobalStaticVariables::physical_time_ = 0.0;
	/** Tag for run particle relaxation for the initial body fitted distribution. */
	system.run_particle_relaxation_ = false;
	/** Tag for reload initially relaxed particles. */
	system.reload_particles_ = true;
	/** Tag for computation from restart files. 0: not from restart files. */
	system.restart_step_ = 0;
// handle command line arguments
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
		SolidBody herat_model(system, makeShared<Heart>("HeartModel"));
		herat_model.defineBodyLevelSetShape()->writeLevelSet(herat_model);
		herat_model.defineParticlesAndMaterial<DiffusionReactionParticles<ElasticSolidParticles>, FiberDirectionDiffusion>();
		herat_model.generateParticles<ParticleGeneratorLattice>();
		/** topology */
		BodyRelationInner herat_model_inner(herat_model);
		/** Random reset the relax solid particle position. */
		RandomizeParticlePosition random_particles(herat_model);
		/** A  Physics relaxation step. */
		relax_dynamics::RelaxationStepInner relaxation_step_inner(herat_model_inner);
		/** Time step for diffusion. */
		GetDiffusionTimeStepSize<SolidBody, ElasticSolidParticles, LocallyOrthotropicMuscle> get_time_step_size(herat_model);
		/** Diffusion process for diffusion body. */
		DiffusionRelaxation diffusion_relaxation(herat_model_inner);
		/** Compute the fiber and sheet after diffusion. */
		ComputeFiberandSheetDirections compute_fiber_sheet(herat_model);
		/** Write the body state to Vtp file. */
		BodyStatesRecordingToVtp write_herat_model_state_to_vtp(in_output, {herat_model});
		/** Write the particle reload files. */
		ReloadParticleIO write_particle_reload_files(in_output, {&herat_model}, {"HeartModel"});
		/** Write material property to xml file. */
		ReloadMaterialParameterIO write_material_property(in_output, herat_model.base_material_, "FiberDirection");
		//----------------------------------------------------------------------
		//	Physics relaxation starts here.
		//----------------------------------------------------------------------
		random_particles.parallel_exec(0.25);
		relaxation_step_inner.surface_bounding_.parallel_exec();
		write_herat_model_state_to_vtp.writeToFile(0.0);
		//----------------------------------------------------------------------
		// From here the time stepping begins.
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
				write_herat_model_state_to_vtp.writeToFile(ite);
			}
		}

		BodySurface surface_part(herat_model);
		/** constraint boundary condition for diffusion. */
		DiffusionBCs impose_diffusion_bc(herat_model, surface_part);
		impose_diffusion_bc.parallel_exec();

		write_herat_model_state_to_vtp.writeToFile(ite);

		Real dt = get_time_step_size.parallel_exec();
		while (ite <= diffusion_step + relax_step)
		{
			diffusion_relaxation.parallel_exec(dt);
			impose_diffusion_bc.parallel_exec();
			if (ite % 10 == 0)
			{
				std::cout << "Diffusion steps N=" << ite - relax_step << "	dt: " << dt << "\n";
				write_herat_model_state_to_vtp.writeToFile(ite);
			}
			ite++;
		}
		compute_fiber_sheet.exec();
		ite++;
		write_herat_model_state_to_vtp.writeToFile(ite);
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
	SolidBody mechanics_heart(system, makeShared<Heart>("MechanicalHeart"));
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
	//----------------------------------------------------------------------
	//	SPH Observation section
	//----------------------------------------------------------------------
	ObserverBody voltage_observer(system, "VoltageObserver");
	voltage_observer.generateParticles<HeartObserverParticleGenerator>();
	ObserverBody myocardium_observer(system, "MyocardiumObserver");
	myocardium_observer.generateParticles<HeartObserverParticleGenerator>();
	//----------------------------------------------------------------------
	//	SPHBody relation (topology) section
	//----------------------------------------------------------------------
	BodyRelationInner physiology_heart_inner(physiology_heart);
	BodyRelationInner mechanics_body_inner(mechanics_heart);
	BodyRelationContact physiology_heart_contact(physiology_heart, {&mechanics_heart});
	BodyRelationContact mechanics_body_contact(mechanics_heart, {&physiology_heart});
	BodyRelationContact voltage_observer_contact(voltage_observer, {&physiology_heart});
	BodyRelationContact myocardium_observer_contact(myocardium_observer, {&mechanics_heart});
	//----------------------------------------------------------------------
	//	SPH Method section
	//----------------------------------------------------------------------
	// Corrected configuration.
	solid_dynamics::CorrectConfiguration correct_configuration_excitation(physiology_heart_inner);
	// Time step size calculation.
	electro_physiology::GetElectroPhysiologyTimeStepSize get_physiology_time_step(physiology_heart);
	// Diffusion process for diffusion body.
	electro_physiology::ElectroPhysiologyDiffusionRelaxationInner diffusion_relaxation(physiology_heart_inner);
	// Solvers for ODE system.
	electro_physiology::ElectroPhysiologyReactionRelaxationForward reaction_relaxation_forward(physiology_heart);
	electro_physiology::ElectroPhysiologyReactionRelaxationBackward reaction_relaxation_backward(physiology_heart);
	//	Apply the Iron stimulus.
	ApplyStimulusCurrentSI apply_stimulus_s1(physiology_heart);
	ApplyStimulusCurrentSII apply_stimulus_s2(physiology_heart);
	// Active mechanics.
	solid_dynamics::CorrectConfiguration correct_configuration_contraction(mechanics_body_inner);
	observer_dynamics::CorrectInterpolationKernelWeights correct_kernel_weights_for_interpolation(mechanics_body_contact);
	/** Interpolate the active contract stress from electrophysiology body. */
	observer_dynamics::InterpolatingAQuantity<Real>
		active_stress_interpolation(mechanics_body_contact, "ActiveContractionStress", "ActiveContractionStress");
	/** Interpolate the particle position in physiology_heart  from mechanics_heart. */
	observer_dynamics::InterpolatingAQuantity<Vecd>
		interpolation_particle_position(physiology_heart_contact, "Position", "Position");
	/** Time step size calculation. */
	solid_dynamics::AcousticTimeStepSize get_mechanics_time_step(mechanics_heart);
	/** active and passive stress relaxation. */
	solid_dynamics::StressRelaxationFirstHalf stress_relaxation_first_half(mechanics_body_inner);
	solid_dynamics::StressRelaxationSecondHalf stress_relaxation_second_half(mechanics_body_inner);
	SimpleDynamics<NormalDirectionFromBodyShape>(mechanics_heart).parallel_exec();
	solid_dynamics::UpdateElasticNormalDirection body_update_normal(mechanics_heart);
	/** Constrain region of the inserted body. */
	MuscleBaseShapeParameters muscle_base_parameters;
	BodyRegionByParticle muscle_base(mechanics_heart, makeShared<TriangleMeshShapeBrick>(muscle_base_parameters, "Holder"));
	solid_dynamics::ConstrainSolidBodyRegion constrain_holder(mechanics_heart, muscle_base);
	BodySurface surface_layer(mechanics_heart);
	IndexVector all_surface_ids = surface_layer.body_part_particles_;
	auto myocardium_particles = dynamic_cast<ElasticSolidParticles*>(mechanics_heart.base_particles_);
	// Boundary conditions
	// Pressure initial conditions and conversion factors
	Real mmHg_to_Pa = 133.3224;
	Real initial_IVC_pressure_LV = 0.0; // 10.0 *  mmHg_to_Pa * 1e-6;
	Real initial_IVC_pressure_RV = 0.0; // 6.0 *  mmHg_to_Pa * 1e-6;
	Real initial_ejection_pressure_LV = 50.0 * mmHg_to_Pa * stress_scale; // 50mmHg initial pressure at start of ejection phase
	Real initial_ejection_pressure_RV = 15.0 * mmHg_to_Pa * stress_scale; // 15mmHg initial pressure at start of ejection phase
	Real duration_IVC = 50.0; // duration isovolumic constraction phase = 50ms
	Real IVC_pressure_slope_LV = (initial_ejection_pressure_LV - initial_IVC_pressure_LV)  / duration_IVC;
	Real IVC_pressure_slope_RV = (initial_ejection_pressure_RV - initial_IVC_pressure_RV) / duration_IVC;
	BodyPartByParticle LV_body(mechanics_heart, "LeftVentricle");
	Vecd source_point_LV = Vec3d(11.94363656084621, -24.4540245, -1.36504656); // Vecd(3.265584517 , -24.4540245 , -1.36504656);
	StdVec<array<Real, 2>> pressure_over_time_LV = {
	 	{0.0, 0.0},
	 	{duration_IVC, initial_ejection_pressure_LV}};
	// BodyPartByParticle RV_body(mechanics_heart, "RightVentricle");
	// Vecd source_point_RV = Vecd(-38.46060618712892, -18.19030649377814, -0.20931284741483314);
	// StdVec<array<Real, 2>> pressure_over_time_RV = {
	// 	{0.0, 0.0},
	// 	{duration_IVC, initial_ejection_pressure_RV}};
	solid_dynamics::SurfacePressureFromSource surface_pressure_LV(mechanics_heart, LV_body, source_point_LV, pressure_over_time_LV);
	StdVec<Vecd> LV_particles_by_sourcepoint_pos_0 = surface_pressure_LV.GetPressureParticles();
	IndexVector pressure_particle_ids = surface_pressure_LV.GetPressureParticlesIDs();
	ObserverBody LV_source_point_observer(system, "LV_source_point_observer");
	LV_source_point_observer.generateParticles<ObserverParticleGenerator>(LV_particles_by_sourcepoint_pos_0);
	// solid_dynamics::SurfacePressureFromSource surface_pressure_RV(mechanics_heart, RV_body, source_point_RV, pressure_over_time_RV);
	// StdVec<Vecd> RV_particles_by_sourcepoint_pos_0 = surface_pressure_RV.GetPressureParticles();
	// ObserverBody RV_source_point_observer(system, "RV_source_point_observer");
	// RV_source_point_observer.generateParticles<ObserverParticleGenerator>(RV_particles_by_sourcepoint_pos_0);
	Real stiffness = 50.0 * 1e3 * stress_scale; // 0.1 kPa/mm --> 100 Pa/mm --> 100 * 1e-6 * [pressure unit] / mm --> 1e-4
	Real damping_ratio = 0.0; // 1 kPa*s/mm --> 1e3 Pa*s/mm --> 1e3 * 1e-6 * [pressure unit] * s / mm --> 1e3 * 1e-6 * 1e3 [pressure unit] * ms / mm --> 1
	Vecd outer_soure_point = Vecd(3.265584517, -75.0, -1.36504656);
	BodyPartByParticle body_part(mechanics_heart, "Pericardium");
	solid_dynamics::SpringNormalOnSurfaceParticles spring_damper_constraint(mechanics_heart, body_part, true, outer_soure_point, stiffness, damping_ratio);
	//----------------------------------------------------------------------
	//	SPH Output section
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp write_states(in_output, system.real_bodies_);
	BodyStatesRecordingToVtp write_observers(in_output, system.observation_bodies_);
	RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Real>>
		write_voltage("Voltage", in_output, voltage_observer_contact);
	RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
		write_displacement("Position", in_output, myocardium_observer_contact);
	//----------------------------------------------------------------------
	//	 Pre-simulation.
	//----------------------------------------------------------------------
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	correct_configuration_excitation.parallel_exec();
	correct_configuration_contraction.parallel_exec();
	correct_kernel_weights_for_interpolation.parallel_exec();
	/** Output initial states and observations */
	write_states.writeToFile(0);
	write_observers.writeToFile(0);
	write_voltage.writeToFile(0);
	write_displacement.writeToFile(0);
	//----------------------------------------------------------------------
	//	 Physical parameters for main loop.
	//----------------------------------------------------------------------
	int screen_output_interval = 10;
	int ite = 0;
	int reaction_step = 2;
	Real End_Time = 100;
	Real Ouput_T = End_Time / 200.0;
	Real Observer_time = 0.01 * Ouput_T;
	Real dt = 0.0;	 /**< Default acoustic time step sizes for physiology. */
	Real dt_s = 0.0; /**< Default acoustic time step sizes for mechanics. */
	/** Statistics for computing time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	std::cout << "Main Loop Starts Here : "
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
					std::cout << std::fixed << std::setprecision(9) << "N=" << ite << "	Time = "
							  << GlobalStaticVariables::physical_time_
							  << "	dt = " << dt
							  << "	dt_s = " << dt_s << "\n";
				}
				/** Apply stimulus excitation. */
				if (0 <= GlobalStaticVariables::physical_time_ && GlobalStaticVariables::physical_time_ <= 0.5)
				{
					apply_stimulus_s1.parallel_exec(dt);
				}
				/** Single spiral wave. */
				// if( 60 <= GlobalStaticVariables::physical_time_
				// 	&&  GlobalStaticVariables::physical_time_ <= 65)
				// {
				// 	apply_stimulus_s2.parallel_exec(dt);
				// }
				/**Strong splitting method. */
				// forward reaction
				int ite_forward = 0;
				while (ite_forward < reaction_step)
				{
					reaction_relaxation_forward.parallel_exec(0.5 * dt / Real(reaction_step));
					ite_forward++;
				}
				/** 2nd Runge-Kutta scheme for diffusion. */
				diffusion_relaxation.parallel_exec(dt);

				// backward reaction
				int ite_backward = 0;
				while (ite_backward < reaction_step)
				{
					reaction_relaxation_backward.parallel_exec(0.5 * dt / Real(reaction_step));
					ite_backward++;
				}

				active_stress_interpolation.parallel_exec();

				Real dt_s_sum = 0.0;
				while (dt_s_sum < dt)
				{
					dt_s = get_mechanics_time_step.parallel_exec();
					if (dt - dt_s_sum < dt_s)
						dt_s = dt - dt_s_sum;
					auto gravity = Gravity({0,0,0});
                    TimeStepInitialization(mechanics_heart, gravity).parallel_exec(dt_s);
					body_update_normal.parallel_exec(dt_s);
					surface_pressure_LV.parallel_exec(dt_s);
					spring_damper_constraint.parallel_exec(dt_s);
					// surface_pressure_RV.parallel_exec(dt_s);
					stress_relaxation_first_half.parallel_exec(dt_s);
					// constrain_holder.parallel_exec(dt_s);
					stress_relaxation_second_half.parallel_exec(dt_s);
					Real current_velocity = 0.0;
					for(size_t i = 0; i < pressure_particle_ids.size(); i++)
					{
						size_t index_j = pressure_particle_ids[i];
						if(std::find(all_surface_ids.begin(), all_surface_ids.end(), index_j) != all_surface_ids.end())
						{
							current_velocity += myocardium_particles->vel_n_[index_j].norm();
						}
					}
					if (ite % 20 == 0) std::cout << "current surface velocities: " << current_velocity << std::endl;
					dt_s_sum += dt_s;
				}

				ite++;
				dt = get_physiology_time_step.parallel_exec();

				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
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
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	return 0;
}
