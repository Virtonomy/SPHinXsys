
// #include "arterial_guidewire.h" /**< Case setup for this example. */
// #include "base_data_type.h"
#include "contact_body_relation.h"
#include "data_type.h"
#include "sphinxsys.h"
#include "type_wrapper.h"
#include <SimTKcommon/internal/UnitVec.h>
#include "level_set.h"

using namespace SPH;

Real get_physical_viscosity_general(Real rho, Real youngs_modulus, Real length_scale, Real shape_constant = 0.4)
{
    // the physical viscosity is defined in the paper of prof. Hu
    // https://arxiv.org/pdf/2103.08932.pdf
    // physical_viscosity = (beta / 4) * sqrt(rho * Young's modulus) * L
    // beta: shape constant (0.4 for beam)
    return shape_constant / 4.0 * std::sqrt(rho * youngs_modulus) * length_scale;
}

template <typename... Ts>
SPH::BoundingBox union_bounding_box(const SPH::BoundingBox &a, const SPH::BoundingBox &b, Ts &&...args)
{
    SPH::BoundingBox out = a;
    out.first_[0] = std::min(a.first_[0], b.first_[0]);
    out.first_[1] = std::min(a.first_[1], b.first_[1]);
    out.first_[2] = std::min(a.first_[2], b.first_[2]);
    out.second_[0] = std::max(a.second_[0], b.second_[0]);
    out.second_[1] = std::max(a.second_[1], b.second_[1]);
    out.second_[2] = std::max(a.second_[2], b.second_[2]);
    if constexpr (sizeof...(args) > 0)
        return union_bounding_box(out, args...);
    else
        return out;
}

struct SolidAlgorithm
{
    InteractionWithUpdate<CorrectedConfigurationInner> corrected_configuration_;
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> stress_relaxation_first_half_;
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half_;
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d>>> position_damping_;
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> time_step_size_;

    SolidAlgorithm(InnerRelation &inner_relation, Real physical_viscosity)
        : corrected_configuration_(inner_relation),
          stress_relaxation_first_half_(inner_relation),
          stress_relaxation_second_half_(inner_relation),
          position_damping_(0.2, inner_relation, "Velocity", physical_viscosity),
          time_step_size_(inner_relation.getSPHBody())
    {
    }

    inline Real get_time_step_size() { return time_step_size_.exec(); }
    inline void correct_config() { corrected_configuration_.exec(); }
    inline void stress_first_half(Real dt) { stress_relaxation_first_half_.exec(dt); }
    inline void stress_second_half(Real dt) { stress_relaxation_second_half_.exec(dt); }
    inline void damping() { position_damping_.exec(); }

    inline void exec_all(Real dt)
    {
        stress_first_half(dt);
        damping();
        stress_second_half(dt);
    }
};

struct SolidContact
{
    SurfaceContactRelation contact_relation_;
    std::unique_ptr<InteractionDynamics<solid_dynamics::ContactDensitySummation>> contact_density_;
    std::unique_ptr<InteractionDynamics<solid_dynamics::ContactForce>> contact_forces_;

    explicit SolidContact(SolidBody &contact_to, SolidBody &contact_from)
        : contact_relation_(contact_to, {&contact_from})
    {
        contact_to.getBaseParticles().registerSharedVariable<Real>("ContactDensity");
        contact_from.getBaseParticles().registerSharedVariable<Real>("ContactDensity");
        contact_density_ = std::make_unique<InteractionDynamics<solid_dynamics::ContactDensitySummation>>(contact_relation_);
        contact_forces_ = std::make_unique<InteractionDynamics<solid_dynamics::ContactForce>>(contact_relation_);
    }
    inline void update_config() { contact_relation_.updateConfiguration(); }
    inline void exec()
    {
        contact_density_->exec();
        contact_forces_->exec();
    }
};



int main(int ac, char *av[])
{
    const double scale = 1e-3; // mm to m
    const double guidewire_diameter = 1.066*scale;
    const double guidewire_length = 100*scale;
    const double guidewire_rho = 695.48;
    const double guidewire_youngs_modulus = 83e9; //83 GPa
    const double guidewire_poisson = 0.3;
    const double guidewire_insert_velocity_magnitude = 5.;
    const double guidewire_diameter_hole_on_aorta = guidewire_diameter + 0.5*scale;


    const double aorta_diameter = 16.*scale;
    const double aorta_length = 100.*scale;
    // const double thickness_aorta = 2.0*scale;
    const double aorta_thickness = 2.0*scale;

    const double aorta_rho0_s = 1265.0;
    const double aorta_poisson = 0.45;
    const double aorta_Youngs_modulus = 5e4;
    const int SimTK_resolution = 10;
    const double res_ref = guidewire_diameter / 4.0;
    const double res_factor_1 = 1.;
    const double res_factor_2 = 1.;
    const double res_aorta = res_ref / double(res_factor_1);
    const double res_guidewide = res_ref / double(res_factor_2);

    //  Define guidewire shape
    //----------------------------------------------------------------------
    Vec3d guidewire_direction = Vec3d(-0.5, -0.5, 0);
    guidewire_direction.normalize();
    Vec3d guidewire_translation = guidewire_length * 0.5 * -1 * guidewire_direction;
    auto guidewire_shape = makeShared<ComplexShape>("GuideWire");
    guidewire_shape->add<TriangleMeshShapeCylinder>(SimTK::UnitVec3(guidewire_direction[0], guidewire_direction[1], guidewire_direction[2]), guidewire_diameter * 0.5,
                                                    guidewire_length * 0.5, SimTK_resolution,
                                                    guidewire_translation);
    //----------------------------------------------------------------------
    // Vec3d aorta_direction = Vec3d::UnitX();
    auto aorta_shape = makeShared<ComplexShape>("Aorta");
    aorta_shape->add<TriangleMeshShapeCylinder>(SimTK::UnitVec3(1., 0., 0.), aorta_diameter * 0.5 + aorta_thickness,
                                                aorta_length * 0.5, SimTK_resolution,
                                                Vec3d::Zero());
    aorta_shape->subtract<TriangleMeshShapeCylinder>(SimTK::UnitVec3(1., 0., 0.), aorta_diameter * 0.5,
                                                     aorta_length * 0.5, SimTK_resolution,
                                                     Vec3d::Zero());

    aorta_shape->subtract<TriangleMeshShapeCylinder>(SimTK::UnitVec3(guidewire_direction[0], guidewire_direction[1], guidewire_direction[2]), guidewire_diameter_hole_on_aorta * 0.5,
                                                    guidewire_length * 0.5, SimTK_resolution,
                                                    guidewire_translation);
    
    ComplexShape aorta_outer_shape("AortaOuter");
    aorta_outer_shape.add<TriangleMeshShapeCylinder>(SimTK::UnitVec3(1., 0., 0.), aorta_diameter * 0.5 + aorta_thickness,
                                                aorta_length * 0.5, SimTK_resolution,
                                                Vec3d::Zero());

    const auto bb_box = union_bounding_box(guidewire_shape->getBounds(), aorta_shape->getBounds());


    /** Setup the system. Please the make sure the global domain bounds are correctly defined. */
    SPHSystem sph_system(bb_box, res_aorta);
    sph_system.setRunParticleRelaxation(false);
    sph_system.setReloadParticles(true);
#ifdef BOOST_AVAILABLE
    sph_system.handleCommandlineOptions(ac, av);
#endif
    IOEnvironment io_environment(sph_system);

    /** create a body with corresponding material, particles and reaction model. */
    SolidBody guidewire_body(sph_system, guidewire_shape);
    guidewire_body.defineAdaptationRatios(1.15, res_aorta/res_guidewide);
    guidewire_body.defineBodyLevelSetShape();
    guidewire_body.defineParticlesAndMaterial<ElasticSolidParticles,SaintVenantKirchhoffSolid>(
        guidewire_rho, guidewire_youngs_modulus, guidewire_poisson);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? guidewire_body.generateParticles<ParticleGeneratorReload>(io_environment,guidewire_body.getName())
        : guidewire_body.generateParticles<ParticleGeneratorLattice>();

    SolidBody aorta_body(sph_system, aorta_shape);
    aorta_body.defineBodyLevelSetShape();
    aorta_body.defineParticlesAndMaterial<ElasticSolidParticles,SaintVenantKirchhoffSolid>(aorta_rho0_s, aorta_Youngs_modulus, aorta_poisson);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? aorta_body.generateParticles<ParticleGeneratorReload>(io_environment, aorta_body.getName())
        : aorta_body.generateParticles<ParticleGeneratorLattice>();


    LevelSetShape aorta_outer_lvl_shape(guidewire_body,aorta_outer_shape);
    /**body relation topology */
    InnerRelation guidewire_inner(guidewire_body);
    InnerRelation aorta_inner(aorta_body);
    ContactRelation guidewire_aorta_contact(guidewire_body, {&aorta_body});
    ContactRelation aorta_guidewire_contact(aorta_body, {&guidewire_body});

    if (sph_system.RunParticleRelaxation())
    {
        
        InnerRelation guidewire_body_inner(guidewire_body);
        InnerRelation aorta_body_shell_inner(aorta_body);
        using namespace relax_dynamics;
        /** Random reset the insert body particle position. */
        SimpleDynamics<RandomizeParticlePosition> guidewire_random_column_particles(guidewire_body);
        SimpleDynamics<RandomizeParticlePosition> aorta_random_column_particles(aorta_body);
        /** Write the body state to Vtp file. */
        BodyStatesRecordingToVtp guidewire_write_column_to_vtp(io_environment,guidewire_body);
        BodyStatesRecordingToVtp aorta_write_column_to_vtp(io_environment,aorta_body);
        /** Write the particle reload files. */

        ReloadParticleIO guidewire_write_particle_reload_files(io_environment,guidewire_body);
        ReloadParticleIO aorta_write_particle_reload_files(io_environment,aorta_body);
        /** A  Physics relaxation step. */
        RelaxationStepInner guidewire_relaxation_step_inner(guidewire_body_inner);
        RelaxationStepInner aorta_relaxation_step_inner(aorta_body_shell_inner);
        /**
         * @brief 	Particle relaxation starts here.
         */
        guidewire_random_column_particles.exec(0.25);
        guidewire_relaxation_step_inner.SurfaceBounding().exec();
        guidewire_write_column_to_vtp.writeToFile(0.0);

        aorta_random_column_particles.exec(0.25);
        aorta_relaxation_step_inner.SurfaceBounding().exec();
        aorta_write_column_to_vtp.writeToFile(0.0);
        /** relax particles of the insert body. */
        int ite_p = 0;
        while (ite_p < 1000)
        {
            guidewire_relaxation_step_inner.exec();
            aorta_relaxation_step_inner.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the column body N = " << ite_p << "\n";
                guidewire_write_column_to_vtp.writeToFile(ite_p);
                aorta_write_column_to_vtp.writeToFile(ite_p);

            }
        }
        std::cout << "The physics relaxation process of cylinder body finish !" << std::endl;
        /** Output results. */
        guidewire_write_column_to_vtp.writeToFile(ite_p);
        aorta_write_column_to_vtp.writeToFile(ite_p);

        guidewire_write_particle_reload_files.writeToFile(0.0);
        aorta_write_particle_reload_files.writeToFile(0.0);
        return 0.;
    }
    //----------------------------------------------------------------------
    //	All numerical methods will be used in this case.
    //----------------------------------------------------------------------
    // Solid solvers
    SolidAlgorithm aorta_algs(aorta_inner, get_physical_viscosity_general(aorta_rho0_s, aorta_Youngs_modulus, aorta_thickness));
    SolidAlgorithm guidewire_algs(guidewire_inner, get_physical_viscosity_general(guidewire_rho, guidewire_youngs_modulus, res_guidewide));
    // Contact
    SolidContact aorta_contact(aorta_body, guidewire_body);
    SolidContact guidewire_contact(guidewire_body, aorta_body);

    // Boundary conditions
    Vec3d guidewire_insert_velocity = guidewire_insert_velocity_magnitude*guidewire_direction;
    auto guidewide_velocity_condition = [&]()
    {
            auto particles_s = dynamic_cast<BaseParticles *>(&guidewire_body.getBaseParticles());

            // auto& particles_s = guidewire_body.getBaseParticles();
             auto vel = particles_s->getVariableByName<Vec3d>("Velocity");
            particle_for(
                execution::ParallelPolicy(),
                 particles_s->total_real_particles_,
                [&](size_t index_i)
                {
                    if(!aorta_outer_lvl_shape.checkContain(particles_s->pos_[index_i]))
                    {
                        (*vel)[index_i] = guidewire_insert_velocity;
                    }
                }
            );
    };    

    guidewire_body.addBodyStateForRecording<Vec3d>("Velocity");
    aorta_body.addBodyStateForRecording<Vec3d>("Velocity");

    //----------------------------------------------------------------------
    // From here the time stepping begins.
    //----------------------------------------------------------------------
    GlobalStaticVariables::physical_time_ = 0.0;
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    aorta_algs.correct_config();
    guidewire_algs.correct_config();

    //----------------------------------------------------------------------
    // Setup time-stepping related simulation parameters.
    //----------------------------------------------------------------------
    Real end_time = 1.0e-4;
    Real output_period = 1.0e-5; // anyway 50 write_states files in total
    Real dt = 0.0;
    /** Statistics for computing time. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    /**define simple data file input and outputs functions. */
    BodyStatesRecordingToVtp write_states(io_environment,sph_system.real_bodies_);
    write_states.writeToFile();
    //----------------------------------------------------------------------
    // Main time-stepping loop.
    //----------------------------------------------------------------------
    size_t output_iteration = 0;

    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        output_iteration = 0;
        const Real dt_solid_acoustic_ref =
            SMIN(aorta_algs.get_time_step_size(), guidewire_algs.get_time_step_size());
        dt = 0;
        size_t ite = 0;

        GlobalStaticVariables::physical_time_ = 0;
        TickCount t1 = TickCount::now();
        TimeInterval interval;
        while (GlobalStaticVariables::physical_time_ < end_time)
        {
            Real integration_time = 0.0;
            while (integration_time < output_period)
            {
                if (ite % 100 == 0)
                    std::cout << "N=" << ite << " Time: " << GlobalStaticVariables::physical_time_ << "	dt: " << dt
                              << "\n";
                aorta_contact.exec();
                guidewire_contact.exec();
                // time step computation
                dt = SMIN(aorta_algs.get_time_step_size(), guidewire_algs.get_time_step_size());
                { // check for extremely decreased dt
                    if (dt < dt_solid_acoustic_ref / 1e2)
                    {
                        std::cout << "dt = " << dt << " <<< dt_solid_acoustic_ref = " << dt_solid_acoustic_ref << "\n";
                        throw std::runtime_error("The time step decreased too much, stopping simulation\n");
                    }
                }
                // stress and damping
                aorta_algs.exec_all(dt);
                guidewire_algs.exec_all(dt);
                guidewide_velocity_condition();
                // update
                aorta_body.updateCellLinkedList();
                guidewire_body.updateCellLinkedList();
                aorta_contact.update_config();
                guidewire_contact.update_config();

                // timestepping
                ++ite;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
            }

            // output
            ++output_iteration;
            TickCount t2 = TickCount::now();
            write_states.writeToFile(output_iteration);
            TickCount t3 = TickCount::now();
            interval += t3 - t2;
        }
        TickCount t4 = TickCount::now();
        TimeInterval tt;
        tt = t4 - t1;
        std::cout << "Total time for computation: " << tt.seconds() << " seconds." << std::endl;
    };
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    return 0;
}
