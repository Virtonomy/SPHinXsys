#include "thin_structure_dynamics.h"

namespace SPH
{
//=====================================================================================================//
namespace thin_structure_dynamics
{
//=================================================================================================//
ShellDynamicsInitialCondition::ShellDynamicsInitialCondition(SPHBody &sph_body)
    : LocalDynamics(sph_body), ShellDataSimple(sph_body),
      n0_(particles_->n0_), n_(particles_->n_), pseudo_n_(particles_->pseudo_n_),
      pos0_(particles_->pos0_), transformation_matrix_(particles_->transformation_matrix_) {}
//=================================================================================================//
ShellAcousticTimeStepSize::ShellAcousticTimeStepSize(SPHBody &sph_body, Real CFL)
    : LocalDynamicsReduce<Real, ReduceMin>(sph_body, Real(MaxRealNumber)),
      ShellDataSimple(sph_body), CFL_(CFL), vel_(particles_->vel_), acc_(particles_->acc_),
      angular_vel_(particles_->angular_vel_), dangular_vel_dt_(particles_->dangular_vel_dt_),
      acc_prior_(particles_->acc_prior_),
      thickness_(particles_->thickness_),
      rho0_(particles_->elastic_solid_.ReferenceDensity()),
      E0_(particles_->elastic_solid_.YoungsModulus()),
      nu_(particles_->elastic_solid_.PoissonRatio()),
      c0_(particles_->elastic_solid_.ReferenceSoundSpeed()),
      smoothing_length_(sph_body.sph_adaptation_->ReferenceSmoothingLength()) {}
//=================================================================================================//
Real ShellAcousticTimeStepSize::reduce(size_t index_i, Real dt)
{
    // Since the particle does not change its configuration in pressure relaxation step,
    // I chose a time-step size according to Eulerian method.
    Real time_setp_0 = SMIN((Real)sqrt(smoothing_length_ / ((acc_[index_i] + acc_prior_[index_i]).norm() + TinyReal)),
                            smoothing_length_ / (c0_ + vel_[index_i].norm()));
    Real time_setp_1 = SMIN((Real)sqrt(1.0 / (dangular_vel_dt_[index_i].norm() + TinyReal)),
                            Real(1.0) / (angular_vel_[index_i].norm() + TinyReal));
    Real time_setp_2 = smoothing_length_ * (Real)sqrt(rho0_ * (1.0 - nu_ * nu_) / E0_ /
                                                      (2.0 + (Pi * Pi / 12.0) * (1.0 - nu_) *
                                                                 (1.0 + 1.5 * pow(smoothing_length_ / thickness_[index_i], 2))));
    return CFL_ * SMIN(time_setp_0, time_setp_1, time_setp_2);
}
//=================================================================================================//
ShellCorrectConfiguration::
    ShellCorrectConfiguration(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), ShellDataInner(inner_relation),
      B_(particles_->B_),
      n0_(particles_->n0_), transformation_matrix_(particles_->transformation_matrix_) {}
//=================================================================================================//
ShellDeformationGradientTensor::
    ShellDeformationGradientTensor(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), ShellDataInner(inner_relation),
      pos_(particles_->pos_), pseudo_n_(particles_->pseudo_n_), n0_(particles_->n0_),
      B_(particles_->B_), F_(particles_->F_), F_bending_(particles_->F_bending_),
      transformation_matrix_(particles_->transformation_matrix_) {}
//=================================================================================================//
BaseShellRelaxation::BaseShellRelaxation(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), ShellDataInner(inner_relation),
      rho_(particles_->rho_),
      thickness_(particles_->thickness_),
      pos_(particles_->pos_), vel_(particles_->vel_),
      acc_(particles_->acc_),
      acc_prior_(particles_->acc_prior_),
      n0_(particles_->n0_), pseudo_n_(particles_->pseudo_n_),
      dpseudo_n_dt_(particles_->dpseudo_n_dt_), dpseudo_n_d2t_(particles_->dpseudo_n_d2t_),
      rotation_(particles_->rotation_), angular_vel_(particles_->angular_vel_),
      dangular_vel_dt_(particles_->dangular_vel_dt_),
      B_(particles_->B_), F_(particles_->F_), dF_dt_(particles_->dF_dt_),
      F_bending_(particles_->F_bending_), dF_bending_dt_(particles_->dF_bending_dt_),
      transformation_matrix_(particles_->transformation_matrix_) {}
//=================================================================================================//
ShellStressRelaxationFirstHalf::
    ShellStressRelaxationFirstHalf(BaseInnerRelation &inner_relation,
                                   int number_of_gaussian_points, bool hourglass_control)
    : BaseShellRelaxation(inner_relation),
      elastic_solid_(particles_->elastic_solid_),
      global_stress_(particles_->global_stress_),
      global_moment_(particles_->global_moment_),
      mid_surface_cauchy_stress_(particles_->mid_surface_cauchy_stress_),
      numerical_damping_scaling_(particles_->numerical_damping_scaling_),
      global_shear_stress_(particles_->global_shear_stress_),
      n_(particles_->n_),
      rho0_(elastic_solid_.ReferenceDensity()),
      inv_rho0_(1.0 / rho0_),
      smoothing_length_(sph_body_.sph_adaptation_->ReferenceSmoothingLength()),
      E0_(elastic_solid_.YoungsModulus()),
      G0_(elastic_solid_.ShearModulus()),
      nu_(elastic_solid_.PoissonRatio()),
      hourglass_control_(hourglass_control),
      number_of_gaussian_points_(number_of_gaussian_points)
{
    /** Note that, only three-point and five-point Gaussian quadrature rules are defined. */
    switch (number_of_gaussian_points)
    {
    case 1:
        gaussian_point_ = one_gaussian_point_;
        gaussian_weight_ = one_gaussian_weight_;
        break;
    case 5:
        gaussian_point_ = five_gaussian_points_;
        gaussian_weight_ = five_gaussian_weights_;
        break;
    default:
        gaussian_point_ = three_gaussian_points_;
        gaussian_weight_ = three_gaussian_weights_;
    }
    /** Define the factor of hourglass control algorithm. */
    hourglass_control_factor_ = 0.01;
}
//=================================================================================================//
void ShellStressRelaxationFirstHalf::initialization(size_t index_i, Real dt)
{
    // Note that F_[index_i], F_bending_[index_i], dF_dt_[index_i], dF_bending_dt_[index_i]
    // and rotation_[index_i], angular_vel_[index_i], dangular_vel_dt_[index_i], B_[index_i]
    // are defined in local coordinates, while others in global coordinates.
    pos_[index_i] += vel_[index_i] * dt * 0.5;
    rotation_[index_i] += angular_vel_[index_i] * dt * 0.5;
    pseudo_n_[index_i] += dpseudo_n_dt_[index_i] * dt * 0.5;

    F_[index_i] += dF_dt_[index_i] * dt * 0.5;
    F_bending_[index_i] += dF_bending_dt_[index_i] * dt * 0.5;

    Real J = F_[index_i].determinant();
    Matd inverse_F = F_[index_i].inverse();

    rho_[index_i] = rho0_ / J;

    /** Calculate the current normal direction of mid-surface. */
    n_[index_i] = transformation_matrix_[index_i].transpose() * getNormalFromDeformationGradientTensor(F_[index_i]);
    /** Get transformation matrix from global coordinates to current local coordinates. */
    Matd current_transformation_matrix = getTransformationMatrix(pseudo_n_[index_i]);
    Matd Q = current_transformation_matrix * transformation_matrix_[index_i].transpose();

    Matd resultant_stress = Matd::Zero();
    Matd resultant_moment = Matd::Zero();
    Vecd resultant_shear_stress = Vecd::Zero();

    for (int i = 0; i != number_of_gaussian_points_; ++i)
    {
        Matd F_gaussian_point = F_[index_i] + gaussian_point_[i] * F_bending_[index_i] * thickness_[index_i] * 0.5;
        Matd dF_gaussian_point_dt = dF_dt_[index_i] + gaussian_point_[i] * dF_bending_dt_[index_i] * thickness_[index_i] * 0.5;
        Matd inverse_F_gaussian_point = F_gaussian_point.inverse();
        Matd e_L = 0.5 * (Matd::Identity() - inverse_F_gaussian_point.transpose() * inverse_F_gaussian_point);
        Matd current_local_almansi_strain = Q * e_L * Q.transpose();

        // current_local_almansi_strain(2, 2) = 0;
        Matd cauchy_stress = elastic_solid_.StressCauchy(current_local_almansi_strain, F_gaussian_point, index_i);
        if (!cauchy_stress.allFinite())
        {
            std::string info = "i = " + std::to_string(index_i) + ", Gaussian point = " + std::to_string(i);
            info += "\nF_gaussian_point: " + std::to_string(F_gaussian_point(0, 0)) + " " + std::to_string(F_gaussian_point(0, 1)) + " " + std::to_string(F_gaussian_point(0, 2)) + "\n" +
                    std::to_string(F_gaussian_point(1, 0)) + " " + std::to_string(F_gaussian_point(1, 1)) + " " + std::to_string(F_gaussian_point(1, 2)) + "\n" +
                    std::to_string(F_gaussian_point(2, 0)) + " " + std::to_string(F_gaussian_point(2, 1)) + " " + std::to_string(F_gaussian_point(2, 2)) + "\n";
            info += "\nQ = " + std::to_string(Q(0, 0)) + " " + std::to_string(Q(0, 1)) + " " + std::to_string(Q(0, 2)) + "\n" +
                    std::to_string(Q(1, 0)) + " " + std::to_string(Q(1, 1)) + " " + std::to_string(Q(1, 2)) + "\n" +
                    std::to_string(Q(2, 0)) + " " + std::to_string(Q(2, 1)) + " " + std::to_string(Q(2, 2)) + "\n";
            info += "\ne_L: " + std::to_string(e_L(0, 0)) + " " + std::to_string(e_L(0, 1)) + " " + std::to_string(e_L(0, 2)) + "\n" +
                    std::to_string(e_L(1, 0)) + " " + std::to_string(e_L(1, 1)) + " " + std::to_string(e_L(1, 2)) + "\n" +
                    std::to_string(e_L(2, 0)) + " " + std::to_string(e_L(2, 1)) + " " + std::to_string(e_L(2, 2)) + "\n";
            info += "\ncurrent_local_almansi_strain: " + std::to_string(current_local_almansi_strain(0, 0)) + " " + std::to_string(current_local_almansi_strain(0, 1)) + " " + std::to_string(current_local_almansi_strain(0, 2)) + "\n" +
                    std::to_string(current_local_almansi_strain(1, 0)) + " " + std::to_string(current_local_almansi_strain(1, 1)) + " " + std::to_string(current_local_almansi_strain(1, 2)) + "\n" +
                    std::to_string(current_local_almansi_strain(2, 0)) + " " + std::to_string(current_local_almansi_strain(2, 1)) + " " + std::to_string(current_local_almansi_strain(2, 2)) + "\n";
            info += "\ndet(e_L): " + std::to_string(e_L.determinant()) + "\n";
            info += "det(e_l): " + std::to_string(current_local_almansi_strain.determinant()) + "\n";
            info += "B: " + std::to_string(((-2.0 * current_local_almansi_strain + SPH::Mat3d::Identity()).inverse())(0, 0)) + " " +
                    std::to_string(((-2.0 * current_local_almansi_strain + SPH::Mat3d::Identity()).inverse())(0, 1)) + " " +
                    std::to_string(((-2.0 * current_local_almansi_strain + SPH::Mat3d::Identity()).inverse())(0, 2)) + "\n" +
                    std::to_string(((-2.0 * current_local_almansi_strain + SPH::Mat3d::Identity()).inverse())(1, 0)) + " " +
                    std::to_string(((-2.0 * current_local_almansi_strain + SPH::Mat3d::Identity()).inverse())(1, 1)) + " " +
                    std::to_string(((-2.0 * current_local_almansi_strain + SPH::Mat3d::Identity()).inverse())(1, 2)) + "\n" +
                    std::to_string(((-2.0 * current_local_almansi_strain + SPH::Mat3d::Identity()).inverse())(2, 0)) + " " +
                    std::to_string(((-2.0 * current_local_almansi_strain + SPH::Mat3d::Identity()).inverse())(2, 1)) + " " +
                    std::to_string(((-2.0 * current_local_almansi_strain + SPH::Mat3d::Identity()).inverse())(2, 2)) + "\n";
            throw std::runtime_error(info);
        }
        { /// Enforce plane stress condition adapting algorithm from Sec. 5.4.1 from http://dx.doi.org/10.18419/opus-14215
          /// Differential geometry and the geometrically non-linear Reissner-Mindlin shell model
          /// @WARN Algorithm is not guaranteed to converge, see discussion in Sec. 5.4.1
          ///       even more so considering we do not reuse analytical derivatives.

            // Based on the Cauchy stress and the inverse left Cauchy-Green tensor c = I-2*almansi_strain
            // Initial secant is based on the linear elastic model
            double slope = -0.5 * (elastic_solid_.BulkModulus() + 4.0 / 3.0 * elastic_solid_.ShearModulus());
            int it = 0;
            constexpr int max_iterations = 20; // @WARN hard-coded maximum number of iterations

            std::vector<double> s_history;
            std::vector<double> e_history;
            std::vector<double> E_history;
            s_history.reserve(max_iterations);
            e_history.reserve(max_iterations);
            E_history.reserve(max_iterations);
            s_history.push_back(cauchy_stress(2, 2));
            e_history.push_back(current_local_almansi_strain(2, 2));
            E_history.push_back(slope);

            constexpr auto infinity = std::numeric_limits<Real>::infinity();
            auto tolerance_sqr = [](const Matd &stress)
            {
                return std::max(Eps, Eps * stress.block<2, 2>(0, 0).colwise().squaredNorm().minCoeff());
            };
            for (double s_next = infinity;
                 s_next * s_next > tolerance_sqr(cauchy_stress) && it < max_iterations;
                 ++it)
            {
                double s_prev = cauchy_stress(2, 2);
                Matd c = -2.0 * current_local_almansi_strain + Matd::Identity();
                double c_prev = c(2, 2);
                double d_c = -s_prev / slope;
                double det_c = c.determinant();
                double m_12 = c.block<2, 2>(0, 0).determinant();
                double c_next = c_prev + det_c * std::expm1(d_c * m_12 / det_c) / m_12;
                double e_next = 0.5 * (1.0 - c_next);
                current_local_almansi_strain(2, 2) = e_next;
                cauchy_stress = elastic_solid_.StressCauchy(current_local_almansi_strain, F_gaussian_point, index_i);
                s_next = cauchy_stress(2, 2);
                slope = (s_next - s_prev) / (c_next - c_prev);

                s_history.push_back(s_next);
                e_history.push_back(e_next);
                E_history.push_back(slope);
            }
            if (cauchy_stress.allFinite() == false || it == max_iterations)
            {
                std::string info = "i = " + std::to_string(index_i) + ", history: ";
                for (size_t i = 0; i < s_history.size(); ++i)
                {
                    info += "\n  " + std::to_string(i) + ": s=" + std::to_string(s_history[i]) +
                            ", e=" + std::to_string(e_history[i]) + ", E=" + std::to_string(E_history[i]);
                }
                for (size_t m = 0; m < 3; ++m)
                {
                    for (size_t n = 0; n < 3; ++n)
                    {
                        info += "\n  almansi strain(" + std::to_string(m) + ", " + std::to_string(n) + ")=" + std::to_string(current_local_almansi_strain(m, n));
                    }
                }
                throw std::runtime_error(info);
            }
        }
        /// Impact of including numerical damping in the algorithm above unclear
        /// Left here in absence of discriminating factors
        Matd damping = current_transformation_matrix * transformation_matrix_[index_i].transpose() * F_gaussian_point * elastic_solid_.NumericalDampingRightCauchy(F_gaussian_point, dF_gaussian_point_dt, numerical_damping_scaling_[index_i], index_i) * F_gaussian_point.transpose() * transformation_matrix_[index_i] * current_transformation_matrix.transpose() / F_gaussian_point.determinant();
        cauchy_stress += damping;

        /** Impose modeling assumptions. */
        cauchy_stress.col(Dimensions - 1) *= shear_correction_factor_;
        cauchy_stress.row(Dimensions - 1) *= shear_correction_factor_;
        cauchy_stress(Dimensions - 1, Dimensions - 1) = 0.0;

        if (i == 0)
        {
            mid_surface_cauchy_stress_[index_i] = cauchy_stress;
        }

        /** Integrate Cauchy stress along thickness. */
        resultant_stress +=
            0.5 * thickness_[index_i] * gaussian_weight_[i] * cauchy_stress;
        resultant_moment +=
            0.5 * thickness_[index_i] * gaussian_weight_[i] * (cauchy_stress * gaussian_point_[i] * thickness_[index_i] * 0.5);
        resultant_shear_stress -=
            0.5 * thickness_[index_i] * gaussian_weight_[i] * cauchy_stress.col(Dimensions - 1);

        resultant_stress.col(Dimensions - 1) = Vecd::Zero();
        resultant_moment.col(Dimensions - 1) = Vecd::Zero();
    }

    /** stress and moment in global coordinates for pair interaction */
    global_stress_[index_i] = J * current_transformation_matrix.transpose() * resultant_stress * current_transformation_matrix * transformation_matrix_[index_i].transpose() * inverse_F.transpose() * transformation_matrix_[index_i];
    global_moment_[index_i] = J * current_transformation_matrix.transpose() * resultant_moment * current_transformation_matrix * transformation_matrix_[index_i].transpose() * inverse_F.transpose() * transformation_matrix_[index_i];
    global_shear_stress_[index_i] = J * current_transformation_matrix.transpose() * resultant_shear_stress;
}
//=================================================================================================//
void ShellStressRelaxationFirstHalf::update(size_t index_i, Real dt)
{
    vel_[index_i] += (acc_prior_[index_i] + acc_[index_i]) * dt;
    angular_vel_[index_i] += dangular_vel_dt_[index_i] * dt;
}
//=================================================================================================//
void ShellStressRelaxationSecondHalf::initialization(size_t index_i, Real dt)
{
    pos_[index_i] += vel_[index_i] * dt * 0.5;
    rotation_[index_i] += angular_vel_[index_i] * dt * 0.5;
    dpseudo_n_dt_[index_i] = transformation_matrix_[index_i].transpose() *
                             getVectorChangeRateAfterThinStructureRotation(local_pseudo_n_0, rotation_[index_i], angular_vel_[index_i]);
    pseudo_n_[index_i] += dpseudo_n_dt_[index_i] * dt * 0.5;
}
//=================================================================================================//
void ShellStressRelaxationSecondHalf::update(size_t index_i, Real dt)
{
    F_[index_i] += dF_dt_[index_i] * dt * 0.5;
    F_bending_[index_i] += dF_bending_dt_[index_i] * dt * 0.5;
}
//=================================================================================================//
ConstrainShellBodyRegion::
    ConstrainShellBodyRegion(BodyPartByParticle &body_part)
    : BaseLocalDynamics<BodyPartByParticle>(body_part), ShellDataSimple(sph_body_),
      vel_(particles_->vel_), angular_vel_(particles_->angular_vel_) {}
//=================================================================================================//
void ConstrainShellBodyRegion::update(size_t index_i, Real dt)
{
    vel_[index_i] = Vecd::Zero();
    angular_vel_[index_i] = Vecd::Zero();
}
//=================================================================================================//
ConstrainShellBodyRegionAlongAxis::ConstrainShellBodyRegionAlongAxis(BodyPartByParticle &body_part, int axis)
    : BaseLocalDynamics<BodyPartByParticle>(body_part), ShellDataSimple(sph_body_),
      axis_(axis), pos_(particles_->pos_), pos0_(particles_->pos0_), vel_(particles_->vel_),
      acc_(particles_->acc_), rotation_(particles_->rotation_), angular_vel_(particles_->angular_vel_),
      dangular_vel_dt_(particles_->dangular_vel_dt_) {}
//=================================================================================================//
void ConstrainShellBodyRegionAlongAxis::update(size_t index_i, Real dt)
{
    vel_[index_i][axis_] = 0.0;
    vel_[index_i][2] = 0.0;
    acc_[index_i][axis_] = 0.0;
    acc_[index_i][2] = 0.0;

    angular_vel_[index_i][1 - axis_] = 0.0;
    dangular_vel_dt_[index_i][1 - axis_] = 0.0;
}
//=================================================================================================//
DistributingPointForcesToShell::
    DistributingPointForcesToShell(SPHBody &sph_body, std::vector<Vecd> point_forces,
                                   std::vector<Vecd> reference_positions, Real time_to_full_external_force,
                                   Real particle_spacing_ref, Real h_spacing_ratio)
    : LocalDynamics(sph_body), ShellDataSimple(sph_body),
      point_forces_(point_forces), reference_positions_(reference_positions),
      time_to_full_external_force_(time_to_full_external_force),
      particle_spacing_ref_(particle_spacing_ref), h_spacing_ratio_(h_spacing_ratio),
      pos0_(particles_->pos0_), acc_prior_(particles_->acc_prior_),
      thickness_(particles_->thickness_)
{
    for (size_t i = 0; i < point_forces_.size(); i++)
    {
        weight_.push_back(StdLargeVec<Real>(0.0));
        time_dependent_point_forces_.push_back(Vecd::Zero());
        sum_of_weight_.push_back(0.0);
        particles_->registerVariable(weight_[i], "Weight_" + std::to_string(i));
    }

    getWeight(); // TODO: should be revised and parallelized, using SimpleDynamics
}
//=================================================================================================//
void DistributingPointForcesToShell::getWeight()
{
    Kernel *kernel_ = sph_body_.sph_adaptation_->getKernel();
    Real reference_smoothing_length = sph_body_.sph_adaptation_->ReferenceSmoothingLength();
    Real smoothing_length = h_spacing_ratio_ * particle_spacing_ref_;
    Real h_ratio = reference_smoothing_length / smoothing_length;
    Real cutoff_radius_sqr = pow(2.0 * smoothing_length, 2);
    for (size_t i = 0; i < point_forces_.size(); ++i)
    {
        sum_of_weight_[i] = 0.0;
        for (size_t index = 0; index < particles_->total_real_particles_; ++index)
        {
            weight_[i][index] = 0.0;
            Vecd displacement = reference_positions_[i] - pos0_[index];
            if (displacement.squaredNorm() <= cutoff_radius_sqr)
            {
                weight_[i][index] = kernel_->W(h_ratio, displacement.norm(), displacement);
                sum_of_weight_[i] += weight_[i][index];
            }
        }
    }
}
//=================================================================================================//
void DistributingPointForcesToShell::setupDynamics(Real dt)
{
    Real current_time = GlobalStaticVariables::physical_time_;
    for (size_t i = 0; i < point_forces_.size(); ++i)
    {
        time_dependent_point_forces_[i] = current_time < time_to_full_external_force_
                                              ? current_time * point_forces_[i] / time_to_full_external_force_
                                              : point_forces_[i];
    }
}
//=================================================================================================//
void DistributingPointForcesToShell::update(size_t index_i, Real dt)
{
    acc_prior_[index_i] = Vecd::Zero();
    for (size_t i = 0; i < point_forces_.size(); ++i)
    {
        Vecd force = weight_[i][index_i] / (sum_of_weight_[i] + TinyReal) * time_dependent_point_forces_[i];
        acc_prior_[index_i] += force / particles_->ParticleMass(index_i);
    }
}
//=================================================================================================//
} // namespace thin_structure_dynamics
} // namespace SPH