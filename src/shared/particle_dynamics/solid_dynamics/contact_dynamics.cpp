#include "contact_dynamics.h"

#ifdef max
#undef max
#endif

namespace SPH
{
//=========================================================================================================//
namespace solid_dynamics
{
//=================================================================================================//
SelfRepulsionFactorSummation::
    SelfRepulsionFactorSummation(SelfSurfaceContactRelation &self_contact_relation)
    : RepulsionFactorAccessor(self_contact_relation.base_particles_, "SelfRepulsionFactor"),
      LocalDynamics(self_contact_relation.getSPHBody()),
      SolidDataInner(self_contact_relation),
      Vol_(particles_->Vol_)
{
    Real dp_1 = self_contact_relation.getSPHBody().sph_adaptation_->ReferenceSpacing();
    offset_W_ij_ = self_contact_relation.getSPHBody().sph_adaptation_->getKernel()->W(dp_1, ZeroVecd);
}
//=================================================================================================//
RepulsionFactorSummation::
    RepulsionFactorSummation(SurfaceContactRelation &solid_body_contact_relation)
    : RepulsionFactorAccessor(solid_body_contact_relation.base_particles_, "RepulsionFactor"),
      LocalDynamics(solid_body_contact_relation.getSPHBody()),
      ContactDynamicsData(solid_body_contact_relation), Vol_(particles_->Vol_),
      offset_W_ij_(StdVec<Real>(contact_configuration_.size(), 0.0))
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_Vol_.push_back(&(contact_particles_[k]->Vol_));
    }

    // we modify the default formulation by an offset, so that exactly touching bodies produce 0 initial force
    // subtract summation of the kernel function of 2 particles at 1 particle distance, and if the result is negative, we take 0
    // different resolution: distance = 0.5 * dp1 + 0.5 * dp2
    // dp1, dp2 half reference spacing
    Real dp_1 = solid_body_contact_relation.getSPHBody().sph_adaptation_->ReferenceSpacing();
    // different resolution: distance = 0.5 * dp1 + 0.5 * dp2
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Real dp_2 = solid_body_contact_relation.contact_bodies_[k]->sph_adaptation_->ReferenceSpacing();
        Real distance = 0.5 * dp_1 + 0.5 * dp_2;
        offset_W_ij_[k] = solid_body_contact_relation.getSPHBody().sph_adaptation_->getKernel()->W(distance, ZeroVecd);
    }
}
//=================================================================================================//
ShellRepulsionFactor::ShellRepulsionFactor(SurfaceContactRelation &solid_body_contact_relation)
    : RepulsionFactorAccessor(solid_body_contact_relation.base_particles_, "RepulsionFactor"),
      LocalDynamics(solid_body_contact_relation.getSPHBody()),
      ContactDynamicsData(solid_body_contact_relation), solid_(particles_->solid_),
      kernel_(solid_body_contact_relation.getSPHBody().sph_adaptation_->getKernel()),
      particle_spacing_(solid_body_contact_relation.getSPHBody().sph_adaptation_->ReferenceSpacing())
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        Real dp_k = solid_body_contact_relation.contact_bodies_[k]->sph_adaptation_->ReferenceSpacing();
        Real average_spacing_k = 0.5 * particle_spacing_ + 0.5 * dp_k;
        Real h_ratio_k = particle_spacing_ / average_spacing_k;
        offset_W_ij_.push_back(kernel_->W(h_ratio_k, average_spacing_k, ZeroVecd));

        Real contact_max(0.0);
        for (int l = 0; l != 3; ++l)
        {
            Real temp = three_gaussian_points_[l] * average_spacing_k * 0.5 + average_spacing_k * 0.5;
            Real contact_temp = 2.0 * (kernel_->W(h_ratio_k, temp, ZeroVecd) - offset_W_ij_[k]) *
                                average_spacing_k * 0.5 * three_gaussian_weights_[l];
            contact_max += Dimensions == 2 ? contact_temp : contact_temp * Pi * temp;
        }
        /** a calibration factor to avoid particle penetration into shell structure */
        calibration_factor_.push_back(1.0 / (contact_max + Eps));

        contact_Vol_.push_back(&(contact_particles_[k]->Vol_));
    }
}
//=================================================================================================//
SelfContactForce::
    SelfContactForce(SelfSurfaceContactRelation &self_contact_relation)
    : LocalDynamics(self_contact_relation.getSPHBody()),
      SolidDataInner(self_contact_relation),
      solid_(particles_->solid_), mass_(particles_->mass_),
      self_repulsion_factor_(*particles_->getVariableByName<Real>("SelfRepulsionFactor")),
      Vol_(particles_->Vol_), acc_prior_(particles_->acc_prior_),
      vel_(particles_->vel_),
      contact_impedance_(solid_.ReferenceDensity() * sqrt(solid_.ContactStiffness())) {}
//=================================================================================================//
ContactForce::ContactForce(SurfaceContactRelation &solid_body_contact_relation)
    : LocalDynamics(solid_body_contact_relation.getSPHBody()),
      ContactDynamicsData(solid_body_contact_relation),
      solid_(particles_->solid_),
      repulsion_factor_(*particles_->getVariableByName<Real>("RepulsionFactor")),
      Vol_(particles_->Vol_), mass_(particles_->mass_),
      acc_prior_(particles_->acc_prior_)
{
    // The Hertzian theory of non-adhesive elastic contact between two balls gives the analytical solution F = 4/3* E* sqrt(R)* delta^(3/2),
    // where the composite Young's modulus of elasticity 1/E* = (1-v1^2)/E1 + (1-v2^2)/E2, the effective radius 1/R = 1/R1 + 1/R2,
    // Ref: https://en.wikipedia.org/wiki/Contact_mechanics#Hertzian_theory_of_non-adhesive_elastic_contact
    // The pinball contact algorithm also defines the contact force with a similar formula F2 = Gi * Gj / (Gi + Gj) * sqrt(Ri * Rj/(Ri+Rj)) * p^{3/2}
    // Ref: https://doi.org/10.1016/0045-7825(93)90064-5, The splitting pinball method for contact-impact problems
    // Inspired by the composite modulus in these formulas, we define the contact stiffness as the harmonic average K = 2 * Ki * Kj / (Ki + Kj), where K1, K2 are the contact stiffness of the two bodies.
    // In comparison of geometric average, the harmonic average is dominant by the softer material.
    // Empirically, the harmonic average is sufficient to prevent penetration, and matches the time-step size of the softer material.
    // This allows us to use a different time-step size for the two materials
    Real K_1 = solid_.ContactStiffness();
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_solids_.push_back(&contact_particles_[k]->solid_);
        contact_repulsion_factor_.push_back(contact_particles_[k]->getVariableByName<Real>("RepulsionFactor"));
        Real K_2 = contact_solids_[k]->ContactStiffness();
        contact_stiffness_.emplace_back(2 * K_1 * K_2 / (K_1 + K_2));
    }
}
//=================================================================================================//
ContactForceFromWall::ContactForceFromWall(SurfaceContactRelation &solid_body_contact_relation)
    : LocalDynamics(solid_body_contact_relation.getSPHBody()),
      ContactWithWallData(solid_body_contact_relation), solid_(particles_->solid_),
      repulsion_factor_(*particles_->getVariableByName<Real>("RepulsionFactor")),
      Vol_(particles_->Vol_), mass_(particles_->mass_),
      acc_prior_(particles_->acc_prior_) {}
//=================================================================================================//
ContactForceToWall::ContactForceToWall(SurfaceContactRelation &solid_body_contact_relation)
    : LocalDynamics(solid_body_contact_relation.getSPHBody()),
      ContactDynamicsData(solid_body_contact_relation),
      Vol_(particles_->Vol_), mass_(particles_->mass_),
      acc_prior_(particles_->acc_prior_)
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_solids_.push_back(&contact_particles_[k]->solid_);
        contact_repulsion_factor_.push_back(contact_particles_[k]->getVariableByName<Real>("RepulsionFactor"));
    }
}
//=================================================================================================//
PairwiseFrictionFromWall::
    PairwiseFrictionFromWall(BaseContactRelation &contact_relation, Real eta)
    : LocalDynamics(contact_relation.getSPHBody()), ContactWithWallData(contact_relation),
      eta_(eta), Vol_(particles_->Vol_), mass_(particles_->mass_),
      vel_(particles_->vel_)
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        wall_vel_n_.push_back(&contact_particles_[k]->vel_);
        wall_n_.push_back(&contact_particles_[k]->n_);
    }
}
//=================================================================================================//
DynamicContactForceWithWall::
    DynamicContactForceWithWall(SurfaceContactRelation &solid_body_contact_relation, Real penalty_strength)
    : LocalDynamics(solid_body_contact_relation.getSPHBody()),
      ContactDynamicsData(solid_body_contact_relation),
      solid_(particles_->solid_),
      Vol_(particles_->Vol_), mass_(particles_->mass_),
      vel_(particles_->vel_), acc_prior_(particles_->acc_prior_),
      penalty_strength_(penalty_strength)
{
    impedance_ = solid_.ReferenceDensity() * sqrt(solid_.ContactStiffness());
    reference_pressure_ = solid_.ReferenceDensity() * solid_.ContactStiffness();
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_vel_.push_back(&(contact_particles_[k]->vel_));
        contact_n_.push_back(&(contact_particles_[k]->n_));
    }
}
//=================================================================================================//
} // namespace solid_dynamics
} // namespace SPH
