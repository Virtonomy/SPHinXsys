#include "non_newtonian_dynamics.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
Oldroyd_BIntegration1stHalf<Inner<>>::
    Oldroyd_BIntegration1stHalf(BaseInnerRelation &inner_relation)
    : Integration1stHalfInnerRiemann(inner_relation)
{
    particles_->registerVariable(tau_, "ElasticStress");
    particles_->registerVariable(dtau_dt_, "ElasticStressChangeRate");
    particles_->addVariableToSort<Matd>("ElasticStress");
    particles_->addVariableToRestart<Matd>("ElasticStress");
}
//=================================================================================================//
void Oldroyd_BIntegration1stHalf<Inner<>>::initialization(size_t index_i, Real dt)
{
    Integration1stHalfInnerRiemann::initialization(index_i, dt);

    tau_[index_i] += dtau_dt_[index_i] * dt * 0.5;
}
//=================================================================================================//
void Oldroyd_BIntegration1stHalf<Inner<>>::interaction(size_t index_i, Real dt)
{
    Integration1stHalfInnerRiemann::interaction(index_i, dt);

    Vecd force = Vecd::Zero();
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Vecd nablaW_ijV_j = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];

        // elastic force
        force += mass_[index_i] * (tau_[index_i] + tau_[index_j]) * nablaW_ijV_j;
    }

    force_[index_i] += force / rho_[index_i];
}
//=================================================================================================//
Oldroyd_BIntegration1stHalf<Contact<Wall>>::
    Oldroyd_BIntegration1stHalf(BaseContactRelation &wall_contact_relation)
    : Integration1stHalfContactWallRiemann(wall_contact_relation),
      tau_(*particles_->getVariableByName<Matd>("ElasticStress")){};
//=================================================================================================//
void Oldroyd_BIntegration1stHalf<Contact<Wall>>::interaction(size_t index_i, Real dt)
{
    Integration1stHalfContactWallRiemann::interaction(index_i, dt);

    Real rho_i = rho_[index_i];
    Matd tau_i = tau_[index_i];

    Vecd force = Vecd::Zero();
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Neighborhood &wall_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
        {
            Vecd nablaW_ijV_j = wall_neighborhood.dW_ijV_j_[n] * wall_neighborhood.e_ij_[n];
            /** stress boundary condition. */
            force += mass_[index_i] * 2.0 * tau_i * nablaW_ijV_j / rho_i;
        }
    }

    force_[index_i] += force;
}
//=================================================================================================//
Oldroyd_BIntegration2ndHalf<Inner<>>::
    Oldroyd_BIntegration2ndHalf(BaseInnerRelation &inner_relation)
    : Integration2ndHalfInnerRiemann(inner_relation),
      oldroyd_b_fluid_(DynamicCast<Oldroyd_B_Fluid>(this, particles_->getBaseMaterial())),
      vel_grad_(*particles_->getVariableByName<Matd>("VelocityGradient")),
      tau_(*particles_->getVariableByName<Matd>("ElasticStress")),
      dtau_dt_(*particles_->getVariableByName<Matd>("ElasticStressChangeRate"))
{
    mu_p_ = oldroyd_b_fluid_.ReferencePolymericViscosity();
    lambda_ = oldroyd_b_fluid_.getReferenceRelaxationTime();
}
//=================================================================================================//
void Oldroyd_BIntegration2ndHalf<Inner<>>::update(size_t index_i, Real dt)
{
    Integration2ndHalfInnerRiemann::update(index_i, dt);

    Matd vel_grad_transpose = vel_grad_[index_i].transpose();
    dtau_dt_[index_i] = vel_grad_transpose * tau_[index_i] + tau_[index_i] * vel_grad_[index_i] -
                        tau_[index_i] / lambda_ +
                        (vel_grad_transpose + vel_grad_[index_i]) * mu_p_ / lambda_;
    tau_[index_i] += dtau_dt_[index_i] * dt * 0.5;
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
