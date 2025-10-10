#pragma once

#include "fluid_dynamics_inner_correction.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class RiemannSolverType>
BaseIntegration1stHalfCorrect<RiemannSolverType>::BaseIntegration1stHalfCorrect(BaseInnerRelation &inner_relation)
    : BaseIntegration1stHalf<RiemannSolverType>(inner_relation),
      B_(*this->particles_->template registerSharedVariable<Matd>("CorrectionMatrix", Matd::Identity()))
{
    this->particles_->registerVariable(p_B_, "CorrectedPressure");
}
//=================================================================================================//
template <class RiemannSolverType>
void BaseIntegration1stHalfCorrect<RiemannSolverType>::initialization(size_t index_i, Real dt)
{
    BaseIntegration1stHalf<RiemannSolverType>::initialization(index_i, dt);
    p_B_[index_i] = this->p_[index_i] * this->B_[index_i];
}
//=================================================================================================//
template <class RiemannSolverType>
void BaseIntegration1stHalfCorrect<RiemannSolverType>::interaction(size_t index_i, Real dt)
{
    Vecd acceleration = Vecd::Zero();
    Real rho_dissipation(0);
    const Neighborhood &inner_neighborhood = this->inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
        const Vecd &e_ij = inner_neighborhood.e_ij_[n];

        acceleration -= (p_B_[index_i] + p_B_[index_j]) * dW_ijV_j * e_ij;
        rho_dissipation += this->riemann_solver_.DissipativeUJump(this->p_[index_i] - this->p_[index_j]) * dW_ijV_j;
    }
    this->acc_[index_i] += acceleration / this->rho_[index_i];
    this->drho_dt_[index_i] = rho_dissipation * this->rho_[index_i];
}
//=================================================================================================//
template <class RiemannSolverType>
BaseIntegration2ndHalfCorrect<RiemannSolverType>::BaseIntegration2ndHalfCorrect(BaseInnerRelation &inner_relation)
    : BaseIntegration2ndHalf<RiemannSolverType>(inner_relation),
      B_(*this->particles_->template registerSharedVariable<Matd>("CorrectionMatrix", Matd::Identity())) {}
//=================================================================================================//
template <class RiemannSolverType>
void BaseIntegration2ndHalfCorrect<RiemannSolverType>::interaction(size_t index_i, Real dt)
{
    Real density_change_rate(0);
    Vecd p_dissipation = Vecd::Zero();
    const Neighborhood &inner_neighborhood = this->inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        const Vecd e_ij_corrected = B_[index_i] * inner_neighborhood.e_ij_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];

        Real u_jump = (this->vel_[index_i] - this->vel_[index_j]).dot(e_ij_corrected);
        density_change_rate += u_jump * dW_ijV_j;
        p_dissipation += this->riemann_solver_.DissipativePJump(u_jump) * dW_ijV_j * e_ij_corrected;
    }
    this->drho_dt_[index_i] += density_change_rate * this->rho_[index_i];
    this->acc_[index_i] = p_dissipation / this->rho_[index_i];
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
