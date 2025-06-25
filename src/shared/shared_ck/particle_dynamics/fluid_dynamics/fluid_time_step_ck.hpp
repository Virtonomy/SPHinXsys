#ifndef VIRTOSIM_FLUID_TIME_STEP_CK_HPP_D55BD858_6531_4A31_A05A_144EFCE2CBA3
#define VIRTOSIM_FLUID_TIME_STEP_CK_HPP_D55BD858_6531_4A31_A05A_144EFCE2CBA3

#include "fluid_time_step_ck.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class FluidType>
AcousticTimeStepCK<FluidType>::AcousticTimeStepCK(SPHBody &sph_body, Real acousticCFL)
    : LocalDynamicsReduce<ReduceMax>(sph_body),
      fluid_(DynamicCast<FluidType>(this, particles_->getBaseMaterial())),
      dv_rho_(particles_->getVariableByName<Real>("Density")),
      dv_p_(particles_->getVariableByName<Real>("Pressure")),
      dv_mass_(particles_->getVariableByName<Real>("Mass")),
      dv_vel_(particles_->getVariableByName<Vecd>("Velocity")),
      dv_force_(particles_->getVariableByName<Vecd>("Force")),
      dv_force_prior_(particles_->getVariableByName<Vecd>("ForcePrior")),
      h_min_(sph_body.getSPHAdaptation().MinimumSmoothingLength()),
      acousticCFL_(acousticCFL) {}
//=================================================================================================//
template <class FluidType>
AcousticTimeStepCK<FluidType>::FinishDynamics::
    FinishDynamics(AcousticTimeStepCK<FluidType> &encloser)
    : h_min_(encloser.h_min_), acousticCFL_(encloser.acousticCFL_) {}
//=================================================================================================//
template <class FluidType>
Real AcousticTimeStepCK<FluidType>::FinishDynamics::Result(Real reduced_value)
{
    // since the particle does not change its configuration in the acoustic time steps
    // I chose a time-step size according to Eulerian method
    return acousticCFL_ * h_min_ / (reduced_value + TinyReal);
}
//=================================================================================================//
template <class FluidType>
template <class ExecutionPolicy>
AcousticTimeStepCK<FluidType>::ReduceKernel::ReduceKernel(
    const ExecutionPolicy &ex_policy, AcousticTimeStepCK<FluidType> &encloser)
    : eos_(encloser.fluid_),
      rho_(encloser.dv_rho_->DelegatedData(ex_policy)),
      p_(encloser.dv_p_->DelegatedData(ex_policy)),
      mass_(encloser.dv_mass_->DelegatedData(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
      force_(encloser.dv_force_->DelegatedData(ex_policy)),
      force_prior_(encloser.dv_force_prior_->DelegatedData(ex_policy)),
      h_min_(encloser.h_min_) {}
//=================================================================================================//
template <class ExecutionPolicy>
AdvectionStepSetup::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, AdvectionStepSetup &encloser)
    : Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      mass_(encloser.dv_mass_->DelegatedData(ex_policy)),
      rho_(encloser.dv_rho_->DelegatedData(ex_policy)),
      dpos_(encloser.dv_dpos_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class ExecutionPolicy>
AdvectionStepClose::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, AdvectionStepClose &encloser)
    : pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      dpos_(encloser.dv_dpos_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class ParticleScopeType, class FluidType>
AcousticTimeStepCK_v2<ParticleScopeType, FluidType>::AcousticTimeStepCK_v2(SPHBody &sph_body, Real acousticCFL)
    : LocalDynamicsReduce<ReduceMax>(sph_body),
      fluid_(DynamicCast<FluidType>(this, particles_->getBaseMaterial())),
      dv_rho_(particles_->getVariableByName<Real>("Density")),
      dv_p_(particles_->getVariableByName<Real>("Pressure")),
      dv_mass_(particles_->getVariableByName<Real>("Mass")),
      dv_vel_(particles_->getVariableByName<Vecd>("Velocity")),
      dv_force_(particles_->getVariableByName<Vecd>("Force")),
      dv_force_prior_(particles_->getVariableByName<Vecd>("ForcePrior")),
      h_min_(sph_body.getSPHAdaptation().MinimumSmoothingLength()),
      acousticCFL_(acousticCFL),
      within_scope_method_(this->particles_) {}
//=================================================================================================//
template <class ParticleScopeType, class FluidType>
AcousticTimeStepCK_v2<ParticleScopeType, FluidType>::FinishDynamics::
    FinishDynamics(AcousticTimeStepCK_v2<ParticleScopeType, FluidType> &encloser)
    : h_min_(encloser.h_min_), acousticCFL_(encloser.acousticCFL_) {}
//=================================================================================================//
template <class ParticleScopeType, class FluidType>
Real AcousticTimeStepCK_v2<ParticleScopeType, FluidType>::FinishDynamics::Result(Real reduced_value)
{
    // since the particle does not change its configuration in the acoustic time steps
    // I chose a time-step size according to Eulerian method
    return acousticCFL_ * h_min_ / (reduced_value + TinyReal);
}
//=================================================================================================//
template <class ParticleScopeType, class FluidType>
template <class ExecutionPolicy>
AcousticTimeStepCK_v2<ParticleScopeType, FluidType>::ReduceKernel::ReduceKernel(
    const ExecutionPolicy &ex_policy, AcousticTimeStepCK_v2<ParticleScopeType, FluidType> &encloser)
    : eos_(encloser.fluid_),
      rho_(encloser.dv_rho_->DelegatedData(ex_policy)),
      p_(encloser.dv_p_->DelegatedData(ex_policy)),
      mass_(encloser.dv_mass_->DelegatedData(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
      force_(encloser.dv_force_->DelegatedData(ex_policy)),
      force_prior_(encloser.dv_force_prior_->DelegatedData(ex_policy)),
      h_min_(encloser.h_min_),
      within_scope_(ex_policy, encloser.within_scope_method_, *this) {}
//=================================================================================================//
template <class ParticleScopeType>
AdvectionTimeStepCK_v2<ParticleScopeType>::
    AdvectionTimeStepCK_v2(SPHBody &sph_body, Real U_ref, Real advectionCFL)
    : LocalDynamicsReduce<ReduceMax>(sph_body),
      h_min_(sph_body.getSPHAdaptation().MinimumSmoothingLength()),
      speed_ref_(U_ref), advectionCFL_(advectionCFL),
      dv_mass_(particles_->getVariableByName<Real>("Mass")),
      dv_vel_(particles_->getVariableByName<Vecd>("Velocity")),
      dv_force_(particles_->getVariableByName<Vecd>("Force")),
      dv_force_prior_(particles_->getVariableByName<Vecd>("ForcePrior")),
      within_scope_method_(this->particles_) {}
//=================================================================================================//
template <class ParticleScopeType>
AdvectionTimeStepCK_v2<ParticleScopeType>::FinishDynamics::FinishDynamics(AdvectionTimeStepCK_v2<ParticleScopeType> &encloser)
    : h_min_(encloser.h_min_), speed_ref_(encloser.speed_ref_),
      advectionCFL_(encloser.advectionCFL_) {}
//=================================================================================================//
template <class ParticleScopeType>
Real AdvectionTimeStepCK_v2<ParticleScopeType>::FinishDynamics::Result(Real reduced_value)
{
    return advectionCFL_ * h_min_ / (SMAX(std::sqrt(reduced_value), speed_ref_) + TinyReal);
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
#endif // VIRTOSIM_FLUID_TIME_STEP_CK_HPP_D55BD858_6531_4A31_A05A_144EFCE2CBA3
