#ifndef FLOWRATE_CK_HPP
#define FLOWRATE_CK_HPP
#include "flowrate_ck.h"

namespace SPH
{
namespace fluid_dynamics
{

//========================================================================================
//   1) Implementation of the Base Class
//========================================================================================
template <class BaseInteractionType>
template <class DynamicsIdentifier>
FlowrateCKBase<BaseInteractionType>::
    FlowrateCKBase(DynamicsIdentifier &identifier)
    : BaseInteractionType(identifier),
      dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")),
      dv_dpos_(this->particles_->template getVariableByName<Vecd>("Displacement")),
      dv_vel_(this->particles_->template getVariableByName<Vecd>("Velocity")),
      dv_volumetric_flux_(this->particles_->template registerStateVariableOnly<Vecd>("VelocityContribute"))
{
}
//========================================================================================
//   2) Partial Specialization:
//      <Inner<WithUpdate, KernelCorrectionType, ResolutionType, ExtraParams...>
//========================================================================================
template <class UpdatePolicy, class KernelCorrectionType, class ResolutionType, typename... Parameters>
FlowrateCK<Inner<UpdatePolicy, KernelCorrectionType, ResolutionType, Parameters...>>::
    FlowrateCK(Relation<Inner<Parameters...>> &inner_relation)
    : FlowrateCKBase<Interaction<Inner<Parameters...>>>(inner_relation),
      kernel_correction_(this->particles_)
{
}

//------------------------------------------------------------------------------
// InteractKernel Implementation
//------------------------------------------------------------------------------
template <class UpdatePolicy, class KernelCorrectionType, class ResolutionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
FlowrateCK<Inner<UpdatePolicy, KernelCorrectionType, ResolutionType, Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseInteraction::InteractKernel(ex_policy, encloser),
      correction_(ex_policy, encloser.kernel_correction_),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      dpos_(encloser.dv_dpos_->DelegatedData(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
      volumetric_flux_(encloser.dv_volumetric_flux_->DelegatedData(ex_policy))
{
}

template <class UpdatePolicy, class KernelCorrectionType, class ResolutionType, typename... Parameters>
void FlowrateCK<Inner<UpdatePolicy, KernelCorrectionType, ResolutionType, Parameters...>>::InteractKernel::
    interact(size_t index_i, Real dt)
{
    {
        Vecd volumetric_flux = Vecd::Zero();
        const Vecd vel_i = this->vel_[index_i];

        for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
        {
            UnsignedInt index_j = this->neighbor_index_[n];
            const Real Vj = Vol_[index_j];

            // kernel gradient: ∇W = (dW/dr) * e_ij
            const Real dW_dr = this->dW_ij(index_i, index_j);
            const Vecd e_ij = this->e_ij(index_i, index_j);
            const Vecd gradW = dW_dr * e_ij;

            // neighbor velocity
            const Vecd vel_j = this->vel_[index_j];

            // acceleration for transport velocity
            volumetric_flux += Vj * (vel_j - vel_i).cwiseProduct(gradW) * (this->correction_(index_i) + this->correction_(index_j));
        }
        this->volumetric_flux_[index_i] = volumetric_flux; // still missing * volume_i
    }
}
//------------------------------------------------------------------------------
// UpdateKernel Implementation
//------------------------------------------------------------------------------
template <class UpdatePolicy,
          class KernelCorrectionType,
          class ResolutionType,
          typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
FlowrateCK<Inner<UpdatePolicy, KernelCorrectionType, ResolutionType, Parameters...>>::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      volumetric_flux_(encloser.dv_volumetric_flux_->DelegatedData(ex_policy))
{
}

template <class UpdatePolicy,
          class KernelCorrectionType,
          class ResolutionType,
          typename... ExtraParams>
void FlowrateCK<
    Inner<UpdatePolicy, KernelCorrectionType, ResolutionType, ExtraParams...>>::UpdateKernel::
    update(size_t index_i, Real dt)
{
    volumetric_flux_[index_i] *= vol_[index_i];
}

//========================================================================================
//   Partial Specialization:
//      <Contact<WithUpdate, KernelCorrectionType, ExtraParams...>>
//========================================================================================

template <class KernelCorrectionType, class ResolutionType, typename... Parameters>
FlowrateCK<Contact<Wall, KernelCorrectionType, ResolutionType, Parameters...>>::
    FlowrateCK(Relation<Contact<Parameters...>> &contact_relation)
    : FlowrateCKBase<Interaction<Contact<Wall, Parameters...>>>(contact_relation),
      kernel_correction_(this->particles_)
{
    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        dv_contact_wall_vel_.push_back(this->contact_particles_[k]->template getVariableByName<Vecd>("Velocity"));
        dv_contact_wall_Vol_.push_back(this->contact_particles_[k]->template getVariableByName<Real>("VolumetricMeasure"));
    }
}
//------------------------------------------------------------------------------
// InteractKernel Implementation
//------------------------------------------------------------------------------
template <class KernelCorrectionType, class ResolutionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
FlowrateCK<Contact<Wall, KernelCorrectionType, ResolutionType, Parameters...>>::
    InteractKernel::InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index)
    : BaseInteraction::InteractKernel(ex_policy, encloser, contact_index),
      correction_(ex_policy, encloser.kernel_correction_),
      vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
      contact_wall_vel_(encloser.dv_contact_wall_vel_[contact_index]->DelegatedData(ex_policy)),
      volumetric_flux_(encloser.dv_volumetric_flux_->DelegatedData(ex_policy)),
      contact_wall_Vol_(encloser.dv_contact_wall_Vol_[contact_index]->DelegatedData(ex_policy))
{
}

template <class KernelCorrectionType, class ResolutionType, typename... Parameters>
void FlowrateCK<Contact<Wall, KernelCorrectionType, ResolutionType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Vecd volumetric_flux = Vecd::Zero();
    const Vecd vel_i = this->vel_[index_i];

    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {

        UnsignedInt index_j = this->neighbor_index_[n];
        const Real Vj = contact_wall_Vol_[index_j];

        // kernel gradient: ∇W = (dW/dr) * e_ij
        const Real dW_dr = this->dW_ij(index_i, index_j);
        const Vecd e_ij = this->e_ij(index_i, index_j);
        const Vecd gradW = dW_dr * e_ij;

        // neighbor velocity
        const Vecd vel_j = this->contact_wall_vel_[index_j];

        // acceleration for transport velocity
        volumetric_flux += Vj * (vel_j - vel_i).cwiseProduct(gradW) * (2 * this->correction_(index_i));
    }
    this->volumetric_flux_[index_i] += volumetric_flux;
}

} // end namespace fluid_dynamics
} // end namespace SPH
#endif // FLOWRATE_CK_HPP
