#ifndef VIRTOSIM_DENSITY_REGULARIZATION_20COPY_HPP_CEB4A1AB_A372_4275_BCE5_A03E8572E18C
#define VIRTOSIM_DENSITY_REGULARIZATION_20COPY_HPP_CEB4A1AB_A372_4275_BCE5_A03E8572E18C
#ifndef VIRTOSIM_DENSITY_REGULARIZATION_HPP_D6998CB9_73EC_414F_9E1A_5FAE9F4CBFDC
#define VIRTOSIM_DENSITY_REGULARIZATION_HPP_D6998CB9_73EC_414F_9E1A_5FAE9F4CBFDC

#include "renormalizationVW.h"

#include "base_particles.hpp"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <template <typename...> class RelationType, typename... Parameters>
template <class DynamicsIdentifier>
RenormalizationVW<Base, RelationType<Parameters...>>::
    RenormalizationVW(DynamicsIdentifier &identifier)
    : Interaction<RelationType<Parameters...>>(identifier),
      dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")),
      dv_renormalizationVW_(this->particles_->template registerStateVariableOnly<Real>("ReNormalizationVW"))
{
}

//=================================================================================================//
template <typename RegularizationType, typename ParticleScopeType, typename... Parameters>
RenormalizationVW<Inner<InteractionOnly, RegularizationType, ParticleScopeType, Parameters...>>::
    RenormalizationVW(Relation<Inner<Parameters...>> &inner_relation)
    : RenormalizationVW<Base, Inner<Parameters...>>(inner_relation),
      within_scope_method_(this->particles_)
{
}
//=================================================================================================//
template <typename RegularizationType, typename ParticleScopeType, typename... Parameters>
template <class ExecutionPolicy>
RenormalizationVW<Inner<InteractionOnly, RegularizationType, ParticleScopeType, Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy,
                   RenormalizationVW<Inner<InteractionOnly, RegularizationType, ParticleScopeType, Parameters...>> &encloser)
    : RenormalizationVW<Base, Inner<Parameters...>>::InteractKernel(ex_policy, encloser), Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      renormalizationVW_(encloser.dv_renormalizationVW_->DelegatedData(ex_policy)), W0_(this->kernel_.W(ZeroData<Vecd>::value))
{
}
//=================================================================================================//
template <typename RegularizationType, typename ParticleScopeType, typename... Parameters>
void RenormalizationVW<Inner<InteractionOnly, RegularizationType, ParticleScopeType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Real renormalize_factor = W0_ * this->Vol_[index_i];
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {

        renormalize_factor += this->W_ij(index_i, this->neighbor_index_[n]) * this->Vol_[index_i];
    }

    this->renormalizationVW_[index_i] = renormalize_factor;
}

//=================================================================================================//
template <typename... Parameters>
RenormalizationVW<Contact<Parameters...>>::
    RenormalizationVW(Relation<Contact<Parameters...>> &contact_relation)
    : RenormalizationVW<Base, Contact<Parameters...>>(contact_relation)
{
    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        Real rho0_k = this->contact_bodies_[k]->getBaseMaterial().ReferenceDensity();
        contact_inv_rho0_.push_back(1.0 / rho0_k);
        dv_contact_mass_.push_back(this->contact_particles_[k]->template getVariableByName<Real>("Mass"));
    }
}
//=================================================================================================//
template <typename... Parameters>
template <class ExecutionPolicy>
RenormalizationVW<Contact<Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy,
                   RenormalizationVW<Contact<Parameters...>> &encloser,
                   size_t contact_index)
    : RenormalizationVW<Base, Contact<Parameters...>>::
          InteractKernel(ex_policy, encloser, contact_index),
      contact_inv_rho0_k_(encloser.contact_inv_rho0_[contact_index]),
      contact_mass_k_(encloser.dv_contact_mass_[contact_index]->DelegatedData(ex_policy)),
      renormalizationVW_(encloser.dv_renormalizationVW_->DelegatedData(ex_policy))
{
}
//=================================================================================================//
template <typename... Parameters>
void RenormalizationVW<Contact<Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{

    Real renormalize_factor = 0.;
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        renormalize_factor += this->W_ij(index_i, index_j) * contact_inv_rho0_k_ * contact_mass_k_[index_j];
        ;
    }

    this->renormalizationVW_[index_i] += renormalize_factor;
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
#endif // VIRTOSIM_DENSITY_REGULARIZATION_HPP_D6998CB9_73EC_414F_9E1A_5FAE9F4CBFDC

#endif // VIRTOSIM_DENSITY_REGULARIZATION_20COPY_HPP_CEB4A1AB_A372_4275_BCE5_A03E8572E18C
