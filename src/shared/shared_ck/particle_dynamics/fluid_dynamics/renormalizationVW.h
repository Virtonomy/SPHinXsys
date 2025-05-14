/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file density_regularization.h
 * @brief Here, we define the algorithm classes for computing
 * the density of a continuum by kernel function summation.
 * @details We are using templates and their explicit or partial specializations
 * to identify variations of the interaction types..
 * @author Xiangyu Hu
 */

#ifndef VIRTOSIM_RENORMALIZATIONVW_H_ABD7A600_6AA3_49EE_BC9C_AF356672CE05
#define VIRTOSIM_RENORMALIZATIONVW_H_ABD7A600_6AA3_49EE_BC9C_AF356672CE05

#include "base_fluid_dynamics.h"
#include "interaction_ck.hpp"
#include "particle_functors_ck.h" // or wherever ParticleScopeTypeCK is defined

namespace SPH
{
namespace fluid_dynamics
{
template <typename... RelationTypes>
class RenormalizationVW;

template <template <typename...> class RelationType, typename... Parameters>
class RenormalizationVW<Base, RelationType<Parameters...>>
    : public Interaction<RelationType<Parameters...>>
{
  public:
    template <class DynamicsIdentifier>
    explicit RenormalizationVW(DynamicsIdentifier &identifier);
    virtual ~RenormalizationVW() {};

  protected:
    DiscreteVariable<Real> *dv_Vol_;
    DiscreteVariable<Real> *dv_renormalizationVW_;
};

template <class FlowType, class ParticleScopeType, typename... Parameters>
class RenormalizationVW<Inner<InteractionOnly, FlowType, ParticleScopeType, Parameters...>>
    : public RenormalizationVW<Base, Inner<Parameters...>>
{
    using ParticleScopeTypeKernel = typename ParticleScopeTypeCK<ParticleScopeType>::ComputingKernel;

  public:
    explicit RenormalizationVW(Relation<Inner<Parameters...>> &inner_relation);
    virtual ~RenormalizationVW() {};

    class InteractKernel
        : public RenormalizationVW<Base, Inner<Parameters...>>::InteractKernel
    {
      public:
        template <class ExecutionPolicy>
        InteractKernel(const ExecutionPolicy &ex_policy,
                       RenormalizationVW<Inner<InteractionOnly, FlowType, ParticleScopeType, Parameters...>> &encloser);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real *Vol_;
        Real *renormalizationVW_;
        Real W0_;
    };

  protected:
    ParticleScopeTypeCK<ParticleScopeType> within_scope_method_;
};

template <typename... Parameters>
class RenormalizationVW<Contact<Parameters...>>
    : public RenormalizationVW<Base, Contact<Parameters...>>
{
  public:
    explicit RenormalizationVW(Relation<Contact<Parameters...>> &contact_relation);
    virtual ~RenormalizationVW() {};

    class InteractKernel
        : public RenormalizationVW<Base, Contact<Parameters...>>::InteractKernel
    {
      public:
        template <class ExecutionPolicy>
        InteractKernel(const ExecutionPolicy &ex_policy,
                       RenormalizationVW<Contact<Parameters...>> &encloser,
                       size_t contact_index);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real contact_inv_rho0_k_;
        Real *contact_mass_k_;
        Real *renormalizationVW_;
    };

  protected:
    StdVec<Real> contact_inv_rho0_;
    StdVec<DiscreteVariable<Real> *> dv_contact_mass_;
};
using RegularizationVWComplex = RenormalizationVW<Inner<InteractionOnly, Internal, AllParticles>, Contact<>>;

} // namespace fluid_dynamics
} // namespace SPH

#endif // VIRTOSIM_RENORMALIZATIONVW_H_ABD7A600_6AA3_49EE_BC9C_AF356672CE05
