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

#ifndef VIRTOSIM_DENSITY_REGULARIZATION_H_FC0C6065_700E_479E_9943_3FBA4C37FCA7
#define VIRTOSIM_DENSITY_REGULARIZATION_H_FC0C6065_700E_479E_9943_3FBA4C37FCA7

#include "base_fluid_dynamics.h"
#include "interaction_ck.hpp"
#include "particle_functors_ck.h" // or wherever ParticleScopeTypeCK is defined

namespace SPH
{
namespace fluid_dynamics
{

//------------------------------------------------------------------
// forward declarations of Regularization<FlowType>
//------------------------------------------------------------------
template <typename...>
class Regularization;

template <>
class Regularization<Internal>
{
  public:
    Regularization(BaseParticles *particles) {};

    class ComputingKernel
    {
      public:
        template <class ExecutionPolicy, class ComputingKernelType>
        ComputingKernel(const ExecutionPolicy &ex_policy,
                        Regularization<Internal> &encloser,
                        ComputingKernelType &computing_kernel){};

        Real operator()(Real &rho_sum) { return rho_sum; };
        Real operator()(Real rho_sum, Real /*rho_old*/, size_t) const
        {
            return rho_sum;
        }
    };
};

template <>
class Regularization<FreeSurface>
{
  public:
    Regularization(BaseParticles *particles) {};

    class ComputingKernel
    {
      public:
        template <class ExecutionPolicy, class ComputingKernelType>
        ComputingKernel(const ExecutionPolicy &ex_policy,
                        Regularization<FreeSurface> &encloser,
                        ComputingKernelType &computing_kernel)
            : rho0_(computing_kernel.InitialDensity()){};

        Real operator()(Real &rho_sum, Real &rho, size_t index_i) { return SMAX(rho_sum, rho0_); };

      protected:
        Real rho0_;
    };
};

template <>
class Regularization<FreeStream>
{
  public:
    Regularization(BaseParticles *particles) {}

    class ComputingKernel
    {
      public:
        template <class ExecutionPolicy, class ComputingKernelType>
        ComputingKernel(const ExecutionPolicy & /*ex_policy*/,
                        Regularization<FreeStream> & /*encloser*/,
                        ComputingKernelType &computing_kernel)
            : rho0_(computing_kernel.InitialDensity()) {}

        Real operator()(Real rho_sum, Real rho, size_t index_i) const
        {
            // if (rho_sum < rho)
            //     return rho_sum + (rho - rho_sum) * rho0_ / rho;
            return rho_sum;
        }

      protected:
        Real rho0_;
    };
};

class renormalizationwithVW;
template <>
class Regularization<renormalizationwithVW>
{
  public:
    Regularization(BaseParticles *particles) : dv_renormalizationVW_(particles->template getVariableByName<Real>("ReNormalizationVW"))
    {
    }

    class ComputingKernel
    {
      public:
        template <class ExecutionPolicy, class ComputingKernelType>
        ComputingKernel(const ExecutionPolicy &ex_policy,
                        Regularization<renormalizationwithVW> &encloser,
                        ComputingKernelType &computing_kernel)
            : rho0_(computing_kernel.InitialDensity()),
              renormalizationVW_(encloser.dv_renormalizationVW_->DelegatedData(ex_policy))
        {
        }
        // spacing ratio  Δx / h   : 1 / 2.6  -> 0.385
        // kernel bias C (3-D Wendland-C2) : 0.22
        // f_max = 1+ c(Δx / h)^2 = 1+ C*(0.385)^2
        // Real operator()(Real rho_sum, Real rho_old, size_t index_i) const
        // {
        //     Real f_max = 1.032;
        //     Real f = renormalizationVW_[index_i];
        //     if (f < 0.0)
        //         f = 0.0;
        //     else if (f > f_max)
        //         f = f_max;
        //     if (f < f_max)
        //     {
        //         rho_sum = rho_sum + rho0_ * (f_max - f);
        //     }
        //     if (rho_sum < rho_old)
        //         return rho_sum + (rho_old - rho_sum) * rho0_ / rho_old;
        //     return rho_sum;
        // }

        Real operator()(Real rho_sum, Real rho_old, size_t i) const
        {
            constexpr Real f_max = 1.032; // Wendland-C2, Δx/h ≈ 0.385
            Real f = renormalizationVW_[i];
            if (f < 0)
                f = 0;
            if (f > f_max)
                f = f_max;
            /* Shepard fill-in term with reference density */
            Real rho_f_sum = rho_sum + rho0_ * (f_max - f);
            Real rho_freestream_sum = 0.;
            if (rho_sum < rho_old)
                rho_freestream_sum = rho_sum + (rho_old - rho_sum) * rho0_ / rho_old;
            else
                rho_freestream_sum = rho_sum;
            Real rho_diff = abs(rho_f_sum - rho_freestream_sum);
            Real rho_max_change = 0.01 * 0.025 * rho0_; // only change 0.01%
            if (rho_diff < rho_max_change)
                return rho_old + rho_diff;
            else
                return rho_old + rho_max_change * ((rho_f_sum - rho_freestream_sum > 0) ? 1 : -1);
        }

        // Real operator()(Real rho_sum, Real rho_old, size_t i) const
        // {
        //     constexpr Real f_max = 1.00937; // Wendland-C2, Δx/h ≈ 0.385

        //     Real f = renormalizationVW_[i];

        //     if (f > 0.95)
        //     {
        //         if (rho_sum < rho_old)
        //             return rho_sum + (rho_old - rho_sum) * rho0_ / rho_old;
        //         return rho_sum;
        //     }
        //     if (f < 0)
        //         f = 0;
        //     if (f > f_max)
        //         f = f_max;

        //     /* Shepard fill-in term with reference density */
        //     Real rho_f_sum = rho_sum + rho0_ * (f_max - f);
        //     Real rho_freestream_sum = 0.;
        //     if (rho_sum < rho_old)
        //         rho_freestream_sum = rho_sum + (rho_old - rho_sum) * rho_old / rho_old;
        //     else
        //         rho_freestream_sum = rho_sum;

        //     Real rho_diff = abs(rho_f_sum - rho_freestream_sum);
        //     Real rho_max_change = 0.01 * 0.025 * rho0_; // only change 0.025%
        //     if (rho_diff < rho_max_change)
        //         return rho_old + rho_diff;
        //     else
        //         return rho_old + rho_max_change * ((rho_f_sum - rho_freestream_sum > 0) ? 1 : -1);
        // }

      protected:
        Real rho0_;
        Real *renormalizationVW_;
    };
    DiscreteVariable<Real> *dv_renormalizationVW_;
};

template <typename... RelationTypes>
class DensityRegularization;

template <template <typename...> class RelationType, typename... Parameters>
class DensityRegularization<Base, RelationType<Parameters...>>
    : public Interaction<RelationType<Parameters...>>
{
  public:
    template <class DynamicsIdentifier>
    explicit DensityRegularization(DynamicsIdentifier &identifier);
    virtual ~DensityRegularization() {};

    class InteractKernel : public Interaction<RelationType<Parameters...>>::InteractKernel
    {
      public:
        template <class ExecutionPolicy, typename... Args>
        InteractKernel(const ExecutionPolicy &ex_policy,
                       DensityRegularization<Base, RelationType<Parameters...>> &encloser,
                       Args &&...args);
        Real InitialDensity() { return rho0_; };

      protected:
        Real *rho_, *mass_, *rho_sum_, *Vol_;
        Real rho0_, inv_sigma0_;
    };

  protected:
    DiscreteVariable<Real> *dv_rho_, *dv_mass_, *dv_rho_sum_, *dv_Vol_;
    Real rho0_, inv_sigma0_;
};

template <class FlowType, class ParticleScopeType, typename... Parameters>
class DensityRegularization<Inner<WithUpdate, FlowType, ParticleScopeType, Parameters...>>
    : public DensityRegularization<Base, Inner<Parameters...>>
{
    using RegularizationKernel = typename Regularization<FlowType>::ComputingKernel;
    using FreeStreamRegularizationKernel = typename Regularization<FreeStream>::ComputingKernel;
    using ParticleScopeTypeKernel = typename ParticleScopeTypeCK<ParticleScopeType>::ComputingKernel;

  public:
    explicit DensityRegularization(Relation<Inner<Parameters...>> &inner_relation);
    virtual ~DensityRegularization() {};

    class InteractKernel
        : public DensityRegularization<Base, Inner<Parameters...>>::InteractKernel
    {
      public:
        template <class ExecutionPolicy>
        InteractKernel(const ExecutionPolicy &ex_policy,
                       DensityRegularization<Inner<WithUpdate, FlowType, ParticleScopeType, Parameters...>> &encloser);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real W0_;
    };

    class UpdateKernel
        : public DensityRegularization<Base, Inner<Parameters...>>::InteractKernel
    {
      public:
        template <class ExecutionPolicy>
        UpdateKernel(const ExecutionPolicy &ex_policy,
                     DensityRegularization<Inner<WithUpdate, FlowType, ParticleScopeType, Parameters...>> &encloser);
        void update(size_t index_i, Real dt = 0.0);

      protected:
        RegularizationKernel regularization_;
        FreeStreamRegularizationKernel freestream_regularization_;
        ParticleScopeTypeKernel particle_scope_;
        Real *Vol_;
        Real *mass_;
    };

  protected:
    Regularization<FlowType> regularization_method_;
    Regularization<FreeStream> freestream_regularization_method_;
    ParticleScopeTypeCK<ParticleScopeType> within_scope_method_;
};

template <typename... Parameters>
class DensityRegularization<Contact<Parameters...>>
    : public DensityRegularization<Base, Contact<Parameters...>>
{
  public:
    explicit DensityRegularization(Relation<Contact<Parameters...>> &contact_relation);
    virtual ~DensityRegularization() {};

    class InteractKernel
        : public DensityRegularization<Base, Contact<Parameters...>>::InteractKernel
    {
      public:
        template <class ExecutionPolicy>
        InteractKernel(const ExecutionPolicy &ex_policy,
                       DensityRegularization<Contact<Parameters...>> &encloser,
                       size_t contact_index);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real contact_inv_rho0_k_;
        Real *contact_mass_k_;
    };

  protected:
    StdVec<Real> contact_inv_rho0_;
    StdVec<DiscreteVariable<Real> *> dv_contact_mass_;
};

using DensityRegularizationComplex = DensityRegularization<Inner<WithUpdate, Internal, AllParticles>, Contact<>>;
using DensityRegularizationComplexFreeSurface = DensityRegularization<Inner<WithUpdate, FreeSurface, AllParticles>, Contact<>>;
using DensityRegularizationComplexFreeStream = DensityRegularization<Inner<WithUpdate, FreeStream, AllParticles>, Contact<>>;
using DensityRegularizationComplexInternalPressureBoundary = DensityRegularization<Inner<WithUpdate, Internal, ExcludeBufferParticles>, Contact<>>;
using DensityRegularizationComplexRenormalization = DensityRegularization<Inner<WithUpdate, renormalizationwithVW, BulkParticles>, Contact<>>;

} // namespace fluid_dynamics
} // namespace SPH

#endif // VIRTOSIM_DENSITY_REGULARIZATION_H_FC0C6065_700E_479E_9943_3FBA4C37FCA7
