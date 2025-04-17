#ifndef FLOWRATE_CK_H
#define FLOWRATE_CK_H

#include "base_fluid_dynamics.h"
#include "interaction_ck.hpp"
#include "kernel_correction_ck.hpp"
#include "particle_functors_ck.h"

namespace SPH
{
namespace fluid_dynamics
{

//--------------------------------------------------------------------------------------
// Base class for Flowrate
//--------------------------------------------------------------------------------------
template <class BaseInteractionType>
class FlowrateCKBase : public BaseInteractionType
{
  public:
    template <class DynamicsIdentifier>
    explicit FlowrateCKBase(DynamicsIdentifier &identifier);
    virtual ~FlowrateCKBase() {}

  protected:
    DiscreteVariable<Real> *dv_Vol_;  ///< "VolumetricMeasure"
    DiscreteVariable<Vecd> *dv_dpos_; ///< "Position"
    DiscreteVariable<Vecd> *dv_vel_;
    DiscreteVariable<Vecd> *dv_volumetric_flux_;
};

//--------------------------------------------------------------------------------------
// Main template declaration for FLOWRATE_CK
//--------------------------------------------------------------------------------------
template <typename...>
class FlowrateCK;

template <class UpdatePolicy, class KernelCorrectionType, class ResolutionType, typename... Parameters>
class FlowrateCK<Inner<UpdatePolicy, KernelCorrectionType, ResolutionType, Parameters...>>
    : public FlowrateCKBase<Interaction<Inner<Parameters...>>>
{
    using BaseInteraction = FlowrateCKBase<Interaction<Inner<Parameters...>>>;
    using CorrectionKernel = typename KernelCorrectionType::ComputingKernel;

  public:
    explicit FlowrateCK(Relation<Inner<Parameters...>> &inner_relation);

    virtual ~FlowrateCK() {}

    //====================== Interact Kernel ======================//
    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        CorrectionKernel correction_;
        Real *Vol_;
        Vecd *dpos_;
        Vecd *vel_;
        Vecd *volumetric_flux_;
    };
    //====================== Update Kernel ======================//
    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(size_t index_i, Real dt = 0.0);

      protected:
        Real *vol_;
        Vecd *volumetric_flux_;
    };

  protected:
    KernelCorrectionType kernel_correction_;
};
//----------------------------------------------
//  2) Partial specialization for Contact<...>
//----------------------------------------------
template <class KernelCorrectionType, class ResolutionType, typename... Parameters>
class FlowrateCK<Contact<Wall, KernelCorrectionType, ResolutionType, Parameters...>>
    : public FlowrateCKBase<Interaction<Contact<Wall, Parameters...>>>
{
    using BaseInteraction = FlowrateCKBase<Interaction<Contact<Wall, Parameters...>>>;
    using CorrectionKernel = typename KernelCorrectionType::ComputingKernel;

  public:
    explicit FlowrateCK(Relation<Contact<Parameters...>> &contact_relation);
    virtual ~FlowrateCK() {}

    //====================== Interact Kernel ======================//
    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        CorrectionKernel correction_;
        Vecd *vel_;
        Vecd *contact_wall_vel_;
        Vecd *volumetric_flux_;
        Real *contact_wall_Vol_;
    };

  protected:
    KernelCorrectionType kernel_correction_;
    StdVec<DiscreteVariable<Vecd> *> dv_contact_wall_vel_;
    StdVec<DiscreteVariable<Real> *> dv_contact_wall_Vol_;
};

//--------------------------------------------------------------------------------------
// Alias Definitions for Specific Configurations
//--------------------------------------------------------------------------------------

using FlowrateCorrectionInnerNoCorrectionCK =
    FlowrateCK<
        Inner<WithUpdate, NoKernelCorrectionCK, SingleResolution>>;

using FlowrateCorrectionWallNoCorrectionCK =
    FlowrateCK<
        Inner<WithUpdate, NoKernelCorrectionCK, SingleResolution>,
        Contact<Wall, NoKernelCorrectionCK, SingleResolution>>;

using FlowrateCorrectionCorrectedComplexCK =
    FlowrateCK<
        Inner<WithUpdate, LinearCorrectionCK, SingleResolution>,
        Contact<Wall, LinearCorrectionCK, SingleResolution>>;

} // namespace fluid_dynamics
} // namespace SPH

#endif