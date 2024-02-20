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
 * @file 	viscous_dynamics.h
 * @brief Here, we define the algorithm classes for computing viscous forces in fluids.
 * @details TBD.
 * @author	Xiangyu Hu
 */

#include "fluid_dynamics_complex.hpp"
#include "fluid_dynamics_inner.hpp"
#include "fluid_surface_inner.hpp"
namespace SPH
{
class NablaWVInner : public LocalDynamics, public fluid_dynamics::FluidDataInner
{
  public:
    explicit NablaWVInner(BaseInnerRelation &inner_relation) : LocalDynamics(inner_relation.getSPHBody()), fluid_dynamics::FluidDataInner(inner_relation),
                                                               kernel_sum_(*this->particles_->template registerSharedVariable<Vecd>("KernelSummation")){};
    virtual ~NablaWVInner() = default;

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        kernel_sum_[index_i] = Vecd::Zero();
        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            kernel_sum_[index_i] += inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
        }
    };

  protected:
    StdLargeVec<Vecd> &kernel_sum_;
};

class NablaWV : public BaseInteractionComplex<NablaWVInner, fluid_dynamics::FluidContactData>
{
  protected:
  public:
    template <typename... Args>
    NablaWV(Args &&...args)
        : BaseInteractionComplex<NablaWVInner, fluid_dynamics::FluidContactData>(
              std::forward<Args>(args)...){};
    virtual ~NablaWV(){};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {

        // kernel_sum_[index_i] = Vecd::Zero();

        // const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        // for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        // {
        //     kernel_sum_[index_i] += inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
        // }
        for (size_t k = 0; k < contact_configuration_.size(); ++k)
        {
            Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                kernel_sum_[index_i] += contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.e_ij_[n];
            }
        }
    };
};

using NablaWVComplex = ComplexInteraction<NablaWVInner, NablaWV>;

} // namespace SPH