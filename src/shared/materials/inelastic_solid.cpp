#include "inelastic_solid.h"
#include "base_particles.hpp"

namespace SPH
{
//=================================================================================================//
void HardeningPlasticSolid::initializeLocalParameters(BaseParticles *base_particles)
{
    PlasticSolid::initializeLocalParameters(base_particles);
    base_particles->registerVariable(inverse_plastic_strain_, "InversePlasticRightCauchyStrain",
                                     [&](size_t i) -> Matd
                                     { return Matd::Identity(); });
    base_particles->registerVariable(hardening_parameter_, "HardeningParameter");
    base_particles->addVariableToRestart<Matd>("InversePlasticRightCauchyStrain");
    base_particles->addVariableToRestart<Real>("HardeningParameter");
}
//=================================================================================================//
Matd HardeningPlasticSolid::PlasticConstitutiveRelation(const Matd &F, size_t index_i, Real dt)
{
    // resulting stress is not rate dependent (does not depend on the speed of deformation)
    // inverse_plastic_strain_[index_i]=inverse plastic right Cauchy-Green based on previous time step
    // elastic trial left Cauchy-Green tensor (stress )
    Matd be = F * inverse_plastic_strain_[index_i] * F.transpose();
    // Remove volumetric part
    Matd normalized_be = be * pow(be.determinant(), -OneOverDimensions);
    // Mean value of the isochoric elastic trial tensor.
    Real normalized_be_isentropic = normalized_be.trace() * OneOverDimensions;
    // deviatoric Kirchhoff stress tensor
    Matd deviatoric_PK = DeviatoricKirchhoff(normalized_be - normalized_be_isentropic * Matd::Identity());
    Real deviatoric_PK_norm = deviatoric_PK.norm();

    // Yield function for isotropic hardening:
    Real trial_function = deviatoric_PK_norm -
                          sqrt_2_over_3_ * (hardening_modulus_ * hardening_parameter_[index_i] + yield_stress_);

    // how much plastic stress there is over the yield stress
    if (trial_function > 0.0)
    {   // Plastic correction: return the trial stress back to the yield surface along the deviatoric
        // stress direction

        // Effective shear modulus scaled by the isochoric elastic deformation state (direction of deformation).
        Real renormalized_shear_modulus = normalized_be_isentropic * G0_;

        // equivalent plastic strain increment
        // The increment is chosen so that the corrected deviatoric stress lies on the isotropically
        // hardened yield surface.
        Real relax_increment = 0.5 * trial_function / (renormalized_shear_modulus + hardening_modulus_ / 3.0);
        // accumulated equivalent plastic strain
        hardening_parameter_[index_i] += sqrt_2_over_3_ * relax_increment;
        // Correct the deviatoric Kirchhoff stress after plastic relaxation.
        deviatoric_PK -= 2.0 * renormalized_shear_modulus * relax_increment * deviatoric_PK / deviatoric_PK_norm;
        /// Reconstruct the corrected isochoric elastic tensor from the relaxed deviatoric stress.
        Matd relaxed_be = deviatoric_PK / G0_ + normalized_be_isentropic * Matd::Identity();
        // Normalize again to remove any volumetric contribution introduced by the correction.
        normalized_be = relaxed_be * pow(relaxed_be.determinant(), -OneOverDimensions);
    }
    Matd inverse_F = F.inverse();
    Matd inverse_F_T = inverse_F.transpose();
    // Update the stored inverse plastic tensor using the corrected elastic tensor.
    inverse_plastic_strain_[index_i] = inverse_F * normalized_be * inverse_F_T;

    // returns the first Piola-Kirchhoff stress tensor
    return (deviatoric_PK + VolumetricKirchhoff(F.determinant()) * Matd::Identity()) * inverse_F_T;
}
//=================================================================================================//
void MultiLinearHardeningPlasticSolid::initializeLocalParameters(BaseParticles *base_particles)
{
    PlasticSolid::initializeLocalParameters(base_particles);
    base_particles->registerVariable(inverse_plastic_strain_, "InversePlasticRightCauchyStrain",
                                     [&](size_t i) -> Matd
                                     { return Matd::Identity(); });
    base_particles->registerVariable(hardening_parameter_, "HardeningParameter");
    base_particles->registerVariable(hardening_modulus_, "HardeningModulus");

    base_particles->addVariableToRestart<Matd>("InversePlasticRightCauchyStrain");
    base_particles->addVariableToRestart<Real>("HardeningParameter");
    base_particles->addVariableToRestart<Real>("HardeningModulus");
}
//=================================================================================================//
Matd MultiLinearHardeningPlasticSolid::PlasticConstitutiveRelation(const Matd &F, size_t index_i, Real dt)
{

    Matd be = F * inverse_plastic_strain_[index_i] * F.transpose();
    Matd normalized_be = be * pow(be.determinant(), -OneOverDimensions);
    Real normalized_be_isentropic = normalized_be.trace() * OneOverDimensions;
    Matd deviatoric_PK = DeviatoricKirchhoff(normalized_be - normalized_be_isentropic * Matd::Identity());
    Real deviatoric_PK_norm = deviatoric_PK.norm();
    // Sigma is the current yield stress based on the spline created from the experimental data
    Real sigma = hardening_stresses_(hardening_parameter_[index_i])(0);
    // the hardening modulus is now a function of the hardening parameter, which comes from the spline
    hardening_modulus_[index_i] = hardening_moduluses_(hardening_parameter_[index_i])(0);

    // instead of using the yield_stress_ parameter, we use the yield stress from the spline because the yield
    // stress is a function of the hardening parameter and changes as the material deforms plastically
    Real trial_function =
        deviatoric_PK_norm - (sqrt_2_over_3_ * sigma);

    // check if the trial function is greater than zero at every yield point not just at the intial yield stress
    if (trial_function > 0.0)
    {
        Real renormalized_shear_modulus = normalized_be_isentropic * G0_;
        Real relax_increment = trial_function / (2.0 * renormalized_shear_modulus + (2.0 / 3.0) * hardening_modulus_[index_i]);
        hardening_parameter_[index_i] += sqrt_2_over_3_ * relax_increment;
        deviatoric_PK -= 2.0 * renormalized_shear_modulus * relax_increment * deviatoric_PK / deviatoric_PK_norm;
        Matd relaxed_be = deviatoric_PK / G0_ + normalized_be_isentropic * Matd::Identity();
        normalized_be = relaxed_be * pow(relaxed_be.determinant(), -OneOverDimensions);
    }
    Matd inverse_F = F.inverse();
    Matd inverse_F_T = inverse_F.transpose();
    inverse_plastic_strain_[index_i] = inverse_F * normalized_be * inverse_F_T;

    return (deviatoric_PK + VolumetricKirchhoff(F.determinant()) * Matd::Identity()) * inverse_F_T;
}

//=================================================================================================//
} // namespace SPH
