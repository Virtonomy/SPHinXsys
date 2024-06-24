/**
 * @file 	base_particles.hpp
 * @brief 	This is the implementation of the template functions in base_particles.h
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef BASE_PARTICLES_HPP
#define BASE_PARTICLES_HPP

#include "base_particles.h"

//=====================================================================================================//
namespace SPH
{
//=================================================================================================//
template <typename DataType>
StdLargeVec<DataType> *BaseParticles::initializeVariable(DiscreteVariable<DataType> *variable,
                                                         DataType initial_value)
{
    if (variable->DataField() == nullptr)
    {
        variable->allocateDataField(particles_bound_, initial_value);
        constexpr int type_index = DataTypeIndex<DataType>::value;
        std::get<type_index>(all_particle_data_).push_back(variable->DataField());
    }
    else
    {
        std::cout << "\n Error: the variable '" << variable->Name() << "' has already been allocated!" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }
    return variable->DataField();
}
//=================================================================================================//
template <typename DataType, class InitializationFunction>
StdLargeVec<DataType> *BaseParticles::initializeVariable(DiscreteVariable<DataType> *variable,
                                                         const InitializationFunction &initialization)
{

    auto &contained_data = *initializeVariable(variable);
    for (size_t i = 0; i != total_real_particles_; ++i)
    {
        contained_data[i] = initialization(i); // Here, lambda function is applied for initialization.
    }
    return &contained_data;
}
//=================================================================================================//
template <typename DataType>
DataType *BaseParticles::registerSingleVariable(const std::string &variable_name, DataType initial_value)
{
    SingleVariable<DataType> *variable = findVariableByName<DataType>(all_single_variables_, variable_name);

    return variable != nullptr
               ? variable->ValueAddress()
               : addVariableToAssemble<DataType>(all_single_variables_,
                                                 all_global_variable_ptrs_, variable_name, initial_value)
                     ->ValueAddress();
}
//=================================================================================================//
template <typename DataType>
DataType *BaseParticles::getSingleVariableByName(const std::string &variable_name)
{
    SingleVariable<DataType> *variable = findVariableByName<DataType>(all_single_variables_, variable_name);

    if (variable != nullptr)
    {
        return variable->ValueAddress();
    }

    std::cout << "\nError: the variable '" << variable_name << "' is not registered!\n";
    std::cout << __FILE__ << ':' << __LINE__ << std::endl;
    return nullptr;

    return nullptr;
}
//=================================================================================================//
template <typename DataType>
DiscreteVariable<DataType> *BaseParticles::addSharedVariable(const std::string &variable_name)
{
    DiscreteVariable<DataType> *variable = findVariableByName<DataType>(all_discrete_variables_, variable_name);
    if (variable == nullptr)
    {
        variable = addVariableToAssemble<DataType>(all_discrete_variables_, all_discrete_variable_ptrs_, variable_name);
    }
    return variable;
}
//=================================================================================================//
template <typename DataType, typename... Args>
StdLargeVec<DataType> *BaseParticles::registerSharedVariable(const std::string &variable_name, Args &&...args)
{

    DiscreteVariable<DataType> *variable = addSharedVariable<DataType>(variable_name);

    if (variable->DataField() == nullptr)
    {
        initializeVariable(variable, std::forward<Args>(args)...);
    }
    return variable->DataField();
}
//=================================================================================================//
template <typename DataType>
StdLargeVec<DataType> *BaseParticles::registerSharedVariableFrom(
    const std::string &new_name, const std::string &old_name)
{
    DiscreteVariable<DataType> *variable = findVariableByName<DataType>(all_discrete_variables_, old_name);

    if (variable == nullptr)
    {
        std::cout << "\nError: the old variable '" << old_name << "' is not registered!\n";
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }

    StdLargeVec<DataType> &old_data = *variable->DataField();
    return registerSharedVariable<DataType>(new_name, [&](size_t index)
                                            { return old_data[index]; });
}
//=================================================================================================//
template <typename DataType>
StdLargeVec<DataType> *BaseParticles::registerSharedVariableFrom(
    const std::string &new_name, const StdLargeVec<DataType> &geometric_data)
{
    return registerSharedVariable<DataType>(new_name, [&](size_t index)
                                            { return geometric_data[index]; });
}
//=================================================================================================//
template <typename DataType>
DiscreteVariable<DataType> *BaseParticles::getVariableByName(const std::string &variable_name)
{
    DiscreteVariable<DataType> *variable = findVariableByName<DataType>(all_discrete_variables_, variable_name);
    if (variable == nullptr)
    {
        std::cout << "\nError: the variable '" << variable_name << "' is not registered!\n";
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }
    return variable;
}
//=================================================================================================//
template <typename DataType>
StdLargeVec<DataType> *BaseParticles::getVariableDataByName(const std::string &variable_name)
{
    DiscreteVariable<DataType> *variable = getVariableByName<DataType>(variable_name);

    if (variable->DataField() == nullptr)
    {
        std::cout << "\nError: the variable '" << variable_name << "' has not been allocated yet!\n";
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }

    return variable->DataField();
}
//=================================================================================================//
template <typename DataType>
DiscreteVariable<DataType> *BaseParticles::
    addVariableToList(ParticleVariables &variable_set, const std::string &variable_name)
{
    DiscreteVariable<DataType> *variable = findVariableByName<DataType>(all_discrete_variables_, variable_name);

    if (variable != nullptr)
    {
        DiscreteVariable<DataType> *listed_variable = findVariableByName<DataType>(variable_set, variable_name);

        if (listed_variable == nullptr)
        {
            constexpr int type_index = DataTypeIndex<DataType>::value;
            std::get<type_index>(variable_set).push_back(variable);
            return variable;
        }
    }
    else
    {
        std::cout << "\n Error: the variable '" << variable_name << "' to write is not particle data!" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }
    return nullptr;
}
//=================================================================================================//
template <typename DataType>
void BaseParticles::addVariableToSort(const std::string &variable_name)
{
    DiscreteVariable<DataType> *new_sortable =
        addVariableToList<DataType>(sortable_variables_, variable_name);
    if (new_sortable != nullptr)
    {
        constexpr int type_index = DataTypeIndex<DataType>::value;
        StdLargeVec<DataType> *variable_data = new_sortable->DataField();
        std::get<type_index>(sortable_data_).push_back(variable_data);
    }
}
//=================================================================================================//
template <typename DataType>
void BaseParticles::addVariableToWrite(const std::string &variable_name)
{
    addVariableToList<DataType>(variables_to_write_, variable_name);
}
//=================================================================================================//
template <typename DataType>
void BaseParticles::addVariableToRestart(const std::string &variable_name)
{
    addVariableToList<DataType>(variables_to_restart_, variable_name);
}
//=================================================================================================//
template <typename DataType>
void BaseParticles::addVariableToReload(const std::string &variable_name)
{
    addSharedVariable<DataType>(variable_name);
    addVariableToList<DataType>(variables_to_reload_, variable_name);
}
//=================================================================================================//
template <typename SequenceMethod>
void BaseParticles::sortParticles(SequenceMethod &sequence_method)
{
    StdLargeVec<size_t> &sequence = sequence_method.computingSequence(*this);
    particle_sorting_.sortingParticleData(sequence.data(), total_real_particles_);
}
//=================================================================================================//
template <typename DataType>
void BaseParticles::ResizeParticles::
operator()(DataContainerAddressKeeper<StdLargeVec<DataType>> &data_keeper, size_t new_size)
{
    for (size_t i = 0; i != data_keeper.size(); ++i)
    {
        data_keeper[i]->resize(new_size, ZeroData<DataType>::value);
    }
}
//=================================================================================================//
template <typename DataType>
void BaseParticles::CopyParticleData::
operator()(DataContainerAddressKeeper<StdLargeVec<DataType>> &data_keeper, size_t index, size_t another_index)
{
    for (size_t i = 0; i != data_keeper.size(); ++i)
    {
        (*data_keeper[i])[index] = (*data_keeper[i])[another_index];
    }
}
//=================================================================================================//
template <typename DataType>
void BaseParticles::WriteAParticleVariableToXml::
operator()(DataContainerAddressKeeper<DiscreteVariable<DataType>> &variables, ParticleData &all_particle_data)
{
    for (size_t i = 0; i != variables.size(); ++i)
    {
        size_t index = 0;
        StdLargeVec<DataType> &variable_data = *variables[i]->DataField();
        for (auto child = xml_parser_.first_element_->FirstChildElement(); child; child = child->NextSiblingElement())
        {
            xml_parser_.setAttributeToElement(child, variables[i]->Name(), variable_data[index]);
            index++;
        }
    }
}
//=================================================================================================//
template <typename DataType>
void BaseParticles::ReadAParticleVariableFromXml::
operator()(DataContainerAddressKeeper<DiscreteVariable<DataType>> &variables, BaseParticles *base_particles)
{
    for (size_t i = 0; i != variables.size(); ++i)
    {
        size_t index = 0;
        StdLargeVec<DataType> &variable_data = variables[i]->DataField() != nullptr
                                                   ? *variables[i]->DataField()
                                                   : *base_particles->initializeVariable<DataType>(variables[i]);
        for (auto child = xml_parser_.first_element_->FirstChildElement(); child; child = child->NextSiblingElement())
        {
            xml_parser_.queryAttributeValue(child, variables[i]->Name(), variable_data[index]);
            index++;
        }
    }
}
//=================================================================================================//
template <typename StreamType>
void BaseParticles::writeParticlesToVtk(StreamType &output_stream)
{
    size_t total_real_particles = total_real_particles_;

    // write sorted particles ID
    output_stream << "    <DataArray Name=\"SortedParticle_ID\" type=\"Int32\" Format=\"ascii\">\n";
    output_stream << "    ";
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        output_stream << i << " ";
    }
    output_stream << std::endl;
    output_stream << "    </DataArray>\n";

    // write unsorted particles ID
    output_stream << "    <DataArray Name=\"UnsortedParticle_ID\" type=\"Int32\" Format=\"ascii\">\n";
    output_stream << "    ";
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        output_stream << unsorted_id_[i] << " ";
    }
    output_stream << std::endl;
    output_stream << "    </DataArray>\n";

    // write integers
    constexpr int type_index_int = DataTypeIndex<int>::value;
    for (DiscreteVariable<int> *variable : std::get<type_index_int>(variables_to_write_))
    {
        StdLargeVec<int> &variable_data = *variable->DataField();
        output_stream << "    <DataArray Name=\"" << variable->Name() << "\" type=\"Int32\" Format=\"ascii\">\n";
        output_stream << "    ";
        for (size_t i = 0; i != total_real_particles; ++i)
        {
            output_stream << std::fixed << std::setprecision(9) << variable_data[i] << " ";
        }
        output_stream << std::endl;
        output_stream << "    </DataArray>\n";
    }

    // write scalars
    constexpr int type_index_Real = DataTypeIndex<Real>::value;
    for (DiscreteVariable<Real> *variable : std::get<type_index_Real>(variables_to_write_))
    {
        StdLargeVec<Real> &variable_data = *variable->DataField();
        output_stream << "    <DataArray Name=\"" << variable->Name() << "\" type=\"Float32\" Format=\"ascii\">\n";
        output_stream << "    ";
        for (size_t i = 0; i != total_real_particles; ++i)
        {
            output_stream << std::fixed << std::setprecision(9) << variable_data[i] << " ";
        }
        output_stream << std::endl;
        output_stream << "    </DataArray>\n";
    }

    // write vectors
    constexpr int type_index_Vecd = DataTypeIndex<Vecd>::value;
    for (DiscreteVariable<Vecd> *variable : std::get<type_index_Vecd>(variables_to_write_))
    {
        StdLargeVec<Vecd> &variable_data = *variable->DataField();
        output_stream << "    <DataArray Name=\"" << variable->Name() << "\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
        output_stream << "    ";
        for (size_t i = 0; i != total_real_particles; ++i)
        {
            Vec3d vector_value = upgradeToVec3d(variable_data[i]);
            output_stream << std::fixed << std::setprecision(9) << vector_value[0] << " " << vector_value[1] << " " << vector_value[2] << " ";
        }
        output_stream << std::endl;
        output_stream << "    </DataArray>\n";
    }

    // write matrices
    constexpr int type_index_Matd = DataTypeIndex<Matd>::value;
    for (DiscreteVariable<Matd> *variable : std::get<type_index_Matd>(variables_to_write_))
    {
        StdLargeVec<Matd> &variable_data = *variable->DataField();
        output_stream << "    <DataArray Name=\"" << variable->Name() << "\" type= \"Float32\"  NumberOfComponents=\"9\" Format=\"ascii\">\n";
        output_stream << "    ";
        for (size_t i = 0; i != total_real_particles; ++i)
        {
            Mat3d matrix_value = upgradeToMat3d(variable_data[i]);
            for (int k = 0; k != 3; ++k)
            {
                Vec3d col_vector = matrix_value.col(k);
                output_stream << std::fixed << std::setprecision(9) << col_vector[0] << " " << col_vector[1] << " " << col_vector[2] << " ";
            }
        }
        output_stream << std::endl;
        output_stream << "    </DataArray>\n";
    }
}
//=================================================================================================//
} // namespace SPH
#endif // BASE_PARTICLES_HPP