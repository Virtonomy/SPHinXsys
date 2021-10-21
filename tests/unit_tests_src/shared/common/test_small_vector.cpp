#include <gtest/gtest.h>
#include "small_vectors.h"

using namespace SPH;

TEST(small_vectors, getCosineOfAngleBetweenTwoVectors)
{
    Vec3d vector_1 = Vec3d(2, -4, 0);
    Vec3d vector_2 = Vec3d(3, 2, 5);

    Vec3d vector_3 = Vec3d(1, 0, 0);
    Vec3d vector_4 = Vec3d(0, 0, 1);

    Vec3d vector_5 = Vec3d(1, 0, 0);
    Vec3d vector_6 = Vec3d(-1, 0, 0);

    Real cos_teta_1 = getCosineOfAngleBetweenTwoVectors(vector_1, vector_2);
    Real cos_teta_ref_1 = -2 / (sqrt(20) * sqrt(38));

    Real cos_teta_2 = getCosineOfAngleBetweenTwoVectors(vector_3, vector_4);
    Real cos_teta_ref_2 = 0;

    Real cos_teta_3 = getCosineOfAngleBetweenTwoVectors(vector_5, vector_6);
    Real cos_teta_ref_3 = -1;
 
	EXPECT_EQ(cos_teta_1, cos_teta_ref_1);
    EXPECT_EQ(cos_teta_2, cos_teta_ref_2);
    EXPECT_EQ(cos_teta_3, cos_teta_ref_3);

}
//=================================================================================================//
TEST(small_vectors, getVectorProjectionOfVector)
{
    Vec3d vector_1 = Vec3d(1, 2, 3);
    Vec3d vector_2 = Vec3d(2, -1, 4);
    

    //angle exactly 90°
    Vec3d vector_3 = Vec3d( 0, 0, 0);
    Vec3d vector_4 = Vec3d(0, 1, 0);
    
    //vectors in opposite direction
    Vec3d vector_5 = Vec3d(1, 1, 0);
    Vec3d vector_6 = Vec3d(-1, -1, 0);

    Vec3d vector_7 = Vec3d(5.44009282066326e-15, -0.000058164135562623, 2.22856150706995e-07);
    Vec3d vector_8 = Vec3d(-9.6283166736093e-11, 0.999698258170137, -0.00394363294839056);

    Vec3d proj_vector_1 = getVectorProjectionOfVector(vector_1, vector_2);
    Vec3d proj_vector_1_ref = Vec3d(1.142857, -0.571428, 2.285714);

    Vec3d proj_vector_2 = getVectorProjectionOfVector(vector_3, vector_4);
    Vec3d proj_vector_2_ref = Vec3d(0, 0, 0);

    Vec3d proj_vector_3 = getVectorProjectionOfVector(vector_5, vector_6);
    Vec3d proj_vector_3_ref = Vec3d(1, 1, 0);

    //values from spring normal to surface
    Vec3d proj_vector_4 = getVectorProjectionOfVector(vector_7, vector_8);
    Vec3d proj_vector_4_ref = Vec3d(5.60191499112963e-15, -0.000058164109562364, 2.29447132681609e-07);

	for (size_t i = 0; i < 3; i++)
	{
        EXPECT_NEAR(proj_vector_1[i], proj_vector_1_ref[i], 1e-6);
        EXPECT_NEAR(proj_vector_2[i], proj_vector_2_ref[i], 1e-6);
        EXPECT_NEAR(proj_vector_3[i], proj_vector_3_ref[i], 1e-6);
        EXPECT_NEAR(proj_vector_4[i], proj_vector_4_ref[i], 1e-8);
	}

}
//=================================================================================================//
TEST(small_vectors, getPrincipalValuesFromMatrix)
{
    Mat3d stress_tensor_1 = Mat3d(
        Vec3d(50, 30, 20),
        Vec3d(30, -20, -10),
        Vec3d(20, -10, 10)
    );
    Vec3d principal_stress_1 = getPrincipalValuesFromMatrix(stress_tensor_1);
    Vec3d principal_stress_1_ref = {65.527, 11.531, -37.058};
    for (size_t i = 0; i < 3; i++) EXPECT_NEAR(principal_stress_1[i], principal_stress_1_ref[i], 1e-3);

    Mat3d stress_tensor_2 = Mat3d(
        Vec3d(12.5, 0, 0),
        Vec3d(0, -50.4, 0),
        Vec3d(0, 0, 15.3)
    );
    Vec3d principal_stress_2 = getPrincipalValuesFromMatrix(stress_tensor_2);
    Vec3d principal_stress_2_ref = {15.3, 12.5, -50.4};
    for (size_t i = 0; i < 3; i++) EXPECT_DOUBLE_EQ(principal_stress_2[i], principal_stress_2_ref[i]);
}

//=================================================================================================//
//=================================================================================================//
int main(int argc, char* argv[])
{	
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}