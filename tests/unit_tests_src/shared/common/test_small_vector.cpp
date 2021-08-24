#include <gtest/gtest.h>
#include "small_vectors.h"

using namespace SPH;

TEST(AngleBetweenTwo3DVectors, getAngleBetweenTwo3DVectors)
{
    Vec3d vector_1 = Vec3d(2, -4, 0);
    Vec3d vector_2 = Vec3d(3, 2, 5);

    Vec3d vector_3 = Vec3d(1, 0, 0);
    Vec3d vector_4 = Vec3d(0, 0, 1);

    Vec3d vector_5 = Vec3d(1, 0, 0);
    Vec3d vector_6 = Vec3d(-1, 0, 0);

    Real cos_teta_1 = getAngleBetweenTwo3DVectors(vector_1, vector_2);
    Real cos_teta_ref_1 = -2 / (sqrt(20) * sqrt(38));

    Real cos_teta_2 = getAngleBetweenTwo3DVectors(vector_3, vector_4);
    Real cos_teta_ref_2 = 0;

    Real cos_teta_3 = getAngleBetweenTwo3DVectors(vector_5, vector_6);
    Real cos_teta_ref_3 = -1;
 
	EXPECT_EQ(cos_teta_1, cos_teta_ref_1);
    EXPECT_EQ(cos_teta_2, cos_teta_ref_2);
    EXPECT_EQ(cos_teta_3, cos_teta_ref_3);

}
//=================================================================================================//
TEST(VectorProjectionOf3DVector, getVectorProjectionOf3DVector)
{
    Vec3d vector_1 = Vec3d(1, 2, 3);
    Vec3d vector_2 = Vec3d(2, -1, 4);
    

    //angle more than 90°
    Vec3d vector_3 = Vec3d(1, 1, 0);
    Vec3d vector_4 = Vec3d(-1, 0, 0);
    
    //vectors in opposite direction
    Vec3d vector_5 = Vec3d(1, 1, 0);
    Vec3d vector_6 = Vec3d(-1, -1, 0);

    Vec3d proj_vector_1 = getVectorProjectionOf3DVector(vector_1, vector_2);
    Vec3d proj_vector_1_ref = Vec3d(1.142857, -0.571428, 2.285714);

    Vec3d proj_vector_2 = getVectorProjectionOf3DVector(vector_3, vector_4);
    Vec3d proj_vector_2_ref = Vec3d(1, 0, 0);

    Vec3d proj_vector_3 = getVectorProjectionOf3DVector(vector_5, vector_6);
    Vec3d proj_vector_3_ref = Vec3d(1, 1, 0);

	for (size_t i = 0; i < 3; i++)
	{
        EXPECT_NEAR(proj_vector_1[i], proj_vector_1_ref[i], 1e-6);
        EXPECT_NEAR(proj_vector_2[i], proj_vector_2_ref[i], 1e-6);
        EXPECT_NEAR(proj_vector_3[i], proj_vector_3_ref[i], 1e-6);
	}

}
//=================================================================================================//
int main(int argc, char* argv[])
{	
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}