#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <iostream>

using namespace Eigen;


int main(){

	/*************** Fixed sized matrix, stored on the stack ***************/

	Matrix3d f; // 3 x 3 matrix <=> Matrix<double, 3, 3>; 

	// Matrix3f f(3,3);
	// Matrix3i f(3,3); 
	
    /* Initilization 1 */
	f << 1, 2, 3, 
	     4, 5, 6, 
	     7, 8, 9;

	/* Initilization 2 */
	// f = Matrix3d::Random();

	/* Initilization 3 */
	// f = Matrix3d::Constant(1.0);  // <=> f = Matrix3d::Ones();

 	/* Initilization 4 */
	// f = Matrix3d::Identity();

 	/* Initilization 5 */
	// f = Matrix3d::Zero();

	std::cout << "The fixed matrix is: " << std::endl;
	std::cout << f << std::endl;
	std::cout << "The size of fixed matrix is: " << f.size() << std::endl << std::endl;


	/*************** Dynamic matrix -> resizable, stored on the heap ***************/

	MatrixXd d1; // MatrixXf, <=> Matrix<double, Dynamic, Dynamic>; 
	// typedef Matrix<double, Dynamic, Dynamic> MatrixXd;

	/* Initilization 1 */
	d1.resize(2,2);
	d1(0,0) = 1.0;
	d1(0,1) = 2.0;
	d1(1,0) = 3.0;
	d1(1,1) = 4.0;

	/* Initilization 2 */
	// d1 = MatrixXd::Random(5, 5);

	/* Initilization 3 */
	// d1 = MatrixXd::Constant(5,5,0.1);

	std::cout << "The dynamic matrix d1 is: " << std::endl;
	std::cout << d1 << std::endl;
	std::cout << "The size of dynamic matrix d1 is: " << d1.size() << std::endl;

	std::cout << std::endl;

	MatrixXd d2(2,5);

	d2.resize(4,3);
	std::cout << "The size of dynamic matrix d2 is: "
              << d2.rows() << " x " << d2.cols() << std::endl;

    std::cout << "It has " << d2.size() << " coefficients" << std::endl << std::endl;


    /*************** Matrix Calculation ***************/

	// matrix.transpose()
	// matrix.inverse()
	// matrix.array().square()
	// matrix.sum()
	// matrix.prod()
	// matrix.mean()
	// matrix.minCoeff()
	// matrix.maxCoeff()
	// matrix.trace()

    Matrix<float, 1, 2> m1; // 1 x 2
    m1 << 1, 2;

    Matrix<float, 2, 2> m2; // 2 x 2
    m2 << 1, 2, 3, 4;

    MatrixXf m3 = m1 * m2;  // 1 x 2

    std::cout << "Matrix m3 is: " << m3 << std::endl << std::endl;


    VectorXf v(2); // 2 x 1 column vector
    v << 5, 6;

    MatrixXf m4 = v.transpose() * m2; 
    std::cout << "Matrix m4 is: " << m4 << std::endl;

    // Add 1 to all the elements of m2
    MatrixXf m5 = m2.array() + 1;
    std::cout << "Matrix m5 is: " << std::endl;
    std::cout << m5 << std::endl << std::endl;

    Vector3f v1 = Vector3f::Random();
	Vector3f v2 = Vector3f::Random();
	std::cout << v1*v2.transpose() << std::endl << std::endl;

	std::cout << v1.dot(v2) << std::endl << std::endl;
	std::cout << v1.normalized() << std::endl << std::endl;
	std::cout << v1.cross(v2) << std::endl << std::endl;

	Vector3f s = Vector3f::Random();
	std::cout << s << std::endl << std::endl;
	Vector4f q = s.homogeneous();
	std::cout << q << std::endl << std::endl;

	std::cout << v1.array() * v2.array() << std::endl << std::endl;
	std::cout << v1.array().sin() << std::endl << std::endl;

	return 0;

}