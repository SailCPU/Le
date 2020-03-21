//
// Created by Sail Yang on 2019-09-21.
//

#include "LinearAlgebra.h"
using namespace le::math;

Eigen::MatrixXd LinearAlgebra::pseudoInverse(const Eigen::MatrixXd& am, double precision) {
    /*
    const Eigen::JacobiSVD<Eigen::MatrixXd> svd = am.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::JacobiSVD<Eigen::MatrixXd>::SingularValuesType sigma= svd.singularValues();
	for ( long i=0; i<sigma.cols(); ++i) {
		if ( sigma(i) > precision )
			sigma(i)=1.0/sigma(i);
		else
			sigma(i)=0;
	}

	return svd.matrixV() * sigma.asDiagonal()* svd.matrixU().transpose();
    */

    if (am.rows() < am.cols()){
        //RW_THROW("pseudoInverse require rows >= to cols!");
        Eigen::MatrixXd a = am.transpose();
        Eigen::JacobiSVD < Eigen::MatrixXd > svd = a.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
        double tolerance = precision * std::max(a.cols(), a.rows()) * svd.singularValues().array().abs().maxCoeff();
        return
                (svd.matrixV()
                 * Eigen::MatrixXd(
                        (svd.singularValues().array().abs() > tolerance).select(
                                svd.singularValues().array().inverse(), 0)).asDiagonal() * svd.matrixU().adjoint()).transpose();

    } else {

        Eigen::JacobiSVD < Eigen::MatrixXd > svd = am.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
        double tolerance = precision * std::max(am.cols(), am.rows()) * svd.singularValues().array().abs().maxCoeff();
        return
                svd.matrixV()
                * Eigen::MatrixXd(
                        (svd.singularValues().array().abs() > tolerance).select(
                                svd.singularValues().array().inverse(), 0)).asDiagonal() * svd.matrixU().adjoint();
    }
}

void LinearAlgebra::svd(const Eigen::MatrixXd& M, Eigen::MatrixXd& U, Eigen::VectorXd& sigma, Eigen::MatrixXd& V) {
    const Eigen::JacobiSVD<Eigen::MatrixXd> svd = M.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV);
    U = svd.matrixU();
    sigma = svd.singularValues();
    V = svd.matrixV();
}


bool LinearAlgebra::checkPenroseConditions(
        const Eigen::MatrixXd& A,
        const Eigen::MatrixXd& X,
        double prec)
{
    const Eigen::MatrixXd AX = A*X;
    const Eigen::MatrixXd XA = X*A;

    if (((AX*A)-A).lpNorm<Eigen::Infinity>() > prec)
        return false;
    if (((XA*X)-X).lpNorm<Eigen::Infinity>() > prec)
        return false;
    if ((AX.transpose() - AX).lpNorm<Eigen::Infinity>() > prec)
        return false;
    if ((XA.transpose() - XA).lpNorm<Eigen::Infinity>() > prec)
        return false;
    return true;
}

template<>
std::pair<LinearAlgebra::EigenMatrix<double>::type, LinearAlgebra::EigenVector<double>::type> LinearAlgebra::eigenDecompositionSymmetric<double>(const Eigen::MatrixXd& Am1) {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver;
    eigenSolver.compute(Am1);
    return std::make_pair(eigenSolver.eigenvectors(), eigenSolver.eigenvalues());
}

template<>
std::pair<LinearAlgebra::EigenMatrix<std::complex<double> >::type, LinearAlgebra::EigenVector<std::complex<double> >::type > LinearAlgebra::eigenDecomposition<double>(const Eigen::MatrixXd& Am1) {
    Eigen::EigenSolver<Eigen::MatrixXd> eigenSolver;
    eigenSolver.compute(Am1);
    return std::make_pair(eigenSolver.eigenvectors(), eigenSolver.eigenvalues());
}


