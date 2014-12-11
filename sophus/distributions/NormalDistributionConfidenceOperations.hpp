/*
* Copyright (c) 2014 Jacobs University Robotics Group
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without modification, 
* are strictly limited to non-commercial academic purposes and are only permitted 
* provided that the following conditions are met:
*     * Redistributions of source code must retain the above copyright
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above copyright
*       notice, this list of conditions and the following disclaimer in the
*       documentation and/or other materials provided with the distribution.
*     * Neither the name of Jacobs University Bremen nor the name jacobs_robotics nor the
*       names of its contributors may be used to endorse or promote products
*       derived from this software without specific prior written permission.
*
* If you are interested in using this code for a commercial or non-academic purpose, 
* please contact us.
*
* THIS SOFTWARE IS PROVIDED BY Jacobs Robotics ``AS IS'' AND ANY
* EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL Jacobs Robotics BE LIABLE FOR ANY
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
* Contact: Max Pfingsthorn, m.pfingsthorn@jacobs-university.de
*          Andreas Birk, a.birk@jacobs-university.de
*
* Paper Mail:
*
*   Jacobs University Bremen gGmbH
*   Andreas Birk
*   Campus Ring 1
*   D-28759 Bremen
*   Germany
*/

#pragma once

#include <cmath>
#include <list>

#include <Eigen/Eigenvalues>

#include <boost/math/distributions/chi_squared.hpp>

#include <sophus/sophus.hpp>

#include "NormalDistributionOn.hpp"



namespace Sophus {
namespace Distributions {

template < typename Scalar, int Dim >
std::list< Eigen::Matrix<Scalar,Dim,1> > sample_sphere(int steps, int d) {
	using namespace std;
	std::list< Eigen::Matrix<Scalar,Dim,1> > out;
	
	double inc = 2*Sophus::SophusConstants<Scalar>::pi()/(steps-1);
	
		int d1 = d+1;
		int d2 = d+2;
		if( d1 >= Dim ) {
			d1-=Dim;
		}
		if( d2 >= Dim ) {
			d2-=Dim;
		}
		
		Eigen::Matrix<Scalar,Dim,1> s = Eigen::Matrix<Scalar,Dim,1>::Zero();
		Scalar a=0;
		for( int i=0; i < steps; i++, a+=inc) {
			s(d1) = cos(a);
			s(d2) = sin(a);
			out.push_back( s );
		}

	return out;
}

template < typename LieGroup >
std::list< LieGroup > confidence_region_contours(const NormalDistributionOn<LieGroup>& N, int dim, typename LieGroup::Scalar confidence) {
	using namespace std;
	
	list< typename LieGroup::Tangent > sphere = sample_sphere< typename LieGroup::Scalar, LieGroup::DoF >(50, dim);
	
	boost::math::chi_squared_distribution<typename LieGroup::Scalar> chi2( LieGroup::DoF );
	double scale = sqrt( boost::math::quantile(chi2, confidence) ) ;
	
	Eigen::SelfAdjointEigenSolver< typename NormalDistributionOn< LieGroup >::Covariance > eigs;
	eigs.compute(N.covariance());
	
	typename NormalDistributionOn<SE2d>::Covariance sqrt_cov = eigs.eigenvectors() * eigs.eigenvalues().array().sqrt().matrix().asDiagonal();
	
	std::list< LieGroup > out;
	
	for( typename list< typename LieGroup::Tangent >::iterator it = sphere.begin(); it != sphere.end(); it++ ) {
		out.push_back( N.mean() * LieGroup::exp( scale* sqrt_cov * (*it) ) );
	}
	
	return out;
	
}

template < typename LieGroup >
bool inside_confidence_region(const NormalDistributionOn<LieGroup>& N, const LieGroup& x, typename LieGroup::Scalar confidence) {
	return squared_mahalanobis_distance(N,x) < boost::math::quantile( boost::math::chi_squared_distribution<typename LieGroup::Scalar>( NormalDistributionOn<LieGroup>::DoF ), confidence);
}

}
}
