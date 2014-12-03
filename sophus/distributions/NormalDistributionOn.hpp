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
#include <vector>
#include <limits>

#include <Eigen/Eigen>
#include <Eigen/Geometry>

#include <Distributions.hpp>

namespace Sophus {
namespace Distributions {

/**  Normal Distribution on any Lie Group implemented in Sophus.
  *  NOTE: The covariance is defined in the exponential coordinates!
  *
  *  See also Timothy D. Barfoot and Paul T. Furgale, "Associating 
  *  Uncertainty With Three-Dimensional Poses for Use in Estimation
  *  Problems", IEEE TRANSACTIONS ON ROBOTICS, VOL. 30, NO. 3, JUNE 2014
  *  and
  *  Yunfeng Wang and Gregory S. Chirikjian, "Nonparametric Second-order
  *  Theory of Error Propagation on Motion Groups", The International 
  *  Journal of Robotics Research November/December 2008 27: 1258-1273
  */
template < typename GroupT >
class NormalDistributionOn {
	
public:
	typedef GroupT Group;
	
	typedef typename Group::Scalar Scalar;
	typedef typename Group::Transformation Transformation;
	typedef typename Group::Tangent Tangent;
	typedef typename Group::Adjoint Adjoint;
	
	typedef Adjoint Covariance;
	
	static const int DoF = Group::DoF;
	
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	
	
	// covariance lives in tangent space!
	NormalDistributionOn( const Tangent& mean, const Covariance& cov ) : mMean(Group::exp(mean)), mCovarianceMatrix(cov) { init(); }
	NormalDistributionOn( const Group& mean, const Covariance& cov ) : mMean(mean), mCovarianceMatrix(cov) { init(); }
	
	inline const Group& mean() const { return mMean; }
	inline const Group& mean_inverse() const { return nMeanInv; }
	
	inline const Covariance& covariance() const { return mCovarianceMatrix; }
	inline const Covariance& information() const { return nInformationMatrix; }
	
	inline const Scalar& normalizingConstant() const { return nNormalizingConstant; }
	inline const Scalar& lnNormalizingConstant() const { return nLnNormalizingConstant; }
	
	
protected:
	Group mMean;
	
	Covariance mCovarianceMatrix;
	
	// derived
	Group nMeanInv;
	Scalar nNormalizingConstant;
	Scalar nLnNormalizingConstant;
	Covariance nInformationMatrix;
	
	
	void init() {
		using namespace std;
		
		nMeanInv = mMean.inverse();
		
		// this an approximation
		nNormalizingConstant = 1. / ( pow( sqrt( 2 * SophusConstants<Scalar>::pi() ), DoF ) * sqrt( mCovarianceMatrix.determinant() ) );
		nLnNormalizingConstant = log(nNormalizingConstant);
		
		// precompute for later, maybe explicitly use Eigen::ColPivHouseholderQR (R.transpose() is cholesky factor)
		nInformationMatrix = mCovarianceMatrix.inverse();
	}
};





template < typename Group >
Group argmax( const NormalDistributionOn<Group>& N ) {
	return N.mean();
}

template < typename Group >
size_t total_number_of_components(const NormalDistributionOn<Group>& N ) {
	return 1;
}


// probability with other group scalar type (e.g. ceres::jet)
template< typename GroupA, typename GroupB >
typename GroupB::Scalar pdf_( const NormalDistributionOn<GroupA>& N, const GroupB& x ) {
	typename NormalDistributionOn<GroupB>::Covariance inf;
	for(int i=0; i<GroupB::DoF; i++)
	for(int j=0; j<GroupB::DoF; j++) {
		inf(i,j) = typename GroupB::Scalar( N.information()(i,j) );
	}
	
	typename GroupB::Tangent mean_inv_vec;
	for(int i=0; i<GroupB::DoF; i++) {
		mean_inv_vec(i) = typename GroupB::Scalar( N.mean_inverse().log()(i) );
	}
	GroupB mean_inv = GroupB::exp(mean_inv_vec);
	
	typename GroupB::Tangent y = GroupB::log( mean_inv * x );
	return typename GroupB::Scalar(N.normalizingConstant()) * std::exp( -.5 * y.transpose() * inf * y );
}

// probability
template< typename Group >
typename Group::Scalar pdf( const NormalDistributionOn<Group>& N, const Group& x ) {
	typename Group::Tangent y = Group::log( N.mean_inverse() * x );
	return N.normalizingConstant() * std::exp( -.5 * y.transpose() * N.information() * y );
}

// log_n probability
template< typename Group >
typename Group::Scalar ln_pdf( const NormalDistributionOn<Group>& N, const Group& x ) {
	typename Group::Tangent y = Group::log( N.mean_inverse() * x );
	return N.lnNormalizingConstant()  -.5 * y.transpose() * N.information() * y ;
}

// squared mahalanobis distance
template< typename GroupA, typename GroupB >
typename GroupB::Scalar squared_mahalanobis_distance( const NormalDistributionOn<GroupA>& N, const GroupB& x ) {
	typename NormalDistributionOn<GroupB>::Covariance inf;
	for(int i=0; i<GroupB::DoF; i++)
	for(int j=0; j<GroupB::DoF; j++) {
		inf(i,j) = typename GroupB::Scalar( N.information()(i,j) );
	}
	
	GroupB mean_inv;
	for(int i=0; i<GroupB::num_parameters; i++) {
		mean_inv.data()[i] = typename GroupB::Scalar( N.mean_inverse().data()[i] );
	}
	
	typename GroupB::Tangent y = GroupB::log( mean_inv * x );
	return  y.transpose() * inf * y ;
}

/// estimate only covariance from samples, must be an stl collection
template< typename GroupCollection >
typename NormalDistributionOn< typename GroupCollection::value_type >::Covariance estimate_covariance(const typename GroupCollection::value_type& mean, const GroupCollection& samples) {
	typedef typename GroupCollection::value_type Group;

	typename NormalDistributionOn<Group>::Covariance cov = NormalDistributionOn<Group>::Covariance::Zero();
	
	Group meanInv = mean.inverse();
	
	for( typename GroupCollection::const_iterator s = samples.begin(); s != samples.end(); s++) {
		typename Group::Tangent y = Group::log( meanInv * (*s) );
		
		cov += y * y.transpose();
	}
	
	cov /= samples.size()-1;
	
	return cov ;
}

/// estimate full normal distribution from samples, must be an stl collection
template< typename GroupCollection >
NormalDistributionOn< typename GroupCollection::value_type> estimate(const GroupCollection& samples, size_t meanSteps = 3) {
	typedef typename GroupCollection::value_type Group;
	
	if( samples.size() < 1 ) {
		// maybe later throw error
		return NormalDistributionOn<Group>( Group(), NormalDistributionOn<Group>::Covariance::Zero() );
	}
	
	Group mean = samples.front();
	Group meanInv = mean.inverse();
	
	for(size_t i=0; i<meanSteps; i++) {
		typename Group::Tangent sum = Group::Tangent::Zero();
		
		for( typename GroupCollection::const_iterator s = samples.begin(); s != samples.end(); s++) {
			sum += Group::log( meanInv * (*s) );
		}
		
		sum /= samples.size();
		
		mean *= Group::exp( sum );
		meanInv = mean.inverse();
	}
	
	typename NormalDistributionOn<Group>::Covariance cov = estimate_covariance(mean, samples);
	
	return NormalDistributionOn<Group>( mean, cov );
}




}
}

