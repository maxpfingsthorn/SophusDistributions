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

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>


#include <Eigen/Dense>

#include "NormalDistributionOn.hpp"

namespace Sophus {
namespace Distributions {


// randomly sample according to distribution
template< typename LieGroup >
class SamplerOfNormalDistributionOn {
public:
	typedef LieGroup Group;
	typedef NormalDistributionOn<LieGroup> Distribution;
	
	const NormalDistributionOn<LieGroup>& N;
	
	boost::normal_distribution<> nd;
	
	
	typename NormalDistributionOn<LieGroup>::Covariance sqrt_cov;
	
	SamplerOfNormalDistributionOn(const NormalDistributionOn<LieGroup>& n) : N(n), nd(0,1) {
		
		Eigen::JacobiSVD< typename Distribution::Covariance, Eigen::NoQRPreconditioner > svd( n.covariance() , Eigen::ComputeFullU);
		
		typename Distribution::Tangent sqrt_eigs = svd.singularValues();
		for(size_t i=0; i<Distribution::DoF; i++) {
			sqrt_eigs(i) = sqrt( sqrt_eigs(i) );
		}
		
		sqrt_cov = svd.matrixU() * sqrt_eigs.asDiagonal();
	}
	
	template < typename RNG >
	Group operator() (RNG& rng) {
		boost::variate_generator< RNG&, boost::normal_distribution<> > rnd_gen(rng, nd);
		
		typename Distribution::Tangent random;
		
		for( size_t i=0; i<Distribution::DoF; i++) {
			random(i) = rnd_gen();
		}
		
		typename Distribution::Tangent t = sqrt_cov * random;
		
		return N.mean() * Group::exp(t);
	}	
	
};

template < typename LieGroup >
struct SampleTraits< NormalDistributionOn<LieGroup> > {
	static const bool supports_sampling = true;
	
	typedef SamplerOfNormalDistributionOn<LieGroup> Sampler;
};



}
}


