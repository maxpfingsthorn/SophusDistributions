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


namespace Sophus {
namespace Distributions {

template < typename Dist >
class MixtureOfSample {
public:
	typedef typename Dist::Group Group;
	typedef typename DistributionSampleTraits<Dist>::Sampler DistSampler;
	
	const MixtureOf<Dist>& distribution;
	
	boost::random::discrete_distribution<size_t, typename Dist::Scalar> rnd_dist;
	std::vector< DistSampler* > comp_samplers;
	
	MixtureOfSample( const MixtureOf<Dist>& d ) : distribution(d), rnd_dist( distribution.weights() ) {
		for(size_t i=0; i<d.numComponents(); i++) {
			comp_samplers.push_back( 
					new DistSampler(d.component(i)) 
			);
		}
	}
	
	~MixtureOfSample() {
		for(size_t i = 0; i < comp_samplers.size(); ++i)
		{
			delete comp_samplers[i];
		}
		comp_samplers.clear();
	}
	
	template < typename RNG >
	Group operator() (RNG& rng) {
		size_t comp = rnd_dist(rng);
		return (*comp_samplers[comp])(rng);
	}
};


template < typename Dist >
struct DistributionSampleTraits< MixtureOf<Dist> > {
	static const bool supports_sampling = DistributionSampleTraits<Dist>::supports_sampling;
	
	typedef MixtureOfSample<Dist> Sampler;
};

template < typename Dist >
bool inside_confidence_region(const MixtureOf<Dist>& M, const typename Dist::Group& x, typename Dist::Scalar confidence) {
	for(size_t i=0; i<M.numComponents(); i++) {
		if( inside_confidence_region(M.component(i), x, confidence ) ) { // approximation
			return true;
		}
	}
	return false;
}


}
}

