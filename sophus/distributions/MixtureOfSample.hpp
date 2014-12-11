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

namespace Sophus {
namespace Distributions {

template < typename Dist >
class MixtureOfSample {
public:
	typedef typename Dist::Group Group;
	typedef typename SampleTraits<Dist>::Sampler ComponentSampler;
	
	const MixtureOf<Dist>& distribution;
	
	boost::uniform_real<typename Dist::Scalar> rnd_real;
	std::vector< ComponentSampler* > comp_samplers;
	
	MixtureOfSample( const MixtureOf<Dist>& d ) : distribution(d), rnd_real(0,1) {
		for(size_t i=0; i<d.numComponents(); i++) {
			comp_samplers.push_back( 
					new ComponentSampler(d.component(i)) 
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
		typename Group::Scalar d = rnd_real(rng);
		typename Group::Scalar prob = 0;
		for(size_t i=0; i<distribution.numComponents(); i++) {
			prob += distribution.weight(i);
			if( d <= prob ) {
				return (*comp_samplers[i])(rng);
			}
		}
		return (*comp_samplers.back())(rng);
	}
};


template < typename Dist >
struct SampleTraits< MixtureOf<Dist> > {
	static const bool supports_sampling = SampleTraits<Dist>::supports_sampling;
	
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

