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

#include <Distributions.hpp>

namespace Sophus {
namespace Distributions {

template < typename Dist >
class MixtureOf {
public:
	static const int DoF = Dist::DoF;

	typedef typename Dist::Scalar Scalar;
	typedef typename Dist::Group Group;
	
	typedef Dist ComponentDistribution;
	
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	
	MixtureOf() : mWeightsDoSumToOne(false) {}
	
	inline size_t numComponents() const { return mComponents.size(); }
	
	inline const Dist& component(size_t i) const { return mComponents[i]; }
	inline double weight(size_t i) const { return mWeights[i]; }
	
	inline const std::vector< typename Dist::Scalar >& weights() const { return mWeights; }
	
	inline void normalizeWeights() {
		typename Dist::Scalar sum = sumOfWeights();
		for(size_t i = 0; i < mWeights.size(); i++) {
			mWeights[i] /= sum;
		}
		mWeightsDoSumToOne = true;
	}

	inline typename Dist::Scalar sumOfWeights() const {
		typename Dist::Scalar sum(0);
		for(size_t i = 0; i < mWeights.size(); i++) {
			sum += mWeights[i];
		}
		return sum;
	}

	inline bool weightsDoSumToOne() const {
		return mWeightsDoSumToOne;
	}
	
	void addComponent(const typename Dist::Scalar& w, const Dist& d) {
		mComponents.push_back(d);
		mWeights.push_back(w);
	}
	
	
protected:
	std::vector< Dist > mComponents;
	std::vector< typename Dist::Scalar > mWeights;
	bool mWeightsDoSumToOne;
};

template< typename Dist >
typename Dist::Group argmax( const MixtureOf<Dist>& N ) {
	typename Dist::Scalar w(0);
	typename Dist::Group g;
	
	for(size_t i = 0; i < N.numComponents(); i++) {
		if( i == 0 || w < N.weight(i) ) {
			w = N.weight(i);
			g = argmax( N.component(i) );
		}
	}
	
	return g;
}

template< typename Dist >
size_t max_component( const MixtureOf<Dist>& N ) {
	typename Dist::Scalar max_v=0;
	size_t max_i=0;
	
	for(size_t i = 0; i < N.numComponents(); i++) {
		typename Dist::Scalar pr = N.weight(i);
		if( i==0 || pr > max_v ) {
			max_v = pr;
			max_i = i;
		}
	}
	
	return max_i;
}

template< typename Dist >
size_t max_component( const MixtureOf<Dist>& N, const typename Dist::Group& x ) {
	typename Dist::Scalar max_v=0;
	size_t max_i=0;
	
	for(size_t i = 0; i < N.numComponents(); i++) {
		typename Dist::Scalar pr = ln_pdf( N.component(i), x);
		if( i==0 || pr > max_v ) {
			max_v = pr;
			max_i = i;
		}
	}
	
	return max_i;
}

template < typename Dist >
size_t total_number_of_components(const MixtureOf<Dist>& N ) {
	size_t n=0;
	for(size_t i = 0; i < N.numComponents(); i++) {
		n += total_number_of_components( N.component(i) );
	}
	if ( !N.weightsDoSumToOne() ) {
		return n+1;
	}
	return n;
}


// probability
template< typename Dist, typename Group >
typename Group::Scalar pdf_( const MixtureOf<Dist>& N, const Group& x ) {
	if( N.numComponents() == 1 ) {
		return pdf_(N.component(0), x);
	}
	
	typename Group::Scalar p(0);
	for(size_t i = 0; i < N.numComponents(); i++) {
		p += N.weight(i) * pdf_( N.component(i), x );
	}
	
	return p;
}

template< typename Dist >
typename Dist::Scalar pdf( const MixtureOf<Dist>& N, const typename Dist::Group& x ) {
	if( N.numComponents() == 1 ) {
		return pdf(N.component(0), x);
	}
	
	typename Dist::Scalar p(0);
	for(size_t i = 0; i < N.numComponents(); i++) {
		p += N.weight(i) * pdf( N.component(i), x );
	}
	
	return p;
}

template< typename Dist >
typename Dist::Scalar ln_pdf( const MixtureOf<Dist>& N, const typename Dist::Group& x ) {
	if( N.numComponents() == 1 ) {
		return ln_pdf(N.component(0), x);
	}
	
	typename Dist::Scalar p = pdf(N, x);
	
	if( p > Constants<typename Dist::Scalar>::epsilon() ) {
		return log(p);
	}

	// here, we are far away from any mode!
	
	// approximate with "max-mixture" model
	p = Constants<typename Dist::Scalar>::lowest();
	for(size_t i = 0; i < N.numComponents(); i++) {
		typename Dist::Scalar pr = ln_pdf( N.component(i), x);
		if( i==0 || pr > p ) {
			p = pr;
		}
	}
	
	return p;
}
	
}
}

