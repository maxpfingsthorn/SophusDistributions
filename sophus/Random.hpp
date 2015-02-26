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

#include <iostream>
#include <vector>
#include <algorithm>
#include <boost/random.hpp>


namespace Sophus {
namespace Random {

template < typename LieGroup >
class uniform_group {
public:
	typedef typename LieGroup::Scalar input_type;
	typedef LieGroup result_type;

	explicit uniform_group(const typename LieGroup::Point& min_v, const typename LieGroup::Point& max_v)
		: _translation_dim(min_v.innerSize()), _rotation_dim(LieGroup::num_parameters - _translation_dim), _rotation_random(_rotation_dim)
	{
		init_translation_random(min_v, max_v);
	}
	explicit uniform_group(const typename LieGroup::Point& abs_v) 
		: _translation_dim(abs_v.innerSize()), _rotation_dim(LieGroup::num_parameters - _translation_dim), _rotation_random(_rotation_dim)
	{
		init_translation_random(-abs_v, abs_v);
	}

	template < typename Engine >
	const result_type& operator()(Engine& eng) {
		typedef  typename LieGroup::Scalar Scalar;

		// construct uniform rotation
		boost::uniform_01< Scalar > uni;
		Scalar _r1(0), _r2(0), _cached_rho(0);
		bool _valid=false;

		using std::sqrt; using std::log; using std::sin; using std::cos;
		const Scalar pi = Scalar(3.14159265358979323846);

		Scalar sqsum = 0;
		for(int i=0; i<_rotation_dim; i++) {
			if(!_valid) {
				_r1 = uni(eng);
				_r2 = uni(eng);
				_cached_rho = sqrt(-Scalar(2) * log(Scalar(1)-_r2));
				_valid = true;
			} else {
				_valid = false;
			}

			Scalar n = _cached_rho * (_valid ? cos(Scalar(2)*pi*_r1) : sin(Scalar(2)*pi*_r1));
 			_group.data()[i] = n;

 			sqsum += n*n;

			//std::cerr << "_g.data()[" << i << "] = " << _group.data()[i] << ", n = " << n << std::endl;
		}

		sqsum = sqrt(sqsum);

		// if we made a zero vector, make it identity
		if( sqsum < std::numeric_limits< Scalar >::epsilon() ) {
			for(int i=0; i<_rotation_dim; i++) {
				_group.data()[i] = Scalar(0);
				if( i == _rotation_dim -1 ) {
					_group.data()[i] = Scalar(1);
				}
			}
		} else { // normalize rotation
			// if we made a quaternion, make sure last element is positive
			if( _rotation_dim == 4 && _group.data()[_rotation_dim-1] < 0 ) {
				sqsum = -sqsum;
			}
			for(int i=0; i<_rotation_dim; i++) {
				_group.data()[i] /= sqsum;
			}
		}

		// construct uniform (within bounds) translation
		for(int i=0; i<_translation_dim; i++) {
			_group.translation()[i] = _translation_random[i](eng);
		}

		// don't have to normalize group representation here, already done above

		return _group;
	}
private:
	int _translation_dim;
	int _rotation_dim;

	result_type _group;

	inline void init_translation_random(const typename LieGroup::Point& min_v, const typename LieGroup::Point& max_v) {
		//std::cerr << "group has " << _translation_dim << " translation and " << _rotation_dim << " rotation dims" << std::endl;

		for(int i=0; i<_translation_dim; i++) {
			_translation_random.push_back( boost::uniform_real< typename LieGroup::Scalar >(min_v[i], max_v[i]) );
			//std::cerr << "dim " << i << ": " << min_v[i] << " < " << max_v[i] << std::endl;
		}
	}

	std::vector< boost::uniform_real< typename LieGroup::Scalar > > _translation_random;
	boost::uniform_on_sphere< typename LieGroup::Scalar > _rotation_random;
};


}}

