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

//#include <ceres/jet.h>
#include <sophus.hpp>

namespace NormalDistributionConversions {
	

	
template < typename Group, typename Derived1, typename Derived2 >
void  project_covariance_to_exponential_coordinates(const Group& g, const Eigen::MatrixBase<Derived1>& in_cov, Eigen::MatrixBase<Derived2>& out_cov) {
	
	
	Eigen::Matrix< typename Group::Scalar, Group::DoF, Group::num_parameters > jacobian;
	
	/*  // autodiff somehow does not work, gives huge values for rotation part
	{
		using namespace ceres;
	
		typedef Jet< typename Group::Scalar, Group::num_parameters> ActiveScalar;
		typedef typename Sophus::cast_traits<Group>::template to< ActiveScalar >::Type ActiveGroup;
	
		ActiveGroup ag = g.template cast< ActiveScalar >();
		// set infinitisimal parts
		for(int i=0; i<Group::num_parameters; i++) {
			ag.data()[i].v[i] = 1;
			//std::cout << "input row " << i << ": " << ag.data()[i].v.transpose() << std::endl;
		}
	
	
	
		typename ActiveGroup::Tangent t = ag.log();
	
		
	
	
		for(int i=0; i<Group::DoF; i++) {
			//std::cout << "auto jacobian row " << i << ": " << t[i].v.transpose() << std::endl;
			jacobian.row(i) = t[i].v.transpose();
		}
	}
	std::cout << "au jacobian: " << std::endl << jacobian << std::endl;
	//*/
	
	
	// for testing: numeric diff
	//*
	{
		typename Group::Scalar h (1e-7);
		Group fd_g(g);
		for(int i=0; i<Group::num_parameters; i++) {
			typename Group::Scalar old_v = fd_g.data()[i];
			fd_g.data()[i]= old_v + h;
			fd_g.normalize();
			typename Group::Tangent t = fd_g.log();
			
			fd_g.data()[i]= old_v - h;
			fd_g.normalize();
			
			t -= fd_g.log();
			t /= 2*h;
			
			fd_g.data()[i] = old_v;
			
			jacobian.col(i) = t;
		}
	}
	//*/
	
	out_cov = jacobian * in_cov * jacobian.transpose();
	
}

template < typename Group, typename Derived1, typename Derived2 >
void  project_covariance_to_logarithmic_coordinates(const Group& g, const Eigen::MatrixBase<Derived1>& in_cov, Eigen::MatrixBase<Derived2>& out_cov) {
	
	typename Group::Tangent t = g.log();
	
	Eigen::Matrix< typename Group::Scalar, Group::num_parameters, Group::DoF > jacobian;
	

	//*
	{
		typename Group::Scalar h (1e-7);

		for(int i=0; i<Group::DoF; i++) {
			Eigen::Matrix< typename Group::Scalar, Group::num_parameters, 1 > fd_g_ph, fd_g_mh, fd_g;
			Eigen::Map< Group > fd_g_ph_map(fd_g_ph.data());
			Eigen::Map< Group > fd_g_mh_map(fd_g_mh.data());
			
			typename Group::Scalar old_v = t[i];
			
			t[i]= old_v + h;
			fd_g_ph_map = Group::exp(t);
			
			t[i]= old_v - h;
			fd_g_mh_map = Group::exp(t);
			
			//fd_g_map = g_t_mh.inverse() * g_t_ph; // this is basically (g_t_ph - g_t_mh)
			fd_g = fd_g_ph - fd_g_mh;
				
			//std::cout << "col " << i << ": " << fd_g.transpose() << std::endl;
			
			fd_g /= 2*h;
			
			t[i] = old_v;
			
			jacobian.col(i) = fd_g;
		}
	}
	//*/
	
	
	out_cov = jacobian * in_cov * jacobian.transpose();
	
}

}

