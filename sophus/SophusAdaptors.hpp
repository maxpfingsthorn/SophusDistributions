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

#include <string>
#include <stream>

#include <Eigen/Eigen>
#include <Eigen/Geometry>

#include <se3.hpp>
#include <se2.hpp>

// this code could be improved by size checks, both runtime and compiletime with enable_if

namespace Sophus {
namespace Adaptors {


template < typename Derived1, typename Derived2 >
void convert_from_x_y_theta(const Eigen::MatrixBase<Derived1>& x_y_theta, Sophus::SE2GroupBase<Derived2>& g, Eigen::Matrix< typename Derived2::Scalar, Derived2::num_parameters, 3>& jacobian) {
	Eigen::Matrix< typename Eigen::internal::traits<Derived1>::Scalar, Derived2::num_parameters, 1 > conv;
	// we know Derived2::num_parameters == 4
	conv[0] = std::cos(x_y_theta[2]);
	conv[1] = std::sin(x_y_theta[2]);
	conv[2] = x_y_theta[0];
	conv[3] = x_y_theta[1];
	
	g = Eigen::Map< Sophus::SE2Group< typename Eigen::internal::traits<Derived1>::Scalar > >(conv.data());
	
	typename Derived2::Scalar zero(0);
	typename Derived2::Scalar one (1);
	typename Derived2::Scalar s(std::sin(x_y_theta[2]));
	typename Derived2::Scalar c(std::cos(x_y_theta[2]));
	
	jacobian << zero, zero,   -s,
	            zero, zero,    c,
				 one, zero, zero,
				zero,  one, zero;
}

template < typename Derived1, typename Derived2 >
void convert_to_x_y_theta(const Sophus::SE2GroupBase<Derived2>& g, Eigen::MatrixBase<Derived1>& x_y_theta, Eigen::Matrix< typename Derived2::Scalar, 3, Derived2::num_parameters>& jacobian) {
	x_y_theta.head(2) = g.translation();
	x_y_theta[2] = atan2(g.so2().unit_complex().y(), g.so2().unit_complex().x());
	
	
	typename Derived2::Scalar zero(0);
	typename Derived2::Scalar one (1);
	typename Derived2::Scalar x(g.so2().unit_complex().x());
	typename Derived2::Scalar y(g.so2().unit_complex().y());
	
	jacobian << zero, zero,  one, zero,
	            zero, zero, zero,  one,
				  -y,    x, zero, zero;
}

template < typename Derived1, typename Derived2 >
void convert_from_x_y_z_qw_qx_qy_qz(const Eigen::MatrixBase<Derived1>& x_y_z_qw_qx_qy_qz, Sophus::SE3GroupBase<Derived2>& g, Eigen::Matrix< typename Derived2::Scalar, Derived2::num_parameters, 7>& jacobian) {
	
	typedef typename Eigen::internal::traits<Derived1>::Scalar Derived1Scalar;
	
	g = Sophus::SE3Group< Derived1Scalar > (
		Eigen::Quaternion< Derived1Scalar >(
			x_y_z_qw_qx_qy_qz[3], x_y_z_qw_qx_qy_qz[4], x_y_z_qw_qx_qy_qz[5], x_y_z_qw_qx_qy_qz[6]
		),
		x_y_z_qw_qx_qy_qz.head(3)
	);
	
	typename Derived2::Scalar zero(0);
	typename Derived2::Scalar one (1);
	
	// its rotation, then translation
	jacobian << zero, zero, zero, zero,  one, zero, zero,
	            zero, zero, zero, zero, zero,  one, zero,
				zero, zero, zero, zero, zero, zero,  one,
				zero, zero, zero,  one, zero, zero, zero,
				 one, zero, zero, zero, zero, zero, zero,
				zero,  one, zero, zero, zero, zero, zero,
				zero, zero,  one, zero, zero, zero, zero;
}


template < typename Derived1, typename Derived2 >
void convert_from_x_y_z_qx_qy_qz_qw(const Eigen::MatrixBase<Derived1>& x_y_z_qx_qy_qz_qw, Sophus::SE3GroupBase<Derived2>& g, Eigen::Matrix< typename Derived2::Scalar, Derived2::num_parameters, 7>& jacobian) {
	
	g = Sophus::SE3Group< typename Eigen::internal::traits<Derived1>::Scalar > (
		Eigen::Quaternion< typename Eigen::internal::traits<Derived1>::Scalar >(
			x_y_z_qx_qy_qz_qw[6], x_y_z_qx_qy_qz_qw[3], x_y_z_qx_qy_qz_qw[4], x_y_z_qx_qy_qz_qw[5]
		),
		x_y_z_qx_qy_qz_qw.head(3)
	);
	
	
	typename Derived2::Scalar zero(0);
	typename Derived2::Scalar one (1);
	
	// its rotation, then translation
	jacobian << zero, zero, zero,  one, zero, zero, zero,
	            zero, zero, zero, zero,  one, zero, zero,
				zero, zero, zero, zero, zero,  one, zero,
				zero, zero, zero, zero, zero, zero,  one,
				 one, zero, zero, zero, zero, zero, zero,
				zero,  one, zero, zero, zero, zero, zero,
				zero, zero,  one, zero, zero, zero, zero;
}


template < typename Derived1, typename Derived2 >
void convert_to_x_y_z_qx_qy_qz_qw(const Sophus::SE3GroupBase<Derived2>& g, Eigen::MatrixBase<Derived1>& x_y_z_qx_qy_qz_qw, Eigen::Matrix< typename Derived2::Scalar, 7, Derived2::num_parameters>& jacobian) {
	
	
	x_y_z_qx_qy_qz_qw.head(3) = g.translation();
	x_y_z_qx_qy_qz_qw.tail(4) = g.unit_quaternion().coeffs();
	
	typename Derived2::Scalar zero(0);
	typename Derived2::Scalar one (1);
	
	// its rotation, then translation
	jacobian << zero, zero, zero, zero,  one, zero, zero,
	            zero, zero, zero, zero, zero,  one, zero,
				zero, zero, zero, zero, zero, zero,  one,
				 one, zero, zero, zero, zero, zero, zero,
				zero,  one, zero, zero, zero, zero, zero,
				zero, zero,  one, zero, zero, zero, zero,
				zero, zero, zero,  one, zero, zero, zero;
}

template < typename Derived1, typename Derived2 >
void convert_from_x_y_z_qx_qy_qz_of_g2o(const Eigen::MatrixBase<Derived1>& x_y_z_qx_qy_qz, Sophus::SE3GroupBase<Derived2>& g, Eigen::Matrix< typename Derived2::Scalar, Derived2::num_parameters, 6>& jacobian) {
	
	Eigen::Matrix< typename Eigen::internal::traits<Derived1>::Scalar, 4, 1> quat_mat;
	quat_mat.head(3) = x_y_z_qx_qy_qz.tail(3);

	typename Eigen::internal::traits<Derived1>::Scalar w2 = typename Eigen::internal::traits<Derived1>::Scalar(1.)-x_y_z_qx_qy_qz.tail(3).squaredNorm();
	if( w2 < 0 ) {
		throw std::invalid_argument("Quaternion not properly scaled!");
	}
	quat_mat[3] = sqrt(w2);
	
	
	g = Sophus::SE3Group< typename Eigen::internal::traits<Derived1>::Scalar >(
		Eigen::Quaternion< typename Eigen::internal::traits<Derived1>::Scalar >(
			quat_mat
		),
		x_y_z_qx_qy_qz.head(3)
	);
	
	typename Derived2::Scalar zero(0);
	typename Derived2::Scalar one (1);
	

	typename Derived2::Scalar dx_dw(0);
	typename Derived2::Scalar dy_dw(0);
	typename Derived2::Scalar dz_dw(0);
	
	if(quat_mat[3]>1e-10) {
		dx_dw = typename Derived2::Scalar(-quat_mat[0]/quat_mat[3]);
		dy_dw = typename Derived2::Scalar(-quat_mat[1]/quat_mat[3]);
		dz_dw = typename Derived2::Scalar(-quat_mat[2]/quat_mat[3]);
	}
	
	// its rotation, then translation
	jacobian << zero, zero, zero,   one,  zero,  zero,
	            zero, zero, zero,  zero,   one,  zero,
				zero, zero, zero,  zero,  zero,   one,
				zero, zero, zero, dx_dw, dy_dw, dz_dw,
				 one, zero, zero,  zero,  zero,  zero,
				zero,  one, zero,  zero,  zero,  zero,
				zero, zero,  one,  zero,  zero,  zero;
}

template < typename Derived1, typename Derived2 >
void convert_from_x_y_z_qx_qy_qz_qw_of_g2o(const Eigen::MatrixBase<Derived1>& x_y_z_qx_qy_qz_qw, Sophus::SE3GroupBase<Derived2>& g, Eigen::Matrix< typename Derived2::Scalar, Derived2::num_parameters, 6>& jacobian) {
	
	Eigen::Matrix< typename Eigen::internal::traits<Derived1>::Scalar, 4, 1> quat_mat;
	quat_mat = x_y_z_qx_qy_qz_qw.tail(4);
	
	g = Sophus::SE3Group< typename Eigen::internal::traits<Derived1>::Scalar >(
		Eigen::Quaternion< typename Eigen::internal::traits<Derived1>::Scalar >( quat_mat ),
		x_y_z_qx_qy_qz_qw.head(3)
	);
	
	typename Derived2::Scalar zero(0);
	typename Derived2::Scalar one (1);
	

	typename Derived2::Scalar dx_dw(0);
	typename Derived2::Scalar dy_dw(0);
	typename Derived2::Scalar dz_dw(0);
	
	// as if we computed qw
	if(quat_mat[3]>1e-10) {
		dx_dw = typename Derived2::Scalar(-quat_mat[0]/quat_mat[3]);
		dy_dw = typename Derived2::Scalar(-quat_mat[1]/quat_mat[3]);
		dz_dw = typename Derived2::Scalar(-quat_mat[2]/quat_mat[3]);
	}
	
	// its rotation, then translation
	jacobian << zero, zero, zero,   one,  zero,  zero,
	            zero, zero, zero,  zero,   one,  zero,
				zero, zero, zero,  zero,  zero,   one,
				zero, zero, zero, dx_dw, dy_dw, dz_dw,
				 one, zero, zero,  zero,  zero,  zero,
				zero,  one, zero,  zero,  zero,  zero,
				zero, zero,  one,  zero,  zero,  zero;
}

}
}


