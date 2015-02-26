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

#include <Eigen/Dense>

namespace Sophus {
namespace Distributions {
namespace math_helper {

template < typename Derived1, typename Derived2 >
inline bool matrix_square_root_svd(const Eigen::MatrixBase<Derived1>& in, Eigen::MatrixBase<Derived2>& out) {
		typedef typename Eigen::JacobiSVD< Derived1, Eigen::NoQRPreconditioner > SVD;

		SVD svd( in , Eigen::ComputeFullU);
		
		typename SVD::SingularValuesType sqrt_eigs = svd.singularValues();
		for(size_t i=0; i<sqrt_eigs.innerSize(); i++) {
			sqrt_eigs(i) = sqrt( sqrt_eigs(i) );
		}
		
		out = svd.matrixU() * sqrt_eigs.asDiagonal();
}

}}}
