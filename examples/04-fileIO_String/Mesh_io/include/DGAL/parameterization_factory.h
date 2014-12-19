// Copyright (c) 2007-2010  Dalian University of Technology (China). 
// All rights reserved.
//
// This file is part of DGAL; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: http://jjcao1231.googlepages.com $
// $Id: Parameterization_factory.h 2007-04-03$
//
//
// Author(s)     : JJCAO <jjcao@hotmail.com>

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//Warning: If you meet unintelligible compile error when you include this file, 
//         you may try to include this file following #include <CGAL/Cartesian.h> immediately!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifndef DGAL_PARAMETERIZATION_FACTORY_H
#define DGAL_PARAMETERIZATION_FACTORY_H
//for cgal
#include <CGAL/Timer.h>
#include <CGAL/parameterize.h>
#include <CGAL/Parameterization_mesh_patch_3.h>
#include <CGAL/Circular_border_parameterizer_3.h>
#include <CGAL/Square_border_parameterizer_3.h>
#include <CGAL/Two_vertices_parameterizer_3.h>
#include <CGAL/Barycentric_mapping_parameterizer_3.h>
#include <CGAL/Discrete_conformal_map_parameterizer_3.h>
#include <CGAL/Discrete_authalic_parameterizer_3.h>
#include <CGAL/Mean_value_coordinates_parameterizer_3.h>
#include <CGAL/LSCM_parameterizer_3.h>
#include <CGAL/Parameterization_mesh_feature_extractor.h>

#include <CGAL/OpenNL/linear_solver.h>
#ifdef CGAL_USE_TAUCS
    #include <CGAL/Taucs_solver_traits.h>
#endif
//#include <taucs.h>

#include <DGAL/config.h>
#include <DGAL/Planar_parameterizer_3.h>
#include <string>
#include <sstream>
#include <map>

DGAL_BEGIN_NAMESPACE

enum Parameterize_type {
	UNKNOW_PARAMETERIZATION = -1,
	NONE_PARAMETERIZED = 0,
	CONFORMAL_CHORD_CIRCLE = 1,
	FLOATER_CHORD_CIRCLE ,
	BARYCENTRIC_CHORD_CIRCLE ,
	AUTHALIC_CHORD_CIRCLE ,

	CONFORMAL_CHORD_SQUARE = 21,
	FLOATER_CHORD_SQUARE ,
	BARYCENTRIC_CHORD_SQUARE ,
	AUTHALIC_CHORD_SQUARE ,
	LSCM_2PTS
};

class Parameterize_option
{
private:
	std::string type_;
	std::string border_;
	std::string solver_ ;
	bool is_parameterized_;
	std::map<std::string, Parameterize_type> type_map_;
public:
	void type(const char* in){type_ = in;}	
	void border(const char* in){border_ = in;}
	std::string type(){return type_;}
	std::string border(){return border_;}
	std::string solver(){return solver_;}

	void is_parameterized(bool in){is_parameterized_ = in;}
	bool is_parameterized(){ return is_parameterized_;}
public:
	Parameterize_option():is_parameterized_(false){init();}
	Parameterize_option(std::string type, std::string border, std::string solver= "opennl")
		:type_(type), border_(border), solver_(solver),is_parameterized_(false)
	{
		init();
	}
public:
	Parameterize_type paremeterize_type()
	{
		std::map<std::string, Parameterize_type>::iterator it = type_map_.find(type_+border_);
		if ( it == type_map_.end())
			return UNKNOW_PARAMETERIZATION;
		
		return it->second;
	}
	bool is_circle_border()
	{
		Parameterize_type pt = paremeterize_type();		
		if ( pt > 0 && pt < 21)
			return true;
		return false;
	}
	bool is_square_border()
	{
		Parameterize_type pt = paremeterize_type();		
		if ( pt > 20 && pt < 41)
			return true;
		return false;
	}
private:
	void init(){
		if ( type_map_.size() == 0)
		{
			{
			std::stringstream ss;
			ss << "conformal" << "_" << "circle";
			type_map_.insert( std::make_pair(ss.str(),CONFORMAL_CHORD_CIRCLE));
			}

			{
			std::stringstream ss;
			ss << "conformal" << "_" << "square";
			type_map_.insert( std::make_pair(ss.str(),CONFORMAL_CHORD_SQUARE));
			}

			{
			std::stringstream ss;
			ss << "authalic" << "_" << "square";
			type_map_.insert( std::make_pair(ss.str(),AUTHALIC_CHORD_SQUARE));
			}
		}
	}
};

template <class ParameterizationMesh_3>
typename CGAL::Parameterizer_traits_3<ParameterizationMesh_3>::Error_code
constrained_parameterize(ParameterizationMesh_3& mesh,   // Mesh parameterization adaptor
             const char *type="dcp",      // type of parameterization (see usage)
             const char *border="circle")    // type of border parameterization (see usage)
{
	typedef typename OpenNL::DefaultLinearSolverTraits<typename ParameterizationMesh_3::NT>
		                                SolverPolicy;
	typedef typename Hard_constrian_policy<ParameterizationMesh_3,
										   SolverPolicy::Vector,
										   SolverPolicy::Matrix
	                                    >
										ConstrainPolicy;

	return parameterize<ParameterizationMesh_3, 
		                SolverPolicy, 
						ConstrainPolicy>(mesh, type, border);
}

template <class ParameterizationMesh_3>
typename CGAL::Parameterizer_traits_3<ParameterizationMesh_3>::Error_code
parameterize(ParameterizationMesh_3& mesh,   // Mesh parameterization adaptor
             const char *type="mvc",      // type of parameterization (see usage)
             const char *border="circle", // type of border parameterization (see usage)
			 const char* solver="opennl")  // type of solver, opennl or taucs  
{
	typedef CGAL::Parameterizer_traits_3<ParameterizationMesh_3> Parameterizer;
    Parameterizer::Error_code err;

	/*typedef typename OpenNL::DefaultLinearSolverTraits<typename ParameterizationMesh_3::NT>
			SolverPolicy;*/
	if (CGAL_CLIB_STD::strcmp(solver,"opennl") == 0)
	{
		typedef typename OpenNL::DefaultLinearSolverTraits<typename ParameterizationMesh_3::NT>
			SolverPolicy;
		typedef typename No_constrian_policy<ParameterizationMesh_3,
										  SolverPolicy::Vector,SolverPolicy::Matrix
	                                    >
										ConstrainPolicy;

		err = parameterize<ParameterizationMesh_3, 
							SolverPolicy, 
							ConstrainPolicy>(mesh, type, border);
	}
	else if (CGAL_CLIB_STD::strcmp(solver,"taucs") == 0)
	{
#ifdef CGAL_USE_TAUCS
		typedef typename CGAL::Taucs_solver_traits<typename ParameterizationMesh_3::NT>
			SolverPolicy;
		typedef typename No_constrian_policy<ParameterizationMesh_3,
										  SolverPolicy::Vector,SolverPolicy::Matrix
	                                    >
										ConstrainPolicy;

		err = parameterize<ParameterizationMesh_3, 
							SolverPolicy, 
							ConstrainPolicy>(mesh, type, border);
#else
		std::cerr << "FATAL ERROR: TAUCS is not installed" << std::endl;
		err = Parameterizer::ERROR_WRONG_PARAMETER;
#endif
	}
	else
	{
		std::cerr << "FATAL ERROR: invalid solver parameter " << solver << std::endl;
		err = Parameterizer::ERROR_WRONG_PARAMETER;
	}

	return err;	
}

// Call appropriate parameterization method based on application parameters
template<class ParameterizationMesh_3,       // 3D surface
         class SolverPolicy  // Traits class to solve a sparse linear system
>                                  
typename CGAL::Parameterizer_traits_3<ParameterizationMesh_3>::Error_code
parameterize(ParameterizationMesh_3& mesh,   // Mesh parameterization adaptor
             const char *type="mvc",      // type of parameterization (see usage)
             const char *border="circle")
{
	typedef typename No_constrian_policy<ParameterizationMesh_3,
										  SolverPolicy::Vector,SolverPolicy::Matrix>
										  ConstrainPolicy;
	return parameterize<ParameterizationMesh_3, 
		                SolverPolicy, 
						ConstrainPolicy>(mesh, type, border);

}

// Call appropriate parameterization method based on application parameters
template<class ParameterizationMesh_3,       // 3D surface
         class SolverPolicy,  // Traits class to solve a sparse linear system
		 class ConstrainPolicy
>                                  
typename CGAL::Parameterizer_traits_3<ParameterizationMesh_3>::Error_code
parameterize(ParameterizationMesh_3& mesh,   // Mesh parameterization adaptor
             const char *type="mvc",      // type of parameterization (see usage)
             const char *border="circle")    // type of border parameterization (see usage)
{    
	typedef typename CGAL::Circular_border_arc_length_parameterizer_3<ParameterizationMesh_3> 
		                                  Circular_border_policy;
	typedef typename CGAL::Square_border_arc_length_parameterizer_3<ParameterizationMesh_3> 
		                                  Square_border_policy;
	
	
	typename CGAL::Parameterizer_traits_3<ParameterizationMesh_3>::Error_code err;

    if ( (CGAL_CLIB_STD::strcmp(type,"mvc") == 0) && (CGAL_CLIB_STD::strcmp(border,"circle") == 0) )
    {
        err = CGAL::parameterize(
								mesh,
								DGAL::Planar_parameterizer_3<
									 ParameterizationMesh_3, 
									 Circular_border_policy,//Square_border_policy,
									 MVC_policy<ParameterizationMesh_3>,
									 SolverPolicy,
									 ConstrainPolicy
									>()
								);
    }
    else if ( (CGAL_CLIB_STD::strcmp(type,"mvc") == 0) && (CGAL_CLIB_STD::strcmp(border,"square") == 0) )
    {
        err = CGAL::parameterize(
								mesh,
								DGAL::Planar_parameterizer_3<
									 ParameterizationMesh_3, 
									 Square_border_policy,
									 MVC_policy<ParameterizationMesh_3>,
									 SolverPolicy,
									 ConstrainPolicy
									>()
								);
    }
    //else if ( (CGAL_CLIB_STD::strcmp(type,"barycentric") == 0) && (CGAL_CLIB_STD::strcmp(border,"circle") == 0) )
    //{
    //    err = CGAL::parameterize(
    //        mesh,
    //        CGAL::Barycentric_mapping_parameterizer_3<
    //            ParameterizationMesh_3,
    //            CGAL::Circular_border_arc_length_parameterizer_3<ParameterizationMesh_3>,
    //            SparseLinearAlgebraTraits_d
    //        >());
    //}
    //else if ( (CGAL_CLIB_STD::strcmp(type,"barycentric") == 0) && (CGAL_CLIB_STD::strcmp(border,"square") == 0) )
    //{
    //    err = CGAL::parameterize(
    //        mesh,
    //        CGAL::Barycentric_mapping_parameterizer_3<
    //            ParameterizationMesh_3,
    //            CGAL::Square_border_arc_length_parameterizer_3<ParameterizationMesh_3>,
    //            SparseLinearAlgebraTraits_d
    //        >());
    //}
    else if ( (CGAL_CLIB_STD::strcmp(type,"dcp") == 0) && (CGAL_CLIB_STD::strcmp(border,"circle") == 0) )
    {
        err = CGAL::parameterize(
								mesh,
								DGAL::Planar_parameterizer_3<
									 ParameterizationMesh_3, 
									 Circular_border_policy,//Square_border_policy,
									 DCP_policy<ParameterizationMesh_3>,//MVC_policy<ParameterizationMesh_3>,
									 SolverPolicy,
									 ConstrainPolicy
									>()
								);
    }
    else if ( (CGAL_CLIB_STD::strcmp(type,"dcp") == 0) && (CGAL_CLIB_STD::strcmp(border,"square") == 0) )
    {
        err = CGAL::parameterize(
								mesh,
								DGAL::Planar_parameterizer_3<
									 ParameterizationMesh_3, 
									 Square_border_policy,
									 DCP_policy<ParameterizationMesh_3>,//MVC_policy<ParameterizationMesh_3>,
									 SolverPolicy,
									 ConstrainPolicy
									>()
								);
    }
 /*   else if ( (CGAL_CLIB_STD::strcmp(type,"authalic") == 0) && (CGAL_CLIB_STD::strcmp(border,"circle") == 0) )
    {
        err = CGAL::parameterize(
            mesh,
            CGAL::Discrete_authalic_parameterizer_3<
                ParameterizationMesh_3,
                CGAL::Circular_border_arc_length_parameterizer_3<ParameterizationMesh_3>,
                SparseLinearAlgebraTraits_d
            >());
    }
    else if ( (CGAL_CLIB_STD::strcmp(type,"authalic") == 0) && (CGAL_CLIB_STD::strcmp(border,"square") == 0) )
    {
        err = CGAL::parameterize(
            mesh,
            CGAL::Discrete_authalic_parameterizer_3<
                ParameterizationMesh_3,
                CGAL::Square_border_arc_length_parameterizer_3<ParameterizationMesh_3>,
                SparseLinearAlgebraTraits_d
            >());
    }
    else if ( (CGAL_CLIB_STD::strcmp(type,"lscm") == 0) && (CGAL_CLIB_STD::strcmp(border,"2pts") == 0) )
    {
        err = CGAL::parameterize(
            mesh,
            CGAL::LSCM_parameterizer_3<
                ParameterizationMesh_3,
                CGAL::Two_vertices_parameterizer_3<ParameterizationMesh_3>,
                SparseLinearAlgebraTraits_d
            >());
    }*/
    else
    {
        std::cerr << "FATAL ERROR: invalid parameters combination " << type << " + " << border << std::endl;
        err = CGAL::Parameterizer_traits_3<ParameterizationMesh_3>::ERROR_WRONG_PARAMETER;
    }

    return err;
}

DGAL_END_NAMESPACE

#endif //DGAL_PARAMETERIZATION_FACTORY_H

