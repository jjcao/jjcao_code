#ifndef DGAL_BORDER_PARAMETERIZATION_POLICIES_H
#define DGAL_BORDER_PARAMETERIZATION_POLICIES_H

#include <CGAL/circular_border_parameterizer_3.h>
#include <CGAL/square_border_parameterizer_3.h>
#include <DGAL/config.h>
#include <fstream>

DGAL_BEGIN_NAMESPACE
template<class ParameterizationMesh_3>           //< 3D surface
class Square_border_fixed_start_arc_length_parameterizer_3
	: public CGAL::Square_border_arc_length_parameterizer_3<ParameterizationMesh_3>
{
public:
	typedef ParameterizationMesh_3          Adaptor;
private:
	// Private types
private:
    // Mesh_Adaptor_3 subtypes:
    typedef typename Adaptor::NT            NT;
    typedef typename Adaptor::Point_2       Point_2;
    typedef typename Adaptor::Point_3       Point_3;
    typedef typename Adaptor::Vector_2      Vector_2;
    typedef typename Adaptor::Vector_3      Vector_3;
    typedef typename Adaptor::Facet         Facet;
    typedef typename Adaptor::Facet_handle  Facet_handle;
    typedef typename Adaptor::Facet_const_handle
                                            Facet_const_handle;
    typedef typename Adaptor::Facet_iterator Facet_iterator;
    typedef typename Adaptor::Facet_const_iterator
                                            Facet_const_iterator;
    typedef typename Adaptor::Vertex        Vertex;
    typedef typename Adaptor::Vertex_handle Vertex_handle;
    typedef typename Adaptor::Vertex_const_handle
                                            Vertex_const_handle;
    typedef typename Adaptor::Vertex_iterator Vertex_iterator;
    typedef typename Adaptor::Vertex_const_iterator
                                            Vertex_const_iterator;
    typedef typename Adaptor::Border_vertex_iterator
                                            Border_vertex_iterator;
    typedef typename Adaptor::Border_vertex_const_iterator
                                            Border_vertex_const_iterator;
    typedef typename Adaptor::Vertex_around_facet_circulator
                                            Vertex_around_facet_circulator;
    typedef typename Adaptor::Vertex_around_facet_const_circulator
                                            Vertex_around_facet_const_circulator;
    typedef typename Adaptor::Vertex_around_vertex_circulator
                                            Vertex_around_vertex_circulator;
    typedef typename Adaptor::Vertex_around_vertex_const_circulator
                                            Vertex_around_vertex_const_circulator;

	typedef typename std::vector<double>    Offset_map;
	typedef typename std::list<Vertex_handle>::iterator BBorder_vertex_iterator;
private:
	std::list<Vertex_handle> border_vertices_;
	BBorder_vertex_iterator border_vertices_begin()
	{
		return border_vertices_.begin();
	}
	BBorder_vertex_iterator border_vertices_end(){
		return border_vertices_.end();}
public:
	/// Assign to mesh's border vertices a 2D position (ie a (u,v) pair)
    /// on border's shape. Mark them as "parameterized".
	typename CGAL::Parameterizer_traits_3<Adaptor>::Error_code 
		           parameterize_border(Adaptor& mesh)
	{
	#ifdef DEBUG_TRACE
		std::cerr << "  map on a square" << std::endl;
	#endif

		// Nothing to do if no border
		if (mesh.mesh_main_border_vertices_begin() == mesh.mesh_main_border_vertices_end())
			return CGAL::Parameterizer_traits_3<Adaptor>::ERROR_BORDER_TOO_SHORT;

		// Compute the total border length
		double total_len = compute_border_length(mesh);
		if (total_len == 0)
			return CGAL::Parameterizer_traits_3<Adaptor>::ERROR_BORDER_TOO_SHORT;

		// map to [0,4[
		double len = 0.0;           // current position on square in [0, total_len[
		Offset_map offset;          // vertex index -> offset map
		offset.resize(mesh.count_mesh_vertices());

		// reform a border with pointed start vertex
	order_border(mesh.mesh_main_border_vertices_begin(), mesh.mesh_main_border_vertices_end());
	
	// parameterize boundary
	//std::ofstream ofile("boundary.txt");
	//ofile << "{";
	BBorder_vertex_iterator it = border_vertices_begin();	BBorder_vertex_iterator end = border_vertices_end();	
    for(; it != end;++it)
    {
		CGAL_surface_mesh_parameterization_assertion(mesh.is_vertex_on_main_border(*it));

		offset[mesh.get_vertex_index(*it)] = 4.0f*len/total_len;
								// current position on square in [0,4[

		// Get next iterator (looping)
		BBorder_vertex_iterator next = it;
		next++;
		if(next == border_vertices_end())
			next = border_vertices_begin();

		// Add edge "length" to 'len'
		len += compute_edge_length(mesh, *it, *next);
	}

		// First square corner is mapped to first vertex.
		// Then find closest points for three other corners.
		BBorder_vertex_iterator it0 = border_vertices_begin();
		BBorder_vertex_iterator it1 = closest_iterator(mesh, offset, 1.0);
		BBorder_vertex_iterator it2 = closest_iterator(mesh, offset, 2.0);
		BBorder_vertex_iterator it3 = closest_iterator(mesh, offset, 3.0);
		//
		// We may get into trouble if the border is too short
		if (it0 == it1 || it1 == it2 || it2 == it3 || it3 == it0)
			return CGAL::Parameterizer_traits_3<Adaptor>::ERROR_BORDER_TOO_SHORT;
		//
		// Snap these vertices to corners
		offset[mesh.get_vertex_index(*it0)] = 0.0;
		offset[mesh.get_vertex_index(*it1)] = 1.0;
		offset[mesh.get_vertex_index(*it2)] = 2.0;
		offset[mesh.get_vertex_index(*it3)] = 3.0;

		// Set vertices along square's sides and mark them as "parameterized"
		for(it = it0; it != it1; it++) // 1st side
		{
			Point_2 uv(offset[mesh.get_vertex_index(*it)], 0.0);
			mesh.set_vertex_uv(*it, uv);
			mesh.set_vertex_parameterized(*it, true);
		}
		for(it = it1; it != it2; it++) // 2nd side
		{
			Point_2 uv(1.0, offset[mesh.get_vertex_index(*it)]-1);
			mesh.set_vertex_uv(*it, uv);
			mesh.set_vertex_parameterized(*it, true);
		}
		for(it = it2; it != it3; it++) // 3rd side
		{
			Point_2 uv(3-offset[mesh.get_vertex_index(*it)], 1.0);
			mesh.set_vertex_uv(*it, uv);
			mesh.set_vertex_parameterized(*it, true);
		}
		for(it = it3; it != border_vertices_end(); it++) // 4th side
		{
			Point_2 uv(0.0, 4-offset[mesh.get_vertex_index(*it)]);
			mesh.set_vertex_uv(*it, uv);
			mesh.set_vertex_parameterized(*it, true);
		}

		return CGAL::Parameterizer_traits_3<Adaptor>::OK;
	}
// Private operations
private:
	// reform a border with pointed start vertex
	void order_border(Border_vertex_iterator begin, Border_vertex_iterator end)
	{
		border_vertices_.clear();

		std::list<Vertex_handle> temp1;
		bool b_2_first(false);
		bool bFound(false);
		for(Border_vertex_iterator it = begin; it != end;++it)
		{
			Vertex_handle vch = it;
			
			if ( !bFound && (SECOND_BORDER_VERTEX == vch->halfedge()->tag())  ){// find second first
				b_2_first = true;
			}

			//original:-2; first:FIRST_BORDER_VERTEX; second: SECOND_BORDER_VERTEX
			if ( FIRST_BORDER_VERTEX == vch->halfedge()->tag())// find first
			{
				bFound = true;
			}
			if ( bFound)
			{
				border_vertices_.push_back(vch);
			}
			else
			{
				temp1.push_back(vch);
			}
		}

		border_vertices_.splice(border_vertices_.end(), temp1);

		if ( b_2_first)
		{
			temp1.clear();
			temp1.push_back(border_vertices_.front());

			BBorder_vertex_iterator it = --border_vertices_end();
			for ( ; it!= border_vertices_begin(); --it)
			{
				temp1.push_back(*it);
			}
			border_vertices_.clear();
			border_vertices_.splice(border_vertices_.end(), temp1);
		}	

		//
		/*std::cout << "param boundary order:" << std::endl;
		std::ofstream ofile("boundary.txt");
		ofile << "{";
		for ( Border_vertex_const_iterator it = border_vertices_begin(); it!= border_vertices_end(); ++it)
		{
			Point_3 p3 = mesh.get_vertex_position( it);
			ofile << "{" << p3.x() << ", " << p3.y() << "},";
		}
		ofile << "}";
		ofile.close();*/
	}
    /// Compute the total length of the border.
    double compute_border_length(const Adaptor& mesh)
	{
		double len = 0.0;
		for(Border_vertex_const_iterator it = mesh.mesh_main_border_vertices_begin();
			it != mesh.mesh_main_border_vertices_end();
			it++)
		{
			CGAL_surface_mesh_parameterization_assertion(mesh.is_vertex_on_main_border(it));

			// Get next iterator (looping)
			Border_vertex_const_iterator next = it;
			next++;
			if(next == mesh.mesh_main_border_vertices_end())
				next = mesh.mesh_main_border_vertices_begin();

			// Add 'length' of it -> next vector to 'len'
			len += compute_edge_length(mesh, it, next);
		}
		return len;
	}

    /// Get mesh iterator whose offset is closest to 'value'.
    BBorder_vertex_iterator closest_iterator(Adaptor& mesh, const Offset_map& offset, double value)
	{
		BBorder_vertex_iterator best;
		double min = DBL_MAX;           // distance for 'best'

		for (BBorder_vertex_iterator it = border_vertices_begin();
			 it != border_vertices_end();
			 it++)
		{
			double d = std::fabs(offset[mesh.get_vertex_index(*it)] - value);
			if (d < min)
			{
				best = it;
				min = d;
			}
		}

		return best;
	}
};


template<class ParameterizationMesh_3>           //< 3D surface
class Circular_border_fixed_start_arc_length_parameterizer_3
	: public CGAL::Circular_border_arc_length_parameterizer_3<ParameterizationMesh_3>
{
public:
	typedef ParameterizationMesh_3          Adaptor;
private:
	// Private types
private:
    // Mesh_Adaptor_3 subtypes:
    typedef typename Adaptor::NT            NT;
    typedef typename Adaptor::Point_2       Point_2;
    typedef typename Adaptor::Point_3       Point_3;
    typedef typename Adaptor::Vector_2      Vector_2;
    typedef typename Adaptor::Vector_3      Vector_3;
    typedef typename Adaptor::Facet         Facet;
    typedef typename Adaptor::Facet_handle  Facet_handle;
    typedef typename Adaptor::Facet_const_handle
                                            Facet_const_handle;
    typedef typename Adaptor::Facet_iterator Facet_iterator;
    typedef typename Adaptor::Facet_const_iterator
                                            Facet_const_iterator;
    typedef typename Adaptor::Vertex        Vertex;
    typedef typename Adaptor::Vertex_handle Vertex_handle;
    typedef typename Adaptor::Vertex_const_handle
                                            Vertex_const_handle;
    typedef typename Adaptor::Vertex_iterator Vertex_iterator;
    typedef typename Adaptor::Vertex_const_iterator
                                            Vertex_const_iterator;
    typedef typename Adaptor::Border_vertex_iterator
                                            Border_vertex_iterator;
    typedef typename Adaptor::Border_vertex_const_iterator
                                            Border_vertex_const_iterator;
    typedef typename Adaptor::Vertex_around_facet_circulator
                                            Vertex_around_facet_circulator;
    typedef typename Adaptor::Vertex_around_facet_const_circulator
                                            Vertex_around_facet_const_circulator;
    typedef typename Adaptor::Vertex_around_vertex_circulator
                                            Vertex_around_vertex_circulator;
    typedef typename Adaptor::Vertex_around_vertex_const_circulator
                                            Vertex_around_vertex_const_circulator;
public:
	/// Assign to mesh's border vertices a 2D position (ie a (u,v) pair)
    /// on border's shape. Mark them as "parameterized".
	typename CGAL::Parameterizer_traits_3<Adaptor>::Error_code
                                        parameterize_border(Adaptor& mesh);
private:
	std::list<Vertex_handle> border_vertices_;
	Border_vertex_iterator border_vertices_begin(){
		return Border_vertex_iterator(border_vertices_.begin());}
	Border_vertex_iterator border_vertices_end(){
		return Border_vertex_iterator(border_vertices_.end());}
	// reform a border with pointed start vertex
	void order_border(Border_vertex_iterator begin, Border_vertex_iterator end)
	{
		border_vertices_.clear();
		std::list<Vertex_handle> temp1;
		bool b_2_first(false);
		bool bFound(false);
		for(Border_vertex_iterator it = begin; it != end;++it)
		{
			Vertex_handle vch = it;
			
			if ( !bFound && (SECOND_BORDER_VERTEX == vch->halfedge()->tag())  ){// find second first
				b_2_first = true;
			}

			//original:-2; first:FIRST_BORDER_VERTEX; second: SECOND_BORDER_VERTEX
			if ( FIRST_BORDER_VERTEX == vch->halfedge()->tag())// find first
			{
				bFound = true;
			}
			if ( bFound)
			{
				border_vertices_.push_back(vch);
			}
			else
			{
				temp1.push_back(vch);
			}
		}

		border_vertices_.splice(border_vertices_.end(), temp1);
		if ( b_2_first)
		{
			temp1.clear();
			temp1.push_back(border_vertices_.front());

			Border_vertex_iterator it = --border_vertices_end();
			for ( ; it!= border_vertices_begin(); --it)
			{
				temp1.push_back(it);
			}
			border_vertices_.clear();
			border_vertices_.splice(border_vertices_.end(), temp1);
		}	
	}
    /// Compute the total length of the border
    double compute_border_length(const Adaptor& mesh)
	{
		double len = 0.0;
		for(Border_vertex_const_iterator it = mesh.mesh_main_border_vertices_begin();
			it != mesh.mesh_main_border_vertices_end();
			it++)
		{
			CGAL_surface_mesh_parameterization_assertion(mesh.is_vertex_on_main_border(it));

			// Get next iterator (looping)
			Border_vertex_const_iterator next = it;
			next++;
			if(next == mesh.mesh_main_border_vertices_end())
				next = mesh.mesh_main_border_vertices_begin();

			// Add 'length' of it -> next vector to 'len'
			len += compute_edge_length(mesh, it, next);
		}
		return len;
	}
};
/// Assign to mesh's border vertices a 2D position (ie a (u,v) pair)
/// on border's shape. Mark them as "parameterized".
template<class Adaptor>
inline 
typename CGAL::Parameterizer_traits_3<Adaptor>::Error_code
Circular_border_fixed_start_arc_length_parameterizer_3<Adaptor>::
parameterize_border(Adaptor& mesh)
{
#ifdef DEBUG_TRACE
    std::cerr << "  map on a circle" << std::endl;
#endif

    // Nothing to do if no border
    if (mesh.mesh_main_border_vertices_begin() == mesh.mesh_main_border_vertices_end())
		return CGAL::Parameterizer_traits_3<Adaptor>::ERROR_BORDER_TOO_SHORT;

    // Compute the total border length
    double total_len = compute_border_length(mesh);
    if (total_len == 0)
        return CGAL::Parameterizer_traits_3<Adaptor>::ERROR_BORDER_TOO_SHORT;

    double pi = 3.14159265359;
    double tmp = 2*pi/total_len;
    double len = 0.0;           // current position on circle in [0, total_len]

	// reform a border with pointed start vertex
	order_border(mesh.mesh_main_border_vertices_begin(), mesh.mesh_main_border_vertices_end());
	
	// parameterize boundary
#ifdef _DEBUG
	std::ofstream ofile("output\\boundary.txt");
	ofile << "{";
#endif
    for(Border_vertex_iterator it = border_vertices_begin(); it != border_vertices_end();++it)
    {
        CGAL_surface_mesh_parameterization_assertion(mesh.is_vertex_on_main_border(it));

        double angle = len*tmp; // current position on the circle in radians

        // map vertex on unit circle
        Point_2 uv;
        uv = Point_2(0.5+0.5*std::cos(-angle),0.5+0.5*std::sin(-angle));
        mesh.set_vertex_uv(it, uv);
#ifdef _DEBUG
		ofile << "{" << uv.x() << ", " << uv.y() << "},";
#endif
        // Mark vertex as "parameterized"
        mesh.set_vertex_parameterized(it, true);

        // Get next iterator (looping)
        Border_vertex_iterator next = it;
        next++;
        if(next == border_vertices_end())
            next = border_vertices_begin();

        // Add 'length' of it -> next vector to 'len'
        len += compute_edge_length(mesh, it, next);
    }
#ifdef _DEBUG
	ofile << "}";
	ofile.close();
#endif

    return CGAL::Parameterizer_traits_3<Adaptor>::OK;
}

//mapping border vertices onto a unit circle centered at (0.5,0.5) by arc length, 
//then rescale them one by one according to their distances, i.e. s(). 
template<class ParameterizationMesh_3>           //< 3D surface
class Measured_border_arc_length_policy
	: public CGAL::Circular_border_arc_length_parameterizer_3<ParameterizationMesh_3>
{
public:
	typedef ParameterizationMesh_3          Adaptor;
private:
	// Private types
private:
    // Mesh_Adaptor_3 subtypes:
    typedef typename Adaptor::NT            NT;
    typedef typename Adaptor::Point_2       Point_2;
    typedef typename Adaptor::Point_3       Point_3;
    typedef typename Adaptor::Vector_2      Vector_2;
    typedef typename Adaptor::Vector_3      Vector_3;
    typedef typename Adaptor::Facet         Facet;
    typedef typename Adaptor::Facet_handle  Facet_handle;
    typedef typename Adaptor::Facet_const_handle
                                            Facet_const_handle;
    typedef typename Adaptor::Facet_iterator Facet_iterator;
    typedef typename Adaptor::Facet_const_iterator
                                            Facet_const_iterator;
    typedef typename Adaptor::Vertex        Vertex;
    typedef typename Adaptor::Vertex_handle Vertex_handle;
    typedef typename Adaptor::Vertex_const_handle
                                            Vertex_const_handle;
    typedef typename Adaptor::Vertex_iterator Vertex_iterator;
    typedef typename Adaptor::Vertex_const_iterator
                                            Vertex_const_iterator;
    typedef typename Adaptor::Border_vertex_iterator
                                            Border_vertex_iterator;
    typedef typename Adaptor::Border_vertex_const_iterator
                                            Border_vertex_const_iterator;
    typedef typename Adaptor::Vertex_around_facet_circulator
                                            Vertex_around_facet_circulator;
    typedef typename Adaptor::Vertex_around_facet_const_circulator
                                            Vertex_around_facet_const_circulator;
    typedef typename Adaptor::Vertex_around_vertex_circulator
                                            Vertex_around_vertex_circulator;
    typedef typename Adaptor::Vertex_around_vertex_const_circulator
                                            Vertex_around_vertex_const_circulator;
public:
	/// Assign to mesh's border vertices a 2D position (ie a (u,v) pair)
    /// on border's shape. Mark them as "parameterized".
	typename CGAL::Parameterizer_traits_3<Adaptor>::Error_code
                                        parameterize_border(Adaptor& mesh)
	{
		#ifdef DEBUG_TRACE
			std::cerr << "  map on a circle, then relocate point by a distance" << std::endl;
		#endif

			// Nothing to do if no border
			if (mesh.mesh_main_border_vertices_begin() == mesh.mesh_main_border_vertices_end())
				return CGAL::Parameterizer_traits_3<Adaptor>::ERROR_BORDER_TOO_SHORT;

			// Compute the total border length
			double total_len = compute_border_length(mesh);
			if (total_len == 0)
				return CGAL::Parameterizer_traits_3<Adaptor>::ERROR_BORDER_TOO_SHORT;

			const double PI = 3.14159265359;
			const double tmp = 2*PI/total_len;
			double len = 0.0;           // current position on circle in [0, total_len]
			Border_vertex_iterator it = mesh.mesh_main_border_vertices_begin();
			double radius = it->s();
			for(;it != mesh.mesh_main_border_vertices_end();++it)
			{
				CGAL_surface_mesh_parameterization_assertion(mesh.is_vertex_on_main_border(it));

				double angle = len*tmp; // current position on the "circle" in radians

				// map vertex on unit circle
				Vertex_handle vh = it;
				double r = vh->s()/(radius*2);
				double u = 0.5 + std::cos(angle) * r;
				double v = 0.5 + std::sin(angle) * r;
				mesh.set_vertex_uv(vh, Point_2(u,v)); 

				// Mark vertex as "parameterized"
				mesh.set_vertex_parameterized(it, true);

				// Get next iterator (looping)
				Border_vertex_iterator next = it;
				next++;
				if(next == mesh.mesh_main_border_vertices_end())
					next = mesh.mesh_main_border_vertices_begin();

				// Add 'length' of it -> next vector to 'len'
				len += compute_edge_length(mesh, it, next);
			}

			return CGAL::Parameterizer_traits_3<Adaptor>::OK;
	}
private:
    /// Compute the total length of the border
    double compute_border_length(const Adaptor& mesh)
	{
		double len = 0.0;
		for(Border_vertex_const_iterator it = mesh.mesh_main_border_vertices_begin();
			it != mesh.mesh_main_border_vertices_end();
			it++)
		{
			CGAL_surface_mesh_parameterization_assertion(mesh.is_vertex_on_main_border(it));

			// Get next iterator (looping)
			Border_vertex_const_iterator next = it;
			next++;
			if(next == mesh.mesh_main_border_vertices_end())
				next = mesh.mesh_main_border_vertices_begin();

			// Add 'length' of it -> next vector to 'len'
			len += compute_edge_length(mesh, it, next);
		}
		return len;
	}
};

DGAL_END_NAMESPACE
#endif //DGAL_BORDER_PARAMETERIZATION_POLICIES_H