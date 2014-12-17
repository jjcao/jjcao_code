#ifndef FM_DISTANCE_FIELD_BUILDER__H
#define FM_DISTANCE_FIELD_BUILDER__H

#include <DGAL/config.h>
using std::string;
using std::cerr;
using std::cout;
using std::endl;
#include <gw/gw_core/GW_MathsWrapper.h>
#include <gw/gw_geodesic/GW_GeodesicMesh.h>
#include <gw/gw_toolkit/GW_OFFLoader.h>
using namespace GW;

DGAL_BEGIN_NAMESPACE
void* fb_instance;

/// Compute a geodesic distance field from specified vertices. 
template<
	class Polyhedron_3
>
class Fm_distance_field_builder
{
public:
	typedef typename Polyhedron_3 Polyhedron;
	typedef Fm_distance_field_builder<Polyhedron>		Self;
	typedef typename Polyhedron::Vertex Vertex;
	typedef typename Polyhedron::Vertex_handle Vertex_handle;
	typedef typename Polyhedron::Vertex_iterator Vertex_iterator;
	typedef typename Polyhedron::Point_3 Point_3;
	typedef typename Polyhedron::Halfedge_around_vertex_circulator Halfedge_vertex_circulator;
	typedef typename Polyhedron::Facet_iterator Facet_iterator;
	typedef typename Polyhedron::Halfedge_around_facet_circulator Halfedge_facet_circulator;
	typedef typename Polyhedron::Vector_3 Vector_3;
	struct Set_scalar
	{
		Set_scalar(std::vector<double>& s):scalars(s){}
		void operator() (Vertex& vh){
			double s = scalars[vh.index()];
			if (abs(s)<1.0e-015)	
				s = 0;
			vh.s(s);
		}
		void operator()(Vertex_handle& elem){
			this->operator ()(*elem);
		}
		std::vector<double>& scalars;
	};
public:
	Fm_distance_field_builder():m_geo_mesh(0),m_mesh(0){
		m_dist_max = 1e9;
	}
	bool compute(Polyhedron* mesh,std::list<Vertex_handle>& start_vertices, 
				std::vector<double> *init_start_distance=0){
		if (mesh!=m_mesh)
			set_mesh(mesh);
#ifdef OUTPUT_INFO
		CGAL::Timer timer;	timer.start();
#endif	
		std::vector<int> start_vertices_id;
		start_vertices_id.reserve(start_vertices.size());
		for(std::list<Vertex_handle>::iterator it = start_vertices.begin(); it!=start_vertices.end();++it)
		{
			start_vertices_id.push_back( (*it)->index());		
		}
		
		// for using memeber function as callback, I've to set all member variables used in callback to self, but not this
		Self* self = (Self*)fb_instance;
		int nverts = m_geo_mesh->GetNbrVertex();
		self->m_nbr_iter = 0;
		self->m_niter_max = std::numeric_limits<int>::max();
		self->m_niter_max = std::min(self->m_niter_max, int(1.2*nverts));

		if(self->m_weight.empty())
			self->m_weight.insert(self->m_weight.begin(), nverts, 1.0);

		// set up fast marching		
		m_geo_mesh->ResetGeodesicMesh();
		for( int i=0; i<start_vertices_id.size(); ++i )
		{
			GW_GeodesicVertex* v = (GW_GeodesicVertex*) m_geo_mesh->GetVertex((GW_U32) start_vertices_id[i]);
			GW_ASSERT( v!=NULL );
			m_geo_mesh->AddStartVertex( *v );
		}
		m_geo_mesh->SetUpFastMarching();
		m_geo_mesh->RegisterWeightCallbackFunction( WeightCallback );
		m_geo_mesh->RegisterForceStopCallbackFunction( StopMarchingCallback );
		m_geo_mesh->RegisterVertexInsersionCallbackFunction( InsersionCallback );

		if( !self->m_H.empty() )
			m_geo_mesh->RegisterHeuristicToGoalCallbackFunction( HeuristicCallback );
			
		// initialize the distance of the starting points
		
		if( init_start_distance )
		{
			std::vector<double>& isd = * init_start_distance;
			for( int i=0; i<start_vertices_id.size(); ++i )
			{
				GW_GeodesicVertex* v = (GW_GeodesicVertex*) m_geo_mesh->GetVertex((GW_U32) start_vertices_id[i]);
				GW_ASSERT( v!=NULL );
				v->SetDistance( isd[i] );
			}
		}
			
		// perform fast marching
		m_geo_mesh->PerformFastMarching();

		// output result
		m_distance.clear();		m_distance.reserve(nverts);
		m_status.clear();		m_status.reserve(nverts);
		m_Q.clear();		    m_Q.reserve(nverts);
		bool result(true);
		for( int i=0; i<nverts; ++i )
		{
			GW_GeodesicVertex* v = (GW_GeodesicVertex*) m_geo_mesh->GetVertex((GW_U32) i);
			GW_ASSERT( v!=NULL );
			m_distance.push_back( v->GetDistance());
			m_status.push_back( v->GetState());
			if ( v->GetState()<2) //kFar 0,kAlive 1, kDead 2
				result = false;
			GW_GeodesicVertex* v1 = v->GetFront();
			if( v1==NULL )
				m_Q.push_back( -1);
			else
				m_Q.push_back( v1->GetID());
		}
		
		std::for_each(mesh->vertices_begin(),mesh->vertices_end(),Set_scalar(m_distance));

#ifdef OUTPUT_INFO
		std::cout << "fm build: " << timer.time() << " seconds." << std::endl;
		timer.reset();
#endif	
		return result;
	}
	void fm_mesh(char* filename)
	{
		// create the mesh
		GW_GeodesicMesh mesh;
		cout << "Loading OFF file " << filename << " ... ";
		GW_I32 nRet = GW_OFFLoader::Load( mesh, filename );
		if( nRet<0 )
		{
			cout << endl << "Can't load file.";
			return;
		}
		cout << "done." << endl;

		/* set up the connectivity */
		cout << "Building connectivity ... ";
		mesh.BuildConnectivity();
		cout << "done." << endl;
	}
	void* matrix(){return 0;}
	void* factor(){return 0;}
private:
	static GW_Float WeightCallback(GW_GeodesicVertex& Vert)
	{
		Self* self = (Self*)fb_instance;
		GW_U32 i = Vert.GetID();
		return self->m_weight[i];
	}
	static GW_Bool StopMarchingCallback( GW_GeodesicVertex& Vert )
	{
		Self* self = (Self*)fb_instance;
		// check if the end point has been reached
		GW_U32 i = Vert.GetID();
		if( Vert.GetDistance()>self->m_dist_max )
			return true;
			
		for(std::vector<int>::iterator it = self->m_end_vertices_id.begin(); it!=self->m_end_vertices_id.end();++it)
			if( *it==i )
				return true;
		return false;
	}
	static GW_Bool InsersionCallback( GW_GeodesicVertex& Vert, GW_Float rNewDist )
	{
		Self* self = (Self*)fb_instance;
		// check if the distance of the new point is less than the given distance
		GW_U32 i = Vert.GetID();
		bool doinsersion = self->m_nbr_iter <= self->m_niter_max;
		if( !self->m_L.empty() )
			doinsersion = doinsersion && (rNewDist<self->m_L[i]);
		++(self->m_nbr_iter);
		return doinsersion;
	}
	static GW_Float HeuristicCallback( GW_GeodesicVertex& Vert )
	{
		Self* self = (Self*)fb_instance;
		// return the heuristic distance
		GW_U32 i = Vert.GetID();
		return self->m_H[i];
	}
	/// Preconditions:
	/// vertices of the mesh are indexed.
	void set_mesh(Polyhedron* mesh)
	{		
		if(m_geo_mesh) delete m_geo_mesh;
		
		m_geo_mesh = new GW_GeodesicMesh;
		int nverts(mesh->size_of_vertices());
		m_geo_mesh->SetNbrVertex(nverts);		
		for( Vertex_iterator it = mesh->vertices_begin(); it!=mesh->vertices_end(); ++it)
		{
			GW_GeodesicVertex& vert = (GW_GeodesicVertex&) m_geo_mesh->CreateNewVertex();
			Point_3 pt = it->point();
			vert.SetPosition( GW_Vector3D(pt.x(),pt.y(),pt.z()) );
			m_geo_mesh->SetVertex(it->index(), &vert);
		}
		int nfaces(mesh->size_of_facets());
		m_geo_mesh->SetNbrFace(nfaces);
		int i(0);
		for(Facet_iterator fi = mesh->facets_begin(); fi!=mesh->facets_end(); ++fi,++i)
		{
			GW_GeodesicFace& face = (GW_GeodesicFace&) m_geo_mesh->CreateNewFace();
			Halfedge_facet_circulator fc = fi->facet_begin();
			Halfedge_facet_circulator fend = fc;
			GW_Vertex* v[3];
			int j(0);
			CGAL_For_all(fc,fend)
			{
				v[j] = m_geo_mesh->GetVertex(fc->vertex()->index()); GW_ASSERT( v[j]!=NULL );
				++j;
			}
			face.SetVertex(*v[0],*v[1],*v[2]);
			
			m_geo_mesh->SetFace(i, &face);
		}
		m_geo_mesh->BuildConnectivity();
		
		m_mesh = mesh;
	}
	void set_end_vertices(std::vector<int>& in){
		Self* self = (Self*)fb_instance;
		self->m_end_vertices_id = in;
	}
	Polyhedron* m_mesh;//just a reference, do not new and del in this class
	GW_GeodesicMesh *m_geo_mesh;
	int m_nbr_iter;//current iteration step
	int m_niter_max;//max iteration steps. stop when a given number of iterations is reached.
	double m_dist_max;// max distance
	std::vector<double> m_weight;	// weight,provide non-uniform speed for each vertex.
	std::vector<int> m_end_vertices_id; //stop when these points are reached
	std::vector<double> m_H;	// heuristic(typically that try to guess the distance 
	//that remains from a given node to a given target). This is an array of same size as m_weight.
	std::vector<double> m_L;	// reduce the set of explored points. 
	// Only points with current distance smaller than L will be expanded. 
	//Set some entries of L to -Inf to avoid any exploration of these points.
	
	//////////////////////////////////////////////////////////////////////////
	/// results
	std::vector<double> m_distance;//the distance function to the set of starting points.
	std::vector<int> m_status;//the final state of the points : -1 for dead (ie the distance
	//	%       has been computed), 0 for open (ie the distance is only a temporary
	//	%       value), 1 for far (ie point not already computed). Distance function
	//	%       for far points is Inf.
	std::vector<int> m_Q;//index of the closest point. Q is set to 0 for far points. Q provide a Voronoi decomposition of the domain. 
};

DGAL_END_NAMESPACE
#endif//FM_DISTANCE_FIELD_BUILDER__H