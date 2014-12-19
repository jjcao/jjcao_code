#include <boost/config.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <vector>
#include "graphAlgorithms.h"

using namespace boost;

int dijkstra(Polyhedron& mesh, int idxSource, double* scalarRange)
{
	int nv = mesh.size_of_vertices();
	std::vector<float> dist(nv);	
	int ne = mesh.size_of_halfedges()/2;

	//////////////////////////////////////////////////////////////
	// from mesh => edge list and weight list
	typedef std::pair<int, int> Edge;
	std::vector<Edge> edges;    edges.reserve(ne);
	std::vector<float> weights; weights.reserve(ne);	
	int i1, i2; float len; 
	Kernel::Point_3 p1, p2; 
	Kernel::Vector_3 edge;
	for (Polyhedron::Edge_iterator edgeIt = mesh.edges_begin(); edgeIt != mesh.edges_end(); ++edgeIt)
	{
		i1 = edgeIt->vertex()->index_;
		i2 = edgeIt->opposite()->vertex()->index_;
		edges.push_back( Edge(i1,i2) );

		p1 = edgeIt->vertex()->point();
		p2 = edgeIt->opposite()->vertex()->point();
		edge = p1 - p2;
		len = std::sqrt(edge*edge);
		weights.push_back(len);
	}

	//////////////////////////////////////////////////////////////
	// edge list & weight list => Graph => dijkstra_shortest_paths
	typedef adjacency_list < listS, vecS, undirectedS, no_property, property < edge_weight_t, float > > Graph;
	typedef graph_traits < Graph >::vertex_descriptor vertex_descriptor;
	typedef graph_traits < Graph >::edge_descriptor edge_descriptor;
	
	Graph g(edges.begin(), edges.end(), weights.begin(), nv);
	property_map<Graph, edge_weight_t>::type weightmap = get(edge_weight, g);

	std::vector<vertex_descriptor> path(nv);
	vertex_descriptor s = vertex(idxSource, g);

	dijkstra_shortest_paths( g, s, predecessor_map(&path[0]).distance_map(&dist[0]) );

	//////////////////////////////////////////////////////////////
	// shortest distance => mesh.vertex.scalar
	int i(0);	
	scalarRange[0] = scalarRange[1] = dist[0];
	for(Polyhedron::Vertex_iterator pv = mesh.vertices_begin(); pv !=  mesh.vertices_end(); ++pv, ++i)
	{
		pv->scalar_ = dist[i];
		if ( pv->scalar_ < scalarRange[0])
			scalarRange[0] = pv->scalar_;
		if ( pv->scalar_ > scalarRange[1])
			scalarRange[1] = pv->scalar_;

		pv->tag_ = path[i];// parent of current vertex
	}

	return 0;
}