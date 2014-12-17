#ifndef _CALCULATOR_H
#define _CALCULATOR_H

#pragma warning(disable:4819 4244 4267 4503 4996)
#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Timer.h>

#include <fstream>
#include <list>
#include <vector>

#include <DGAL/Taucs_solver_traits.h>
#include <OpenNL/linear_solver.h>

//#define USE_OPENNL
#ifdef USE_OPENNL
typedef OpenNL::DefaultLinearSolverTraits<double> Solver;
#else
typedef DGAL::Taucs_solver_traits<double>     Solver;
#endif
typedef DGAL::Taucs_symmetric_solver_traits<double>     SymmetricSolver;

typedef Solver::Vector      Vector;
typedef Solver::Matrix      Matrix;
typedef SymmetricSolver::Matrix      SymmetricMatrix;

template <class Refs, class T, class P, class Norm>
class MyVertex : public CGAL::HalfedgeDS_vertex_base<Refs, T, P>
{
private:
	bool is_constrained_;
	bool is_parameterized_;
	double u_;   // parameter
	double v_;	
	double s_;
	int index_;
	// tag
	int tag_;
	// normal
	Norm normal_;
public:
	// life cycle
	// -1 means uninitialized
	MyVertex():
	  is_constrained_(false),is_parameterized_(false),
		  u_(-1.0), v_(-1.0),s_(-1.0),
		  index_(-1),tag_(-1)
	  {
	  }
	  // repeat mandatory constructors
	  MyVertex(const P& pt):
	  CGAL::HalfedgeDS_vertex_base<Refs, T, P>(pt),
		  is_constrained_(false),is_parameterized_(false),
		  u_(-1.0), v_(-1.0),s_(-1.0),
		  index_(-1),tag_(-1)
	  {		
	  }

	  bool is_parameterized() const {	return is_parameterized_;}
	  void is_parameterized(bool in){is_parameterized_ = in;}

	  // normal
	  typedef Norm Normal_3;
	  Normal_3& normal() { return normal_; }
	  const Normal_3& normal() const { return normal_; }

	  // u,v texture
	  const double u() const {return u_;}
	  const double v() const {  return v_; }
	  const double s() const {  return s_; }
	  void uv(double u, double v)  { u_=u; v_=v; }
	  void u(double in){u_=in;}
	  void v(double in){v_=in;}
	  void s(double in){s_=in;}

	  bool is_constrained() const {
		  return is_constrained_; 
	  }
	  void is_constrained(bool in) { 
		  is_constrained_ = in; 
	  }

	  int index() const {return index_;}
	  void index(int in){ index_ = in;}
	  int tag() const {return tag_;}
	  void tag(int in){ tag_ = in;}
};
struct My_items : public CGAL::Polyhedron_items_3 {
	// wrap vertex
	template <class Refs, class Traits>
	struct Vertex_wrapper
	{
		typedef typename Traits::Point_3  Point;
		typedef typename Traits::Vector_3 Normal;
		typedef MyVertex<Refs, CGAL::Tag_true,Point, Normal> Vertex;
	};
};


typedef CGAL::Cartesian<double>             Kernel;
typedef CGAL::Polyhedron_3<Kernel,My_items> Polyhedron;
typedef Polyhedron::Vertex_handle Vertex_handle;
typedef Polyhedron::Vertex_iterator Vertex_iterator;
typedef Polyhedron::Halfedge_around_vertex_circulator Halfedge_vertex_circulator;
typedef Kernel::Point_2                 Point_2;
typedef Kernel::Point_3					Point_3;
typedef Kernel::Vector_2				Vector_2;
typedef Kernel::Vector_3				Vector_3;
class ConstrianedVertex
{
public:
	ConstrianedVertex(Vertex_handle vh, double s):m_scalar(s),m_vh(vh){}
	double    m_scalar;
	Vertex_handle m_vh;
};
// Compute w_ij = (i,j) coefficient of matrix A for j neighbor vertex of i.
class WeightOption
{
public:
	virtual double operator() (Polyhedron* mesh,
                        Vertex_handle main_vertex_v_i,
                        Halfedge_vertex_circulator neighbor_vertex_v_j) = 0;
protected:
    //                                                    -> ->
    /// Return cotangent of (P,Q,R) corner (i.e. cotan of QP,QR angle).
    double cotangent(const Point_3& P,
                     const Point_3& Q,
                     const Point_3& R)
    {
        Vector_3 u = P - Q;
        Vector_3 v = R - Q;
        // (u . v)/((u x v).len)
        double dot = (u*v);
        Vector_3 cross_vector = CGAL::cross_product(u,v);
        double cross_norm = std::sqrt(cross_vector*cross_vector);
        if(cross_norm != 0.0)
            return (dot/cross_norm);
        else
            return 0.0; // undefined
    }

    //                                                    -> ->
    /// Return tangent of (P,Q,R) corner (i.e. tangent of QP,QR angle).
    double tangent(const Point_3& P,
                   const Point_3& Q,
                   const Point_3& R)
    {
        Vector_3 u = P - Q;
        Vector_3 v = R - Q;
        // (u . v)/((u x v).len)
        double dot = (u*v);
        Vector_3 cross_vector = CGAL::cross_product(u,v);
        double cross_norm = std::sqrt(cross_vector*cross_vector);
        if(dot != 0.0)
            return (cross_norm/dot);
        else
            return 0.0; // undefined
    }

    //                                                       -> ->
    /// Return angle (in radians) of of (P,Q,R) corner (i.e. QP,QR angle).
    double compute_angle_rad(const Point_3& P,
                                    const Point_3& Q,
                                    const Point_3& R)
    {
        static const double PI = 3.14159265359;

        Vector_3 u = P - Q;
        Vector_3 v = R - Q;

        // check
        double product = std::sqrt(u*u) * std::sqrt(v*v);
        if(product == 0)
            return 0.0;

        // cosine
        double dot = (u*v);
        double cosine = dot / product;

        // sine
        Vector_3 w = CGAL::cross_product(u,v);
        double AbsSine = std::sqrt(w*w) / product;

        if(cosine >= 0)
            return std::asin(fix_sine(AbsSine));
        else
            return PI-std::asin(fix_sine(AbsSine));
    }

// Private operations
private:
    /// Fix sine.
    static double fix_sine(double sine)
    {
        if(sine >= 1)
            return 1;
        else if(sine <= -1)
            return -1;
        else
            return sine;
    }
};
class DcpOption : public WeightOption
{
public:
	double operator() (Polyhedron* mesh,
                        Vertex_handle main_vertex_v_i,
                        Halfedge_vertex_circulator neighbor_vertex_v_j)
	{
		Point_3 position_v_i = main_vertex_v_i->point();
		Point_3 position_v_j = neighbor_vertex_v_j->opposite()->vertex()->point();

		// Compute cotangent of (v_i,v_k,v_j) corner (i.e. cotan of v_k corner)
		// if v_k is the vertex before v_j when circulating around v_i
		Halfedge_vertex_circulator previous_vertex_v_k = neighbor_vertex_v_j;
		previous_vertex_v_k --;
		Point_3 position_v_k = previous_vertex_v_k->opposite()->vertex()->point();
		double cotg_beta_ij  = cotangent(position_v_i, position_v_k, position_v_j);

		// Compute cotangent of (v_j,v_l,v_i) corner (i.e. cotan of v_l corner)
		// if v_l is the vertex after v_j when circulating around v_i
		Halfedge_vertex_circulator next_vertex_v_l = neighbor_vertex_v_j;
		next_vertex_v_l ++;
		Point_3 position_v_l = next_vertex_v_l->opposite()->vertex()->point();
		double cotg_alpha_ij = cotangent(position_v_j, position_v_l, position_v_i);

		double weight = cotg_beta_ij+cotg_alpha_ij;
		return weight;
	}
};
class MvcOption : public WeightOption
{
public:
	double operator() (Polyhedron* mesh,
                        Vertex_handle main_vertex_v_i,
                        Halfedge_vertex_circulator neighbor_vertex_v_j)
	{
		Point_3 position_v_i = main_vertex_v_i->point();
		Point_3 position_v_j = neighbor_vertex_v_j->opposite()->vertex()->point();

        // Compute the norm of v_j -> v_i vector
        Vector_3 edge = position_v_i - position_v_j;
        double len = std::sqrt(edge*edge);

        // Compute angle of (v_j,v_i,v_k) corner (ie angle of v_i corner)
        // if v_k is the vertex before v_j when circulating around v_i
        Halfedge_vertex_circulator previous_vertex_v_k = neighbor_vertex_v_j;
        previous_vertex_v_k --;
        Point_3 position_v_k = previous_vertex_v_k->opposite()->vertex()->point();
        double gamma_ij  = compute_angle_rad(position_v_j, position_v_i, position_v_k);

        // Compute angle of (v_l,v_i,v_j) corner (ie angle of v_i corner)
        // if v_l is the vertex after v_j when circulating around v_i
        Halfedge_vertex_circulator next_vertex_v_l = neighbor_vertex_v_j;
        next_vertex_v_l ++;
        Point_3 position_v_l = next_vertex_v_l->opposite()->vertex()->point();
        double delta_ij = compute_angle_rad(position_v_l, position_v_i, position_v_j);

        double weight = 0.0;
        assert(len != 0.0);    // two points are identical!
        if(len != 0.0)
            weight = (std::tan(0.5*gamma_ij) + std::tan(0.5*delta_ij)) / len;
		assert(weight > 0);

        return weight;
	}
};
class TutteOption : public WeightOption
{
public:
	double operator() (Polyhedron* mesh,
                        Vertex_handle main_vertex_v_i,
                        Halfedge_vertex_circulator neighbor_vertex_v_j)
	{        
  //      Halfedge_vertex_circulator vvc = neighbor_vertex_v_j;
  //      int i (0);

		//do {//访问一环顶点
		//	++i;
		//} while ( ++vvc != neighbor_vertex_v_j);

  //      return 1.0/i;
		return 1;
	}
};
class SpringOption1 : public WeightOption
{
public:
	double operator() (Polyhedron* mesh,
                        Vertex_handle main_vertex_v_i,
                        Halfedge_vertex_circulator neighbor_vertex_v_j)
	{        
        Point_3 position_v_i = main_vertex_v_i->point();
		Point_3 position_v_j = neighbor_vertex_v_j->opposite()->vertex()->point();

        // Compute the norm of v_j -> v_i vector
        Vector_3 edge = position_v_i - position_v_j;
        double len = std::sqrt(edge*edge);//edge*edge;//

        double weight = 0.0;
        if(len != 0.0)
            weight = 1.0/len;

        return weight;
	}
};
class SpringOption2 : public WeightOption
{
public:
	double operator() (Polyhedron* mesh,
                        Vertex_handle main_vertex_v_i,
                        Halfedge_vertex_circulator neighbor_vertex_v_j)
	{        
        Point_3 position_v_i = main_vertex_v_i->point();
		Point_3 position_v_j = neighbor_vertex_v_j->opposite()->vertex()->point();

        // Compute the norm of v_j -> v_i vector
        Vector_3 edge = position_v_i - position_v_j;
        double len = edge*edge;

        double weight = 0.0;
        if(len != 0.0)
            weight = 1.0/len;

        return weight;
	}
};

class Calculator
{
public:
	Calculator(int weightOption = 0);
	~Calculator();
	int compute(Polyhedron* mesh,std::list<ConstrianedVertex*>& vcMap,
						std::vector<double> &cValues, char* solverType, std::ofstream& logger);
private:
	int index_mesh_vertices(Polyhedron* mesh);
	void setScalarForVc(Polyhedron* mesh, std::list<ConstrianedVertex*>& vcMap);
	void initialize_system_from_mesh_vc(Matrix& A, Vector& Bu, std::list<ConstrianedVertex*>& vcMap);
	int setup_inner_vertex_relations(Matrix& A, Vector& Bu,
		                             Polyhedron* mesh, Vertex_handle vh, std::vector<double> &cValues);
	int setup_inner_vertex_relations(SymmetricMatrix& A, Vector& Bu,
		                             Polyhedron* mesh, Vertex_handle vh, std::vector<double> &cValues);

	void set_mesh_scalar_from_system(Polyhedron* mesh, Vector& Xu);

	int lu_symbolic_factor_solve(Matrix& A, Matrix& MA, Vector& B, Vector& X, std::ofstream& logger);
	int lu_symbolic_factor_solve1(Matrix& A, Matrix& MA, Vector& B, Vector& X, std::ofstream& logger);
	int ooc_lu_factor_solve(Matrix& A,Vector& B, Vector& X, std::ofstream& logger);
	int ll_symbolic_factor_solve(SymmetricMatrix &SA, SymmetricMatrix &SMA, Vector& B, Vector& X, std::ofstream& logger);

	void setupSystem(Polyhedron* mesh, Matrix&A, Vector&Bu, std::list<ConstrianedVertex*>& vcMap, std::vector<double> &cValues);
	void setupSymmetricSystem(Polyhedron* mesh, SymmetricMatrix&A, Vector&Bu, std::list<ConstrianedVertex*>& vcMap, std::vector<double> &cValues);
private:
	WeightOption* m_weightOption;
};

#endif//_CALCULATOR_H