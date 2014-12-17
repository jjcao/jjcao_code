#include "SgpProp.h"
#include <CGAL/Cartesian.h>
#include <CGAL/Timer.h>
#include <vtkRungeKutta4.h>
#include <vtkStreamTracer.h>
#include <vtkTubeFilter.h>
#include <vtkPointData.h>
#include <vtkSphereSource.h>
#include <vtkIntArray.h>
#include <vtkCellData.h>
SgpProp::SgpProp(vtkRenderer* render, QObject *parent)
	: QObject(parent),m_render(render)
{
	m_lines=new c2vPolyLines(render,this);
	m_basePlaneActor = 0;
}

SgpProp::~SgpProp()
{
	if (m_lines)
	{
		delete m_lines;
		m_lines=0;
	}

	if(m_basePlaneActor)
	{
		m_render->RemoveViewProp(m_basePlaneActor);
		m_basePlaneActor->Delete();
	}
}

//create a plane passing the center, with the average of some normals as normal
//the plane created is stored in m_basePlane
//用：CGAL::Plane_3<Kernel>我们这里默认了所得到的边界点是有顺序的！！！
void SgpProp::createPlane(Vertex_handle &center, Polyhedron* mesh)
{
	std::cout << "  SgpProp::createPlane() begin!" <<std::endl;
	CGAL::Timer timer;
	timer.start();	

	std::list<Vertex_handle>& main_border = mesh->main_border();
	if (main_border.empty())
	{
		main_border = mesh->extract_longest_border();
	}

	std::list<Vector_3 > spokes;//a circular linked list
	for (std::list<Vertex_handle>::iterator it=main_border.begin(); it!=main_border.end(); ++it)
	{
		Vertex_handle vh = *it;
		spokes.push_back(vh->point() - center->point());
	}
	spokes.push_back(*(spokes.begin()) );

	std::list<Vector_3 >::iterator bVector = spokes.begin();
	std::list<Vector_3 >::iterator nVector = bVector;
	Vector_3 sum_norm(0,0,0);
	for(++nVector; nVector != spokes.end(); ++nVector,++bVector)
	{
		sum_norm = sum_norm + CGAL::cross_product(*bVector,*nVector);
	} 

	m_basePlaneNormal = sum_norm/(spokes.size()+1);
	
	std::cout << "....Time: " << timer.time() << " seconds." << std::endl<<std::endl;
}
//compute the cross point
Point_3 compute_cross_point(Plane_3 plane, Point_3 start, Point_3 end)
{
	Vector_3 normal = plane.orthogonal_vector();
	Vector_3 line_direction = end - start;
	Point_3  p= plane.point();
	double t;
	double a = (start.x() - p.x()) * normal.x() + (start.y() - p.y()) * normal.y() + (start.z() - p.z()) * normal.z();
	double b = line_direction.x() * normal.x() + line_direction.y() * normal.y() + line_direction.z() * normal.z();

	assert(b != 0);
	t = -a / b;

	return start + t * line_direction;
}
// compute one cross point between 1-ring of c_vh and cutPlane
bool vertexTo(Vertex_handle &c_vh, Halfedge_handle &c_hh, Vertex_handle center, 
			  Halfedge_vertex_circulator &optimal_start_spoke, 
			  const Plane_3 &cutPlane, Polyhedron* mesh, std::list<Point_3> &cross_points,
			  int nthTarget)
{
	bool result(false);
	Halfedge_vertex_circulator  hc;
	if ( c_vh == center )
		hc =  optimal_start_spoke;
	else
		hc =  c_vh->vertex_begin();

	Halfedge_vertex_circulator end  =  hc;			
	CGAL_For_all(hc,end)
	{
		Vertex_handle lvh = hc->opposite()->vertex();
		Vertex_handle rvh = hc->next()->vertex();
		Point_3 lp = lvh->point();
		Point_3 rp = rvh->point();
		int il = cutPlane.oriented_side(lp);
		int ir = cutPlane.oriented_side(rp);
		int tmp = il*ir;
		if ( tmp<0)//异侧
		{
			Halfedge_handle lrhh = hc->next()->next();
			if (lrhh->tag()==nthTarget)//在找这次的target的过程中已经用过了
			{
				continue;
			}
			else
			{
				lrhh->tag(nthTarget);
				lrhh->opposite()->tag(nthTarget);
			}

			if ( c_vh == center )
			{
				optimal_start_spoke = hc;// ++optimal_start_spoke;
			}

			Point_3 cp = compute_cross_point(cutPlane,lp,rp);
			cross_points.push_back(cp);
			
			c_hh = hc->next()->next()->opposite();
			c_vh = 0;
			result = true;
			break;
		}
		else if(tmp>0)//同侧
		{
			continue;
		}

		if (ir)//lp is on the cut plane
		{		
			if (lvh->tag()==nthTarget)//在找这次的target的过程中已经用过了
			{
				continue;
			}
			else
			{
				lvh->tag(nthTarget);
			}
			if ( c_vh == center)
			{
				optimal_start_spoke = hc;// ++optimal_start_spoke;
			}
			c_vh = lvh;			
		}
		else
		{	
			if (rvh->tag()==nthTarget)//在找这次的target的过程中已经用过了
			{
				continue;
			}
			else
			{
				rvh->tag(nthTarget);
			}
			if ( c_vh == center)
			{
				optimal_start_spoke = hc; ++optimal_start_spoke; //++optimal_start_spoke;
			}
			c_vh = rvh;			
		}				
		cross_points.push_back(c_vh->point());
		c_hh = 0;
		result = true;
		break;
	}
	
	return result;
}
//tvh: target vertex
//optimal_start_spoke: estimated start spoke vertex, which point to the center; 
//根据上一个tvh所用到的spoke，估计本次的出发spoke
//return length: -1:死循环
double straightest_geodesic_path(Vertex_handle &center, Vertex_handle &tvh,Halfedge_vertex_circulator &optimal_start_spoke,
								 Vector_3 norm, Polyhedron* mesh, std::list<Point_3> &cross_points, int nthTarget)
{
	double len(0);
	cross_points.clear();
	std::list<Point_3> bakCps;

	//Vector_3 c2tVector = tvh->point()-center->point();
	Vector_3 cNorm=CGAL::cross_product(center->point()-tvh->point(),norm);
	Plane_3 cutPlane(center->point(),cNorm);

	Vertex_handle c_vh = center;
	cross_points.push_back(c_vh->point());
	c_vh->tag(nthTarget);

	Halfedge_handle c_hh = 0;	
	if ( optimal_start_spoke==0)
	{
		optimal_start_spoke = c_vh->vertex_begin();
	}
	int bBreak(0);

	do
	{		
		if(c_vh==0)
		{
			Vector_3 v = cross_points.back()-tvh->point();
			double dot = std::sqrt(v*v);
			if(dot<0.000000001)//考虑涉入误差的影响
			{
				cross_points.pop_back();
				cross_points.push_back(tvh->point());
				break;
			}

			if (c_hh->is_border())
			{
				if (bBreak==1)//两边都没有交到，取距离近的
				{		
					Vector_3 v0 = cross_points.back() - tvh->point();					
					Vector_3 v1 = bakCps.back() - tvh->point();
					double d0 = std::sqrt(v0*v0);
					double d1 = std::sqrt(v1*v1);
					if(d0>d1)
					{
						cross_points.clear();
						cross_points.insert(cross_points.begin(), bakCps.begin(), bakCps.end());
					}
					break;
				}
				++bBreak;
				c_vh = center;
				c_hh = 0;
				bakCps.insert(bakCps.begin(), cross_points.begin(), cross_points.end());
				cross_points.erase(++cross_points.begin(), cross_points.end());
				continue;
			}

			Vertex_handle nvh = c_hh->next()->vertex();			
			Point_3 np = nvh->point();
			int in = cutPlane.oriented_side(np);
			if (in==0)
			{
				if ( nvh->tag() == nthTarget)
				{
					break;//return -1.0;
				}

				c_vh = nvh;
				c_vh->tag(nthTarget);				
				cross_points.push_back(np);
				c_hh = 0;
			}
			else
			{
				Point_3 cp;
				//left			
				Vertex_handle lvh = c_hh->vertex();
				Point_3 lp = lvh->point();
				int il = cutPlane.oriented_side(lp);			
				int tmp = il*in;
				if ( tmp<0)//切面的两侧，即和左边有交
				{
					Halfedge_handle lhh = c_hh->next();
					if ( lhh->tag() == nthTarget)
					{
						break;//return -1.0;
					}
					else
					{
						lhh->tag(nthTarget);
						lhh->opposite()->tag(nthTarget);
					}

					cp = compute_cross_point(cutPlane,lp,np);
					
					c_hh = c_hh->next()->opposite();				
				}
				else//和右边有交
				{
					Halfedge_handle rhh = c_hh->next()->next();
					if ( rhh->tag() == nthTarget)
					{
						break;//return -1.0;
					}
					else
					{
						rhh->tag(nthTarget);
						rhh->opposite()->tag(nthTarget);
					}

					c_hh = c_hh->next()->next();//
					cp = compute_cross_point(cutPlane,c_hh->vertex()->point(),np);
					c_hh = c_hh->opposite();				
				}	
				
				Vector_3 v = cp - tvh->point();
				double dot = std::sqrt(v*v);
				if(dot<0.000000001)//考虑涉入误差的影响
				{
					cross_points.push_back(tvh->point());
					c_hh =0;
					c_vh = tvh;
				}
				else
				{
					cross_points.push_back(cp);
					c_vh = 0;
				}
			}
		}
		else
		{
			if (mesh->is_border(c_vh))
			{
				c_vh = center;
				c_hh = 0;
				cross_points.erase(++cross_points.begin(), cross_points.end());
			}
			else
			{
				bool result = vertexTo(c_vh, c_hh, center, optimal_start_spoke, cutPlane, mesh,cross_points, nthTarget);
				if (!result) 
					break;
			}
		}
	}while(c_vh != tvh);
	
	std::list<Point_3>::iterator itl = cross_points.begin();
	std::list<Point_3>::iterator itr = itl; ++itr;
	for(; itr!=cross_points.end();++itl,++itr)
	{		
		len += std::sqrt( (*itr-*itl)*(*itr-*itl));
	}
	return len;
}
//计算poisson侧地距离，梯度方向
bool isVertex(Polyhedron *mesh,Point_3 &p,Vertex_handle &tvh)//ok
{
	for(Vertex_iterator pVertex = mesh->vertices_begin(); pVertex != mesh->vertices_end(); ++pVertex)
	{
		Vector_3 v_3=pVertex->point()-p;
		double len=v_3*v_3;
		//if (pVertex->point()==p)//精确
		if(len<0.0000000000001)//带误差
		{
			tvh=pVertex;
			return true;
		}
	}
	return false;
}
Vector_3 current_trangle_gradient(Halfedge_handle &h)//ok
{  
	Halfedge_handle hh=h;
	Vector_3 gradient(0.0,0.0,0.0);
	double area=1;
	for(int i=0;i<3;i++)//
	{
		Vector_3 vec=hh->vertex()->point()- hh->opposite()->vertex()->point();
		Line_3 l_3(hh->vertex()->point(),hh->opposite()->vertex()->point());
		double leng=std::sqrt(vec*vec);
		hh=hh->next();
		Vector_3 next_vec=hh->vertex()->point()- hh->opposite()->vertex()->point();
		Vector_3 vvv=CGAL::cross_product(vec,next_vec);
		area=0.5*std::sqrt(vvv*vvv);
		double next_u=hh->vertex()->s();     
		Plane_3 plan(hh->vertex()->point(),vec);
		Plane_3 plan1(l_3,hh->vertex()->point());

		CGAL::Object obj=CGAL::intersection(plan,plan1);
		Line_3 li_3=CGAL::object_cast<Line_3>(obj);
		Vector_3 v=li_3.to_vector();
		if (v*next_vec<0)
		{
			v=-v;
		}
		gradient=gradient+next_u*(leng/std::sqrt(v*v))*v;	
	}	
	return gradient/area;
}

Point_3 compute_cross_points(Point_3 &start,Vector_3 &gradient, Facet_handle &current_face)//?
{
	Point_3 s=start;
	Line_3 l(start, gradient);	
	Point_2 p_2(start.x(),start.y());
	Point_2 po_2(0,0);//start所在边对应的vertex_points
	Point_2 poi_2(0,0);//start对面的vertex_points,or (0,0)
	Vector_2 gradient_2(gradient.x(),gradient.y());
	Halfedge_handle intersect_h=current_face->facet_begin();
	Halfedge_handle start_at_h=current_face->facet_begin();
	Halfedge_facet_circulator h= current_face->facet_begin();
	Halfedge_facet_circulator hh= current_face->facet_begin();
	do {   
		//Line_3 line(h->opposite()->vertex()->point(),h->vertex()->point());
		Vector_3 start_h_o=h->opposite()->vertex()->point()-start;
		Vector_3 start_h=h->vertex()->point()-start;
		Vector_3 cross=CGAL::cross_product(start_h,start_h_o);
		double lengths=std::sqrt(cross*cross);

		// if(line.has_on(start)) 太精确了
		if(lengths<0.00000000001)//带误差
		{
			start_at_h=h;
			if(start==h->vertex()->point())
			{
				h++;
				hh=h;
				po_2=Point_2(h->vertex()->point().x(),h->vertex()->point().y());
				h++;
				poi_2=Point_2(h->vertex()->point().x(),h->vertex()->point().y());
				break;			
			}
			else{
				hh=h;
				po_2=Point_2(h->vertex()->point().x(),h->vertex()->point().y());
				h++;
				poi_2=Point_2(h->vertex()->point().x(),h->vertex()->point().y());
				break;
			}			
		}

	} while (++h!=current_face->facet_begin());

	Line_2 l_2(p_2,gradient_2);
	if (l_2.has_on_positive_side(po_2)+l_2.has_on_positive_side(poi_2)==1)
	{
		intersect_h=++hh;
	}
	else
	{
		intersect_h=--hh;

	}

	current_face=intersect_h->opposite()->facet();	

	Point_3 ren_point(start.x(),start.y()+100,start.z()+200);
	Plane_3 plan(intersect_h->vertex()->point(),intersect_h->opposite()->vertex()->point(),ren_point);
	CGAL::Object obj=CGAL::intersection(plan,l);
	start=CGAL::object_cast<Point_3>(obj);

	/*if (CGAL::assign(start, obj))
	{
	int i=start.x();
	}
	else if (CGAL::assign(l, obj))
	{
	int i=start.x();
	}
	else
	{
	int i=start.x();
	}*/
	Vector_3 v_d=start-s;
	double k=v_d*gradient;
	if (k<-0.000000001&&intersect_h->vertex()->s()<-1)//防止出现迂回，打折线 from big to small???
	{
		if (start_at_h->vertex()->s() > start_at_h->opposite()->vertex()->s())
		{
			return start_at_h->opposite()->vertex()->point();
		}
		else
			return start_at_h->vertex()->point();
	}
	if (k<-0.000000001)//from small to big
	{
		if (start_at_h->vertex()->s() > start_at_h->opposite()->vertex()->s())
		{
			return start_at_h->vertex()->point();
		}
		else
			return start_at_h->opposite()->vertex()->point();
	}

	return start ;
}
double  Poisson_geodesic_path(Vertex_handle &center, Vertex_handle &tvh,
							  Polyhedron* mesh, std::list<Point_3> &geodesic_path)
{

	Point_3 sourcepoint = center->point();//源点	
	Point_3 p=tvh->point();
	geodesic_path.push_back(p);	

	Vector_3 gradient(0.0,0.0,0.0);	
	Vector_3 con_gradient=-gradient;
	Facet_handle current_face;//=tvh->vertex_begin()->opposite()->facet();
	double currentdistance=100000;
	double distance=1000000;
	double edge_length=0;
	int step=0;	
	if (tvh->s()<center->s())//from small to big
	{
		do
		{  	
			if(isVertex(mesh,p,tvh))//当前点是否是网格顶点
			{     
				double n_v_u=1000;
				double n_v_n_u=1000;
				double v_u=tvh->s();
				double max_u=-1000;
				Point_3 max_p;
				Halfedge_handle max_h;//

				Halfedge_vertex_circulator h_v=tvh->vertex_begin();
				Halfedge_vertex_circulator end=h_v;
				CGAL_For_all(h_v, end)
				{
					Halfedge_handle h=h_v->opposite();

					gradient=current_trangle_gradient(h);

					//判断
					Vector_2 v_2(gradient.x(),gradient.y());
					Point_2 p_2(p.x(),p.y());
					Line_2 l_2(p_2,v_2);
					Halfedge_handle h_h=h;
					n_v_u=h_h->vertex()->s();

					Vector_3 vv=h_h->vertex()->point()-h_h->opposite()->vertex()->point();
					edge_length=std::sqrt(vv*vv);
					n_v_u=(n_v_u-v_u)/edge_length;//单位化

					h_h=h_h->next();
					n_v_n_u=h_h->vertex()->s();
                    Vector_3 vvv=h_h->vertex()->point()-h_h->opposite()->vertex()->point();
					double edge_length1=std::sqrt(vvv*vvv);
					n_v_n_u=(n_v_n_u-v_u)/edge_length1;//单位化

					Point_2 po_2(h->vertex()->point().x(),h->vertex()->point().y());
					Point_2 poi_2(h_h->vertex()->point().x(),h_h->vertex()->point().y());
					/*Vector_3 vv=h_h->vertex()->point()-h_h->opposite()->vertex()->point();
					edge_length=std::sqrt(vv*vv);*/

					if (max_u < n_v_u)//判断u值最大的点和对应的边，出现误差时用
					{
						max_u=n_v_u;
						max_p=h->vertex()->point();
						max_h=h;
					}

					if (l_2.has_on_positive_side(po_2)+l_2.has_on_positive_side(poi_2)==1&&l_2.has_on(po_2)+l_2.has_on(poi_2)==0)//前进方向，可能有问题???
					{
						/*if (n_v_u>v_u)
						{
							v_u=n_v_u;	
                            con_gradient=gradient;
							current_face=h->facet();
						}	*/	
						if (0.5*(n_v_u+n_v_n_u)>v_u)
						{
							v_u=0.5*(n_v_u+n_v_n_u);	
							con_gradient=gradient;
							current_face=h->facet();
						}	
					}				
				}



				if (v_u==tvh->s())//如果出现误差，梯度不经过任何三角形
				{

					p=max_p;
					geodesic_path.push_back(p);

					currentdistance=std::sqrt((p.x()-sourcepoint.x())*(p.x()-sourcepoint.x()))+
						std::sqrt((p.y()-sourcepoint.y())*(p.y()-sourcepoint.y()))+
						std::sqrt((p.z()-sourcepoint.z())*(p.z()-sourcepoint.z()));//当前点距离源点距离

					if(step>200||currentdistance<edge_length) 
						break;
					step++;
					continue;
				}
			}
			else
			{

				Halfedge_handle h_c=current_face->facet_begin();
				gradient=current_trangle_gradient(h_c);
				con_gradient=gradient;
				
			}

			p=compute_cross_points(p,con_gradient,current_face);//下个路径点和下个三角面（非网格点）


			currentdistance=std::sqrt((p.x()-sourcepoint.x())*(p.x()-sourcepoint.x()))+
				std::sqrt((p.y()-sourcepoint.y())*(p.y()-sourcepoint.y()))+
				std::sqrt((p.z()-sourcepoint.z())*(p.z()-sourcepoint.z()));

			if((currentdistance==distance||currentdistance>distance)&&distance<edge_length*2)//如果计算离源点足够近可以停止（因为有误差，无法到达源点）
				break;
			distance=currentdistance;
			geodesic_path.push_back(p);

			if(step>200) 
				break;
			step++;
		}while(currentdistance>0.001);


	}
	else{
		do
		{  	
			if(isVertex(mesh,p,tvh))//当前点是否是网格顶点
			{     
				double n_v_u=0;
				double n_v_n_u=0;
				double v_u=tvh->s();
				double min_u=10000;
				Point_3 min_p;
				Halfedge_handle min_h;//

				Halfedge_vertex_circulator h_v=tvh->vertex_begin();
				Halfedge_vertex_circulator end=h_v;
				CGAL_For_all(h_v, end)
				{
					Halfedge_handle h=h_v->opposite();

					gradient=current_trangle_gradient(h);

					//判断
					Vector_2 v_2(gradient.x(),gradient.y());
					Point_2 p_2(p.x(),p.y());
					Line_2 l_2(p_2,v_2);
					Halfedge_handle h_h=h;
					n_v_u=h_h->vertex()->s();

					Vector_3 vv=h_h->vertex()->point()-h_h->opposite()->vertex()->point();
					edge_length=std::sqrt(vv*vv);
					n_v_u=(n_v_u-v_u)/edge_length;//单位化

					h_h=h_h->next();
					n_v_n_u=h_h->vertex()->s();
                    Vector_3 vvvv=h_h->vertex()->point()-h_h->opposite()->vertex()->point();
					double edge_length2=std::sqrt(vvvv*vvvv);
					n_v_n_u=(n_v_n_u-v_u)/edge_length2;//单位化


					Point_2 po_2(h->vertex()->point().x(),h->vertex()->point().y());
					Point_2 poi_2(h_h->vertex()->point().x(),h_h->vertex()->point().y());
					//Vector_3 vv=h_h->vertex()->point()-h_h->opposite()->vertex()->point();
					//edge_length=std::sqrt(vv*vv);

					if (min_u>n_v_u)//判断u值最小的点和对应的边，出现误差时用
					{
						min_u=n_v_u;
						min_p=h->vertex()->point();
						min_h=h;
					}

					if (l_2.has_on_positive_side(po_2)+l_2.has_on_positive_side(poi_2)==1&&l_2.has_on(po_2)+l_2.has_on(poi_2)==0)//前进方向，可能有问题???
					{
						/*if (n_v_u<v_u)
						{
							v_u=n_v_u;	
							con_gradient=-gradient;
							current_face=h->facet();
						}	*/	
						if (0.5*(n_v_u+n_v_n_u)<v_u)
						{
							v_u=0.5*(n_v_u+n_v_n_u);	
							con_gradient=-gradient;
							current_face=h->facet();
						}
						
					}				
				}



				if (v_u==tvh->s())//如果出现误差，梯度不经过任何三角形
				{

					p=min_p;
					geodesic_path.push_back(p);

					currentdistance=std::sqrt((p.x()-sourcepoint.x())*(p.x()-sourcepoint.x()))+
						std::sqrt((p.y()-sourcepoint.y())*(p.y()-sourcepoint.y()))+
						std::sqrt((p.z()-sourcepoint.z())*(p.z()-sourcepoint.z()));//当前点距离源点距离

					if(step>200||currentdistance<edge_length) 
						break;
					step++;
					continue;
				}
			}
			else
			{

				Halfedge_handle h_c=current_face->facet_begin();
				gradient=current_trangle_gradient(h_c);
				con_gradient=-gradient;
			}

			p=compute_cross_points(p,con_gradient,current_face);//下个路径点和下个三角面（非网格点）


			currentdistance=std::sqrt((p.x()-sourcepoint.x())*(p.x()-sourcepoint.x()))+
				std::sqrt((p.y()-sourcepoint.y())*(p.y()-sourcepoint.y()))+
				std::sqrt((p.z()-sourcepoint.z())*(p.z()-sourcepoint.z()));

			if((currentdistance==distance||currentdistance>distance)&&distance<edge_length*2)//如果计算离源点足够近可以停止（因为有误差，无法到达源点）
				break;
			distance=currentdistance;
			geodesic_path.push_back(p);

			if(step>200) 
				break;
			step++;
		}while(currentdistance>0.001);

	}
	

	geodesic_path.push_back(sourcepoint);

	double len(0);
	std::list<Point_3>::iterator itl = geodesic_path.begin();
	std::list<Point_3>::iterator itr = itl; ++itr;
	for(; itr!=geodesic_path.end();++itl,++itr)
	{		
		len += std::sqrt( (*itr-*itl)*(*itr-*itl));
	}
	return len;
}

//create all straightest path actor from center to every vertices in vhs by cutting plane,
//and save relative sgp distance to them
void SgpProp::createSgp(Vertex_handle &center, Polyhedron* mesh, std::list<Vertex_handle>& vhs)
{
	//init for label
	for(Vertex_iterator it = mesh->vertices_begin(); it != mesh->vertices_end(); ++it)
	{		
		it->tag(-1);
	}
	for(Halfedge_iterator it = mesh->halfedges_begin(); it != mesh->halfedges_end(); ++it)
	{		
		it->tag(-1);
	}	
    m_sgps.clear();

	//
	std::list<double> lens;
	std::list<Point_3> cross_points;    
	Halfedge_vertex_circulator  optimal_start_spoke;
	
	int i(0);
	for(std::list<Vertex_handle>::iterator pVertex = vhs.begin(); pVertex != vhs.end(); ++pVertex,++i)
	{		
		Vertex_handle vh = *pVertex;
		double len = straightest_geodesic_path(center, vh, optimal_start_spoke, 
											   m_basePlaneNormal, mesh, cross_points, i);
		vh->s(len);
		lens.push_back(len);
		m_sgps.push_back(cross_points);
	}
}
void SgpProp::createPgp(Vertex_handle &center, Polyhedron* mesh, std::list<Vertex_handle>& vhs)
{
	//init for label
	for(Vertex_iterator it = mesh->vertices_begin(); it != mesh->vertices_end(); ++it)
	{		
		it->tag(-1);
	}
	for(Halfedge_iterator it = mesh->halfedges_begin(); it != mesh->halfedges_end(); ++it)
	{		
		it->tag(-1);
	}	
	m_sgps.clear();

	//
	std::vector<double> lengs;
	lengs.reserve(vhs.size());
	   
	for(std::list<Vertex_handle>::iterator pVertex = vhs.begin(); pVertex != vhs.end(); ++pVertex)
	{	
		std::list<Point_3> P_geodesic_path; 
		Vertex_handle vh = *pVertex;
		double len = Poisson_geodesic_path(center, vh, mesh, P_geodesic_path);
		vh->u(len);
		lengs.push_back(len);
		m_sgps.push_back(P_geodesic_path);
	}

	/*int i = 0;
	for(std::list<Vertex_handle>::iterator pVertex = vhs.begin(); pVertex != vhs.end(); ++pVertex,++i)
	{	
		Vertex_handle vh = *pVertex;
		vh->s( lengs[i] );
	}*/
}
//create all straightest path actor from center to every border vertices by cutting plane,
//and save relative sgp distance to them. 
void SgpProp::createBorderSgp(Vertex_handle &center, Polyhedron* mesh)
{
	clearActors();
	
	std::cout << "  SgpProp::createBorderSgp() begin!" <<std::endl;
	CGAL::Timer timer;
	timer.start();
	createSgp(center, mesh, mesh->main_border());
	std::cout << "....Time: " << timer.time() << " seconds." << std::endl<<std::endl;

	//
	m_lines->loadData(m_sgps);
	m_lines->setPickable(false);
	emit updateVTK();
}
void SgpProp::createBorderPgp(Vertex_handle &center, Polyhedron* mesh)
{
	clearActors();
	createPgp(center, mesh, mesh->main_border());
	//
	m_lines->loadData(m_sgps);
	m_lines->setPickable(false);
	emit updateVTK();
}
//create all straightest path actor from center to every vertices by cutting plane,
//and save relative sgp distance to them, if idx =-1.
//If idx>0, create a straightest path from center to the vertex, whose index = idx;
void SgpProp::createSgp(Vertex_handle &center, Polyhedron* mesh, int idx)
{
	clearActors();

	std::list<Vertex_handle> vhs;
	if ( idx<0)
	{
		for(Vertex_iterator vi = mesh->vertices_begin(); vi != mesh->vertices_end(); ++vi)
		{
			vhs.push_back(vi);
		}
	}
	else
	{
		vhs.push_back(mesh->getVertexHandle(idx));
	}
	createSgp(center, mesh, vhs);

	if ( idx<0)
		return;

	m_lines->loadData(m_sgps);
	m_lines->setPickable(false);
	emit updateVTK();
}
void SgpProp::createPgp(Vertex_handle &center, Polyhedron* mesh, int idx)
{
	clearActors();

	std::list<Vertex_handle> vhs;
	if ( idx<0)
	{
		for(Vertex_iterator vi = mesh->vertices_begin(); vi != mesh->vertices_end(); ++vi)
		{
			vhs.push_back(vi);
		}
	}
	else
	{
		vhs.push_back(mesh->getVertexHandle(idx));
	}
	createPgp(center, mesh, vhs);

	if ( idx<0)
		return;

	m_lines->loadData(m_sgps);
	m_lines->setPickable(false);
	emit updateVTK();
}
void SgpProp::clearActors()
{
	delete m_lines;
	m_lines=new c2vPolyLines(m_render,this);

	for(std::list<vtkActor*>::iterator itor = m_actors.begin();itor!=m_actors.end(); ++itor)
	{
		vtkActor* a = *itor;
		m_render->RemoveViewProp(a);
		a->Delete();
	}
	m_actors.clear();
}
void SgpProp::setVisibility(bool in)
{
	m_lines->setVisibility(in);
}
void SgpProp::setBasePlaneVisibility(bool in)
{
}


void SgpProp::createRK(Vertex_handle &center, Polyhedron* mesh, vtkPolyData* meshData, int idx)
{
	clearActors();

	std::list<Vertex_handle> vhs;
	if ( idx<0)
	{
		for(Vertex_iterator vi = mesh->vertices_begin(); vi != mesh->vertices_end(); ++vi)
		{
			vhs.push_back(vi);
		}
	}
	else
	{
		vhs.push_back(mesh->getVertexHandle(idx));
	}
	createRK(center, mesh, meshData, vhs);

	emit updateVTK();
}
void SgpProp::createRK(Vertex_handle &center, Polyhedron* mesh, vtkPolyData* meshData, std::list<Vertex_handle>& vhs)
{
	//vtkRungeKutta4* integ = vtkRungeKutta4::New();
	for(std::list<Vertex_handle>::iterator pVertex = vhs.begin(); pVertex != vhs.end(); ++pVertex)
	{		
		Vertex_handle vh = *pVertex;		
		Point_3 p = vh->point();

		vtkStreamTracer* streamer = vtkStreamTracer::New();
		streamer->SetInput(meshData);
		streamer->SetStartPosition(p.x(), p.y(), p.z());
		//streamer->SetSourceConnection(sphere->GetOutputPort());
		streamer->SetMaximumPropagation(vtkStreamTracer::CELL_LENGTH_UNIT,500);
		streamer->SetInitialIntegrationStep(vtkStreamTracer::CELL_LENGTH_UNIT,0.1);
		streamer->SetMinimumIntegrationStep(vtkStreamTracer::CELL_LENGTH_UNIT,0.1);
		streamer->SetMaximumIntegrationStep(vtkStreamTracer::CELL_LENGTH_UNIT,1);
		streamer->SetMaximumError(1.0);
		streamer->SetIntegrationDirectionToBoth();
		//streamer->SetIntegrationDirectionToBackward();
		//streamer->SetIntegrator(integ);
		streamer->Modified();
		streamer->Update();

		vtkPolyData* sline = streamer->GetOutput();
		sline->Modified();
		sline->Update();

		vtkStreamTracer::ReasonForTermination RFT;
vtkIntArray *Int = (vtkIntArray *)streamer->GetOutput()->GetCellData()->GetArray("ReasonForTermination");
QString str;
if(Int && Int->GetSize()>0) {
   RFT = (vtkStreamTracer::ReasonForTermination)Int->GetValue(0);
   switch (RFT)
   {
   case vtkStreamTracer::OUT_OF_DOMAIN:
	   str = "OUT_OF_DOMAIN";
	   break;
   case vtkStreamTracer::NOT_INITIALIZED:
	   str = "NOT_INITIALIZED";
	   break;
   case vtkStreamTracer::UNEXPECTED_VALUE:
	   str = "UNEXPECTED_VALUE";
	   break;
   case vtkStreamTracer::OUT_OF_TIME:
	   str = "OUT_OF_TIME";
	   break;
   case vtkStreamTracer::OUT_OF_STEPS:
	   str = "OUT_OF_STEPS";
	   break;
   case vtkStreamTracer::STAGNATION:
	   str = "STAGNATION";
	   break;
   default:
	   str = "UNKNOWN";
   }   
} else 
{	
   str = "Error getting streamline properties!";
}
emit setStatusBarInfo(str);

		//vtkTubeFilter* streamTube = vtkTubeFilter::New();
		//streamTube->SetInputConnection(streamer->GetOutputPort());
		//streamTube->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_CELLS,vtkDataSetAttributes::VECTORS);
		//streamTube->SetRadius( 0.02);
		//streamTube->SetNumberOfSides( 12);
		//streamTube->SetVaryRadiusToVaryRadiusByVector();

		vtkPolyDataMapper *mapStreamTube = vtkPolyDataMapper::New();
		mapStreamTube->SetInputConnection(streamer->GetOutputPort());
		mapStreamTube->SetScalarRange(meshData->GetPointData()->GetScalars()->GetRange());
		vtkActor *streamTubeActor = vtkActor::New();
		streamTubeActor->SetMapper(mapStreamTube);
		streamTubeActor->GetProperty()->BackfaceCullingOn();

		m_actors.push_back(streamTubeActor);
		m_render->AddActor(streamTubeActor);
	}
}
void SgpProp::createBorderRK(Vertex_handle &center, Polyhedron* mesh, vtkPolyData* meshData)
{
	clearActors();
	createRK(center, mesh, meshData, mesh->main_border());
	emit updateVTK();
}