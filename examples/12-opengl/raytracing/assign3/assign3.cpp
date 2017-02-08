/*
CSCI 480
Assignment 3 Raytracer

Name: <Your name here>
*/

#define cimg_use_jpeg /* see <your jpeg-9 directory>/jpeglib.h for using libjpeg 9 */
#include "CImg.h" 
using namespace cimg_library;

#include <windows.h>
#include <stdlib.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <vector>
#include "Utils.h"

#define MAX_TRIANGLES 2000
#define MAX_SPHERES 10
#define MAX_LIGHTS 10

char *filename="screen.jpg";

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2
int mode=MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
float WIDTH(320*2);//640
float HEIGHT(240*2);//480
float ar = WIDTH / HEIGHT; // aspect ratio
//the field of view of the camera
//#define fov 60.0
#define fov 1.05 // 60 degreen in radians

//unsigned char buffer[HEIGHT][WIDTH][3];
CImg<unsigned char> buffer(WIDTH, HEIGHT, 1, 3);

struct Vertex
{
  float position[3];
  float color_diffuse[3];
  float color_specular[3];
  float normal[3];
  float shininess;
};

typedef struct _Triangle
{
  struct Vertex v[3];
} Triangle;

typedef struct _Sphere
{
  float position[3];
  float color_diffuse[3];
  float color_specular[3];
  float shininess;
  float radius;
} Sphere;

typedef struct _Light
{
  float position[3];
  float color[3];
} Light;

enum ObjType { NONE, SPHERE, TRIANGLE};
Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
float ambient_light[3];
//Material coefficients for each color:
float ka = 0.1;
float kd = 0.6;
float ks = 0.3;
// attenuation constant: 1/(a+bq+cq^2), q is the distance
float aa(1.0);
float ab(0.0);
float ac(1.0);


int num_triangles=0;
int num_spheres=0;
int num_lights=0;

void save_jpg()
{
	printf("Saving JPEG file: %s\n", filename);
	buffer.save(filename);
}

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

// does the ray hit the sphere?
// ray: p(t) = origin + dt;
// d is the normalized direction;
// t returns t for the nearest intersection point to origin
// t == -1 if it does not hit anything
float hit_ray(Vector& origin, Vector& d, Sphere& sphere)
{
	float t(-1.0);

	Vector p0(origin);
	Vector ce(sphere.position);// center
	float r = sphere.radius;
	Vector tmpv = p0 - ce;
	float b = 2 * d.Dot(tmpv);
	float c = tmpv.Dot(tmpv) - r * r;
	float tmpf = b*b - 4 * c;

	if (tmpf < 0) return t;

	tmpf = std::sqrt(tmpf);
	float t0 = 0.5*(-b + tmpf);
	float t1 = 0.5*(-b - tmpf);

	t = std::numeric_limits<float>::max();
	if (t0 > 0 && t0 < t) t = t0;
	if (t1 > 0 && t1 < t) t = t1;

	if (t > 0.9*std::numeric_limits<float>::max()) t = -1;
	return t;
}
// does the line from origin to target intersect any object?
// ray: p(t) = origin + dt;
// d is the normalized direction determined by: target - origin;
// t returns t for the nearest intersection point to origin
// if distance(origin, target) < distance(origin, p(t)), return -1;
// return t and object (type & id) of the hitted object
// t == -1 if it does not hit anything
// if bnearst == true, return the object (type & id) of the nearest hitted object also.
float intersect_line(Vector& origin, Vector& target, int & objid, ObjType& objtype, bool bnearst=false)
{
	float t(-1.0);	
	Vector d(target.x - origin.x, target.y - origin.y, target.z - origin.z);
	float distThre = d.Magnitude();
	if (distThre < 20 * std::numeric_limits<float>::epsilon()) return t;
	d.Normalize(); // direction

	float mint = std::numeric_limits<float>::max();
	for (int i = 0; i < num_spheres; ++i)
	{
		t = hit_ray(origin, d, spheres[i]);
		if (t < 0 || t > distThre) continue;
		if (t > mint) continue;

		mint = t;
		objid = i;
		objtype = SPHERE;
		t = mint;

		if (!bnearst) break;
	}

	// todo triangles

	return t;
}

std::vector<float> trace_pixel(unsigned int x, unsigned int y)
{
	std::vector<float> color(3, 0.0f);

	// level0: scale
	float fy = 2.0f * tan(fov*0.5) * y / HEIGHT; // assume that near plane: z = -1
	float fx = 2.0f * tan(fov*0.5) * ar * x / WIDTH;
	fy = fy - tan(fov*0.5); // assume that near plane: z = -1
	fx = fx - ar * tan(fov*0.5);
	
	// level1: sent out rays
	Vector origin(0,0,0); // origin
	Vector dirPrim(fx - origin.x, fy - origin.y, -1 - origin.z);
	dirPrim.Normalize(); // direction of the primary ray

	float mint = std::numeric_limits<float>::max();
	for (int i = 0; i < num_spheres; ++i)
	{
		// level2: intersection
		float t = hit_ray(origin, dirPrim, spheres[i]);
		if (t< 0 || t > mint) continue;
		mint = t;		

		//color[0] = 255.0; 
		Vector ip = origin + dirPrim*mint; //intersection point		

		// level3: shadow rays
		int objid(-1);
		ObjType objtype(NONE);
		Vector illu(0.0, 0.0, 0.0); // accumulate illumination
		for (int j = 0; j < num_lights; ++j)
		{
			Vector target(lights[j].position);

			float t = intersect_line(ip, target, objid, objtype);
			if (t < 0 || (objtype==SPHERE && objid == i))
			{
				// todo compute its shading by local phone model				
				Vector L(lights[j].color);
				Vector sc(spheres[i].position);// center
				Vector n = (ip - sc).Normalize();
				Vector dirLight = (target - ip).Normalize(); // actually the negative direction of lighting
				Vector v = dirPrim*(-1); // unit vector to camera
				Vector r = n*dirLight.Dot(n)*2-dirLight;// unit reflected vector;
				float tmpd = std::max<float>(0.0f, n.Dot(dirLight));
				float tmps = std::max<float>(0.0f, r.Dot(v));
				Vector I = L*kd*tmpd+ L*ks*tmps;
				float dist = (target - ip).Magnitude(); ;
				illu = illu + I / (aa + ab*t+ac*dist * dist);
			}
			// else the pixel is blocked, i.e. no light hit the pixel
		}
		
		Vector La(ambient_light);
		illu = illu + La*ka;
		if (illu.x > 1.0) illu.x = 1.0;
		if (illu.y > 1.0) illu.y = 1.0;
		if (illu.z > 1.0) illu.z = 1.0;

		color[0] = 255.0*illu.x;
		color[1] = 255.0*illu.y;
		color[2] = 255.0*illu.z;
		// level4: nonzero specular component todo
	}

	for (int i = 0; i < num_triangles; ++i)
	{
		// level2: intersection
		Vector p(origin);
		Vector p0(triangles[i].v[0].position);
		Vector p1(triangles[i].v[1].position);
		Vector p2(triangles[i].v[2].position);

		CImg<float> A(3, 3, 1, 1); CImg<float> b(1, 3, 1, 1);
		Vector tmpv(dirPrim*(-1.0));
		//tmpv = origin;
		int j = 0;
		A(j, 0) = tmpv.x; A(j, 1) = tmpv.y; A(j, 2) = tmpv.z;
		
		tmpv = p1 - p0; 
		//tmpv = origin;
		j = 1;
		A(j, 0) = tmpv.x; A(j, 1) = tmpv.y; A(j, 2) = tmpv.z;
		
		tmpv = p2 - p0;
		//tmpv = origin;
		j = 2;
		A(j, 0) = tmpv.x; A(j, 1) = tmpv.y; A(j, 2) = tmpv.z;
		//A(0, 0) = 1; A(1, 1) = 1; A(2, 2) = 1;


		tmpv = p - p0;
		j = 0;
		b(j,0) = tmpv.x; b(j, 1) = tmpv.y; b(j, 2) = tmpv.z;


		CImg<float> sol = b.solve(A);
		float t = sol(j, 0);
		float u = sol(j, 1);
		float v = sol(j, 2);

		if (t < 0 || t > mint) continue;
		if (u < 0 || v < 0) continue;
		if (u + v > 1) continue;

		mint = t;
		//color[0] = 255.0; //color[1] = 1.0; color[2] = 1.0;
		Vector ip = p + dirPrim*mint; //intersection point	

		// level3: shadow rays todo
		// level4: nonzero specular component todo
	}
	//color[0] = 255.0;
	return color;
}

//MODIFY THIS FUNCTION
void draw_scene()
{
  unsigned int x,y;
  //simple output
  for(x=0; x<WIDTH; x++)
  {
    glPointSize(2.0);  
    glBegin(GL_POINTS);
    for(y=0;y < HEIGHT;y++)
    {
		//plot_pixel(x,y,x%256,y%256,(x+y)%256);
		std::vector<float> color = trace_pixel(x, y);      
		plot_pixel(x, y, color[0], color[1], color[2]);
    }
    glEnd();
    glFlush();
  }
  printf("Done!\n"); fflush(stdout);
}

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
  glColor3f(((float)r)/256.f,((float)g)/256.f,((float)b)/256.f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
	//buffer[HEIGHT-y-1][x][0]=r;
	//buffer[HEIGHT-y-1][x][1]=g;
	//buffer[HEIGHT-y-1][x][2]=b;
	buffer(x, y, 0, 0) = r;
	buffer(x, y, 0, 1) = g;
	buffer(x, y, 0, 2) = b;
}

void plot_pixel(int x,int y,unsigned char r,unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  //if(mode == MODE_JPEG)
      plot_pixel_jpeg(x,y,r,g,b);
}



void parse_check(char *expected,char *found)
{
  if(_stricmp(expected,found))
    {
      printf("Expected '%s ' found '%s '\n",expected,found);
      printf("Parse error, abnormal abortion\n");
      exit(0);
    }

}

void parse_floats(FILE*file, char *check, float p[3])
{
  char str[100];
  fscanf_s(file,"%s",str,100);
  parse_check(check,str);
  fscanf_s(file,"%f %f %f",&p[0],&p[1],&p[2]);
  printf("%s %f %f %f\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE*file,float *r)
{
  char str[100];
  fscanf_s(file,"%s",str, 100);
  parse_check("rad:",str);
  fscanf_s(file,"%f",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE*file,float *shi)
{
  char s[100];
  fscanf_s(file,"%s",s, 100);
  parse_check("shi:",s);
  fscanf_s(file,"%f",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
	FILE *file;
	fopen_s(&file, argv, "r");
  int number_of_objects;
  char type[50];
  int i;
  Triangle t;
  Sphere s;
  Light l;
  fscanf_s(file,"%i",&number_of_objects);

  printf("number of objects: %i\n",number_of_objects);

  parse_floats(file,"amb:",ambient_light);

  for(i=0;i < number_of_objects;i++)
    {
      fscanf_s(file,"%s\n",type, 50);
      printf("%s\n",type);
      if(_stricmp(type,"triangle")==0)
	{

	  printf("found triangle\n");
	  int j;

	  for(j=0;j < 3;j++)
	    {
	      parse_floats(file,"pos:",t.v[j].position);
	      parse_floats(file,"nor:",t.v[j].normal);
	      parse_floats(file,"dif:",t.v[j].color_diffuse);
	      parse_floats(file,"spe:",t.v[j].color_specular);
	      parse_shi(file,&t.v[j].shininess);
	    }

	  if(num_triangles == MAX_TRIANGLES)
	    {
	      printf("too many triangles, you should increase MAX_TRIANGLES!\n");
	      exit(0);
	    }
	  triangles[num_triangles++] = t;
	}
      else if(_stricmp(type,"sphere")==0)
	{
	  printf("found sphere\n");

	  parse_floats(file,"pos:",s.position);
	  parse_rad(file,&s.radius);
	  parse_floats(file,"dif:",s.color_diffuse);
	  parse_floats(file,"spe:",s.color_specular);
	  parse_shi(file,&s.shininess);

	  if(num_spheres == MAX_SPHERES)
	    {
	      printf("too many spheres, you should increase MAX_SPHERES!\n");
	      exit(0);
	    }
	  spheres[num_spheres++] = s;
	}
      else if(_stricmp(type,"light")==0)
	{
	  printf("found light\n");
	  parse_floats(file,"pos:",l.position);
	  parse_floats(file,"col:",l.color);

	  if(num_lights == MAX_LIGHTS)
	    {
	      printf("too many lights, you should increase MAX_LIGHTS!\n");
	      exit(0);
	    }
	  lights[num_lights++] = l;
	}
      else
	{
	  printf("unknown type in scene description:\n%s\n",type);
	  exit(0);
	}
    }
  return 0;
}

void display()
{

}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  //hack to make it only draw once
  static int once=0;
  if(!once)
  {
      draw_scene();
      //if(mode == MODE_JPEG)
		save_jpg();
    }
  once=1;
}

int main (int argc, char ** argv)
{
  if (argc<2 || argc > 3)
  {  
    printf ("usage: %s <scenefile> [jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
    {
      mode = MODE_JPEG;
      filename = argv[2];
    }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}
