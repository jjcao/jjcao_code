====BUILD ON MAC OS X====

cc main.cc -o test -framework OpenGL -framework GLUT

====BUILD ON LINUX (e.g., Ubuntu)====

1) install graphics driver, should be by default
2) install C/C++ compiler, should be by default
3) install OpenGL and freeglut libraries
3.1) OpenGL: sudo apt-get install mesa-common-dev
3.2) freeglut: sudo apt-get install freeglut3-dev
3.3) use these #include "GL/freeglut.h" and #include "GL/gl.h"
3.4) gcc main.cc -o test -lGL -lGLU -lglut

====BUILD ON WINDOWS====

http://user.xmission.com/~nate/glut.html
or via Visual Studio