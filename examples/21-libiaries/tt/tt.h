//#ifndef ADD_H
//#define ADD_H
//int add(int a, int b);
//#endif  // ADD_H

#ifndef ADD_H
#define ADD_H
#ifdef BUILD_DLL
#define PORT_DLL __declspec(dllexport)
#else
#define PORT_DLL __declspec(dllimport)
#endif
int PORT_DLL add(int a, int b);
#endif  // ADD_H