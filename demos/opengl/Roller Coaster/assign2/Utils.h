/*
  15-462 Computer Graphics I
  Assignment 3: Ray Tracer
  C++ Utility Classes and Functions
  Modified from rtark Oct 2007

  NOTE: You do not need to edit this file for this assignment 
  but may do so, especially to the Camera class

  This file defines the following:
  Vector Class
  Matrix Class
  Camera Class
*/

#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
  Vector Class - A float triplet class with several vector-operation functions

  This class can be used to simplify vector math programming
*/
class Vector
{
 public:
    float x, y, z, w;

    // -- Constructors & Destructors --
    // - Default Constructor - Initializes to Vector <0, 0, 0, [1]>
    Vector (void) : x (0.0f), y (0.0f), z (0.0f), w (1.0f) {}

    // - Parameter Constructor - Initializes to Vector <a, b, c, 1>
    Vector (float a, float b, float c) : x (a), y (b), z (c), w (1.0f) {}

    // - Parameter Constructor - Initializes to Vector <a, b, c, d>
    Vector (float a, float b, float c, float d) : x (a), y (b), z (c), w (d) {}

	Vector(float* a) : x(a[0]), y(a[1]), z(a[2]), w(1.0f) {}

    // - Default Destructor -
    ~Vector ()	{}

    // -- Utility Functions --
    // - Magnitude - Returns the Magnitude of the current Vector
    float Magnitude (void);

    // - Normalize - Normalizes to a Unit Vector (Scales to magnitude of 1)
    Vector Normalize (void);

    // - Scale - Scales the Vector by a factor
    Vector Scale (float scaleFactor);

    // - Dot - Calculates the Dot-Product between this and another Vector
    float Dot (Vector vec2);

    // - Cross - Returns the Cross-Product between this and another Vector
    Vector Cross (Vector vec2);

    // -- Operator Overloads to the class --
    // - Assignment Operator - Allows you to simply write "vec1 = vec2"
    Vector operator = (const Vector vec2);

    // NOTE: The following arithmetic operator overloads DO NOT change 
    // the value of the current vector
    // - Add Operator - Returns the sum of vectors
    Vector operator + (const Vector vec2);

    // - Subtract Operator - Returns the difference of vectors
    Vector operator - (const Vector vec2);

    // - Multiply Operator - Returns the vector scaled by a factor
    Vector operator * (const float scaleFactor);

    // - Divide Operator - Returns the vectors scaled by a factor
    Vector operator / (const float scaleFactor);

    // - Vector Multiply Operator -
    Vector operator * (const Vector scaleVector);
};

/*
  Matrix Class - A 4x4 Matrix class with several matrix-matrix and 
  matrix-vector operation functions

  This class can be used to simplify matrix & vector math programming
*/
class Matrix
{
 public:
    // Values are defined with a naming convention of row_column
    //    e.g. _23 is the 2nd Row and 3rd Column element
    float _11; float _12; float _13; float _14;
    float _21; float _22; float _23; float _24;
    float _31; float _32; float _33; float _34;
    float _41; float _42; float _43; float _44;

    // -- Constructors & Destructors --
    // - Default Constructor - Initializes to Identity Matrix
    Matrix (void)
	{
	    _12 = _13 = _14 = 0;
	    _21 = _23 = _24 = 0;
	    _31 = _32 = _34 = 0;
	    _41 = _42 = _43 = 0;
	    _11 = _22 = _33 = _44 = 1;
	}
    // - Default Destructor -
    ~Matrix (void)
	{}

    // -- Utitility Functions --
    // - Identity - Initializes the matrix back to the Identity
    Matrix Identity (void);

    // - Transpose - Transposes the Matrix
    Matrix Transpose (void);

    // - Inverse - Calculates the Inverse of the 4x4 Matrix
    Matrix Inverse (void);

    // -- Operator Overloads to the class --
    // - Assignment Operator - Allows you to simply write "mat1 = mat2"
    Matrix operator = (Matrix b);

    // - Matrix-Vector Multiplication - 
    // Returns the Vector multiplied by this matrix
    Vector operator * (const Vector vec);

    // - Matrix-Scale Multiplication - 
    // Returns the Matrix multiplied by the scalar value
    Matrix operator * (const float scalar);

    // - Matrix-Matrix Multiplication - 
    // Returns the Matrix multiplied by the 2nd Matrix
    Matrix operator * (const Matrix mat); // Multiply Operator
};


/*
  Camera Class - Class to help keep track of camera movements

  This class keeps the camera vectors together and 
  adds a bit of functionality around it
  NOTE: This is unoptimized code for better readability
*/
class Camera
{
 public:
    Vector position;
    Vector target;
    Vector up;
    float fieldOfView;
    float nearClip, farClip;

    // -- Constructors & Destructors --
    // - Default Constructor - Initializes to Vector <0, 0, 0>
    Camera (void) {}

    // - Parameter Constructor - Initializes to Vector <a, b, c>
    Camera (Vector p, Vector t, Vector u) { position = p; target = t; up = u; }

    // - Default Destructor -
    ~Camera ()	{}

    // -- Accessor Functions --
    // - GetPosition - Returns the position vector of the camera
    Vector GetPosition (void);

    // - SetPosition - Sets the position vector of the camera
    void SetPosition (const Vector vec);

    // - GetTarget - Returns the target vector of the camera
    Vector GetTarget (void);

    // - SetTarget - Sets the target vector of the camera
    void SetTarget (const Vector vec);

    // - GetUp - Returns the up vector of the camera
    Vector GetUp (void);

    // - SetUp - Sets the up vector of the camera
    void SetUp (const Vector vec);

    // - GetFOV - Returns the Field of View
    float GetFOV (void);

    // - SetFOV - Sets of Field of View
    void SetFOV (float fov);

    // - GetNearClip - Returns the Near-Clip Distance
    float GetNearClip (void);

    // - SetNearClip - Sets the Near-Clip Distance
    void SetNearClip (float nc);

    // - GetFarClip - Returns the Far-Clip Distance
    float GetFarClip (void);

    // - SetNearClip - Sets the Near-Clip Distance
    void SetFarClip (float fc);
};

#endif // UTILS_H
