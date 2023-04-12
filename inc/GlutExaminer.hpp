
//=============================================================================
//
//  CLASS GlutExaminer
//
//=============================================================================


#ifndef GLUTEXAMINER_HH
#define GLUTEXAMINER_HH


//== INCLUDES =================================================================

#include "GlutViewer.hpp"
#include <OpenMesh/Core/Math/VectorT.hh>

#include <string>
#include <vector>

using namespace OpenMesh;


//== CLASS DEFINITION =========================================================

	      

class GlutExaminer : public GlutViewer
{
public:
   
  GlutExaminer(const char* _title, int _width, int _height);

  void   set_scene(const Vec3f& _center, float _radius);
  void   view_all();
  double measure_fps();


protected:

  virtual void init();
  virtual void draw(const string& _draw_mode);


  // overloaded glut functions
  virtual void motion(int x, int y);
  virtual void mouse(int button, int state, int x, int y);
  virtual void reshape(int w, int h); 
  virtual void keyboard(int key, int x, int y);


  // updates projection matrix
  void update_projection_matrix();
  // translate the scene and update modelview matrix
  void translate(const Vec3f& _trans);
  // rotate the scene (around its center) and update modelview matrix
  void rotate(const Vec3f& _axis, float _angle);


  // virtual trackball: map 2D screen point to unit sphere
  bool map_to_sphere(const Vec2i& _point, Vec3f& _result);


  // mouse processing functions
  void rotation(int x, int y);
  void translation(int x, int y);
  void zoom(int x, int y);


protected:

  // scene position and dimension
  Vec3f    center_;
  float    radius_;


  // projection parameters
  float    near_, far_, fovy_;


  // OpenGL matrices
  double   projection_matrix_[16],
           modelview_matrix_[16];

  
  // trackball helpers
  Vec2i    last_point_2D_;
  Vec3f    last_point_3D_;
  bool     last_point_ok_;
  bool     button_down_[10];
};


//=============================================================================
#endif // GLUTEXAMINER_HH defined
//=============================================================================

