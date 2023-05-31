
//=============================================================================
//
//  CLASS CurvatureViewer
//
//=============================================================================


#ifndef VIEWERWIDGET_HH
#define VIEWERWIDGET_HH


//== INCLUDES =================================================================


#include "MeshViewer.hpp"
#include "ComputingTools.hpp"
#include "MeshTools.h"
#include <Windows.h>
#include <commdlg.h>


//== CLASS DEFINITION =========================================================

	      

class Viewer : public MeshViewer
{
public:
   
  /// default constructor
  Viewer(const char* _title, int _width, int _height);

  // destructor
  ~Viewer();

  /// open mesh
  virtual bool open_mesh(const char* _filename);

  //Getter
  OpenMesh::TriMesh_ArrayKernelT<> getMesh();

  
protected:

  virtual void init();

  virtual void draw(const std::string& _draw_mode);

  void store_initial_points();

  void reset_mesh();

  /// calculate vertex and edge weights
  void calc_mean_weights();

  /// calculate curvature per vertex
  void calc_mean_curvature();

  // calculate gaussian curvature per vertex
  void calc_gauss_curvature();

  /// set vertex color from vertex curvature
  void color_coding(OpenMesh::VPropHandleT<Mesh::Scalar> _curv);

  // color 1-ring of a given vertex
  void color_1_ring(Mesh::VertexHandle _vh);

  // easier access to vertex weights
  Mesh::Scalar& weight(Mesh::VertexHandle _vh) 
  { return mesh.property(vweight_, _vh); }

  // easier access to vertex curvature  
  Mesh::Scalar& curvature_me(Mesh::VertexHandle _vh)
  {
	  return mesh.property(v_mean_curvature_, _vh);
  }

  Mesh::Scalar& gauss_curvature(Mesh::VertexHandle _vh)
  {
	  return mesh.property(v_gauss_curvature_, _vh);
  }

  Vec3f& initial_coord(Mesh::VertexHandle _vh)
  {
	  return mesh.property(initial_coords, _vh);
  }


  // easier access to edge weights
  Mesh::Scalar& weight(Mesh::EdgeHandle _eh) 
  { return mesh.property(eweight_, _eh); }


  /// keyboard interaction
  void keyboard(int key, int x, int y);

  /// open file explorer
  void browse_meshes();

private:

  OpenMesh::VPropHandleT<Mesh::Scalar>  vweight_, v_mean_curvature_, v_gauss_curvature_;
  OpenMesh::EPropHandleT<Mesh::Scalar>  eweight_;
  OpenMesh::VPropHandleT<Vec3f> initial_coords;

  GLuint  textureID_;
  ComputingTools ctools;
  MeshTools mtools;

  bool isModified{ false };
  bool calledSmoothing{ false };
  int v_id{ 0 };
  int neighbour_offset{ 0 };

  TopologyGraph* graph{ nullptr };

  void draw_1_ring(const VertexHandle vh, Vec3f color);
};


//=============================================================================
#endif // VIEWERWIDGET_HH defined
//=============================================================================

