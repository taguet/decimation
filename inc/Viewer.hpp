
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


  // easier access to edge weights
  Mesh::Scalar& weight(Mesh::EdgeHandle _eh) 
  { return mesh.property(eweight_, _eh); }



private:

  OpenMesh::VPropHandleT<Mesh::Scalar>  vweight_, v_mean_curvature_, v_gauss_curvature_;
  OpenMesh::EPropHandleT<Mesh::Scalar>  eweight_;

  GLuint  textureID_;
  ComputingTools ctools;
};


//=============================================================================
#endif // VIEWERWIDGET_HH defined
//=============================================================================

