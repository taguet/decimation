
//=============================================================================
//
//  CLASS MeshViewerWidget
//
//=============================================================================


#ifndef MESH_VIEWER_WIDGET_HH
#define MESH_VIEWER_WIDGET_HH


//== INCLUDES =================================================================


#include "GlutExaminer.hpp"
#include <OpenMesh/Core/Mesh/Types/TriMesh_ArrayKernelT.hh>

//== CLASS DEFINITION =========================================================

#define EPS 10e-4

class MeshViewer : public GlutExaminer
{

protected:

  typedef OpenMesh::TriMesh_ArrayKernelT<>  Mesh;

public:
  
  /// default constructor
  MeshViewer(const char* _title, int _width, int _height);

  /// open mesh
  virtual bool open_mesh(const char* _filename);

  /// update buffer with face indices
  void update_face_indices(Mesh &mesh_, vector<unsigned int> &indices_);

  /// draw the scene
  virtual void draw(const string& _draw_mode);

  /// get the indices from the viewer
  vector<unsigned int> getIndicesMesh();

  /// keyboard interaction
  void keyboard(int key, int x, int y);

  void map_normals();
 
protected:

  Mesh mesh;

  Mesh::Point bbMin, bbMax;
  
  vector<unsigned int> indices;
  
  unsigned int prev_id_draw_mode;
  unsigned int default_id_draw_mode;
  string prev_draw_mode;
  string current_draw_mode;

};


//=============================================================================
#endif // MESH_VIEWER_WIDGET_HH defined
//=============================================================================

