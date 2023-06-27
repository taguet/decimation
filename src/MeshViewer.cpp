
//=============================================================================
//
//  CLASS MeshViewer - IMPLEMENTATION
//
//=============================================================================


//== INCLUDES =================================================================


#include <OpenMesh/Core/IO/MeshIO.hh>
#include "MeshViewer.hpp"
#include "gl.h"
#include <iostream>
#include <fstream>
#include "ComputingTools.hpp"


//== IMPLEMENTATION ========================================================== 


MeshViewer::
MeshViewer(const char* _title, int _width, int _height)
  : GlutExaminer(_title, _width, _height)
{

  //Requiert les normales des faces et des sommets pour la couche de base
  mesh.request_face_normals();
  mesh.request_vertex_normals();

  clear_draw_modes();
  add_draw_mode("Points");
  add_draw_mode("Wireframe");
  add_draw_mode("Hidden Line");
  add_draw_mode("Solid Flat");
  add_draw_mode("Solid Smooth");
  add_draw_mode("Normal Vectors");

  prev_id_draw_mode=4;
  default_id_draw_mode = prev_id_draw_mode;

  set_draw_mode(prev_id_draw_mode);
}


//-----------------------------------------------------------------------------

//Ouverture d'un fichier de CAO pour stocker les éléments du maillage initial dans l'un des deux maillages de surface.
bool
MeshViewer::
open_mesh(const char* _filename)
{
  // load mesh
  if (OpenMesh::IO::read_mesh(mesh, _filename))
  {
    // set center and radius
    Mesh::ConstVertexIter  v_it(mesh.vertices_begin()), 
                           v_end(mesh.vertices_end());
    //Mesh::Point            bbMin, bbMax;

    bbMin = bbMax = mesh.point(v_it);
    for (; v_it!=v_end; ++v_it)
    {
      bbMin.minimize(mesh.point(v_it));
      bbMax.maximize(mesh.point(v_it));
    }
    set_scene( (Vec3f)(bbMin + bbMax)*0.5, 0.5*(bbMin - bbMax).norm());

    // compute face & vertex normals
    mesh.update_normals();

    // update face indices for faster rendering
    update_face_indices(mesh,indices);

    // info
    cerr << mesh.n_vertices() << " vertices, "
	      << mesh.n_faces()    << " faces\n";

    return true;
  }

  else{
	cerr << "Error loading mesh from file " << _filename << endl;
  }

  return false;
}


bool MeshViewer::write_mesh(const char* _filename) {
	if (OpenMesh::IO::write_mesh(mesh, _filename)) {
		std::cerr << "Succesfully wrote mesh to file " << _filename;
		return true;
	}
	else {
		std::cerr << "Error writing mesh to file " << _filename;
	}
	return false;
}


//-----------------------------------------------------------------------------


void
MeshViewer::
update_face_indices(Mesh &mesh_, vector<unsigned int> &indices_)
{
  Mesh::ConstFaceIter f_it(mesh_.faces_sbegin()), f_end(mesh_.faces_end());
  Mesh::ConstFaceVertexIter  fv_it;

  indices_.clear();
  indices_.reserve(mesh_.n_faces()*3);

  for (; f_it!=f_end; ++f_it)
  {
    indices_.push_back((fv_it=mesh_.cfv_iter(f_it)).handle().idx());
    indices_.push_back((++fv_it).handle().idx());
    indices_.push_back((++fv_it).handle().idx());
  }
}


//-----------------------------------------------------------------------------

 
vector<unsigned int> MeshViewer::getIndicesMesh()
{
	return indices;
}


//-----------------------------------------------------------------------------

void 
MeshViewer::
draw(const string& _draw_mode)
{

  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT); 

  if (indices.empty())
  {
    GlutExaminer::draw(_draw_mode);
    return;
  }


  //Incompatible displays

  current_draw_mode = _draw_mode;


  if(current_draw_mode == "Wireframe")
  {
   
	glColor3f(1.0, 1.0, 1.0);

	glDisable(GL_LIGHTING);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	glEnableClientState(GL_VERTEX_ARRAY);
	GL::glVertexPointer(mesh.points());
	glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, &indices[0]);
	
	glDisableClientState(GL_VERTEX_ARRAY);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	prev_draw_mode = current_draw_mode;
	prev_id_draw_mode=get_draw_mode();
  }


  else if (current_draw_mode == "Hidden Line")
  {
    glDisable(GL_LIGHTING);

    glColor3f(1.0, 1.0, 1.0);
    glEnableClientState(GL_VERTEX_ARRAY);
    GL::glVertexPointer(mesh.points());


	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, &indices[0]);

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1.0, 1.0);
	glColor3f(0.0, 0.0, 0.0);
	glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, &indices[0]);
	glDisable(GL_POLYGON_OFFSET_FILL);

	prev_draw_mode = current_draw_mode;
	prev_id_draw_mode=get_draw_mode();
  }


  else if (current_draw_mode == "Solid Flat")
  {
    Mesh::ConstFaceIter f_it(mesh.faces_begin()), f_end(mesh.faces_end());
    Mesh::ConstFaceVertexIter fv_it;

	glEnable(GL_LIGHTING);
	glShadeModel(GL_FLAT);

	glColor3f(0.3, 0.3, 0.3);
	glShadeModel(GL_FLAT);
	glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT);
	glEnable(GL_COLOR_MATERIAL);

	glBegin(GL_TRIANGLES);
	for (; f_it!=f_end; ++f_it)
	{
	  GL::glNormal(mesh.normal(f_it));
	  fv_it = mesh.cfv_iter(f_it.handle()); 
	  GL::glVertex(mesh.point(fv_it));
	  ++fv_it;
	  GL::glVertex(mesh.point(fv_it));
	  ++fv_it;
	  GL::glVertex(mesh.point(fv_it));
	}
	glEnd();

	glDisable(GL_COLOR_MATERIAL);


	prev_draw_mode = current_draw_mode;
	prev_id_draw_mode=get_draw_mode();
  }


  else if (current_draw_mode == "Solid Smooth")
  {

	glEnable(GL_LIGHTING);

	glShadeModel(GL_SMOOTH);
	glColor3f(0.3, 0.3, 0.3);
	glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT);
	glEnable(GL_COLOR_MATERIAL);

	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);
	GL::glVertexPointer(mesh.points());
	GL::glNormalPointer(mesh.vertex_normals());

	glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, &indices[0]);

	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);


	prev_draw_mode = current_draw_mode;
	prev_id_draw_mode=get_draw_mode();
}


	//Optionnal display
	else if (current_draw_mode == "Normal Vectors")
	{
		map_normals();
		glDisable(GL_LIGHTING);

		glShadeModel(GL_SMOOTH);
		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_NORMAL_ARRAY);
		glEnableClientState(GL_COLOR_ARRAY);

		//Base surface
		GL::glVertexPointer(mesh.points());
		GL::glNormalPointer(mesh.vertex_normals());
		GL::glColorPointer(mesh.vertex_colors());
		glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, &indices[0]);


		glDisableClientState(GL_VERTEX_ARRAY);
		glDisableClientState(GL_NORMAL_ARRAY);
		glDisableClientState(GL_COLOR_ARRAY);


		prev_draw_mode = current_draw_mode;
		prev_id_draw_mode = get_draw_mode();
	}


   else if(current_draw_mode == "Points" ) // -----------------------------------------
   {
	 glDisable(GL_LIGHTING);

	 glColor3f(1.0f, 1.0f, 1.0f);

     glEnableClientState(GL_VERTEX_ARRAY);
     glVertexPointer(3, GL_FLOAT, 0, mesh.points());
     glDrawArrays( GL_POINTS, 0, mesh.n_vertices() );
     glDisableClientState(GL_VERTEX_ARRAY);
	 
	 glEnable(GL_LIGHTING);

	prev_draw_mode = current_draw_mode;
	prev_id_draw_mode=get_draw_mode();
  }
}

// Maps vertices normals to colors
void MeshViewer::map_normals() {
	for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
		Vec3f normal = mesh.normal(v_it);
		Vec3f to_2d = (normal + Vec3f(1.0, 1.0, 1.0)) / 2.0;
		Mesh::Color color = Mesh::Color(to_2d*255.0);
		mesh.set_color(v_it, color);
	}
}

//=============================================================================


void MeshViewer::keyboard(int key, int x, int y) 
{
  switch (key)
  {
    default:
    {
      GlutExaminer::keyboard(key, x, y);
      break;
    }
  }
}


