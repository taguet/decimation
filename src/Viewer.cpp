
//=============================================================================
//
//  CLASS Viewer - IMPLEMENTATION
//
//=============================================================================


//== INCLUDES =================================================================

#include "Viewer.hpp"
#include <vector>
#include <float.h>


//== IMPLEMENTATION ========================================================== 


Viewer::Viewer(const char* _title, int _width, int _height): MeshViewer(_title, _width, _height) {

  mesh.request_vertex_colors();
	
  mesh.add_property(v_mean_curvature_);
  mesh.add_property(v_gauss_curvature_);

  mesh.add_property(vweight_);
  mesh.add_property(eweight_);

  add_draw_mode("Mean Curvature");
  add_draw_mode("Gaussian Curvature");
  add_draw_mode("Reflection Lines");
  add_draw_mode("Vertex 1-Ring");

  init();
}

//-----------------------------------------------------------------------------

Viewer::~Viewer() {
  if (glIsTexture(textureID_))  
    glDeleteTextures( 1, &textureID_);
}

//-----------------------------------------------------------------------------

OpenMesh::TriMesh_ArrayKernelT<> Viewer::getMesh() {
  return mesh;
}

//-----------------------------------------------------------------------------

void Viewer::init() {

  // base class first
  MeshViewer::init();
  ctools = ComputingTools();
	
	
  // generate checkerboard-like image
  GLubyte tex[256*256*3], *tp=tex;
  for (int x=0; x<256; ++x)
    for (int y=0; y<256; ++y)
      if (((x+2)/4 % 10) == 0 || ((y+2)/4 % 10) == 0)
      {
				*(tp++) = 0;
				*(tp++) = 0;
				*(tp++) = 0;
      }
			else
			{
				*(tp++) = 255;
				*(tp++) = 255;
				*(tp++) = 255;
			}
				
				
	// generate texture
	if (!glIsTexture(textureID_))
		glGenTextures(1, &textureID_);
  glBindTexture(GL_TEXTURE_2D, textureID_);
	
	
  // copy texture to GL
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexImage2D(GL_TEXTURE_2D, 0, 3, 256, 256, 0, GL_RGB, GL_UNSIGNED_BYTE, tex);

}



//-----------------------------------------------------------------------------


bool Viewer::open_mesh(const char* _filename) {

	  // load mesh
	  if (MeshViewer::open_mesh(_filename))
	  {
		ctools.set_mesh(mesh);
		ctools.calc_princ_curvatures();
		
		//Compute all the information about curvature
		calc_mean_curvature();
		calc_gauss_curvature();

		glutPostRedisplay();
		return true;
	  }
	  return false;
}


//-----------------------------------------------------------------------------


void Viewer::calc_mean_weights() {

	  Mesh::VertexIter        v_it, v_end(mesh.vertices_end());
	  Mesh::EdgeIter          e_it, e_end(mesh.edges_end());
	  Mesh::VertexFaceIter    vf_it;
	  Mesh::FaceVertexIter    fv_it;
	  Mesh::HalfedgeHandle    h0, h1, h2;
	  Mesh::VertexHandle      v0, v1;
	  Mesh::Point             p0, p1, p2, d0, d1;
	  Mesh::Scalar            w, area, b(0.99);


	  for (e_it=mesh.edges_begin(); e_it!=e_end; ++e_it)
	  {
		w  = 0.0;
		
		h0 = mesh.halfedge_handle(e_it.handle(), 0);
		v0 = mesh.to_vertex_handle(h0);
		p0 = mesh.point(v0);
		
		h1 = mesh.halfedge_handle(e_it.handle(), 1);
		v1 = mesh.to_vertex_handle(h1);
		p1 = mesh.point(v1);
		
			if (!mesh.is_boundary(h0))
			{
				h2 = mesh.next_halfedge_handle(h0);
				p2 = mesh.point(mesh.to_vertex_handle(h2));
				d0 = (p0 - p2).normalize();
				d1 = (p1 - p2).normalize();
				w += 1.0 / tan(acos(std::max(-b, std::min(b, (d0|d1)))));
			}
		
			if (!mesh.is_boundary(h1))
			{
				h2 = mesh.next_halfedge_handle(h1);
				p2 = mesh.point(mesh.to_vertex_handle(h2));
				d0 = (p0 - p2).normalize();
				d1 = (p1 - p2).normalize();
				w += 1.0 / tan(acos(std::max(-b, std::min(b, (d0|d1)))));
			}
		
			// force weights to be non-negative for higher robustness
			w = std::max(w, 0.0f);
		
		weight(e_it) = w;

	  }
	
	  //int k = 0; cout << endl;

	  for (v_it=mesh.vertices_begin(); v_it!=v_end; ++v_it)
	  {
		area = 0.0;
		
		for (vf_it=mesh.vf_iter(v_it); vf_it; ++vf_it)
		{
		  fv_it = mesh.fv_iter(vf_it);

		  const Mesh::Point& P = mesh.point(fv_it);  ++fv_it;
		  const Mesh::Point& Q = mesh.point(fv_it);  ++fv_it;
		  const Mesh::Point& R = mesh.point(fv_it);
			
		  area += ((Q-P)%(R-P)).norm() * 0.5f * 0.3333f;
		}
		
		//if (k<10)
		//	cout << "Aire pour courbure moyenne : " << area << endl;
		//k++;

		weight(v_it) = (fabs(area)>FLT_MIN ? 1.0 / (2.0 * area) : 0.0);
	  }
}


//-----------------------------------------------------------------------------


void Viewer::calc_mean_curvature() {

	  Mesh::VertexIter        v_it, v_end(mesh.vertices_end());
	  Mesh::HalfedgeHandle    h;
	  Mesh::EdgeHandle        e;
	  Mesh::VertexVertexIter  vv_it;
	  Mesh::Point             laplace(0.0, 0.0, 0.0);
	  Mesh::EdgeIter          e_it, e_end(mesh.edges_end());
	
	  calc_mean_weights();

	  //float summ = 0.0f;

	  for (v_it=mesh.vertices_begin(); v_it!=v_end; ++v_it)
	  {

		//summ = 0.0f;
		
		curvature_me(v_it) = 0.0;
		laplace = Mesh::Point(0,0,0);
		
		if (!mesh.is_boundary(v_it.handle()))
		{
		  for (vv_it=mesh.vv_iter(v_it); vv_it; ++vv_it)
		  {
					h = vv_it.current_halfedge_handle();
					e = mesh.edge_handle(h);
				
					laplace += weight(e) * (mesh.point(vv_it) - mesh.point(v_it));

					//summ += weight(e);
		  }
		  laplace *= weight(v_it);
		
		  //summ*= weight(v_it);
		  curvature_me(v_it) = laplace.norm();
		}

		//cout << " summ : " << summ << endl;
	  }
}


void Viewer::calc_gauss_curvature() {
	std::vector<std::vector<float>> princ_curvatures = ctools.principal_curvatures();
	int i = 0;
	for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it, ++i) {
		float k1 = princ_curvatures[i][0];
		float k2 = princ_curvatures[i][1];
		gauss_curvature(v_it) = k1 * k2;
	}
}


//-----------------------------------------------------------------------------


void Viewer::color_coding(OpenMesh::VPropHandleT<Mesh::Scalar> _curv) {
	

	  Mesh::VertexIter  v_it, v_end(mesh.vertices_end());
	  Mesh::Scalar      curv, min_curv(FLT_MAX), max_curv(-FLT_MAX);
	  Mesh::Color       col;
	
	  // put all curvature values into one array
	  std::vector<Mesh::Scalar> curv_values;
	  curv_values.reserve(mesh.n_vertices());
	  for (v_it=mesh.vertices_begin(); v_it!=v_end; ++v_it) {
		  curv_values.push_back(
			  mesh.property(_curv, v_it)
				//curvature(v_it)
			);
	  }
	
	  // discard upper and lower 5%
	  unsigned int n = curv_values.size()-1;
	  unsigned int i = n /20;
	  std::sort(curv_values.begin(), curv_values.end());
	  min_curv = curv_values[i];
	  max_curv = curv_values[n-1-i];

	  // define uniform color intervalls [v0,v1,v2,v3,v4]
	  Mesh::Scalar v0, v1, v2, v3, v4;
	  v0 = min_curv + 0.0/4.0 * (max_curv - min_curv);
	  v1 = min_curv + 1.0/4.0 * (max_curv - min_curv);
	  v2 = min_curv + 2.0/4.0 * (max_curv - min_curv);
	  v3 = min_curv + 3.0/4.0 * (max_curv - min_curv);
	  v4 = min_curv + 4.0/4.0 * (max_curv - min_curv);
	
	
	
	  // map curvatures to colors
	  for (v_it=mesh.vertices_begin(); v_it!=v_end; ++v_it) {

			//curv = curvature(v_it);
			curv = mesh.property(_curv, v_it);

			col = Mesh::Color(255,255,255);
    
			unsigned char u;
		
			if (curv < v0)
			{
			  col = Mesh::Color(0, 0, 255);
			}
			else if (curv > v4) 
			{
			  col = Mesh::Color(255, 0, 0);
			}
		
			else if (curv <= v2) 
			{
			  if (curv <= v1) // [v0, v1]
			  {
				u = (unsigned char) (255.0 * (curv - v0) / (v1 - v0));
				col = Mesh::Color(0, u, 255);
			  }      
			  else // ]v1, v2]
			  {
				u = (unsigned char) (255.0 * (curv - v1) / (v2 - v1));
				col = Mesh::Color(0, 255, 255-u);
			  }
			}
			else 
			{
			  if (curv <= v3) // ]v2, v3]
			  {
				u = (unsigned char) (255.0 * (curv - v2) / (v3 - v2));
				col = Mesh::Color(u, 255, 0);
			  }
			  else // ]v3, v4]
			  {
				u = (unsigned char) (255.0 * (curv - v3) / (v4 - v3));
				col = Mesh::Color(255, 255-u, 0);
			  }
			}
		
			mesh.set_color(v_it, col);
	  }
}


//-----------------------------------------------------------------------------


void Viewer::draw(const std::string& _draw_mode) {


 if (_draw_mode == "Mean Curvature") {

		color_coding(v_mean_curvature_);
		
		glDisable(GL_LIGHTING);
	  //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	
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
		prev_id_draw_mode=get_draw_mode();
  }

	else if (_draw_mode == "Gaussian Curvature") {

		color_coding(v_gauss_curvature_);

		glDisable(GL_LIGHTING);
		//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);


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


  else if (_draw_mode == "Reflection Lines") {

		glColor3f(0.3, 0.3, 0.3);
		glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT);
		glEnable(GL_COLOR_MATERIAL);

		glTexGeni( GL_S, GL_TEXTURE_GEN_MODE, GL_SPHERE_MAP );
		glTexGeni( GL_T, GL_TEXTURE_GEN_MODE, GL_SPHERE_MAP );
		glEnable( GL_TEXTURE_GEN_S );
		glEnable( GL_TEXTURE_GEN_T );
		glEnable( GL_TEXTURE_2D );    
		glEnable(GL_LIGHTING);
		glShadeModel(GL_SMOOTH);

		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_NORMAL_ARRAY);

		//Base surface
		GL::glVertexPointer(mesh.points());
		GL::glNormalPointer(mesh.vertex_normals());
		glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, &indices[0]);


		glDisableClientState(GL_VERTEX_ARRAY);
		glDisableClientState(GL_NORMAL_ARRAY);
		
		glDisable( GL_TEXTURE_GEN_S );
		glDisable( GL_TEXTURE_GEN_T );
		glDisable( GL_TEXTURE_2D );

		prev_draw_mode = current_draw_mode;
		prev_id_draw_mode=get_draw_mode();
  }
	
	else if (_draw_mode == "Vertex 1-Ring")
	{
		glDisable(GL_LIGHTING);

		VertexHandle _vh = mesh.vertex_handle(0); //TODO
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glEnable(GL_POLYGON_OFFSET_LINE);
		glPolygonOffset(1.0, -1.0);

		glBegin(GL_LINE_LOOP);
		glColor3f(1.0, 0.0, 0.0);
		for (auto vv_it = mesh.vv_iter(_vh); vv_it; ++vv_it) {
			GL::glVertex(mesh.point(vv_it));
		}
		glEnd();
		glDisable(GL_POLYGON_OFFSET_LINE);

		glColor3f(.2, .2, .2);
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
		prev_id_draw_mode = get_draw_mode();
	}
	
  else 
	  MeshViewer::draw(_draw_mode);
}

// Colors faces of the 1-ring of a given vertex
void Viewer::color_1_ring(Mesh::VertexHandle _vh) {
	Mesh::Color neighbor_color = Mesh::Color(255, 0, 0);
	for (auto vv_it = mesh.vv_iter(_vh); vv_it; ++vv_it) {
		mesh.set_color(vv_it, neighbor_color);
	}
}


//=============================================================================
