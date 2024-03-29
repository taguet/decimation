
//=============================================================================
//
//  CLASS Viewer - IMPLEMENTATION
//
//=============================================================================


//== INCLUDES =================================================================

#include "Viewer.hpp"
#include <vector>
#include <set>
#include <map>
#include <float.h>


//== IMPLEMENTATION ========================================================== 

Viewer::Viewer(const char* _title, int _width, int _height): MeshViewer(_title, _width, _height) {

  mesh.request_vertex_colors();
	
  mesh.add_property(initial_coords);

  mesh.add_property(v_mean_curvature_);
  mesh.add_property(v_gauss_curvature_);

  mesh.add_property(vweight_);
  mesh.add_property(eweight_);

  add_draw_mode("Mean Curvature");
  add_draw_mode("Gaussian Curvature");
  add_draw_mode("Reflection Lines");
  add_draw_mode("Vertex 1-Ring");
  add_draw_mode("Uniform Laplacian");
  add_draw_mode("Cotangent Laplacian");
  add_draw_mode("Anisotropic Laplacian");
  add_draw_mode("Planar Region Extraction");
  add_draw_mode("Mesh Simplification");
  add_draw_mode("Erreur plan");

  add_draw_mode("Debug Mesh Simplification");
  add_draw_mode("Debug Planar Region Extraction");
  add_draw_mode("Debug plane-region fitting");
  add_draw_mode("Debug contour");
  add_draw_mode("Debug contour projection");
  add_draw_mode("Debug opposite angles");
  add_draw_mode("Debug laplacien cotangente");
  add_draw_mode("Debug lissage");
  add_draw_mode("Debug centro�de");
  add_draw_mode("Debug normal angles");
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
  mtools = MeshTools();
	
	
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
		  isModified = false;
		  v_id = 0;
		  neighbour_offset = 0;
		ctools.set_mesh(mesh);
		mtools.setMesh(mesh);
		ctools.calc_princ_curvatures();
		
		//Store initial coordinates for all vertices
		store_initial_points();
		
		//Compute all the information about curvature
		calc_mean_curvature();
		calc_gauss_curvature();

		glutPostRedisplay();
		return true;
	  }
	  return false;
}

bool Viewer::write_mesh(const char* _filename) {
	if (MeshViewer::write_mesh(_filename)) {
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------

void Viewer::store_initial_points() {
	for (auto v_it{ mesh.vertices_begin() }; v_it != mesh.vertices_end(); ++v_it) {
		initial_coord(v_it) = mesh.point(v_it);
	}
}


void Viewer::reset_mesh() {
	for (auto v_it{ mesh.vertices_begin() }; v_it != mesh.vertices_end(); ++v_it) {
		mesh.set_point(v_it, initial_coord(v_it));
	}
	isModified = false;
	mesh.update_normals();
	set_draw_mode(default_id_draw_mode);
}


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
	current_draw_mode = _draw_mode;

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

		VertexHandle _vh = mesh.vertex_handle(v_id);
		draw_1_ring(_vh, { 1.0, 0.0, 0.0 });

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

	else if (_draw_mode == "Uniform Laplacian") {
		if (!isModified) {
			//mtools.taubinSmoothing<UniformLaplacian>(20, 1, -.5);
			mtools.smoothMesh<UniformLaplacian>(20, 1);
			isModified = true;
		}
		prev_draw_mode = current_draw_mode;
		prev_id_draw_mode = get_draw_mode();
		draw("Solid Smooth");
	}
	else if (_draw_mode == "Cotangent Laplacian") {
		if (!isModified) {
			mtools.taubinSmoothing<CotangentLaplacian>(50, .005, -.005);
			//mtools.smoothMesh<CotangentLaplacian>(20, .005);
			isModified = true;
		}
		prev_draw_mode = current_draw_mode;
		prev_id_draw_mode = get_draw_mode();
		draw("Solid Smooth");
	}
	else if (_draw_mode == "Anisotropic Laplacian") {
		if (!isModified) {
			auto start{ std::chrono::steady_clock::now() };
			//mtools.taubinSmoothing<UniformLaplacian>(20, 1, -.5);
			mtools.smoothMesh<AnisotropicLaplacian>(20, 1);
			isModified = true; 
			auto end{ std::chrono::steady_clock::now() };
			auto duration{ std::chrono::duration_cast<std::chrono::milliseconds>(end - start) };
			std::cerr << "Regions extracted in " << duration.count() << " milliseconds.\n\n";
		}
		prev_draw_mode = current_draw_mode;
		prev_id_draw_mode = get_draw_mode();
		draw("Solid Smooth");
	}
	else if (_draw_mode == "Planar Region Extraction") {
		if (!isModified) {
			auto start{ std::chrono::steady_clock::now() };
			isModified = true;
			this->graph = new TopologyGraph(mesh, 1.0f, 1.0f);
			mtools.extractRegions(*graph);
			this->graph->projectContourVertices();
			auto end{ std::chrono::steady_clock::now() };
			auto duration{ std::chrono::duration_cast<std::chrono::milliseconds>(end - start) };

			std::cerr << "Executed in " << duration.count() << " milliseconds.\n\n";
		}
		draw("Solid Smooth"); 
		prev_draw_mode = current_draw_mode;
		prev_id_draw_mode = get_draw_mode();
	}
	else if (_draw_mode == "Mesh Simplification") {
		if (!isModified) {
			auto start{ std::chrono::steady_clock::now() };
			isModified = true;
			this->graph = new TopologyGraph(mesh, 1.0f, 1.0f);
			mtools.extractRegions(*graph);
			this->graph->projectContourVertices();
			mtools.simplifyMesh(*graph, graph->size()*2);
			update_face_indices(mesh, indices);
			auto end{ std::chrono::steady_clock::now() };
			auto duration{ std::chrono::duration_cast<std::chrono::milliseconds>(end - start) };

			std::cerr << "Executed in " << duration.count() << " milliseconds.\n\n";
		}
		draw("Solid Smooth");
		prev_draw_mode = current_draw_mode;
		prev_id_draw_mode = get_draw_mode();
	}
	else if (_draw_mode == "Debug Mesh Simplification") {
		if (!isModified) {
			isModified = true;
			this->graph = new TopologyGraph(mesh, 1.0f, 1.0f);
			mtools.extractRegions(*graph);
			this->graph->projectContourVertices();
			mesh.request_vertex_status();
			mesh.request_edge_status();
			mesh.request_face_status();
			ec = new EdgeCollapse{ mesh, *graph };
		}
		if (calledCollapse && ec != nullptr) {
			ec->collapse();
			//mesh.garbage_collection();
			update_face_indices(mesh, indices);
		}
		draw("Hidden Line"); 
		prev_draw_mode = current_draw_mode;
		prev_id_draw_mode = get_draw_mode();
	}
	else if (_draw_mode == "Debug Planar Region Extraction") {
		static std::map<int, std::unique_ptr<int[]>> rgb{};
		if (!isModified) {
			isModified = true;
			rgb.clear();
			this->graph = new TopologyGraph{ mesh, 1.0f, 1.0f };

			auto start{ std::chrono::steady_clock::now() };
			mtools.extractRegions(*graph);
			auto end{ std::chrono::steady_clock::now() };
			auto duration{ std::chrono::duration_cast<std::chrono::milliseconds>(end - start) }; 
			std::cerr << "Regions extracted in " << duration.count() << " milliseconds.\n\n";
			srand(0);
			std::set<int> groups_found{graph->getRegionIDs()};
			for (int groupID : groups_found) {
				//Generate color
				std::unique_ptr<int[]> color{ new int[3] {rand() % 256, rand() % 256, rand() % 256} };
				rgb[groupID] = std::move(color);
				Equation::Plane plane{ graph->getRegion(groupID).plane};
				//std::cout << "z = " << params[1] << "x + " << params[2] << "y + " << params[0] << std::endl;
			}
			const Mesh::Color* c{ mesh.vertex_colors() };
			std::unique_ptr<int[]> color{ new int[3] { c->data()[0], c->data()[1], c->data()[2]} };
			rgb[-1] = std::move(color);
			std::cerr << "Found " << groups_found.size() << " out of " << graph->size() << " groups.\n";
			for (int i : groups_found) {
				//for (auto face : graph->getRegion(i).getFaceHandles())
				//std::cerr << i << ": " << face.idx() << "\n";
			}
		}
		glDisable(GL_LIGHTING);
		for (auto f_iter{ mesh.faces_begin() }; f_iter != mesh.faces_end(); ++f_iter) {
			const int id{ graph->faceGroup(f_iter)};
			glColor3ub(rgb.at(id)[0], rgb.at(id)[1], rgb.at(id)[2]);
			glBegin(GL_TRIANGLES);
			Mesh::HalfedgeHandle heh{ mesh.halfedge_handle(f_iter) };
			GL::glVertex(mesh.point(mesh.from_vertex_handle(heh)));
			GL::glVertex(mesh.point(mesh.to_vertex_handle(heh)));
			GL::glVertex(mesh.point(mesh.to_vertex_handle(mesh.next_halfedge_handle(heh))));
			glEnd();
		}

		prev_draw_mode = current_draw_mode;
		prev_id_draw_mode = get_draw_mode();
	}
	else if (_draw_mode == "Debug plane-region fitting") {
		if (!isModified) {
			isModified = true;
			this->graph = new TopologyGraph(mesh, 1.0f, 1.0f);
			mtools.extractRegions(*graph);
		}

		//draw("Solid Smooth");
		set<int> ids{ graph->getRegionIDs() };
		auto it{ ids.begin() }; 
		std::advance(it, region_id% graph->size());
		int regionID{ *it };
		glEnable(GL_LIGHTING);
		for (auto fh : graph->getRegion(regionID).getFaceHandles()) {
			glColor3f(0.0f, 0.0f, 1.0f);
			glBegin(GL_TRIANGLES);
			Mesh::HalfedgeHandle heh{ mesh.halfedge_handle(fh) };
			GL::glNormal(mesh.normal(fh));
			GL::glVertex(mesh.point(mesh.from_vertex_handle(heh)));
			GL::glVertex(mesh.point(mesh.to_vertex_handle(heh)));
			GL::glVertex(mesh.point(mesh.opposite_vh(heh)));
			glEnd();
		}

		const Equation::Plane& plane{ graph->getRegion(regionID).plane };
		std::cerr << "Viewing region " << regionID << "\nEquation: "<< plane.a() << "x + " << plane.b() << "y + " << plane.c() << " z + " << plane.d() << " = 0\n";
		glColor3f(1.0f, 0.0f, 0.0f);
		std::vector<Eigen::Vector3f> points{ computePlane(plane) };
		Vector3f& p0{ points[0]};
		Vector3f& p1{ points[1] };
		Vector3f& p2{ points[2] };
		std::cerr << "P0(" << p0 << ")\tnP1(" << p1 << "\nP2" << p2 << "\n";
		glBegin(GL_TRIANGLES);
			GL::glVertex(OpenMesh::Vec3f{ p0[0], p0[1], p0[2] });
			GL::glVertex(OpenMesh::Vec3f{ p1[0], p1[1], p1[2] });
			GL::glVertex(OpenMesh::Vec3f{ p2[0], p2[1], p2[2] });
		glEnd();
		VertexHandle vh{ *graph->getRegion(regionID).getVertexHandles().begin()};
		draw_vector(plane.getNormal(), MeshUtils::toEigen(mesh.point(vh)));

		prev_draw_mode = current_draw_mode;
		prev_id_draw_mode = get_draw_mode();
	}
	else if (_draw_mode == "Debug contour") {
		if (!isModified) {
			contour_edges = {};
			contour_vertices = {};
			lines = {};
			isModified = true;
			this->graph = new TopologyGraph(mesh, 1.0f, 1.0f);
			mtools.extractRegions(*graph);
			contour_edges = graph->extractContour();
			for (auto edge : contour_edges) {
				Mesh::HalfedgeHandle heh{ mesh.halfedge_handle(edge, 0) };
				contour_vertices.insert(mesh.from_vertex_handle(heh));
				contour_vertices.insert(mesh.to_vertex_handle(heh));
			}
			lines = this->graph->findPlanePlaneIntersections();
		}
		glEnable(GL_LIGHTING);

		draw("Solid Smooth");

		glEnable(GL_DEPTH_TEST);
		glEnable(GL_POLYGON_OFFSET_LINE);
		glPolygonOffset(-1.0, 1.0);
		glLineWidth(3.0f);
		glBegin(GL_LINES);
			glColor3f(0.0f, 0.0f, 1.0f);	// draw contour on mesh
			for (auto edge : contour_edges) {
				Mesh::HalfedgeHandle heh{ mesh.halfedge_handle(edge, 0) };
				GL::glVertex(mesh.point(mesh.from_vertex_handle(heh)));
				GL::glVertex(mesh.point(mesh.to_vertex_handle(heh)));
			}
			glColor3f(1.0f, 0.0f, 0.0f);	// draw plane plane intersections
			for (auto& line : lines) {
				Vector3f p0{ line.evaluate(-2.0f) };
				auto t = p0[0];
				Vector3f p1{ line.evaluate(2.0f) };
				std::cerr << "Line: (" << line.origin.transpose() << ") + t*(" << line.direction.transpose() << ")\n";
				std::cerr << "p0=" << p0.transpose() << "\tp1=" << p1.transpose() << '\n';
				GL::glVertex(OpenMesh::Vec3f{ p0(0), p0(1), p0(2) });
				GL::glVertex(OpenMesh::Vec3f{ p1(0), p1(1), p1(2) });
			}
		glEnd();
		glLineWidth(1.0f);
		glDisable(GL_POLYGON_OFFSET_LINE);

		prev_draw_mode = current_draw_mode;
		prev_id_draw_mode = get_draw_mode();
	}
	else if (_draw_mode == "Debug contour projection") {
		if (!isModified) {
			contour_edges = {};
			contour_vertices = {};
			lines = {};
			isModified = true;
			this->graph = new TopologyGraph(mesh, 1.0f, 1.0f);
			mtools.extractRegions(*graph);
			auto start{ std::chrono::steady_clock::now() };
			contour_edges = graph->extractContour();
			for (auto edge : contour_edges) {
				Mesh::HalfedgeHandle heh{ mesh.halfedge_handle(edge, 0) };
				contour_vertices.insert(mesh.from_vertex_handle(heh));
				contour_vertices.insert(mesh.to_vertex_handle(heh));
			}
			this->graph->projectContourVertices();
			auto end{ std::chrono::steady_clock::now() };
			auto duration{ std::chrono::duration_cast<std::chrono::milliseconds>(end - start) };

			std::cerr << "Executed in " << duration.count() << " milliseconds.\n\n";
		}
		glEnable(GL_LIGHTING);

		draw("Solid Smooth");

		glEnable(GL_DEPTH_TEST);
		glEnable(GL_POLYGON_OFFSET_LINE);
		glPolygonOffset(-1.0, 1.0);
		glLineWidth(3.0f);
		glBegin(GL_LINES);
		glColor3f(0.0f, 0.0f, 1.0f);	// draw contour on mesh
		for (auto edge : contour_edges) {
			Mesh::HalfedgeHandle heh{ mesh.halfedge_handle(edge, 0) };
			GL::glVertex(mesh.point(mesh.from_vertex_handle(heh)));
			GL::glVertex(mesh.point(mesh.to_vertex_handle(heh)));
		}
		glColor3f(1.0f, 0.0f, 0.0f);	// draw plane plane intersections
		for (auto& line : lines) {
			Vector3f p0{ line.evaluate(-2.0f) };
			auto t = p0[0];
			Vector3f p1{ line.evaluate(2.0f) };
			std::cerr << "Line: (" << line.origin.transpose() << ") + t*(" << line.direction.transpose() << ")\n";
			std::cerr << "p0=" << p0.transpose() << "\tp1=" << p1.transpose() << '\n';
			GL::glVertex(OpenMesh::Vec3f{ p0(0), p0(1), p0(2) });
			GL::glVertex(OpenMesh::Vec3f{ p1(0), p1(1), p1(2) });
		}
		glEnd();
		glLineWidth(1.0f);
		glDisable(GL_POLYGON_OFFSET_LINE);

		prev_draw_mode = current_draw_mode;
		prev_id_draw_mode = get_draw_mode();
	}
	else if (_draw_mode == "Erreur plan") {
		if (!isModified) {
			contour_edges = {};
			contour_vertices = {};
			lines = {};
			isModified = true;
			this->graph = new TopologyGraph(mesh, 1.0f, 1.0f);
			mtools.extractRegions(*graph);
			computeFittingError();
		}
		draw("Solid Smooth");
		glEnable(GL_DEPTH_TEST);
		glEnable(GL_POLYGON_OFFSET_LINE);
		glPolygonOffset(-1.0, 1.0);
		glPointSize(10.0f);
		glDisable(GL_LIGHTING);
		glBegin(GL_POINTS);
		for (auto vh{ mesh.vertices_begin() }; vh != mesh.vertices_end(); ++vh) {
			const float error{ mesh.property(fit_error, vh) };
			if (error <= 0.1f)
				glColor3f(0.0f, 1.0f, 0.0f);
			else if (error > 0.1f && error <= 0.2f)
				glColor3f(0.0f, 0.0f, 1.0f);
			else if (error > 0.2f)
				glColor3f(1.0f, 0.0f, 0.0f);
			else
				continue;
			GL::glVertex(mesh.point(vh));
		}
		glEnd();
		glDisable(GL_POLYGON_OFFSET_LINE);		
		prev_draw_mode = current_draw_mode;
		prev_id_draw_mode = get_draw_mode();
	}
	else if (_draw_mode == "Debug opposite angles") {
		static int cur_v_id = -1;
		static std::vector<HalfedgeHandle> ohs{};
		if (cur_v_id != v_id) {
			ohs.clear();
			for (auto voh_it{ mesh.voh_iter(mesh.vertex_handle(v_id)) }; voh_it; ++voh_it)
				ohs.push_back(voh_it);
			cur_v_id = v_id;
		}
		static HalfedgeHandle cur_oh{};

		HalfedgeHandle oh{ ohs[neighbour_offset % ohs.size()] };
		glDisable(GL_LIGHTING);

		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glEnable(GL_POLYGON_OFFSET_LINE);
		glPolygonOffset(.5, 1.0);

		OpenMesh::HalfedgeHandle opposite{ mesh.opposite_halfedge_handle(oh) };
		OpenMesh::HalfedgeHandle left{ mesh.next_halfedge_handle(oh) };
		OpenMesh::HalfedgeHandle right{ mesh.next_halfedge_handle(opposite) };

		if (cur_oh != oh) {
			cur_oh = oh;
			std::cerr << "v_i=" << mesh.point(mesh.from_vertex_handle(oh)) << '\n'
					  << "left hh: " << left << '\n'
					  << "v_j=" << mesh.point(mesh.to_vertex_handle(oh)) << '\n'
					  << "right hh: " << right << '\n';
			Vec3f a;
			mesh.calc_edge_vector(left, a);
			std::cerr << "next()=" << a << '\n';
			std::pair<float, float> angles{ MeshUtils::getOppositeAngles(mesh, oh) };
			Vec3f e1, e2;
			mesh.calc_sector_vectors(left, e1, e2);
			std::cerr << "alpha=" << angles.first << "\te1=" << e1 << "\te2=" << e2 << '\n'
					  << "cotan alpha=" << MeshUtils::cotan(angles.first) << '\n';
			mesh.calc_sector_vectors(right, e1, e2);
			std::cerr << "beta=" << angles.second << "\te1=" << e1 << "\te2=" << e2 << '\n'
					  << "cotan beta=" << MeshUtils::cotan(angles.second) << '\n';
		}

		draw_1_ring(mesh.vertex_handle(v_id), {0.5, 0.5, 0.5});

		glPolygonOffset(1.0, -1.0);
		glBegin(GL_LINES);
			glColor3f(0.0, 0.0, 1.0);
			GL::glVertex(mesh.point(mesh.from_vertex_handle(oh)));
			GL::glVertex(mesh.point(mesh.to_vertex_handle(oh)));
		glEnd();

		glBegin(GL_LINE_STRIP);
			glColor3f(1.0, 0.0, 0.0);
			GL::glVertex(mesh.point(mesh.from_vertex_handle(left)));
			GL::glVertex(mesh.point(mesh.to_vertex_handle(left)));
			GL::glVertex(mesh.point(mesh.to_vertex_handle(mesh.next_halfedge_handle(left))));
		glEnd();

		glBegin(GL_LINE_STRIP);
			glColor3f(0.0, 1.0, 0.0);
			GL::glVertex(mesh.point(mesh.from_vertex_handle(right)));
			GL::glVertex(mesh.point(mesh.to_vertex_handle(right)));
			GL::glVertex(mesh.point(mesh.to_vertex_handle(mesh.next_halfedge_handle(right))));
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
	else if (_draw_mode == "Debug laplacien cotangente") {
		VertexHandle _vh = mesh.vertex_handle(v_id);
		static int cur_v_id = -1;
		static CotangentLaplacian debug_clapl{ CotangentLaplacian(mesh) };
		if (cur_v_id != v_id) {
			std::cerr << "v_i=" << mesh.point(_vh) << '\n'
					  << "L(v_i)=" << debug_clapl.computeLaplacian(_vh) << '\n'
					  << "A(v_i)=" << MeshUtils::computeVertexArea(mesh, _vh) << '\n';
			cur_v_id = v_id;
		}
		glDisable(GL_LIGHTING);
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
	else if (_draw_mode == "Debug lissage") {
		if (calledSmoothing) {
			mtools.smoothMesh<AnisotropicLaplacian>(1, 1.0);
		}
		glEnable(GL_LIGHTING);

		glShadeModel(GL_SMOOTH);
		glColor3f(0.3, 0.3, 0.3);
		glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT);
		glEnable(GL_COLOR_MATERIAL);

		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_NORMAL_ARRAY);
		GL::glVertexPointer(mesh.points());
		GL::glNormalPointer(mesh.vertex_normals());

		glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, &indices[0]);

		glDisableClientState(GL_VERTEX_ARRAY);
		glDisableClientState(GL_NORMAL_ARRAY);


		prev_draw_mode = current_draw_mode;
		prev_id_draw_mode = get_draw_mode();
	}
	else if (_draw_mode == "Debug centro�de") {
		VertexHandle _vh = mesh.vertex_handle(v_id);
		static int cur_v_id = -1;
		static std::vector<HalfedgeHandle> ohs{};
		if (cur_v_id != v_id) {
			cur_v_id = v_id;
			ohs.clear();
			for (auto voh_it{ mesh.voh_iter(mesh.vertex_handle(v_id)) }; voh_it; ++voh_it)
				ohs.push_back(voh_it);
		}
		static HalfedgeHandle cur_oh{};		
		HalfedgeHandle oh{ ohs[neighbour_offset % ohs.size()] };

		OpenMesh::VertexHandle v0, v1, v2;
		v0 = mesh.from_vertex_handle(oh);
		v1 = mesh.to_vertex_handle(oh);
		v2 = mesh.to_vertex_handle(mesh.next_halfedge_handle(oh));
		Mesh::Point p0, p1, p2;
		p0 = mesh.point(v0);
		p1 = mesh.point(v1);
		p2 = mesh.point(v2);

		if (cur_oh != oh) {
			cur_oh = oh;
			std::cerr << "v_i=" << mesh.point(mesh.from_vertex_handle(oh)) << '\n'
					  << "v_j=" << mesh.point(mesh.to_vertex_handle(oh)) << '\n'
					  << "v0=" << p0 << '\n'
					  << "v1=" << p1 << '\n'
					  << "v2=" << p2 << '\n'
					  << "centroid=" << MeshUtils::computeCentroid(mesh, oh) << '\n'
					  << (p0 + p1 + p2) / 3.0f << '\n';
		}

		glDisable(GL_LIGHTING);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glEnable(GL_POLYGON_OFFSET_LINE);
		glPolygonOffset(1.0, -1.0);

		glBegin(GL_LINE_LOOP);
		glColor3f(1.0, 0.0, 0.0);
			GL::glVertex(p0);
		glColor3f(0.0, 1.0, 0.0);
			GL::glVertex(p1);
		glColor3f(0.0, 0.0, 1.0);
			GL::glVertex(p2);
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
	else if (_draw_mode == "Debug normal angles") {
		static int cur_e_id = -1;
		if (cur_e_id != e_id) {
			cur_e_id = e_id;
			EdgeHandle eh{ mesh.edge_handle(e_id) };
			HalfedgeHandle hh0{ mesh.halfedge_handle(eh, 0) };
			HalfedgeHandle hh1{ mesh.halfedge_handle(eh, 1) };
			if (hh0 != Mesh::InvalidHalfedgeHandle && hh1 != Mesh::InvalidHalfedgeHandle) {
				FaceHandle fh0{ mesh.face_handle(hh0) };
				FaceHandle fh1{ mesh.face_handle(hh1) };
				if (fh0.is_valid() && fh1.is_valid()) {
					Mesh::Normal n0{ mesh.normal(fh0) };
					Mesh::Normal n1{ mesh.normal(fh1) };
					float cos{ dot(n0, n1) };
					float sin{ cross(n0, n1).length() };
					float angle{ acosf(cos) };
					std::cerr << "Normal 1 (" << n0 << ")\tNormal 2(" << n1 << ")\n";
					std::cerr << "cos = " << cos << "\tsin = " << sin << "\nAngle = " << angle << "(" << rad_to_deg(angle) << ")\n";
					std::cerr << "Normals are close ? (> cos(" << deg_to_rad(20.0f) << ")=" << cosf(deg_to_rad(20.0f)) << "\t" << (dot(n0, n1) > cosf(deg_to_rad(20.0f))) << "\n\n";
				}
			}
		}
		EdgeHandle eh{ mesh.edge_handle(e_id) };
		Mesh::VertexHandle vh_0, vh_1;
		MeshUtils::vertex_handles(mesh, eh, vh_0, vh_1); 
		glDisable(GL_LIGHTING);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glEnable(GL_POLYGON_OFFSET_LINE);

		glPolygonOffset(1.0, -1.0);
		glBegin(GL_LINES);
			glColor3f(0.0, 0.0, 1.0);
			GL::glVertex(mesh.point(vh_0));
			GL::glVertex(mesh.point(vh_1));
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

		HalfedgeHandle hh{ MeshUtils::findHalfedge(mesh, vh_0, vh_1) };
		draw_normal(mesh.face_handle(hh), 1.0f);
		draw_normal(mesh.face_handle(mesh.opposite_halfedge_handle(hh)), 1.0f);

		prev_draw_mode = current_draw_mode;
		prev_id_draw_mode = get_draw_mode();
	}

  else 
	  MeshViewer::draw(_draw_mode);
	calledSmoothing = false;
	calledCollapse = false;
}

void Viewer::draw_1_ring(const VertexHandle vh, Vec3f color) {
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glEnable(GL_POLYGON_OFFSET_LINE);
	glPolygonOffset(1.0, -1.0);

	glBegin(GL_LINE_LOOP);
	glColor3f(color[0], color[1], color[2]);
	for (auto vv_it = mesh.vv_iter(vh); vv_it; ++vv_it) {
		GL::glVertex(mesh.point(vv_it));
	}
	glEnd();
	glDisable(GL_POLYGON_OFFSET_LINE);
}

void Viewer::draw_normal(const FaceHandle fh, float length) {
	if (!fh.is_valid())	return;
	Vec3f centroid{ MeshUtils::computeCentroid(mesh, mesh.halfedge_handle(fh)) };
	Vec3f dest{ centroid + mesh.normal(fh) * length };

	glBegin(GL_LINES);
		glColor3f(1.0f, 0.0f, 0.0f);
		glLineWidth(1.0f);
		GL::glVertex(centroid);
		GL::glVertex(dest);
	glEnd();
}


void Viewer::draw_vector(const Vector3f vec, const Vector3f position) {
	Vector3f dest{ position + vec };
	glBegin(GL_LINES);
		glColor3f(1.0f, 0.0f, 0.0f);
		glLineWidth(1.0f);
		GL::glVertex(Vec3f(position[0], position[1], position[2]));
		GL::glVertex(Vec3f(dest[0], dest[1], dest[2]));
	glEnd();
}


// Colors faces of the 1-ring of a given vertex
void Viewer::color_1_ring(Mesh::VertexHandle _vh) {
	Mesh::Color neighbor_color = Mesh::Color(255, 0, 0);
	for (auto vv_it = mesh.vv_iter(_vh); vv_it; ++vv_it) {
		mesh.set_color(vv_it, neighbor_color);
	}
}

// Opens file explorer to select a new mesh file
void Viewer::browse_meshes() {
	OPENFILENAME ofn;
	wchar_t szPath[100];
	ZeroMemory(&ofn, sizeof(ofn));

	ofn.lStructSize = sizeof(ofn);
	ofn.lpstrFilter = L"Mesh\0*.off;*.stl;*.obj\0\0";
	ofn.lpstrFile = szPath;
	ofn.lpstrFile[0] = '\0';
	ofn.nMaxFile = 100;
	ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;

	if (GetOpenFileName(&ofn)) {
		wstring ws(szPath);
		string str(ws.begin(), ws.end());
		open_mesh(str.c_str());
	}
}


void Viewer::save_mesh() {
	OPENFILENAME ofn;
	wchar_t szPath[100];
	ZeroMemory(&ofn, sizeof(ofn));

	ofn.lStructSize = sizeof(ofn);
	ofn.lpstrFilter = L"Mesh\0*.off;*.stl;*.obj\0\0";
	ofn.lpstrFile = szPath;
	ofn.lpstrFile[0] = '\0';
	ofn.lpstrDefExt = L"off";
	ofn.nMaxFile = 100;
	ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;

	if (GetSaveFileName(&ofn)) {
		wstring ws(szPath);
		string str(ws.begin(), ws.end());
		write_mesh(str.c_str());
	}
}


void Viewer::keyboard(int key, int x, int y) {
	switch (key)
	{
		case 'o':
		{
			browse_meshes();
			break;
		}
		case 'p':
		{
			++region_id;
			glutPostRedisplay();
			break;
		}
		case 32:
		{
			calledSmoothing = true;
			calledCollapse = true;
			break;
		}
		case 'r':
		{
			reset_mesh();
			break;
		}
		case 's':
		{
			save_mesh();
			break;
		}
		default:
		{
			MeshViewer::keyboard(key, x, y);
			break;
		}
	}
}


void Viewer::special(int key, int x, int y) {
	switch (key) {
		case GLUT_KEY_LEFT:
		{			
			if (current_draw_mode == "Debug normal angles") {
				e_id = (e_id - 1) % mesh.n_edges();
				std::cerr << "e_id=" << e_id << '\n';
			}
			else {
				v_id = (v_id - 1) % mesh.n_vertices();
				std::cerr << "v_id=" << v_id << '\n';
			}
			glutPostRedisplay();
			break;
		}
		case GLUT_KEY_RIGHT:
		{
			if (current_draw_mode == "Debug normal angles") {
				e_id = (e_id + 1) % mesh.n_edges();
				std::cerr << "e_id=" << e_id << '\n';
			}
			else {
				v_id = (v_id + 1) % mesh.n_vertices();
				std::cerr << "v_id=" << v_id << '\n';
			}
			glutPostRedisplay();
			break;
		}
		case GLUT_KEY_UP:
		{
			if (current_draw_mode == "Debug normal angles") {
				e_id = (e_id + 100) % mesh.n_edges();
				std::cerr << "e_id=" << e_id << '\n';
			}
			else {
				++neighbour_offset;
				std::cerr << "neighbour_offset=" << neighbour_offset << '\n';
			}
			glutPostRedisplay();
			break;
		}
		case GLUT_KEY_DOWN:
		{
			if (current_draw_mode == "Debug normal angles") {
				e_id = (e_id - 100) % mesh.n_edges();
				std::cerr << "e_id=" << e_id << '\n';
			}
			else {
				--neighbour_offset;
				std::cerr << "neighbour_offset=" << neighbour_offset << '\n';
			}
			glutPostRedisplay();
			break;
		}
		default:
		{
			MeshViewer::special(key, x, y);
			break;
		}
	}
}


/// @brief Compute 3 points of a plane
/// @param plane_params A vector (c, a, b)
/// @return A 3x3 matrix with each column representing a point.
std::vector<Eigen::Vector3f> Viewer::computePlane(const Equation::Plane& plane) const {
	//ax + by + c - z = 0
	//On trouve 3 points sur le plan
	float a{ plane.a()};
	float b{ plane.b()};
	float c{ plane.c() };
	float d{ plane.d() };
	std::vector<Eigen::Vector3f> points{ 3 };
	points.at(0) = Vector3f{0, 0, -d / c};	//x=0 and y=0
	points.at(1) = Vector3f{0, -d/b, 0};	//x=0 and z=0
	points.at(2) = Vector3f{-d/a, 0, 0};	//y=0 and z=0
	return points;
}


void Viewer::computeFittingError() {
	mesh.add_property(fit_error);
	for (auto vh{ mesh.vertices_begin() }; vh != mesh.vertices_end(); ++vh) {
		mesh.property(fit_error, vh) = -1;
	}
	for (auto fh{ mesh.faces_begin() }; fh != mesh.faces_end(); ++fh) {
		const int id{ this->graph->getFaceRegion(fh) };
		if (id == -1) continue;
		const Plane& plane{ this->graph->getPlane(id) };
		const HalfedgeHandle heh{ mesh.halfedge_handle(fh) };
		const VertexHandle vh_0{ mesh.to_vertex_handle(heh) };
		const VertexHandle vh_1{ mesh.from_vertex_handle(heh) };
		const VertexHandle vh_2{ mesh.opposite_vh(heh)};

		mesh.property(fit_error, vh_0) = plane.distToPoint(mesh.point(vh_0));
		mesh.property(fit_error, vh_1) = plane.distToPoint(mesh.point(vh_1));
		mesh.property(fit_error, vh_2) = plane.distToPoint(mesh.point(vh_2));
	}
}

//=============================================================================
