

#include <iostream>
#include "ComputingTools.hpp"

#include <minmax.h>

//#define mini(a,b)  ((a < b) ? a : b)
//#define maxi(a,b)  ((a > b) ? a : b)



ComputingTools::ComputingTools(void)
{
	volume_ = -1;
	area_ = -1;

	principal_curvatures_.clear();
}


ComputingTools::ComputingTools(Mesh _mesh)
{
	principal_curvatures_.clear();
}


ComputingTools::~ComputingTools(void)
{
}


double ComputingTools::volume()
{
	return volume_;
}


void ComputingTools::compute_volume()
{
	volume_ = 0.0f;
	Mesh::Point pointA, pointB, pointC;
	Mesh::ConstFaceVertexIter cfvIt;

	if(!mesh_.faces_empty()){

		//Iterator over all faces
		for (Mesh::FaceIter f_it=mesh_.faces_begin(); f_it!=mesh_.faces_end(); ++f_it){
			
			//Circulate around the vertex belonging to the current face
			cfvIt = mesh_.cfv_iter(f_it);

			pointA = mesh_.point(cfvIt.handle());
			pointB = mesh_.point((++cfvIt).handle());
			pointC = mesh_.point((++cfvIt).handle());

			volume_ += ((pointA+pointB+pointC)/3.0f)|((pointB-pointA)%(pointC-pointA));
			
		}
		volume_ /= 6.0f;
	}
	else{
		volume_ = -1;
		std::cerr<<"Warning: the mesh is empty !!"<<std::endl;
	}

	
	std::cout<<"Computation of the volume done: "<<volume_<<std::endl;
}


double ComputingTools::area()
{
	return area_;
}


void ComputingTools::compute_area()
{
	area_ = 0.0f;
	Mesh::Point pointA, pointB, pointC;
	Mesh::ConstFaceVertexIter cfvIt;

	if(!mesh_.faces_empty()){

		//Iterator over all faces
		for (Mesh::FaceIter f_it=mesh_.faces_begin(); f_it!=mesh_.faces_end(); ++f_it){
			
			//Circulate around the vertex belonging to the current face
			cfvIt = mesh_.cfv_iter(f_it);

			pointA = mesh_.point(cfvIt.handle());
			pointB = mesh_.point((++cfvIt).handle());
			pointC = mesh_.point((++cfvIt).handle());
			area_ += ((pointB-pointA)%(pointC-pointA)).norm();		
		}
		area_ *= 0.5f;
	}
	else{
		area_=-1;
		std::cerr<<"Warning: the mesh is empty !!"<<std::endl;
	}

	std::cout<<"Computation of the area done: "<<area_<<std::endl;
}


OpenMesh::Vec3f ComputingTools::halfedgeVector(Mesh::HalfedgeHandle const &h) const {
	return mesh_.point(mesh_.to_vertex_handle(h)) - 
		mesh_.point(mesh_.from_vertex_handle(h));
}


void ComputingTools::calc_princ_curvatures(){

	float k1,k2;

	for(Mesh::ConstVertexIter i = mesh_.vertices_begin(), ie = mesh_.vertices_end(); i != ie; ++i) {
	
		compute_local_principalCurvatures(k1,k2,i.handle());
		
		float max_k = max(k1,k2);
		float min_k = min(k1,k2);

		std::vector<float> main_curv; 
		main_curv.push_back(min_k);
		main_curv.push_back(max_k);
		principal_curvatures_.push_back(main_curv);

		//std::cout<<"(min,max) = ("<<min_k<<","<<max_k<<")"<<std::endl;	
	}

}


void ComputingTools::compute_local_principalCurvatures(float &_k1, float &_k2, const Mesh::VertexHandle _vh){

	float H_loc = 0.0f;
	float K_loc = 0.0f;
	float area_i=0.0f;

	//Pour vérification :
	//Théoriquement pour notre sphère avec normales à l'intérieur :
	// k1, k2 = 0.70915 : ok
	// K = 0.50289 : ok
	// H = 0.70915 : ok
	//Avec les normales à l'extérieur, les courbures principales devraient être négatives => delta_max et rho_max seraient infinies
	//Rayon théorique : 1.414214

	const Mesh::Point pointA = mesh_.point(_vh);
	Mesh::Point pointB, pointC;	


	//Compute the local gaussian curvature ************************

	//"Classical" formulation

	/*K_loc = float(TWOPI);
	if (mesh_.is_boundary(_vh)) K_loc = float(PI);

	Mesh::VertexOHalfedgeIter voh_it = mesh_.voh_iter(_vh);
	
	//Outgoing halfedge iterator 
	for(Mesh::VertexOHalfedgeIter voh_it = mesh_.voh_iter(_vh); voh_it; ++voh_it){
		
		//For the local Area
		pointB = mesh_.point(mesh_.to_vertex_handle(voh_it)); //Handle of the vertex pointed by the current outgoing halfedge
		Mesh::HalfedgeHandle opposite_heh = mesh_.opposite_halfedge_handle(voh_it);//Switch to opposite halfedge
		Mesh::HalfedgeHandle next_heh = mesh_.next_halfedge_handle(voh_it);//Next halfedge points to neighbouring vertex
		pointC = mesh_.point(mesh_.to_vertex_handle(next_heh));			
		area_i += ((pointB-pointA)%(pointC-pointA)).norm();
		

	    K_loc -= acos(OpenMesh::sane_aarg( ((pointB-pointA).normalize() | (pointC-pointA).normalize()) ));
	}
	area_i*=0.5f;
	K_loc *=3.0f/area_i;*/

	//Formulation taking into account the lengths of the opposite edges (Boix)

	float denom;
	float opp_lengths = 0.0f;
	float alpha_i;
	float l_i_sqr;

	K_loc = float(TWOPI);
	if (mesh_.is_boundary(_vh)) K_loc = float(PI);

	Mesh::VertexOHalfedgeIter voh_it = mesh_.voh_iter(_vh);

	//Outgoing halfedge iterator
	for(Mesh::VertexOHalfedgeIter voh_it = mesh_.voh_iter(_vh); voh_it; ++voh_it){

	//For the local Area
	pointB = mesh_.point(mesh_.to_vertex_handle(voh_it)); //Handle of the vertex pointed by the current outgoing halfedge
	Mesh::HalfedgeHandle opposite_heh = mesh_.opposite_halfedge_handle(voh_it);//Switch to opposite halfedge
	Mesh::HalfedgeHandle next_heh = mesh_.next_halfedge_handle(voh_it);//Next halfedge points to neighbouring vertex
	pointC = mesh_.point(mesh_.to_vertex_handle(next_heh));
	area_i += cross(pointB-pointA,pointC-pointA).norm();

	alpha_i = acos(OpenMesh::sane_aarg( ((pointB-pointA).normalize() | (pointC-pointA).normalize()) ));
	l_i_sqr = (pointC-pointB).sqrnorm();

	opp_lengths -= cot(alpha_i)*l_i_sqr;
	K_loc -= alpha_i;
	}

	area_i*=0.25f;
	opp_lengths *= 0.125f;
	denom = area_i + opp_lengths;
	K_loc /= denom;


	//Compute the local mean curvature ********************************
	
	//"Classical" formulation
	for(Mesh::ConstVertexEdgeIter j = mesh_.ve_iter(_vh); j; ++j) {
		
		double angle;
		Mesh::HalfedgeHandle h1 = mesh_.halfedge_handle(j.handle(), 0);
		Mesh::HalfedgeHandle h2 = mesh_.halfedge_handle(j.handle(), 1);
		OpenMesh::Vec3f v = halfedgeVector(h1);

		if(mesh_.is_boundary(h1) || mesh_.is_boundary(h2))
			angle = 0.0;
		else {
			OpenMesh::Vec3f n1,n2;
			n1 = mesh_.normal(mesh_.face_handle(h1));
			n2 = mesh_.normal(mesh_.face_handle(h2));
		
			angle = acos(min(max(n1 | n2, -1.0f), 1.0f));
			angle *= ((n1 % n2) | v) >= 0.0 ? 1.0 : -1.0;
		}
		H_loc += float(angle) * v.norm();
	}

	H_loc *=0.25f*3.0f/area_i;
	
	std::cout << H_loc << " ";
	
	H_loc = 0;

	//Formulation taking into account the lengths of the opposite ed)ges (Boix)
	for(Mesh::ConstVertexEdgeIter j = mesh_.ve_iter(_vh); j; ++j) {

		double angle;
		Mesh::HalfedgeHandle h1 = mesh_.halfedge_handle(j.handle(), 0);
		Mesh::HalfedgeHandle h2 = mesh_.halfedge_handle(j.handle(), 1);
		OpenMesh::Vec3f v = halfedgeVector(h1);

		if(mesh_.is_boundary(h1) || mesh_.is_boundary(h2))
		angle = 0.0;
		else {
			OpenMesh::Vec3f n1,n2;

			n2 = mesh_.normal(mesh_.face_handle(h1));
			n1 = mesh_.normal(mesh_.face_handle(h2));
		
			angle = acos(min(max(n1 | n2, -1.0f), 1.0f));

			angle *= ((n1 % n2) | v) >= 0.0 ? 1.0 : -1.0;
		}
		H_loc += float(angle) * v.norm();
	}

	H_loc *= 0.25f;
	H_loc /= denom;

	std::cout << H_loc <<std::endl;




	//Compute the principal curvatures
	_k1 = H_loc;// +sqrt(std::abs(H_loc*H_loc - K_loc) < EPS ? 0.0f : std::abs(H_loc*H_loc - K_loc));
	_k2 = K_loc;// -sqrt(std::abs(H_loc*H_loc - K_loc) < EPS ? 0.0f : std::abs(H_loc*H_loc - K_loc));

	//if(_k1>100 || _k2>100)
	//{
	//	std::cout<<"Les courbures principales sont (k1,k2) = ("<<_k1<<","<<_k2<<")"<<std::endl;
	//}

	//std::cout<<H_loc*H_loc-K_loc<<std::endl;
	//std::cout<<"Les courbures principales sont (k1,k2) = ("<<_k1<<","<<_k2<<")"<<std::endl;
}	


void ComputingTools::set_mesh(Mesh _mesh)
{
	mesh_ = _mesh;
}

std::vector<std::vector<float>> ComputingTools::principal_curvatures(){

	return principal_curvatures_;
}

float ComputingTools::cot(float _t) {

	return std::cos(_t) / std::sin(_t);
}