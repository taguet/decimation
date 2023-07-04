#include "MeshTools.h"
#include <list>

void MeshTools::extractRegions(TopologyGraph& graph) {
	//initialization
	std::list<OpenMesh::FaceHandle> ungrouped_faces{};
	for (auto f_iter{ mesh_->faces_begin() }; f_iter != mesh_->faces_end(); ++f_iter) {
		graph.faceGroup(f_iter) = -1;
		ungrouped_faces.push_back(f_iter);
	}
	std::cerr << "Executing region growing algorithm...\n";
	growRegions(ungrouped_faces, graph);
	std::cerr << "Done\nBuilding topology graph...\n";
	buildTopologyGraph(graph);
}


void MeshTools::growRegions(std::list<Mesh::FaceHandle>& ungrouped_faces, TopologyGraph& graph) {
	int current_region_id{-1};
	while (!ungrouped_faces.empty()) {
		Mesh::FaceHandle fh{ ungrouped_faces.front() };
		ungrouped_faces.pop_front();
		graph.faceGroup(fh) = ++current_region_id;
		graph.addFaceToRegion(current_region_id, fh);
		std::list<Mesh::FaceHandle>  neighbors;
		MeshUtils::getFaceNeighbors(*mesh_, fh, neighbors);
		std::cerr << "Region " << current_region_id << '\n';

		for (auto f_neighbor{ neighbors.begin() }; f_neighbor != neighbors.end(); ) {
			Mesh::Normal f_normal_neighbor{ mesh_->normal(*f_neighbor)};
			Mesh::Normal f_normal{ mesh_->normal(fh)};

			if (faceIsGrouped(*f_neighbor, graph)) {
				f_neighbor = neighbors.erase(f_neighbor);
				continue;
			}
			else if (normalsAreCloseEnough(f_normal_neighbor, f_normal, 0.349066f)) {
				graph.addFaceToRegion(current_region_id, *f_neighbor);
				extendNeighborhood(*f_neighbor, neighbors);
				ungrouped_faces.remove(*f_neighbor);
			}
			++f_neighbor;
		}
	}
}


void MeshTools::buildTopologyGraph(TopologyGraph& graph) {
	for (auto f_iter{ mesh_->faces_begin() }; f_iter != mesh_->faces_end(); ++f_iter) {
		for (auto neighbor{ mesh_->ff_iter(f_iter) }; neighbor; ++neighbor) {
			if (graph.faceGroup(f_iter) != graph.faceGroup(neighbor)) {
				graph.connectRegions(graph.faceGroup(f_iter), graph.faceGroup(neighbor));
			}
		}
	}
	std::cerr << "Done\nFitting planes...\n";
	graph.fitPlanes();
	std::cerr << "Done\nSimplifying graph...\n";
	while (graph.simplifyGraph());	//simplify graph until no changes are made
	std::cerr << "Done\n";
}


bool MeshTools::faceIsGrouped(const Mesh::FaceHandle fh, const TopologyGraph& graph) const {
	return graph.faceGroup(fh) != -1;
}


bool MeshTools::normalsAreCloseEnough(const Mesh::Normal& n_1, const Mesh::Normal& n_2, float threshold) const {
	return dot(n_1, n_2) > sin(threshold);
}


void MeshTools::extendNeighborhood(const Mesh::FaceHandle fh, std::list<Mesh::FaceHandle>& neighbors) {
	std::list<Mesh::FaceHandle> extended_neighborhood;
	MeshUtils::getFaceNeighbors(*mesh_, fh, extended_neighborhood);
	neighbors.splice(neighbors.end(), extended_neighborhood);
}


Quadric EdgeCollapse::computeQuadric(const Equation::Plane& plane) {
	return plane.parameters * plane.parameters.transpose();
}


Quadric EdgeCollapse::computeVertexQuadric(const Mesh::VertexHandle vh) {
	Quadric quadric{};
	for (auto& f_it{ mesh->cvf_iter(vh) }; f_it; ++f_it) {
		Equation::Plane face_plane{ MeshUtils::computeFacePlane(*mesh, f_it) };
		quadric += computeQuadric(face_plane);
	}
	return quadric;
}


Quadric EdgeCollapse::computeEdgeQuadric(const Mesh::VertexHandle vh_0, const Mesh::VertexHandle vh_1) {
	return vertexQuadric(vh_0) + vertexQuadric(vh_1);
}


void EdgeCollapse::computeVerticesQuadrics() {
	for (auto v_it{ mesh->vertices_begin() }; v_it != mesh->vertices_end(); ++v_it) {
		if (graph->allNeighborFacesAreInSameRegion(v_it)) {
			Mesh::FaceHandle fh{ mesh->vf_iter(v_it) };
			const RegionID id{ graph->faceGroup(fh) };
			vertexQuadric(v_it) = region_quadrics.at(id);
		}
		else {
			vertexQuadric(v_it) = computeVertexQuadric(v_it);
		}
	}
}


float EdgeCollapse::computeCollapseError(const Vector3f& v, const Quadric& e_quadric) {
	Vector4f extended_v;
	extended_v << v, 1;
	return (extended_v.transpose() * e_quadric * extended_v).value();
}


EdgeCollapse::CollapseResult EdgeCollapse::computeCollapseResult(const Mesh::VertexHandle vh_0, const Mesh::VertexHandle vh_1)
{
	Quadric e_quadric{ computeEdgeQuadric(vh_0, vh_1) };
	Eigen::Matrix4f coeff_mat{ e_quadric };
	coeff_mat.row(3) << 0, 0, 0, 1;
	Eigen::FullPivLU<Eigen::Matrix4f> lu{ coeff_mat };
	Vector3f new_vertex;
	float cost;

	if (lu.isInvertible()) {	// We pick either of the two vertices or their middle
		const Vector3f p0{ MeshUtils::toEigen(mesh->point(vh_0)) };
		const Vector3f p1{ MeshUtils::toEigen(mesh->point(vh_1)) };
		const Vector3f middle{ (p0 + p1) / 2 };
		const std::vector<float> errors{ computeCollapseError(p0, e_quadric), computeCollapseError(p1, e_quadric), computeCollapseError(middle, e_quadric) };
		const float* error_for_v0{ &errors[0]};
		const float* error_for_v1{ &errors[1]};
		const float* error_for_middle{ &errors[2]};

		const float* min_error{ &(*std::min_element(errors.begin(), errors.end())) };
		if (min_error == error_for_v0) {
			new_vertex = p0;
		}
		else if (min_error == error_for_v1) {
			new_vertex = p1;
		}
		else {
			new_vertex = middle;
		}
		cost = *min_error;
	}
	else {	// In this case we have to compute a new vertex
		const Vector4f vertex_4dim = lu.solve(Vector4f{ 0, 0, 0, 1 });
		new_vertex = vertex_4dim.head(2);
		cost = computeCollapseError(new_vertex, e_quadric);
	}
	return { new_vertex, cost };
}


	}
	return new_vertex;
}


void EdgeCollapse::computeRegionQuadrics() {
	this->region_quadrics = {};
	std::set<RegionID> regionIDs{ graph->getRegionIDs() };
	for (auto const& id : regionIDs) {
		region_quadrics[id] = computeQuadric(graph->getPlane(id));
	}
}