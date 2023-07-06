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


void MeshTools::simplifyMesh(TopologyGraph& graph, int targetVerticesAmount) {
	EdgeCollapse ec{ *mesh_, graph };
	mesh_->request_vertex_status();
	mesh_->request_edge_status();
	mesh_->request_face_status();

	while (mesh_->n_vertices() > targetVerticesAmount) {
		std::cerr << mesh_->n_vertices() << '\n';
		ec.collapse();
	}	
	mesh_->update_normals();

	mesh_->release_vertex_status();
	mesh_->release_edge_status();
	mesh_->release_face_status();
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
	if (graph->allNeighborFacesAreInSameRegion(vh)) {
		Mesh::FaceHandle fh{ mesh->vf_iter(vh) };
		const RegionID id{ graph->faceGroup(fh) };
		return region_quadrics.at(id);
	}
	else {
		Quadric quadric{};
		for (auto& f_it{ mesh->cvf_iter(vh) }; f_it; ++f_it) {
			Equation::Plane face_plane{ MeshUtils::computeFacePlane(*mesh, f_it) };
			quadric += computeQuadric(face_plane);
		}
		return quadric;
	}
}


Quadric EdgeCollapse::computeEdgeQuadric(const Mesh::VertexHandle vh_0, const Mesh::VertexHandle vh_1) {
	return vertexQuadric(vh_0) + vertexQuadric(vh_1);
}


void EdgeCollapse::computeVerticesQuadrics() {
	for (auto v_it{ mesh->vertices_begin() }; v_it != mesh->vertices_end(); ++v_it) {
		vertexQuadric(v_it) = computeVertexQuadric(v_it);
	}
}


float EdgeCollapse::computeCollapseError(const Vector3f& v, const Quadric& e_quadric) {
	Vector4f extended_v;
	extended_v << v[0], v[1], v[2], 1;
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
		new_vertex = vertex_4dim.head(3);
		cost = computeCollapseError(new_vertex, e_quadric); //TODO check if removing last row messes with the cost computation
	}
	return { new_vertex, cost };
}


void EdgeCollapse::computePotentialCollapses() {
	for (auto& v_it{ mesh->vertices_begin() }; v_it != mesh->vertices_end(); ++v_it) {
		for (auto& vv_it{ mesh->vv_iter(v_it) }; vv_it; ++vv_it) {
			potentialCollapse(v_it, vv_it) = {v_it, vv_it, computeCollapseResult(v_it, vv_it)};
			collapses.insert(&potentialCollapse(v_it, vv_it));
		}
	}
}


void EdgeCollapse::collapse() {
	Collapse* collapse{ *collapses.begin()};
	collapses.erase(collapses.begin());
	if (!collapse->vh_0.is_valid() || !collapse->vh_1.is_valid())	return;
	Vector3f result{ collapse->result.vertex };
	Mesh::HalfedgeHandle hh{ MeshUtils::findHalfedge(*mesh, collapse->vh_0, collapse->vh_1)};
	if (!hh.is_valid())	return;
	Mesh::Point created_vertex{ result[0], result[1], result[2]};
	Mesh::VertexHandle vh{ mesh->to_vertex_handle(hh) };

	mesh->collapse(hh);
	mesh->set_point(vh, created_vertex);

	mesh->garbage_collection();

	std::set<Mesh::VertexHandle> modified_vertices{ MeshUtils::getNeighboringVertices(*mesh, vh) };
	updateVertices(modified_vertices);
	updatePotentialCollaspes(vh);
}


void EdgeCollapse::updateVertex(const Mesh::VertexHandle vh) {
	vertexQuadric(vh) = computeVertexQuadric(vh);
}


void EdgeCollapse::updateVertices(const std::set<Mesh::VertexHandle>& vhs) {
	for (auto const& vh : vhs) {
		updateVertex(vh);
	}
}


void EdgeCollapse::updatePotentialCollaspes(const Mesh::VertexHandle vh) {
	for (auto& vih_it{ mesh->vih_iter(vh) }; vih_it; ++vih_it) {
		Mesh::HalfedgeHandle hh{ mesh->halfedge_handle(vih_it) };
		Mesh::VertexHandle vh_n{ mesh->from_vertex_handle(hh) };
		Mesh::EdgeHandle eh{ mesh->edge_handle(hh) };

		auto collapse_node{ collapses.extract(&potentialCollapse(eh)) };
		potentialCollapse(eh) = Collapse{ vh, vh_n, computeCollapseResult(vh, vh_n) };
		collapse_node.value() = &potentialCollapse(eh);
		collapses.insert(std::move(collapse_node));
	}
}


void EdgeCollapse::computeRegionQuadrics() {
	this->region_quadrics = {};
	std::set<RegionID> regionIDs{ graph->getRegionIDs() };
	for (auto const& id : regionIDs) {
		region_quadrics[id] = computeQuadric(graph->getPlane(id));
	}
}

