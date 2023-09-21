#include "MeshTools.h"
#include <list>

void MeshTools::extractRegions(TopologyGraph& graph) {
	//initialization
	std::set<OpenMesh::FaceHandle> ungrouped_faces{};
	for (auto f_iter{ mesh_->faces_begin() }; f_iter != mesh_->faces_end(); ++f_iter) {
		graph.faceGroup(f_iter) = -1;
		ungrouped_faces.insert(f_iter);
	}
	std::cerr << "Executing region growing algorithm...\n";
	growRegions(ungrouped_faces, graph);
	std::cerr << "Done\nBuilding topology graph...\n";
	buildTopologyGraph(graph);
}


void MeshTools::simplifyMesh(TopologyGraph& graph, int targetVerticesAmount) {
	mesh_->request_vertex_status();
	mesh_->request_edge_status();
	mesh_->request_face_status();
	

	for (EdgeCollapse ec{ *mesh_, graph }; ec.n_vertices() > targetVerticesAmount; ) {
		ec.collapse();
	}
	for (auto& v_it{ mesh_->vertices_begin() }; v_it != mesh_->vertices_end(); ++v_it) {
		if (mesh_->status(v_it).deleted()) {
			continue;
		}
		mesh_->update_normal(v_it);
	}
	mesh_->garbage_collection();


	mesh_->release_face_status();
	mesh_->release_edge_status();
	mesh_->release_vertex_status();
}


void MeshTools::growRegions(std::set<Mesh::FaceHandle>& ungrouped_faces, TopologyGraph& graph) {
	int current_region_id{-1};
	while (!ungrouped_faces.empty()) {
		Mesh::FaceHandle fh{ ungrouped_faces.extract(ungrouped_faces.begin()).value()};	// Pop first ungrouped face
		graph.addFaceToRegion(++current_region_id, fh);

		std::list<Mesh::FaceHandle>  neighbors;
		MeshUtils::getFaceNeighbors(*mesh_, fh, neighbors);
		std::cerr << "Region " << current_region_id << '\n';

		for (auto f_neighbor{ neighbors.begin() }; f_neighbor != neighbors.end(); ) {
			Mesh::Normal f_neighbor_normal{ mesh_->normal(*f_neighbor)};
			Mesh::Normal f_normal{ mesh_->normal(fh)};

			if (faceIsGrouped(*f_neighbor, graph)) {
				f_neighbor = neighbors.erase(f_neighbor);
				continue;
			}
			else if (normalsAreCloseEnough(f_neighbor_normal, f_normal, 0.349066f)) {
				graph.addFaceToRegion(current_region_id, *f_neighbor);
				extendNeighborhood(*f_neighbor, neighbors);
				ungrouped_faces.erase(*f_neighbor);
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
	graph.computeFaceAreas();
	graph.computeVertexProjectedDistances();
	std::cerr << "\tSimplifying small regions...\t";
	while (graph.simplifySmallRegions());	//simplify graph until no changes are made
	std::cerr << "Done\n";
	std::cerr << "\tSimplifying large regions...\t";
	while (graph.simplifyLargeRegions());
	std::cerr << "Done\n";
}


bool MeshTools::faceIsGrouped(const Mesh::FaceHandle fh, const TopologyGraph& graph) const {
	return graph.faceGroup(fh) != -1;
}


bool MeshTools::normalsAreCloseEnough(const Mesh::Normal& n_1, const Mesh::Normal& n_2, float threshold) const {
	return dot(n_1, n_2) > cos(threshold);	//The higher the threshold, the closer the normals have to be?
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
		Quadric quadric{Eigen::Matrix4f::Zero()};
		for (auto& f_it{ mesh->cvf_iter(vh) }; f_it; ++f_it) {
			if (mesh->status(f_it).deleted())
				continue;
			Equation::Plane face_plane{ MeshUtils::computeFacePlane(*mesh, f_it) };
			quadric += computeQuadric(face_plane);
		}
		return quadric;
	}
}


Quadric EdgeCollapse::computeEdgeQuadric(const Mesh::VertexHandle vh_0, const Mesh::VertexHandle vh_1) {
	auto v0{ vertexQuadric(vh_0) };
	auto v1{ vertexQuadric(vh_1) };
	return vertexQuadric(vh_0) + vertexQuadric(vh_1);
}


void EdgeCollapse::computeVerticesQuadrics() {
	for (auto v_it{ mesh->vertices_begin() }; v_it != mesh->vertices_end(); ++v_it) {
		if (mesh->status(v_it).deleted() || mesh->is_isolated(v_it))
			continue;
		vertexQuadric(v_it) = computeVertexQuadric(v_it);
	}
}


float EdgeCollapse::computeCollapseError(const Vector3f& v, const Quadric& e_quadric) {
	Vector4f extended_v{};
	extended_v << v[0], v[1], v[2], 1;
	return (extended_v.transpose() * e_quadric * extended_v).value();
}


float EdgeCollapse::computeCollapseError(const Vector4f& v, const Quadric& e_quadric) {
	//Vector4f extended_v;
	//extended_v << v[0], v[1], v[2], 1;
	return (v.transpose() * e_quadric * v).value();
}


EdgeCollapse::CollapseResult EdgeCollapse::computeCollapseResult(const Mesh::VertexHandle vh_0, const Mesh::VertexHandle vh_1)
{
	Quadric e_quadric{ computeEdgeQuadric(vh_0, vh_1) };
	Eigen::Matrix4f coeff_mat{ e_quadric };
	coeff_mat.row(3) = Vector4f{ 0, 0, 0, 1 };
	Eigen::FullPivLU<Eigen::Matrix4f> lu{ coeff_mat };
	Vector3f new_vertex;
	float cost;

	if (!lu.isInvertible()) {	// We pick either of the two vertices or their middle
		const Vector3f p0{ MeshUtils::toEigen(mesh->point(vh_0)) };
		const Vector3f p1{ MeshUtils::toEigen(mesh->point(vh_1)) };
		const Vector3f middle{ (p0 + p1) / 2 };
		const std::vector<float> errors{ computeCollapseError(p0, e_quadric), computeCollapseError(p1, e_quadric), computeCollapseError(middle, e_quadric) };
		const float* error_for_v0{ &errors[0]};
		const float* error_for_v1{ &errors[1]};
		const float* error_for_middle{ &errors[2]};

		const float* min_error{ findMinError(errors, 0.000001) };
		if (min_error == error_for_v0) {
			new_vertex = p0;
		}
		else if (min_error == error_for_v1) {
			new_vertex = p1;
		}
		else {
			new_vertex = middle;
		}
		//new_vertex = middle;
		cost = *min_error;
		//cost = *error_for_middle;
	}
	else {	// In this case we have to compute a new vertex
		const Vector4f vertex_4dim = lu.solve(Vector4f{ 0, 0, 0, 1 });
		cost = computeCollapseError(vertex_4dim, e_quadric); 
		new_vertex = vertex_4dim.head(3);
	}
	return { new_vertex, cost };
}


void EdgeCollapse::computePotentialCollapses() {
	for (auto& v_it{ mesh->vertices_begin() }; v_it != mesh->vertices_end(); ++v_it) {
		if (mesh->status(v_it).deleted())
			continue;
		for (auto& vv_it{ mesh->vv_iter(v_it) }; vv_it; ++vv_it) {
			if (mesh->status(vv_it).deleted())
				continue;
			computePotentialCollapse(v_it, vv_it);
		}
	}
}


void EdgeCollapse::computePotentialCollapse(const Mesh::VertexHandle vh_0, const Mesh::VertexHandle vh_1) {
	potentialCollapse(vh_0, vh_1) = { vh_0, vh_1, computeCollapseResult(vh_0, vh_1) };
	float cost{ potentialCollapse(vh_0, vh_1).result.cost };
	if (collapses.contains(cost)) {
		collapses.at(potentialCollapse(vh_0, vh_1).result.cost).push_back(potentialCollapse(vh_0, vh_1));
	}
	else {
		collapses.insert({ cost, std::deque<Collapse>{potentialCollapse(vh_0, vh_1) } });
	}
}


void EdgeCollapse::collapse() {
	std::deque<Collapse>& c_arr{ collapses.begin()->second };
	Collapse collapse{ c_arr[0]};
	c_arr.pop_front();
	if (c_arr.empty()) {
		collapses.erase(collapses.begin());
	}

	//std::cerr << "Cost: " << collapse.result.cost << '\n';

	Mesh::VertexHandle vh_0{ collapse.edge.vh_0};
	Mesh::VertexHandle vh_1{ collapse.edge.vh_1 };

	if (mesh->status(vh_0).deleted() || mesh->status(vh_1).deleted()) {
		//std:cerr << "Potential collapse contained deleted vertex. Skipping.\n";
		return;
	}

	if (!collapse.edge.vh_0.is_valid() || !collapse.edge.vh_1.is_valid()) {
		return;
	}

	if (vh_0 == vh_1) {
		//std::cerr << "Potential collapse between two same vertices. Skipping.\n";
		return;
	}

	Vector3f result{ collapse.result.vertex };
	Mesh::HalfedgeHandle hh{ MeshUtils::findHalfedge(*mesh, vh_0, vh_1)};
	if (!hh.is_valid())	
		return;

	Mesh::VertexHandle del_vh{ mesh->from_vertex_handle(hh) };
	Mesh::Point created_vertex{ result[0], result[1], result[2]};
	Mesh::VertexHandle vh{ mesh->to_vertex_handle(hh) };

	if (!mesh->is_collapse_ok(hh)) {
		std::cerr << "Collapse not OK. Skipping\n";
		skipped_vertices.insert(vh_0);
		skipped_vertices.insert(vh_1);
		return;
	}

	mesh->set_point(vh, created_vertex);
	mesh->set_point(del_vh, created_vertex);

	mesh->collapse(hh);
	++removed_vertices;

	std::set<Mesh::VertexHandle> modified_vertices{ MeshUtils::getNeighboringVertices(*mesh, vh)};
	modified_vertices.insert(vh);
	modified_vertices.merge(skipped_vertices);

	std::cerr << "Merged " << del_vh.idx() << " into " << vh.idx() << '\n';
	std::cerr << "New vertex : " << created_vertex << '\n';

	updateVertices(modified_vertices);
	updatePotentialCollapses(modified_vertices);
}


void EdgeCollapse::updateVertex(const Mesh::VertexHandle vh) {
	vertexQuadric(vh) = computeVertexQuadric(vh);
}


void EdgeCollapse::updateVertices(const std::set<Mesh::VertexHandle>& vhs) {
	for (auto const& vh : vhs) {
		updateVertex(vh);
	}
}


unsigned int EdgeCollapse::n_vertices() {
	return mesh->n_vertices() - removed_vertices;
}


void EdgeCollapse::updatePotentialCollapses(const std::set<Mesh::VertexHandle>& vhs) {
	removeOutdatedCollapses(vhs);

	for (auto& v_it : vhs) {	// Recompute all operations on potentially modified edges		
		for (auto& vv_it{ mesh->vv_iter(v_it) }; vv_it; ++vv_it) {
			computePotentialCollapse(v_it, vv_it);
		}
	}
}


void EdgeCollapse::removeOutdatedCollapses(const std::set<Mesh::VertexHandle>& modified_vertices) {
	for (auto& collapses_it{ collapses.begin() }; collapses_it != collapses.end(); ) {
		std::deque<Collapse> collapse_list{ collapses_it->second };
		for (auto& collapse_it{ collapse_list.begin() }; collapse_it != collapse_list.end(); ) {
			const Edge& edge{ collapse_it->edge };
			Mesh::EdgeHandle eh{ MeshUtils::findEdge(*mesh, edge.vh_0, edge.vh_1) };

			if (!eh.is_valid() || mesh->status(eh).deleted() || mesh->status(edge.vh_0).deleted() || mesh->status(edge.vh_1).deleted()) {	// If operation on deleted edge
				collapse_it = collapse_list.erase(collapse_it);
				continue;
			}
			else if (modified_vertices.contains(edge.vh_0) || modified_vertices.contains(edge.vh_1)) {	// Remove operations on obsolete edges
				collapse_it = collapse_list.erase(collapse_it);
				continue;
			}
			++collapse_it;
		}
		if (collapse_list.empty()) {
			collapses_it = collapses.erase(collapses_it);
			continue;
		}
		++collapses_it;
	}
}


void EdgeCollapse::computeRegionQuadrics() {
	this->region_quadrics = {};
	std::set<RegionID> regionIDs{ graph->getRegionIDs() };
	for (auto const& id : regionIDs) {
		region_quadrics[id] = computeQuadric(graph->getPlane(id));
	}
}


bool EdgeCollapse::lessThan(float a, float b, float epsilon) {
	return (b - a) > ((fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}


const float* EdgeCollapse::findMinError(const std::vector<float>& errors, const float epsilon) {
	assert(!errors.empty());
	const float* min{ &errors.at(0)};
	if (errors.size() == 1) {
		return min;
	}
	else {
		for (auto& it{ errors.cbegin()+1 }; it != errors.cend(); ++it) {
			if (lessThan(*it, *min, epsilon)) {
				min = &(*it);
			}
		}
	}
	return min;
}