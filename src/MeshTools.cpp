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



		}
	}
}


void computeEdgeErrors(const Mesh& mesh) {

}


void EdgeCollapse::computeRegionQuadrics() {
	this->region_quadrics = {};
	std::set<RegionID> regionIDs{ graph->getRegionIDs() };
	for (auto const& id : regionIDs) {
		region_quadrics[id] = computeQuadric(graph->getPlane(id));
	}
}