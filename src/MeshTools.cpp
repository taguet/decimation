#include "MeshTools.h"
#include <math.h>


MeshTools::MeshTools(void) {
	
}

MeshTools::MeshTools(Mesh& mesh) {
	this->mesh_ = &mesh;
}

void MeshTools::setMesh(Mesh& mesh) {
	this->mesh_ = &mesh;
}


