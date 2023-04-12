
//=============================================================================
//
// Often used abbreviations for OpenMesh entitites.
//
//=============================================================================

#ifndef OMAbbreviations_HH
#define OMAbbreviations_HH

#define import_om_abbreviations( Mesh )    \
\
  typedef Mesh::Point Point;   \
  typedef Mesh::Scalar Scalar; \
\
  typedef Mesh::VertexHandle   VH; \
  typedef Mesh::EdgeHandle     EH; \
  typedef Mesh::HalfedgeHandle HH; \
  typedef Mesh::FaceHandle     FH; \
\
  typedef Mesh::VertexIter     VI; \
  typedef Mesh::HalfedgeIter   HI; \
  typedef Mesh::EdgeIter       EI; \
  typedef Mesh::FaceIter       FI; \
\
  typedef Mesh::Vertex    Vertex;   \
  typedef Mesh::Halfedge  Halfedge; \
  typedef Mesh::Edge      Edge;     \
  typedef Mesh::Face      Face;     \
\
  typedef Mesh::ConstVertexIter    CVI; \
  typedef Mesh::ConstHalfedgeIter  CHI; \
  typedef Mesh::ConstEdgeIter      CEI; \
  typedef Mesh::ConstFaceIter      CFI; \
\
  typedef Mesh::VertexVertexIter    VVI;  \
  typedef Mesh::VertexOHalfedgeIter VOHI; \
  typedef Mesh::VertexIHalfedgeIter VIHI; \
  typedef Mesh::VertexEdgeIter      VEI;  \
  typedef Mesh::VertexFaceIter      VFI;  \
  typedef Mesh::FaceVertexIter      FVI;  \
  typedef Mesh::FaceHalfedgeIter    FHI;  \
  typedef Mesh::FaceEdgeIter        FEI;  \
  typedef Mesh::FaceFaceIter        FFI;  \
\
  typedef Mesh::ConstVertexVertexIter    CVVI;  \
  typedef Mesh::ConstVertexOHalfedgeIter CVOHI; \
  typedef Mesh::ConstVertexIHalfedgeIter CVIHI; \
  typedef Mesh::ConstVertexEdgeIter      CVEI;  \
  typedef Mesh::ConstVertexFaceIter      CVFI;  \
  typedef Mesh::ConstFaceVertexIter      CFVI;  \
  typedef Mesh::ConstFaceHalfedgeIter    CFHI;  \
  typedef Mesh::ConstFaceEdgeIter        CFEI;  \
  typedef Mesh::ConstFaceFaceIter        CFFI;  \

#endif
