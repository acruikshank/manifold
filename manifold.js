// renderer: face -> [geometry model, three.js, csg.js, stl, etc.]
// manifold: [0-1] -> [face]
// vertexGen: [0-1] -> vertex


/*
vertexHandler: vertex ->
transformer([0-1]) -> 3dTransform
facer(renderer) -> vertexHandler // vertices have row data
manifold: indexedVertexHandler -> vertexHandler
rib: vertexHandler ->
*/


XYPointRib([[1,1],[1,-1],[-1,-1],[-1,1]])(
  linearZ(0,4), )
) 