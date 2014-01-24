/*
Step: [0-1]
StepSink: Step ->
StepSource: StepSink ->
  step(n): stepSource // produce n linear steps between 0 and 1 inclusive
VertexSink: Vertex, Step, Step* ->
VertexGenerator: VertexSink -> StepSink
  vertices(vertices) // take the given vertices to create a VertexGenerator
  vertices(f): (Step -> Vertex) -> VertexGenerator // convert a vertex producer into a vertex generator 
VertexMap: VertexSink -> VertexSink
  zTranslate(start,end): translate vertices along z axis
  zRotate(start*,end*): rotate vertices around z axis
  tag(edges): mark vertices with tags
Facer: FaceSink -> VertexSink
  face(options): Produce a facer with the given options for tesselating and joining
Renderer: * -> FaceSink
  ThreeJSRenderer
  CSGRenerer
  STLRenderer
  NativeRenderer
*/

// cube
var renderer = ThreeJSRenderer();

// cube
step(2)(
  vertices([[1,1,0],[-1,1,0],[-1,-1,0],[1,-1,0]])(
  transform(zTranslate(1), label) (
  face( edgeLoop, tesselate('bottom'), tesselate('top',true) )( renderer ))))

var geometry = renderer.geometry;

// OR
extrudeClosed = face( edgeLoop, tesselate('bottom'), tesselate('top',true) )

// sphere
var R = 40, Q = 20, 2PI = 2*Math.Pi, hPI = Math.PI/2;
function simiCircle(s) {return [R*Math.cos(tween(s,-hPI,hPI)),0,R*Math.sin(tween(s,-hPI,hPI))]}

step(Q)(
  step(Q)(
    parametric(semiCircle))(
      zRotate()(
        face(connect:'tb', singular:'lr')( renderer ))))

// STEP
function step(iterations) {
  return function(stepSink) {
    for (var i=0; i < iterations; i++) stepSink( i / (iterations - 1) );
  }
}

// VERTICES
function vertices(points) {
  return function(vertexSink) {
    return function(step) {
      for (var i=0,p; p=points[i]; i++)
        vertexSink(p, step, i / (points.length-1));
    }
  }
}

function parametric(f) {
  return function(vertexSink) {
    return function( ribStep ) {
      return function( transformStep ) {
        vertexSink( f(ribStep, transformStep), transformStep, ribStep );
      }
    }
  }
}

/// TESSELATE

function vadd(a,v) { return [a[0]+v[0], a[1]+v[1], a[2]+v[2]]; }
function vsub(a,v) { return [a[0]-v[0], a[1]-v[1], a[2]-v[2]]; }
function vscale(a,c) { return [a[0]*c, a[1]*c, a[2]*c]; }
function vdot(a,v) { return a[0]*v[0] + a[1]*v[1] + a[2]*v[2]; }
function vcross(v) { return [a[1]*v[2] - a[2]*v[1], a[2]*v[0] - a[0]*v[2], a[0]*v[1] - a[1]*v[0]]; }
function vlength(v) { return Math.sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]); }
function vnorm(v) { var l=vlength(v); return v > 0 ? [v[0]/l, v[1]/l, v[2]/l] : v; }

function vectorAverage(vs) { return vs.length ? vscale( vs.reduce(vadd,[0,0,0]), 1/vs.length ) : [0,0,0]; }

function loopMeanNormal(vs) { 
  function v(i) { return vs[i%vs.length] }
  return vectorAverage(vs.map(function(v0,i){return vcross(vsub(v0,v(i+1)), vsub(v(i+2), v(i+1)))})) 
}

function intersect(v1, v2, v3, v4) {
  var det = (v1.x-v2.x)*(v3.y-v4.y) - (v1.y-v2.y)*(v3.x-v4.x);
  return !det ? null : {
    x: ((v1.x*v2.y-v2.x*v1.y)*(v3.x-v4.x) - (v3.x*v4.y-v4.x*v3.y)*(v1.x-v2.x))/det,
    y: ((v1.x*v2.y-v2.x*v1.y)*(v3.y-v4.y) - (v3.x*v4.y-v4.x*v3.y)*(v1.y-v2.y))/det
  };
}

function doesIntersect(v1, v2, v3, v4) {
  var p = intersect(v1,v2,v3,v4);
  if (!p) return false;
  return p.x >= Math.min(v1.x,v2.x) && p.x >= Math.min(v3.x,v4.x)
    && p.x <= Math.max(v1.x,v2.x) && p.x <= Math.max(v3.x,v4.x)
    && p.y >= Math.min(v1.y,v2.y) && p.y >= Math.min(v3.y,v4.y)
    && p.y <= Math.max(v1.y,v2.y) && p.y <= Math.max(v3.y,v4.y);
}

function isSelfIntersecting( p1, p2, vertices ) {
  for (var i=0,v1,v2; v1=vertices[i], v2=vertices[i+1]; i++)
    if ( doesIntersect(p1,p2,v1,v2) ) return true;
  return false;
}

function projectIntoPlane( p, planePoint, planeNormal ) {
  // p - pp - (((p-pp) . pn) * pn)
  var pTranslated = vsub( p, planePoint );
  return vsub( pTranslated, vscale( planeNormal, vdot( pTranslated, planeNormal )))
}

function tesselate3d( points3d, faceSink ) {
  var planeNormal = vnorm(loopMeanNormal(points3d));
  var planePoint = vectorAverage(points3d);

  var offNormal = vadd( planePoint, Math.abs(planeNormal[0] < .5) ? [1,0,0] : [0,1,0] );
  var axis1 = vnorm(projectIntoPlane(offNormal ,planePoint, planeNormal));
  var axis2 = vcross( axis1, planeNormal );

  var vertices2d = points3d.map( function(p) {
    var pointInPlane = projectIntoPlane(p, planePoint, planeNormal);
    return [vdot(pointInPlane,axis1), vdot(pointInPlane,axis2)];
  })

  tesslateVertices( vertices2d, function(face) {
    faceSink(face);
  })
}

function tesselateVertices(vertices, faceSink) {
  if ( vertices.length < 3 ) return;
  if ( vertices.length < 4 ) return tesselateConvex(vertices);
    
  // traveling in either direction, find the first negative vertex followed by a positive one.
  // If there are no negative vertices, the shape is convex. Tesselate and return.
  // Starting with the vertex after that, find the first negative vertex or the last vertext s.t.
  //   the angle formed from it and the first two vertices is < π.
  // Create two new sets of vertices, one containing all the vertices found above, and the
  //   other containing the rest plus the first and last from above.
  // The first set is convex, so pass it on to be tesselated.
  // The second list is "less concave", because the points at the ends now have less negative angles.
  // Recurse with the second list. 
  function index(i) { return (i+vertices.length)%vertices.length;  }
  function center(i) { return vertices[ i%vertices.length ] }
  function left(i) { return vertices[ (i+vertices.length-1)%vertices.length ] }
  function right(i) { return vertices[ (i+1)%vertices.length ] }
  function angle(l,p,r) { return Math.atan2( 
    (r.x-p.x)*(l.y-p.y) - (r.y-p.y)*(l.x-p.x),
    (r.x-p.x)*(l.x-p.x) + (r.y-p.y)*(l.y-p.y)
  ) }

  var angles = vertices.map( function(p,i) { return angle(left(i),p,right(i)); })
  
  var startPoint;
  for (var i=vertices.length-1; i>=0; i--) {
    if ( angles[i] <= 0 && angles[index(i+1)] > 0 ) {
      startPoint = i;

      for (var j=2, endPoint; j<vertices.length-1; j++) {
      endPoint = j + startPoint;
        if ( angle(left(endPoint),center(endPoint),right(endPoint)) <= 0 )
          break;

        if ( angle(center(endPoint),center(startPoint), right(startPoint)) <= 0 ) {
          endPoint = (endPoint+vertices.length-1) % vertices.length;
          break;
        }
      }
      
      endPoint = index(endPoint);

      var newShape, rest;
      if (endPoint < startPoint) {
        newShape = vertices.slice(startPoint).concat(vertices.slice(0,endPoint+1));
        rest = vertices.slice(endPoint, startPoint+1);
      } else {
        newShape = vertices.slice(startPoint, endPoint+1);
        rest = vertices.slice(endPoint).concat(vertices.slice(0,startPoint+1));
      }

      if ( ! isSelfIntersecting( vertices[startPoint], vertices[endPoint], rest.slice(1,rest.length-1) ) ) {
        if ( Math.abs(endPoint- startPoint) <= 1 || Math.abs(endPoint- startPoint) >= vertices.length -1 ) {
          console.log("wrap around", startPoint, endPoint, vertices.length);
          return;
        }
        if (endPoint < startPoint) {
          tesselate( vertices.slice(startPoint).concat(vertices.slice(0,endPoint+1)) );
          tesselate( vertices.slice(endPoint, startPoint+1 ) );
        } else {
          tesselate( vertices.slice(startPoint, endPoint+1 ) );
          tesselate( vertices.slice(endPoint).concat(vertices.slice(0,startPoint+1)) );
        }

        ctx.strokeStyle = '#fff';
        ctx.lineWidth = 2;
        var p = center(startPoint);
        ctx.beginPath()
        ctx.moveTo( p.x, p.y );
        p = center(endPoint);
        ctx.lineTo( p.x, p.y );
        ctx.stroke();  
        return;
      }
    }
  }
  
  return tesselateConvex(vertices);
  
}

function tesselateConvex(vertices) {
  var ai = 0, bi=(vertices.length/3)|0, ci=(2*vertices.length/3)|0;
  var a = vertices[ai], b=vertices[bi], c=vertices[ci];

  ctx.fillStyle = randomColor();
  ctx.strokeStyle = '#000';
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(a.x,a.y);
  ctx.lineTo(b.x,b.y);
  ctx.lineTo(c.x,c.y);
  ctx.closePath();
  ctx.fill();
  ctx.stroke();
  
  if (bi-ai > 1) tesselateConvex(vertices.slice(ai,bi+1));
  if (ci-bi > 1) tesselateConvex(vertices.slice(bi,ci+1));
  if (vertices.length-ci > 1) 
    tesselateConvex(vertices.slice(ci,vertices.length).concat([a]));
}
