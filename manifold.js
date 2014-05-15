/*
Step: [0-1]
StepSink: Step ->
StepSource: StepSink ->
  step(n): stepSource // produce n linear steps between 0 and 1 inclusive
VertexSink: Vertex, Step, Step* ->
VertexGenerator: VertexSink -> StepSink
  vertices(vertices) // take the given vertices to create a VertexGenerator
  parametric(f): (Step -> Vertex) -> VertexGenerator // convert a vertex producer into a vertex generator 
VertexMap: VertexSink -> VertexSink
  zTranslate(start,end): translate vertices along z axis
  zRotate(start*,end*): rotate vertices around z axis
  tag(edges): mark vertices with tags
Facer: FaceSink -> VertexSink
  face(options): Produce a facer with the given options for tesselating and joining
Renderer: FaceSink
  ThreeJSRenderer
  CSGRenerer
  STLRenderer
  NativeRenderer
*/

(function(context) {
  /*
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


  */

  function Vertex( point, transformStep, ribStep, id) {
    this[0] = point[0];
    this[1] = point[1];
    this[2] = point[2];
    this[3] = transformStep;
    this[4] = ribStep;
    this[5] = id;
  }
  Vertex.prototype = new Float64Array(6);
  Object.defineProperty(Vertex.prototype, 'x', {get:function() {return this[0]}, set:function(x) {this[0]=x}})
  Object.defineProperty(Vertex.prototype, 'y', {get:function() {return this[1]}, set:function(y) {this[1]=y}})
  Object.defineProperty(Vertex.prototype, 'z', {get:function() {return this[2]}, set:function(z) {this[2]=z}})
  Object.defineProperty(Vertex.prototype, 'transformStep', {
    get:function() {return this[3]}, 
    set:function(transformStep) {this[3]=transformStep}})
  Object.defineProperty(Vertex.prototype, 'ribStep', {
    get:function() {return this[4]}, 
    set:function(ribStep) {this[4]=ribStep}})
  Object.defineProperty(Vertex.prototype, 'id', {
    get:function() {return this[5]}, 
    set:function(id) {this[5]=id}})


  function resolveCurve( points, s, start, size ) {
    start = start || 0;
    size = size || points.length;
    if ( size < 2 ) return points[start];

    var p1 = resolveCurve( points, s, start, size-1 ), p2 = resolveCurve( points, s, start+1, size-1 );
    return vadd( p1, vscale( vsub(p2, p1), s) );
  }

  function curveTangent( points, s, start, size ) {
    start = start || 0;
    size = size || points.length;
    if ( size < 2 ) return [0,0,0];

    var p1 = resolveCurve( points, s, start, size-1 ), p2 = resolveCurve( points, s, start+1, size-1 );
    return vnorm( vsub( p2, p1 ) );
  }

  function curveCurl( points, s, start, size ) {
    start = start || 0;
    size = size || points.length;
    if ( size < 3 ) return [0,0,0];

    var p1 = resolveCurve( points, s, start, size-2 ), 
        p2 = resolveCurve( points, s, start+1, size-2 ),
        p3 = resolveCurve( points, s, start+2, size-2 );
    return vcross( vsub( p3, p2 ), vsub( p1, p2 ) );
  }

  var pathIndex = 0;
  function Path( start ) {
    var segments = [], points = [start], path = {}, totalWeight = 0;

    path.curve = function(ps, weight) {
      if (weight == null) weight = 1;

      for (var i=0,point;point = ps[i];i++) points.push(point);

      segments.push({s:points.length-ps.length-1, o:ps.length+1, w:weight});

      totalWeight += weight;
      return path;
    }

    path.line = function(point, weight) {
      return path.curve([point], weight);
    }

    path.spline = function(ps, weight) {
      if (points.length < 2) path.curve(ps);

      var p1 = points[points.length-1];
      var p2 = points[points.length-2];
      var ctl = vadd( p1, vscale( vsub(p2,p1), -1) );

      return path.curve( [ctl].concat(ps), weight);
    }

    path.vertices = function(vertexSink, step, numVertices) {
      var allDivisions = Math.max(numVertices, segments.length + 1);
      var remainingDivisions = allDivisions - 1 - segments.length;
      var remainingWeight = totalWeight;
      var index = 0;

      for (var i=0,s; s = segments[i]; i++) {
        var divisions = i >= segments.length - 1 ? remainingDivisions + 1 : (s.w * remainingDivisions / remainingWeight + 1)|0;
        for (var j=0; j<divisions; j++)
          vertexSink( new Vertex(resolveCurve(points,j/divisions, s.s, s.o), step, (index++)/(allDivisions-1), pathIndex++) );

        remainingDivisions -= divisions - 1;
        remainingWeight -= s.w;
      }

      vertexSink( new Vertex(points[points.length-1], step, 1, pathIndex++) );
    }
    
    return path;
  }

  // STEP
  function step(iterations) {
    return function(stepSink) {
      for (var i=0; i < iterations; i++) stepSink( i / (iterations - 1) );
    }
  }

  // VERTICES
  function vertices(points) {
    return function(vertexSink) {
      var index = 0;
      return function(step) {
        for (var i=0,p; p=points[i]; i++)
          vertexSink( new Vertex(p, step, i / (points.length-1), index++) );
      }
    }
  }

  function parametric(f, stepper) {
    return function(vertexSink) {
      var index = 0;
      return function( transformStep ) {
        stepper( function(ribStep) {
          vertexSink( new Vertex(f(ribStep, transformStep), transformStep, ribStep, index++) );
        })
      }
    }
  }

  // TRANSFORM
  function translate(translations) {
    return function( vertexSink ) {
      return function( vertex ) {
        console.log(vertex)
        var tIndex = parseInt(vertex.transformStep * (translations.length-1))
        var translate = vertex.transformStep==1 ? 
          translations[tIndex] : 
          vinterp( translations[tIndex], translations[tIndex+1], vertex.transformStep-tIndex)
        vertexSink( new Vertex(vadd(translate,vertex), vertex.transformStep, vertex.ribStep, vertex.id ) );
      }
    }
  }

  // TESSELATE
  function skin( faceSink ) {
    var lastRib;
    var nextRib = [];
    var bIndex = 0;
    var lastTransformStep = 0;
    return function skinVertexSink( vertex ) {
      if ( vertex.transformStep > lastTransformStep ) {
        lastRib = nextRib;
        nextRib = [];
        bIndex=0;
        lastTransformStep = vertex.transformStep;
      }

      if (lastRib) {
        var blVertex = lastRib[bIndex];
        var brVertex = lastRib[bIndex+1];

        if (nextRib.length) {
          console.log('top face', blVertex.id, nextRib[nextRib.length-1].id, vertex.id)
          faceSink( [blVertex, nextRib[nextRib.length-1], vertex] );
        }

        while ( brVertex &&  vertex.ribStep > blVertex.ribStep + (brVertex.ribStep-blVertex.ribStep)/2 ) {
          faceSink( [blVertex, vertex, brVertex] );
          console.log("face", blVertex.id, vertex.id, brVertex.id)
          blVertex = brVertex;
          brVertex = lastRib[++bIndex + 1];
        }
      }

      nextRib.push(vertex)
    }
  }

  function facers(facer1, facer2, etc) {
    var facers = Array.prototype.slice.call(arguments,0);
    return function facersFacer( faceSink ) {
      var vertexSinks = facers.map( function(f) { return f(faceSink); } );
      return function facersVertexSink( vertex ) {
        vertexSinks.forEach(function(vs) {vs(vertex)});
      }
    }
  }

  function closeEdge( faceSink ) {
    var bottomFirstInRib, bottomLastInRib;
    var topFirstInRib, topLastInRib;
    return function wrapEdgeVertexSink( vertex ) {
      if ( ! topFirstInRib )
        topFirstInRib = vertex;
      else if ( vertex.transformStep == topFirstInRib.transformStep )
        topLastInRib = vertex;

      if ( vertex.ribStep == 1 || vertex.transformStep > topFirstInRib.transformStep ) {
        if ( bottomFirstInRib ) {
          if ( topLastInRib && bottomLastInRib )
            faceSink( [bottomLastInRib, topLastInRib, topFirstInRib] );

          if ( bottomLastInRib ) 
            faceSink( [bottomLastInRib, topFirstInRib, bottomFirstInRib] );
          else if ( topLastInRib )
            faceSink( [bottomFirstInRib, topLastInRib, topFirstInRib] );
        }

        bottomFirstInRib = topFirstInRib;
        bottomLastInRib = topLastInRib;
        topFirstInRib = vertex.transformStep > topFirstInRib.transformStep ? vertex : null;
        topLastInRib = null;
      }
    }
  }

  function capBottom( faceSink ) {
    var rib = [];
    return function capBottomVertexSink( vertex ) {
      if ( vertex.transformStep == 0 )
        return rib.push(vertex);

      if (rib) {
        tesselate(rib, reverseFaceSink(faceSink) );
        rib = null;
      }
    }
  }

  function capTop( faceSink ) {
    var rib = [];
    return function capBottomVertexSink( vertex ) {
      if ( vertex.transformStep < 1 ) return;

      rib.push(vertex);

      if (vertex.ribStep == 1)
        tesselate(rib, faceSink );
    }
  }

  // FACER TRANSFORMS

  // Facer -> Facer
  function reverse(facer) {
    return function( faceSink ) {
      return facer( reverseFaceSink(faceSink) )
    }
  }

  // FaceSink -> FaceSink
  function reverseFaceSink( faceSink ) {
    return function(face) {
      faceSink( [face[2],face[1],face[0]] )
    }
  }

  function vertexString(v) {return v?'{'+v.x+','+v.y+','+v.z+'} ':'null'}
  function faceString(face) { return face.map(vertexString); }

  // RENDER
  function ThreeJSRenderer() {
    var geometry = new THREE.Geometry();
    return {
      renderer : function( face ) {
        saveThreeJSVertex( geometry, face[0] );
        saveThreeJSVertex( geometry, face[1] );
        saveThreeJSVertex( geometry, face[2] );
        geometry.faces.push( new THREE.Face3(face[0].id, face[1].id, face[2].id) );
      },
      geometry: geometry
    }
  }

  function saveThreeJSVertex(geometry, vertex) {
    if (! geometry.vertices[vertex.id])
      geometry.vertices[vertex.id] = new THREE.Vector3( vertex[0], vertex[1], vertex[2] );
  }

  function STLRenderer(){
    var doc = 'solid pixel';
    return {
      renderer : function( face ) {
        var normal = vnorm( vcross( vsub(face[1],face[0]), vsub(face[2],face[0])) );
        doc += "facet normal " + normal.x + " " + normal.y + " " + normal.z + " \n";
        doc += "outer loop \n";
        doc += "vertex " + face[0].x + " " + face[0].y + " " + face[0].z + " \n";
        doc += "vertex " + face[1].x + " " + face[1].y + " " + face[1].z + " \n";
        doc += "vertex " + face[2].x + " " + face[2].y + " " + face[2].z + " \n";
        doc += "endloop \n";
        doc += "endfacet \n";
      },
      doc : function() { return doc + 'endsolid' }
    }
  }


  // MATH

  function vadd(a,v) { return [a[0]+v[0], a[1]+v[1], a[2]+v[2]]; }
  function vsub(a,v) { return [a[0]-v[0], a[1]-v[1], a[2]-v[2]]; }
  function vscale(a,c) { return [a[0]*c, a[1]*c, a[2]*c]; }
  function vdot(a,v) { return a[0]*v[0] + a[1]*v[1] + a[2]*v[2]; }
  function vcross(a,v) { return [a[1]*v[2] - a[2]*v[1], a[2]*v[0] - a[0]*v[2], a[0]*v[1] - a[1]*v[0]]; }
  function vlength(v) { return Math.sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]); }
  function vnorm(v) { var l=vlength(v); return l > 0 ? [v[0]/l, v[1]/l, v[2]/l] : v; }
  function vinterp(a,b,c) { return vadd(a,vscale(vsub(b,a),c)); }

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

  function loop( vertices, index ) {
    return vertices[ (index + vertices.length) % vertices.length ];
  }

  // TESSLATION

  // Assume 2d with l2.x >= v.x >= l1.x
  function vertexIsAboveLine( v, l1, l2 ) {
    return (v[0] == l1[0] && v[1] > l1[1] && (l2[0] > l1[0] || v[1] > l2[1])) 
          || (v[1]-l1[1])/(v[0]-l1[0]) > (l2[1]-l1[1])/(l2[0]-l1[0]);
  }

  function isConvex(v1, v2, v3) {
    return vcross(vsub(v1,v2),vsub(v3,v2))[2] > 0;
  }

  function MonotonePolygon( lowerVertices ) {
    var upper = [];
    var lower = lowerVertices || [];

    function addUpper( v, faceSink ) {
      while ( upper.length > 1 && isConvex(upper[upper.length-2],upper[upper.length-1],v) )
        faceSink( [upper[upper.length-2], upper.pop(), v] );

      for (var i=0, vl1, vl2; vl1 = lower[i], vl2=lower[i+1]; i++)
        faceSink( [vl1, v, vl2] );

      upper.push(v);
      lower = lower.slice(lower.length-1,lower.length);
    }

    function addLower( v, faceSink ) {
      while ( lower.length > 1 && isConvex(v,lower[lower.length-1],lower[lower.length-2]) )
        faceSink( [lower[lower.length-2], lower.pop(), v] );

      for (var i=0, vl1, vl2; vl1 = upper[i], vl2=upper[i+1]; i++)
        faceSink( [vl2, v, vl1] );

      lower.push(v);
      upper = upper.slice(upper.length-1,upper.length);
    }

    function attemptAdd( v, vertices, faceSink ) {
      var vu1 = upper[upper.length-1];
      var vu2 = loop( vertices, vu1.index+1 );
      if ( v === vu2 ) {
        addUpper( v, faceSink );
        return "HANDLED";
      }

      var vl1 = lower[upper.length-1] || vu1;
      var vl2 = loop( vertices, vl1.index - 1 );
      if ( v === vl2 ) {
        addLower( v, faceSink );
        return "HANDLED";
      }

      if ( vertexIsAboveLine(v, vu1, vu2 ) )
        return "ABOVE";
      else if ( vertexIsAboveLine(v, vl1, vl2 ) )
        return "INSIDE";
      else
        return "BELOW";
    }

    function split( v, faceSink ) {
      var other = MonotonePolygon( lower );
      addUpper( v, faceSink );

      lower = [];
      addLower( v, faceSink );

      return other;
    }

    return { addUpper:addUpper, addlower:addLower, attemptAdd:attemptAdd, split:split };
  }


  function tesselate( points3d, faceSink ) {
    // determine approximate plane through loop and map points onto it.
    var planeNormal = vnorm(loopMeanNormal(points3d));
    var planePoint = vectorAverage(points3d);

    var offNormal = vadd( planePoint, Math.abs(planeNormal[0] < .5) ? [1,0,0] : [0,1,0] );
    var axis1 = vnorm(projectIntoPlane(offNormal ,planePoint, planeNormal));
    var axis2 = vcross( planeNormal, axis1 );

    var vertices2d = points3d.map( function(p,index) {
      var pointInPlane = projectIntoPlane(p, planePoint, planeNormal);
      return new Vertex([vdot(pointInPlane,axis1), vdot(pointInPlane,axis2),0],0,0,index);
    })

    // tesselate 2d loop
    tesselateVertices( vertices2d, function(face) {
      faceSink([ points3d[face[0].id], points3d[face[1].id], points3d[face[2].id] ]);
    })
  }


  /*
   Better triangulation algorithm:
   Monotone polygon - a polygon for which no vertical line passes through it more than twice.
      Implemented as an upper and lower list of vertices.
      Assume vertices are added from left to right.
      addUpper v
        while lower.length > 1
          triangle(lower[0],v,lower[1])
          lower = lower.slice(1)
        if upper.length == 1 && angle(lower[0],upper[0],v) is convex
          triangle(lower[0],upper[0],p)
          upper = []
        upper.push(v)

      addLower v - similar but opposite

      attemptAdd v
        lastUpper = upper[upper.length-1]
        nextUpper = vertices[lastUpper.index + 1]
        lastUpper = upper[upper.length-1]
        nextLower = vertices[lastLower.index + 1]
        if v.index == nextUpper
          addUpper v
          return HANDLED
        else if v.index == nextLower
          addLower v
          return HANDLED
        else if aboveLine(v,lastUpper,nextUpper)
          return ABOVE
        else if belowLine(v,lastLower,nextLower)
          return BELOW
        return INSIDE

      split v
        next = MonotonePolygon()
        for (vl in lower) next.addLower vl
        next.addUpper v
        lower = []
        addLower v
        return next


    Triangulation algorithm:
      Sort vertices from left to right.
      monotonePolygons = []
      for v in vertices
        for poly in monotonePolygons
          position = poly.attemptAdd v
          if position != BELOW
            break
        if position == ABOVE
          insert new monontonePolygon above
        if position == BELOW
          insert new monotonePolygon below
        if position == INSIDE
          insert poly.split(v) below

   */

  function tesselateVertices(vertices, faceSink) {
    console.log('tesslate', vertices.map(function(v) { return [v.x,v.y].join(',')}))
    if ( vertices.length < 3 ) return;
    if ( vertices.length < 4 ) return tesselateConvex(vertices, faceSink);
      
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
        console.log('start end', startPoint,endPoint)

        var newShape, rest;
        if (endPoint < startPoint) {
          newShape = vertices.slice(startPoint).concat(vertices.slice(0,endPoint+1));
          rest = vertices.slice(endPoint, startPoint+1);
        } else {
          newShape = vertices.slice(startPoint, endPoint+1);
          rest = vertices.slice(endPoint).concat(vertices.slice(0,startPoint+1));
        }

        if ( ! isSelfIntersecting( vertices[startPoint], vertices[endPoint], rest.slice(1,rest.length-1) ) ) {
          console.log('slice', [vertices[startPoint], vertices[endPoint]].map(function(v) { return [v.x,v.y].join(',')}))
          if ( Math.abs(endPoint- startPoint) <= 1 || Math.abs(endPoint- startPoint) >= vertices.length -1 ) {
            console.log("wrap around", startPoint, endPoint, vertices.length);
            return;
          }
          if (endPoint < startPoint) {
            tesselateVertices( vertices.slice(startPoint).concat(vertices.slice(0,endPoint+1)), faceSink );
            tesselateVertices( vertices.slice(endPoint, startPoint+1 ), faceSink );
          } else {
            tesselateVertices( vertices.slice(startPoint, endPoint+1 ), faceSink );
            tesselateVertices( vertices.slice(endPoint).concat(vertices.slice(0,startPoint+1)), faceSink );
          }

          return;
        }
      }
    }
    
    return tesselateConvex(vertices, faceSink);
    
  }

  function tesselateConvex(vertices, faceSink) {
    console.log('tesselateConvex', vertices.map(function(v) { return [v.x,v.y].join(',')}))
    var ai = 0, bi=(vertices.length/3)|0, ci=(2*vertices.length/3)|0;
    var a = vertices[ai], b=vertices[bi], c=vertices[ci];

    faceSink( [c,b,a] );
    
    if (bi-ai > 1) tesselateConvex(vertices.slice(ai,bi+1), faceSink);
    if (ci-bi > 1) tesselateConvex(vertices.slice(bi,ci+1), faceSink);
    if (vertices.length-ci > 1) 
      tesselateConvex(vertices.slice(ci,vertices.length).concat([a]), faceSink);
  }

  var all = {
      vadd:vadd, vsub:vsub, vscale:vscale, vdot:vdot, vcross: vcross, vlength:vlength, vnorm:vnorm,
      step:step,
      Path:Path,
      Vertex:Vertex, vertices:vertices, parametric:parametric,
      translate:translate,
      skin:skin, facers:facers, closeEdge:closeEdge, capBottom:capBottom, capTop:capTop,
      reverse:reverse,
      ThreeJSRenderer:ThreeJSRenderer, STLRenderer:STLRenderer
  };
  for (var k in all) context[k] = all[k];

})(typeof window != 'undefined' ? (window.manifold={}) : exports);
