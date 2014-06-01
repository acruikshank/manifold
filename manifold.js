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
  var debug = false;
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
  var vertexId = 1;

  function Vertex( point, transformStep, ribStep, id) {
    this[0] = point[0];
    this[1] = point[1];
    this[2] = point[2];
    this[3] = transformStep;
    this[4] = ribStep;
    this[5] = id == null ? ++vertexId : id;
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

    path.vertices = function(numVertices, step) {
      return function(vertexSink) {
        var allDivisions = Math.max(numVertices, segments.length + 1);
        var remainingDivisions = allDivisions - 1 - segments.length;
        var remainingWeight = totalWeight;
        var index = 0;

        for (var i=0,s; s = segments[i]; i++) {
          var divisions = i >= segments.length - 1 ? remainingDivisions + 1 : (s.w * remainingDivisions / remainingWeight + 1)|0;
          for (var j=0; j<divisions; j++) {
            vertexSink( new Vertex(resolveCurve(points,j/divisions, s.s, s.o), step, (index++)/(allDivisions-1)) );
          }

          remainingDivisions -= divisions - 1;
          remainingWeight -= s.w;
        }

        vertexSink( new Vertex(points[points.length-1], step, 1) );
      }
    }
    
    return path;
  }

  // STEP
  function step(iterations) {
    return function(stepSink) {
      for (var i=0; i < iterations; i++) stepSink( i / (iterations - 1) );
    }
  }

  // VERTEX GENERATORS

  // Given a list of 3d points, generated a vertex for each point per step 
  // (points need to be transformed)
  function vertices(points) {
    return function(vertexSink) {
      return function(step) {
        for (var i=0,p; p=points[i]; i++)
          vertexSink( new Vertex(p, step, i / (points.length-1)) );
      }
    }
  }

  // Convert steps into vertices using a given function.
  function parametric(f, stepper) {
    return function(vertexSink) {
      return function( transformStep ) {
        stepper( function(ribStep) {
          vertexSink( new Vertex(f(ribStep, transformStep), transformStep, ribStep) );
        })
      }
    }
  }

  // Simplifies step generators
  function vertexGenerator( f ) {
    return function(vertexSink) {
      return function(step) {
        f(step,vertexSink)
      }
    }
  }

  // Convert a path into a manifold using a path generator function that converts vertices
  // on the path into new paths.
  function PathParameterized(path, transformSteps, ribSteps) {
    return function(generator) {
      return function(vertexSink) {
        path.vertices(transformSteps, 0)(function(vertex) {
          generator(vertex).vertices(ribSteps, vertex.ribStep)(vertexSink);
        });
      }
    }
  }

  // Create a vertext generator that converts a single vertex into a circle of vertices
  // such that the circle passes through the vertex and its center lies on the center normal.
  function CircleRib( steps, centerNormal, centerOffset ) {
    centerNormal = vnorm(centerNormal);
    centerOffset = centerOffset || [0,0,0];
    centerOffset = vsub(centerOffset, vscale(centerNormal, vdot(centerOffset,centerNormal)));
    return function(vertexSink) {
      return function(vertex) {
        var center = vadd(vscale(centerNormal, vdot(centerNormal,vertex)), centerOffset);        
        var a = vsub(vertex,center);
        var radius = vlength(a);
        if (radius === 0)
          return vertexSink( new Vertex(center, vertex.ribStep, 1) );

        var b = vcross(centerNormal,a);
        for (var i=0; i<steps; i++) {
          var point = vadd(vadd(vscale(a,Math.cos(2*i*Math.PI/(steps))),vscale(b,-Math.sin(2*i*Math.PI/(steps)))),center);
          vertexSink( new Vertex(point, vertex.ribStep, i/(steps-1) ) );
        }
      }
    }
  }

  // TRANSFORM
  function translate(translations) {
    return function( vertexSink ) {
      return function( vertex ) {
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
          faceSink( [blVertex, nextRib[nextRib.length-1], vertex] );
        }

        while ( brVertex &&  vertex.ribStep > blVertex.ribStep + (brVertex.ribStep-blVertex.ribStep)/2 ) {
          faceSink( [blVertex, vertex, brVertex] );
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
      if ( vertex.transformStep === 0 )
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
    geometry._manifoldIdMap = {};
    return {
      renderer : function( face ) {
        var index0 = saveThreeJSVertex( geometry, face[0] );
        var index1 = saveThreeJSVertex( geometry, face[1] );
        var index2 = saveThreeJSVertex( geometry, face[2] );
        geometry.faces.push( new THREE.Face3(index0, index1, index2) );
      },
      geometry: geometry
    }
  }

  function CSGRenderer() {
    var polygons = [];
    return {
      renderer : function( face ) {
        var faceNormal = vnorm(vcross(vsub(face[0],face[1]),vsub(face[2],face[1])));
        var vertices = [ 
          new CSG.Vertex( face[0], faceNormal ), 
          new CSG.Vertex( face[1], faceNormal ), 
          new CSG.Vertex( face[2], faceNormal )];
        polygons.push( new CSG.Polygon( vertices ) );        
      },
      csgObject: function() { return CSG.fromPolygons(polygons); }
    }
  }

  function saveThreeJSVertex(geometry, vertex) {
    var location = geometry._manifoldIdMap[vertex.id];
    if ( location == null ) {
      location = geometry._manifoldIdMap[vertex.id] = geometry.vertices.length;
      geometry.vertices.push( new THREE.Vector3( vertex[0], vertex[1], vertex[2] ) );
    }
    return location;
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

  function printVertex(v) { return '('+v.id+':'+v.x.toFixed(2)+','+v.y.toFixed(2)+','+v.z.toFixed(2)+')'}

  function MonotonePolygon( first, lowerVertices ) {
    var upper = first ? [first] : [];    
    var lower = lowerVertices || [];
    var mergePolygon;

    function addUpper( v, faceSink ) {
      while ( upper.length > 1 && isConvex(upper[upper.length-2],upper[upper.length-1],v) ) {
        faceSink( [upper[upper.length-1], v, upper[upper.length-2]] );
        upper.pop();
      }

      if (upper.length == 1 && lower.length && isConvex(lower[0],upper[0],v))
        faceSink( [lower[0], upper.pop(), v] );

      for (var i=0, vl1, vl2; vl1 = lower[i], vl2=lower[i+1]; i++)
        faceSink( [vl1, v, vl2] );

      upper.push(v);
      lower = lower.slice(lower.length-1,lower.length);
    }

    function addLower( v, faceSink ) {
      while ( lower.length > 1 && isConvex(v,lower[lower.length-1],lower[lower.length-2]) ) {
        faceSink( [v, lower[lower.length-1], lower[lower.length-2]] );
        lower.pop();
      }

      if (upper.length && lower.length == 1 && isConvex(v,lower[0],upper[0]))
        faceSink( [lower.pop(), upper[0], v] );

      for (var i=0, vl1, vl2; vl1 = upper[i], vl2=upper[i+1]; i++)
        faceSink( [vl2, v, vl1] );

      lower.push(v);
      upper = upper.slice(upper.length-1,upper.length);
    }


    function lastLower() {
      return lower[lower.length-1] || upper[0]
    }

    function aboveBottom(v, vertices) {
      var vl1 = lastLower();
      var vl2 = loop( vertices, vl1.id - 1 );
      return vertexIsAboveLine(v, vl1, vl2 )
    }

    function isBottom(v, vertices) {
      return v === loop( vertices, lastLower().id - 1 );
    }

    function getUpper() { return upper; }
    function getLower() { return lower; }

    function attemptAdd( v, vertices, faceSink ) {
      if (debug)
        console.log('vertex',v, upper.map(printVertex).join(','), lower.map(printVertex).join(','));
      var vu1 = upper[upper.length-1];
      if ( ! vu1 ) {
        addUpper( v, faceSink );
        return "TOP";
      }

      var vu2 = loop( vertices, vu1.id+1 );
      if ( v === vu2 ) {
        addUpper( v, faceSink );
        if (mergePolygon) {
          mergePolygon.addUpper(v, faceSink);
          upper = mergePolygon.getUpper();
          lower = mergePolygon.getLower();
          mergePolygon = null;
        }
        return isBottom(v,vertices) ? "DONE" : "TOP";
      }

      if ( mergePolygon ? mergePolygon.isBottom(v,vertices) : isBottom(v,vertices) ) {
        addLower( v, faceSink );
        if (mergePolygon) {
          mergePolygon.addLower( v, faceSink );
          mergePolygon = null;
        }
        var vl1 = lastLower();
        var vl2 = loop( vertices, vl1.id - 1 );
        return vl2.x >= vl1.x ? "BOTTOM" : "MERGE";
      }

      if ( vertexIsAboveLine(v, vu1, vu2 ) )
        return "ABOVE";
      else if ( mergePolygon ? mergePolygon.aboveBottom(v,vertices) : aboveBottom(v,vertices) )
        return "INSIDE";
      else
        return "BELOW";
    }

    function merge( polygon ) {
      mergePolygon = polygon;
    }

    function split( v, faceSink ) {
      var other = mergePolygon || MonotonePolygon( null, lower.length ? [lastLower()] : [upper[0]] );
      mergePolygon = null;

      if (debug) {
        console.log('SPLIT', upper.map(printVertex).join(','), lower.map(printVertex).join(','));
        console.log('OTHER', other.getUpper().map(printVertex).join(','), other.getLower().map(printVertex).join(','))
      }
      other.addUpper( v, faceSink );
      addLower( v, faceSink );

      return other;
    }

    return { addUpper:addUpper, addLower:addLower, attemptAdd:attemptAdd, merge:merge, split:split,
             isBottom:isBottom, aboveBottom:aboveBottom, getUpper:getUpper, getLower:getLower };
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
    tesselate2d( vertices2d, function(face) {
      faceSink([ points3d[face[2].id], points3d[face[1].id], points3d[face[0].id] ]);
    })
  }

  function tesselate2d( vertices, faceSink ) {
    var sortedVertices = vertices.slice().sort(function(a,b) { return a.x - b.x || b.y - a.y; });
    var monotones = [MonotonePolygon(sortedVertices[0])];

    for (var i=1, v; v = sortedVertices[i]; i++) {
      var result;
      for (var j=0,mp; mp = monotones[j]; j++) {
        result = mp.attemptAdd(v, vertices, faceSink);
        if (result !== 'BELOW') break;
      }
      if (debug)
        console.log( "TESSELATE", j, result )

      if ( result === 'ABOVE' )
        monotones.splice( j, 0, MonotonePolygon(v) );
      else if ( result === 'INSIDE' )
        monotones.splice( j+1, 0, mp.split(v, faceSink) );
      else if ( result === 'BELOW' )
        monotones.push( MonotonePolygon(v) );
      else if ( result === 'MERGE' ) {
        monotones[j+1].attemptAdd( v, vertices, faceSink );
        monotones[j].merge(monotones[j+1])
        monotones.splice(j+1,1)
      } else if (result === 'DONE' )
        monotones.splice(j, 1);
    } 
  }



  var all = {
      vadd:vadd, vsub:vsub, vscale:vscale, vdot:vdot, vcross: vcross, vlength:vlength, vnorm:vnorm,
      step:step,
      Path:Path, PathParameterized:PathParameterized, CircleRib:CircleRib,
      Vertex:Vertex, vertices:vertices, parametric:parametric, vertexGenerator:vertexGenerator,
      translate:translate,MonotonePolygon:MonotonePolygon, tesselate2d:tesselate2d,
      skin:skin, facers:facers, closeEdge:closeEdge, capBottom:capBottom, capTop:capTop,
      reverse:reverse,
      ThreeJSRenderer:ThreeJSRenderer, CSGRenderer:CSGRenderer, STLRenderer:STLRenderer
  };
  for (var k in all) context[k] = all[k];

})(typeof window != 'undefined' ? (window.manifold={}) : exports);
