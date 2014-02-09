var expect = require('expect.js')
var M = require('../manifold')

function prettyVertex(v) { return '{x:'+(v[0]).toFixed(2)+',y:'+(v[1]).toFixed(2)+',z:'+(v[2]).toFixed(2)+'}'; }
function prettyFace(f) { return 'face:' + f.map(prettyVertex).join(', '); }
function faceNormal(face) { return M.vnorm( M.vcross(M.vsub(face[1],face[0]), M.vsub(face[2],face[0])) ) }

function faceContainsPoint( face, p ) {
  for (var i=0,v; v = face[i]; i++ )
    if ( M.vlength(M.vsub(p,v)) < .00001) return true;
  return false;
}

// create vertex from each of a 2d array of coordinates using the index for rib and transform steps
// and feed it into facer.
function applyVertices( facer, vertices ) {
  vertices.forEach(function(rib, transformStep) { rib.forEach(function(v, ribStep) {
    var ts = transformStep/(vertices.length-1), rs = rib.length > 1 ? ribStep/(rib.length - 1) : 1;
    facer( new M.Vertex(v, ts, rs) );
  }) })
}

function facesContainsFace( faces, points ) {
  for (var i=0,face; face = faces[i]; i++)
    if (points.reduce(function(m,p) { return m && faceContainsPoint(face,p) }, true))
      return true;
  return false;
}

expect.Assertion.prototype.containPoint = function( p ) {
  this.assert( faceContainsPoint(this.obj, p), 
      function() { return 'expected face to contain ' + prettyVertex(p) }, 
      function() { return 'expected face to not contain ' + prettyVertex(p) })
}

expect.Assertion.prototype.containFace = function( points ) {
  this.assert( facesContainsFace(this.obj, points), 
      function() { return 'expected faces to contain ' + prettyFace(points) }, 
      function() { return 'expected faces to not contain ' + prettyFace(points) })
}

var DEFAULT_DISTANCE = .001;
expect.Assertion.prototype.nearTo = function( vec, difference ) {
  if (difference == undefined) difference = DEFAULT_DISTANCE;

  var expectStr = Array.prototype.map.call(this.obj, function(x) {return x.toFixed(2)})
  var vecStr = vec.map(function(x) {return x.toFixed(2)})
  this.assert(
    M.vlength(M.vsub(this.obj, vec)) <= difference
    , function(){ return 'expected ' + expectStr + ' to be close to ' + vecStr }
    , function(){ return 'expected ' + expectStr + ' to be far from ' + vecSter }
  )
}

describe( 'math:', function() {
  it('adds correctly', function() {
    expect( M.vadd([1,2,3], [4,-3,2]) ).to.be.nearTo([5,-1,5]);
    expect( M.vadd([-.1,.4,-.6], [.01,.3,.5]) ).to.be.nearTo([-.09,.7,-.1]);
  })

  it('subtracts correctly', function() {
    expect( M.vsub([1,2,3], [4,-3,2]) ).to.be.nearTo([-3,5,1]);
    expect( M.vsub([-.1,.4,-.6], [.01,.3,.5]) ).to.be.nearTo([-.11,.1,-1.1]);    
  })

  it('dot products correctly', function() {
    expect( M.vdot([1,2,3], [4,-3,2]) ).to.equal(4);
    expect( M.vdot([-.1,.4,-.6], [.01,.3,.5]) ).to.equal(-.181);
  })

  it('crosses axis normals', function() {
    expect( M.vcross([1,0,0], [0,1,0]) ).to.be.nearTo( [0,0,1] )
    expect( M.vcross([0,1,0], [0,0,1]) ).to.be.nearTo( [1,0,0] )
    expect( M.vcross([1,0,0], [0,0,1]) ).to.be.nearTo( [0,-1,0] )
  })

  it('cross product is orthogonal to both vectors', function() {
    var a = [3.2389, 4.23342, 984394];
    var b = [-22.3829, 6.34834, 6.43834];
    var c = M.vcross(a,b);
    expect( M.vdot(a,c) ).to.be.within(-.0001, .0001);
    expect( M.vdot(b,c) ).to.be.within(-.0001, .0001);
  })

  it('The length of a normal is one', function() {
    var a = [3.2389, 4.23342, 9.84394];
    var b = [-22.3829, 6.34834, 6.43834];
    expect( M.vlength( M.vnorm(a) ) ).to.within(.999999, 1.000001)    
    expect( M.vlength( M.vnorm(b) ) ).to.within(.999999, 1.000001)    
  })

  it('norm(a) * |a| == a', function() {
    var a = [3.2389, 4.23342, 9.84394];
    var aNorm = M.vnorm(a)
    var b = [-22.3829, 6.34834, 6.43834];
    var bNorm = M.vnorm(b)
    expect( M.vscale(aNorm, M.vlength(a)) ).to.be.nearTo(a);
    expect( M.vscale(bNorm, M.vlength(b)) ).to.be.nearTo(b);
  })
})

describe( 'vertices:', function() {
  it('emits a rib of three vertices when given 3 points', function() {
    var vertices = [];    
    M.vertices([[17,11,5],[7,13,3],[19,2,23]])(function(vertex) { vertices.push(vertex); })(0)

    expect(vertices.length).to.equal(3);
    expect(vertices[0].x).to.equal(17);
    expect(vertices[0].y).to.equal(11);
    expect(vertices[0].z).to.equal(5);

    expect(vertices[1].x).to.equal(7);
    expect(vertices[1].y).to.equal(13);
    expect(vertices[1].z).to.equal(3);

    expect(vertices[2].x).to.equal(19);
    expect(vertices[2].y).to.equal(2);
    expect(vertices[2].z).to.equal(23);
  })
})

describe( 'skin:', function() {
  it('emits two faces when given a square on two ribs', function() {
    var faces = [];
    var facer = M.skin(function(face) { faces.push(face); });
    facer( new M.Vertex([0,0,0],0,0) )
    facer( new M.Vertex([1,0,0],0,1) )
    facer( new M.Vertex([0,0,1],1,0) )
    facer( new M.Vertex([1,0,1],1,1) )

    expect(faces.length).to.equal(2);
    expect(faces[0]).to.containPoint( [0,0,0] );
    expect(faces[0]).to.containPoint( [1,0,1] );
    expect(faces[0]).to.containPoint( [0,0,1] );

    expect(faces[1]).to.containPoint( [0,0,0] );
    expect(faces[1]).to.containPoint( [1,0,1] );
    expect(faces[1]).to.containPoint( [1,0,0] );
  })

  it('emits 4 faces for 3 ribs of 3,1 and 3 vertices', function() {
    var faces = [];
    var facer = M.skin(function(face) { faces.push(face); });
    facer( new M.Vertex([-1,0,-1],0,0) )
    facer( new M.Vertex([0,0,-1],0,.5) )
    facer( new M.Vertex([1,0,-1],0,1) )
    facer( new M.Vertex([0,0,0],.5,1) )
    facer( new M.Vertex([-1,0,1],1,0) )
    facer( new M.Vertex([0,0,1],1,.5) )
    facer( new M.Vertex([1,0,1],1,1) )

    expect(faces.length).to.equal(4);
  })
})

describe( 'reverse:', function() { 
  it('reverses face normal', function() {
    var vertices = [new M.Vertex([0,0,0],0,0), new M.Vertex([1,0,0],0,1), new M.Vertex([0,0,1],1,1)]
    var face;
    var facer = M.skin(function(f) { face = f; });
    for (var i=0,v; v=vertices[i]; i++) facer( v );

    // compute face normal
    var forward = faceNormal(face);

    face = null;
    facer = M.reverse(M.skin) (function(f) { face = f; } );
    for (var i=0,v; v=vertices[i]; i++) facer( v );

    // compute face normal
    var reverse = faceNormal(face);

    expect( M.vdot(forward,reverse) ).to.be.within( -1.00001, -.99999 );
  }) 
})

describe( 'closeEdge', function() {
  var faces;

  describe("with a simple cube", function() {
    beforeEach( function() {
      faces = [];
      var vertices = [  // a pyramid atop a cube atop an inverted pyramid.
        [[1,1,2],[1,-1,2],[-1,-1,2],[-1,1,2]],
        [[1,1,4],[1,-1,4],[-1,-1,4],[-1,1,4]]];

      var faceSink = function(face) { faces.push(face); };
      var facer = M.facers( M.skin, M.closeEdge )( faceSink );

      applyVertices( facer, vertices );
    })

    it ('produces all skin faces plus 2 faces to close', function() {
      expect( faces.length ).to.equal( 8 );
    })

    it ('closes the cube', function() {
      expect( faces ).to.containFace( [[-1,1,2],[-1,1,4],[1,1,4]] )
      expect( faces ).to.containFace( [[-1,1,2],[1,1,4],[-1,1,2]] )
    })

  })

  describe("with degenerate ribs", function() { 
    beforeEach( function() {
      faces = [];
      var vertices = [  // a pyramid atop a cube atop an inverted pyramid.
        [[0,0,0]],
        [[1,1,2],[1,-1,2],[-1,-1,2],[-1,1,2]],
        [[1,1,4],[1,-1,4],[-1,-1,4],[-1,1,4]],
        [[0,0,6]]];

      var faceSink = function(face) { faces.push(face); };
      var facer = M.facers( M.skin, M.closeEdge )( faceSink );

      applyVertices( facer, vertices );
    })

    it ('produces all skin faces plus 4 faces to close', function() {
      expect( faces.length ).to.equal( 16 );
    })

    it ('closes the bottom pyramid', function() {
      expect( faces ).to.containFace( [[0,0,0],[-1,1,2],[1,1,2]] )
    })

    it ('closes the middle cube', function() {
      expect( faces ).to.containFace( [[-1,1,2],[-1,1,4],[1,1,4]] )
      expect( faces ).to.containFace( [[-1,1,2],[1,1,4],[-1,1,2]] )
    })

    it ('closes the top pyramid', function() {
      expect( faces ).to.containFace( [[0,0,6],[-1,1,4],[1,1,4]] )
    })
  })
})

describe('parametric', function() {
  var vertices;

  it ('emits vertices', function() {
    vertices = [];
    function linear(r,t) { return [r,t,5] }
    M.step(2) (M.parametric(linear, M.step(2)) (function(v) {vertices.push(v)}))

    expect( vertices.length ).to.be(4);
    expect( vertices[0] ).to.be.nearTo([0,0,5])
    expect( vertices[1] ).to.be.nearTo([1,0,5])
    expect( vertices[2] ).to.be.nearTo([0,1,5])
    expect( vertices[3] ).to.be.nearTo([1,1,5])
  })
})

describe( 'cap', function() {
  var faces;
  var vertices = [  // a pyramid atop a cube atop an inverted pyramid.
    [[1,1,2],[1,-1,2],[-1,-1,2],[-1,1,2]],
    [[1,1,4],[1,-1,4],[-1,-1,4],[-1,1,4]]];

  describe("bottom", function() {
    beforeEach( function() {
      faces = [];

      var faceSink = function(face) { faces.push(face); };
      var facer = M.capBottom( faceSink );

      applyVertices( facer, vertices );
    })

    it ('produces 2 vertices for the bottom', function() {
      expect( faces.length ).to.equal( 2 );
    })

    it ('faces point down', function() {
      expect( faceNormal(faces[0]) ).to.be.nearTo( [0,0,-1] )
      expect( faceNormal(faces[1]) ).to.be.nearTo( [0,0,-1] )
    })
  })

  describe("top", function() {
    beforeEach( function() {
      faces = [];

      var faceSink = function(face) { faces.push(face); };
      var facer = M.capTop( faceSink );

      applyVertices( facer, vertices );
    })

    it ('produces 2 vertices for the top', function() {
      expect( faces.length ).to.equal( 2 );
    })

    it ('faces point up', function() {
      expect( faceNormal(faces[0]) ).to.be.nearTo( [0,0,1] )
      expect( faceNormal(faces[1]) ).to.be.nearTo( [0,0,1] )
    })

  })
})

