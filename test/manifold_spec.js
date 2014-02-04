var expect = require('expect.js')
var M = require('../manifold')


var DEFAULT_DISTANCE = .001;
expect.Assertion.prototype.nearTo = function( vec, difference ) {
  if (difference == undefined) difference = DEFAULT_DISTANCE;

  var expectStr = this.obj.map(function(x) {return x.toFixed(2)})
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

describe( 'facer:', function() {
  it('emits two faces when given a square on two ribs', function() {
    var faces = [];
    var facer = M.facer()(function(face) { faces.push(face); });
    facer( new M.Vertex([0,0,0],0,0) )
    facer( new M.Vertex([1,0,0],0,1) )
    facer( new M.Vertex([0,0,1],1,0) )
    facer( new M.Vertex([1,0,1],1,1) )

    expect(faces.length).to.equal(2);
  })
})
