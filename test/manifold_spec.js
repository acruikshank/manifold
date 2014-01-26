expect = require('expect.js')
manifold = require('../manifold')

var M = manifold.math;

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
})
