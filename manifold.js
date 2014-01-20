/*
Step: [0-1]
StepSink: Step ->
StepSource: StepSink ->
  step(n): stepSource // produce n linear steps between 0 and 1 inclusive
VertexSink: Vertex, Step, Step* ->
VertexGenerator: VertexSink -> StepSink
  vertices(vertices) // take the given vertices to create a VertexGenerator
  vertices(f): (Step -> Vertex) -> VertexGenerator // convert a vertex producer into a vertex generator 
Transformer: VertexSink -> VertexSink
  zTranslate(start,end) -> translate vertices along z axis
  zRotate(start*,end*) -> rotate vertices around z axis
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
  vertices([[1,1],[-1,1],[-1,-1],[1,-1]])(
    zTranslate(1)(
      face({connect:'lr',tesselate:'tb'})( renderer )))))

var geometry = renderer.geometry;

// sphere
var R = 40, Q = 20, 2PI = 2*Math.Pi, hPI = Math.PI/2;
function simiCircle(s) {return [R*Math.cos(tween(s,-hPI,hPI)),0,R*Math.sin(tween(s,-hPI,hPI))]}

step(Q)(
  step(Q)(
    vertices(semiCircle))(
      zRotate()(
        face(connect:'tb', singular:'lr')( renderer ))))

