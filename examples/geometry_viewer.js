function geometryViewer(geometrySource) {
  var camera, scene, renderer;
  var particleLight1, pointLight1;
  var particleLight2, pointLight2;
  var meshes = {};
  var dragging = false, startMouse;

  var M = manifold;

  initScene();

  function onWindowResize() {
    camera.aspect = window.innerWidth / window.innerHeight;
    camera.updateProjectionMatrix();

    renderer.setSize( window.innerWidth, window.innerHeight );
  }

  var threeGeometry = M.ThreeJSRenderer();
  var geometry = threeGeometry.geometry;

  geometrySource(threeGeometry.renderer);

  geometry.computeBoundingSphere();
  geometry.computeFaceNormals();


  mesh = new THREE.Mesh(geometry, new THREE.MeshPhongMaterial( { color: '#ffffff' } ) );
  scene.add(mesh);

  function render( ) {
    renderer.render( scene, camera );
  }

  render();

  function initScene() {
    renderer = new THREE.WebGLRenderer();
    renderer.setSize( window.innerWidth, window.innerHeight );
    renderer.shadowMapEnabled = true;
    document.body.appendChild( renderer.domElement );

    camera = new THREE.PerspectiveCamera( 70, window.innerWidth / window.innerHeight, 1, 500 );
    camera.position.y = -60;
    camera.position.z = 60;
    camera.lookAt( new THREE.Vector3(0,0,0))


    scene = new THREE.Scene();

    scene.fog = new THREE.Fog( 0xffffff, 1, 1000 );
    scene.fog.color.setHSL( 0.6, 0, .2 );

    var light = new THREE.DirectionalLight(0xbbbbbb);
    light.position.set(50,-50,100);
    light.rotation.x = light.rotation.y = light.rotation.z = 0;
    light.target.position.set(0.0,0.0,0.0);
    light.target.updateMatrixWorld();
    light.castShadow = true;
    light.shadowMapWidth = 4096;
    light.shadowMapHeight = 4096;

    var d = 200;
    light.shadowCameraLeft = -d;
    light.shadowCameraRight = d;
    light.shadowCameraTop = d;
    light.shadowCameraBottom = -d;

    scene.add(light);

    var hemiLight = new THREE.HemisphereLight( 0xffffff, 0xffffff, 0.6 );
    hemiLight.color.setHSL( 0.095, 0, 0.6 );
    hemiLight.groundColor.setHSL( 0.095, .05, 0.8 );
    hemiLight.position.set( 0, 500, 0 );
    scene.add( hemiLight );


    // GROUND
    var groundGeo = new THREE.PlaneGeometry( 10000, 10000 );
    var groundMat = new THREE.MeshPhongMaterial( { ambient: 0x050505, color: 0x050505, specular: 0x050505 } );
    groundMat.color.setHSL( 0.095, .05, 0.8 );

    ground = new THREE.Mesh( groundGeo, groundMat );
    ground.position.z = -20;
    scene.add( ground );

    ground.receiveShadow = true;

    renderer.shadowMapEnabled = true;
    renderer.shadowMapSoft = true;
    renderer.shadowMapCullFace = THREE.CullFaceBack;

    renderer.render( scene, camera );


    window.addEventListener( 'resize', onWindowResize, false );
  }

  function mousemove(e) {
    if (! dragging) return;
    e.preventDefault();

    var xRotation = -(e.clientY) * Math.PI * 2 / document.body.offsetHeight;
    var yRotation = (e.clientX) * Math.PI * 2 / document.body.offsetWidth;

    scene.remove(mesh);
    mesh = new THREE.Mesh(geometry, new THREE.MeshPhongMaterial( { color: '#ffffff' } ) );
    scene.add(mesh);

    mesh.rotateOnAxis(new THREE.Vector3(0,1,0), yRotation);
    mesh.rotateOnAxis(new THREE.Vector3(1,0,0), xRotation);

    render();
    startMouse = {x:e.clientX,y:e.clientY};
  }
  document.body.onmousemove = mousemove;
  document.body.onmousedown = function(e) { e.preventDefault(); startMouse = {x:e.clientX,y:e.clientY}; dragging = true; }
  document.body.onmouseup = function(e) { mousemove(e); dragging = false; }
}