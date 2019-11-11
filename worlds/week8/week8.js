"use strict"


////////////////////////////// USEFUL VECTOR OPERATIONS

let dot = (a, b) => {
   let value = 0;
   for (let i = 0 ; i < a.length ; i++)
      value += a[i] * b[i];
   return value;
}

let subtract = (a,b) => {
   let c = [];
   for (let i = 0 ; i < a.length ; i++)
      c.push(a[i] - b[i]);
   return c;
}

let normalize = a => {
   let s = Math.sqrt(dot(a, a)), b = [];
   for (let i = 0 ; i < a.length ; i++)
      b.push(a[i] / s);
   return b;
}

let cross = (a, b) => [ a[1] * b[2] - a[2] * b[1],
                        a[2] * b[0] - a[0] * b[2],
                        a[0] * b[1] - a[1] * b[0] ];
                        
////////////////////////////// MATRIX SUPPORT

let cos = t => Math.cos(t);
let sin = t => Math.sin(t);
let identity = ()       => [1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1];
let rotateX = t         => [1,0,0,0, 0,cos(t),sin(t),0, 0,-sin(t),cos(t),0, 0,0,0,1];
let rotateY = t         => [cos(t),0,-sin(t),0, 0,1,0,0, sin(t),0,cos(t),0, 0,0,0,1];
let rotateZ = t         => [cos(t),sin(t),0,0, -sin(t),cos(t),0,0, 0,0,1,0, 0,0,0,1];
let scale = (x,y,z)     => [x,0,0,0, 0,y,0,0, 0,0,z,0, 0,0,0,1];
let translate = (x,y,z) => [1,0,0,0, 0,1,0,0, 0,0,1,0, x,y,z,1];
let inverse = src => {
  let dst = [], det = 0, cofactor = (c, r) => {
     let s = (i, j) => src[c+i & 3 | (r+j & 3) << 2];
     return (c+r & 1 ? -1 : 1) * ( (s(1,1) * (s(2,2) * s(3,3) - s(3,2) * s(2,3)))
                                 - (s(2,1) * (s(1,2) * s(3,3) - s(3,2) * s(1,3)))
                                 + (s(3,1) * (s(1,2) * s(2,3) - s(2,2) * s(1,3))) );
  }
  for (let n = 0 ; n < 16 ; n++) dst.push(cofactor(n >> 2, n & 3));
  for (let n = 0 ; n <  4 ; n++) det += src[n] * dst[n << 2];
  for (let n = 0 ; n < 16 ; n++) dst[n] /= det;
  return dst;
}
let multiply = (a, b)   => {
   let c = [];
   for (let n = 0 ; n < 16 ; n++)
      c.push( a[n&3     ] * b[    n&12] +
              a[n&3 |  4] * b[1 | n&12] +
              a[n&3 |  8] * b[2 | n&12] +
              a[n&3 | 12] * b[3 | n&12] );
   return c;
}
let transpose = m => [ m[0],m[4],m[ 8],m[12],
                       m[1],m[5],m[ 9],m[13],
                       m[2],m[6],m[10],m[14],
                       m[3],m[7],m[11],m[15] ];

let transform = (m, v) => [
   m[0] * v[0] + m[4] * v[1] + m[ 8] * v[2] + m[12] * v[3],
   m[1] * v[0] + m[5] * v[1] + m[ 9] * v[2] + m[13] * v[3],
   m[2] * v[0] + m[6] * v[1] + m[10] * v[2] + m[14] * v[3],
   m[3] * v[0] + m[7] * v[1] + m[11] * v[2] + m[15] * v[3]
];

let Matrix = function() {
   let topIndex = 0,
       stack = [ identity() ],
       getVal = () => stack[topIndex],
       setVal = m => stack[topIndex] = m;

   this.identity  = ()      => setVal(identity());
   this.restore   = ()      => --topIndex;
   this.rotateX   = t       => setVal(multiply(getVal(), rotateX(t)));
   this.rotateY   = t       => setVal(multiply(getVal(), rotateY(t)));
   this.rotateZ   = t       => setVal(multiply(getVal(), rotateZ(t)));
   this.save      = ()      => stack[++topIndex] = stack[topIndex-1].slice();
   this.scale     = (x,y,z) => setVal(multiply(getVal(), scale(x,y,z)));
   this.translate = (x,y,z) => setVal(multiply(getVal(), translate(x,y,z)));
   this.value     = ()      => getVal();
}

////////////////////////////// SUPPORT FOR CREATING 3D SHAPES

const VERTEX_SIZE = 8;

let createCubeVertices = () => {
   let V = [], P = [ -1,-1, 1, 0,0, 1, 0,0,   1, 1, 1, 0,0, 1, 1,1,  -1, 1, 1, 0,1, 1, 0,1,
                      1, 1, 1, 0,0, 1, 1,1,  -1,-1, 1, 0,0, 1, 0,0,   1,-1, 1, 0,0, 1, 1,0,
                      1, 1,-1, 0,0,-1, 0,0,  -1,-1,-1, 0,0,-1, 1,1,  -1, 1,-1, 0,0,-1, 1,0,
                     -1,-1,-1, 0,0,-1, 1,1,   1, 1,-1, 0,0,-1, 0,0,   1,-1,-1, 0,0,-1, 0,1 ];
   for (let n = 0 ; n < 3 ; n++)
      for (let i = 0 ; i < P.length ; i += 8) {
         let p0 = [P[i],P[i+1],P[i+2]], p1 = [P[i+3],P[i+4],P[i+5]], uv = [P[i+6],P[i+7]];
	 V = V.concat(p0).concat(p1).concat(uv);
	 for (let j = 0 ; j < 3 ; j++) {
	    P[i   + j] = p0[(j+1) % 3];
	    P[i+3 + j] = p1[(j+1) % 3];
         }
      }
   return V;
}

let cubeVertices = createCubeVertices();


////////////////////////////// SUPPORT FOR SPLINES


let HermiteBasisMatrix = [
    2,-3, 0, 1,
   -2, 3, 0, 0,
    1,-2, 1, 0,
    1,-1, 0, 0
];

let BezierBasisMatrix = [
   -1,  3, -3,  1,
    3, -6,  3,  0,
   -3,  3,  0,  0,
    1,  0,  0,  0
];

let toCubicCurveCoefficients = (basisMatrix, M) => {
   let C = [];
   for (let i = 0 ; i < M.length ; i++)
      C.push(transform(basisMatrix, M[i]));
   return C;
}

let toCubicPatchCoefficients = (basisMatrix, M) => {
   let C = [];
   for (let i = 0 ; i < M.length ; i++)
      C.push(multiply(basisMatrix, multiply(M[i], transpose(basisMatrix))));
   return C;
}


////////////////////////////// SUPPORT FOR CREATING 3D SHAPES


// const VERTEX_SIZE = 8;    // EACH VERTEX IS: [ x,y,z, nx,ny,nz, u,v ]


// FUNCTION createMeshVertices() REPEATEDLY CALLS uvToShape(u, v, args).
// EACH CALL ADDS ANOTHER VERTEX TO THE MESH, IN THE FORM: [x,y,z, nx,ny,nz, u,v]

function createMeshVertices(M, N, uvToShape, arg) {

    // IMPLEMENTATION NOTES:

    // THIS IS ESSENTIALLY WHAT YOU HAVE ALREADY IMPLEMENTED.
    // THE ONLY SIGNIFICANT DIFFERENCE IS THAT YOU NEED TO PASS IN
    // arg AS A THIRD ARGUMENT WHEN YOU CALL uvToShape().

    // M column, N row
    if (M == 1 || N == 1) throw "Wrong column or row!";
    let vertices = [];
    let addVertex = (u, v) => {
        let a = uvToShape(u, v, arg);
        if (a){
        for (let i = 0 ; i < a.length ; i++)
            vertices.push(a[i]);}
    }

    let du = 1. / (M - 1);
    let dv = 1. / (N - 1);
    for (let row = 0; row < N - 1; row++) {
        let u0 = row % 2 == 0 ? 1 : 0;
        let sign = row % 2 == 0 ? -1 : 1;
        let vBot = row * dv;
        let vTop = (row + 1) * dv;
        if (row == 0) addVertex(u0, vBot);
        addVertex(u0, vTop);
        // let numSteps = M - 1;
        for (let col = 1; col < M; col++) {
            let u = u0 + sign * col * du;
            addVertex(u, vBot);
            addVertex(u, vTop);
        }
    }
    return vertices;
}

// FOR uvCubicCurvesRibbon(), arg IS IN THE BELOW FORM:
//
// {
//    width: width,
//    data: [
//       [ [a0x,b0x,c0x,d0x], [a0y,b0y,c0y,d0y], [a0z,b0z,c0z,d0z] ], // CURVE 0
//       [ [a1x,b1x,c1x,d1x], [a1y,b1y,c1y,d1y], [a1z,b1z,c1z,d1z] ], // CURVE 1
//       ...                                                         // ...
//    ]
// }


let uvToCubicCurvesRibbon = (u, v, arg) => {

    // IMPLEMENTATION NOTES:

    // THE MULTIPLE CURVES TOGETHER SPAN THE LENGTH OF THE RIBBON,
    // FROM u == 0.0 TO u == 1.0 FROM ONE END OF THE RIBBON TO THE OTHER.

    // arg.width SPECIFIES THE WIDTH OF THE RIBBON. THIS IS THE DIRECTION
    // THAT SPANS FROM v == 0.0 TO v == 1.0.

    // EACH ELEMENT OF arg.data PROVIDES CUBIC COEFFICIENTS FOR THE X,Y AND Z
    // COORDINATES, RESPECTIVELY, OF ONE CUBIC SEGMENT ALONG THE CURVE.

    // THE KEY TO IMPLEMENTATION IS TO EVAL THE CUBIC AT TWO SLIGHTLY
    // DIFFERENT VALUES OF u, SUCH AS u AND u+0.001.
    // THE DIFFERENCE VECTOR [dx,dy] IN X,Y CAN THEN BE USED TO
    // CREATE THE VECTOR THAT VARIES ALONG THE WIDTH OF THE RIBBON,
    // FROM p + (dy,-dx) WHEN v == 0.0 TO p + (-dy,dx) WHEN v == 1.0. 

    // VARYING IN Z IS TRICKY, BECAUSE YOU NEED TO FIGURE OUT HOW TO
    // COMPUTE A CORRECT VALUE FOR THE SURFACE NORMAL AT EACH VERTEX.
    // IF YOU CAN'T FIGURE OUT HOW TO PRODUCE A RIBBON THAT VARIES IN Z,
    // IT IS OK TO CREATE A RIBBON WHERE ALL THE Z VALUES ARE THE SAME.
    
    let numCurves = arg.data.length;
    let curveIndex = u == 1 ? numCurves - 1 : Math.floor(u * numCurves);
    let t = u * numCurves - curveIndex;
    let curveCoefficients = arg.data[curveIndex]; // [ [ax,bx,cx,dx], [ay,by,cy,dy], [az,bz,cz,dz] ]
    let calculateXYZ = (coef, t) => {
        let xyz = [];
        for (let i = 0; i < 3; i++) {
            xyz.push(coef[i][0]*t*t*t +
                     coef[i][1]*t*t +
                     coef[i][2]*t +
                     coef[i][3]);
        }
        return xyz;
    }
    let p0 = calculateXYZ(curveCoefficients, t);
    let p1 = calculateXYZ(curveCoefficients, t+0.001);
    let dp = subtract(p1, p0); // [dx, dy, 0]
    let vDirNorm = normalize([-dp[1], dp[0]]); // [-dy, dx]
    let width = arg.width;
    
    let x = p0[0] + vDirNorm[0] * (v - 0.5) * width;
    let y = p0[1] + vDirNorm[1] * (v - 0.5) * width;
    let z = p0[2];
    let nx = 0;
    let ny = 0;
    let nz = 1;
    return [x, y ,z, nx, ny, nz, u, v];
}


// For uvToCubicPatch, arg consists of bicubic coefficents in the form:
//
// [
//    [x0,x1, ... x15],  // Bicubic coefficients in x
//    [y0,y1, ... y15],  // Bicubic coefficients in y
//    [z0,z1, ... z15]   // Bicubic coefficients in z
// ]

let uvToCubicPatch = (u, v, arg) => {

    // IMPLEMENTATION NOTES:

    // THE THREE 4x4 MATRICES IN VARIABLE arg ARE VALUES OF Cx, Cy AND Cz.
    // THESE ARE THE BICUBIC COEFFICIENTS FOR X, Y AND Z, RESPECTIVELY.

    // TO EVAL THE X,Y AND Z COORDS AT ANY [u,v] YOU NEED TO MULTIPLY THREE TERMS:

    //   x = U * Cx * transpose( V )
    //   y = U * Cy * transpose( V )
    //   z = U * Cz * transpose( V )

    // WHERE U = [ u*u*u , u*u , u , 1 ] AND V = [ v*v*v , v*v , v , 1 ]

    // NOW YOU HAVE THE SURFACE POINT p = [x,y,z].

    // TO COMPUTE THE SURFACE NORMAL, YOU NEED TO EVALUATE AT SLIGHTLY
    // DIFFERENT PARAMETRIC LOCATIONS pu = [u+.001,v] AND pv = [u,v+.001].

    // THEN YOU CAN TAKE THE DIFFERENCE VECTORS (pu - p) AND (pv - p).

    // THE CROSS PRODUCT OF THOSE TWO VECTORS IS IN A DIRECTION PERPENDICULAR
    // TO THE SURFACE. YOU CAN NORMALIZE THIS VECTOR TO GET THE SURFACE NORMAL.

    // FINALLY, RETURN [ x, y, z,  nx, ny, nz,  u, v ]
    
    let U = [u*u*u, u*u, u, 1];
    let V = [v*v*v, v*v, v, 1];
    let calculateVMV = (vec1, mat, vec2) => {
        // this function calculates the result of vec1*mat*transpose(vec2),
        // where vec1 and vec2 are a 1x4 vectors and mat is a 4x4 matrix.
        let ret = 0;
        for (let i = 0; i < 4; i++) {
            ret += vec1[i] * (mat[i]*vec2[0] + mat[i+4]*vec2[1] + mat[i+8]*vec2[2] + mat[i+12]*vec2[3]);
        }
        return ret;
    }
    let calculateXYZ = (vec1, vec2) => {
        let XYZ = [];
        for (let i = 0; i < 3; i++) {
            XYZ.push(calculateVMV(vec1, arg[i], vec2));
        }
        return XYZ;
    }
    let xyz = calculateXYZ(U, V);
    
    let dU = [3*u*u, 2*u, 1, 0];
    let dV = [3*v*v, 2*v, 1, 0];
    let Su = calculateXYZ(dU, V);
    let Sv = calculateXYZ(U, dV);
    let norm = cross(Su, Sv);
    return [xyz[0], xyz[1], xyz[2], norm[0], norm[1], norm[2], u, v];
}


////////////////////////////// SCENE SPECIFIC CODE

async function setup(state) {
    hotReloadFile(getPath('week8.js'));

    const images = await imgutil.loadImagesPromise([
       getPath("textures/crate.jpg"),
       getPath("textures/grass-2-large.jpg"),
       getPath("textures/golf_ball.jpg"),
       getPath("textures/blood.jpg"),
       getPath("textures/night_sky.jpg"),
       getPath("textures/teeth.png")
    ]);

    let libSources = await MREditor.loadAndRegisterShaderLibrariesForLiveEditing(gl, "libs", [
        { key : "pnoise"    , path : "shaders/noise.glsl"     , foldDefault : true },
        { key : "sharedlib1", path : "shaders/sharedlib1.glsl", foldDefault : true },      
    ]);
    if (! libSources)
        throw new Error("Could not load shader library");

    // load vertex and fragment shaders from the server, register with the editor
    let shaderSource = await MREditor.loadAndRegisterShaderForLiveEditing(
        gl,
        "mainShader",
        { 
            onNeedsCompilation : (args, libMap, userData) => {
                const stages = [args.vertex, args.fragment];
                const output = [args.vertex, args.fragment];
                const implicitNoiseInclude = true;
                if (implicitNoiseInclude) {
                    let libCode = MREditor.libMap.get('pnoise');
                    for (let i = 0; i < 2; i++) {
                        const stageCode = stages[i];
                        const hdrEndIdx = stageCode.indexOf(';');
                        const hdr = stageCode.substring(0, hdrEndIdx + 1);
                        output[i] = hdr + '\n#line 2 1\n' + 
                                    '#include<pnoise>\n#line ' + (hdr.split('\n').length + 1) + ' 0' + 
                                    stageCode.substring(hdrEndIdx + 1);
                    }
                }
                MREditor.preprocessAndCreateShaderProgramFromStringsAndHandleErrors(
                    output[0],
                    output[1],
                    libMap
                );
            },
            onAfterCompilation : (program) => {
                gl.useProgram(state.program = program);
                state.uColorLoc    = gl.getUniformLocation(program, 'uColor');
                state.uCursorLoc   = gl.getUniformLocation(program, 'uCursor');
                state.uModelLoc    = gl.getUniformLocation(program, 'uModel');
                state.uProjLoc     = gl.getUniformLocation(program, 'uProj');
                state.uTexIndexLoc = gl.getUniformLocation(program, 'uTexIndex');
                state.uTimeLoc     = gl.getUniformLocation(program, 'uTime');
                state.uViewLoc     = gl.getUniformLocation(program, 'uView');
		            state.uTexLoc = [];
            		for (let n = 0 ; n < 8 ; n++) {
            		   state.uTexLoc[n] = gl.getUniformLocation(program, 'uTex' + n);
                               gl.uniform1i(state.uTexLoc[n], n);
            		}
                
                var NL = 2;
                state.uLightsLoc = [];
                for (var i = 0; i < NL; i++) {
                    state.uLightsLoc[i] = {};
                    state.uLightsLoc[i].src = gl.getUniformLocation(program, 'uLights['+i+'].src');
                    state.uLightsLoc[i].col = gl.getUniformLocation(program, 'uLights['+i+'].col');
                }

                var NS = 1;
                state.uMaterialsLoc = [];
                for (var i = 0; i < NS; i++) {
                    state.uMaterialsLoc[i] = {};
                    state.uMaterialsLoc[i].ambient    = gl.getUniformLocation(program, 'uMaterials['+i+'].ambient');
                    state.uMaterialsLoc[i].diffuse    = gl.getUniformLocation(program, 'uMaterials['+i+'].diffuse');
                    state.uMaterialsLoc[i].specular   = gl.getUniformLocation(program, 'uMaterials['+i+'].specular');
                    state.uMaterialsLoc[i].power      = gl.getUniformLocation(program, 'uMaterials['+i+'].power');
                }
            } 
        },
        {
            paths : {
                vertex   : "shaders/vertex.vert.glsl",
                fragment : "shaders/fragment.frag.glsl"
            },
            foldDefault : {
                vertex   : true,
                fragment : false
            }
        }
    );
    if (! shaderSource)
        throw new Error("Could not load shader");

    state.cursor = ScreenCursor.trackCursor(MR.getCanvas());

    state.buffer = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, state.buffer);

    let bpe = Float32Array.BYTES_PER_ELEMENT;

    let aPos = gl.getAttribLocation(state.program, 'aPos');
    gl.enableVertexAttribArray(aPos);
    gl.vertexAttribPointer(aPos, 3, gl.FLOAT, false, bpe * VERTEX_SIZE, bpe * 0);

    let aNor = gl.getAttribLocation(state.program, 'aNor');
    gl.enableVertexAttribArray(aNor);
    gl.vertexAttribPointer(aNor, 3, gl.FLOAT, false, bpe * VERTEX_SIZE, bpe * 3);

    let aUV  = gl.getAttribLocation(state.program, 'aUV');
    gl.enableVertexAttribArray(aUV);
    gl.vertexAttribPointer(aUV , 2, gl.FLOAT, false, bpe * VERTEX_SIZE, bpe * 6);

    for (let i = 0 ; i < images.length ; i++) {
        gl.activeTexture (gl.TEXTURE0 + i);
        gl.bindTexture   (gl.TEXTURE_2D, gl.createTexture());
        gl.texParameteri (gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.REPEAT);
        gl.texParameteri (gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.REPEAT);
        gl.texParameteri (gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR_MIPMAP_NEAREST);
        gl.texParameteri (gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
        gl.texImage2D    (gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, images[i]);
        gl.generateMipmap(gl.TEXTURE_2D);
    }
}

let noise = new ImprovedNoise();
let m = new Matrix();
let turnAngle = 0, cursorPrev = [0,0,0];

function onStartFrame(t, state) {
    if (! state.tStart)
        state.tStart = t;
    state.time = (t - state.tStart) / 1000;

    let cursorValue = () => {
       let p = state.cursor.position(), canvas = MR.getCanvas();
       return [ p[0] / canvas.clientWidth * 2 - 1, 1 - p[1] / canvas.clientHeight * 2, p[2] ];
    }

    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

    let cursorXYZ = cursorValue();
    if (cursorXYZ[2] && cursorPrev[2])
        turnAngle += 2 * (cursorXYZ[0] - cursorPrev[0]);
    cursorPrev = cursorXYZ;

    gl.uniform3fv(state.uCursorLoc     , cursorXYZ);
    gl.uniform1f (state.uTimeLoc       , state.time);

    gl.enable(gl.DEPTH_TEST);
    gl.enable(gl.CULL_FACE);
        
    // camera
    gl.uniform3fv(state.cameraLoc, [0., 0., 5.]);
    
    // Lights
    gl.uniform3fv(state.uLightsLoc[0].src, [-1.,1.,1.]);
    gl.uniform3fv(state.uLightsLoc[0].col, [.9,.9,.9]);

    gl.uniform3fv(state.uLightsLoc[1].src, [1.,1.,-2.]);
    gl.uniform3fv(state.uLightsLoc[1].col, [.4,.4,.4]);
}

function onDraw(t, projMat, viewMat, state, eyeIdx) {
    gl.uniformMatrix4fv(state.uViewLoc, false, new Float32Array(viewMat));
    gl.uniformMatrix4fv(state.uProjLoc, false, new Float32Array(projMat));

    let drawShape = (color, type, vertices, texture) => {
       gl.uniform3fv(state.uColorLoc, color);
       gl.uniformMatrix4fv(state.uModelLoc, false, m.value());
       gl.uniform1i(state.uTexIndexLoc, texture === undefined ? -1 : texture);
       gl.bufferData(gl.ARRAY_BUFFER, new Float32Array( vertices ), gl.STATIC_DRAW);
       gl.drawArrays(type, 0, vertices.length / VERTEX_SIZE);
    }

    m.identity();
    m.rotateY(turnAngle);

    // grass
    m.save();
       m.translate(0,-2,0);
       m.scale(6,.01,6);
       drawShape([1,1,1], gl.TRIANGLES, cubeVertices, 1);
    m.restore();

    /// crates and lasers
    m.save();
      m.translate(0,noise.noise(6, 7, 100 * -1 + state.time / 3) / 7 - 0.1, 0);
      // crates
      for (let z = -4 ; z <= 4 ; z += 2)
      for (let x = -4 ; x <= 4 ; x += 2) {
         m.save();
            let y = -1.5;
            m.translate(x, y, z);
            m.scale(.3,.3,.3);
            drawShape([1,1,1], gl.TRIANGLES, cubeVertices, 0);
         m.restore();
      }
      
      // laser grid
      for (let x = -4; x <= 4; x += 2) {
        m.save();
            let y = -1.5;
            let z = 0;
            let c = 0.5*Math.abs(Math.cos(state.time));
            m.translate(x, y, z);
            m.rotateY(Math.PI/2);
            m.scale(4, .01, .01);
            drawShape([1,c,c], gl.TRIANGLES, cubeVertices);
         m.restore();
      }
      for (let z = -4; z <= 4; z += 2) {
        m.save();
            let y = -1.5;
            let x = 0;
            let c = 0.5*Math.abs(Math.cos(state.time));
            m.translate(x, y, z);
            m.scale(4, .01, .01);
            drawShape([1,c,c], gl.TRIANGLES, cubeVertices);
         m.restore();
      }
    m.restore();
    
    /////// background ///////
    let z = -.3;
    let backgroundP = [
      [
        -1,-1/3, 1/3, 1,
        -1,-1/3, 1/3, 1,
        -1,-1/3, 1/3, 1,
        -1,-1/3, 1/3, 1
      ],
      [
        -1  ,-1.  ,-1  ,-1,
        -1/3,-1/3,-1/3,-1/3,
         1/3, 1/3, 1/3, 1/3,
         1  , 1  , 1  , 1
      ],
      [
        0,  z,  z,  0,
        0,  z,  z,  0,
        0,  z,  z,  0,
        0,  z,  z,  0
      ]
    ];
    let background = createMeshVertices(16, 16, uvToCubicPatch,
       toCubicPatchCoefficients(BezierBasisMatrix, backgroundP)
    );
    m.save();
    m.translate(0,0,-4);
      for (let i = 0; i < 2*Math.PI; i += Math.PI/3){
        m.save();
        m.rotateY(i);
        m.translate(0,0,-10.35);
        m.scale(6,6,6);
        gl.uniform3fv(state.uMaterialsLoc[0].ambient , [1,1,1]);
        gl.uniform3fv(state.uMaterialsLoc[0].diffuse , [.1,.1,.1]);
        gl.uniform3fv(state.uMaterialsLoc[0].specular, [.1,.1,.1]);
        gl.uniform1f (state.uMaterialsLoc[0].power   , 15.);
        drawShape([1,1,1], gl.TRIANGLE_STRIP, background, 4);
        m.restore();
      }
    m.restore();
    
    /////// Eye of Cthulhu ///////
    m.save();
    m.rotateY(-0.2 * state.time);
    m.translate(0,noise.noise(6, 7, 100 * -1 + state.time / 3) / 7 - 0.1, -1.2);
    m.rotateY(-0.1+noise.noise(6, 7, 100 * -2 + state.time / 3) / 7);
    m.rotateZ(noise.noise(3, 2, 100 * -1 + state.time / 3) / 7);
    m.scale(.6,.6,.6);

    // helper function
    let frontToBack = fP => {
      let bP = [];
      let flip = M => {
        let ret = [];
        for (let i = 0; i < 16; i+=4) {
          ret[i] = M[12 - i];
          ret[i+1] = M[13 - i];
          ret[i+2] = M[14 - i];
          ret[i+3] = M[15 - i];
        }
        return ret;
      }
      let reverse = M => {
        let ret = [];
        for (let i = 0 ; i < M.length ; i++)
          ret.push(-M[i]);
        return ret;
      }
      bP = [flip(fP[0]), flip(fP[1]), reverse(flip(fP[2]))];
      return bP;
    }
      // eye ball //
      m.save();
      let eyeFrontP = [
         [
           -1, -1, .2, .5,
           -1,-.9,-.3,-.1,
           -1,-.9,-.3,-.1,
           -1, -1, .2, .5
         ],
         [
            0, -1, -1,-.5,
            0,-.6,-.5,-.5,
            0, .6, .5, .5,
            0,  1,  1, .5
         ],
         [
            0,  0,   0,  0,
            0,  1, 1.3, .7,
            0,  1, 1.3, .7,
            0,  0,   0,  0
         ]
      ];
      let eyeBackP = frontToBack(eyeFrontP);
      let eyeFront = createMeshVertices(16, 16, uvToCubicPatch,
         toCubicPatchCoefficients(BezierBasisMatrix, eyeFrontP)
      );
      let eyeBack = createMeshVertices(16, 16, uvToCubicPatch,
         toCubicPatchCoefficients(BezierBasisMatrix, eyeBackP)
      );

      m.scale(.6,.6,.6);
      gl.uniform3fv(state.uMaterialsLoc[0].ambient , [.6,.6,.6]);
      gl.uniform3fv(state.uMaterialsLoc[0].diffuse , [.9,.9,.9]);
      gl.uniform3fv(state.uMaterialsLoc[0].specular, [.8,.8,.8]);
      gl.uniform1f (state.uMaterialsLoc[0].power   , 2.);
      drawShape([1,1,1], gl.TRIANGLE_STRIP, eyeFront, 2);
      drawShape([1,1,1], gl.TRIANGLE_STRIP, eyeBack, 2);
      
      // mouth //
      let mouthFrontP = [
         [
          .07, .1, .3, .5,
          .07,  0,-.1,-.4,
          .07,  0,-.1,-.4,
          .07, .1, .3, .5
         ],
         [
           0,-.3,-.5,-.5,
           0,-.3,-.5,-.5,
           0, .3, .5, .5,
           0, .3, .5, .5
         ],
         [
          .54, .54, .3,  0,
          .54, .54, .3,  0,
          .54, .54, .3,  0,
          .54, .54, .3,  0
         ]
      ];
      let mouthBackP = frontToBack(mouthFrontP);
      let mouthFront = createMeshVertices(16, 16, uvToCubicPatch,
         toCubicPatchCoefficients(BezierBasisMatrix, mouthFrontP)
      );
      let mouthBack = createMeshVertices(16, 16, uvToCubicPatch,
         toCubicPatchCoefficients(BezierBasisMatrix, mouthBackP)
      );

      let mouthColor = [.8,.07,.08]
      gl.uniform3fv(state.uMaterialsLoc[0].ambient , [.16,.007,.008]);
      gl.uniform3fv(state.uMaterialsLoc[0].diffuse , [.4,.035,.04]);
      gl.uniform3fv(state.uMaterialsLoc[0].specular, [.6,.6,.6]);
      gl.uniform1f (state.uMaterialsLoc[0].power   , 30.);
      drawShape(mouthColor, gl.TRIANGLE_STRIP, mouthFront);
      drawShape(mouthColor, gl.TRIANGLE_STRIP, mouthBack);
      m.restore();
    
      // tooth //
      let toothFrontP = [
         [
           0, 0, 0, 0,
           -.05,-.05/3,.05/3, .05,
           -.1,-.1/3,.1/3, .1,
           -.15, -.05, .05, .15
         ],
         [
           -1, -1, -1, -1,
          -1/3,-1/3,-1/3,-1/3,
           1/3, 1/3, 1/3, 1/3,
            1,  1,  1,  1
         ],
         [
            0,  0,   0,  0,
            0, .2/3,.2/3,  0,
            0, .4/3,.4/3,  0,
            0, .2,  .2,  0
         ]
      ];
      let toothBackP = frontToBack(toothFrontP);
      let toothFront = createMeshVertices(16, 16, uvToCubicPatch,
         toCubicPatchCoefficients(BezierBasisMatrix, toothFrontP)
      );
      let toothBack = createMeshVertices(16, 16, uvToCubicPatch,
         toCubicPatchCoefficients(BezierBasisMatrix, toothBackP)
      );
      
      let toothColor = [.48, .43, .35];
      gl.uniform3fv(state.uMaterialsLoc[0].ambient , [.09,.065,.025]);
      gl.uniform3fv(state.uMaterialsLoc[0].diffuse , [.18, .13, .05]);
      gl.uniform3fv(state.uMaterialsLoc[0].specular, [.6,.6,.6]);
      gl.uniform1f (state.uMaterialsLoc[0].power   , 3.);
      for (let updown = -1; updown <= 1; updown += 2) {
        m.save();
        m.translate(.27,updown*.21,0);
        m.rotateZ(updown == -1 ? (Math.PI - 0.43) : 0.43);
        m.scale(.2,.09,.2);
        drawShape(toothColor, gl.TRIANGLE_STRIP, toothFront, 5);
        drawShape(toothColor, gl.TRIANGLE_STRIP, toothBack, 5);
        m.restore();
        for (let side = -1; side <= 1; side += 2) {
          m.save();
          m.translate(.18,updown*.2,side * 0.1);
          m.rotateZ(updown == -1 ? (Math.PI - 0.45) : 0.45);
          m.scale(.2,.1,.2);
          drawShape(toothColor, gl.TRIANGLE_STRIP, toothFront, 5);
          drawShape(toothColor, gl.TRIANGLE_STRIP, toothBack, 5);
          m.restore();
          m.save();
          m.translate(.11,updown*.14,side * 0.18);
          m.rotateZ(updown == -1 ? (Math.PI - 0.5) : 0.5);
          m.scale(.2,.1,.2);
          drawShape(toothColor, gl.TRIANGLE_STRIP, toothFront, 5);
          drawShape(toothColor, gl.TRIANGLE_STRIP, toothBack, 5);
          m.restore();
          m.save();
          m.translate(.06,updown*.06,side * 0.26);
          m.rotateZ(updown == -1 ? (Math.PI - 0.6) : 0.6);
          m.scale(.2,.05,.2);
          drawShape(toothColor, gl.TRIANGLE_STRIP, toothFront, 5);
          drawShape(toothColor, gl.TRIANGLE_STRIP, toothBack, 5);
          m.restore();
        }
      }
      
      
      /// blood vines and tails
      let bloodVineColor = [.5, .02, .03];
      let bloodVinesAndTailsPFrontToBack = fP => {
        let bP = [];
        let reverseX = X => {
          let ret = [];
          for (let i = 0; i < X.length; i++) {
            ret.push(-X[i]);
          }
          return ret;
        }
        bP.push(reverseX(fP[0]));
        bP.push(fP[1]);
        bP.push(fP[2]);
        return bP;
      }
      
      // blood vines //
      let bloodVinesFrontP = [];
      bloodVinesFrontP[0] = [
        [ -.6, -.5, 0,  0.09], // A.x B.x C.x D.x
        [  0,  .2,  0,  0.2], // A.y B.y C.y D.y
        [  0,  .6, .5, .25]  // A.z B.z C.z D.z
      ];
      bloodVinesFrontP[1] = [
        [ -.58, -.5, 0,  0.03], // A.x B.x C.x D.x
        [ -.2,  -.2,  0, -0.05], // A.y B.y C.y D.y
        [  0,  .55, .55, .31]  // A.z B.z C.z D.z
      ];
      bloodVinesFrontP[2] = [
                  [ -.45,-.4, -.2,  0], // A.x B.x C.x D.x
                  [ -.4, -.4, -.3, -0.3], // A.y B.y C.y D.y
                  [   0,  .25, .45, .28]  // A.z B.z C.z D.z
      ];
      bloodVinesFrontP[3] = [
                  [ -.23,-.1, -.15,  0.2], // A.x B.x C.x D.x
                  [ -.43,-.4, -.45, -0.27], // A.y B.y C.y D.y
                  [  .15, .2, .3, .18]  // A.z B.z C.z D.z
      ];
      bloodVinesFrontP[4] = [
                  [ -.41,-.26, -.15,0.2], // A.x B.x C.x D.x
                  [ .35,.4, .45, 0.27], // A.y B.y C.y D.y
                  [  .15, .2, .3, .18]  // A.z B.z C.z D.z
      ];
      bloodVinesFrontP[5] = [
                  [-.1,-.08, -.06, 0], // A.x B.x C.x D.x
                  [ .1, .15,.1, .2], // A.y B.y C.y D.y
                  [ .43, .43, .4, .35]  // A.z B.z C.z D.z
      ];
      bloodVinesFrontP[6] = [
                  [-.3,-.25, -.2, -.15], // A.x B.x C.x D.x
                  [ .1, .17,.2, .2], // A.y B.y C.y D.y
                  [ .43, .41, .42, .4]  // A.z B.z C.z D.z
      ];
      bloodVinesFrontP[7] = [
                  [-.2,-.17, -.15, -.05], // A.x B.x C.x D.x
                  [ .1, .05,-.05, -.0], // A.y B.y C.y D.y
                  [ .45, .47, .46, .42]  // A.z B.z C.z D.z
      ];
      bloodVinesFrontP[8] = [
                  [-.43,-.25, -.2, -.13], // A.x B.x C.x D.x
                  [-.16, 0,.05, -.005], // A.y B.y C.y D.y
                  [ .35, .5, .46, .45]  // A.z B.z C.z D.z
      ];
      bloodVinesFrontP[9] = [
                  [-.3, -.25, -.2, -.13], // A.x B.x C.x D.x
                  [-.12,-.18,-.24, -.24], // A.y B.y C.y D.y
                  [ .43, .41, .4, .39]  // A.z B.z C.z D.z
      ];
      bloodVinesFrontP[10] = [
                  [-.13, -.13, -.1, -.03], // A.x B.x C.x D.x
                  [-.06,-.06,-.13, -.14], // A.y B.y C.y D.y
                  [ .45, .45, .45, .39]  // A.z B.z C.z D.z
      ];
      bloodVinesFrontP[11] = [
                  [-.45, -.23, -.1, .08], // A.x B.x C.x D.x
                  [-.4,-.45,-.46, -.18], // A.y B.y C.y D.y
                  [ 0, .22, .3, .3]  // A.z B.z C.z D.z
      ];
      bloodVinesFrontP[12] = [
                  [-.35, -.3, -.25, -.25], // A.x B.x C.x D.x
                  [-.37,-.25,-.24, -.24], // A.y B.y C.y D.y
                  [ .22, .35, .4, .35]  // A.z B.z C.z D.z
      ];
      bloodVinesFrontP[13] = [
                  [.03,.03, .02, .02], // A.x B.x C.x D.x
                  [ -.27, -.25,-.15, -.15], // A.y B.y C.y D.y
                  [ .25, .34, .33, .35]  // A.z B.z C.z D.z
      ];
      bloodVinesFrontP[14] = [
                  [-.5,-.4, -.2, -.15], // A.x B.x C.x D.x
                  [ .2, .2,.37, .3], // A.y B.y C.y D.y
                  [ .23, .35, .33, .34]  // A.z B.z C.z D.z
      ];
      bloodVinesFrontP[15] = [
                  [-.15,-.05, .05, .1], // A.x B.x C.x D.x
                  [ .3, .2,.2, .31], // A.y B.y C.y D.y
                  [ .34, .4, .4,.21]  // A.z B.z C.z D.z
      ];
      bloodVinesFrontP[16] = [
                  [-.4,-.36, -.3, -.3], // A.x B.x C.x D.x
                  [ .23, .3,.321, .35], // A.y B.y C.y D.y
                  [ .3, .29, .28, .25]  // A.z B.z C.z D.z
      ];
      bloodVinesFrontP[17] = [
                  [-.43,-.36, -.3, -.3], // A.x B.x C.x D.x
                  [ .22, .18,.18, .14], // A.y B.y C.y D.y
                  [ .28, .37, .41, .41]  // A.z B.z C.z D.z
      ];
      bloodVinesFrontP[18] = [
                  [-.25,-.2,-.1,-.1], // A.x B.x C.x D.x
                  [ .4, .44, .5, .5], // A.y B.y C.y D.y
                  [ .2, .15,.04,  0]  // A.z B.z C.z D.z
      ];
      bloodVinesFrontP[19] = [
                  [   0, .1, .2,  .25], // A.x B.x C.x D.x
                  [ .36, .4, .4,  .35], // A.y B.y C.y D.y
                  [ .23, .13,.04, .02]  // A.z B.z C.z D.z
      ];
      bloodVinesFrontP[20] = [
                  [   0, .1, .2,  .25], // A.x B.x C.x D.x
                  [-.36,-.4,-.4, -.35], // A.y B.y C.y D.y
                  [ .23, .13,.04, .02]  // A.z B.z C.z D.z
      ];
      bloodVinesFrontP[21] = [
                  [ -.2, -.1, 0,  .05], // A.x B.x C.x D.x
                  [-.42,-.5,-.5, -.47], // A.y B.y C.y D.y
                  [ .18, .05,.04, .02]  // A.z B.z C.z D.z
      ];
      bloodVinesFrontP[22] = [
                  [ -.5, -.4,-.34, -.34], // A.x B.x C.x D.x
                  [ .08,   0,  0, -.05], // A.y B.y C.y D.y
                  [ .3, .38, .45,  .4]  // A.z B.z C.z D.z
      ];
      bloodVinesFrontP[23] = [
                  [-.15, -.1,-.05,  .1], // A.x B.x C.x D.x
                  [ .44, .42, .42, .46], // A.y B.y C.y D.y
                  [ .14, .18,  .2, .01]  // A.z B.z C.z D.z
      ];
      
      gl.uniform3fv(state.uMaterialsLoc[0].ambient , [.25, .01, .015]);
      gl.uniform3fv(state.uMaterialsLoc[0].diffuse , [.5, .02, .03]);
      gl.uniform3fv(state.uMaterialsLoc[0].specular, [.6,.6,.6]);
      gl.uniform1f (state.uMaterialsLoc[0].power   , 10.);
      m.save();
      let bloodVinesFront = [];
      let bloodVinesBack = [];
      for (let i = 0; i < bloodVinesFrontP.length; i++) {
        m.save();
        bloodVinesFront[i] = createMeshVertices(16, 2, uvToCubicCurvesRibbon,
          {
            width: i <= 4 ? 0.02 : 0.01,
            data: [
               toCubicCurveCoefficients(BezierBasisMatrix, bloodVinesFrontP[i])
            ]
          }
        );
        drawShape(bloodVineColor, gl.TRIANGLE_STRIP, bloodVinesFront[i], 3);
        m.restore();
        m.save();
        m.rotateY(Math.PI);
        bloodVinesBack[i] = createMeshVertices(16, 2, uvToCubicCurvesRibbon,
          {
            width: i <= 4 ? 0.02 : 0.01,
            data: [
               toCubicCurveCoefficients(BezierBasisMatrix, bloodVinesAndTailsPFrontToBack(bloodVinesFrontP[i]))
            ]
          }
        );
        drawShape(bloodVineColor, gl.TRIANGLE_STRIP, bloodVinesBack[i], 3);
        m.restore();
      }
      m.restore();
      
      // tails //
      let tailsFrontP = [];
      let c = Math.cos(7 * state.time);
      let s = c > .5 ? .8 : (c < -.5 ? -.4 : .4);
      tailsFrontP[0] = [
        [-.53,  -.62,  -.65,  -.7], // A.x B.x C.x D.x
        [   0, .05*s, .05*s,    0], // A.y B.y C.y D.y
        [   0,     0,     0,    0]  // A.z B.z C.z D.z
      ];
      tailsFrontP[1] = [
        [ -.7,  -.75,  -.75, -.85], // A.x B.x C.x D.x
        [   0,-.05*s,  .1*s,    0], // A.y B.y C.y D.y
        [   0,     0,     0,    0]  // A.z B.z C.z D.z
      ];
      tailsFrontP[2] = [
        [-.85,  -.95,  -.95,   -1], // A.x B.x C.x D.x
        [   0, -.1*s,  .1*s,.05*s], // A.y B.y C.y D.y
        [   0,     0,     0,    0]  // A.z B.z C.z D.z
      ];
      tailsFrontP[3] = [
        [-.53,  -.62,  -.65,  -.7], // A.x B.x C.x D.x
        [ .18, .18-.05*s,  .18-.05*s,  .18+.05*s], // A.y B.y C.y D.y
        [  0,      0,     0,    0]  // A.z B.z C.z D.z
      ];
      
      let tailColor = bloodVineColor;
      for (let theta = - Math.PI/3; theta < Math.PI/3; theta += Math.PI/7) {
        m.save();
        m.rotateZ(theta);
        m.translate(0, -0.14, 0);
        let tailsFront = [];
        let tailsBack = [];
        for (let i = 0; i < tailsFrontP.length; i++) {
          m.save();
          tailsFront[i] = createMeshVertices(16, 2, uvToCubicCurvesRibbon,
            {
              width: i <= 2 ? 0.03 : 0.025,
              data: [
                 toCubicCurveCoefficients(BezierBasisMatrix, tailsFrontP[i])
              ]
            }
          );
          drawShape(tailColor, gl.TRIANGLE_STRIP, tailsFront[i], 3);
          m.restore();
          m.save();
          m.rotateY(Math.PI);
          tailsBack[i] = createMeshVertices(16, 2, uvToCubicCurvesRibbon,
            {
              width: i <= 2 ? 0.03 : 0.025,
              data: [
                 toCubicCurveCoefficients(BezierBasisMatrix, bloodVinesAndTailsPFrontToBack(tailsFrontP[i]))
              ]
            }
          );
          drawShape(tailColor, gl.TRIANGLE_STRIP, tailsBack[i], 3);
          m.restore();
        }
        m.restore();
      }
      
    m.restore();
}

function onEndFrame(t, state) {}

export default function main() {
    const def = {
        name         : 'week8',
        setup        : setup,
        onStartFrame : onStartFrame,
        onEndFrame   : onEndFrame,
        onDraw       : onDraw,
    };

    return def;
}

