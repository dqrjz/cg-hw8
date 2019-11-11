"use strict"

////////////////////////////// MATRIX SUPPORT

let cos = t => Math.cos(t);
let sin = t => Math.sin(t);
let identity = ()       => [1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1];
let rotateX = t         => [1,0,0,0, 0,cos(t),sin(t),0, 0,-sin(t),cos(t),0, 0,0,0,1];
let rotateY = t         => [cos(t),0,-sin(t),0, 0,1,0,0, sin(t),0,cos(t),0, 0,0,0,1];
let rotateZ = t         => [cos(t),sin(t),0,0, -sin(t),cos(t),0,0, 0,0,1,0, 0,0,0,1];
let scale = (x,y,z)     => [x,0,0,0, 0,y,0,0, 0,0,z,0, 0,0,0,1];
let translate = (x,y,z) => [1,0,0,0, 0,1,0,0, 0,0,1,0, x,y,z,1];
let multiply = (a, b)   => {
   let c = [];
   for (let n = 0 ; n < 16 ; n++)
      c.push( a[n&3     ] * b[    n&12] +
              a[n&3 |  4] * b[1 | n&12] +
              a[n&3 |  8] * b[2 | n&12] +
              a[n&3 | 12] * b[3 | n&12] );
   return c;
}

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
       getPath("textures/brick.png"),
       getPath("textures/tiles.jpg"),
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

    m.save();
       m.translate(0,-2,0);
       m.scale(6,.01,6);
       drawShape([1,1,1], gl.TRIANGLES, cubeVertices, 1);
    m.restore();

    for (let z = -3 ; z <= 3 ; z += 2)
    for (let x = -3 ; x <= 3 ; x += 2) {
       m.save();
          let y = Math.max(Math.abs(x),Math.abs(z)) / 3 - 1 +
	          noise.noise(x, 0, 100 * z + state.time / 2) / 5;
          m.translate(x, y, z);
          m.scale(.3,.3,.3);
          drawShape([1,1,1], gl.TRIANGLES, cubeVertices, 0);
       m.restore();
    }
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

