
rooms.raytrace5 = function () {

    lib3D();
 
    description = `
    <style>
    input{
       font-family: 'Roboto', sans-serif;
    }
    </style>
 <small>
 
 
     <br> <input type=range id=refract value=10> Refraction
     <br>
     <br>
     
      <input type=range id=transZ  value=0> 
      Try Different View :D
      <div id=iorInfo></div>
      <br>
      <div>
       <button id=nav2 onclick="handleZPlus()" onkeydown="handleZPlus()">Z +</button>
       <br>
       <button id=nav1 onclick="handleXPlus()" onkeydown="handleXPlus()"
       > X +</button>
       <button id=nav4 onclick="handleXMinus()" onkeydown="handleXMinus()"
       >X -</button>
       <br>
       <button id=nav3 onclick="handleZMinus()" onkeydown="handleZMinus()">Z -</button>
      </div>
      <br>
       <button id=light onclick="handleLight()">BG</button>
     
 </small>
 `;
 
    code = {
       'init': `
 
    // DEFINE MATERIALS TO BE RENDERED VIA PHONG REFLECTANCE MODEL
 
    S.redPlastic    = [.2,.1,.1,0,  .5,.2,.2,0,  2,2,2,20,  0,0,0,0];
    S.greenPlastic  = [.1,.2,.1,0,  .2,.5,.2,0,  2,2,2,20,  0,0,0,0];
    S.bluePlastic   = [.1,.1,.2,0,  .2,.2,.5,0,  2,2,2,20,  0,0,0,0];
    S.whitePlastic  = [.2,.2,.2,0,  .5,.5,.5,0,  2,2,2,20,  0,0,0,0];
    S.blackGlass    = [.1,.1,.1,0,  .1,.1,.1,0,  2,2,2,200,  0,0,0,0];
    S.matteGray     = [.0,.0,.0,0,  .5,.5,.5,0,  .2,.2,.2,20,  0,0,0,0];
    S.translucentWhite   = [.0,.0,.0,0,  .5,.5,.5,0,  .2,.2,.2,20,  0,0,0,.5];
    S.gold = [.25,.15,.025,0, .5,.3,.05,0, 1,.6,.1,600,  0,0,0,0]; // GOLD
    S.whitebright = [.0,.0,.0,0,  .5,.5,.5,0,  2,2,2,20,  0,0,0,0];
 
    // Dark
 
       S.bg= 10;
       S.lightOn = false;
       S.nav_x = 0;
       S.nav_z = 0;
 
 
       // Light
 `,
 
       fragment: `
 S.setFragmentShader(\`
 
    // DECLARE CONSTANTS, UNIFORMS, VARYING VARIABLES
 
    const int nQ = \` + S.nQ + \`;
    const int nL = \` + S.nL + \`;
    uniform vec3 uBgColor;
    uniform vec3 uLd[nL];
    uniform vec3 uLpos[nL];
    uniform vec3 uLc[nL];
    uniform mat4 uQ[nQ];
    uniform mat4 uPhong[nQ];
    uniform int  uShape[nQ];
    uniform float uIor;
    varying vec3 vPos;
    uniform float uFl;
    // float fl = 3.0;
 
    // <p>  <input type=range id=red   > red
    // <br> <input type=range id=green > green
    // <br> <input type=range id=blue  > blue
 
    
 
    vec3 normalQ(vec3 P, mat4 Q){
       // calculate normal to quadric Q at point P
       float fx = 2.*Q[0][0]*P.x + (Q[0][1]+Q[1][0])*P.y + (Q[0][2]+Q[2][0])*P.z + (Q[0][3]+Q[3][0]);
       float fy = (Q[0][1]+Q[1][0])*P.x + 2.*Q[1][1]*P.y + (Q[1][2]+Q[2][1])*P.z + (Q[1][3]+Q[3][1]);
       float fz = (Q[0][2]+Q[2][0])*P.x + (Q[1][2]+Q[2][1])*P.y + 2.*Q[2][2]*P.z + (Q[2][3]+Q[3][2]);
       
       return normalize(vec3(fx,fy,fz));
    }
 
 
    vec2 rayQ(vec3 V, vec3 W, mat4 Q){
       vec4 W0 = vec4(W,0.);
       vec4 V1 = vec4(V,1.);
       // calculate intersection of ray VW with quadric 
 
          float a = dot(W0, Q*W0);
          float b = dot(V1, Q*W0) + dot(W0, Q*V1);
          float c = dot(V1, Q*V1);
          float d = b*b - 4.*a*c;
          if (d < 0.) return vec2(-1.,-1.);
          float t1 = (-b - sqrt(d)) / (2.*a);
          float t2 = (-b + sqrt(d)) / (2.*a);
          // if (t1 < 0.) 
          //    return vec2(-1.,-1.);
          return vec2(t1,t2);
          // Return both roots as a vec2
       }
    
    vec3 Q1(vec3 T, int n, vec2 t){
          float tIn  = t.x;
          float tOut = t.y;
          if (tIn > 0. && tIn < tOut && tIn < T.y)
             T = vec3(n, tIn, tOut);
          return T;
       }
 
    vec3 Q2(vec3 T, int n, vec2 t0, vec2 t1){
          float tIn  = max(t0.x, t1.x);
          float tOut = min(t0.y, t1.y);
          if (tIn > 0. && tIn < tOut && tIn < T.y){
             int i = t0.x==tIn ? 0 : 1;
             T = vec3(n+i,tIn,tOut);
          }
          return T;
       }
    
    vec3 Q3(vec3 T, int n, vec2 t0, vec2 t1, vec2 t2){
          float tIn  = max(t0.x, max(t1.x, t2.x));
          float tOut = min(t0.y, min(t1.y, t2.y));
          if (tIn > 0. && tIn < tOut && tIn < T.y){
             int i = t0.x==tIn ? 0 : t1.x==tIn ? 1 : 2;
             T = vec3(n+i,tIn,tOut);
          }
          return T;
       }
 
 
    vec3 Q4(vec3 T, int n, vec2 t0, vec2 t1, vec2 t2, vec2 t3){
          float tIn  = max(t0.x,max(max(t1.x,t2.x),t3.x));
          float tOut = min(t0.y,min(min(t1.y,t2.y),t3.y));
          if (tIn > 0. && tIn < tOut && tIn < T.y){
             int i = t0.x==tIn ? 0 : t1.x==tIn ? 1 : t2.x==tIn ? 2 : 3;
             T = vec3(n+i,tIn,tOut);
          }
          return T;
       }
 
    vec3 rayScene(vec3 V, vec3 W){
       vec3 T = vec3(-1.,1000., 0.);
       for (int i = 0; i < nQ; i++){
          if (uShape[i] == 1) T = Q1(T, i, rayQ(V,W,uQ[i]));
          if (uShape[i] == 2) T = Q2(T, i, rayQ(V,W,uQ[i]),rayQ(V,W,uQ[i+1]));
          if (uShape[i] == 3) T = Q3(T, i, rayQ(V,W,uQ[i]),rayQ(V,W,uQ[i+1]),rayQ(V,W,uQ[i+2]));
          if (uShape[i] == 4) T = Q4(T, i, rayQ(V,W,uQ[i]),rayQ(V,W,uQ[i+1]),rayQ(V,W,uQ[i+2]),rayQ(V,W,uQ[i+3]));
       }
       return T;
    }
 
    vec3 shadeSurface(vec3 P, vec3 N, mat4 phong){
       vec3  ambient  = phong[0].rgb;
       vec3  diffuse  = phong[1].rgb;
       vec3  specular = phong[2].rgb;
       float p        = phong[2].a;
 
       // vec3 lightPos = vec3(0.,0.,0.);
       // vec3 lihgtN = lightPos - P;
       // vec3 lightDir = normalize(lightPos - P);
 
       // vec3 N = N;
       vec3 c = mix(ambient, uBgColor, .3);
       vec3 E = vec3(0.,0.,1.);
 
       for (int i = 0; i < nL; i++){
          vec3 lightPos = uLpos[i];
          // vec3 lightPos = vec3(uLpos[i][0], uLpos[i][1], uLpos[i][2]);
          // vec3 lightPos = vec3(-.5,.5,.5);
          vec3 lightN = lightPos - P;
          vec3 lightDir = normalize(lightPos - P);
          float dist = length(lightN);
 
          // float t = -1.;
          // float tIn = rayScene(P, uLd[i]).y;
          // if (tIn != 1000.){
          //    t = max(t, tIn);
          // }
 
          float t = -1.;
          float tIn = rayScene(P, lightPos).y;
          if (tIn != 1000.){
             t = max(t, tIn);
          }
 
          if (dist < 0.1) dist = 0.01;
          float attenuation = 1. / (1. + .2 * dist + .3 * dist * dist);
 
           // COMPUTE DIFFUSE AND SPECULAR FOR THIS LIGHT SOURCE
 
          // if (t < 0.) {
          //     vec3 R = 1. * dot(N, uLd[i]) * N - uLd[i];
          //     c += uLc[i] * (diffuse * max(0.,dot(N, uLd[i]))
          //                 + specular * pow(max(0., dot(R, E)), p));
          // }
 
          if (t < 0.) {
             vec3 R = 1. * dot(N, lightDir) * N - lightDir;
             vec3 lightColor = vec3(.4,.15,.1);
             c += lightColor * (attenuation * diffuse * max(0.,dot(N, lightDir))
                         + attenuation * specular * pow(max(0., dot(R, E)), p));
          }
 
 // Spot Light
          lightPos = vec3(-.5,.6,.1);
          vec3 lightTarget = vec3(-.1,-.01,0.);
          float intensity = .2;
          vec3 lightColor = vec3(.2,.3,.3) * intensity;
          lightN = lightPos - P;
          vec3 lightStraight = normalize(lightPos - lightTarget);
          lightDir = normalize(lightPos - P);
          dist = length(lightN);
          float theta = dot(lightDir, lightStraight);
          if(theta > 0.76) {
          t = -1.;
          tIn = rayScene(P, lightPos).y;
          if (tIn != 1000.){
             t = max(t, tIn);
          }
          if (dist < 0.1) dist = 0.01;
          attenuation = 1. / (1. + .1 * dist + .3 * dist * dist);
          if (t < 0.) {
             vec3 R = 1. * dot(N, lightDir) * N - lightDir;
             
             c += lightColor * (attenuation * diffuse * max(0.,dot(N, lightDir))
                          + attenuation * .01 * specular * pow(max(0., dot(R, E)), p));
          }
          }else{
          }
       }
       return c;
    }
 
    vec3 refractRay(vec3 W, vec3 N, float n){
       vec3 C = dot(W,N) * N;
       vec3 S = W - C;   
       float sin1 = length(S);
       float sin2 = sin1 * n;
       float cos1 = sqrt(1. - sin1*sin1);
       float cos2 = sqrt(1. - sin2*sin2);
       vec3 C2 = C * cos2 / cos1;
       vec3 S2 = S * sin2 / sin1;
       vec3 W2 = C2 + S2;
       return W2;
    }
 
 
 // /********* CODE IMPLEMENTATION OF FRAGMENT SHADER ********** 
 
 
    void main() {
       
       vec3 color = uBgColor;
 
       vec3 V = vec3(0.,0.,uFl);
       vec3 W = normalize(vec3(vPos.xy, -uFl));
 
       vec3 T = rayScene(V, W);
 
       for (int n = 0 ; n < nQ ; n++){
          if (n == int(T.x)) {
             mat4 Q = uQ[n]; 
             // if(n >=0) {
                vec3 P = V + T.y * W;
                vec3 N = normalQ(P, Q);
                vec3 R = 2. * dot(N, -W) * N + W;
 
                color += shadeSurface(P, N, uPhong[n]);
 
                // Do reflection
                T = rayScene(P, R);
 
                for (int i = 0; i < nQ; i++){
                   if (i == int(T.x)) {
                      // if (i >= 0) {
                         vec3 M = P + T.y * R;
                         // REFLECTION INTENSITY
                         color += shadeSurface(M, normalQ(M,uQ[i]), uPhong[i])/2.;
                      // }
                   }
                }
                
                // DO refraction
                // (1) SHOOT RAY TO FIND REAR OF THIS OBJECT (USE 2nd ROOT):
 
                W = refractRay(W, N, uIor);
                T = rayScene( P-.01*W , W);
                for (int j = 0; j < nQ; j++){
                   if (j == int(T.x)) {
                      P = P + T.z * W;
                      N = normalQ(P, uQ[j]);
                   }
                }
                
                // (2) SHOOT RAY FROM REAR OF THIS OBJECT TO FIND ANOTHER OBJECT:
 
                W = refractRay(W, N, 1./uIor);
                T = rayScene( P , W);
 
                for (int k = 0; k < nQ; k++){
                   if (k == int(T.x)) {
                      if (k >= 0) {
                         vec3 M = P + T.y * W;
                         color += uPhong[k][1].xyz * shadeSurface(M, normalQ(M,uQ[k]), uPhong[k]);
                      }
                   }
                }
             // }
          }
       }
       
 
       gl_FragColor = vec4(color, 1);
       vec4 BrightColor;
       float brightness = dot(color, vec3(0.2126, 0.7152, 0.0722));
       if (brightness > 1.) {
          BrightColor = gl_FragColor;
       } else{
          BrightColor = vec4(0.0, 0.0, 0.0, 1.0);
       }
    }
 \`);
 `,
       vertex: `
 S.setVertexShader(\`
 
    attribute vec3 aPos;
    varying   vec3 vPos;
 
    void main() {
       vPos = aPos;
       gl_Position = vec4(aPos, 1.);
    }
 
 \`)
 
 `,
       render: `
 
    // USEFUL VECTOR FUNCTIONS
 
    let add = (a,b) => [ a[0]+b[0], a[1]+b[1], a[2]+b[2] ];
    let dot = (a,b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let norm = v => Math.sqrt(dot(v,v));
    let normalize = v => { let s = norm(v); return [ v[0]/s, v[1]/s, v[2]/s ]; }
    let scale = (v,s) => [ s * v[0], s * v[1], s * v[2] ];
    let subtract = (a,b) => [ a[0]-b[0], a[1]-b[1], a[2]-b[2] ];
    let z_value = -transZ.value / 50 * 3 + 3;
 
    // SEND LIGHT SOURCE DATA TO GPU
 
    let ldData = [ normalize([.1,.1,.1]),
                   normalize([-1,-1,-1]) ];
    // let lposData = [ [ -.5,.5,.5 ], [ 0+ S.nav_x, 0, 0 + S.nav_z ] ];
    let lposData = [ [ -.5,.5,.5 ], [ 0, 0, 0 ] ];
    // S.setUniform('3fv', 'uLpos', lposData.flat());
    S.setUniform('3fv', 'uLpos', [ -.5,.5,.5,  0, 0, 0 ]);
    S.setUniform('3fv', 'uLd', ldData.flat());
    // light color abs value: intensity and color is based on the proportion of the light color
    S.setUniform('3fv', 'uLc', [ .1,.1,.1, .5,.3,.1 ]);
    // S.setUniform('3fv', 'uLc', [ .1,.1,.1, 1,.6,.05 ]);
 
    // DEFINE NUMBER OF LIGHTS FOR GPU
 
    S.nL = ldData.length;
 
    // SEND BACKGROUND COLOR TO GPU
 
    // S.setUniform('3fv', 'uBgColor', [ red.value   / 100,
    //                                   green.value / 100,
    //                                   blue.value  / 100 ]);
    S.setUniform('3fv', 'uBgColor', [ S.bg  / 100,
                                      S.bg / 100,
                                      S.bg / 100 ]);
 
    // SEND INDEX OF REFRACTION TO GPU
 
    let ior = refract.value / 100 + 1;
    S.setUniform('1f', 'uIor', ior);
 
    // DIFFERENT QUADRIC SURFACES
 
 //                xx        yy         zz           c
 
    let qSlabX  = [1,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,-1]; // x*x - 1 <= 0
    let qSlabY  = [0,0,0,0, 0,1,0,0, 0,0,0,0, 0,0,0,-1]; // y*y - 1 <= 0
    let qSlabZ  = [0,0,0,0, 0,0,0,0, 0,0,1,0, 0,0,0,-1]; // z*z - 1 <= 0
    let qSphere = [1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,-1]; // x*x + y*y + z*z - 1 <= 0
    let qTubeX  = [0,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,-1]; // y*y + z*z - 1 <= 0
    let qTubeY  = [1,0,0,0, 0,0,0,0, 0,0,1,0, 0,0,0,-1]; // x*x + z*z - 1 <= 0
    let qTubeZ  = [1,0,0,0, 0,1,0,0, 0,0,0,0, 0,0,0,-1]; // x*x + y*y - 1 <= 0
 
    // SHAPES ARE INTERSECTIONS OF QUADRIC SURFACES
 
    let shape = [], coefs = [], xform = [], phong = [], M;
 
    let sphere = (m, M) => {
       shape.push(1);
       phong.push(m);
       xform.push(M);
       coefs.push(qSphere);
    }
 
    let tubeX = (m, M) => {
       shape.push(2, 0);
       phong.push(m, m);
       xform.push(M, M);
       coefs.push(qTubeX, qSlabX);
    }
 
    let tubeY = (m, M) => {
       shape.push(2, 0);
       phong.push(m, m);
       xform.push(M, M);
       coefs.push(qTubeY, qSlabY);
    }
 
    let tubeZ = (m, M) => {
       shape.push(2, 0);
       phong.push(m, m);
       xform.push(M, M);
       coefs.push(qTubeZ, qSlabZ);
    }
 
    let cube = (m, M) => {
       shape.push(3, 0, 0);
       phong.push(m, m, m);
       xform.push(M, M, M);
       coefs.push(qSlabX, qSlabY, qSlabZ);
    }
 
    let octahedron = (m, M) => {
       shape.push(4, 0, 0, 0);
       phong.push(m, m, m, m);
       xform.push(M, M, M, M);
       coefs.push([1, 2, 2, 0,  0, 1, 2, 0,  0,0,1,0,  0,0,0,-1]);
       coefs.push([1,-2,-2, 0,  0, 1, 2, 0,  0,0,1,0,  0,0,0,-1]);
       coefs.push([1,-2, 2, 0,  0, 1,-2, 0,  0,0,1,0,  0,0,0,-1]);
       coefs.push([1, 2,-2, 0,  0, 1,-2, 0,  0,0,1,0,  0,0,0,-1]);
    }
 
    // CREATE THE SCENE
    // top
    
    // tubeY(S.redPlastic,
    //       mScale(.4,.027,.12,
    //       mRoty(0,
    //       mRotz(0,
    //       mRotx(0,
    //       matrixTranslate(0,.6,0))))));
    // middle tube
    // tubeY(S.redPlastic,
    //    mScale(.6,.02,.2,
    //    mRoty(0,
    //    mRotz(0,
    //    mRotx(0,
    //    matrixTranslate(0,.65,0))))));
    // tubeX(S.redPlastic,
    //       mScale(.6,.04,.04,
    //       mRoty(1,
    //       mRotz(1.55,
    //       mRotx(1.2,
    //       matrixTranslate(0,-.01,0))))));
    // btm
    // tubeY(S.redPlastic,
    //       mScale(.8,.04,.2,
    //       mRotx(0,
    //       mRoty(0,
    //       mRotz(0,
    //       matrixTranslate(0,-.6,0))))));
    
          // tubeJ(S.redPlastic,
          //    mScale(.4,.02,.2,
          //    mRoty(.1,
          //    mRotz(.6,
          //    mRotx(1,
          //    matrixTranslate(0,.1,0))))));
 
 
    // octahedron(S.greenPlastic,
    //     mScale(.18,.18,.18,
    //     mRoty(time * 1.2 /5,
    //     mRotz(time * 1.3 /5,
    //     mRotx(time * 1.1 /5,
    //     matrixTranslate(0,-Math.cos(time)*.2,Math.sin(time)*.4))))));
 
       // LOOP
       // sin(2PI/n*i + time)
       // cos(2PI/n*i + time)
 // Decoration
    // for (let i = 0; i < 10; i++) {
    //    let x = Math.sin(2*Math.PI/10*i + time/5);
    //    let z = Math.cos(2*Math.PI/10*i + time/5);
    //    let y = .55;
    //    let r = .3;
    //    let s = .03;
    //    octahedron(S.blackGlass,
    //       mScale(s,s,s,
    //       mRoty(time * 1.2 /5,
    //       mRotz(0,
    //       mRotx(0,
    //       matrixTranslate(x*r + S.nav_x,y-0.5 + Math.cos(time)*.02,z*r + S.nav_z))))));
 
    // }
    // octahedron(S.greenPlastic,
    //    mScale(.03,.03,.03,
    //    mRoty(time * 1.2 /5,
    //    mRotz(0,
    //    mRotx(time * 1.1 /5,
    //    matrixTranslate(Math.sin(time)*.3,.55,-Math.cos(time)*.3))))));
    
       // octahedron(S.greenPlastic,
       //    mScale(.03,.03,.03,
       //    mRoty(time * 1.2 /5,
       //    mRotz(0,
       //    mRotx(time * 1.1 /5,
       //    matrixTranslate(Math.sin(time)*.3 + .2,.55,-Math.cos(time)*.3))))));
 
       //    octahedron(S.greenPlastic,
       //       mScale(.03,.03,.03,
       //       mRoty(time * 1.2 /5,
       //       mRotz(0,
       //       mRotx(time * 1.1 /5,
       //       matrixTranslate(Math.sin(time)*.3 + .1,.55,-Math.cos(time)*.3))))));
 
    // cube(S.whitePlastic,
    //     mScale(.18,.03,.12,
    //     mRoty(time * 1.2 /5,
    //     mRotz(time * 1.1 /5,
    //     mRotx(time * 1.3 /5,
    //     matrixTranslate(0,Math.cos(time)*.2,0))))));
 
    // cube(S.whitePlastic,
    //       mScale(1.2,1.2,1.2,
    //       mRoty(time * 1.2 /5,
    //       mRotz(time * 1.1 /5,
    //       mRotx(time * 1.3 /5,
    //       matrixTranslate(0,Math.cos(time)*.2,-3))))));
 
    // sphere(S.blackGlass,
    //       mScale(.12,.04,.04,
    //       matrixTranslate(
    //          0 + S.nav_x,
    //          .2 + Math.cos(time)*.02, 
    //          1 + S.nav_z)));
 
    
    // sphere(S.blackGlass,
    //       mScale(.15,.15,.15,
    //       matrixTranslate(
    //          -0.7 + S.nav_x,
    //          0.2  + Math.cos(time)*.02,
    //          0 + S.nav_z)));
    
    // sphere(S.blackGlass,
    //       mScale(.2,.15,.18,
    //       matrixTranslate(
    //          0 + S.nav_x,
    //          .1 + Math.cos(time)*.02,
    //          -.3 + S.nav_z)));
 
    // sphere(S.gold,
    //       mScale(.2,.15,.18,
    //       matrixTranslate(S.nav_x + 0,0+Math.cos(time)*.02,S.nav_z+0)));
 
    // sphere(S.gold,
    //       mScale(.2 + .15 * Math.cos(time/10),.2 + .15 * Math.cos(time/10),.2 + .15 * Math.cos(time/10),
    //       matrixTranslate(S.nav_x + 0,0,S.nav_z + 0)));
 
    sphere(S.gold,
          mScale(.1,.1,.1,
          matrixTranslate(-.5,.6,.1)));
 // sphere(S.gold,
 //          mScale(.2,.15,.18,
 //          matrixTranslate(-.5 ,.5,.5)));
 
    // sphere(S.blackGlass,
    //       mScale(.2,.15,.18,
    //       matrixTranslate(
    //          4.3 + S.nav_x,
    //          .1 + Math.cos(time)*.02,
    //          -.3 + S.nav_z)));
 
    sphere(S.whitebright,
          mScale(.18,.18,.18,
          matrixTranslate(
             0.15 + S.nav_x,
             .2 + Math.cos(time)*.02,
             -.3 + S.nav_z)));
 
    cube(S.matteGray,
          mScale(2.9,.5,3.2,
          mRoty(.3,
          mRotz(.0,
          mRotx(0.,
          matrixTranslate(
             0 + S.nav_x,
             -.9,
             .5 + S.nav_z))))));
    
    cube(S.translucentWhite,
 
          mScale(.9,2,2.2,
          mRoty(.3,
          mRotz(.0,
          mRotx(.0,
          matrixTranslate(
             1.7 + S.nav_x,
             -.2,
             .1 + S.nav_z))))));
 
 
    // SEND SCENE DATA TO GPU
 
    for (let n = 0 ; n < coefs.length ; n++) {
       let IM = matrixInverse(xform[n]);
       coefs[n] = matrixMultiply(matrixTranspose(IM), matrixMultiply(coefs[n], IM));
    }
    S.setUniform('1iv', 'uShape', shape);
    S.setUniform('1f', 'uFl', z_value);
    S.setUniform('Matrix4fv', 'uQ', false, coefs.flat());
    S.setUniform('Matrix4fv', 'uPhong', false, phong.flat());
 
    // DEFINE NUMBER OF QUADRIC SURFACES FOR GPU
 
    S.nQ = coefs.length;
 
    // RENDER THIS ANIMATION FRAME
 
    S.gl.drawArrays(S.gl.TRIANGLE_STRIP, 0, 4);
 
    // SET ANY HTML INFO
 
    // iorInfo.innerHTML = 'index of refraction = ' + (ior * 100 >> 0) / 100;
 `,
       events: `
 
       handleXPlus=(e)=>{
          // console.log('event',e.key);
          if(S.nav_x < 1.5)
             S.nav_x += .05;
       }
 
       handleXMinus=()=>{
          if(S.nav_x > -1.5)
             S.nav_x -= .05;
       }
 
       handleZPlus=()=>{
          if(S.nav_z < 3.4)
             S.nav_z += .05;
       }
 
       handleZMinus=()=>{
          if(S.nav_z > -2)
             S.nav_z -= .05;
       }
 
       handleLight=()=>{
          if(S.lightOn){
             S.bg=12;
             S.lightOn=false;
          }else{
             S.bg=20;
             S.lightOn=true;
          }
       }
 
    ;
 `
    };
 
 }
 
 
 