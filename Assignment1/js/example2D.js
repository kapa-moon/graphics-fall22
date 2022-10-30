
rooms.example2D = function () {

  lib2D();

  description = 'Simple 2:D example';

  code = {
    'explanation': `
  S.html(\`
     A 2D canvas lets you create paths.
     <p>
     You can either
     draw <i>strokes</i> along those paths or else
     create solid shapes by <i>filling</i> those paths.
  \`);
`,
    init: `
  S.x = 400;
  S.y = 400;
`,
    assets: `
  S.line = (ax,ay,bx,by) => {
     S.context.beginPath();
     S.context.moveTo(ax,ay);
     S.context.lineTo(bx,by);
     S.context.stroke();
  }

  S.rect = (x,y,w,h) => {
     S.context.beginPath();
     S.context.rect(x,y,w,h);

     S.context.strokeStyle = 'white';
     S.context.stroke();

     if (S.isSpace) {
        S.context.fillStyle = 'gray';
        S.context.fill();
     }
  }
`,
    render: `

  let t = 3 * Math.sin(2 * Math.PI * time * time/2);

  let c = S.context;

  c.lineWidth = 10;
  c.lineCap = 'round'; 

  let wx = 300;
  let wy = 46 +t;
  
  if(time < 5){
    S.rect(wx,wy, 5,5);
  }else if(time >= 5 && time < 10){
    S.rect(wx,-325+15 *time*time, 5,5);
  }


  c.beginPath();
  c.moveTo(100,100);
  c.bezierCurveTo(0,100, 0,300, 300,50);
  c.stroke();
`,
    events: `
  onDrag = (x,y) => {
     S.x = x;
     S.y = y;
  }
  onKeyPress   = key => S.isSpace = key == 32;
  onKeyRelease = key => S.isSpace = false;
`
  };

}

