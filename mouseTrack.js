/* Â© 2009 ROBO Design
 * http://www.robodesign.ro
 */
var showControl = false;
var curves = new Array();
var context;
var canvas;
function updateShowControl()
{
    var elem = document.getElementById("button1");
    var txt = elem.textContent || elem.innerText;
    if (txt == "Show Control Point") txt = "Hide Control Point";
    else txt = "Show Control Point";
    elem.textContent = txt;
    elem.innerText = txt;
    showControl = !showControl;
    drawCurves();
}
function drawCircle(center_x, center_y, radius)
{
    context.beginPath();
    context.arc(center_x,center_y, radius, 0, 2 * Math.PI, false);
    context.fillStyle = "yellow";
    context.fill();
    context.lineWidth = 1;
    context.strokeStyle = "black";
    context.stroke();
}

function drawLine(start_x, start_y, end_x, end_y)
{
    context.lineWidth = 1;
    context.beginPath();
    context.moveTo(start_x,start_y);
    context.lineTo(end_x,end_y);
    context.stroke();
}

function drawCurves()
{
    context.clearRect(0, 0, canvas.width, canvas.height);
    for(var s = 0; s < curves.length; ++s)
    {
        var bPoints = curves[s];
        for(var i = 0; i < bPoints.length-1; i += 3)
        {
            context.beginPath();
            context.moveTo(bPoints[i].x,bPoints[i].y);
            context.bezierCurveTo(bPoints[i+1].x,bPoints[i+1].y, bPoints[i+2].x,bPoints[i+2].y, bPoints[i+3].x,bPoints[i+3].y);
            context.stroke();
            if(showControl){
                //draw control point edge
                drawLine(bPoints[i].x, bPoints[i].y, bPoints[i+1].x, bPoints[i+1].y);
                drawLine(bPoints[i+3].x, bPoints[i+3].y, bPoints[i+2].x, bPoints[i+2].y);
            }
        }
        if(showControl){
            for(var i = 0; i < bPoints.length; i++)
            {
                //draw control point
                drawCircle(bPoints[i].x, bPoints[i].y, 3);
            }
        }
    }
}

if(window.addEventListener) {
    window.addEventListener('load', function () {
                                var tool;
                                var stroke = new Array();

                                function init () {
                                    // Find the canvas element.
                                    canvas = document.getElementById('imageView');
                                    if (!canvas) {
                                        alert('Error: I cannot find the canvas element!');
                                        return;
                                    }

                                    if (!canvas.getContext) {
                                        alert('Error: no canvas.getContext!');
                                        return;
                                    }
                                    canvas.width  = 800;
                                    canvas.height = 600;
                                    // Get the 2D canvas context.
                                    context = canvas.getContext('2d');
                                    if (!context) {
                                        alert('Error: failed to getContext!');
                                        return;
                                    }

                                    // Pencil tool instance.
                                    tool = new tool_pencil();

                                    // Attach the mousedown, mousemove and mouseup event listeners.
                                    canvas.addEventListener('mousedown', ev_canvas, false);
                                    canvas.addEventListener('mousemove', ev_canvas, false);
                                    canvas.addEventListener('mouseup',   ev_canvas, false);
                                }

                                // This painting tool works like a drawing pencil which tracks the mouse
                                // movements.
                                function tool_pencil () {
                                    var tool = this;
                                    this.started = false;

                                    // This is called when you start holding down the mouse button.
                                    // This starts the pencil drawing.
                                    this.mousedown = function (ev) {
                                             context.beginPath();
                                             context.moveTo(ev._x, ev._y);
                                             stroke.push(new Point(ev._x, ev._y));
                                             tool.started = true;
                                         };

                                    // This function is called every time you move the mouse. Obviously, it only
                                    // draws if the tool.started state is set to true (when you are holding down
                                    // the mouse button).
                                    this.mousemove = function (ev) {
                                             if (tool.started) {
                                                 context.lineTo(ev._x, ev._y);
                                                 context.stroke();
                                                 stroke.push(new Point(ev._x, ev._y));
                                             }
                                         };

                                    // This is called when you release the mouse button.
                                    this.mouseup = function (ev) {
                                             if (tool.started) {
                                                 tool.mousemove(ev);
                                                 tool.started = false;
                                                 context.font = '18pt Calibri';
                                                 context.fillStyle = 'black';
                                                 //document.write(p.x + "," + p.y+"<br>");
                                                 var bezier_points = new Array();
                                                 var tolerance=10;
                                                 if(stroke.length>2)
                                                     fitBezier(stroke, bezier_points, tolerance);
                                                 curves.push(bezier_points.slice(0));

                                                 drawCurves();
                                                 stroke.length = 0;
                                                 bezier_points.length = 0;
                                             }
                                         }
                                }
                                    // The general-purpose event handler. This function just determines the mouse
                                    // position relative to the canvas element.
                                    function ev_canvas (ev) {
                                        if (ev.layerX || ev.layerX == 0) { // Firefox
                                            ev._x = ev.layerX;
                                            ev._y = ev.layerY;
                                        } else if (ev.offsetX || ev.offsetX == 0) { // Opera
                                            ev._x = ev.offsetX;
                                            ev._y = ev.offsetY;
                                        }

                                        // Call the event handler of the tool.
                                        var func = tool[ev.type];
                                        if (func) {
                                            func(ev);
                                        }
                                    }

                                    init();

                                }, false);
                            }

                            // vim:set spell spl=en fo=wan1croql tw=80 ts=2 sw=2 sts=2 sta et ai cin fenc=utf-8 ff=unix:
