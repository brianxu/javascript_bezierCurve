/* Â© 2009 ROBO Design
 * http://www.robodesign.ro
 */
if(window.addEventListener) {
    window.addEventListener('load', function () {
                            var canvas, context, tool;
                            var curves = new Array();
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
                            }
                            }
                            
                            
                            bezier_points.length = 0;
                            stroke.length = 0;
                            }
                            };
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
