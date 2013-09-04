var pixelChain = new Array();
var bezierControlPoints = new Array();
var lon = new Array();
var polLine = new Array();
var INFTY = 10000000
/* ---------------------------------------------------------------------- */
/*
 *  B0, B1, B2, B3 : Bezier multipliers
 */
function B0(u) { return ( 1.0 - u )  *  ( 1.0 - u )  *  ( 1.0 - u ); }
function B1(u) { return 3 * u  *  ( 1.0 - u )  *  ( 1.0 - u ); }
function B2(u) { return 3 * u * u  *  ( 1.0 - u ); }
function B3(u) { return u * u * u; }

function Point(x, y) {
    this.x = x;
    this.y = y;
    this.norm = function() { return Math.sqrt(this.x*this.x+this.y*this.y);}
    this.plus = function(p) { return new Point(this.x+p.x, this.y+p.y);}
    this.minus = function(p) { return new Point(this.x-p.x, this.y-p.y);}
    this.is_zero = function() { return this.x == 0 && this.y == 0;}
    this.neg = function() {return new Point(-this.x, -this.y);}
    this.equal = function(p){return this.x == p.x && this.y == p.y;}
    this.time = function(a){return new Point(this.x*a, this.y*a);}
    this.normalize = function(){var normV = this.norm();this.x /= normV; this.y /= normV;}
    this.copy = function(){return new Point(this.x, this.y);}
    this.set = function(p){this.x == p.x;this.y == p.y;}
}

function fitBezier(freeHandDrawnPath, bezier_points, tolerance)
{
    pixelChain.length = 0;
    bezierControlPoints.length = 0;
    lon.length = 0;
    polLine.length = 0;
    //document.write("---Debug inside fitBezier")
    var points_size = freeHandDrawnPath.length
    //document.write(" size:" + points_size + "<br>")
    
    //----
    
    var p0, p;
	var closedCurve = false;
	var dist_x, dist_y;
	var len = freeHandDrawnPath.length;
	var i;
    
	p0 = freeHandDrawnPath[0];

    //document.write("len:" + len + "<br>");
	for(i=1;i<len;i++)
	{
		p = freeHandDrawnPath[i];
		bresenham(p0, p);
		p0 = p;
	}
	//check for duplicate points and eliminate them
	for(var i=1; i<pixelChain.length; )
    {
        if(pixelChain[i].x == pixelChain[i-1].x && pixelChain[i].y == pixelChain[i-1].y)
        {
            pixelChain.splice(i,1);
        } else {
            i++;
        }
    }
    
    
	//check to see whether it's a closed chain
    len = pixelChain.length;
	dist_x = Math.abs(pixelChain[0].x - pixelChain[len-1].x);
	dist_y = Math.abs(pixelChain[0].y - pixelChain[len-1].y);
    
	if(dist_x<=2 && dist_y<=2)
	{
		var p = pixelChain[0];
		pixelChain[len-1].x = p.x;
		pixelChain[len-1].y = p.y;
		closedCurve = true;
	}
//    for(var i = 0; i < len; i++)
//    {
//        document.write("pixelChain:(" + pixelChain[i].x + "," + pixelChain[i].y + ")<br>");
//    }
    
    calc_lon();
	bestpolygon();
    smooth(1, tolerance);
    
    len = bezierControlPoints.length;
    //document.write("len:(" + len + ")<br>");
	for(i=0;i<len;i++)
	{
        var p = new Point(bezierControlPoints[i].x, bezierControlPoints[i].y);
		bezier_points.push(p);
        //document.write("bezier_points:[" +i + "]=("+ bezier_points[i].x + "," + bezier_points[i].y + ")<br>");
	}
    
	return closedCurve;
}

function bresenham(p0, p1)
{
    
    var x0 = p0.x;
	var y0 = p0.y;
	var x1 = p1.x;
	var y1 = p1.y;
	var computed_pixels = new Array();
	var computed_pixels_size;
    
	var steep = Math.abs(y1 - y0) > Math.abs(x1 - x0);
	var swap_points = false;
	if (steep)
	{
        y0 = [x0, x0 = y0][0];//swap x0, y0
		y1 = [x1, x1 = y1][0];//swap x1, y1
	}
	if (x0 > x1)
	{
		swap_points = true;
        x1 = [x0, x0 = x1][0];//swap x0, x1
		y1 = [y0, y0 = y1][0];//swap y0, y1
	}
    
	var deltax = x1 - x0;
	var deltay = Math.abs(y1 - y0);
	var error = parseInt(-(deltax + 1) / 2, 10);
	var ystep;
	var y = y0;
	if (y0 < y1)
		ystep = 1;
	else
		ystep = -1;
    
//    document.write("x0y0:(" + x0 + "," + y0 + ")<br>");
//    document.write("x1y1:(" + x1 + "," + y1 + ")<br>");
//    document.write("steep:(" + steep + ")<br>");

	for (var x = x0; x <= x1; x++)
	{
        var p = new Point(0, 0);
		if (steep){
			p.x = y;
			p.y = x;
			computed_pixels.push(p);
		}
        else {
			p.x = x;
			p.y = y;
			computed_pixels.push(p);
		}
		error = error + deltay;
		if (error >= 0){
			y = y + ystep;
			error = error - deltax;
		}
        computed_pixels_size = computed_pixels.length;
	}
    computed_pixels_size = computed_pixels.length;
	
	if(swap_points)
	{
		//invert the pixels
		for (var i=computed_pixels_size-1; i>=0;i--)
		{
			pixelChain.push(computed_pixels[i]);
            //document.write("pixelChain1:(" + computed_pixels[i].x + "," + computed_pixels[i].y + ")<br>");
		}
	}
	else {
		for (var i=0; i<computed_pixels_size;i++)
		{
			pixelChain.push(computed_pixels[i]);
		}
	}
}

//POTRACE
/* ---------------------------------------------------------------------- */
/* Stage 1: determine the straight subpaths (Sec. 2.2.1). Fill in the
 "lon" component of a path object (based on pt/len).	For each i,
 lon[i] is the furthest index such that a straight line can be drawn
 from i to lon[i]. Return 1 on error with errno set, else 0. */

/* this algorithm depends on the fact that the existence of straight
 subpaths is a triplewise property. I.e., there exists a straight
 line through squares i0,...,in iff there exists a straight line
 through i,j,k, for all i0<=i<j<k<=in. (Proof?) */

/* this implementation of calc_lon is O(n^2). It replaces an older
 O(n^3) version. A "constraint" means that future points must
 satisfy xprod(constraint[0], cur) >= 0 and xprod(constraint[1],
 cur) <= 0. */

/* Remark for Potrace 1.1: the current implementation of calc_lon is
 more complex than the implementation found in Potrace 1.0, but it
 is considerably faster. The introduction of the "nc" data structure
 means that we only have to test the constraints for "corner"
 points. On a typical input file, this speeds up the calc_lon
 function by a factor of 31.2, thereby decreasing its time share
 within the overall Potrace algorithm from 72.6% to 7.82%, and
 speeding up the overall algorithm by a factor of 3.36. On another
 input file, calc_lon was sped up by a factor of 6.7, decreasing its
 time share from 51.4% to 13.61%, and speeding up the overall
 algorithm by a factor of 1.78. In any case, the savings are
 substantial. */

function calc_lon()
{
	var len = pixelChain.length;
    
	var i;
	var j;
	var k, k1, dir, dir_index;
	var ct = new Array();
	var constraint = new Array();
	var cur;
	var off;
	var p_dir;
	var pivk = new Array();
    var nc = new Array(); /* nc[len]: next corner */
	var dk;  /* direction of k-k1 */
	var a, b, c, d;
    
	/* initialize the nc data structure. Point from each point to the
     furthest future point to which it is connected by a vertical or
     horizontal segment or is the next neighbour  */
	k = len-1;
	for (i=len-2; i>=0; i--) {
		if (pixelChain[i].x != pixelChain[k].x && pixelChain[i].y != pixelChain[k].y) {
			k = i+1;  /* necessarily i<n-1 in this case */
		}
		nc[i] = k;
        //document.write("nc:" + nc[i] + "<br>");
	}
	nc[len-1] = len-1;
	pivk[len-1] = len-1;
	/* determine pivot points: for each i, let pivk[i] be the furthest k
     such that all j with i<j<k lie on a line connecting i,k. */
	var direction = new Array();
	direction[0] = 1; direction[1] = 0; direction[2] = 7;
	direction[3] = 2; direction[4] = 8; direction[5] = 6;
	direction[6] = 3; direction[7] = 4; direction[8] = 5;
    
	for (i=len-2; i>=0; i--) {
        
		ct[0] = ct[1] = ct[2] = ct[3] = 0;
        
        p_dir = pixelChain[i+1].minus(pixelChain[i]);
		/* keep track of "directions" that have occurred */
		p_dir = ddir(p_dir);
		dir_index= (4+3*(p_dir.x)+p_dir.y);
        //document.write("p_dir:" + p_dir.x + " " + p_dir.y +"<br>");
		dir = parseInt((direction[dir_index]+0.5)/2) - 1;
        //document.write("dir:" + dir +"<br>");
		ct[dir]=ct[dir]+1;
        
		if(direction[dir_index]%2>0)
		{
			ct[mod(dir+1,4)]=ct[mod(dir+1,4)]+1;
		}
        
        constraint[0] = new Point(0,0);
        constraint[1] = new Point(0,0);
        
		/* find the next k such that no straight line from i to k */
        
		k = nc[i];
		k1 = i;
        var gotoLabel = -1;
		while (1) {
            p_dir = pixelChain[i+1].minus(pixelChain[i]);
			p_dir = ddir(p_dir);
			dir_index= (4+3*p_dir.x+p_dir.y);
			dir = parseInt((direction[dir_index]+0.5)/2) - 1;
            
            //document.write("dir2:" + dir +"<br>");
			ct[dir]++;
			if(direction[dir_index]%2>0)
			{
				ct[mod(dir+1,4)]++;
			}
			/* if all four "directions" have occurred, cut this path */
			if (ct[0] && ct[1] && ct[2] && ct[3]) {
				pivk[i] = k1;
				//goto foundk;
                gotoLabel = 2;
                break;
			}
            
            cur = pixelChain[k].minus(pixelChain[i]);
            
			/* see if current constraint is violated */
			if (xprod(constraint[0], cur) < 0 || xprod(constraint[1], cur) > 0) {
				//goto constraint_viol;
                gotoLabel = 1;
                break;
			}
            
			/* else, update constraint */
			if (Math.abs(cur.x) <= 1 && Math.abs(cur.y) <= 1) {
				/* no constraint */
			} else {
                off = new Point();
				off.x = cur.x + ((cur.y>=0 && (cur.y>0 || cur.x<0)) ? 1 : -1);
				off.y = cur.y + ((cur.x<=0 && (cur.x<0 || cur.y<0)) ? 1 : -1);
				if (xprod(constraint[0], off) >= 0) {
					constraint[0] = off;
				}
                off = new Point();
				off.x = cur.x + ((cur.y<=0 && (cur.y<0 || cur.x<0)) ? 1 : -1);
				off.y = cur.y + ((cur.x>=0 && (cur.x>0 || cur.y<0)) ? 1 : -1);
				if (xprod(constraint[1], off) <= 0) {
					constraint[1] = off;
				}
			}
			k1 = k;
			k = nc[k1];
			if (!cyclic(k,i,k1)) {
                gotoLabel = 1;
				break;
			}
		} //end while
        //document.write("gotoLabel:" + gotoLabel +"<br>");
        switch(gotoLabel)
        {
            case 1://constraint_viol
            {
                /* k1 was the last "corner" satisfying the current constraint, and
                 k is the first one violating it. We now need to find the last
                 point along k1..k which satisfied the constraint. */
                dk = new Point();
                cur = new Point();
                dk.x = sign(pixelChain[k].x-pixelChain[k1].x);
                dk.y = sign(pixelChain[k].y-pixelChain[k1].y);
                cur.x = pixelChain[k1].x - pixelChain[i].x;
                cur.y = pixelChain[k1].y - pixelChain[i].y;
                
                /* find largest integer j such that xprod(constraint[0], cur+j*dk)
                 >= 0 and xprod(constraint[1], cur+j*dk) <= 0. Use bilinearity
                 of xprod. */
                a = xprod(constraint[0], cur);
                b = xprod(constraint[0], dk);
                c = xprod(constraint[1], cur);
                d = xprod(constraint[1], dk);
                /* find largest integer j such that a+j*b>=0 and c+j*d<=0. This
                 can be solved with integer arithmetic. */
                j = INFTY;
                if (b<0) {
                    j = floordiv(a,-b);
                }
                if (d>0) {
                    j = Math.min(j, floordiv(-c,d));
                }
                pivk[i] = mod(k1+j,len);
                break;
            }
            case 2://foundk
                break;
        }
        
        //document.write("pivk: " + pivk[i] + "<br>");
        
	} 
    
	j=pivk[len-1];
	lon[len-1]=j;
	i=len-2;
	for (i=len-2; i>=0; i--) {
		if (cyclic(i+1,pivk[i],j)) {
			j=pivk[i];
		}
		lon[i]=j;
	}
    
    for(var i=0; i<lon.length;++i)
    {
        //document.write("lon " + i + " :" + lon[i] + "<br>");
    }
    
}
/*********************************************************************************************/
/* ---------------------------------------------------------------------- */
/* Stage 2: calculate the optimal polygon (Sec. 2.2.2-2.2.4). */

/* Auxiliary function: calculate the penalty of an edge from i to j in
 the given path. This needs the "lon" and "sum*" data. */

function penalty3(i, j) {
    
	/* assume 0<=i<j<=n  */
	var len = pixelChain.length;
	var x, y, x0, y0, x1, y1;
	var AB, sum;
    
    
	if(i==j) return 0;
    
	if(i>j)
	{
		var tmp = j;
		j = i; i = tmp;
	}
	if(j>len-1)
		j = len-1;
    
	var A, B, C;
	x0 = pixelChain[i].x;
    y0 = pixelChain[i].y;
	x1 = pixelChain[j].x;
    y1 = pixelChain[j].y;
	A = y0 - y1;
	B = x1 - x0;
	C = x0*y1 - x1*y0;
    
	AB = Math.sqrt(A*A+B*B);
	sum = 0;
	for(var k=i;k<=j;k++)
	{
		x = pixelChain[k].x;
        y = pixelChain[k].y;
		sum +=Math.abs(A*x+B*y+C)/AB;
	}
	sum = AB*Math.sqrt(sum/(j-i+1));
	return sum;
}

/* find the optimal polygon. Fill in the m and po components. Return 1
 on failure with errno set, else 0. Non-cyclic version: assumes i=0
 is in the polygon. Fixme: ### implement cyclic version. */
function bestpolygon()
{
    var i, j, m, k;
	var n = pixelChain.length;
	var pen = new Array(); // pen[n+1]: penalty vector
	var prev = new Array();   // prev[n+1]: best path pointer vector
	var clip0 = new Array();  // clip0[n]: longest segment pointer, non-cyclic
	var clip1 = new Array();  // clip1[n+1]: backwards segment pointer, non-cyclic
	var seg0 = new Array();    // seg0[m+1]: forward segment bounds, m<=n
	var seg1 = new Array();   // seg1[m+1]: backward segment bounds, m<=n
	var thispen;
	var best;
	var c;
    
	// calculate clipped paths
	for (i=0; i<n; i++) {
		c = lon[i];
		if (c == i) {
			c = mod(i+1,n);
		}
		if (c < i) {
			clip0[i] = n;
		} else {
			clip0[i] = c;
		}
	}
    
	// calculate backwards path clipping, non-cyclic. j <= clip0[i] if
	//clip1[j] <= i, for i,j=0..n.
	j = 1;
	for (i=0; i<n; i++) {
		while (j <= clip0[i]) {
			clip1[j] = i;
			j++;
		}
	}
    
	// calculate seg0[j] = longest path from 0 with j segments
	i = 0;
	for (j=0; i<n; j++) {
		seg0[j] = i;
		i = clip0[i];
	}
	seg0[j] = n;
	m = j;
	// calculate seg1[j] = longest path to n with m-j segments
	i = n;
	for (j=m; j>0; j--) {
		seg1[j] = i;
		i = clip1[i];
	}
	seg1[0] = 0;
    
	// now find the shortest path with m segments, based on penalty3
	// note: the outer 2 loops jointly have at most n interations, thus
	//the worst-case behavior here is quadratic. In practice, it is
	//close to linear since the inner loop tends to be short.
	pen[0]=0;
	for (j=1; j<=m; j++) {
		for (i=seg1[j]; i<=seg0[j]; i++) {
			best = -1;
			for (k=seg0[j-1]; k>=clip1[i]; k--) {
				thispen = penalty3(k, i) + pen[k];
				if (best < 0 || thispen < best) {
					prev[i] = k;
					best = thispen;
				}
			}
			pen[i] = best;
		}
	}
    
	// read off shortest path
	for (i=n, j=m-1; i>0; j--) {
		i = prev[i];
		polLine[j] = i;
        //document.write("polline:" + i + "<br>");
	}
}

function smooth(zoom, tolerance)
{
    var alphamax = 2;
    var m = polLine.length;
    
	var i, j, k;
	var dd, denom, alpha;
    
	var continuousCurve = new Array();
	var isFirst = true;
	var isLast = false;
	//for (i=0; i<m; i++) {
	//	polLineCorners.push_back(CURVE);
	//}
    
	/* examine each vertex and find its best fit */
	continuousCurve.push(pixelChain[polLine[0]]);
	for (i=0; i<m-2; i++) {
		j = mod(i+1, m);
		k = mod(i+2, m);
        
		denom = ddenom(pixelChain[polLine[i]], pixelChain[polLine[k]]);
		if (denom != 0.0) {
			dd = dpara(pixelChain[polLine[i]], pixelChain[polLine[j]], pixelChain[polLine[k]]) / denom;
			dd = Math.abs(dd);
			alpha = dd>1 ? (1 - 1.0/dd) : 0;
			alpha = alpha / 0.75;
		}
		else {
			alpha = 4/3.0;
		}
		//curve->alpha0[j] = alpha;	 /* remember "original" value of alpha */
        
		if (alpha >= alphamax) {  /* pointed corner */
			//polLineCorners[j] = LINE;
			continuousCurve.push(pixelChain[polLine[j]]);
			fillBezierPoints(continuousCurve, zoom, tolerance, isFirst, isLast);
			isFirst = false;
			continuousCurve.length = 0;
			continuousCurve.push(pixelChain[polLine[j]]);
		}
		else {
			//polLineCorners[j] = CURVE;
			continuousCurve.push(pixelChain[polLine[j]]);
		}
	}//end for(i<m-2)
	//add the last point
	if(continuousCurve.length == 0) //does this really happens?
	{
		continuousCurve.push(pixelChain[polLine[m-2]]);
		continuousCurve.push(pixelChain[polLine[m-1]]);
	}
	else
		continuousCurve.push(pixelChain[polLine[m-1]]);
	isLast = true;
	fillBezierPoints(continuousCurve, zoom, tolerance, isFirst, isLast);
	return 0;
}

//fit Beziers that pass through all points of the polyline (twice as many points)
function fillBezierPoints(continuousCurve, zoom, tolerance, isFirst, isLast)
{
	var curveLen = continuousCurve.length;
    for(var i = 0; i < curveLen; ++i)
    {
        //document.write("continuousCurve[" + i + "]=", continuousCurve[i].x + "," + continuousCurve[i].y + "<br>");
    }
    
	//int Max_Beziers = (curveLen+3)/4;
	var Max_Beziers = curveLen; //one bezier for each couple of points
	var bezier = new Array();
	var tolerance_sq = Math.abs(zoom)*tolerance;
	var tHat1 = new Point(0,0);
	var tHat2 = new Point(0,0);
	//int split_points[];
	var nsegs;
    
	nsegs = bezier_fit_cubic(bezier, null, continuousCurve, curveLen, tHat1, tHat2, tolerance_sq, Max_Beziers);
    //document.write("bezier length="+bezier.length +"<br>");
    //document.write("nsegs ="+nsegs +"<br>");
//    for(var i=0;i<bezier.length;i++)
//    {
//        document.write("bezier(" + bezier[i].x + "," + bezier[i].y + ")<br>");
//    }
	for(var i=0;i<nsegs;i++)
	{
		var start = i*4;
		/*std::cerr<<"bezier["<<start<<"] = ("<<bezier[start].getX()<<"; "<<bezier[start].getY()<<"); ";
         std::cerr<<"bezier["<<start+1<<"] = ("<<bezier[start+1].getX()<<"; "<<bezier[start+1].getY()<<"); ";
         std::cerr<<"bezier["<<start+2<<"] = ("<<bezier[start+2].getX()<<"; "<<bezier[start+2].getY()<<"); ";
         std::cerr<<"bezier["<<start+3<<"] = ("<<bezier[start+3].getX()<<"; "<<bezier[start+3].getY()<<"); \n";*/
        
		bezierControlPoints.push(bezier[start].copy());
        //document.write("bezier[" + bezier[start].x + "," + bezier[start].y + "]<br>");
		bezierControlPoints.push(bezier[start+1].copy());
        //document.write("bezier[" + bezier[start+1].x + "," + bezier[start+1].y + "]<br>");
		bezierControlPoints.push(bezier[start+2].copy());
        //document.write("bezier[" + bezier[start+2].x + "," + bezier[start+2].y + "]<br>");
	}
	if(isLast && nsegs>0)
	{
		bezierControlPoints.push(bezier[nsegs*4-1]);
	}
    
}

/**
 * Fit a multi-segment Bezier curve to a set of digitized points, without
 * possible weedout of identical points and NaNs.
 *
 * \pre data is uniqued, i.e. not exist i: data[i] == data[i + 1].
 * \param max_beziers Maximum number of generated segments
 * \param Result array, must be large enough for n. segments * 4 elements.
 */
function bezier_fit_cubic(bezier, split_points, data, len, tHat1, tHat2, error, max_beziers)
{
//    document.write("bezier len1:" + bezier.length + "<br>")
//    for(var i = 0; i < data.length; ++i)
//    {
//        document.write("data("+data[i].x+ ","+data[i].y + ")<br>");
//    }
//    document.write("len:"+len+"<br>");
	//Point bezier[4];
	//QVector<Point> data;
	var u = new Array();
	//Point tHat1 = unconstrained_tangent; Point tHat2 = unconstrained_tangent;
	var dist;
	var difPoint;
	var maxIterations = 4;   /* Max times to try iterating */
	var nsegs1, nsegs2;
	
	if ( len < 2 ) {
		return 0;}
    
	if ( len == 2 ) {
        bezier[1] = new Point(0,0);
        bezier[2] = new Point(0,0);
        
		/* We have 2 points, which can be fitted by a line segment. */
		bezier[0] = data[0].copy();
		bezier[3] = data[len - 1].copy();
		/* Straight line segment. */
        var dPx, dPy;
        dPx = bezier[3].x - bezier[0].x;
        dPy = bezier[3].y - bezier[0].y;
		difPoint = new Point(dPx, dPy);
        
        bezier[1].x = bezier[2].x = (bezier[0].x + bezier[3].x)/2;
        bezier[1].y = bezier[2].y = bezier[0].y + (difPoint.y/(difPoint.x != 0 ? difPoint.x : 1))*(bezier[1].x - bezier[0].x);
        
		//curved Bezier
		/*
         dist = (difPoint.L2_norm2D() / 3.0 );
         bezier[1] = ( tHat1.is_zero()
         ? ( bezier[0] * 2 + bezier[3] ) / 3.
         : bezier[0] + tHat1 * dist);
         bezier[2] = ( tHat2.is_zero()
         ? ( bezier[0] + bezier[3] * 2 ) / 3.
         : bezier[3] + tHat2 * dist );*/
        
		return 1;
	}
    
	/*  Parameterize points, and attempt to fit curve */
	var splitPoint = new Array();   /* Point to split point set at. */
	var is_corner;
	{
		chord_length_parameterize(data, u, len);
		if ( u[len - 1] == 0.0 ) {
			/* Zero-length path: every point in data[] is the same.
             *
             * (Clients aren't allowed to pass such data; handling the case is defensive
             * programming.)
             */
			u.length = 0;
			return 0;
		}
        
        
		//std::cerr<<"\n        tHat1 = ("<<tHat1.x<<", "<<tHat1.y<<");\n";
		//std::cerr<<"        tHat2 = ("<<tHat2.x<<", "<<tHat2.y<<");\n\n";
		generate_bezier(bezier, data, u, len, tHat1, tHat2, error);
		reparameterize(data, len, u, bezier);
        
		/* Find max deviation of points to fitted curve. */
		var tolerance = 1;//sqrt(error + 1e-9);
		var maxErrorRatio = compute_max_error_ratio(data, u, len, bezier, tolerance, splitPoint);
        //document.write("bezier_fit_cubic splitPoint = " + splitPoint[0] + "<br>");
		if ( Math.abs(maxErrorRatio) <= error ) {
			u.length = 0;
			return 1;
		}
        
		/* If error not too large, then try some reparameterization and iteration. */
		if ( 0.0 <= maxErrorRatio && maxErrorRatio <= 3.0*tolerance ) {
			for (var i = 0; i < maxIterations; i++) {
				generate_bezier(bezier, data, u, len, tHat1, tHat2, error);
				reparameterize(data, len, u, bezier);
				maxErrorRatio = compute_max_error_ratio(data, u, len, bezier, tolerance, splitPoint);
				if ( Math.abs(maxErrorRatio) <= error ) {
					u.length=0;
					return 1;
				}
			}
		}
        
		u.length = 0;
		is_corner = (maxErrorRatio < 0);
	}
    //document.write("is_corner:"+is_corner+"<br>");
	if (is_corner) {
		//std::cerr<<"\nis corner\n";
		if (splitPoint[0] == 0) {
			if (tHat1.is_zero()) {
				/* Got spike even with unconstrained initial tangent. */
				++splitPoint[0];
			} else {
				var val = bezier_fit_cubic(bezier, split_points, data, len, {x:0,y:0}, tHat2, error, max_beziers);
				return val;
			}
		} else if (splitPoint[0] == len - 1) {
			if (tHat2.is_zero()) {
				/* Got spike even with unconstrained final tangent. */
				--splitPoint[0];
			} else {
				var val = bezier_fit_cubic(bezier, split_points, data, len, tHat1, {x:0,y:0},
                                       error, max_beziers);
				return val;
			}
		}
	}
	//return 1;
	if ( 1 < max_beziers ) {
		//std::cerr<<"RECURSIVE\n";
		/*
         *  Fitting failed -- split at max error point and fit recursively
         */
		var rec_max_beziers1 = splitPoint[0];
        
		var recTHat2, recTHat1;
		if (is_corner) {
			if(!(0 < splitPoint[0] && splitPoint[0] < (len - 1)))
			{
				return -1;
			}
			recTHat1 = new Point(0,0);
            recTHat2 = new Point(0,0);
		}
		else {
			/* Unit tangent vector at splitPoint. */
			recTHat2 = estimate_center_tangent(data, splitPoint[0], len);
			recTHat1 = recTHat2.neg();
		}
        var bezier1 = new Array();
		nsegs1 = bezier_fit_cubic(bezier1, split_points, data, splitPoint[0] + 1,
                                  tHat1, recTHat2, error, rec_max_beziers1);
		/*std::cerr<<"        nsegs1 = "<<nsegs1<<"; rec_max_beziers1 = "<<rec_max_beziers1<<"\n";*/
        //document.write("bezier nsegs1 len:" + bezier1.length + "<br>")
		if ( nsegs1 < 0 ) {
			return -1;
		}
		//assert( nsegs1 != 0 );
		if (split_points != null) {
			split_points[nsegs1 - 1] = splitPoint[0];
		}
		var rec_max_beziers2 = max_beziers - splitPoint[0] -1;
        
		var data_temp = new Array();
        for(var i=0; i<data.length;++i)
        {
            data_temp.push(data[i].copy());
        }
        data_temp.splice(0, splitPoint[0]);
        //document.write("nsegs2:"+nsegs1*4 + "<br>");
        var bezier2 = new Array();
		var nsegs2 = bezier_fit_cubic(bezier2,
                                      ( split_points == null ? null : split_points + nsegs1 ),
                                      data_temp, len - splitPoint[0],
                                      recTHat1, tHat2, error, rec_max_beziers2);
        
        //document.write("bezier nsegs2 len:" + bezier.length + "<br>")
        for(var i = 0; i < bezier1.length; ++i)
        {
            bezier[i] = new Point(bezier1[i].x, bezier1[i].y);
        }
        for(var i = 0; i < bezier2.length; ++i)
        {
            bezier[bezier1.length+i] = new Point(bezier2[i].x, bezier2[i].y);
        }
        //bezier = bezier1.concat(bezier2);
//        for(var i = 0; i < bezier2.length; ++i)
//        {
//            document.write("->>bezier2:" + bezier2[i].x +","+bezier2[i].y+"<br>");
//        }
		//std::cerr<<"nsegs2 = "<<nsegs2<<"; rec_max_beziers2 = "<<rec_max_beziers2<<"\n";
		if ( nsegs2 < 0 ) {
			return -1;
		}
		return nsegs1 + nsegs2;
        
	}
	else {
		return 1;
	}
}

/**
 * Given set of points and their parameterization, try to find a better assignment of parameter
 * values for the points on a 4-point Bezier.
 *
 *  \param d  Array of digitized points.
 *  \param u  Current parameter values.
 *  \param bezCurve  Current 4-point fitted curve.
 *  \param len  Number of values in both d and u arrays.
 *              Also the size of the array that is allocated for return.
 */
function reparameterize(d, len, u, bezCurve)
{
	//assert( 2 <= len );
    
	var last = len - 1;
	//assert( u[0] == 0.0 );
	//assert( u[last] == 1.0 );
	/* Otherwise, consider including 0 and last in the below loop. */
    
	for (var i = 1; i < last; i++) {
		u[i] = NewtonRaphsonRootFind(bezCurve, d[i], u[i]);
	}
}

/**
 *  Use Newton-Raphson iteration to find better root.
 *
 *  \param Q  Current fitted curve
 *  \param P  Digitized point
 *  \param u  Parameter value for "P"
 *
 *  \return Improved u
 */
function NewtonRaphsonRootFind(Q, P, u)
{
	//assert( 0.0 <= u );
	//assert( u <= 1.0 );
    
	/* Generate control vertices for Q'. */
	var Q1 = new Array();
	for (var i = 0; i < 3; i++) {
        Q1[i] = new Point((Q[i+1].x-Q[i].x)*3,(Q[i+1].y-Q[i].y)*3);
		//Q1[i] = Q[i+1].minus(Q[i]).time(3.0);
	}
    
	/* Generate control vertices for Q''. */
	var Q2 = new Array();;
	for (var i = 0; i < 2; i++) {
        Q2[i] = new Point(0,0);
		Q2[i] = Q1[i+1].minus(Q1[i]).time(2.0);
	}
    
	/* Compute Q(u), Q'(u) and Q''(u). */
	var Q_u  = bezier_pt(3, Q, u);
	var Q1_u = bezier_pt(2, Q1, u);
	var Q2_u = bezier_pt(1, Q2, u);
    
	/* Compute f(u)/f'(u), where f is the derivative wrt u of distsq(u) = 0.5 * the square of the
     distance from P to Q(u).  Here we're using Newton-Raphson to find a stationary point in the
     distsq(u), hopefully corresponding to a local minimum in distsq (and hence a local minimum
     distance from P to Q(u)). */
	var diff = Q_u.minus(P);
	var numerator = dot(diff, Q1_u);
	var denominator = dot(Q1_u, Q1_u) + dot(diff, Q2_u);
    
	var improved_u;
	if ( denominator > 0. ) {
		/* One iteration of Newton-Raphson:
         improved_u = u - f(u)/f'(u) */
		improved_u = u - ( numerator / denominator );
	} else {
		/* Using Newton-Raphson would move in the wrong direction (towards a local maximum rather
         than local minimum), so we move an arbitrary amount in the right direction. */
		if ( numerator > 0. ) {
			improved_u = u * .98 - .01;
		} else if ( numerator < 0. ) {
			/* Deliberately asymmetrical, to reduce the chance of cycling. */
			improved_u = .031 + u * .98;
		} else {
			improved_u = u;
		}
	}
    
    if (isNaN(improved_u) || improved_u>=INFTY) {
		improved_u = u;
	} else if ( improved_u < 0.0 ) {
		improved_u = 0.0;
	} else if ( improved_u > 1.0 ) {
		improved_u = 1.0;
	}
    
	/* Ensure that improved_u isn't actually worse. */
	{
		var diff_lensq = dot(diff,diff);
		for (var proportion = .125; ; proportion += .125) {
			var temp = bezier_pt(3, Q, improved_u).minus(P) ;
			if ( dot(temp, temp)  > diff_lensq ) {
				if ( proportion > 1.0 ) {
					//g_warning("found proportion %g", proportion);
					improved_u = u;
					break;
				}
				improved_u = ( ( 1 - proportion ) * improved_u  +
                              proportion         * u            );
			} else {
				break;
			}
		}
	}
    
    
	return improved_u;
}

//////////////////////////////////////////////////////////////////////////////////
/**
 * Fill in bezier[] based on the given data and tangent requirements, using
 * a least-squares fit.
 *
 * Each of tHat1 and tHat2 should be either a zero vector or a unit vector.
 * If it is zero, then bezier[1 or 2] is estimated without constraint; otherwise,
 * it bezier[1 or 2] is placed in the specified direction from bezier[0 or 3].
 *
 * \param tolerance_sq Used only for an initial guess as to tangent directions
 *   when tHat1 or tHat2 is zero.
 */
function generate_bezier(bezier, data, u, len, tHat1, tHat2, tolerance_sq)
{
	var est1 = (tHat1.x == 0 && tHat1.y == 0);
	var est2 = (tHat2.x == 0 && tHat2.y == 0);
	var est_tHat1, est_tHat2;
	est_tHat1 =  est1 ? estimate_left_tangent(data, len, tolerance_sq):tHat1;
	est_tHat2 =  est2 ? estimate_right_tangent(data, len, tolerance_sq):tHat2;
    //document.write("est_tHat1:" + est_tHat1.x +","+est_tHat1.y +"<br>");
    //document.write("est_tHat2:" + est_tHat2.x +","+est_tHat2.y +"<br>");
	estimate_lengths(bezier, data, u, len, est_tHat1, est_tHat2);
    
    
    
	/* We find that sp_darray_right_tangent tends to produce better results
     for our current freehand tool than full estimation. */
	if (est1) {
		estimate_bi(bezier, 1, data, u, len);
		if (!bezier[1].equal(bezier[0])) {
			est_tHat1 = new Point((bezier[1].x - bezier[0].x), (bezier[1].y - bezier[0].y));
            est_tHat1.normalize();
		}
		estimate_lengths(bezier, data, u,  len, est_tHat1, est_tHat2);
	}
    
//    for(var i = 0; i < 4; ++i)
//    {
//        //document.write("generate_bezier 4 ("+bezier[i].x+ ","+bezier[i].y + ")<br> ");
//    }
}

/**
 * Estimates the (backward) tangent at d[center], by averaging the two
 * segments connected to d[center] (and then normalizing the result).
 *
 * \note The tangent is "backwards", i.e. it is with respect to
 * decreasing index rather than increasing index.
 *
 * \pre (0 \< center \< len - 1) and d is uniqued (at least in
 * the immediate vicinity of \a center).
 */
function estimate_center_tangent(d, center, len)
{
	//assert( center != 0 );
	//assert( center < len - 1 );
	var ret;
	var P1 = d[center + 1].copy();
	var P2 = d[center - 1].copy();
    
	if( P1.equal(P2))
	{
		/* Rotate 90 degrees in an arbitrary direction. */
		var diff = d[center].minus(d[center - 1]);
		ret.x = -diff.y;
		ret.y = diff.x;
	}
	else
	{
		ret = d[center - 1].minus(d[center + 1]);
	}
	ret.normalize();
	return ret;
}

/**
 * Estimate the (forward) tangent at point d[0].
 *
 * Unlike the center and right versions, this calculates the tangent in
 * the way one might expect, i.e., wrt increasing index into d.
 *
 * \pre 2 \<= len.
 * \pre d[0] != d[1].
 * \pre all[p in d] in_svg_plane(p).
 * \post is_unit_vector(ret).
 **/

function estimate_left_tangent(d, len, tolerance_sq)
{
	//assert( 2 <= len );
	//assert( 0 <= tolerance_sq );
    
	//tolerance_sq is always 0!!!)
	//return estimate_left_tangent(d, len);
	//std::cerr<<"estimate_left_tangent COMPUTE TANGENT LEFT: len = "<<len<<"; tolerance_sq = "<<tolerance_sq<<"\n";
    var p, pi, t;
	var distsq;
	for (var i = 1;;) {
		pi = d[i].copy();
		t = new Point(pi.x - d[0].x, pi.y - d[0].y);
		distsq = dot(t, t);
		if ( tolerance_sq < distsq ) {
			p = t.copy();
            p.normalize();
			return p;
		}
		++i;
		if (i == len) {
			p = t.copy();
            p.normalize();
			return ( distsq == 0
                    ? estimate_left_tangent2(d, len)
                    : p );
		}
	}
}

/**
 * Estimates the (backward) tangent at d[last].
 *
 * \note The tangent is "backwards", i.e. it is with respect to
 * decreasing index rather than increasing index.
 *
 * \pre 2 \<= len.
 * \pre d[len - 1] != d[len - 2].
 * \pre all[p in d] in_svg_plane(p).
 */
function estimate_right_tangent(d, len, tolerance_sq)

{
	//return estimate_right_tangent(d, len);
	//tolerance_sq =0 !
	//assert( 2 <= len );
	//assert( 0 <= tolerance_sq );
	var last = len - 1;
	var p, pi, t;
	var distsq;
    
	for (var i = last - 1;; i--) {
		pi = d[i].copy();
		t = new Point(pi.x - d[last].x, pi.y - d[last].y);
		distsq = dot(t, t);
		if ( tolerance_sq < distsq ) {
			p = t.copy();
            p.normalize();
			return p;
		}
		if (i == 0) {
			p = t.copy();
            p.normalize();
			return ( distsq == 0
                    ? estimate_right_tangent2(d, len)
                    : p );
		}
	}
}

/**
 *  Assign parameter values to digitized points using relative distances between points.
 *  use for data represented by one cubic Bezier
 */
function chord_length_parameterize(d, u, len )
{
	//int len = d->size();
	//assert( 2 <= len );
    
	/* First let u[i] equal the distance travelled along the path from d[0] to d[i]. */
	u[0] = 0.0;
	for (var i = 1; i < len; i++) {
        var p = d[i].minus(d[i-1]);
		u[i] = u[i-1] + p.norm();
	}
    
	/* Then scale to [0.0 .. 1.0]. */
	var tot_len = u[len - 1];
	//assert( tot_len != 0 );
	for (var i = 1; i < len; ++i) {
		u[i] /= tot_len;
	}
	u[len - 1] = 1;
    
//    for(var i = 0; i < u.length; ++i)
//    {
//        document.write("chord_length_parameterize "+u[i]+"<br>");
//    }
}

/*Functions used for computing the Bezier control points;*/

/**
 * Estimate the (forward) tangent at point d[first + 0.5].
 *
 * Unlike the center and right versions, this calculates the tangent in
 * the way one might expect, i.e., wrt increasing index into d.
 * \pre (2 \<= len) and (d[0] != d[1]).
 **/
function estimate_left_tangent2(d, len)
{
	//compute the left tangent using apportioned chords
	//the tangent ensures that the Bezier passes throught the second data point, also
	var p;
	var x0, x1, x2, x3;
	var tHat1;
	var L1, L2, L3;
	var t1, t2, m1, n1, m2, n2;
	//assert( len >= 2 );
	
	switch (len)
	{
        case 2:
            p = d[1].minus(d[0]);
            p.normalize();
            break;
        case 3:
            x0 = d[0].copy();
            x1 = d[1].copy();
            x2 = d[1].copy();
            x3 = d[2].copy();
            L1 = x1.minus(x0);
            L2 = x2.minus(x1);
            t1 = L1.norm()/(L1.norm() + L2.norm());
            tHat1 = (x1.minus(x0.time(Math.pow(1-t1,3))).minus(x3.time(Math.pow(t1, 3)))).time(1.0/(3*t1*pow(1-t1,2))*(1-t1));
            p = tHat1.minus(x0);
            p.normalize();
        default:
            x0 = d[0].copy();
            x1 = d[1].copy();
            x2 = d[2].copy();
            x3 = d[3].copy();
            L1 = x1.minus(x0);
            L2 = x2.minus(x1);
            L3 = x3.minus(x2);
            t1 = ( L1.norm() )/(L1.norm() + L2.norm() + L3.norm() );
            t2 = ( L1.norm() + L2.norm() )/(L1.norm() + L2.norm() + L3.norm() );
            
            m1 = B1(t1); m2 = B1(t2);
            n1 = B2(t1); n2 = B2(t2);
            L1 = x1.minus(x0.time(B0(t1))).minus(x3.time(B3(t1)));
            L2 = x2.minus(x0.time(B0(t2))).minus(x3.time(B3(t2)));
            tHat1 = L2.time(1.0/(m2-m1/n1)).minus(L1.time(1.0/(n1*m2-m1)));
            p = tHat1.minus(x0);
            p.normalize();
	}
    
	return p;
    
}

/**
 * Estimates the (backward) tangent at d[last - 0.5].
 *
 * \note The tangent is "backwards", i.e. it is with respect to
 * decreasing index rather than increasing index.
 *
 * \pre 2 \<= len.
 * \pre d[len - 1] != d[len - 2].
 * \pre all[p in d] in_svg_plane(p).
 */
function estimate_right_tangent2(d, len)
{
	//compute the right tangent using apportioned chords
	//the tangent ensures that the Bezier passes throught the third data point, also
	var p;
	var x0, x1, x2, x3;
	var tHat2;
	var L1, L2, L3;
	var t1, t2, m1, n1, m2, n2;
	var last = len - 1;
	//assert( len >= 2 );
	
	switch (len)
	{
        case 2:
            p = d[last-1].minus(d[last]);
            p.normalize();
            break;
        case 3:
            x0 = d[last-2].copy();
            x1 = d[last-1].copy();
            x2 = d[last-1].copy();
            x3 = d[last].copy();
            L1 = x1.minus(x0);
            L2 = x2.minus(x1);
            t1 = L1.norm()/(L1.norm() + L2.norm());
            tHat2 = (x1.minus(x0.time(Math.pow(1-t1,3))).minus(x3.time(Math.pow(t1, 3)))).time(1.0/(3*t1*pow(1-t1,2))*(1-t1));
            p = tHat2-x3;
            p.normalize();
        default:
            x0 = d[last-3].copy();
            x1 = d[last-2].copy();
            x2 = d[last-1].copy();
            x3 = d[last].copy();
            L1 = x1.minus(x0);
            L2 = x2.minus(x1);
            L3 = x3.minus(x2);
            t1 = ( L1.norm() )/(L1.norm() + L2.norm() + L3.norm() );
            t2 = ( L1.norm() + L2.norm() )/(L1.norm() + L2.norm() + L3.norm() );
            
            m1 = B1(t1); m2 = B1(t2);
            n1 = B2(t1); n2 = B2(t2);
            
            L1 = x1.minus(x0.time(B0(t1))).minus(x3.time(B3(t1)));
            L2 = x2.minus(x0.time(B0(t2))).minus(x3.time(B3(t2)));
            tHat2 = L1.time(m2).minus(L2.time(m1)).time(1/(m2*n1-m1*n2));
            
            p = tHat2.minus(x3);
            
            p.normalize();
	}
	return p;
    
    
	//assert( 2 <= len );
	//unsigned const last = len - 1;
	//unsigned const prev = last - 1;
	//assert( d->at(last) != d->at(prev) );
	//Point p = (d->at(prev) - d->at(last)); p.normalize();
	//return p;
}

function estimate_bi(bezier, ei, data, u, len)
{
	if(!(1 <= ei && ei <= 2))
		return;
	var oi = 3 - ei;
	var num = new Point(0.,0.);
	var den = 0.;
	for (var i = 0; i < len; ++i) {
        var ui = u[i];
		var b = [
			B0(ui),
			B1(ui),
			B2(ui),
			B3(ui)
		];
        
		{
			num.x += b[ei] * (b[0]  * bezier[0].x +
                              b[oi] * bezier[oi].x +
                              b[3]  * bezier[3].x +
                              - data[i].x);
            
			num.y += b[ei] * (b[0]  * bezier[0].y +
                              b[oi] * bezier[oi].y +
                              b[3]  * bezier[3].y +
                              - data[i].y);
		}
		den -= b[ei] * b[ei];
	}
    
	if (den != 0.) {
        
		bezier[ei] = num.time(1/den);
	}
	else {
		bezier[ei] = bezier[0].time(oi).plus(bezier[3].time(ei)).time(1.0/3.);
	}
}

/**
 *  Estimates length of tangent vectors P1, P2, when direction is given
 *  fills in bezier with the correct values for P1, P2
 *  Point bezier[4]
 */
function estimate_lengths(bezier, data, uPrime, len, tHat1, tHat2)
{
    //document.write("estimate_lengths len:"+len+"<br>");bezier
    //document.write("estimate_lengths bezier len:"+bezier.length+"<br>");
//    for(var i = 0; i < len; ++i)
//    {
//        document.write("estimate_lengths data("+data[i].x+ ","+data[i].y + ")<br> ");
//    }
//    for(var i = 0; i < bezier.length; ++i)
//    {
//        document.write("estimate_lengths bezier("+bezier[i].x+ ","+bezier[i].y + ")<br> ");
//    }
    
	var C = [];/* Matrix C. */
    for(var x = 0; x < 2; x++){
        C[x] = [];
        for(var y = 0; y < 2; y++){
            C[x][y] = 0.0;
        }
    }
    
	var X = [0.0, 0.0];      /* Matrix X. */
    
	/* First and last control points of the Bezier curve are positioned exactly at the first and
     last data points. */
	bezier[0] = new Point(data[0].x, data[0].y);
	bezier[3] = new Point(data[len - 1].x, data[len - 1].y);
    
	for (var i = 0; i < len; i++) {
		/* Bezier control point coefficients. */
		var b0 = B0(uPrime[i]);
		var b1 = B1(uPrime[i]);
		var b2 = B2(uPrime[i]);
		var b3 = B3(uPrime[i]);
        
		/* rhs for eqn */
		var a1 = {x:tHat1.x *b1, y:tHat1.y*b1};
		var a2 = {x:tHat2.x *b2, y:tHat2.y*b2};
        
		C[0][0] += dot(a1, a1);
		C[0][1] += dot(a1, a2);
		C[1][0] = C[0][1];
		C[1][1] += dot(a2, a2);
        
		/* Additional offset to the data point from the predicted point if we were to set bezier[1]
         to bezier[0] and bezier[2] to bezier[3]. */
		var shortfall = {x:(data[i].x-(bezier[0].x *(b0 + b1))-(bezier[3].x*(b2+b3))), y:(data[i].y-(bezier[0].y*(b0 + b1))-(bezier[3].y*(b2+b3)))};
		X[0] += dot(a1, shortfall);
		X[1] += dot(a2, shortfall);
	}
    
	/* We've constructed a pair of equations in the form of a matrix product C * alpha = X.
     Now solve for alpha. */
	var alpha_l, alpha_r;
    
	/* Compute the determinants of C and X. */
	var det_C0_C1 = C[0][0] * C[1][1] - C[1][0] * C[0][1];
	if ( det_C0_C1 != 0 ) {
		/* Apparently Kramer's rule. */
		var det_C0_X  = C[0][0] * X[1]    - C[0][1] * X[0];
		var det_X_C1  = X[0]    * C[1][1] - X[1]    * C[0][1];
		alpha_l = det_X_C1 / det_C0_C1;
		alpha_r = det_C0_X / det_C0_C1;
	} else {
		/* The matrix is under-determined.  Try requiring alpha_l == alpha_r.
         *
         * One way of implementing the constraint alpha_l == alpha_r is to treat them as the same
         * variable in the equations.  We can do this by adding the columns of C to form a single
         * column, to be multiplied by alpha to give the column vector X.
         *
         * We try each row in turn.
         */
		var c0 = C[0][0] + C[0][1];
		if (c0 != 0) {
			alpha_l = alpha_r = X[0] / c0;
		} else {
			var c1 = C[1][0] + C[1][1];
			if (c1 != 0) {
				alpha_l = alpha_r = X[1] / c1;
			} else {
				/* Let the below code handle this. */
				alpha_l = alpha_r = 0.;
			}
		}
	}
    
	/* If alpha negative, use the Wu/Barsky heuristic (see text).  (If alpha is 0, you get
     coincident control points that lead to divide by zero in any subsequent
     NewtonRaphsonRootFind() call.) */
	/// \todo Check whether this special-casing is necessary now that
	/// NewtonRaphsonRootFind handles non-positive denominator.
	if ( alpha_l < 1.0e-6 ||
		alpha_r < 1.0e-6   )
	{
		var p = {x:data[len - 1].x- data[0].x, y:data[len - 1].y- data[0].y};
		alpha_l = alpha_r = ( Math.sqrt(p.x*p.x +p.y*p.y)/3.0 );
	}
    
	/* Control points 1 and 2 are positioned an alpha distance out on the tangent vectors, left and
     right, respectively. */
	//alpha min is 4 ..so the points are far enough away
	/*bezier[1] = tHat1 * max(alpha_l,4.) + bezier[0];
     bezier[2] = tHat2 * max(alpha_r,4.) + bezier[3];*/
    
	bezier[1] = new Point(tHat1.x * alpha_l + bezier[0].x,tHat1.y * alpha_l + bezier[0].y);
	bezier[2] = new Point(tHat2.x * alpha_r + bezier[3].x,tHat2.y * alpha_r + bezier[3].y);
    
	return;
}

/* integer arithmetic */
function mod(a, n) {
	return a>=n ? a%n : a>=0 ? a : n-1-(-1-a)%n;
}
/*for two pixels that are not neighbours: P0, P1
 determine a pixel neighbouring P0 that is on the line uniting P0 to P1
 returns direction with center in P0: (-1,0,1)
 */
function ddir(a)
{
	var p = new Point();
	p.x = 0; p.y = 0;
	if(a.y==a.x==0)
		return p;
	//atan2 gives an error if both arguments are 0;
	var theta = Math.atan2(a.y,a.x)* 180 / Math.PI;
	if(theta<0)
	{
		if((a.x>=0)&&(a.y>=0))
			theta +=180; //
		else if((a.x<=0)&&(a.y<=0))
			theta +=360; //
		else
			document.write("oops!");
	}
	else
	{
		if((a.x>=0)&&(a.y<=0)&&theta<270)
			theta = 360-theta; //
	}
    //document.write("angle:" + theta + "<br>");
	//0/360
	if(theta<22.5 || theta>=337.5)
	{
		p.x = 1; p.y = 0;
		return p;
	}
	//45
	if(theta>=22.5 && theta<67.5)
	{
		p.x = 1; p.y = 1;
		return p;
	}
	//90
	if(theta>=67.5 && theta<112.5)
	{
		p.x = 0; p.y = 1;
		return p;
	}
	//135
	if(theta>=112.5 && theta<157.5)
	{
		p.x = -1; p.y = 1;
		return p;
	}
	//180
	if(theta>=157.5 && theta<202.5)
	{
		p.x = -1; p.y = 0;
		return p;
	}
	//225
	if(theta>=202.5 && theta<247.5)
	{
		p.x = -1; p.y = -1;
		return p;
	}
	//270
	if(theta>=247.5 && theta<292.5)
	{
		p.x = 0; p.y = -1;
		return p;
	}
	//315
	if(theta>=292.5 && theta<337.5)
	{
		p.x = 1; p.y = -1;
		return p;
	}
	return p;
}

/* return 1 if a <= b < c < a, in a cyclic sense (mod n) */
function cyclic(a, b, c) {
    //document.write("a, b, c:" + a + "," + b + "," + c + "<br>");
	if (a<=c) {
		return (a<=b && b<c);
	} else {
		return (a<=b || b<c);
	}
}

function xprod(p1, p2)
{
	return p1.x*p2.y - p1.y*p2.x;
}

function sign(x)
{
	return ((x)>0 ? 1 : (x)<0 ? -1 : 0);
}

/* return (p1-p0)x(p2-p0), the area of the parallelogram */
function dpara(p0, p1, p2) {
	var x1, y1, x2, y2;
    
	x1 = p1.x-p0.x;
	y1 = p1.y-p0.y;
	x2 = p2.x-p0.x;
	y2 = p2.y-p0.y;
    
	return x1*y2 - x2*y1;
}

/* ddenom/dpara have the property that the square of radius 1 centered
 at p1 intersects the line p0p2 iff |dpara(p0,p1,p2)| <= ddenom(p0,p2) */
function ddenom(p0, p2) {
	var r = new Point();
	r.y = sign(p2.x-p0.x);
	r.x = -sign(p2.y-p0.y);
	return r.y*(p2.x-p0.x) - r.x*(p2.y-p0.y);
}

function  dot(P1, P2) {return P1.x*P2.x + P1.y*P2.y;}

function floordiv(a, n) {
	return (a>=0 ? parseInt(a/n) : -1-parseInt((-1-a)/n));
}

/**
 * Evaluate a Bezier curve at parameter value \a t.
 *
 * \param degree The degree of the Bezier curve: 3 for cubic, 2 for quadratic etc.
 * \param V The control points for the Bezier curve.  Must have (\a degree+1)
 *    elements.
 * \param t The "parameter" value, specifying whereabouts along the curve to
 *    evaluate.  Typically in the range [0.0, 1.0].
 *
 * Let s = 1 - t.
 * BezierII(1, V) gives (s, t) * V, i.e. t of the way
 * from V[0] to V[1].
 * BezierII(2, V) gives (s**2, 2*s*t, t**2) * V.
 * BezierII(3, V) gives (s**3, 3 s**2 t, 3s t**2, t**3) * V.
 *
 * The derivative of BezierII(i, V) with respect to t
 * is i * BezierII(i-1, V'), where for all j, V'[j] =
 * V[j + 1] - V[j].
 */
function bezier_pt(degree, V, t)
{
	/** Pascal's triangle. */
    var pascal = new Array();
    pascal[0] = [1];
    pascal[1] = [1,1];
    pascal[2] = [1,2,1];
    pascal[3] = [1,3,3,1];
    
	//assert( degree < 4 );
	var s = 1.0 - t;
    
	/* Calculate powers of t and s. */
	var spow = new Array();
	var tpow = new Array();
	spow[0] = 1.0; spow[1] = s;
	tpow[0] = 1.0; tpow[1] = t;
	for (var i = 1; i < degree; ++i) {
		spow[i + 1] = spow[i] * s;
		tpow[i + 1] = tpow[i] * t;
	}
    
	var ret = new Point(V[0].x * spow[degree], V[0].y * spow[degree]);
	for (var i = 1; i <= degree; ++i) {
        var p = new Point(V[i].x, V[i].y);
		ret = ret.plus(p.time(pascal[degree][i] * spow[degree - i] * tpow[i]));
	}
    //document.write("ret:"+ret.x + "," + ret.y+"<br>");
	return ret;
}

/**
 * Find the maximum squared distance of digitized points to fitted curve, and (if this maximum
 * error is non-zero) set \a *splitPoint to the corresponding index.
 *
 * \pre 2 \<= len.
 * \pre u[0] == 0.
 * \pre u[len - 1] == 1.0.
 * \post ((ret == 0.0)
 *        || ((*splitPoint \< len - 1)
 *            \&\& (*splitPoint != 0 || ret \< 0.0))).
 */
function compute_max_error_ratio(d, u, len, bezCurve, tolerance, splitPoint)

{
//    for(var i = 0; i < d.length; ++i)
//    {
//        document.write("("+d[i].x+ ","+d[i].y + ")<br>");
//    }
//    
//    for(var i = 0; i < u.length; ++i)
//    {
//        document.write(u[i]+"<br>");
//    }
//    
//    for(var i = 0; i < 4; ++i)
//    {
//        document.write("("+bezCurve[i].x+ ","+bezCurve[i].y + ")<br>");
//    }
//    document.write("<br>");
    
	//assert( 2 <= len );
	var last = len - 1;
	//assert( u[0] == 0.0 );
	//assert( u[last] == 1.0 );
	/* I.e. assert that the error for the first & last points is zero.
     * Otherwise we should include those points in the below loop.
     * The assertion is also necessary to ensure 0 < splitPoint < last.
     */
    
	var maxDistsq = 0.0; /* Maximum error */
	var max_hook_ratio = 0.0;
	var snap_end = 0;
	var prev = new Point(bezCurve[0].x, bezCurve[0].y);
	for (var i = 1; i <= last; i++) {
		var curr = bezier_pt(3, bezCurve, u[i]);
		var distsq = dot( curr.minus(d[i]),  curr.minus(d[i]) );
        //document.write("distsq:"+distsq+"<br>");
		if ( distsq > maxDistsq ) {
			maxDistsq = distsq;
			//*splitPoint = i;
            splitPoint[0] = i;
		}
		var hook_ratio = compute_hook(prev, curr, .5 * (u[i - 1] + u[i]), bezCurve, tolerance);
		
		if (max_hook_ratio < hook_ratio) {
			max_hook_ratio = hook_ratio;
			snap_end = i;
		}
		prev = curr.copy();
	}
	var dist_ratio = Math.sqrt(maxDistsq) / tolerance;
	var ret;
	if (max_hook_ratio <= dist_ratio) {
		ret = dist_ratio;
	} else {
		//assert(0 < snap_end);
		ret = -max_hook_ratio;
		splitPoint[0] = snap_end - 1;
	}
//    document.write("splitPoint:"+splitPoint[0]+"<br>");
//	//assert( ret == 0.0 || ( ( splitPoint < last ) && ( splitPoint != 0 || ret < 0. ) ) );
//    document.write("maxDistsq:"+maxDistsq+"<br>");
//    document.write("ret:"+ret+"<br>");
	return ret;
}


/**
 * Whereas compute_max_error_ratio() checks for itself that each data point
 * is near some point on the curve, this function checks that each point on
 * the curve is near some data point (or near some point on the polyline
 * defined by the data points, or something like that: we allow for a
 * "reasonable curviness" from such a polyline).  "Reasonable curviness"
 * means we draw a circle centred at the midpoint of a..b, of radius
 * proportional to the length |a - b|, and require that each point on the
 * segment of bezCurve between the parameters of a and b be within that circle.
 * If any point P on the bezCurve segment is outside of that allowable
 * region (circle), then we return some metric that increases with the
 * distance from P to the circle.
 *
 *  Given that this is a fairly arbitrary criterion for finding appropriate
 *  places for sharp corners, we test only one point on bezCurve, namely
 *  the point on bezCurve with parameter halfway between our estimated
 *  parameters for a and b.  (Alternatives are taking the farthest of a
 *  few parameters between those of a and b, or even using a variant of
 *  NewtonRaphsonFindRoot() for finding the maximum rather than minimum
 *  distance.)
 */
function compute_hook(a, b, u, bezCurve, tolerance)
{
	var P = bezier_pt(3, bezCurve, u);
	var diff = a.plus(b).time(0.5).minus(P);
	var dist = diff.norm();
	if (dist < tolerance) {
		return 0;
	}
	P = b.minus(a);
	var allowed = P.norm() + tolerance;
	return dist / allowed;
	/** \todo
     * effic: Hooks are very rare.  We could start by comparing
     * distsq, only resorting to the more expensive L2 in cases of
     * uncertainty.
     */
}


