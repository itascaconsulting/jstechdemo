



// input
nu = 0.3;
E = 1.8e9;
c = 1.5e6;
phi = 23.0;
p0 = 7e6;
r0 = 5.0;

// derived
phi *= Math.PI/180.0;
k = (1+Math.sin(phi))/(1-Math.sin(phi));
sigma_cm = 2*c*Math.cos(phi)/(1-Math.sin(phi)); // ucs
p_cr = (2 * p0 - sigma_cm)/(1+k); // critical pressure

var r_p = function (p_i) {
    if (p_i > p_cr) {
        return r0;
    }
    return r0*(2*(p0*(k-1)+sigma_cm)/((1+k)*((k-1)*p_i+sigma_cm)))**(1/(k-1));
}
var u = function (p_i) {  // tunnel wall radial displacement, inward is positive.
    if (p_i > p_cr) {
        return r0*(1+nu)*(p0-p_i)/E; // elastic, this is the Kirsch solution
    }
    return r0*(1+nu)*(2*(1-nu)*(p0-p_cr)*(r_p(p_i)/r0)**2-(1-2*nu)*(p0-p_i))/E;
}


function show(data) {
    console.log(math.format(data,  {notation: 'engineering'}));
}

// show(sigma_cm);
// show(k);
// show(p_cr);
show(r_p(0.0));
show(r_p(1e6));
show(r_p(6e6));
show(r_p(7e6));
show(u(0.0));
show(u(1e6));
show(u(6e6));
show(u(7e6));


var u_r = function (p_i,r) {
    if (p_i > p_cr) {
        return {R:r0, u:r0**2*(1+nu)*(p0-p_i)/E/r}; // Kirsch solution
    }
    var delta = Math.PI/4.0 + phi/4.0; // fracture angle, a free parameter
    var Q = Math.tan(delta)/(Math.tan(delta-phi))-1;
    var tan2phi = Math.tan(Math.PI/4.0 + phi/2.0)**2
    var cotphi = math.cot(phi);
    var denom = (1+tan2phi)*(p_i + c*cotphi);
    var R = r0 * ((2*p0-sigma_cm+(1+tan2phi)*c*cotphi)/denom)**(1/Q);
    var b = ((tan2phi-1)*p0+sigma_cm)/(tan2phi+1)*R**2;
    var t0 = p0+c*cotphi-(p_i + c*cotphi*(R/r0)**Q);
    var t = (1-nu)/E*R**2*t0 + b*(1+nu)/E;
    var u_r = (1-nu)/E*(p_i*r**(Q+1)/r0**Q - p0*r) + t/r;
    if (r>=R) {
        return {R:R, u: r0**2*(1+nu)*(p0-p_i)/E/r};
    }
    return {R:R, u: u_r};
}



show(u_r(0,r0));
show(u_r(1e6,r0));
show(u_r(6e6,r0));
show(u_r(7e6,r0));
// show(u_r(0,r0+0.1));
// show(u_r(0,r0+0.2));
// show(u_r(0,r0+0.3));
// show(u_r(0,7));

var np_linspace = function (start,end,num,endpoint) {
    var local_endpoint = endpoint || false;
    var delta = end-start;
    return math.range(start,end,delta/num,local_endpoint)._data;
}


support_pressures = np_linspace(0,p0,25,true);
reaction = support_pressures.map(u);
console.log(reaction);
reaction2 = support_pressures.map(function (d) { return u_r(d,r0).u; });
console.log(reaction2);


var depth = np_linspace(r0,r0+4.0,25,true);
console.log(depth);
u_depth = depth.map(function (d) { return u_r(0,d).u;});
console.log(u_depth);

var plot_xy = function (xarray, yarray) {
    console.log(xarray);
    console.log(yarray);

    var margin = {top: 30, right: 20, bottom: 30, left: 50},
        width = 400 - margin.left - margin.right,
        height = 220 - margin.top - margin.bottom;
    var     x = d3.scale.linear().range([0, width]);
    var     y = d3.scale.linear().range([height, 0]);
    //
    //d3.extent(yarray)
    var xAxis = d3.svg.axis().scale(x).orient("bottom").ticks(5);
    var yAxis = d3.svg.axis().scale(y).orient("left").ticks(5);

    var data = {x: xarray, y: yarray};

    var valueline = function(xa, ya){
        return d3.svg.line()
            .x(function(d,i) { return x(xa[i]); })
            .y(function(d,i) { return y(ya[i]); })
        (Array(xa.length));
    }

    // var valueline = d3.svg.line()
    //     .x(function(d,i) { return x(d.x[i]);})
    //     .y(function(d,i) { return y(d.y[i]);})

    var chart1 = d3.select("body")
        .append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom)
        .append("g")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    x.domain(d3.extent(d3.extent(xarray)));
    y.domain([0,d3.max(d3.extent(yarray))]);

    // Add the valueline path.
    chart1.append("path")
        .attr("class", "line")
        .attr("d", valueline(xarray, yarray));

    // Add the X Axis
    chart1.append("g")
        .attr("class", "x axis")
        .attr("transform", "translate(0," + height + ")")
        .call(xAxis);

    // Add the Y Axis
    chart1.append("g")
        .attr("class", "y axis")
        .call(yAxis);
}

plot_xy(depth, u_depth);
