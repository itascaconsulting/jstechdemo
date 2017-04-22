



// input
var nu = 0.3,
    E = 1.8e9,
    c = 1.5e6,
    phi = 23.0,
    p0 = 7e6,
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

var u_re = function (p_i,r) {
    return {R:r0, u:r0**2*(1+nu)*(p0-p_i)/E/r}; // Kirsch solution
}


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
    //var sigma_R = p0-b/R**2; // which is the critical stress..
    //var sigma_Rp = (p_i+c*cotphi)*(R/r0)**Q-c*cotphi; // which is the critical stress..
    var sigma_re = p0-b/r**2;
    var sigma_rp = (p_i+c*cotphi)*(r/r0)**Q-c*cotphi;



    if (r>R) {
        var u_R = (1-nu)/E*(p_i*R**(Q+1)/r0**Q - p0*R) + t/R; // disp at R
        var u_re = R**2*(1+nu)*(p0-p_cr)/E/R;
        var dadj = R**2*(1+nu)*(p0-p_cr)/E/u_R;

        var du = u_R-u_re;
        // stress match to Kirsch solution? not sure if this is correct
        return {R:R,
                u:R**2*(1+nu)*(p0-p_cr)/E/(dadj+(r-R)),
                sigma: sigma_re};
    }

    return {R:R, u: u_r, sigma: sigma_rp};
}



show(u_r(0,r0));
show(u_r(1e6,r0));
show(u_r(6e6,r0));
show(u_r(7e6,r0));
show("disp at R: " + u_r(0,u_r(0,r0).R).u);


var np_linspace = function (start,end,num,endpoint) {
    var local_endpoint = endpoint || false;
    var delta = end-start;
    return math.range(start,end,delta/num,local_endpoint)._data;
}


support_pressures = np_linspace(0,p0,25,true);
reaction = support_pressures.map(u);
reaction2 = support_pressures.map(function (d) { return u_r(d,r0).u; });


var depth = np_linspace(r0,r0+50.0,100,true);
u_depth = depth.map(function (d) { return u_r(0,d).u;});
u_depthe = depth.map(function (d) { return u_re(0,d).u;});
u_stress = depth.map(function (d) { return u_r(0,d).sigma;});


var plot_xy = function (datasets) {
    var margin = {top: 30, right: 20, bottom: 30, left: 70},
        width = 400 - margin.left - margin.right,
        height = 220 - margin.top - margin.bottom;
    var     x = d3.scale.linear().range([0, width]);
    var     y = d3.scale.linear().range([height, 0]);

    var xAxis = d3.svg.axis().scale(x)
        .orient("bottom")
        .ticks(5);
    var yAxis = d3.svg.axis().scale(y)
        .orient("left")
        .ticks(5)
        .tickFormat(d3.format("2e"));

    var valueline = function(xa, ya){
        return d3.svg.line()
            .x(function(d,i) { return x(xa[i]); })
            .y(function(d,i) { return y(ya[i]); })
        (Array(xa.length));
    }

    var chart1 = d3.select("body")
        .append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom)
        .append("g")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    var xmin = math.max(datasets[0][0]);
    var xmax = math.min(datasets[0][0]);
    var ymin = math.max(datasets[0][1]);
    var ymax = math.min(datasets[0][1]);
    datasets.forEach( function (d) {
        xarray = d[0];
        yarray = d[1];
        if (math.max(xarray) > xmax) { xmax = math.max(xarray); }
        if (math.min(xarray) < xmin) { xmin = math.min(xarray); }
        if (math.max(yarray) > ymax) { ymax = math.max(yarray); }
        if (math.min(yarray) < ymin) { ymin = math.min(yarray); }
    });

    x.domain([xmin,xmax]);
    y.domain([0,ymax]);

    var colors = d3.scale.category10().domain([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]);
    datasets.forEach(function (d, i) {
        xarray = d[0];
        yarray = d[1];

        // Add the valueline path.
        chart1.append("path")
            .attr("class", "line")
            .attr("stroke", colors(i%10))
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
    });
}

plot_xy([[depth, u_stress]])
plot_xy([[depth, u_depth], [depth, u_depthe]]);
plot_xy([[reaction, support_pressures], [reaction2, support_pressures]]);
