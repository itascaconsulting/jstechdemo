onload = function () {


// default inputs
A = 201e-6;
E = 1.4e11;
slen = 0.1;
Kbond = 8;
Sbond = 2.8e5;
n = 5;
Yield = 180e3;
Rupture = 0.2;
ds = 1e7;
dy = 62.8e3;
dr = 0.41;
theta = 40.0;
disp = 4e-3;


var controls = [["Cable area [m^2]",           100e-6,  A,       100e-5],
                ["Cable Young's modulus [Pa]",   1e10, E,       1e12],
                ["Node spacing [m]",           0.05,    slen,    0.4],
                ["Grout stiffness [N/m/m]",          7,     Kbond,   11],
                ["Grout strength [N/m]",          1e4,     Sbond,   1e6],
                ["Number of nodes [half]",              5,       n,       25],
                ["Cable yield strength [N]",          50e3,    Yield,   250e3],
                ["Cable rupture strain []",        0.05,    Rupture, 0.8],
                ["Cable shear stiffness [N/m]",    1e6,     ds,      4e7],
                ["Cable shear strength [N]",    1e4,     dy,      2e5],
                ["Shear rupture strain []",  1e-1,    dr,      6e-1],
                ["Joint friction angle [deg]", 0.0, theta, 50.0],
                ["Displacement [m]",           1e-5,    disp,    5e-2]];

sliders = [];

for (i=0; i<controls.length; i++) {
    name = controls[i][0];
    min_val = controls[i][1];
    start_val = controls[i][2];
    max_val = controls[i][3];

    var show_sliders = document.getElementById("div__sliders");

    // show_sliders.appendChild(h1);

    var h1 = document.createElement("div");
    h1.setAttribute('style', 'white-space: pre;');
    h1.textContent += name + "\r\n";
    show_sliders.appendChild(h1);


    var temp_slider = document.createElement("div");
    temp_slider.setAttribute("id", "slider"+i);
    show_sliders.appendChild(temp_slider);
    var tmp_d3slider = d3.slider().min(min_val).max(max_val).ticks(10).showRange(true).value(start_val).callback(scb);

    d3.select('#slider'+i).call(tmp_d3slider);

    sliders.push(tmp_d3slider)
}


function displacements(A,E,slen,Kbond,Sbond,n,Yield,Rupture,ds,dy,dr,theta,disp) {
    //derived properties
    n = 2*parseInt(n);
    Kbond = Math.pow(10,Kbond);
    c1 = E*A/slen;
    c2 = Kbond*slen;
    c3 = Sbond*slen;
    c2i = math.multiply(math.ones(n), c2);
    ab = math.sparse();
    b = math.matrix();
    b.resize([n]);
    db0 = 0.0;
    db1 = disp;

    // construct matrix
    for(i = 0; i < n; i++) {
        if (i > Math.floor(n/2-1)) {
            b.set([i], -db1*c2i.get([i]));
        } else {
            b.set([i], -db0*c2i.get([i]));
        }
        if (i==0) {
            ab.set([i,i+1], c1);
            ab.set([i, i], -c1-c2i.get([i]));
        } else if (i==n-1) {
            ab.set([i,i-1], c1);
            ab.set([i, i], -c1-c2i.get([i]));
        } else {
            ab.set([i, i-1], c1);
            ab.set([i, i+1], c1);
            ab.set([i, i], -2*c1-c2i.get([i]));
        }
    }
    ndis = math.lusolve(ab,b);
    x = math.multiply(math.range(0, n), slen).toArray();
    ddata = math.flatten(ndis).toArray();
    return {x:x, ddata:ddata};
}

solution = displacements(A,E,slen,Kbond,Sbond,n,Yield,Rupture,ds,dy,dr,theta,disp)

// plotting

var data = [ { label: "Displacement",
               x: x.map(function(i){return i*2}),
               y: ddata },
             ] ;

var xy_chart = d3_xy_chart()
    // .width(960)
    // .height(500)
    .xlabel("Length along cable [m]")
    .ylabel("Displacement [m]") ;

var svg = d3.select("#section__graph").append("svg")
    .datum(data)
    .call(xy_chart) ;

function d3_xy_chart() {
    var width = 640,
        height = 480,
        xlabel = "X Axis Label",
        ylabel = "Y Axis Label" ;

    function chart(selection) {
        selection.each(function(datasets) {
            //
            // Create the plot.
            //
            var margin = {top: 20, right: 80, bottom: 30, left: 50},
                innerwidth = width - margin.left - margin.right,
                innerheight = height - margin.top - margin.bottom ;

            var x_scale = d3.scale.linear()
                .range([0, innerwidth])
                .domain([ d3.min(datasets, function(d) { return d3.min(d.x)-0.05*d3.max(d.x); }),
                          d3.max(datasets, function(d) { return 1.05*d3.max(d.x); }) ]) ;

            var y_scale = d3.scale.linear()
                .range([innerheight, 0])
                .domain([ d3.min(datasets, function(d) { return d3.min(d.y)-.05*d3.max(d.y); }),
                          d3.max(datasets, function(d) { return 1.05*d3.max(d.y); }) ]) ;

            var color_scale = d3.scale.category10()
                .domain(d3.range(datasets.length)) ;

            var x_axis = d3.svg.axis()
                .scale(x_scale)
                .orient("bottom") ;

            var y_axis = d3.svg.axis()
                .scale(y_scale)
                .orient("left") ;

            var x_grid = d3.svg.axis()
                .scale(x_scale)
                .orient("bottom")
                .tickSize(-innerheight)
                .tickFormat("") ;

            var y_grid = d3.svg.axis()
                .scale(y_scale)
                .orient("left")
                .tickSize(-innerwidth)
                .tickFormat("") ;

            var draw_line = d3.svg.line()
                .interpolate("basis")
                .x(function(d) { return x_scale(d[0]); })
                .y(function(d) { return y_scale(d[1]); }) ;

            var svg = d3.select(this)
                .attr("width", width)
                .attr("height", height)
                .append("g")
                .attr("transform", "translate(" + margin.left + "," + margin.top + ")") ;

            svg.append("g")
                .attr("class", "x grid")
                .attr("transform", "translate(0," + innerheight + ")")
                .call(x_grid) ;

            svg.append("g")
                .attr("class", "y grid")
                .call(y_grid) ;

            svg.append("g")
                .attr("class", "x axis")
                .attr("transform", "translate(0," + innerheight + ")")
                .call(x_axis)
                .append("text")
                .attr("dy", "-.71em")
                .attr("x", innerwidth)
                .style("text-anchor", "end")
                .text(xlabel) ;

            svg.append("g")
                .attr("class", "y axis")
                .call(y_axis)
                .append("text")
                .attr("transform", "rotate(-90)")
                .attr("y", 6)
                .attr("dy", "0.71em")
                .style("text-anchor", "end")
                .text(ylabel) ;

            var data_lines = svg.selectAll(".d3_xy_chart_line")
                .data(datasets.map(function(d) {return d3.zip(d.x, d.y);}))
                .enter().append("g")
                .attr("class", "d3_xy_chart_line") ;

            data_lines.append("path")
                .attr("class", "line")
                .attr("d", function(d) {return draw_line(d); })
                .attr("stroke", function(_, i) {return color_scale(i);}) ;

            data_lines.append("text")
                .datum(function(d, i) { return {name: datasets[i].label, final: d[d.length-1]}; })
                .attr("transform", function(d) {
                    return ( "translate(" + x_scale(d.final[0]) + "," +
                             y_scale(d.final[1]) + ")" ) ; })
                .attr("x", 3)
                .attr("dy", ".35em")
                .attr("fill", function(_, i) { return color_scale(i); })
                .text(function(d) { return d.name; }) ;

        }) ;
    }

    chart.width = function(value) {
        if (!arguments.length) return width;
        width = value;
        return chart;
    };

    chart.height = function(value) {
        if (!arguments.length) return height;
        height = value;
        return chart;
    };

    chart.xlabel = function(value) {
        if(!arguments.length) return xlabel ;
        xlabel = value ;
        return chart ;
    } ;

    chart.ylabel = function(value) {
        if(!arguments.length) return ylabel ;
        ylabel = value ;
        return chart ;
    } ;

    return chart;
}

var update_plots = function() {
  // console.log("wow");
  // ok, this happens every time the slider is clicked or clicked and dragged.

  data = [ { label: "Displacement",
                 x: x,
                 y: ddata },] ;

  xy_chart = d3_xy_chart()
      // .width(960)
      // .height(500)
      .xlabel("Length along cable [m]")
      .ylabel("Displacement [m]") ;

  document.querySelectorAll("#section__graph")[0].innerHTML ="<h1>Output Graph</h1>";
  // simply clears the previous graph and re-adds the h1 title.

  d3.select("#section__graph").append("svg")
      .datum(data)
      .call(xy_chart) ;
};



function scb() {
    local_values = [];
    data = "";
    for (i=0; i<sliders.length; i++) {
        data += sliders[i].value() + ", ";
        local_values.push(sliders[i].value());
    }
    displacements.apply(null, local_values);
    update_plots();
}

scb()





};
