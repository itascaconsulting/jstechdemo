



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


r = math.range(0,p0,p0/25,true);
console.log(r.map(u)._data);

var u_r = function (p_i,r) {
    if (p_i > p_cr) {
        return {R:r0, u:r0**2*(1+nu)*(p0-p_i)/E/r}; // Kirsch solution
    }
    var delta = Math.PI/4.0 + phi/2.0;
    var Q = Math.tan(delta)/(Math.tan(delta-phi))-1;
    var tan2phi = Math.tan(Math.PI/4.0 + phi/2.0)**2
    var cotphi = math.cot(phi);
    var denom = (1+tan2phi)*(p_i + c*cotphi);
    var R = r0 * ((2*p0-sigma_cm+(1+tan2phi)*c*cotphi)/denom)**(1/Q);
    var b = ((tan2phi-1)*p0+sigma_cm)/(tan2phi+1)*R**2;
    var t0 = p0+c*cotphi-(p_i + c*cotphi*(R/r0)**Q);
    var t = (1-nu)/E*R**2*t0 + b*(1+nu)/E;
    var u_r = (1-nu)/E*(p_i*r**(Q+1)/r0**Q - p0*r) + t/r;
    return {R:R, u: u_r};
}

// plastic radius only matches when p_i === 0

show(u_r(0,r0));
show(u_r(1e6,r0));
show(u_r(6e6,r0));
show(u_r(7e6,r0));
// show(u_r(0,r0+0.1));
// show(u_r(0,r0+0.2));
// show(u_r(0,r0+0.3));
// show(u_r(0,7));

console.log(r.map(function (d) { return u_r(d,r0).u; })._data);
