/* The arithmetic-geometric mean of two non-negative numbers */
Math.agm = function(a0,g0)
{
    var maxIter = 50;
    var an = (a0+g0)/2;
    var gn = Math.sqrt(a0*g0);

    for (var iter = 0; (iter < maxIter) && (Math.abs(an-gn) > 1e-15); iter++)
    {
        a0 = 0.5 * (an + gn);
        g0 = Math.sqrt(an*gn);
        an = a0;
        gn = g0;
    };
    if (iter == maxIter)
        console.warn("Math.agm hit iteration limit of " + maxIter + ", probably didn't converge.");
    return an;
};

/* EllipticK(m) - The complete elliptic integral of the first type.
 * The argument is the *parameter* m = k^2, where k is the *modulus*.
 * The parameter must satisfy m < 1.
 * In terms of the integral definition, we have
 * K(m) = \int_0^{\pi/2} \frac{d\theta}{\sqrt{1 - m (\sin\theta)^2}}
 * See http://dlmf.nist.gov/19.8.E5 for this method.
 */
Math.EllipticK = function( m )
{
    var kprime = Math.sqrt(1 - m);
    return 0.5 * Math.PI / Math.agm(1, kprime);
};

Math.Ellipj = function(u, m) {
    // Constants
    const MACHEP = 1.1102230246251565e-16;
    const PIO2 = Math.PI / 2.0;

    // Variables
    let ai, b, phi, t, twon;
    const a = [1.0];
    const c = [Math.sqrt(m)];
    let i = 0;

    // Check for special cases
    if (m < 0.0 || m > 1.0) {
        console.error("ellpj: DOMAIN error");
        return { sn: 0.0, cn: 0.0, dn: 0.0, ph: 0.0 };
    }

    if (m < 1.0e-9) {
        t = Math.sin(u);
        b = Math.cos(u);
        ai = 0.25 * m * (u - t * b);
        return {
            sn: t - ai * b,
            cn: b + ai * t,
            ph: u - ai,
            dn: 1.0 - 0.5 * m * t * t
        };
    }

    if (m >= 0.9999999999) {
        ai = 0.25 * (1.0 - m);
        b = Math.cosh(u);
        t = Math.tanh(u);
        phi = 1.0 / b;
        twon = b * Math.sinh(u);
        ai *= t * phi;
        return {
            sn: t + ai * (twon - u) / (b * b),
            ph: 2.0 * Math.atan(Math.exp(u)) - PIO2 + ai * (twon - u) / b,
            cn: phi - ai * (twon - u),
            dn: phi + ai * (twon + u)
        };
    }

    // A. G. M. scale
    b = Math.sqrt(1.0 - m);
    c[0] = b;
    twon = 1.0;

    while (Math.abs(c[i] / a[i]) > MACHEP) {
        if (i > 7) {
            console.error("ellpj: OVERFLOW error");
            return { sn: 0.0, cn: 0.0, dn: 0.0, ph: 0.0 };
        }
        ai = a[i];
        ++i;
        c[i] = (ai - b) / 2.0;
        t = Math.sqrt(ai * b);
        a[i] = (ai + b) / 2.0;
        b = t;
        twon *= 2.0;
    }

    // Backward recurrence
    phi = twon * a[i] * u;

    do {
        t = c[i] * Math.sin(phi) / a[i];
        b = phi;
        phi = (Math.asin(t) + phi) / 2.0;
    } while (--i);

    t = Math.sin(phi);

    return {
        sn: t,
        cn: Math.cos(phi),
        dn: Math.sqrt(1.0 - m * t * t),
        ph: phi
    };
}

Math.JacobiSn = function(kOrU, time, timeperiod) {
    let a_elliptic_integral, sn, cn, angle_amplitude;
    // Based on the article "https://doi.org/10.1016/j.jfluidstructs.2019.01.020"
    if (0 <= kOrU && kOrU < 1) {
        a_elliptic_integral = Math.EllipticK(Math.sqrt(kOrU));
        sn = Math.Ellipj(4 * a_elliptic_integral * time / timeperiod, Math.sqrt(kOrU)).sn;
        angle_amplitude = sn;
    } else if (-1 < kOrU && kOrU < 0) {
        a_elliptic_integral = Math.EllipticK(Math.sqrt(Math.abs(kOrU)));
        cn = Math.Ellipj(4 * a_elliptic_integral * time / timeperiod - a_elliptic_integral, Math.sqrt(Math.abs(kOrU))).cn;
        angle_amplitude = cn;
    } else {
        angle_amplitude = "Error";
    }
    return angle_amplitude;
};
