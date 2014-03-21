#include "hartreefock.h"

#include <QTimer>
#include <QDebug>
#include <iostream>
#include <math.h>

// qreal rk(){
//     double *k1_x, *k2_x, *k2_v, *k3_v, k3_x, *k4_v, k4_x, *k1_v, *v, *x;
//     qreal dt;
//     int n_cicli;
//
//     for(unsigned int i = 0; i < n_cicli; ++i){
//         k1_x = v[i]*dt;
//         k1_v = - x[i]*dt;
//         /// v[i+1] = v[i] + k1_v*0.5;
//         /// x[i+1] = x[i] + k1_x*0.5;
//         k2_x = (v[i] + k1_v*0.5)*dt;
//         k2_v = -(x[i] + k1_x*0.5)*dt;
//         ///v[i+1] = v[i] + k2_v*0.5;
//         ///x[i+1] = x[i] + k2_x*0.5;
//         k3_x =(v[i]+k2_v*0.5)*dt;
//         k3_v= -(x[i]+k2_x*0.5)*dt;
//         ///v[i+1] = v[i] + k3_v*0.5;
//         ///x[i+1] = x[i] + k3_x*0.5;
//         k4_x =(v[i]+k3_v)*dt;
//         k4_v =-(x[i]+k3_x)*dt;
//         x[i+1] = x[i] + (k1_x+2.0*k2_x+2.0*k3_x+k4_x)/(6.0);
//         v[i+1] = v[i] + (k1_v+2.0*k2_v+2.0*k3_v+k4_v)/(6.0);
//         //scrittura dei dati sul file
// //         fprintf(fp, "%.10lf %.10lf %.10lf\n",t[i], x[i], v[i]);
//     }
// }

qreal R0(qreal r)
{
    qreal Z = 1;
    qreal Zs = Z-5.0/16;
    qreal a = 1; // r di bohr
    return 2*pow(Zs/a, 0.5)*Zs*r*exp(-Zs*r/a)/a;
}

hartreefock::hartreefock()
 : dx(0.1),
   xmax(100),
   m(1),
   hbar(1),
   Z(1),
   e(1)
{
    int steps = xmax/dx;

    // mesh r[i] = [dx -> xmax, delta=dx]

    m_rho.resize(steps);
    m_rho.fill(0.01);
    m_phi.resize(steps);
    m_phi.fill(0);
    m_chisq.resize(steps);
    m_chisq.fill(0);
    m_R.resize(steps);
    m_R.fill(0);

    stabilizeE();
}

hartreefock::~hartreefock()
{}

qreal hartreefock::phiIntegrand(int step, qreal x)
{
//     qDebug() << 4*3.141*x*x*m_rho.at(step)/abs(x-m_rho.at(step));
    return 4*3.141*x*x*m_rho.at(step)/fabs(x-m_rho.at(step));
}

QVector< qreal > hartreefock::updatePhi()
{
    QVector<qreal> phi;
    phi.resize(xmax/dx);

    qreal omega1 = ((qreal)128.)/225;
    qreal omega23 = ((322+13*sqrt((qreal)70))/900);
    qreal omega45 = ((322-13*sqrt((qreal)70))/900);

    qreal xi23 = ( ((qreal)1) /3)*sqrt(5-2*sqrt( (qreal)(10./7) ));
    qreal xi45 = ( ((qreal)1) /3)*sqrt(5+2*sqrt( (qreal)(10./7) ));

    for (int rstep = 0; rstep < m_rho.size(); rstep++) {

        qreal integral = 0;

        // integrate for every step
        for (int i = 0; i < 1000; i += 1) {
            qreal x_i = i;
            qreal x_i1 = i+1;
            qreal c = (x_i1+x_i)/2.;
            qreal m = (x_i1-x_i)/2.;

            // http://mathworld.wolfram.com/Legendre-GaussQuadrature.html
    //         qreal new
            integral += ( m*omega1*phiIntegrand(rstep, c) ); // root = 0
//             qDebug() << m*omega1*phiIntegrand(rstep, c);

            integral += ( m*omega23*phiIntegrand(rstep, c - m*xi23) );
            integral += ( m*omega23*phiIntegrand(rstep, c + m*xi23) );

            integral += ( m*omega45*phiIntegrand(rstep, c - m*xi45) );
            integral += ( m*omega45*phiIntegrand(rstep, c + m*xi45) );

//             if (x > 1) {
//                 return QVector<qreal>();
//             }
        }
//         qDebug() << integral;

        phi[rstep] = integral;
    }
    return phi;
}

qreal hartreefock::Veff(int i)
{
    return -0.5*m*(2*Z*e*e/((i+1)*dx)-m_phi(i))/hbar*hbar;
}

qreal hartreefock::doNumerov(qreal E, bool setR)
{
    QVector<qreal> r;
    for (qreal i = dx; i <= xmax; i+= dx) {
        r.append(i);
    }
    qreal N = r.size();

    QVector<qreal> forward; // ufb (:1)
    QVector<qreal> backward; // ufb (:2)
    forward.resize(r.size());
    backward.resize(r.size());

//     ufb(1,1) = m_R.first(); // R(0)
//     ufb(2,1) = m_R.first(); // R(dx)
//     ufb(N,2) = m_R.last(); // R(xmax)
//     ufb(N-1,2) = m_R.last(); // R(xmax-dx)

    forward[0] = m_R.first(); // R(dx)
    forward[1] = m_R.first(); // R(2dx)
    backward[r.size()-1] = exp(sqrt(2*E)*xmax); // R(xmax)
    backward[r.size()-2] = exp(sqrt(2*E)*(xmax-dx)); // R(xmax-dx)

    for (int i=2; i < N; i++) {
        qreal k3 = -2*( Veff(i) - E);
        qreal k2 = -2*( Veff(i-1) - E);
        qreal k1 = -2*( Veff(i-2) - E);
        forward[i] = 2*forward[i-1]*(1-(5.0/12)*k2*dx*dx)/(1+(1.0/12)*k3*dx*dx)
                     - forward[i-2]*(1+(1.0/12)*k1*dx*dx)/(1+(1.0/12)*k3*dx*dx);

        k3 = -2*( Veff(N-(i+1)) -E);
        k2 = -2*( Veff(N-(i+1) +1) -E);
        k1 = -2*( Veff(N-(i+1) +2) -E);
        backward[N-(i+1)] = 2*backward[N-(i+1) +1]*(1-(5.0/12)*k2*dx*dx)/(1+(1.0/12)*k3*dx*dx)
                            - backward[N-(i+1) +2]*(1+(1.0/12)*k1*dx*dx)/(1+(1.0/12)*k3*dx*dx);
    }

    qreal rc = xmax/(dx*2);
    qreal C = forward[rc]/backward[rc];

    if (!setR) {
        return C*(forward[rc-1]-backward[rc-1])/forward[rc-1];
    } else {
        for (int i = 0; i < m_R.size(); i++) {
            if (i < rc) {
                m_R[i] = forward.at(i)*C;
            } else {
                m_R[i] = backward.at(i)*C;
            }
        }
    }
    return 0;
}

void hartreefock::stabilizeE()
{
    m_phi = updatePhi();

    for (qreal e = 0.1; e < 10; e+=0.1) {
        qDebug() << doNumerov(e, false);
    }
//     qDebug() << m_phi.at(1);
//     for (int i = 0; i < m_R.size(); i++) {
//         m_chisq[i] = m_R.at(i)/i*dx;
//     }
}

qreal hartreefock::rho0(qreal r)
{
    return 1/r;
}

qreal hartreefock::phi(qreal r)
{
//     integrate(rho0())
}

qreal hartreefock::integrate(void* f)
{
    return 1;
}

void hartreefock::output()
{
//     std::cout << "Hello World!" << std::endl;
}

#include "hartreefock.moc"
