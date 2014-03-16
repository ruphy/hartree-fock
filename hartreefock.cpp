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

qreal R(qreal r)
{
    qreal Z = 1;
    qreal Zs = Z-5.0/16;
    qreal a = 1; // r di bohr
    return 2*pow(Zs/a, 0.5)*Zs*r*exp(-Zs*r/a)/a;
}

hartreefock::hartreefock()
 : dx(0.1),
   xmax(100)
{
    int steps = xmax/dx;

    m_rho.resize(steps); // 10000 is max radius
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

    for (int rstep = 0; rstep < m_rho.size()-1; rstep++) {

        qreal integral = 0;

//         qDebug() << rstep;
        for (qreal x = 0; x < xmax-dx; x += dx) {
            qreal x_i = x;
            qreal x_i1 = x+dx;
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
        qDebug() << integral;

        phi[rstep] = integral;
    }
    return phi;
}

void hartreefock::stabilizeE()
{
    m_phi = updatePhi();
    qDebug() << m_phi.at(1);
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
