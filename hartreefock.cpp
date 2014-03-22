#include "hartreefock.h"

#include <QTimer>
#include <QDebug>
#include <iostream>
#include <math.h>

qreal hartreefock::R0(qreal r)
{
    qreal Zs = Z-5.0/16;
    qreal a = 0.5299; // r di bohr
    return 2*pow(Zs/a, 0.5)*Zs*r*exp(-Zs*r/a)/a;
}

hartreefock::hartreefock()
{
    m = 1;
    l = 0; // ang. momentum
    hbar = sqrt(7.6359)/100;
    Z = 2;
    xmax = 10/Z;
    dx = 0.001/Z;
    e = sqrt(14.409)/10;
    m_steps = xmax/dx;

    qDebug() << xmax << dx << m_steps;
    // mesh r[i] = [dx -> xmax, delta=dx]


    m_rho.resize(m_steps);
    m_ri.resize(m_steps);

    for (int i = 0; i < m_steps; i++) {
        m_ri[i] = dx*(i+1);
    }

    m_phi.resize(m_steps);
    m_phi.fill(0);
    m_chisq.resize(m_steps);
    m_chisq.fill(0);
    m_R.resize(m_steps);
    m_R.fill(0);

    stabilizeE();
}

hartreefock::~hartreefock()
{}

qreal hartreefock::phiIntegrand(int step, qreal x) const
{
//     qDebug() << 4*3.141*x*x*m_rho.at(step)/abs(x-m_rho.at(step));
//     qDebug() << m_ri.at(step) << x;
//     return e*e*4*3.141*x*x*m_rho.at(step)/fabs(x-m_ri.at(step));
    if (fabs(x - m_ri.at(step)) < dx*dx) { // suppress too small values
        return 0;
    }

    return e*e*2*m_R.at(step)*m_R.at(step)/fabs(x-m_ri.at(step));
}



QVector< qreal > hartreefock::updatePhi() const
{
    QVector<qreal> phi;
    phi.resize(m_steps);

    qreal omega1 = ((qreal)128.)/225;
    qreal omega23 = ((322+13*sqrt((qreal)70))/900);
    qreal omega45 = ((322-13*sqrt((qreal)70))/900);

    qreal xi23 = ( ((qreal)1) /3)*sqrt(5-2*sqrt( (qreal)(10./7) ));
    qreal xi45 = ( ((qreal)1) /3)*sqrt(5+2*sqrt( (qreal)(10./7) ));

    for (int rstep = 0; rstep < m_rho.size(); rstep++) {

        qreal integral = 0;
//         qDebug()<< "Integrating phi(" << rstep << "); m_R=" << m_R[rstep];

        // integrate for every step
        for (qreal i = 0; i < xmax; i += 0.1) {
            qreal x_i = i;
            qreal x_i1 = i+1;
            qreal c = (x_i1+x_i)/2.;
            qreal m = (x_i1-x_i)/2.;

            // http://mathworld.wolfram.com/Legendre-GaussQuadrature.html
            integral += ( m*omega1*phiIntegrand(rstep, c) );

            integral += ( m*omega23*phiIntegrand(rstep, c - m*xi23) );
            integral += ( m*omega23*phiIntegrand(rstep, c + m*xi23) );

            integral += ( m*omega45*phiIntegrand(rstep, c - m*xi45) );
            integral += ( m*omega45*phiIntegrand(rstep, c + m*xi45) );
        }

        phi[rstep] = integral;
    }
    return phi;
}

qreal hartreefock::Veff(int i) const
{
//     qDebug() << "QEFF"  </*<*/ i;
    if (i > 1) {
        return l*(l+1)/(2*pow(m_ri[i],2)) - Z/m_ri[i];
    } else {
        return 0;
    }
//     return -0.5*m*(2*Z*e*e/m_ri[i]-m_phi[i])/hbar*hbar;
//     return -0.5*m*(2*Z*e*e/m_ri[i]-m_phi[i])/hbar*hbar;
}

qreal hartreefock::doNumerov(qreal E, bool setR)
{
    qreal N = m_ri.size();

    QVector<qreal> forward; // ufb (:1)
    QVector<qreal> backward; // ufb (:2)
    forward.resize(m_ri.size());
    backward.resize(m_ri.size());

    forward[0] = m_R.first(); // R(dx)
    forward[1] = 2*dx; // R(2dx)
    backward[m_ri.size()-1] = exp(sqrt(2*E)*xmax); // R(xmax)
    backward[m_ri.size()-2] = exp(sqrt(2*E)*(xmax-dx)); // R(xmax-dx)

    for (int i=2; i < N; i++) {
//         qDebug() << i;
        qreal k3 = -2*( Veff(i) - E);
        qreal k2 = -2*( Veff(i-1) - E);
        qreal k1 = -2*( Veff(i-2) - E);
        forward[i] = 2*forward[i-1]*(1-(5.0/12)*k2*dx*dx)/(1+(1.0/12)*k3*dx*dx)
                     - forward[i-2]*(1+(1.0/12)*k1*dx*dx)/(1+(1.0/12)*k3*dx*dx);


        k3 = -2*( Veff( (N-i) -1) -E);
        k2 = -2*( Veff((N-i) ) -E);
        k1 = -2*( Veff((N-i)+1) -E);
        backward[(N-i) -1] = 2.*backward[(N-i)]*(1.-(5.0/12)*k2*dx*dx)/(1.+(1.0/12)*k3*dx*dx)
                              - backward[(N-i) + 1]*(1.+(1.0/12)*k1*dx*dx)/(1.+(1.0/12)*k3*dx*dx);
    }

    qreal rc = 100;
    qreal C = forward[rc]/backward[rc];
//         qDebug() << "o forse allora qui?";

//     qDebug() << forward[rc];

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

qreal hartreefock::integratedRdr2()
{
    double k1_x, k2_x, k2_v, k3_v, k3_x, k4_v, k4_x, k1_v;
    QVector<qreal> der;
    qreal dt = dx;
    int n_cicli = m_R.size()-1;
    der.resize(m_R.size());
    der[0] = (m_R[1]-m_R[0])/dx; // first order approximation

    for(unsigned int i = 0; i < n_cicli; ++i){
        k1_x = der[i]*dt;
        k1_v = - m_R[i]*dt;
        /// der[i+1] = der[i] + k1_v*0.5;
        /// m_R[i+1] = m_R[i] + k1_x*0.5;
        k2_x = (der[i] + k1_v*0.5)*dt;
        k2_v = -(m_R[i] + k1_x*0.5)*dt;
        ///der[i+1] = der[i] + k2_v*0.5;
        ///m_R[i+1] = m_R[i] + k2_x*0.5;
        k3_x =(der[i]+k2_v*0.5)*dt;
        k3_v= -(m_R[i]+k2_x*0.5)*dt;
        ///der[i+1] = der[i] + k3_v*0.5;
        ///m_R[i+1] = m_R[i] + k3_x*0.5;
        k4_x =(der[i]+k3_v)*dt;
        k4_v =-(m_R[i]+k3_x)*dt;
        m_R[i+1] = m_R[i] + (k1_x+2.0*k2_x+2.0*k3_x+k4_x)/(6.0);
        der[i+1] = der[i] + (k1_v+2.0*k2_v+2.0*k3_v+k4_v)/(6.0);
    }


    qreal result = 0;
    foreach (const qreal el, der) {
        result += el*el;
    }
    result *= dx;
    return result;
}

qreal hartreefock::calcNewE()
{
    qreal Energy = hbar*integratedRdr2()/m;
    qDebug() << Energy;
    qreal potentialPart = 0;
    for (int i= 0; i < m_R.size(); i++) {
        qreal uno = -Z*e*e/m_ri[i] + m_phi[i]/4.;
        potentialPart += uno*m_rho[i]*4*3.414*m_ri[i]*m_ri[i];
    }
    Energy += potentialPart*dx;

    return Energy;
}

QVector<qreal> hartreefock::normalize(const QVector< qreal > &vector) const
{
    qDebug() << "normalizing";
    QVector<qreal> v = vector;
    qreal acc = 0;
    foreach (const qreal el, v) {
        acc += el*el;
    }
    acc *= dx;
    acc = sqrt(acc);

    qDebug() << "norm" << acc;
    for (int i = 0; i < v.size(); i++) {
        v[i] = v.at(i)/acc;
    }

    return v;
}

void hartreefock::updateRho()
{
    for (int i = 0; i < m_R.size(); i++) {
        m_rho[i] = pow(m_R.at(i)/m_ri[i], 2)/(4*3.14);
    }
    // Normalize rho
    normalize(m_rho);
}

void hartreefock::printVector(const QVector< qreal >& vector) const
{
    qDebug() << "debug started";
    for (int rstep = 0; rstep < m_ri.size(); rstep++) {
        std::cout << m_ri[rstep]<< ","<<m_R[rstep] << std::endl;
    }
    qDebug() << "--- debug ended";
}

void hartreefock::stabilizeE()
{
    // Calculate R
    for (int i = 0; i < m_R.size(); i++) {
        m_R[i] = R0(m_ri[i]);
    }
    m_R = normalize(m_R);
    // Calculate rho
    updateRho();
    m_phi = updatePhi();
    qDebug() << calcNewE();

    QList<qreal> eigenvalues;

    qreal increment = 0.01;
    for (qreal e = 0.1; e < 2; e+=increment) {
        qreal giusto = doNumerov(e, false);
        if (fabs(giusto) < 1) {
            for (qreal ne = e; ne < e+increment; ne += increment*0.01) {
                qreal ngiusto = doNumerov(ne, false);
                if (fabs(ngiusto) < 0.052) {
                    qreal min = 10;
                    qreal mine = 0;
                    for (qreal nne = ne; nne < ne+increment*0.01; nne += increment*0.0001) {
                        qreal nngiusto = doNumerov(nne, false);
                        if (nngiusto < min) {
                            min = nngiusto;
                            mine = nne;
                        }
//                         std::cout << nne << ',' << nngiusto <<std::endl;
                    }
                    if (!eigenvalues.size() or fabs(eigenvalues.last() - mine) > 0.01) {
                        eigenvalues.append(mine);
                    }
                }
//                 std::cout << ne << ',' << ngiusto <<std::endl;
            }
        }
//         std::cout << e << ',' << giusto <<std::endl;
    }
    qDebug() << eigenvalues;

    foreach (qreal eigen, eigenvalues) {
        doNumerov(eigen,true);
        m_R = normalize(m_R);
        // Calculate rho
        updateRho();
        m_R = normalize(m_R);
        printVector(m_R);
        m_phi = updatePhi();
//         qDebug()<< m_phi;
        qDebug() << "eigen" << eigen<< "calc" << calcNewE();
    }
    qDebug() << "end";
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
