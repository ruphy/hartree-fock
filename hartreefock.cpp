#include "hartreefock.h"

#include <QTimer>
#include <QDebug>
#include <iostream>
#include <math.h>

#include <root/TFile.h>
#include <root/TTree.h>
// #include <root/TGraph2D.h>
// #include <root/TGraph.h>

#define PI 3.14159265359

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
    hbar = sqrt(7.6359);
    Z = 1;
    xmax = 10/Z;
    dx = 0.001/Z;
    e = sqrt(14.409);
    m_steps = xmax/dx;

    debugFile = new TFile("out.root", "RECREATE", "An Example ROOT file");
    m_tree = new TTree("bTree", "tree title");

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

void hartreefock::stabilizeE()
{
    // Nota: R ~= wavefunction

    for (int i = 0; i < m_R.size(); i++) {
        m_R[i] = R0(m_ri[i]); // Inizializazione con wavefunction dell'orbitale idrogenico
    }

    m_R = normalize(m_R); // questa riga non dovrebbe servire, ma male non fa
    m_rho = updateRho(); // Calculate rho

    m_phi = updatePhi(); // Phi(r) = \int phiIntegrand(r') dr'

    printVector(m_phi);

    qDebug() << m_phi[0];
    qDebug() << calcNewE();

//     m_phi = normalize(m_phi);

//     for (int i = 0; i < 100; i++) {
//         qDebug() << iterateE();
//     }
}

QVector< qreal > hartreefock::updateRho() const
{
    QVector<qreal> rho;
    rho.resize(m_R.size());
    for (int i = 0; i < m_R.size(); i++) {
        rho[i] = 2*pow(m_R.at(i)/m_ri.at(i), 2)/(4*PI);
//         qDebug() << "rho(" <<m_ri.at(i) << ")=" << rho.at(i);
    }
    return rho;
}

// r' <= stepRPrime, qreal r = r
qreal hartreefock::phiIntegrand(int stepRPrime, qreal r) const
{
    if (!fabs(r - m_ri.at(stepRPrime))) { // suppress too small values
        return 0;
    }

    return e*e*4*PI*
           pow(m_ri.at(stepRPrime),2)*
           m_rho.at(stepRPrime)/fabs(r-m_ri.at(stepRPrime));
}

QVector< qreal > hartreefock::updatePhi() const
{
    QVector<qreal> phi;
    phi.resize(m_steps);

    for (int rstep = 0; rstep < m_rho.size(); rstep++) {
        qreal integral = 0;
        // integrate for every step
        for (int i = 0; i < m_rho.size(); i++) {
            integral += phiIntegrand(i, m_ri[rstep])*dx;
        }
        phi[rstep] = integral;
    }
    return phi;
}

qreal hartreefock::calcNewE()
{
    qreal Energy = hbar*hbar*integratedRdr2()/m;
    qDebug() << "KE" << Energy;

//     qDebug() << m_rho;
    qreal p1 = 0;
    qreal p2 = 0;

    for (int i = 0; i < m_R.size(); i++) {
        qreal uno = -Z*e*e*m_rho[i]*4*PI*m_ri[i];
        qreal due = m_rho[i]*PI*m_ri[i]*m_ri[i]*m_phi[i];
        p1 += uno*dx;
        p2 += due*dx;
    }
    qDebug() << "Electrostatic" << p1;
    qDebug() << "Centrifugal" << p2;
    Energy += p1+p2;

    return Energy;
}


qreal hartreefock::Veff(int i) const
{
//     qDebug() << "QEFF"  </*<*/ i;
//     return m_phi[i]/2;

//     if (i > 1) {
        return l*(l+1)/(2*pow(m_ri[i],2)) - Z/m_ri[i];
//     } else {
//         return 0;
//     }
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
    forward[1] = (m_R.at(2)-m_R.at(1))/dx; // R(2dx)
    backward[m_ri.size()-1] = exp(-sqrt(fabs(2*E))*xmax); // R(xmax)
    backward[m_ri.size()-2] = exp(-sqrt(fabs(2*E))*(xmax-dx)); // R(xmax-dx)

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

    if (!setR) {
        return (-forward[rc-1]+C*backward[rc-1])/forward[rc-1];
    } else {
        for (int i = 0; i < m_R.size(); i++) {
            if (i < rc) {
                m_R[i] = forward.at(i);
            } else {
                m_R[i] = backward.at(i)*C;
            }
        }
    }
    return 0;
}

QVector< qreal > hartreefock::differenciate(const QVector< qreal > &in) const
{
    QVector<qreal> der;
    der.resize(in.size());
    for (int i = 0; i < in.size()-1; i++) {
        der[i] = (in.at(i+1) - in.at(i))/dx;
    }

    der.append(der.last());
    return der;
}

qreal hartreefock::integratedRdr2()
{
    QVector<qreal> der = differenciate(m_R);
    qreal result = 0;
    foreach (const qreal el, der) {
        result += el*el;
    }
    result *= dx;
    return result;
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

void hartreefock::printVector(const QVector< qreal >& vector) const
{

    QString name = "vec";
    double x, y;
//     TGraph g(m_ri.size());

    m_tree->Branch("x", &x);
    m_tree->Branch("y", &y);
//     m_tree->Branch(name.toAscii(), &g);

    qDebug() << "debug started";
    for (int rstep = 0; rstep < m_ri.size(); rstep++) {
//         std::cout << m_ri[rstep]<< ","<<m_R[rstep] << std::endl;
//         g.SetPoint(rstep, m_ri.at(rstep), vector.at(rstep));
        x = m_ri.at(rstep);
        y = vector.at(rstep)*m_rho.at(rstep);
        m_tree->Fill();
    }
//     g.Draw("AC*");
    qDebug() << "--- debug ended";
    debugFile->Write();
}

qreal hartreefock::iterateE()
{
    // Calculate R

//     qDebug() << calcNewE();

    QList<qreal> eigenvalues;

    qreal increment = 0.01;
    for (qreal e = -10; e < -9.8; e+=increment) {
        qreal giusto = doNumerov(e, false);
        if (fabs(giusto) < 0.1) {
            for (qreal ne = e; ne < e+increment; ne += increment*0.01) {
                qreal ngiusto = doNumerov(ne, false);
                if (fabs(ngiusto) < 0.0052) {
                    qreal min = 10;
                    qreal mine = 0;
                    for (qreal nne = ne; nne < ne+increment*0.01; nne += increment*0.0001) {
                        qreal nngiusto = doNumerov(nne, false);
                        if (nngiusto < min) {
                            min = nngiusto;
                            mine = nne;
                        }
                    }
                    if (!eigenvalues.size() or fabs(eigenvalues.last() - mine) > 0.01) {
                        eigenvalues.append(mine);
                    }
                }
            }
        }
        std::cout << e << ',' << giusto <<std::endl;
    }
    qDebug() << eigenvalues;
//     printVector(m_R);

    // WARNING this will do only one cycle
    foreach (qreal eigen, eigenvalues) {
        doNumerov(eigen,true);
        m_R = normalize(m_R);
        updateRho();

//         printVector(m_R);

        m_phi = updatePhi();
        m_phi = normalize(m_phi);
        return calcNewE();
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
