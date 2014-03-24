#include "hartreefock.h"

#include <QTimer>
#include <QDebug>
#include <iostream>
#include <math.h>

#ifdef ROOT
#include <root/TFile.h>
#include <root/TTree.h>
#endif
// #include <root/TGraph2D.h>
// #include <root/TGraph.h>

//15943.pts-2.omero

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
    dx = 0.00005/Z;
    e = sqrt(14.409);
    m_steps = xmax/dx;

#ifdef ROOT
    debugFile = new TFile("out.root", "RECREATE", "An Example ROOT file");
    m_tree = new TTree("bTree", "tree title");
#endif

    qDebug() << xmax << dx << m_steps;
    // mesh r[i] = [dx -> xmax, delta=dx]

    m_rho.resize(m_steps);
    m_ri.resize(m_steps);
    m_dx.resize(m_steps+1);

    for (int i = 0; i < m_steps+1; i++) {
        m_dx[i] = dx;//*log(i+2);
//         std::cout << dx*log(i) << std::endl;
    }

    for (int i = 0; i < m_steps; i++) {
        m_ri[i] = m_dx.at(i)*(i+1);
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

    qDebug() << "Phi != 0";

    qDebug() << calcNewE();


//     for (int i = 0; i < m_R.size(); i++) {
//         m_phi[i] = 0; // Inizializazione con wavefunction dell'orbitale idrogenico
//     }
//     qDebug() << "Phi = 0";
//     qDebug() << calcNewE();

#ifdef ROOT
    m_tree->Branch("energia", &m_energy);
    m_tree->Branch("eigenvalue", &m_eigenvalue);
    m_tree->Branch("iterazione", &m_iteration);
#endif
    for (int i = 0; i < 500; i++) {
        m_iteration = i;
        qDebug() << "--------";
        qDebug() << "Iterando sugli autovalori... iterazione" << i;
        qDebug() << "ENERGIA: " << iterateE(i);

#ifdef ROOT
        m_tree->Fill();
        debugFile->Write();
#endif
    }
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
    if (r == m_ri.at(stepRPrime)) { // suppress too small values
        return 0;
    }

    // FIXME FANCULO MI SONO DIMENTICATO UN 4PI
    return e*e*pow(m_ri.at(stepRPrime),2)*
           m_rho.at(stepRPrime)/fabs(r-m_ri.at(stepRPrime));
}

QVector< qreal > hartreefock::updatePhi() const
{
    QVector<qreal> phi;
    phi.resize(m_steps);

    for (int rstep = 0; rstep < m_rho.size(); rstep++) {
        qreal integral = 0;
        // integrate for every step
        for (int i = 0; i < m_rho.size()-2; i++) {
            qreal tempint = 0;
            // FIXME calculate the right b-a
            tempint += phiIntegrand(i, m_ri[rstep])*m_dx.at(i);
            tempint += phiIntegrand(i+1, m_ri[rstep])*4*m_dx.at(i+1);
            tempint += phiIntegrand(i+2, m_ri[rstep])*m_dx.at(i+2);
            integral += tempint/6.;
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
    qreal Zs = Z;//-5.0/16;

    for (int i = 0; i < m_R.size()-2; i++) {

        qreal uno = -Zs*e*e*m_rho[i]*4*PI*m_ri[i]*m_dx.at(i);
        uno += -Zs*e*e*m_rho[i+1]*4*PI*m_ri[i+1]*4*m_dx.at(i+1);
        uno += -Zs*e*e*m_rho[i+2]*4*PI*m_ri[i+2]*m_dx.at(i+2);

        qreal due = m_rho[i]*PI*m_ri[i]*m_ri[i]*m_phi[i]*m_dx.at(i);
        due += m_rho[i+1]*PI*m_ri[i+1]*m_ri[i+1]*m_phi[i+1]*4*m_dx.at(i+1);
        due += m_rho[i+2]*PI*m_ri[i+2]*m_ri[i+2]*m_phi[i+2]*m_dx.at(i+2);

        p1 += uno/6.;
        p2 += due/6.;
    }

    qDebug() << "Electrostatic" << p1;
    qDebug() << "Centrifugal" << p2;
    Energy += p1+p2;

    return Energy;
}

qreal hartreefock::findLowestEigenValue()
{

    QList<qreal> eigenvalues;

    qreal increment = 0.1;
    for (qreal e = -15; e < 0; e+=increment) {
        qreal giusto = doNumerov(e, false);
        if (fabs(giusto) < 0.1) {
            for (qreal ne = e; ne < e+increment; ne += increment*0.01) {
                qreal ngiusto = doNumerov(ne, false);
                if (fabs(ngiusto) < 0.0001) {
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
                        return mine;
                        eigenvalues.append(mine);
                    }

                }

//                 std::cout << ne << ',' << ngiusto <<std::endl;
            }
        }
    }
    qDebug() << eigenvalues;
}

qreal hartreefock::iterateE(int iteration)
{
    // WARNING this will do only one cycle
    m_eigenvalue = findLowestEigenValue();

    doNumerov(m_eigenvalue,true);
    m_R = normalize(m_R);
    m_rho = updateRho();
    m_phi = updatePhi();

    m_energy = calcNewE();

    return m_energy;
}

qreal hartreefock::Veff(int i) const
{
//     qDebug() << "QEFF"  </*<*/ i;
//     return m_phi[i]/2;

//     if (i > 1) {
//         return l*(l+1)/(2*pow(m_ri[i],2)) - Z/m_ri[i];
//     } else {
//         return 0;
//     }
    return -0.5*m*(2*Z*e*e/m_ri[i]-m_phi[i])/(hbar*hbar);
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
    backward[m_ri.size()-2] = exp(-sqrt(fabs(2*E))*(xmax-m_dx.last())); // R(xmax-dx)

    for (int i=2; i < N; i++) {
//         qDebug() << i;

        qreal k3 = -2*( Veff(i) - E)*pow(m_dx.at(i),2);
        qreal k2 = -2*( Veff(i-1) - E)*pow(m_dx.at(i-1),2);
        qreal k1 = -2*( Veff(i-2) - E)*pow(m_dx.at(i-2),2);

        forward[i] = 2*forward[i-1]*(1-(5.0/12)*k2)/(1+(1.0/12)*k3)
                     - forward[i-2]*(1+(1.0/12)*k1)/(1+(1.0/12)*k3);


        k3 = -2*( Veff( (N-i) -1) -E)*pow(m_dx.at((N-i) -1),2);
        k2 = -2*( Veff((N-i) ) -E)*pow(m_dx.at(N-i),2);
        k1 = -2*( Veff((N-i)+1) -E)*pow(m_dx.at((N-i) +1),2);
        backward[(N-i) -1] = 2.*backward[(N-i)]*(1.-(5.0/12)*k2)/(1.+(1.0/12)*k3)
                              - backward[(N-i) + 1]*(1.+(1.0/12)*k1)/(1.+(1.0/12)*k3);
    }

    int rc = m_steps/2;
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
        der[i] = (in.at(i+1) - in.at(i))/m_dx.at(i);
    }

    der.append(der.last());
    return der;
}

qreal hartreefock::integratedRdr2()
{
    QVector<qreal> der = differenciate(m_R);
    qreal result = 0;
    for (int i = 0; i < der.size(); i++) {
        result += pow(der.at(i), 2)*m_dx.at(i);
    }
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
#ifdef ROOT
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

#endif
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
