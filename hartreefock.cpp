#include "hartreefock.h"

#include <QTimer>
#include <QDebug>
#include <iostream>
#include <math.h>

#include <root/TFile.h>
#include <root/TTree.h>

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
    n_elett = 2;
    meshseed = 0.025/Z;
    e = sqrt(14.409);
    m_steps = xmax/meshseed;

    debugFile = new TFile("out.root", "RECREATE", "An Example ROOT file");
    m_tree = new TTree("bTree", "tree title");

    qDebug() << xmax << meshseed << m_steps;

    m_rho.resize(m_steps);
    m_ri.resize(m_steps);
    m_dx.resize(m_steps);

    for (int i = 0; i < m_steps; i++) {
        qreal x = meshseed + meshseed * i;
        m_dx[i] = meshseed*exp(x)/Z;
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
    for (int i = 0; i < m_R.size(); i++) {
        m_R[i] = R0(m_ri[i]); // Inizializazione con wavefunction dell'orbitale idrogenico
    }
    m_rho = updateRho(false); // Calculate rho

    m_phi = updatePhi();

    printVector(m_phi);

    qDebug() << "Phi != 0";
    qDebug() << calcNewE();

    for (int i = 0; i < m_R.size(); i++) {
        m_phi[i] = 0; // Inizializazione con wavefunction dell'orbitale idrogenico
    }
    qDebug() << "Phi = 0";
    qDebug() << calcNewE();

    m_tree->Branch("energia", &m_energy);
    m_tree->Branch("eigenvalue", &m_eigenvalue);
    m_tree->Branch("iterazione", &m_iteration);

    for (int i = 0; i < 25; i++) {
        m_iteration = i;
        qDebug() << "--------";
        qDebug() << "Iterando sugli autovalori... iterazione" << i;
        qDebug() << "ENERGIA: " << iterateE(i);

        m_tree->Fill();
        debugFile->Write();
    }
}

QVector< qreal > hartreefock::updateRho(bool average = false)
{
    QVector<qreal> rho;
    rho.resize(m_R.size());

    for (int i = 0; i < m_R.size(); i++) {
        qreal newterm = n_elett*pow(m_R.at(i)/m_ri.at(i), 2)/(4*PI);
        if (average) {
            rho[i] = m_rho.at(i)*0.5+newterm*0.5;
        } else {
            rho[i] = newterm;
        }
    }
    return rho;
}

// r' <= stepRPrime, qreal r = r
qreal hartreefock::phiIntegrand(int stepRPrime, qreal r) const
{
    if (r == m_ri.at(stepRPrime)) { // suppress inf values
        return 0;
    }

    return e*e*pow(m_ri.at(stepRPrime),2)*
           m_rho.at(stepRPrime)/fabs(r-m_ri.at(stepRPrime));
}

qreal hartreefock::simpsonIntegrate(const QVector< qreal >& in, int step) const
{
    int maxstep = step;
    if (maxstep == -1 || maxstep > in.size()) {
        maxstep = in.size();
    }
    qreal risultato = 0;
    for (int i = 0; i < maxstep-5; i+=4) {
        qreal tempint = 0;
        qreal dx = (m_ri.at(i+4)-m_ri.at(i))*2/45.;
        tempint += in.at(i)*7;
        tempint += in.at(i+1)*32;
        tempint += in.at(i+2)*12;
        tempint += in.at(i+3)*32;
        tempint += in.at(i+4)*7;
        tempint *= (dx/5);
        risultato += tempint;
    }

    return risultato;
}

QVector< qreal > hartreefock::updatePhi() const
{
    QVector<qreal> phi;
    phi.resize(m_steps);

    for (int rstep = 0; rstep < m_rho.size(); rstep++) {
        qreal integral = 0;
        // integrate for every step
        for (int i = 0; i < m_rho.size()-4; i+=4) {

            qreal tempint = 0;
            qreal dx = (m_ri.at(i+4)-m_ri.at(i))*2/45.;

            tempint += phiIntegrand(i, m_ri[rstep])*7;
            tempint += phiIntegrand(i+1, m_ri[rstep])*32;
            tempint += phiIntegrand(i+2, m_ri[rstep])*12;
            tempint += phiIntegrand(i+3, m_ri[rstep])*32;
            tempint += phiIntegrand(i+4, m_ri[rstep])*7;

            integral += tempint*(dx/5);
        }
        phi[rstep] = integral;
    }
    return phi;
}

qreal hartreefock::calcNewE()
{
    qreal Energy = hbar*hbar*integratedRdr2()/m;
    qDebug() << "KE" << Energy;

    qreal p1 = 0;
    qreal p2 = 0;
    qreal Zs = Z;

    QVector<qreal> centrifugal(m_steps);
    QVector<qreal> electrostatical(m_steps);

    for (int i = 0; i < m_R.size(); i++) {
        electrostatical[i] = -Zs*e*e*m_rho[i]*4*PI*m_ri[i];
        centrifugal[i] = m_rho[i]*PI*m_ri[i]*m_ri[i]*m_phi[i];
    }

    p1 = simpsonIntegrate(electrostatical);
    p2 = simpsonIntegrate(centrifugal);

    qDebug() << "Electrostatic" << p1;
    qDebug() << "Centrifugal" << p2;
    Energy += p1+p2;

    return Energy;
}

qreal hartreefock::energyForNL(int n, int l)
{
    qreal energy = 0;
    QVector< qreal > Rnl = m_Rnl[n][l];

    // E1
    QVector<qreal> one = differenciate(Rnl);
    for (int i = 0; i < one.size(); i++) {
        one[i] = pow(one.at(i), 2);
        one[i] += l*(l+1)*pow(Rnl.at(i),2)/pow(m_ri.at(i),2);
    }
    energy += simpsonIntegrate(one)/2.;

    // E2
    QVector<qreal> integrand(Rnl.size());
    for (int i = 0; i < Rnl.size(); i++) {
        integrand[i] = pow(Rnl.at(i),2);
        integrand[i] *= (-Z/m_ri.at(i)+m_phi.at(i));
    }
    energy += simpsonIntegrate(integrand);



    for (int i = 0; i < Rnl.size(); i++) {

        integrand[i] = Rnl.at(i);
        qreal sum = 0;

        for (int np = 0; np < n; n++) {
            for (int lp = 0; lp < l; l++) {

                qreal Nnpnl = pow(np,2)/(2*lp+1);
                qreal A = 0;

                QVector<qreal> Rnplp = m_Rnl.at(np).at(lp);

                // (A)
                for (int lambda = abs(l-lp); lambda <= l+lp; lambda++) {
                    qreal P = (l+lambda+lp)/2.;
                    // A1
                    qreal A1 = factorial(-l+lp+lambda)*factorial(l-lp+lambda)*factorial(l+lp-lambda);
                    A1 /= (factorial(l+lp+lambda+1)*factorial(P-l)*factorial(P-lp)*factorial(P-lambda));
                    A1 *= A1;

                    QVector<qreal> partIntegrand(i); // will contain only values 0 -> r
                    for (int j = 0; j < i; j++) { // j -> rp
                        partIntegrand[j] = Rnplp.at(j)*
                                           Rnl.at(j)*pow(m_ri.at(j), lambda);
                    }
                    qreal A2 = simpsonIntegrate(partIntegrand);

                    qreal A3 = pow(m_ri.at(i), lambda);

                    QVector<qreal> partIntegrand2(m_Rnl.size());
                    for (int j = 0; j < m_Rnl.size(); j++) { // j -> rp
                        partIntegrand2[j] = Rnplp.at(i)*
                                            Rnl.at(j)/pow(m_ri.at(j), lambda+1);
                    }

                    A3 *= simpsonIntegrate(partIntegrand2);

                    A += A1*A2/pow(m_ri.at(i), lambda+1) + A3;
                }

                sum += -0.5*Rnl.at(i)*Nnpnl*A;
            }
        }

        integrand[i] *= sum;
    }

    energy += simpsonIntegrate(integrand);

    return energy;
}

qreal hartreefock::findLowestEigenValue()
{
    QList<qreal> eigenvalues;

    qreal increment = 0.01;
    for (qreal e = -105; e < 10; e+=increment) {
        qreal giusto = doNumerov(e, false);
        if (fabs(giusto) < 0.1) {
            for (qreal ne = e; ne < e+increment; ne += increment*0.01) {
                qreal ngiusto = doNumerov(ne, false);
                if (fabs(ngiusto) < 0.01) {
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
            }
        }
    }
    qDebug() << eigenvalues;
    return eigenvalues.first();
}

qreal hartreefock::iterateE(int iteration)
{
    m_eigenvalue = findLowestEigenValue();

    doNumerov(m_eigenvalue,true);
    m_R = normalize(m_R);
    m_rho = updateRho(true);
    m_phi = updatePhi();

    m_energy = calcNewE();

    qDebug() << "EigenValue : " << m_eigenvalue;
    return m_energy;
}

qreal hartreefock::Veff(int i) const
{
    if (i > 1) {
        return l*(l+1)/(2*pow(m_ri[i],2)) - Z/m_ri[i];
    } else {
        return 0;
    }
    return -0.5*m*(2*Z*e*e/m_ri[i]-m_phi[i])/(hbar*hbar);
}

qreal hartreefock::doNumerov(qreal E, bool setR)
{
    qreal N = m_ri.size();

    QVector<qreal> forward;
    QVector<qreal> backward;
    forward.resize(m_ri.size());
    backward.resize(m_ri.size());

    forward[0] = m_R.first(); // R(dx)
    forward[1] = m_dx.at(0); // R(2dx)
    backward[m_ri.size()-1] = exp(-sqrt(fabs(2*E))*xmax); // R(xmax)
    backward[m_ri.size()-2] = exp(-sqrt(fabs(2*E))*(xmax-m_dx.last())); // R(xmax-dx)

    for (int i=2; i < N; i++) {
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

void hartreefock::manyElectrons()
{
    int n = pow(n_elett,2)/(2*l+1);
    for (int i = 0; i < n; i++) {
        QList<QVector<qreal> > newL;
        for (int j = 0; i < l; i++) {
            QVector<qreal> newRNL = m_R;
            for (int i = 0; i < m_R.size(); i++) {
                newRNL[i] = R0(m_ri[i]);
            }
            newL.append(m_R);
        }
        m_Rnl.append(newL);
    }
    m_rho = updateManyRho(false);
    m_phi = updatePhi();
    qreal eigen1=0;
    for (int i = 0; i < 100; i++) {
        m_eigenvalue = findLowestEigenValue();

        doNumerov(m_eigenvalue,true);
        m_R = normalize(m_R);
        m_rho = updateManyRho(true);
        m_phi = updatePhi();

        m_energy = 0;
        for (int ni = 0; ni < n; ni++) {
            for (int li=0; li<l; li++) {
                m_energy += energyForNL(ni, li);
            }
        }

        qDebug() << "EigenValue : " << m_eigenvalue;
        if (fabs(eigen1 - m_eigenvalue) < 1e-6) {
            qDebug() << "energy converged to" << m_energy;
        }
    }
}


QVector< qreal > hartreefock::updateManyRho(bool average = false)
{
    QVector<qreal> rho;
    rho.resize(m_R.size());

    for (int ni = 0; ni < pow(n_elett,2)/(2*l+1); ni++) {
        for (int li=0; li<l; li++) {
            for (int i = 0; i < m_R.size(); i++) {
                qreal newterm = 2*pow(m_Rnl.at(ni).at(li).at(i)/m_ri.at(i), 2)/(4*PI);
                if (average) {
                    rho[i] = m_rho.at(i)*0.5+newterm*0.5;
                } else {
                    rho[i] = newterm;
                }
            }
        }
    }
    return rho;
}

QVector< qreal > hartreefock::differenciate(const QVector< qreal > &in) const
{
    QVector<qreal> der;
    der.resize(in.size());
    for (int i = 0; i < in.size()-1; i++) {
        der[i] = (in.at(i+1) - in.at(i))/m_dx.at(i);
    }

    return der;
}

qreal hartreefock::integratedRdr2()
{
    QVector<qreal> der = differenciate(m_R);
    for (int i = 0; i < der.size(); i++) {
        der[i] = pow(der.at(i), 2);
    }
    return simpsonIntegrate(der);
}

QVector<qreal> hartreefock::normalize(const QVector< qreal > &vector) const
{
    qDebug() << "normalizing";
    QVector<qreal> v = vector;
    qreal acc = 0;

    for (int i = 0; i < vector.size(); i++) {
        acc += vector.at(i)*vector.at(i)*(m_dx.at(i)+m_dx.at(i))/2;
    }

//     acc = sqrt(simpsonIntegrate(v));

//     foreach (const qreal el, v) {
//         acc += el*el;
//     }
//     acc *= dx;

    acc = sqrt(acc);

    qDebug() << "norm" << acc;
    for (int i = 0; i < v.size(); i++) {
        v[i] = vector.at(i)/acc;
    }

    return v;
}

void hartreefock::printVector(const QVector< qreal >& vector) const
{
    QString name = "vec";
    double x, y;

    m_tree->Branch("x", &x);
    m_tree->Branch("y", &y);
//     m_tree->Branch(name.toAscii(), &g);

    for (int rstep = 0; rstep < m_ri.size(); rstep++) {
        x = m_ri.at(rstep);
        y = vector.at(rstep)*m_rho.at(rstep);
        m_tree->Fill();
    }
    debugFile->Write();
}


#include "hartreefock.moc"
