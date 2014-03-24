#ifndef hartreefock_H
#define hartreefock_H

#include <QtCore>

class TFile;
class TTree;

class hartreefock : public QObject
{
Q_OBJECT
public:
    hartreefock();
    virtual ~hartreefock();

    qreal rho0(qreal r);
    qreal phi(qreal r);

    qreal integrate(void *f);


private slots:
    void output();

private:
    qreal R0(qreal );
    qreal doNumerov(qreal E, bool setR);
    qreal Veff(int i) const;
    inline qreal phiIntegrand(int stepRPrime, qreal r) const;
    qreal calcNewE();
    qreal integratedRdr2();

    qreal findLowestEigenValue();

    QVector<qreal> differenciate(const QVector< qreal >& in) const;
    QVector< qreal > normalize(const QVector< qreal >& vector) const;
    QVector< qreal > updateRho() const;

    void printVector(const QVector< qreal >& vector) const;

    qreal iterateE(int iteration);

    void stabilizeE();

    QVector< qreal > updatePhi() const;

    qreal m_energy, m_eigenvalue, m_iteration;

    qreal dx;
    int xmax, Z,l; // steps = xmax/dx
    int m_steps;
    qreal m,hbar,e;

#ifdef ROOT
    TFile *debugFile;
    TTree *m_tree;
#endif
    
    QVector<qreal> m_rho;
    QVector<qreal> m_phi;
    QVector<qreal> m_chisq;
    QVector<qreal> m_R;
    QVector<qreal> m_ri;
    QVector<qreal> m_dx;

};

#endif // hartreefock_H
