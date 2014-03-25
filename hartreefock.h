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

private:
    qreal R0(qreal );
    qreal doNumerov(qreal E, bool setR);
    qreal Veff(int i) const;
    inline qreal phiIntegrand(int stepRPrime, qreal r) const;
    qreal calcNewE();
    qreal integratedRdr2();

    int factorial(int n) { return (n <= 1 ) ? 1 : factorial(n - 1) * n; }

    qreal findLowestEigenValue();

    qreal simpsonIntegrate(const QVector< qreal >& in, int step = -1) const;
    
    QVector<qreal> differenciate(const QVector< qreal >& in) const;
    QVector< qreal > normalize(const QVector< qreal >& vector) const;
    QVector< qreal > updateRho(bool average);
    QVector< qreal > updateManyRho(bool average);

    void printVector(const QVector< qreal >& vector) const;

    qreal iterateE(int iteration);

    void manyElectrons();

    void stabilizeE();

    QVector< qreal > updatePhi() const;

    QList< QList< QVector<qreal> > > m_Rnl;

    qreal energyForNL(int n, int l);

    qreal m_energy, m_eigenvalue, m_iteration;

    qreal meshseed;

    int n_elett;

    int xmax, Z,l;
    int m_steps;
    qreal m,hbar,e;

    TFile *debugFile;
    TTree *m_tree;

    QVector<qreal> m_rho;
    QVector<qreal> m_phi;
    QVector<qreal> m_chisq;
    QVector<qreal> m_R;
    QVector<qreal> m_ri;
    QVector<qreal> m_dx;

};

#endif // hartreefock_H
