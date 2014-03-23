#ifndef hartreefock_H
#define hartreefock_H

#include <QtCore>

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
    inline qreal phiIntegrand(int step, qreal x) const;
    qreal calcNewE();
    qreal integratedRdr2();

    QVector< qreal > normalize(const QVector< qreal >& vector) const;
    void updateRho();

    void printVector(const QVector< qreal >& vector) const;

    qreal iterateE();

    void stabilizeE();

    QVector< qreal > updatePhi() const;

    qreal dx;
    int xmax, Z,l; // steps = xmax/dx
    int m_steps;
    qreal m,hbar,e;

    QVector<qreal> m_rho;
    QVector<qreal> m_phi;
    QVector<qreal> m_chisq;
    QVector<qreal> m_R;
    QVector<qreal> m_ri;

};

#endif // hartreefock_H
