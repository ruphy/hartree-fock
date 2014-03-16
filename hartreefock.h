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
    qreal phiIntegrand(int step, qreal x);

    void stabilizeE();
    QVector< qreal > updatePhi();

    const qreal dx;
    const int xmax; // steps = xmax/dx

    QVector<qreal> m_rho;
    QVector<qreal> m_phi;
    QVector<qreal> m_chisq;
    QVector<qreal> m_R;

};

#endif // hartreefock_H
