#include <QCoreApplication>
#include "hartreefock.h"


int main(int argc, char** argv)
{
    QCoreApplication app(argc, argv);
    hartreefock foo;
    return app.exec();
}
