#include "mainwindow.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    a.setApplicationDisplayName("OpenCL image transformation");
    a.setApplicationName("OpenCL in Qt");
    MainWindow w;
    w.setWindowIcon(QIcon(":/images/resources/system icon.png"));
    w.setWindowTitle("OpenCL in Qt");
    w.show();

    return a.exec();
}
