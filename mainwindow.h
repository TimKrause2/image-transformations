#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QPixmap>
#include <QGraphicsScene>
#include <QGraphicsSimpleTextItem>
#include <QFile>
#include <QKeyEvent>
#include <QThread>
#include <CL/cl.h>

namespace Ui {
class MainWindow;
}

class WorkerThread : public QThread
{
    Q_OBJECT
    void run() override {
        cl_int r = clFinish(queue);
        if(r!=CL_SUCCESS){
            qDebug("WorkerThread clFinish failed error:%d",r);
        }
    }
public:
    WorkerThread(cl_command_queue queue) : queue(queue){}

private:
    cl_command_queue queue;
};

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private:
    Ui::MainWindow *ui;
    QPixmap srcPixmap;
    QPixmap dstPixmap;
    QImage srcImage;
    QImage srcImage_argb;
    char *dstImage_data;
    QFile clFile;
    char *clData;
    QGraphicsScene *scene;
    QGraphicsPixmapItem *leftPixmapItem;
    QGraphicsPixmapItem *rightPixmapItem;
    QTimer *timer;

    cl_uint num_platforms;
    cl_platform_id *platforms;
    cl_uint num_devices;
    cl_device_id *devices;
    cl_context context;
    cl_command_queue queue;
    cl_program program;
    cl_kernel kernel;
    // kernel arguments
    cl_mem src_buffer;
    cl_mem dst_buffer;
    cl_int image_width;
    cl_int image_height;
    cl_int image_stride;
    cl_mem M_inv_buffer;

    float theta;
    WorkerThread *workerThread;

    int initialize_opencl(void);

    void build_info(cl_program program, cl_device_id device);
    void rotSliderInc(int delta);

public slots:
    void timer_func(void);
    void rotationReturnPressed(void);
    void rotSliderPressed(void);
    void rotSliderValueChanged(int value);
    void rotpClicked(bool checked);
    void rotmClicked(bool checked);


protected:

};

#endif // MAINWINDOW_H
