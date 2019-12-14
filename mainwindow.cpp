#define GLM_ENABLE_EXPERIMENTAL
#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QtMath>
#include <QTimer>
#include <QMatrix>
#include <glm/glm.hpp>
#include <glm/gtx/matrix_transform_2d.hpp>
#include <glm/gtc/matrix_access.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <math.h>

#define RPS_MIN 0.00001
#define RPS_MAX 10.0

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    srcPixmap(":/images/resources/hello.png"),
    clFile(":/images/resources/kernel.c", parent)
{
    ui->setupUi(this);

    theta = 0.0;
    rps = 0.0001;
    dtheta = rps*2*M_PI/60;
    setFocusPolicy(Qt::ClickFocus);

    scene = new QGraphicsScene(0,0,srcPixmap.width()*2,srcPixmap.height(), this);
    ui->graphicsView->setScene(scene);
    scene->setBackgroundBrush(QBrush(QPixmap(":/images/resources/background.png")));

    leftPixmapItem = scene->addPixmap(srcPixmap);
    leftPixmapItem->setOffset(0.0,0.0);
    rightPixmapItem = scene->addPixmap(srcPixmap);
    rightPixmapItem->setOffset((qreal)srcPixmap.width(),0.0);

    srcImage = srcPixmap.toImage();
    QImage::Format image_format = srcImage.format();
    srcImage_argb = srcImage.convertToFormat(QImage::Format_ARGB32);
    image_width = srcImage_argb.width();
    image_height = srcImage_argb.height();
    image_stride = srcImage_argb.bytesPerLine()/sizeof(QRgb);
    qDebug("width:%d height:%d stride:%d\n",image_width,image_height,image_stride);
    dstImage_data = new char[srcImage_argb.byteCount()];


    QImage::Format image_argb_format = srcImage_argb.format();
    QString msg = QString("format:%1\n argb_format:%2").arg(image_format).arg(image_argb_format);
    QGraphicsSimpleTextItem *msgTextItem = scene->addSimpleText(msg);
    msgTextItem->setBrush(QColor(128,128,128,255));

    clFile.open(QIODevice::ReadOnly);
    clData = new char[clFile.size()+1];
    clFile.read(clData,clFile.size());
    clData[clFile.size()]=0;

    QString clString(clData);
    QGraphicsSimpleTextItem *clTextItem = scene->addSimpleText(clString);
    clTextItem->setBrush(QColor(255,255,255,255));
    clTextItem->setPos(64,64);
    clTextItem->setFont(QFont("Ubuntu mono",8));
    QRectF rect = clTextItem->boundingRect();
    qreal x1,y1,x2,y2;
    rect.getCoords(&x1,&y1,&x2,&y2);
    scene->setSceneRect(0,0,srcPixmap.width()*2,y2+64);

    if(initialize_opencl()){
        timer = new QTimer(this);

        connect(timer, &QTimer::timeout, this, &MainWindow::timer_func);

        timer->start(1000/60);
    }

    workerThread = NULL;
}

MainWindow::~MainWindow()
{
    delete ui;
}

int MainWindow::initialize_opencl(void)
{
    clGetPlatformIDs(0, NULL, &num_platforms);
    if(num_platforms==0){
        qDebug("No OpenCL platforms found");
        return 0;
    }
    platforms = new cl_platform_id[num_platforms];
    clGetPlatformIDs(num_platforms,platforms,NULL);
    clGetDeviceIDs(platforms[0],CL_DEVICE_TYPE_GPU,0,NULL,&num_devices);
    if(num_devices==0){
        qDebug("No GPUs found for platform.");
        return 0;
    }
    devices = new cl_device_id[num_devices];
    clGetDeviceIDs(platforms[0],CL_DEVICE_TYPE_GPU,num_devices,devices,NULL);
    cl_context_properties context_properties[]={
        CL_CONTEXT_PLATFORM,
        (cl_context_properties)platforms[0],
        0
    };
    cl_int r;
    context = clCreateContext(
                context_properties,
                1,devices,NULL,
                NULL,&r);
    if(r!=CL_SUCCESS){
        qDebug("clCreateContext failed");
        return 0;
    }
    qDebug("Context created.");
    queue = clCreateCommandQueue(
                context,
                devices[0],
            0,
            NULL);
    size_t srcLengths[]={clFile.size()};
    program = clCreateProgramWithSource(
                context,
                1,
                (const char**)&clData,
                srcLengths,
                &r);
    if(r!=CL_SUCCESS){
        qDebug("clCreateProgramWithSource failed");
        return 0;
    }
    r = clBuildProgram(program, 1, devices, NULL, NULL, NULL);
    if(r!=CL_SUCCESS){
        qDebug("clBuildProgram failed");
        build_info(program,devices[0]);
        return 0;
    }
    kernel = clCreateKernel(program,"affine_transform_aa",&r);
    if(r!=CL_SUCCESS){
        qDebug("clCreateKernel failed %d",r);
        return 0;
    }
    src_buffer = clCreateBuffer(
                context,
                CL_MEM_READ_ONLY,
                srcImage_argb.byteCount(),
                NULL,
                NULL);
    dst_buffer = clCreateBuffer(
                context,
                CL_MEM_WRITE_ONLY,
                srcImage_argb.byteCount(),
                NULL,
                NULL);
    M_inv_buffer = clCreateBuffer(
                context,
                CL_MEM_READ_ONLY,
                sizeof(glm::mat3),
                NULL,
                NULL);

    return 1;
}






void MainWindow::build_info(cl_program program, cl_device_id device){
    cl_int r;
    cl_build_status status;
    r = clGetProgramBuildInfo(program, device,
                              CL_PROGRAM_BUILD_STATUS,
                           sizeof(cl_build_status),
                              &status, NULL);
    if(r!=CL_SUCCESS){
        qDebug("clGetProgramBuildInfo error:%d\n",r);
        return;
    }
    qDebug("Build status:");
    switch(status){
        case CL_BUILD_NONE:
            qDebug("CL_BUILD_NONE");
            break;
        case CL_BUILD_ERROR:
            qDebug("CL_BUILD_ERROR");
            break;
        case CL_BUILD_SUCCESS:
            qDebug("CL_BUILD_SUCCESS");
            break;
        case CL_BUILD_IN_PROGRESS:
            qDebug("CL_BUILD_IN_PROGRESS");
            break;
    }
    qDebug("\n");
    if(status==CL_BUILD_ERROR){
        size_t log_size;
        r = clGetProgramBuildInfo(program,device,
                                  CL_PROGRAM_BUILD_LOG,
                            0,NULL,&log_size);
        if(r!=CL_SUCCESS){
            qDebug("couldn't get log size. error:%d\n",r);
            return;
        }
        char *build_log = new char[log_size];
        r = clGetProgramBuildInfo(program,device,
                                  CL_PROGRAM_BUILD_LOG,
                            log_size, build_log, NULL);
        if(r!=CL_SUCCESS){
            qDebug("couldn't get build log. error:%d\n",r);
            return;
        }
        qDebug("Build log:%s\n",build_log);
        delete build_log;
    }
}

void MainWindow::timer_func(void)
{
    // test if the worker thread is finished
    if(workerThread){
        if(!workerThread->isFinished()){
            return;
        }
    }

    //glm::vec2 A((float)srcImage_argb.width()/2,-(float)srcImage_argb.height()/2);
    glm::vec2 A((float)(srcImage_argb.width()-1)/2,
                -(float)(srcImage_argb.height()-1)/2);
    //glm::vec2 A(0.5f,-0.5f);
    glm::mat3 M = glm::translate(glm::mat3(1.0f),A);
    M = glm::rotate(M, theta);
    //glm::vec2 scales(1.0f,1.0f);
    //M = glm::scale(M,scales);
    M = glm::translate(M, -A);
    glm::mat3 M_inv = glm::inverse(M);
    //qDebug("a11:%.9f a12:%.9f a21:%.9f a22:%.9f",M_inv[0][0],M_inv[1][0],M_inv[0][1],M_inv[1][1]);
    cl_int r = clEnqueueWriteBuffer(
                queue,
                src_buffer,
                CL_FALSE,
                0,
                srcImage_argb.byteCount(),
                srcImage_argb.bits(),
                0,
                NULL,
                NULL);
    if(r!=CL_SUCCESS){
        qDebug("Write Buffer src buffer failed");
        return;
    }
    r = clEnqueueWriteBuffer(
                queue,
                M_inv_buffer,
                CL_FALSE,
                0,
                sizeof(glm::mat3),
                glm::value_ptr(M_inv),
                0,NULL,NULL);
    if(r!=CL_SUCCESS){
        qDebug("Write Buffer M_inv_buffer failed");
        return;
    }
    clSetKernelArg(kernel,0,sizeof(cl_mem),&src_buffer);
    clSetKernelArg(kernel,1,sizeof(cl_mem),&dst_buffer);
    clSetKernelArg(kernel,2,sizeof(cl_int),&image_stride);
    clSetKernelArg(kernel,3,sizeof(cl_int),&image_width);
    clSetKernelArg(kernel,4,sizeof(cl_int),&image_height);
    clSetKernelArg(kernel,5,sizeof(cl_mem),&M_inv_buffer);

    size_t global_work_size[]={image_height};
    r = clEnqueueNDRangeKernel(
                queue,
                kernel,
                1,
                NULL,
                global_work_size,
                NULL,
                0,NULL,NULL);
    if(r!=CL_SUCCESS){
        qDebug("clEnqueueNDRangeKernel failed r:%d",r);
        return;
    }

    r = clEnqueueReadBuffer(
                queue,
                dst_buffer,
                CL_FALSE,
                0,
                srcImage_argb.byteCount(),
                dstImage_data,
                0,NULL,NULL);
    if(r!=CL_SUCCESS){
        qDebug("clEnqueueReadBuffer error:%d",r);
        return;
    }

    if(workerThread){
        delete workerThread;
    }
    workerThread = new WorkerThread(queue);
    workerThread->start();





    //r = clFinish(queue);
    //if(r!=CL_SUCCESS){
    //    qDebug("clFinish error:%d",r);
    //    return;
    //}

    QImage dstImage((uchar*)dstImage_data,image_width,
                    image_height,image_stride*sizeof(QRgb),
                    QImage::Format_ARGB32);

    dstPixmap = QPixmap::fromImage(dstImage);

    rightPixmapItem->setPixmap(dstPixmap);

    theta+=dtheta;
    if(theta>2*M_PI)theta-=2*M_PI;
}

void MainWindow::update_rps(float m){
    rps*=m;
    if(rps<RPS_MIN)rps = RPS_MIN;
    if(rps>RPS_MAX)rps = RPS_MAX;
    dtheta = rps*2*M_PI/60;
}


bool MainWindow::event(QEvent *event)
{
    if(event->type() == QEvent::KeyPress){
        QKeyEvent *keyEvent = static_cast<QKeyEvent*>(event);
        int key = keyEvent->key();
        qDebug("keyPressEvent key:0x%X",key);
        if(key==Qt::Key_Up){
            if(keyEvent->modifiers()&Qt::ControlModifier){
                update_rps(1.259);
            }else{
                update_rps(10.0);
            }
            return true;
        }
        else if(key==Qt::Key_Down){
            if(keyEvent->modifiers()&Qt::ControlModifier){
                update_rps(0.794);
            }else{
                update_rps(0.1);
            }
            return true;
        }else{
            return false;
        }
    }
    return QMainWindow::event(event);
}
