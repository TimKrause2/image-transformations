#define GLM_ENABLE_EXPERIMENTAL
#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QtMath>
#include <QTimer>
#include <QMatrix>
#include <QDoubleValidator>
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
    QValidator *validator = new QDoubleValidator(0.0,360.0,10,this);
    ui->rotLineEdit->setValidator(validator);
    ui->rotLineEdit->setText("0.0");
    ui->rotpPushButton->setIcon(QIcon(":/images/resources/rotp.png"));
    ui->rotmPushButton->setIcon(QIcon(":/images/resources/rotm.png"));

    connect(ui->rotLineEdit,
            &QLineEdit::returnPressed,
            this,
            &MainWindow::rotationReturnPressed);
    connect(ui->rotSlider,
            &QSlider::sliderPressed,
            this,
            &MainWindow::rotSliderPressed);
    connect(ui->rotSlider,
            &QSlider::valueChanged,
            this,
            &MainWindow::rotSliderValueChanged);
    connect(ui->rotpPushButton,
            &QAbstractButton::clicked,
            this,
            &MainWindow::rotpClicked);
    connect(ui->rotmPushButton,
            &QAbstractButton::clicked,
            this,
            &MainWindow::rotmClicked);

    theta = 0.0;
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

    qDebug("image format:%d image argb format:%d",image_format,image_argb_format);

    clFile.open(QIODevice::ReadOnly);
    clData = new char[clFile.size()+1];
    clFile.read(clData,clFile.size());
    clData[clFile.size()]=0;

    workerThread = NULL;

    if(initialize_opencl()){
        timer = new QTimer(this);

        connect(timer, &QTimer::timeout, this, &MainWindow::timer_func);

        timer->start(1000/60);
    }

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
    size_t srcLengths[]={(size_t)clFile.size()};
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
    kernel = clCreateKernel(program,"affine_transform_aa_exp",&r);
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

void MainWindow::rotSliderInc(int delta)
{
    int value = ui->rotSlider->value();
    value += delta;
    if(value>=0 && value<=360){
        ui->rotSlider->setValue(value);
    }
}

void MainWindow::timer_func(void)
{
    // test if the worker thread is finished
    if(workerThread){
        if(!workerThread->isFinished()){
            return;
        }
        QImage dstImage((uchar*)dstImage_data,image_width,
                        image_height,image_stride*sizeof(QRgb),
                        QImage::Format_ARGB32);

        dstPixmap = QPixmap::fromImage(dstImage);
        dstPixmap.save("/tmp/images.png");

        rightPixmapItem->setPixmap(dstPixmap);


    }

    //glm::vec2 A((float)srcImage_argb.width()/2,-(float)srcImage_argb.height()/2);
    glm::vec2 A((float)(srcImage_argb.width()-1)/2,
                -(float)(srcImage_argb.height()-1)/2);
    //glm::vec2 A(0.5f,-0.5f);
    glm::mat3 M = glm::translate(glm::mat3(1.0f),A);
    M = glm::rotate(M, theta);
    //M = glm::shearX(M,-tan(10.0f/180.0f*(float)M_PI));
    M = glm::scale(M,glm::vec2(-1.0f,1.0f));
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

    size_t global_work_size[]={(size_t)image_height};
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

}

void MainWindow::rotationReturnPressed(void)
{
    QString text = ui->rotLineEdit->text();
    float angle = (float)text.toDouble();
    theta = angle*M_PI/180.0;

}

void MainWindow::rotSliderPressed()
{
    int value = ui->rotSlider->value();
    theta = (float)(value)*M_PI/180.0f;
    QString rot_str= QString("%1.0").arg(value);
    ui->rotLineEdit->setText(rot_str);
}

void MainWindow::rotSliderValueChanged(int value)
{
    theta = (float)(value)*M_PI/180.0f;
    QString rot_str= QString("%1.0").arg(value);
    ui->rotLineEdit->setText(rot_str);
}

void MainWindow::rotpClicked(bool checked)
{
    rotSliderInc(1);
}

void MainWindow::rotmClicked(bool checked)
{
    rotSliderInc(-1);
}
