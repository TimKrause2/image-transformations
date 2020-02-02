# image-transformations
OpenCL image transformations with anti-aliasing
This is a Qt Creator project. Load with Qt Creator and build for your machine.
On Ubuntu the OpenCL development files and glm can be installed with:

$ sudo apt install opencl-headers
$ sudo apt install libglm-dev
$ sudo apt install ocl-icd-opencl-dev

And then you'll also need an ICD(Installable Client Driver). This depends on
your hardware. You can try https://packages.ubuntu.com/search?keywords=opencl-icd&searchon=names&suite=all&section=all
and check for a driver that matches your hardware. You can try the mesa opencl ICD first.

$ sudo apt install mesa-opencl-icd

As you can see you may have to add a repository to apt. This can
be done with:

$ sudo apt-add-repository universe

or if you want to add the multiverse repository

$ sudo apt-add-repository multiverse

If you want to check your installation use clinfo. It can be installed with:

$ sudo apt install clinfo

To use it just type:

$ clinfo -l

And that will show you all of the platforms and and all of the devices for each platform. You can have
multiple platforms and each platform can have multiple devices. If you want to see all of the information
available for all of the platforms and devices just type:

$ clinfo | less

and then use the page-up, page-down, arrow keys and more to view the result. Press 'q' to exit less.

