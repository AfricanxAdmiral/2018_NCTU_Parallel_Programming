#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include <fstream>
#include <iostream>
#include <vector>

#define CL_USE_DEPRECATED_OPENCL_1_2_APIS
#ifdef __APPLE__   // for runing the code on MacOS
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

int main(int argc, char const *argv[])
{
        // get platform number
        cl_int err;
        cl_uint num;
        err = clGetPlatformIDs(0, 0, &num);
        if (err != CL_SUCCESS)
        {
                std::cerr << "Unable to get platforms\n";
                return 0;
        }

        // accutally get the platform
        std::vector<cl_platform_id> platforms(num);
        err = clGetPlatformIDs(num, &platforms[0], &num);
        if (err != CL_SUCCESS)
        {
                std::cerr << "Unable to get platform ID\n";
                return 0;
        }

        // print out platform info (test)
        // for (unsigned int i = 0; i < num; ++i)
        // {
        //      char pbuff[100];
        //      err = clGetPlatformInfo(platforms[i], CL_PLATFORM_VENDOR, sizeof(pbuff), pbuff, NULL);
        //      std::cout << pbuff << std::endl;
        // }

        // build a OpenCl context
        cl_context_properties prop[] = { CL_CONTEXT_PLATFORM, reinterpret_cast<cl_context_properties>(platforms[0]), 0 };
        cl_context context = clCreateContextFromType(prop, CL_DEVICE_TYPE_DEFAULT, NULL, NULL, NULL);
        if (!context)
        {
                std::cerr << "Can't create OpenCL context\n";
                return 0;
        }

        // get info of context
        size_t cb;
        clGetContextInfo(context, CL_CONTEXT_DEVICES, 0, NULL, &cb);
        std::vector<cl_device_id> devices(cb / sizeof(cl_device_id));
        clGetContextInfo(context, CL_CONTEXT_DEVICES, cb, &devices[0], 0);

        // print out device info (test)
        // clGetDeviceInfo(devices[0], CL_DEVICE_VERSION, 0, NULL, &cb);
        // for (unsigned int i = 0; i < devices.size(); i++)
        // {
        //      char dbuff[100];
        //      clGetDeviceInfo(devices[i], CL_DEVICE_VERSION, sizeof(dbuff), dbuff, 0);
        //      std::cout << "Device : " << dbuff << std::endl;
        // }

        // create Command Queue
        cl_command_queue queue = clCreateCommandQueue(context, devices[0], 0, 0);
        if (!queue)
        {
                std::cerr << "Can't Create command queue\n";
                clReleaseContext(context);
                return 0;
        }

        unsigned int* histogram_results = (unsigned int*)malloc(256*3*sizeof(unsigned int));
        unsigned int i=0, a, input_size;
        std::fstream inFile("input", std::ios_base::in);
        std::ofstream outFile("0516069.out", std::ios_base::out);

        inFile >> input_size;
        unsigned int *image = new unsigned int[input_size];
        while( inFile >> a ) {
                image[i++] = a;
        }

        // create memory object on device
        cl_mem cl_img = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_uint)*input_size, image, NULL);
        cl_mem cl_result = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(cl_uint) * 256 * 3, NULL, NULL);
        if (!cl_img || !cl_result)
        {
                std::cerr << "Can't create OpenCL buffer\n";
                clReleaseMemObject(cl_result);
                clReleaseMemObject(cl_img);
                clReleaseCommandQueue(queue);
                clReleaseContext(context);
                return 0;
        }

        // build program
        std::ifstream in("histogram.cl", std::ios_base::binary);
        if(!in.good()) return 0;

        // get file length
        in.seekg(0, std::ios_base::end);
        size_t length = in.tellg();
        in.seekg(0, std::ios_base::beg);

        // read program source
        std::vector<char> data(length + 1);
        in.read(&data[0], length);
        data[length] = 0;

        // create and build program
        const char* source = &data[0];
        cl_program program = clCreateProgramWithSource(context, 1, &source, 0, 0);
        if(program == 0)
        {
                std::cerr << "Failed create program object\n";
                clReleaseMemObject(cl_result);
                clReleaseMemObject(cl_img);
                clReleaseCommandQueue(queue);
                clReleaseContext(context);
                return 0;
        }
        if(clBuildProgram(program, 0, 0, 0, 0, 0) != CL_SUCCESS)
        {
                std::cerr << "Failed building program\n";
                clReleaseMemObject(cl_result);
                clReleaseMemObject(cl_img);
                clReleaseCommandQueue(queue);
                clReleaseContext(context);
                return 0;
        }
        // std::cout << "Finish build program\n";

        // build kernel
        cl_kernel adder = clCreateKernel(program, "adder", 0);
        if(!adder) {
                std::cerr << "Can't load kernel\n";
                clReleaseProgram(program);
                clReleaseMemObject(cl_result);
                clReleaseMemObject(cl_img);
                clReleaseCommandQueue(queue);
                clReleaseContext(context);
                return 0;
        }
        // std::cout << "Finish build kernel\n";

        // set the argument of kernel function
        clSetKernelArg(adder, 0, sizeof(cl_mem), &cl_img);
        clSetKernelArg(adder, 1, sizeof(cl_mem), &cl_result);

        size_t work_size = input_size;
        err = clEnqueueNDRangeKernel(queue, adder, 1, 0, &work_size, 0, 0, 0, 0);

        if(err == CL_SUCCESS)
        {
                err = clEnqueueReadBuffer(queue, cl_result, CL_TRUE, 0, sizeof(unsigned int) * 256 * 3, histogram_results, 0, 0, 0);
                //std::cout << "Program executed : " << err << std::endl;
        }

        if(err != CL_SUCCESS)
        {
                std::cerr << "Failed executing program\n";
                clReleaseProgram(program);
                clReleaseMemObject(cl_result);
                clReleaseMemObject(cl_img);
                clReleaseCommandQueue(queue);
                clReleaseContext(context);
                return 0;
        }

        //histogram_results = histogram(image, input_size);
        for(unsigned int i = 0; i < 256 * 3; ++i) {
                if (i % 256 == 0 && i != 0)
                        outFile << std::endl;
                outFile << histogram_results[i]<< ' ';
        }

        inFile.close();
        outFile.close();

        // release resource in OpenCL
        clReleaseMemObject(cl_result);
        clReleaseMemObject(cl_img);
        clReleaseKernel(adder);
        clReleaseProgram(program);
        clReleaseCommandQueue(queue);
        clReleaseContext(context);

        return 0;
}