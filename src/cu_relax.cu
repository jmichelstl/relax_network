/*
This code performs an energy relaxation of a network using the FIRE algorithm.
*/

#include <iostream>
#include <stdio.h>
#include <fstream>
#include "cu_network_utils.h"
#include <map>
#include <vector>
#include <sstream>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <unistd.h>
#include <cuda_runtime.h>
#include "cublas_v2.h"
#include "math.h"
#include "float.h"

#define TENSILE 0
#define SHEAR 1
#define FX_MASK 0x1
#define FY_MASK 0x2
#define FX_FY 0x3

using std::map;
using std::ostringstream;
using std::cerr;
using std::make_pair;

enum incrform {constant, logarithmic};

struct EdgeDatum{

    EdgeDatum(int ival, double len) : index(ival), l0(len){}

    int index;
    double l0;
};

//Given a planar line graph specified as pairs of floating-point x and y
//values, create a list of unique points, assign each point an index, create
//a list of index pairs for points that share bonds, and list the relaxed 
//distances between pairs of points that share bonds.
void read_edges(ifstream &datfile, int *num_dof, int *npairs, double **pos, int **pair_list, unsigned char** fcodes, vector<int> &top, double **lengths, bool topSlide){

    map<int, vector<EdgeDatum>> neighbor_map;
    vector<double> point_list;
    map<Point, int> pmap;
    string nextline;
    Point p1, p2;
    int pindex = 0, index1, index2, mindex, maxdex, point_count = 0, iter;
    vector<double> edge_data;
    double length, miny, maxy;

    miny = FLT_MAX;
    maxy = FLT_MIN;
    *npairs = 0;

    while(! datfile.eof()){
        getline(datfile, nextline);
        edge_data = parse_doubles(split(nextline, ' '));

        if(edge_data.size() >= 2){
            if(point_count == 0){
                p1 = Point(edge_data[0], edge_data[1]);
                point_count ++;
            }

            else{
                p2 = Point(edge_data[0], edge_data[1]);
                point_count ++;
            }

            if(point_count == 2){
                point_count = 0;
                *npairs = *npairs + 1;

                if(pmap.find(p1) == pmap.end()){
                    pmap.insert(make_pair(p1, pindex));
                    point_list.push_back(p1.x);
                    point_list.push_back(p1.y);
                    neighbor_map.insert(make_pair(pindex, vector<EdgeDatum>()));
                    pindex++;
                    miny = p1.y < miny - FLOAT_TOL ? p1.y : miny;
                    maxy = p1.y > maxy + FLOAT_TOL ? p1.y : maxy;
                }
                if(pmap.find(p2) == pmap.end()){
                    pmap.insert(make_pair(p2, pindex));
                    point_list.push_back(p2.x);
                    point_list.push_back(p2.y);
                    neighbor_map.insert(make_pair(pindex, vector<EdgeDatum>()));
                    pindex++;
                    miny = p2.y < miny - FLOAT_TOL ? p2.y : miny;
                    maxy = p2.y > maxy + FLOAT_TOL ? p2.y : maxy;
                }

                index1 = pmap[p1];
                index2 = pmap[p2];
                length = sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y));
                mindex = index1 < index2 ? index1 : index2;
                maxdex = index1 == mindex ? index2 : index1;

                neighbor_map[mindex].push_back(EdgeDatum(maxdex, length));
            }
        }
    }

    //Establish degrees of freedom and the number of bonds, and allocate
    //arrays used in simulating network mechanics
    *num_dof = point_list.size();
    *pos = (double *) calloc(*num_dof, sizeof(double));
    *fcodes = (unsigned char*) calloc(*num_dof / 2, sizeof(char));
    for(iter = 0; iter < *num_dof; iter++){
        (*pos)[iter] = point_list[iter];
        //Determine what mask to use to indicate whether a point may relax in
        //the x and y directions
        if(iter % 2 == 1){
            if(point_list[iter] <= miny + FLOAT_TOL){
                if(topSlide) (*fcodes)[iter / 2] = FX_MASK;
                else (*fcodes)[iter / 2] = 0;
            }
            else if(point_list[iter] >= maxy - FLOAT_TOL){
                top.push_back(iter / 2);
                if(topSlide) (*fcodes)[iter / 2] = FX_MASK;
                else (*fcodes)[iter / 2] = 0;
            }
            else (*fcodes)[iter / 2] = FX_FY;
        }
    }

    //Create a list of index pairs. If the topSlide parameter is set, then
    //allow the x coordinates of points along the top and bottom to relax
    *pair_list = (int *) calloc(2 * *npairs, sizeof(int));
    *lengths = (double *) calloc(*npairs, sizeof(double));
    iter = 0;
    for(auto map_iter = neighbor_map.begin(); map_iter != neighbor_map.end(); map_iter++){
        index1 = map_iter->first;

        for(EdgeDatum edat : map_iter->second){
            index2 = edat.index;
            (*pair_list)[iter * 2] = index1;
            (*pair_list)[iter * 2 + 1] = index2;
            (*lengths)[iter] = edat.l0;
            iter ++;
        }
    }
}

int get_dec_digits(int value){
    int num_digits = 0;

    do{
        num_digits++;
        value /= 10;
    }while(value > 0);

    return num_digits;
}

string report_name(string base, int num_digits, int count){
    int padding, iter;
    ostringstream oss;

    padding = num_digits - get_dec_digits(count);
    oss << base << "_";
    for(iter = 0; iter < padding; iter++) oss << 0;
    oss << count << ".dat";

    return oss.str();
}

void get_report_params(string& base, int& rep_freq){

    string response;
    int num_read;

    do{
        cout << "Enter base name: ";
        getline(cin, response);
    }while(response.compare("") == 0);

    base = response;

    do{
        cout << "Enter reporting frequency: ";
        getline(cin, response);
        num_read = sscanf(response.c_str(), "%d", &rep_freq);
    }while(num_read < 1 || rep_freq < 1);
}

//I found this hack on the web for atomic addition of double precision floats.
//This is needed because, with compute capability less than 6, there is no
//standard support for such an operation
__device__ double atomicSum(double* address, double val){
    unsigned long long int* address_as_ull =
                          (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
            old = atomicCAS(address_as_ull, assumed,__double_as_longlong(val +
                               __longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
}

__global__ void cu_get_pe(int npairs, double *pos, int *pair_list, double *lengths, double *pe){
    int pair_index, index1, index2;
    double x1, x2, y1, y2, l0, dist;

    pair_index = blockIdx.x * blockDim.x * blockDim.y;
    pair_index += threadIdx.y * blockDim.x + threadIdx.x;
    if(pair_index >= npairs) return;

    index1 = pair_list[pair_index * 2];
    index2 = pair_list[pair_index * 2 + 1];

    x1 = pos[index1*2];
    y1 = pos[index1*2 + 1];
    x2 = pos[2*index2];
    y2 = pos[2*index2 + 1];
    l0 = lengths[pair_index];
    dist = hypot(x1 - x2, y1 - y2);
    pe[pair_index] = .5 * (dist - l0) * (dist - l0);
}

__global__ void get_forces(int npairs, double *pos, unsigned char *fcodes, int *pair_list, double *lengths, double *forces){
    int pair_index, index1, index2;
    double x1, x2, y1, y2, l0, inv_dist, fmult, fx, fy;

    pair_index = blockIdx.x * blockDim.x * blockDim.y;
    pair_index += threadIdx.y * blockDim.x + threadIdx.x;
    if(pair_index >= npairs)  return;

    index1 = pair_list[pair_index * 2];
    index2 = pair_list[pair_index * 2 + 1];

    x1 = pos[index1*2];
    y1 = pos[index1*2 + 1];
    x2 = pos[2*index2];
    y2 = pos[2*index2 + 1];
    l0 = lengths[pair_index];
    inv_dist = rhypot(x1 - x2, y1 - y2);
    fmult = 1 - l0 * inv_dist;
    fx = fmult * (x2 - x1);
    fy = fmult * (y2 - y1);

    //If a node is not at the bottom or top of the network, relax
    //the x and y coordinates. If a node is on the top or bottom, determine
    //whether the node may relax in the x direction, and if so, relax in the x
    //direction.
    atomicAdd(forces + index1*2, fx*(fcodes[index1] & FX_MASK));
    atomicAdd(forces + index1*2 + 1, fy*((fcodes[index1] & FY_MASK) / 2));

    atomicAdd(forces + index2*2, -fx*(fcodes[index2] & FX_MASK));
    atomicAdd(forces + index2*2 + 1, -fy*((fcodes[index2] & FY_MASK) / 2));
}

double  get_pe(int npairs, int *pair_list, double *pos, double *lengths){
    double pe = 0;
    int iter, index1, index2;
    double distance, diff, x1, y1, x2, y2;

    for(iter = 0; iter < npairs; iter ++){
        index1 = pair_list[iter * 2];
        index2 = pair_list[iter * 2 + 1];
        x1 = pos[2*index1];
        y1 = pos[2*index1 + 1];
        x2 = pos[2*index2];
        y2 = pos[2*index2 + 1];
        distance = hypot(x1 - x2, y1 - y2);
        diff = distance - lengths[iter];
        pe += .5 * diff * diff;
    }

    return pe;
}

void report_deformed(int npairs, int *pair_list, double *pos, string filename){
    int iter, index1, index2;
    FILE *datfile;

    datfile = fopen(filename.c_str(), "w");

    for(iter = 0; iter < npairs; iter ++){
        index1 = pair_list[iter * 2];
        index2 = pair_list[iter * 2 + 1];
        fprintf(datfile, "%2.12lf %2.12lf\n", pos[index1*2], pos[2*index1+1]);
        fprintf(datfile, "%2.12lf %2.12lf\n\n",pos[index2*2],pos[2*index2+1]);
    }

    fclose(datfile);
}

/*
void simple_md(int dof, double *pos, double *vel, int *pair_list, double mass, double dt, int num_steps){

    double *forces;
    int step_count, index;
    double dt_sq = dt*dt;

    for(step_count = 0; step_count < num_steps; step_count ++){
        get_forces(neighbor_map, pos, forces, bottom, top, false);
        for(index = 0; index < pos.size(); index++){
            if(bottom.find(index/2) == bottom.end() && top.find(index/2) == top.end()){
                pos[index] += dt*vel[index] + .5 * forces[index]/mass * dt_sq;
                vel[index] += .5 * dt * forces[index] / mass;
            }
        }

        get_forces(neighbor_map, pos, forces, bottom, top, false);
        for(index = 0; index < pos.size(); index++){
            if(bottom.find(index/2) == bottom.end() && top.find(index/2) == top.end()){
                vel[index] += .5 * dt * forces[index] / mass;
            }
        }

    }
}
*/

void check(cudaError_t error){
    if(error != cudaSuccess){
        printf("code: %d, reason: %s\n", error, cudaGetErrorString(error));
        exit(1);
    }
}

int run_fire_md(int num_dof, int npairs, double *pos, double *vel, int *pair_list, double *lengths, unsigned char* fcodes, double mass, double fcut, int max_steps, int nmin, double fire_params[5], bool& flag, dim3 grid, dim3 block, FILE *report_file, int rep_freq, double &frms){

    double *d_forces, *d_pos, *d_vel, *d_lengths, *d_pe;
    int *d_pair_list;
    unsigned char *d_fcodes;
    double alpha, finc, fdec, alpha_start, falpha, fire_dt_max, sqrt_dof, pe;
    double power, dt, dt_mult, dt_sq_mult, vmag, fmag, a_vm_o_fm;
    double neg_alpha; 
    int step_count = 0, since_leq_0 = 0, device  = 0;
    string base, full_name;
    cublasHandle_t handle;

    cudaSetDevice(device);

    //Initialize device vectors
    check(cudaMalloc((void **) &d_forces, sizeof(double) * num_dof));
    cudaDeviceSynchronize();

    check(cudaMalloc((void **) &d_pos, sizeof(double) * num_dof));
    cudaDeviceSynchronize();
    cudaMemcpy(d_pos, pos, sizeof(double) * num_dof, cudaMemcpyHostToDevice);

    check(cudaMalloc((void **) &d_fcodes, sizeof(char) * num_dof));
    cudaDeviceSynchronize();
    cudaMemcpy(d_fcodes, fcodes, sizeof(char)*num_dof, cudaMemcpyHostToDevice);

    check(cudaMalloc((void **) &d_vel, sizeof(double) * num_dof));
    cudaDeviceSynchronize();
    cudaMemcpy(d_vel, vel, sizeof(double) * num_dof, cudaMemcpyHostToDevice);

    check(cudaMalloc((void **) &d_pair_list, sizeof(int) * 2 * npairs));
    cudaDeviceSynchronize();
    cudaMemcpy(d_pair_list, pair_list, 2*sizeof(int)*npairs, cudaMemcpyHostToDevice);

    check(cudaMalloc((void **) &d_lengths, sizeof(double) * npairs));
    cudaDeviceSynchronize();
    cudaMemcpy(d_lengths, lengths, sizeof(double)*npairs, cudaMemcpyHostToDevice);

    check(cudaMalloc((void **) &d_pe, sizeof(double) * npairs));
    cudaDeviceSynchronize();

    //Initialize handle for CUDA-accelerated blas calculations
    cublasCreate(&handle);

    //Unpack parameters of FIRE minimization scheme and initialize values
    alpha_start = fire_params[0];
    falpha = fire_params[1];
    fire_dt_max = fire_params[2];
    finc = fire_params[3];
    fdec = fire_params[4];

    alpha = alpha_start;
    dt = fire_dt_max;
    dt_sq_mult = .5 * dt*dt / mass;
    dt_mult = .5 * dt / mass;
    sqrt_dof = sqrt(num_dof);

    //Find the forces at the outset
    cudaMemset(d_forces, 0, sizeof(double) * num_dof);
    cudaDeviceSynchronize();

    get_forces<<<grid, block>>>(npairs, d_pos, d_fcodes, d_pair_list, d_lengths, d_forces);
    cudaDeviceSynchronize();
    //cublasDnrm2(handle, num_dof, d_vel, 1, &vmag);

    cublasDnrm2(handle, num_dof, d_forces, 1, &fmag);
    cudaDeviceSynchronize();
    cout << "Starting RMS Force: " << fmag / sqrt_dof << "\n";
 
    //Perform molecular dynamics steps using velocity verlet method until the
    //kinetic energy cutoff is reached, or the maximum number of steps have
    //taken place
    while(step_count < max_steps){
        step_count ++;        

        //Update positions and velocities
        cublasDaxpy(handle, num_dof, &dt, d_vel, 1, d_pos, 1);
        cudaDeviceSynchronize();
        cublasDaxpy(handle, num_dof, &dt_sq_mult, d_forces, 1, d_pos, 1);
        cudaDeviceSynchronize();
        cublasDaxpy(handle, num_dof, &dt_mult, d_forces, 1, d_vel, 1);
        cudaDeviceSynchronize();

        //Calculate forces
        cudaMemset(d_forces, 0, sizeof(double) * num_dof);
        cudaDeviceSynchronize();
        get_forces<<<grid, block>>>(npairs, d_pos, d_fcodes, d_pair_list, d_lengths, d_forces);
        cudaDeviceSynchronize();

        //Update velocities and calculate power
        cublasDaxpy(handle, num_dof, &dt_mult, d_forces, 1, d_vel, 1);
        cudaDeviceSynchronize();
        cublasDdot(handle, num_dof, d_vel, 1, d_forces, 1, &power);
        cudaDeviceSynchronize();

        //Adjust velocities according to FIRE algorithm
        cublasDnrm2(handle, num_dof, d_vel, 1, &vmag);
        cublasDnrm2(handle, num_dof, d_forces, 1, &fmag);
        cudaDeviceSynchronize();
        a_vm_o_fm = alpha * vmag / fmag;

        cublasDaxpy(handle, num_dof, &a_vm_o_fm, d_forces, 1, d_vel, 1);
        cudaDeviceSynchronize();
        cublasDaxpy(handle, num_dof, &neg_alpha, d_vel, 1, d_vel, 1);
        cudaDeviceSynchronize();
        
        //Adjust FIRE parameters according to current power
        if(power > 0){
            since_leq_0 ++;
            if(since_leq_0 > nmin){
                dt = min(dt*finc, fire_dt_max);
                dt_mult = .5 * dt / mass;
                dt_sq_mult = .5*dt*dt / mass;
                alpha *= falpha;
                neg_alpha = -alpha;
            }
        }
        else{
            since_leq_0 = 0;
            dt *= fdec;
            dt_mult = .5 * dt / mass;
            dt_sq_mult = dt * dt * .5 / mass;
            alpha = alpha_start;
            neg_alpha = -alpha;
            cudaMemset(vel, 0, sizeof(double) * num_dof);
            cudaDeviceSynchronize();
        }

        if(report_file != NULL && (step_count % rep_freq == 0)){
            cu_get_pe<<<grid, block>>>(npairs, d_pos, d_pair_list, d_lengths, d_pe);
            cudaDeviceSynchronize();
            cublasDasum(handle, npairs, d_pe, 1, &pe);
            cudaDeviceSynchronize();
            fprintf(report_file, "%d\t%4.12lf\t%4.12lf\n", step_count,pe,fmag);
        }

        //Check for kinetic energy convergence
        if(fmag / sqrt_dof < fcut){
            flag = true;
            //cout << "RMS Force: " << fmag / sqrt_dof << "\n";
	    frms = fmag  / sqrt_dof;
            cudaMemcpy(pos, d_pos, sizeof(double) * num_dof, cudaMemcpyDeviceToHost);
            cudaMemcpy(vel, d_vel, sizeof(double) * num_dof, cudaMemcpyDeviceToHost);
            cudaDeviceSynchronize();
            cudaFree(d_forces);
            cudaFree(d_pos);
            cudaFree(d_fcodes);
            cudaFree(d_vel);
            cudaFree(d_pair_list);
            cudaFree(d_lengths);
            cublasDestroy(handle);
            cudaDeviceReset();
            return step_count;
        }
    }

    flag = false;
    //cout << "RMS Force: " << fmag / sqrt_dof << "\n";
    frms = fmag / sqrt_dof;
    cout << "Ending energy: " << .5 * mass * vmag * vmag + get_pe(npairs, pair_list, pos, lengths) << "\n";

    cudaMemcpy(pos, d_pos, sizeof(double) * num_dof, cudaMemcpyDeviceToHost);
    cudaMemcpy(vel, d_vel, sizeof(double) * num_dof, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    cudaFree(d_forces);
    cudaFree(d_pos);
    cudaFree(d_fcodes);
    cudaFree(d_vel);
    cudaFree(d_pair_list);
    cudaFree(d_lengths);
    cublasDestroy(handle);
    cudaDeviceReset();

    return step_count;
}

void displace(double *pos_vec, vector<int> top, double disp, int dmode){

    if(dmode == TENSILE){
        for(int index : top){
            pos_vec[index*2 + 1] += disp;
        }
    }
    else{
        for(int index : top){
            pos_vec[index*2] += disp;
        }
    }
}

void do_relaxation_run(bool abort, int dmode, bool topSlide, dim3 block, incrform iform){

    ifstream netfile;
    double *pos, *vel, *lengths;
    vector<double> sim_results;
    vector<int> top;
    int *pair_list;
    unsigned char *fcodes;
    double mass, fcut, dmin, dmax, dinc, disp, prev_disp, d_inc_mult, pe, frms;
    int i, num_dof, num_read, max_steps, nmin, npairs, digits;
    double fparams[5];
    string response, report_base, report_file, log_name, i_base, full_name;
    bool converged, batch_report, make_log, i_report;
    int scount, rep_count, numd, rep_freq;
    ofstream log_file;
    ostringstream log_line;
    FILE *ereport = NULL, *stream, *intermediate = NULL;

    //Obtain file with information about the network
    open_dat_file("Enter the file containing network data: ", netfile);
    if(! netfile.is_open()){
        cerr << "The file stream was not open for reading.\n";
        return;
    }

    //Obtain information about vertex locations and connections in the
    //undeformed lattice and initialize vectors used in simulating the lattice
    read_edges(netfile, &num_dof, &npairs, &pos, &pair_list, &fcodes, top, &lengths, topSlide);
    cout << "There are " << npairs << " pairs\n";
    netfile.close();
    vel = (double *) calloc(num_dof, sizeof(double));

    //Finalize the grid layout
    dim3 grid(npairs  / (block.x * block.y) + 1, 1);

    //Get FIRE parameters for adjusting velocity
    do{
        cout << "Enter alpha, falpha, dt, finc, fdec, and nmin: ";
        getline(cin, response);
        num_read = sscanf(response.c_str(),"%lf %lf %lf %lf %lf %d", fparams, fparams+1, fparams+2, fparams+3, fparams+4, &nmin);
    }while(num_read < 6);

    do{
        cout << "Enter particle mass, cutoff force, and maximum md steps: ";
        getline(cin, response);
        num_read = sscanf(response.c_str(), "%lf %lf %d", &mass, &fcut, &max_steps);
    }while(num_read < 3);

    //Get range of displacements
    do{
        cout << "Enter minimum and maximum displacement and increment: ";
        getline(cin, response);
        num_read = sscanf(response.c_str(), "%lf %lf %lf", &dmin, &dmax, &dinc);
    }while(num_read < 3);
    d_inc_mult = pow(10, dinc);

    /*
    //Get preconditioning parameters
    if((precondition = yesno("Do a preconditioning run?"))){
        do{
            cout << "Enter preconditioning time step and number of steps: ";
            getline(cin, response);
            num_read = sscanf(response.c_str(), "%lf %d", &pre_dt, &pre_steps);
        }while(num_read < 2);
    }
    */

    //Give the option to do batch reporting of final networks
    if((batch_report = yesno("Do batch reporting of final networks?"))){
        do{
            cout << "Enter a base name for reporting: ";
            getline(cin, report_base);
        }while(report_base.compare("") == 0);
        numd = get_dec_digits((int) (1 + (dmax - dmin)/dinc));
    }

    //Offer the opportunity to write convergence outcomes to a log file
    if((make_log = yesno("Use a log file?"))){
	cout << "Enter the log file name: ";
        getline(cin, log_name);
    }

    i_report = false;
    if(yesno("Report intermediate energies and forces?")){
        digits = get_dec_digits((int) ((dmax - dmin) / dinc));
        i_report = true;
        get_report_params(i_base, rep_freq);
    }

    disp = dmin;
    displace(pos, top, dmin, dmode);
    rep_count = 1;
    while(disp <= dmax + FLOAT_TOL){
        for(i = 0; i < num_dof; i ++) vel[i] = 0;
        scount = 0;

        if(i_report){
            full_name = report_name(i_base, digits, rep_count);
            intermediate = fopen(full_name.c_str(), "w");
        }

        do{
            /*
            //Do optional MD run without FIRE minimization
            if(precondition){
                simple_md(pos_copy, vel, neighbor_map, bottom, top, mass, pre_dt, pre_steps);
            }
            */
            scount += run_fire_md(num_dof, npairs, pos, vel, pair_list, lengths, fcodes, mass, fcut, max_steps, nmin, fparams, converged, grid, block, intermediate, rep_freq, frms);
            if(! converged && !abort){
                cerr << "Convergence not reached in " << scount << " steps.\n";
                if(! yesno("Continue?")) break;
            }
            else{
                log_line.str("");
                if(converged) log_line << "Convergence reached after " << scount << " steps.\n";
                else log_line << "Covergence not reached after " << scount << " steps.\n";
		log_line << "Ending RMS Force: " << frms << "\n";
                if(make_log){
                    log_file.open(log_name, ofstream::app);
                    log_file << log_line.str();
                    log_file.close();
                }
                else cout << log_line.str();
            }
        }while(! converged && ! abort);

        if(batch_report){
            report_file = report_name(report_base, numd, rep_count);
            report_deformed(npairs, pair_list, pos, report_file);
        }

        pe = get_pe(npairs, pair_list, pos, lengths);
        sim_results.push_back(disp);
        sim_results.push_back(pe);

        if(iform == incrform::constant){
            displace(pos, top, dinc, dmode);
            disp += dinc;
        }
        else{
            prev_disp = disp;
            disp *= d_inc_mult;
            displace(pos, top, disp - prev_disp, dmode);
        }
        rep_count ++;

        if(i_report){
            fclose(intermediate);
        }
    }

    response = enter_decline("Enter a file name for reporting energies");
    if(! response.compare("") == 0){
        ereport = fopen(response.c_str(), "w");
    }

    stream = ereport != NULL ? ereport : stdout;

    for(i = 0; i < sim_results.size() / 2; i++){
        fprintf(stream,"%1.12le\t%1.12le\n",sim_results[2*i],sim_results[2*i+1]);
    }

    if(ereport != NULL) fclose(ereport);

    free(pos);
    free(vel);
    free(lengths);
    free(pair_list);
    free(fcodes);
}

int main(int argc, char **argv){

    string filename, nextline;
    ifstream netfile;
    char c;
    bool abort = false, topSlide = false;
    int dmode = TENSILE, num_read;
    int block_x = 128, block_y = 1;
    incrform iform = incrform::constant;

    while((c = getopt(argc, argv, "alstx:y:")) != -1){
        switch(c) {
            case 'a':
                abort = true;
                break;
            case 'l':
                iform = incrform::logarithmic;
                break;
            case 'x':
                num_read = sscanf(optarg,"%d",&block_x);
                if(num_read < 1){
                    cerr << "Option \"x\" requires an integer.\n";
                }
                else if(block_x < 1){
                    cerr << "Block dimensions must be positive.\n";
                    block_x = 128;
                }
                else{
                    printf("Block x dimension: %d\n", block_x);
                }
                break;
            case 'y':
                num_read = sscanf(optarg,"%d",&block_y);
                if(num_read < 1){
                    cerr << "Option \"y\" requires an integer.\n";
                }
                else if(block_y < 1){
                    cerr << "Block dimensions must be positive.\n";
                    block_y = 1;
                }
                else{
                    printf("Block y dimension: %d\n", block_y);
                }
                break;
            case 's':
                dmode = SHEAR;
                break;
            case 't':
                topSlide = true;
                break;
            case '?':
                if(optopt == 'x'){
                    cerr << "Option \"x\" requires a block dimension.\n";
                }
                if(optopt == 'y'){
                    cerr << "Option \"y\" requires a block dimension.\n";
                }
                else if(isprint(optopt)){
                    fprintf(stderr, "Unknown option: -%c.\n", optopt);
                }
                else{
                    fprintf(stderr, "Unknown option character.\n");
                }
            default:
                break;
         }
    }

    while(yesno("Perform a relaxation run?")){
        do_relaxation_run(abort, dmode, topSlide, dim3(block_x,block_y,1), iform);
    }
}
