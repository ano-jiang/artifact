#include <cassert>
#include "rlibm-fast.h"
#include <omp.h>
#include <algorithm>
#include <execution>
#include <mpi.h>
#include <sys/time.h>
#include <unordered_map>
#include <cstring>
//#include <tbb/concurrent_unordered_map.h>

static unsigned int NUM_THREADS = omp_get_max_threads();

polynomial *rlibm_solve_with_soplex(sample_data *sintervals,
                                    size_t ssize, int *power, int termsize) {

    SoPlex mysoplex;
    mysoplex.setBoolParam(SoPlex::RATFACJUMP, true);
    mysoplex.setIntParam(SoPlex::SOLVEMODE, 2);
    mysoplex.setIntParam(SoPlex::CHECKMODE, 2);
    mysoplex.setIntParam(SoPlex::SYNCMODE, 1);
    mysoplex.setIntParam(SoPlex::READMODE, 1);
    mysoplex.setRealParam(SoPlex::FEASTOL, 0.0);
    mysoplex.setRealParam(SoPlex::OPTTOL, 0.0);
    mysoplex.setRealParam(SoPlex::EPSILON_ZERO, 0.0);
    mysoplex.setRealParam(SoPlex::EPSILON_FACTORIZATION, 0.0);
    mysoplex.setRealParam(SoPlex::EPSILON_UPDATE, 0.0);
    mysoplex.setRealParam(SoPlex::EPSILON_PIVOT, 0.0);
    mysoplex.setIntParam(SoPlex::VERBOSITY, 0);
    mysoplex.setRealParam(SoPlex::TIMELIMIT, 5 * 60);

    /* we first add objective variables */
    DSVectorRational dummycol(0);
    for (int i = 0; i < termsize; i++) {
        auto column = LPColRational(1.0, dummycol, infinity, -infinity);
        mysoplex.addColRational(column);
    }

    /* then add constraints one by one */
//#pragma omp parallel for firstprivate(power, sintervals, ssize, termsize) shared(mysoplex) default(none)
// #pragma omp parallel for shared(mysoplex)
    for (int i = 0; i < ssize; i++) {
        Rational xValR(sintervals[i].x);
        DSVectorRational row1(termsize);

        for (int j = 0; j < termsize; j++) {
            Rational toAdd(1.0);
            for (int k = 0; k < power[j]; k++) toAdd *= xValR;

            row1.add(j, toAdd);
        }

        // LPRow: low bound, row, upper bound
        double lbnd = sintervals[i].lb;
        double ubnd = sintervals[i].ub;
//#pragma omp critical
        mysoplex.addRowRational(LPRowRational(lbnd, row1, ubnd));
    }

    /* solve LP */
    SPxSolver::Status stat;
    stat = mysoplex.optimize();

    /* get solution */
    if (stat == SPxSolver::OPTIMAL) {
        DVectorRational prim(termsize);
        mysoplex.getPrimalRational(prim);

        /* generate the polynomial as a plain structure */
        polynomial *p = (polynomial *) calloc(1, sizeof(polynomial));
        p->termsize = termsize;
        p->power = power;
//        p->coeffs = (double *) calloc(termsize, sizeof(double));

        for (int i = 0; i < termsize; i++)
            p->coeffs[i] = mpq_get_d(*(prim[i].getMpqPtr_w()));

        return p;
    } else if (stat == SPxSolver::UNBOUNDED) {

        polynomial *p = (polynomial *) calloc(1, sizeof(polynomial));
        p->termsize = termsize;
        p->power = power;
//        p->coeffs = (double *) calloc(termsize, sizeof(double));

        for (int i = 0; i < termsize; i++)
            p->coeffs[i] = 0.0;

        return p;
    }

    return nullptr;
}

void check_sorted(sample_info *sampled_indices, size_t ssize) {
    double min = sampled_indices[0].key;

    for (size_t i = 0; i < ssize; i++) {
        assert(min <= sampled_indices[i].key);
        min = sampled_indices[i].key;
    }

}

double rlibm_horner_evaluation(double x, polynomial *poly) {

    double ret_val = 0.0;
    // simulated Horner's method
    for (int i = poly->termsize - 1; i > 0; i--) {
        ret_val = ret_val + poly->coeffs[i];
        double xmul = 1.0;
        for (int j = 0; j < (poly->power[i] - poly->power[i - 1]); j++) {
            xmul = xmul * x;
        }
        ret_val = ret_val * xmul;
    }
    ret_val = ret_val + poly->coeffs[0];

    for (int j = 0; j < poly->power[0]; j++) {
        ret_val = ret_val * x;
    }
    return ret_val;
}


double rlibm_poly_evaluation(double x, polynomial *poly) {

    if (poly->power[0] == 0) {

        //for exp with powers starting from 0

        if (poly->termsize == 1) {
            return poly->coeffs[0];
        } else if (poly->termsize == 2) {
            return fma(x, poly->coeffs[1], poly->coeffs[0]);
        } else if (poly->termsize == 3) {
            double temp1 = fma(x, poly->coeffs[2], poly->coeffs[1]);
            return fma(x, temp1, poly->coeffs[0]);
        } else if (poly->termsize == 4) {
            double xsquare = x * x;
            double temp1 = fma(x, poly->coeffs[3], poly->coeffs[2]);
            double temp2 = fma(x, poly->coeffs[1], poly->coeffs[0]);
            return fma(xsquare, temp1, temp2);
        } else if (poly->termsize == 5) {
            double xsquare = x * x;
            double temp1 = fma(x, poly->coeffs[1], poly->coeffs[0]);
            double temp2 = fma(x, poly->coeffs[3], poly->coeffs[2]);
            double temp3 = fma(xsquare, poly->coeffs[4], temp2);
            return fma(xsquare, temp3, temp1);

            //return (poly->coeffs[0] + x*poly->coeffs[1]) + x*x*(poly->coeffs[2] + x*poly->coeffs[3] + x*x*poly->coeffs[4]);

        } else if (poly->termsize == 6) {
            double xsquare = x * x;
            double temp1 = fma(x, poly->coeffs[1], poly->coeffs[0]);
            double temp2 = fma(x, poly->coeffs[5], poly->coeffs[4]);
            double temp3 = fma(x, poly->coeffs[3], poly->coeffs[2]);
            double temp4 = fma(xsquare, temp2, temp3);
            return fma(xsquare, temp4, temp1);
        }
        return rlibm_horner_evaluation(x, poly);
    }

    if (poly->power[0] == 1) {
        // log with powers {1, 2, 3,4,5}
        if (poly->termsize == 1) {
            return x * poly->coeffs[0];
        } else if (poly->termsize == 2) {
            double temp = x * x * poly->coeffs[1];
            return fma(x, poly->coeffs[0], temp);

            //      return x*poly->coeffs[0] + x*x*poly->coeffs[1];
        } else if (poly->termsize == 3) {
            double temp = x * x * fma(x, poly->coeffs[2], poly->coeffs[1]);
            return fma(x, poly->coeffs[0], temp);

            //     return x*poly->coeffs[0] + x*x*(poly->coeffs[1] + x*poly->coeffs[2]);
        } else if (poly->termsize == 4) {

            double temp1 = fma(x, poly->coeffs[2], poly->coeffs[1]);
            double xsquare = x * x;
            double temp2 = fma(xsquare, poly->coeffs[3], temp1);
            double temp3 = xsquare * temp2;
            return fma(x, poly->coeffs[0], temp3);

            //      return x*poly->coeffs[0] + x*x*(poly->coeffs[1] + x*poly->coeffs[2] + x*x*poly->coeffs[3]);
        } else if (poly->termsize == 5) {
            double xsquare = x * x;
            double temp1 = fma(x, poly->coeffs[4], poly->coeffs[3]);
            double temp2 = fma(x, poly->coeffs[2], poly->coeffs[1]);
            double temp3 = fma(xsquare, temp1, temp2);
            double temp4 = xsquare * temp3;
            return fma(x, poly->coeffs[0], temp4);

            //     return x*poly->coeffs[0] + x*x*(poly->coeffs[1] + x*poly->coeffs[2] + x*x*(poly->coeffs[3] + x*poly->coeffs[4]));
        } else if (poly->termsize == 6) {
            double xsquare = x * x;
            double temp1 = fma(x, poly->coeffs[4], poly->coeffs[3]);
            double temp2 = fma(xsquare, poly->coeffs[5], temp1);

            double temp3 = fma(x, poly->coeffs[2], poly->coeffs[1]);
            double temp4 = fma(xsquare, temp2, temp3);
            double temp5 = xsquare * temp4;
            return fma(x, poly->coeffs[0], temp5);


            //   return x*poly->coeffs[0] + x*x*(poly->coeffs[1] + x*poly->coeffs[2] + x*x*(poly->coeffs[3] + x*poly->coeffs[4] + x*x*poly->coeffs[5]));
        } else if (poly->termsize == 7) {
            double xsquare = x * x;
            double temp1 = fma(x, poly->coeffs[6], poly->coeffs[5]);
            double temp2 = fma(x, poly->coeffs[4], poly->coeffs[3]);
            double temp3 = fma(xsquare, temp1, temp2);
            double temp4 = fma(x, poly->coeffs[2], poly->coeffs[1]);

            double temp5 = xsquare * (temp4 + temp3);
            return fma(x, poly->coeffs[0], temp5);

            //    return x*poly->coeffs[0] + x*x*(poly->coeffs[1] + x*poly->coeffs[2] + x*x*(poly->coeffs[3] + x*poly->coeffs[4] + x*x*(poly->coeffs[5] + x*poly->coeffs[6])));
        }

        return rlibm_horner_evaluation(x, poly);
    }

    //default horner evaluation
    return rlibm_horner_evaluation(x, poly);

}

bool rlibm_validate_and_fix_intervals(sample_data *sintervals,
                                      size_t ssize, polynomial *poly) {

    bool return_val = true;
//#pragma omp parallel for reduction(&&:return_val) firstprivate(ssize, sintervals, poly) default(none)
#pragma omp parallel for
    for (size_t i = 0; i < ssize; i++) {
        double y = rlibm_poly_evaluation(sintervals[i].x, poly);

        if (y < sintervals[i].orig_lb) {
#pragma omp critical
            return_val = false;
            double_x lbx;
            lbx.d = sintervals[i].lb;
            if (lbx.d >= 0) {
                lbx.x = lbx.x + 1;
            } else {
                lbx.x = lbx.x - 1;
            }
            sintervals[i].lb = lbx.d;
        } else if (y > sintervals[i].orig_ub) {
#pragma omp critical
            return_val = false;
            double_x ubx;
            ubx.d = sintervals[i].ub;
            if (ubx.d >= 0) {
                ubx.x = ubx.x - 1;
            } else {
                ubx.x = ubx.x + 1;
            }
            sintervals[i].ub = ubx.d;
        }
    }
    return return_val;
}

polynomial *
rlibm_generate_polynomial(sample_data *sintervals, size_t ssize,
                          int *power, int power_size, int max_tries, int *prev_successful_degree) {

    for (int i = *prev_successful_degree; i < power_size; i++) {
        printf("Trying to generate a polynomial  with %d terms \n", i + 1);

        int count = 0;
        while (count < max_tries) {
            polynomial *p = rlibm_solve_with_soplex(sintervals, ssize, power, i + 1);
            //printf("count is %d ",count);
            if (p && rlibm_validate_and_fix_intervals(sintervals, ssize, p)) {
                *prev_successful_degree = i;
                return p;
            }
            if (p != nullptr) {
                free(p);
            }
            count++;
        }
    }

    return nullptr;

}

int sample_compare_ori(const void *s1, const void *s2) {

    sample_info *e1 = (sample_info *) s1;
    sample_info *e2 = (sample_info *) s2;
    return e1->key > e2->key;
}

int sample_compare(const sample_info &s1, const sample_info &s2) {

    return s1.key > s2.key;
}

int sample_compare_rev(const sample_info &s1, const sample_info &s2) {

    return s1.key < s2.key;
}

void rlibm_print_sample(sample_info *sampled_indices, size_t size) {

    double prev = 0.0;
    for (size_t i = 0; i < size; i++) {
        assert(sampled_indices[i].key >= prev);
        prev = sampled_indices[i].key;
        printf("counter=%lu, key=%e, sample_index=%lu\n", i, sampled_indices[i].key,
               sampled_indices[i].index);
    }
}

void rlibm_weighted_random_sample(sample_info* sampled_indices, size_t ssize,
                                  interval_data* intervals, size_t nentries){

#pragma omp parallel for
    for(size_t i = 0; i < ssize; i++){
        sampled_indices[i].key = pow(intervals[i].u, 1./(intervals[i].w));
        sampled_indices[i].index = i;
    }

    qsort(sampled_indices, ssize, sizeof(sample_info), sample_compare_ori);
    //  check_sorted (sampled_indices, ssize);

    /* select the top ssize indices from the entire interval data */

#pragma omp parallel for
    for(size_t i = ssize; i < nentries; i++){

        /* if the ith element is smaller than the least element in the
           sample, then nothing to do */
        size_t j= 0;

        double interval_key = pow(intervals[i].u, 1./(intervals[i].w));

        if(interval_key > sampled_indices[0].key){
            /* do insertion sort of the data */
            while(interval_key > sampled_indices[j].key && j < ssize) j++;

            for(size_t t=1; t < j; t++) {
                sampled_indices[t-1] = sampled_indices[t];
            }
            sampled_indices[j-1].key = interval_key;
            sampled_indices[j-1].index = i;
        }
    }
    //  check_sorted(sampled_indices, ssize);
}


size_t
rlibm_compute_violated_indices(size_t *violated_indices, interval_data *intervals, size_t nentries, polynomial *poly,
                               double *extra_err) {

    size_t num_violated_indices = 0;
#pragma omp parallel for shared(num_violated_indices)
    for (size_t i = 0; i < nentries; ++i) {
        double y = rlibm_poly_evaluation(intervals[i].x, poly);
        if (y < intervals[i].lb || y > intervals[i].ub) {
            violated_indices[num_violated_indices] = i;
#pragma omp atomic
            num_violated_indices++;
        }
    }
    return num_violated_indices;
}

void rlibm_evaluate_and_update_weights(size_t *violated_indices, size_t num_violated_indices,
                                       interval_data *intervals, size_t nentries, size_t d, double *ulps) {
    double w_v = 0.0;
    double w_s = 0.0;

    // this can be optimized to only change the updated weights. For now
    // using an inefficient one
#pragma omp parallel for reduction(+:w_s)
    for (size_t i = 0; i < nentries; i++) {
        w_s = w_s + intervals[i].w;
    }

#pragma omp parallel for reduction(+:w_v)
    for (size_t i = 0; i < num_violated_indices; i++) {
        w_v = w_v + intervals[violated_indices[i]].w;
    }

    //doubling the weights for a lucky iteration
    if (w_v <= 2 * w_s / (9 * d - 1)) {
#pragma omp parallel for
        for (size_t i = 0; i < num_violated_indices; i++) {
            size_t vindex = violated_indices[i];
            intervals[vindex].w = intervals[vindex].w * 2;
        }
    }
}

void
rlibm_regenerate_random_values_and_reset_weights(interval_data *intervals,
                                                 size_t nentries) {

    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

//#pragma omp parallel for firstprivate(nentries, intervals) shared(distribution, generator) default(none)
#pragma omp parallel for
    for (size_t i = 0; i < nentries; i++) {
        intervals[i].u = distribution(generator);
        intervals[i].w = 1.0;
    }
}

bool check_sampled_indices(sample_info *sample, sample_info *prev_sample, size_t size) {

    for (size_t i = 0; i < size; i++) {
        if (sample[i].index != prev_sample[i].index) {
            return false;
        }
    }
    return true;
}

void rlibm_print_polyinfo(polynomial *p) {

    if (p->termsize == 0) {
        printf("Polynomial has no terms!\n");
        exit(0);
    }

    printf("Polynomial: y=%a x^(%d)", p->coeffs[0], p->power[0]);
    for (int j = 1; j < p->termsize; j++) {
        printf(" + %a x^(%d)", p->coeffs[j], p->power[j]);
    }
    printf("\n");

}

interval_data *rlibm_read_interval_file(char *filename, size_t *nentries) {

    FILE *fp = fopen(filename, "r");
    assert(fp != nullptr);

    /* count the number of entries */

    fseek(fp, 0, SEEK_END);
    size_t nentries_total = ftell(fp);
    nentries_total = nentries_total / (3 * sizeof(double));
    printf("number of intervals = %lu\n", nentries_total);
    fseek(fp, 0, SEEK_SET);

    interval_data *intervals = (interval_data *) calloc(nentries_total, sizeof(interval_data));

    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    //#pragma omp parallel for
    for (unsigned long i = 0; i < nentries_total; i++) {
        double data_entry[3];
        size_t bytes = fread(data_entry, sizeof(double), 3, fp);
        intervals[i].w = 1.0;
        intervals[i].u = distribution(generator);
        intervals[i].x = data_entry[0];
        intervals[i].lb = data_entry[1];
        intervals[i].ub = data_entry[2];
    }

    fclose(fp);
    *nentries = nentries_total;
    return intervals;
}


int main(int argc, char **argv) {
    struct timeval start, end, func_beg, func_end;
    double execution_time_sec, execution_time2_sec, execution_time3_sec, exec_time[4]={0.0, 0.0, 0.0, 0.0};
    int itr_best = -1, itr=-1;

    gettimeofday(&start, NULL);
    if (argc != 3) {
        printf("Usage: %s <config_file> <interval file> \n", argv[0]);
        exit(0);
    }

//    printf("max number of thread is %d\n", omp_get_max_threads());
    printf("EXIT_ON_THRESHOLD is %d\n", RLIBM_EXIT_ON_THRESHOLD);
    printf("NUM_THREADS is %d\n", NUM_THREADS);

    FILE *config_file = fopen(argv[1], "r");

    int *powers = NULL;
    int MY_RLIBM_PIECES = 1;
    int read_value = fscanf(config_file, "%d\n", &MY_RLIBM_PIECES);
    assert(read_value == 1);
    int MY_TERMS = 1;
    read_value = fscanf(config_file, "%d ", &MY_TERMS);
    assert(read_value == 1);

    powers = (int *) malloc(MY_TERMS * sizeof(int));
    for (int i = 0; i < MY_TERMS; i++) {
        read_value = fscanf(config_file, "%d", &powers[i]);
        assert(read_value == 1);
    }

    int min_violate[MY_RLIBM_PIECES];
    for (int i = 0; i < MY_RLIBM_PIECES; ++i) {
        read_value = fscanf(config_file, "%d", &min_violate[i]);
        assert(read_value == 1);
    }

    fclose(config_file);
    size_t nentries_total = 0;
    interval_data *intervals_full = rlibm_read_interval_file(argv[2], &nentries_total);

    for (int i = 0; i < MY_TERMS; i++) {
        printf("%d ", powers[i]);
    }

    // int powers[]={1, 3, 5, 7, 9};
    int powers_size = MY_TERMS;

    //log2, log10
    //  int powers[]={0, 1,2,3,4,5,6, 7, 8, 9}; //coshx, exp2, exp10, exp
    // int powers[] = {1,3,5,7}; //sinhx
    // int powers[] = {0,2 ,4,6}; // coshx
    //int powers[] = {1,3 ,5,7}; // sinpi
    // int powers[] = {0,2,4,6}; // cospi
    size_t nentries_pieces = nentries_total / MY_RLIBM_PIECES;
    int set_threads = NUM_THREADS / MY_RLIBM_PIECES;
    if (set_threads == 0) {
        set_threads = 1;
    }
    printf("\nTotal Threads is %d\n", set_threads);
    NUM_THREADS = set_threads;
    omp_set_num_threads(set_threads);

    int prev_successful_degree = 0;

    // Initialize the MPI environment
    MPI_Init(&argc, &argv);
    // Find out rank, size
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    // for (int i = 0; i < MY_RLIBM_PIECES; i++) {

    prev_successful_degree = 0;
    size_t start_index = world_rank * nentries_pieces;

    size_t even_split = (world_rank + 1) * nentries_pieces;
    size_t end_index = even_split > nentries_total ? nentries_total : even_split;

    size_t nentries = end_index - start_index;
    interval_data *intervals = &intervals_full[start_index];


    printf("piece = %d\n", world_rank);
    printf("start_index = %lu\n", start_index);
    printf("end_index = %lu\n", end_index);
    printf("nentries=%lu\n", nentries);

    /* sample size */

    size_t cd = 9 * powers_size * powers_size;
    size_t samplesize = cd;

    size_t n_violated_indices = 0;
    size_t *violated_indices = (size_t *) calloc(nentries, sizeof(size_t));
    double *ulps = (double *) calloc(nentries, sizeof(double));

    sample_info *sampled_indices = (sample_info *) calloc(cd, sizeof(sample_info));

    size_t prev_violated_indices = 0;
    size_t matched_violated_indices = 0;

    sample_data *sampled_intervals = (sample_data *) calloc(cd, sizeof(sample_data));

    polynomial *p = nullptr;
    size_t total_iterations = 0;

    bool is_print=false;
    bool is_record=false;
    do {
        if (p != nullptr) free(p);

        n_violated_indices = 0;

        gettimeofday(&func_beg, NULL);
        rlibm_weighted_random_sample(sampled_indices, cd, intervals, nentries);
        gettimeofday(&func_end, NULL);
        exec_time[0]+=((func_end.tv_sec - func_beg.tv_sec) + (func_end.tv_usec - func_beg.tv_usec) / 1000000.0);

        total_iterations++;
#pragma omp parallel for
        for (size_t i = 0; i < cd; i++) {
            size_t iindex = sampled_indices[i].index;

            sampled_intervals[i].x = intervals[iindex].x;
            sampled_intervals[i].lb = intervals[iindex].lb;
            sampled_intervals[i].ub = intervals[iindex].ub;
            sampled_intervals[i].orig_lb = sampled_intervals[i].lb;
            sampled_intervals[i].orig_ub = sampled_intervals[i].ub;
            sampled_intervals[i].w = intervals[iindex].w;
            sampled_intervals[i].u = intervals[iindex].u;
            sampled_intervals[i].k = sampled_indices[i].key;
        }

        gettimeofday(&func_beg, NULL);
        /* need to implement these functions */
        p = rlibm_generate_polynomial(sampled_intervals, samplesize, powers, powers_size, MAX_TRIES,
                                      &prev_successful_degree);
        gettimeofday(&func_end, NULL);
        exec_time[1]+=((func_end.tv_sec - func_beg.tv_sec) + (func_end.tv_usec - func_beg.tv_usec) / 1000000.0);


        if (p) {
            gettimeofday(&func_beg, NULL);
            n_violated_indices = rlibm_compute_violated_indices(violated_indices, intervals, nentries, p, ulps);
            gettimeofday(&func_end, NULL);
            exec_time[2]+=((func_end.tv_sec - func_beg.tv_sec) + (func_end.tv_usec - func_beg.tv_usec) / 1000000.0);

            printf("p%d: number of violated intervals: %lu, total iterations=%lu \n", world_rank,
                   n_violated_indices, total_iterations);

            if (n_violated_indices <= VIOLATE_THRESHOLD) {
                printf("p%d: VIOLATING INPUTS BELOW THRESHOLD:\n", world_rank);
                printf("p%d: starting input is %a\n", world_rank, intervals[0].x);

                for (size_t m = 0; m < n_violated_indices; m++) {
                    printf("p%d: violated_input is %a, lb is %a, ub is %a\n", world_rank,
                           intervals[violated_indices[m]].x,
                           intervals[violated_indices[m]].lb, intervals[violated_indices[m]].ub);
                }
                rlibm_print_polyinfo(p);

                if(!is_print){
                    gettimeofday(&end, NULL); // 记录结束时间
                    execution_time2_sec = ((end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0);
                    itr_best=total_iterations;
                    is_print=true;
                }
//                break;
//                if(!is_record && n_violated_indices<min_violate[world_rank]){
//                    gettimeofday(&end, NULL);
//                    execution_time3_sec=((end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0);
//                    is_record=true;
//                    itr=total_iterations;
//                }


                if (RLIBM_EXIT_ON_THRESHOLD) {
                    break;
                }
            }

            gettimeofday(&func_beg, NULL);
            rlibm_evaluate_and_update_weights(violated_indices, n_violated_indices, intervals, nentries,
                                              powers_size, ulps);
            gettimeofday(&func_end, NULL);
            exec_time[3]+=((func_end.tv_sec - func_beg.tv_sec) + (func_end.tv_usec - func_beg.tv_usec) / 1000000.0);


        } else {
            if (total_iterations > MAX_ITERATIONS) {
                printf("total iterations exceeded %d, terminating the polynomial geenerator\n", MAX_ITERATIONS);
                if (p != nullptr) {
                    free(p);
                    p = nullptr;
                }
                break;
            }
            printf("p%d: failed to generate polynomial, resetting weights, total_iterations=%lu\n", world_rank,
                   total_iterations);
            prev_successful_degree = 0;
            rlibm_regenerate_random_values_and_reset_weights(intervals, nentries);
        }

        /* debugging feature to reset weights for the sample if not making progress*/
        if (n_violated_indices != 0 && (prev_violated_indices == n_violated_indices)) {
            matched_violated_indices++;
            if (matched_violated_indices > SAMPLE_MATCH_THRESHOLD) {
                matched_violated_indices = 0;
                n_violated_indices = 0;

                printf("p%d: not making progress, same number of violated indices, resetting weights, total_iterations=%lu\n",
                       world_rank, total_iterations);
                prev_successful_degree = 0;
                rlibm_regenerate_random_values_and_reset_weights(intervals, nentries);
                if (p != nullptr) {
                    free(p);
                    p = nullptr;
                }
                continue;
            }
        } else {
            matched_violated_indices = 0;
            prev_violated_indices = n_violated_indices;
        }
        //} while (n_violated_indices > min_violate[world_rank] || !p);
    } while (n_violated_indices > 0 || !p);

    if (p) {
        rlibm_print_polyinfo(p);
    } else {
        printf("p%d: Could not generate the polynomial that satisifies all intervals, check for partial results with a few violated intervals\n",
               world_rank);
    }
    free(p);
    free(sampled_intervals);
    free(sampled_indices);
    free(violated_indices);
    free(ulps);
    // } // ends the outer piecewise loop

    free(intervals_full);

    gettimeofday(&end, NULL);
    execution_time_sec = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
    printf("p%d: program reached 7 conflicts in round %d, execution time：%lf seconds\n", world_rank, itr_best, execution_time2_sec);
    for (int i = 0; i < 4; ++i) {
        printf("p%d: function index %d, Cumulative execution time：%lf seconds\n", world_rank, i, exec_time[i]);
    }
    printf("p%d: program execution time: %lf seconds\n", world_rank, execution_time_sec);

    MPI_Finalize();
    return 0;
}
