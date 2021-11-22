/*
by Osnat Lifshitz, Chemi Gotlibovski (2009)
*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <limits.h>
#include <assert.h>

#ifndef LINE_MAX
    #define LINE_MAX 256
#endif

// max macro with index
#define get_max_idx(max_val,idx_of_max,val,idx) if ((val) > (max_val)) { max_val = (val); idx_of_max = (short)(idx); }
// max macro without index
#define get_max(max_val,val) max_val = (((val) > (max_val)) ? (val) : (max_val))
// travel cost matrix
#define travel_cost(rg1,rg2) travel_cost_arr[(rg1)][(rg2)]
// wage discretization macro
#define get_discrete_index(w) (unsigned short)(((w) < P_W_ERROR_RNG[1]) ?  1 : (((w) < P_W_ERROR_RNG[2]) ? 2 : (((w) < P_W_ERROR_RNG[3]) ? 3 : (((w) <P_W_ERROR_RNG[4]) ? 4 : 5))))

/* 
Implements the Polar form of the Box-Muller Transformation
(c) Copyright 1994, Everett F. Carter Jr.
Permission is granted by the author to use
this software for any application provided this
copyright notice is preserved.
*/
inline float rand01()
{
    return (float)rand()/(float)RAND_MAX;
}

inline float randn01()
{                       
    float x1, x2, w, y1;
    static float y2;
    static int use_last = 0;
    /* use value from previous call */
    if (use_last)
    {
        y1 = y2;
        use_last = 0;
    } else
    {
        do
        {
            x1 = 2.0f * rand01() - 1.0f;
            x2 = 2.0f * rand01() - 1.0f;
            w = x1 * x1 + x2 * x2;
        } while ( w >= 1.0f );

        w = sqrtf( (-2.0f * logf( w ) ) / w );
        y1 = x1 * w;
        y2 = x2 * w;
        use_last = 1;
    }

    return y1;
}

// program constant
const unsigned int TYPE_SIZE    = 3;                // # of types
const unsigned short T          = 20;               // time
#ifdef WIFE_MODE
const unsigned int OBSR         = 529;              // individual 
#else // regular mode
const unsigned int OBSR         = 697;              // individual 
#endif
const unsigned int OBS          = OBSR*TYPE_SIZE;   // individual multiplay by number of types
const unsigned int DRAWS        = 30;               // draws for emax
const unsigned int RG_SIZE      = 7;                // # of regions
const unsigned int TC_SIZE      = 3;                // travel cost
const unsigned int STATE_SIZE   = 3;                // # of states: w,b,ue
#ifdef SIMULATION
const unsigned int BASE_DRAWS_F = 1000;
const unsigned int TYPE_0_OF_1000 = BASE_DRAWS_F*0.085;
const unsigned int TYPE_1_OF_1000 = BASE_DRAWS_F*0.605;
const unsigned int TYPE_2_OF_1000 = BASE_DRAWS_F*0.310;
const unsigned int DRAWS_F      = TYPE_1_OF_1000;   // max of: T0 = 85, T1 = 605, T2 = 310
#else
#ifdef WAGE_SELECTION
const unsigned int DRAWS_F      = 333;              // draws for forward solving
#else
const unsigned int DRAWS_F      = 333;              // draws for forward solving
#endif // WAGE_SELECTION
#endif // SIMULATION
const unsigned int D_WAGE       = 6;                //
const unsigned int TAU          = 50000;
const unsigned int STATE_VECTOR_SIZE = RG_SIZE*RG_SIZE + 2*RG_SIZE; // 63 states

// random draws
float randn_b_arr[DRAWS][OBS][T][RG_SIZE][STATE_SIZE];
#define epsilon_b(draw, I, t, rg, state) randn_b_arr[(draw)][(I)][(t)][(rg)][(state)]

#ifdef SIMULATION
float* randn_f_arr = new float[DRAWS_F*OBS*T*RG_SIZE*STATE_SIZE];
#define epsilon_f(draw, I, t, rg, state) randn_f_arr[((((draw)*OBS + (I))*T + (t))*RG_SIZE + (rg))*STATE_SIZE + (state)]
#else
float randn_f_arr[DRAWS_F][OBS][T][RG_SIZE][STATE_SIZE];
#define epsilon_f(draw, I, t, rg, state) randn_f_arr[(draw)][(I)][(t)][(rg)][(state)]
#endif // SIMULATION

#define draw_wage(wage,prob) (rand01() < (prob)) ? (wage) : -INFINITY
#define draw_wage_f(wage,prob) (rand01() < (prob)) ? (wage) : -INFINITY

inline unsigned short draw_type(float prob_type1, float prob_type2)
{
    float p = rand01();
    if (p < prob_type1)
    {
        return 1;
    }
    else if (p < prob_type1 + prob_type2)
    {
        return 2;
    }
    else
    {
        return 0;
    }
}

static void init_rand()
{
    srand(12345);
    
    for (unsigned short i = 0; i < OBS; ++i)
    {
        for (unsigned short t = 0; t < T; ++t)
        {
            for (unsigned short rg = 0; rg < RG_SIZE; ++rg)
            {
                for (unsigned short s = 0; s < STATE_SIZE; ++s)
                {
                    for (unsigned short d = 0; d < DRAWS; ++d)
                    {
                        epsilon_b(d, i, t, rg, s) = randn01();
                    }
                    for (unsigned short d = 0; d < DRAWS_F; ++d)
                    {
                        epsilon_f(d, i, t, rg, s) = randn01();
                    }
                }
            }
        }
    }
}

//////////////////////// profiling - start - counting time of runing

inline timeval tic()
{
    timeval tv;
    gettimeofday(&tv, NULL);
    return tv;
}

inline float toc(timeval tv1)
{
    timeval tv2;
    gettimeofday(&tv2, NULL);
    return (float)(tv2.tv_sec-tv1.tv_sec) + ((float)(tv2.tv_usec - tv1.tv_usec))/1000000.0f;
}

#ifdef PERF_TRACE
const float ITER_COUNT = (float)(OBS*T*(T+1))/2.0f;
#endif

//////////////////////// profiling - end

// global variables loaded from individual's data file
unsigned short M_arr[OBS];         // married / single
unsigned short KIDS_arr[OBS];      //# of children
unsigned short EXP_U_arr[OBS];     // experience in former USSR
unsigned short SCHOOL_arr[OBS];    // years of education
float          WAGE_arr[OBS];      // wage of ind
unsigned short AGE_arr[OBS];       // age at arrival
unsigned short PERIODS_arr[OBS];   // number of periods (<T)
unsigned long  RENT_MORT_arr[OBS]; // cost of housing
unsigned short D_MORT_arr[OBS];    // dummy
unsigned short REP1_arr[OBS];      // republic 1 at USSR
unsigned short REP2_arr[OBS];      // republic 2 at USSR
unsigned short REP3_arr[OBS];      // republic 3 at USSR
unsigned short TYPE2_arr[OBS];     // dummy that get 1 of individual is type 2
unsigned short TYPE3_arr[OBS];     // dummy that get 1 of individual is type 2
unsigned short WIFE_EDU_arr[OBS];  // wife education array

// loading from file to global variables
static bool load_individuals(const char* filename)
{
    const int COLUMN_NUMBER = 14;

    FILE* fp = fopen(filename,"r");
    if (fp)
    {
        char line[LINE_MAX];
        unsigned short I = 0;
        int col_num = 0;
#ifdef TRACE_LOAD
        printf("Loading %s\n", filename);
        printf("============================================================\n");
#endif
        while (I < OBS)
        {
            if (fgets(line, LINE_MAX, fp) != 0)
            {
                // formatting the line according to:
                // M KIDS EXP_U SCHOOL WAGE AGE PERIODS RENT_MORT D_MORT REP1 REP2 REP3 TYPE2 TYPE3
                col_num = sscanf(line, "%hu%hu%hu%hu%f%hu%hu%lu%hu%hu%hu%hu%hu%hu",
                                 &(M_arr[I]), &(KIDS_arr[I]), &(EXP_U_arr[I]), &(SCHOOL_arr[I]),
                                 &(WAGE_arr[I]), &(AGE_arr[I]), &(PERIODS_arr[I]), &(RENT_MORT_arr[I]), &(D_MORT_arr[I]), &(REP1_arr[I]), &(REP2_arr[I]),
                                 &(REP3_arr[I]), &(TYPE3_arr[I]), &TYPE2_arr[I]);
                if (col_num != COLUMN_NUMBER)
                {
                    printf("wrong format in file %s number of columns: %d \n", filename, col_num);
                    printf("line[%hu]: %s\n", I, line);
                    printf("format: M\tKIDS\tEXP_U\tSCHOOL\tWAGE\tAG\tPERIODS\tRENT_MO\tD_MO\tREP1\tREP2\tREP3\tTYPE3\tTYPE2\n");
                    fclose(fp);
                    return false;
                } 
                else
                {
#ifdef TRACE_LOAD
                    printf("line[%hu]: %hu\t%hu\t%hu\t%hu\t%f\t%hu\t%hu\t%lu\t%hu\t%hu\t%hu\t%hu\t%hu\t%hu\n", I, M_arr[I], KIDS_arr[I],
                         EXP_U_arr[I], SCHOOL_arr[I],
                         WAGE_arr[I], AGE_arr[I], PERIODS_arr[I], RENT_MORT_arr[I], D_MORT_arr[I], REP1_arr[I], REP2_arr[I],
                         REP3_arr[I], TYPE3_arr[I], TYPE2_arr[I]);
#endif
                    ++I;
                }
            } 
            else
            {
                printf("erorr [%s] reading file [%s] at line [%hu]\n", strerror(errno), filename, I);
                fclose(fp);
                return false;
            }
        }
        fclose(fp);
        return true;
    } 
    else
    {
        printf("failed to open individual file %s\n", filename);
        return false;
    }
}

#ifdef SIMULATION
unsigned short IND_FILTER_arr[OBS];  // filter out individuals for summary

static bool load_individuals_filter(const char* filename)
{
    const int COLUMN_NUMBER = 1;

    FILE* fp = fopen(filename,"r");
    if (fp)
    {
        char line[LINE_MAX];
        unsigned short I = 0;
        int col_num = 0;
#ifdef TRACE_LOAD
        printf("Loading %s\n", filename);
        printf("============================================================\n");
#endif
        while (I < OBS)
        {
            if (fgets(line, LINE_MAX, fp) != 0)
            {
                // formatting the line according to:
                col_num = sscanf(line, "%hu", &(IND_FILTER_arr[I]));
                if (col_num != COLUMN_NUMBER)
                {
                    printf("wrong format in file %s number of columns: %d \n", filename, col_num);
                    printf("line[%hu]: %s\n", I, line);
                    fclose(fp);
                    return false;
                } 
                else
                {
                    if (IND_FILTER_arr[I] != 0 && IND_FILTER_arr[I] != 1)
                    {
                        printf("line[%hu]: expected values for filter are 0/1 only (%hu)\n", I, IND_FILTER_arr[I]);
                        fclose(fp);
                        return false;
                    }
#ifdef TRACE_LOAD
                    printf("line[%hu]: %hu\n", I, IND_FILTER_arr[I]);
#endif
                    ++I;
                }
            } 
            else
            {
                printf("erorr [%s] reading file [%s] at line [%hu]\n", strerror(errno), filename, I);
                fclose(fp);
                return false;
            }
        }
        fclose(fp);
        return true;
    } 
    else
    {
        printf("failed to open individual filter file %s\n", filename);
        return false;
    }
}

float FR[RG_SIZE];

static bool load_future_intrest(const char* filename)
{
    const int COLUMN_NUMBER = 7;
    FILE* fp = fopen(filename, "r");
    if (fp)
    {
        char line[LINE_MAX];
        int col_num = 0;
#ifdef TRACE_LOAD
        printf("Loading %s\n", filename);
        printf("============================================================\n");
#endif
        if (fgets(line, LINE_MAX, fp) != 0)
        {
            col_num = sscanf(line, "%f%f%f%f%f%f%f", &FR[0], &FR[1], &FR[2], &FR[3], &FR[4], &FR[5], &FR[6]);
            if (col_num != COLUMN_NUMBER)
            {
                printf("wrong format in file %s number of columns: %d \n", filename, col_num);
                printf("line: %s\n", line);
                fclose(fp);
                return false;
            }
            else
            {
#ifdef TRACE_LOAD
                    printf("line: %s", line);
#endif
            }
        }
        fclose(fp);
        return true;
    }
    else
    {

        printf("failed to open future intrest file %s\n", filename);
        return false;
    }
}

#endif // SIMULATION

// loading from file to global variables
static bool load_wife_edu(const char* filename)
{
    const int COLUMN_NUMBER = 1;

    FILE* fp = fopen(filename,"r");
    if (fp)
    {
        char line[LINE_MAX];
        unsigned short I = 0;
        int col_num = 0;
#ifdef TRACE_LOAD
        printf("Loading %s\n", filename);
        printf("============================================================\n");
#endif
        while (I < OBS)
        {
            if (fgets(line, LINE_MAX, fp) != 0)
            {
                col_num = sscanf(line, "%hu", &(WIFE_EDU_arr[I]));
                if (col_num != COLUMN_NUMBER)
                {
                    printf("wrong format in file %s number of columns: %d \n", filename, col_num);
                    printf("line[%hu]: %s\n", I, line);
                    printf("format: WIFE_EDU\n");
                    fclose(fp);
                    return false;
                } 
                else
                {
#ifdef TRACE_LOAD
                    printf("line[%hu]: %hu\n", I, WIFE_EDU_arr[I]);
#endif
                    ++I;
                }
            } 
            else
            {
                printf("erorr [%s] reading file [%s] at line [%hu]\n", strerror(errno), filename, I);
                fclose(fp);
                return false;
            }
        }
        fclose(fp);
        return true;
    } 
    else
    {
        printf("failed to open wife education file %s\n", filename);
        return false;
    }
}

const unsigned short  UE    = 0;
const unsigned short  BLUE  = 2;
const unsigned short  WHITE = 1;

#ifdef WIFE_MODE
const unsigned short MOMENTS_PERIODS = 12;
#else // regular mode
const unsigned short MOMENTS_PERIODS = 13;
#endif

short occupation_arr[OBS][MOMENTS_PERIODS]; // real occ
short live_arr[OBS][MOMENTS_PERIODS];       // real housing region
short work_arr[OBS][MOMENTS_PERIODS];       // real work region
short sample_arr[OBS][MOMENTS_PERIODS];     // real state 0-63
#define occupation(i,j) occupation_arr[(i)][(j)]
#define live(i,j) live_arr[(i)][(j)]
#define work(i,j) work_arr[(i)][(j)]
#define sample(i,j) sample_arr[(i)][(j)]

static bool load_moments(const char* filename)
{
#ifdef WIFE_MODE
    const int COLUMN_NUMBER = 49; 
#else // regular mode
    const int COLUMN_NUMBER = 53;
#endif
    FILE* fp = fopen(filename,"r");
    if (fp)
    {
        char line[LINE_MAX];
        unsigned short I = 0;
        int col_num = 0;
        unsigned int familyid;
#ifdef TRACE_LOAD
        printf("Loading %s\n", filename);
        printf("============================================================\n");
#endif
        while (I < OBS)
        {
            if (fgets(line, LINE_MAX, fp) != 0)
            {
#ifdef WIFE_MODE
                col_num = sscanf(line, 
                     "%u%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd",
                     &(familyid), 
                     &(occupation(I,0)), &(live(I,0)), &(work(I,0)), &(sample(I,0)),
                     &(occupation(I,1)), &(live(I,1)), &(work(I,1)), &(sample(I,1)),
                     &(occupation(I,2)), &(live(I,2)), &(work(I,2)), &(sample(I,2)),
                     &(occupation(I,3)), &(live(I,3)), &(work(I,3)), &(sample(I,3)),
                     &(occupation(I,4)), &(live(I,4)), &(work(I,4)), &(sample(I,4)),
                     &(occupation(I,5)), &(live(I,5)), &(work(I,5)), &(sample(I,5)),
                     &(occupation(I,6)), &(live(I,6)), &(work(I,6)), &(sample(I,6)),
                     &(occupation(I,7)), &(live(I,7)), &(work(I,7)), &(sample(I,7)),
                     &(occupation(I,8)), &(live(I,8)), &(work(I,8)), &(sample(I,8)),
                     &(occupation(I,9)), &(live(I,9)), &(work(I,9)), &(sample(I,9)),
                     &(occupation(I,10)), &(live(I,10)), &(work(I,10)), &(sample(I,10)),
                     &(occupation(I,11)), &(live(I,11)), &(work(I,11)), &(sample(I,11)));
#else // regular mode
                col_num = sscanf(line, 
                     "%u%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd%hd",
                     &(familyid), 
                     &(occupation(I,0)), &(live(I,0)), &(work(I,0)), &(sample(I,0)),
                     &(occupation(I,1)), &(live(I,1)), &(work(I,1)), &(sample(I,1)),
                     &(occupation(I,2)), &(live(I,2)), &(work(I,2)), &(sample(I,2)),
                     &(occupation(I,3)), &(live(I,3)), &(work(I,3)), &(sample(I,3)),
                     &(occupation(I,4)), &(live(I,4)), &(work(I,4)), &(sample(I,4)),
                     &(occupation(I,5)), &(live(I,5)), &(work(I,5)), &(sample(I,5)),
                     &(occupation(I,6)), &(live(I,6)), &(work(I,6)), &(sample(I,6)),
                     &(occupation(I,7)), &(live(I,7)), &(work(I,7)), &(sample(I,7)),
                     &(occupation(I,8)), &(live(I,8)), &(work(I,8)), &(sample(I,8)),
                     &(occupation(I,9)), &(live(I,9)), &(work(I,9)), &(sample(I,9)),
                     &(occupation(I,10)), &(live(I,10)), &(work(I,10)), &(sample(I,10)),
                     &(occupation(I,11)), &(live(I,11)), &(work(I,11)), &(sample(I,11)),
                     &(occupation(I,12)), &(live(I,12)), &(work(I,12)), &(sample(I,12)));
#endif
                if (col_num != COLUMN_NUMBER)
                {
                    printf("wrong format in file %s number of columns: %d \n", filename, col_num);
                    printf("line[%hu]: %s\n", I, line);
                    fclose(fp);
                    return false;
                } 
                else
                {
#ifdef TRACE_LOAD
                    printf("line[%hu]: %s", I, line);
#endif
                    ++I;
                }
            } 
            else
            {
                printf("erorr [%s] reading file [%s] at line [%hu]\n", strerror(errno), filename, I);
                fclose(fp);
                return false;
            }
        }
        fclose(fp);
        return true;
    } 
    else
    {
        printf("failed to open moments file %s\n", filename);
        return false;
    }
}

const unsigned short int MAX_PARAM_LEN = 158;   //# of parameters
#define set_param_array(param_name,size) float param_name[(size)]; for (j = 0; j < (size); ++i, ++j) param_name[j] = params[i]; 
#define set_param(param_name) float param_name = params[i]; ++i;

static bool load_dynamic_params(const char* filename, float* params, unsigned short* dynamic_param_idx, unsigned short dynamic_param_size)
{
    const int COLUMN_NUMBER = 1;
    FILE* fp = fopen(filename,"r");
    bool use_dynamic = (dynamic_param_size != 0);
    if (!use_dynamic)
    {
        dynamic_param_size = MAX_PARAM_LEN;
    }
    if (fp)
    {
        char line[LINE_MAX];
        unsigned short int line_idx = 0;
        int col_num = 0;
        float tmp_param;
#ifdef TRACE_LOAD
        printf("Loading %s\n", filename);
        printf("============================================================\n");
#endif
        while (line_idx <= dynamic_param_size)
        {
            if (fgets(line, LINE_MAX, fp) != 0)
            {
                col_num = sscanf(line, "%f", &tmp_param);
                if (col_num != COLUMN_NUMBER)
                {
                    printf("wrong format in file %s number of columns: %d \n", filename, col_num);
                    printf("line[%hu]: %s\n", line_idx, line);
                    fclose(fp);
                    return false;
                }
                else
                {
                    if(line_idx == 0)
                    {
                        // first param is the number of params
                        if ((unsigned short)tmp_param != dynamic_param_size)
                        {
                            printf("wrong number of params in file %s is %hu, should be %hu\n", filename, (unsigned short)tmp_param, dynamic_param_size);
                            fclose(fp);
                            return false; 
                        }
                    }
                    else
                    {
                        if (use_dynamic)
                        {
                            params[dynamic_param_idx[line_idx-1]] = tmp_param;
                        }
                        else
                        {
                            params[line_idx-1] = tmp_param;
                        }
#ifdef TRACE_LOAD
                        printf("line[%d]: param[%d]=%f\n", line_idx, dynamic_param_idx[line_idx-1], tmp_param);
#endif
                    }
                    ++line_idx;
                }
            } else
            {
                printf("erorr [%s] reading file [%s] at line [%d]\n", strerror(errno), filename, line_idx);
                fclose(fp);
                return false;
            }

        }
        fclose(fp);
        return true;
    } else
    {

        printf("failed to open starting point file %s\n", filename);
        return false;
    }
}

static bool load_const_params(const char* filename, float* x)
{
    const int COLUMN_NUMBER = 2;
    FILE* fp = fopen(filename,"r");
    if (fp)
    {
        char line[LINE_MAX];
        unsigned short line_idx = 0;
        int col_num = 0;
        float tmp_param;
        unsigned short param_idx;
#ifdef TRACE_LOAD
        printf("Loading %s\n", filename);
        printf("============================================================\n");
#endif
        while (line_idx < MAX_PARAM_LEN)
        {
            if (fgets(line, LINE_MAX, fp) != 0)
            {
                col_num = sscanf(line, "%hu%f", &param_idx, &tmp_param);
                if (col_num != COLUMN_NUMBER)
                {
                    printf("wrong format in file %s number of columns: %d \n", filename, col_num);
                    printf("line[%hu]: %s\n", line_idx, line);
                    fclose(fp);
                    return false;
                } 
                else if (param_idx >= MAX_PARAM_LEN)
                {
                    printf("line[%hu]: invalid param index: %hu\n", line_idx, param_idx);
                    fclose(fp);
                    return false;
                } 
                else
                {
                    x[param_idx] = tmp_param;
#ifdef TRACE_LOAD
                    printf("line[%hu]: %s", line_idx, line);
#endif
                    ++line_idx;
                }
            } 
            else
            {
                printf("erorr [%s] reading file [%s] at line [%d]\n", strerror(errno), filename, line_idx);
                fclose(fp);
                return false;
            }

        }
        fclose(fp);
        return true;
    } 
    else
    {

        printf("failed to open starting point file %s\n", filename);
        return false;
    }
}

// EMAX matrix
// Note: the loop of T goes from 1 to T, C++ arrays are 0 to T-1, so the size must be T+1 
#ifndef WAGE_SLECTION
float* EMAX_mat = (float*)malloc((T+1)*(T+2)*RG_SIZE*RG_SIZE*STATE_SIZE*D_WAGE*TYPE_SIZE*sizeof(float));
#define EMAX(t,k,h_rg,w_rg,state,dwage,type) EMAX_mat[(((((((t)*(T+2) + (k))*RG_SIZE + (h_rg))*RG_SIZE + (w_rg))*STATE_SIZE + (state))*D_WAGE + (dwage))*TYPE_SIZE + (type))]
#else
#define EMAX(t,k,h_rg,w_rg,state,dwage,type) 0
#endif

#ifdef SIMULATION
static const unsigned int RENT_SIM = 1;
static const unsigned int WAGE_SIM = 2;
static const unsigned int TC_SIM = 3;
static const unsigned int FR_SIM = 5;
static const unsigned int MARRIED_SIM = 6;
#endif

// estimation function used inside the optimization process to find the params the find minimum likelihood
// input: array of MAX_PARAM_LEN (157) parameters
// output: likelihood of these params in respect to the individuals' params and the moments

#define RENT_REF_PARAM 1.0
#define WAGE_REF_PARAM 1.0

#ifdef CALC_STDEV
static double estimation(float* params, FILE *fp)
#elif SIMULATION
static double estimation(float* params, unsigned int sim_type, float sim_percent)
#elif REF_PARAM
#undef RENT_REF_PARAM
#undef WAGE_REF_PARAM
#define RENT_REF_PARAM rent_param
#define WAGE_REF_PARAM wage_param
static double estimation(float* params, float rent_param, float wage_param)
#else
static double estimation(float* params)
#endif
{
    const float beta = 0.985;  // discount rate
    const float one_by_beta = 1.0f/beta;
    unsigned short i = 0; // index in the input array v
    unsigned short j = 0; // temp index in output arrays
    // convert input vector v to arguments with meaningful names
    float alfa1[RG_SIZE]; // constant of unemployment in first period for every region [0..6] 
    for (; j < RG_SIZE; ++j, ++i)
    {
        alfa1[j] = expf(params[i]);
    }

    float alfa2[TYPE_SIZE]; // Moving costs by type [7..9]
    for (j = 0; j < TYPE_SIZE; ++j, ++i)
    {
        alfa2[j] = expf(params[i]);
    }

    //Taste for Residential Location
    float teta0[RG_SIZE];   // constant for every region [10..14] - only 5 (at the first and last region =0)
    teta0[0] = 0.0;
    for (j = 1; j < 6; ++j, ++i)
    {
        teta0[j] = params[i]*1000.0f;
    }
    teta0[6] = 0.0;

    float teta1[RG_SIZE];   //  republic 1 for every region [15..19] - only 5 (at the first and last region =0)
    teta1[0] = 0.0;
    for (j = 1; j < 6; ++j, ++i)
    {
        teta1[j] = params[i]*1000.0f;
    }
    teta1[6] = 0.0;

    float teta2[RG_SIZE];   //  republic 2 for every region [20..24] - only 5 (at the first and last region =0)
    teta2[0] = 0.0;
    for (j = 1; j < 6; ++j, ++i)
    {
        teta2[j] = params[i]*1000.0f;
    }
    teta2[6] = 0.0;

    float teta3[RG_SIZE];   //  republic 1 for every region [25..29] - only 5 (at the first and last region =0)
    teta3[0] = 0.0;
    for (j = 1; j < 6; ++j, ++i)
    {
        teta3[j] = params[i]*1000.0f;
    }
    teta3[6] = 0.0;

    // Housing Cost
    float gama0[RG_SIZE]; // constant for every region [30..36]
    for (j = 0; j < RG_SIZE; ++j, ++i)
    {
        gama0[j] = params[i];
    }
    // marital status [37]
    float gama1 = params[i]/10000.0f; ++i;
    // number of children [38]
    float gama2 = params[i]/100.0f; ++i;
    // type1 [39]
    float gama3 = params[i]/10.0f; ++i;
    // type2 [40]
    float gama4 = params[i]/10.0f; ++i;
    // Wage parameters in white collar
    float beta20[RG_SIZE];// constant for every region [41..47]
    for (j = 0; j < RG_SIZE; ++j, ++i)
    {
        beta20[j] = params[i]/10.0f;
    }
    float beta21_1 = params[i]/100.0f; ++i;// schooling [48]
    float beta21_2 = params[i]/100.0f; ++i;// schooling [49]
    float beta21_3 = params[i]/100.0f; ++i;// schooling [50]
    float beta22 = params[i]/1000.0f; ++i;// experience in USSR [51]
    float beta23 = params[i]/1000.0f; ++i;// exp^2  in USSR [52]
    float beta24 = params[i]/100.0f; ++i;// experience in Israel [53]
    float beta25 = params[i]/1000.0f; ++i;// exp^2 in Israel [54]
    float beta26 = params[i]/10.0f; ++i;// age over 40 [55]
    float beta27 = params[i]/10.0f; ++i;// type1 [56]
    float beta28 = params[i]/10.0f; ++i;// type2 [57]
    // Wage parameters in blue collar
    float beta30[RG_SIZE]; // constant for every region[58...64] 
    for (j = 0; j < RG_SIZE; ++j, ++i)
    {
         beta30[j] = params[i]/10.0f;
    }
    float beta31_1 = params[i] / 100.0f; ++i; // schooling[65]
    float beta31_2 = params[i] / 100.0f; ++i; // schooling[66]
    float beta31_3 = params[i] / 100.0f; ++i; //schooling[67]
    float beta32 = beta22; // experience in USSR - as in white collar
    float beta33 = beta23; // exp^2  in USSR - as in white collar
    set_param(beta34)// experience in Israel[68]
    beta34 = beta34 / 100.0f;
    set_param(beta35)// exp^2 in Israel[69]
    beta35 = beta35 / 1000.0f;
    float beta36 = beta26; // age over 40 - as in white collar
    set_param(beta37)// type1[70]
    beta37 = beta37 / 10.0f;
    set_param(beta38)// type2[71]
    beta38 = beta38 / 10.0f;
    // Travelling Costs
    float tc[TC_SIZE]; //[72..74]
    for (j = 0; j < TC_SIZE; ++j, ++i)
    {
        tc[j] = expf(params[i]);
    }
    // tc1 - tc12,tc13,tc45
    // tc2 - tc15,tc14,tc16,tc17,tc25,tc26,tc27,tc35,tc36,tc37,yc46,tc47,tc56,tc57,tc67
    // tc3 - tc23,tc24,tc34
    // tc4 - 0
    // probability of Losing Job, by Type
    set_param_array(aw, STATE_SIZE) // Probability of Losing Job in white color-type 0,1,2[75..77]
    float ab[STATE_SIZE]; // Probability of Losing Job in blue color-type 0 -  types 1 and 2 are the same as in white[78]
    ab[0] = params[i]; ++i;
    ab[1] = aw[1];
    ab[2] = aw[2];
    // Job offer parameters - white collar
    set_param_array(lamda20, RG_SIZE) // constant for every region[79...85]
    float lamda21_1 = params[i]/10.0f; ++i; // schooling[86]
    float lamda21_2 = params[i]/10.0f; ++i; //schooling [87]
    float lamda21_3 = params[i]/10.0f; ++i; //schooling[88]
    set_param(lamda22) // unemployment[89]
    set_param(lamda23) // age at arrival[90]
    lamda23 = lamda23/100.0f;
    set_param(lamda25) // time[91]
    set_param(lamda26) // time^2[92]
    set_param(lamda27) // type1[93]
    set_param(lamda28) // type2[94]
    set_param(psai_w) // t==0 - a new one[95] 
    set_param(psai_b) // t==0 - a new one[96] 

    // Job  offer parameters - blue collar
    set_param_array(lamda30, RG_SIZE) // constant for every region[97...103]
    //float lamda31 = lamda21; // schooling - the same as in white collar
    float lamda32 = lamda22; // unemployment - the same as in white collar
    float lamda33 = lamda23; // age at arrival - the same as in white collar
    set_param(lamda35) // time[104]
    set_param(lamda36) // time^2[105]
    set_param(lamda37) // type1[106]
    set_param(lamda38) // type2[107]
    // type probability
    set_param(type1) // type1 probability[108]
    set_param(type2) // type2 probability[109]
    // Standard deviation of measurement errors
    set_param(error_w) // Standard deviation of measurement errors for wages[110]
    error_w = error_w / 10.0f;
    set_param(error_h) // Standard deviation of measurement errors for houseing[111]
    error_h = expf(error_h / 10.0f);
    set_param(error_c) // Base classification error rate[112]
    error_c = expf(error_c);
    error_c = error_c/(1.0f+error_c);
    // sigma
    set_param_array(sgma, 3) // sigma for white and blue collar and unemployment [113...115]
    sgma[0] = expf(sgma[0]);
    sgma[1] = expf(sgma[1]);
    set_param(wage_error) //punishment for miss the region or occupation in wage [116]
    wage_error = expf(wage_error);      
    set_param(row_w) // serial corellation for white collar [117]
    set_param(row_b) // serial corellation for blue collar  [118]

    set_param_array(R, RG_SIZE) //rent by area [119...125]
    
    set_param_array(psi1, RG_SIZE) //wpmen Married by region[126...132]
    set_param_array(psi2, RG_SIZE) //women education by region[133...139]
    set_param_array(psi3, RG_SIZE) //married*kids by region[140...146]
    set_param_array(psi4, RG_SIZE) //married*women age by region[147...153]
    
    set_param(type1_edu_15) // ??? [154]
    set_param(type1_edu_17) // ??? [155]
    set_param(type2_edu_15) // ??? [156]
    set_param(type2_edu_17) // ??? [157]


    //float PROB_T1=expf(type1)/(1.0f+(expf(type1)+expf(type2)));  
    //float PROB_T2=expf(type2)/(1.0f+(expf(type1)+expf(type2)));
    //float PROB_T0=1.0f-PROB_T1-PROB_T2;
    // P_W_ERROR[0] = 0, P(x<P_W_ERROR[1]) = 10%, P(x<P_W_ERROR[2]) = 30%, P(x<P_W_ERROR[3]) = 50%, P(x<P_W_ERROR[4]) = 70%, P(x<P_W_ERROR[5]) = 90%
    static const float P_W_ERROR[D_WAGE] = {0.0f, -1.281551f, -0.524401f, 0.0f, 0.524401f, 1.281551f};
    // P(x<P_W_ERROR_RNG[1]) = 20%, P(x<P_W_ERROR_RNG[2]) = 40%, P(x<P_W_ERROR_RNG[3]) = 60%, P(x<P_W_ERROR_RNG[4]) = 80%
    static const float P_W_ERROR_RNG[D_WAGE] = {0.0f, -0.841621f, -0.253347f, 0.253347f, 0.841621f, 0.0f};
    // tc1 - tc12,tc13,tc45 // tc0 - tc01,tc02,tc34
    // tc2 - tc15,tc14,tc16,tc17,tc25,tc26,tc27,tc34,tc35,tc36,tc37,yc46,tc47,tc56,tc57,tc67
    // tc3 - tc23,tc24  // tc2 - tc12,tc13,tc23
    // same index - 0
    // note: actual indexes to array are the above minus one
    float travel_cost_arr[RG_SIZE][RG_SIZE];
#ifdef SIMULATION
    
    travel_cost_arr[0][1] = tc[0]*((sim_type == TC_SIM) ? (1.0f - sim_percent) : 1.0f);
    travel_cost_arr[0][2] = tc[0]*((sim_type == TC_SIM) ? (1.0f - sim_percent) : 1.0f);
    travel_cost_arr[3][4] = tc[2]*((sim_type == TC_SIM) ? (1.0f - sim_percent) : 1.0f);
    travel_cost_arr[0][4] = tc[1]*((sim_type == TC_SIM) ? (1.0f - sim_percent) : 1.0f);
    travel_cost_arr[0][3] = tc[1];
    travel_cost_arr[0][5] = tc[1]*((sim_type == TC_SIM) ? (1.0f - sim_percent) : 1.0f);
    travel_cost_arr[0][6] = tc[1];
    travel_cost_arr[1][4] = tc[1]*((sim_type == TC_SIM) ? (1.0f - sim_percent) : 1.0f);
    travel_cost_arr[1][6] = tc[1];
    travel_cost_arr[1][5] = tc[1]*((sim_type == TC_SIM) ? (1.0f - sim_percent) : 1.0f);
    travel_cost_arr[2][4] = tc[1]*((sim_type == TC_SIM) ? (1.0f - sim_percent) : 1.0f);
    travel_cost_arr[2][5] = tc[1]*((sim_type == TC_SIM) ? (1.0f - sim_percent) : 1.0f);
    travel_cost_arr[2][6] = tc[1];
    travel_cost_arr[3][5] = tc[1]*((sim_type == TC_SIM) ? (1.0f - sim_percent) : 1.0f);
    travel_cost_arr[3][6] = tc[1];
    travel_cost_arr[4][5] = tc[1]*((sim_type == TC_SIM) ? (1.0f - sim_percent) : 1.0f);
    travel_cost_arr[4][6] = tc[1];
    travel_cost_arr[5][6] = tc[1];
    travel_cost_arr[1][2] = tc[2]*((sim_type == TC_SIM) ? (1.0f - sim_percent) : 1.0f);
    travel_cost_arr[1][3] = tc[2];
    travel_cost_arr[2][3] = tc[1];

    travel_cost_arr[1][0] = tc[0];
    travel_cost_arr[2][0] = tc[0];
    travel_cost_arr[4][3] = tc[2];
    travel_cost_arr[4][0] = tc[1];
    travel_cost_arr[3][0] = tc[1];
    travel_cost_arr[5][0] = tc[1];
    travel_cost_arr[6][0] = tc[1];
    travel_cost_arr[4][1] = tc[1]*((sim_type == TC_SIM) ? (1.0f - sim_percent) : 1.0f);
    travel_cost_arr[6][1] = tc[1]*((sim_type == TC_SIM) ? (1.0f - sim_percent) : 1.0f);
    travel_cost_arr[5][1] = tc[1]*((sim_type == TC_SIM) ? (1.0f - sim_percent) : 1.0f);
    travel_cost_arr[4][2] = tc[1]*((sim_type == TC_SIM) ? (1.0f - sim_percent) : 1.0f);
    travel_cost_arr[5][2] = tc[1]*((sim_type == TC_SIM) ? (1.0f - sim_percent) : 1.0f);
    travel_cost_arr[6][2] = tc[1]*((sim_type == TC_SIM) ? (1.0f - sim_percent) : 1.0f);
    travel_cost_arr[5][3] = tc[1];
    travel_cost_arr[6][3] = tc[1];
    travel_cost_arr[5][4] = tc[1]*((sim_type == TC_SIM) ? (1.0f - sim_percent) : 1.0f);
    travel_cost_arr[6][4] = tc[1]*((sim_type == TC_SIM) ? (1.0f - sim_percent) : 1.0f);
    travel_cost_arr[6][5] = tc[1]*((sim_type == TC_SIM) ? (1.0f - sim_percent) : 1.0f);
    travel_cost_arr[2][1] = tc[2]*((sim_type == TC_SIM) ? (1.0f - sim_percent) : 1.0f);
    travel_cost_arr[3][1] = tc[2]*((sim_type == TC_SIM) ? (1.0f - sim_percent) : 1.0f);
    travel_cost_arr[3][2] = tc[1]*((sim_type == TC_SIM) ? (1.0f - sim_percent) : 1.0f);

#else

    travel_cost_arr[0][1] = tc[0];
    travel_cost_arr[0][2] = tc[0];
    travel_cost_arr[3][4] = tc[2];
    travel_cost_arr[0][4] = tc[1];
    travel_cost_arr[0][3] = tc[1];
    travel_cost_arr[0][5] = tc[1];
    travel_cost_arr[0][6] = tc[1];
    travel_cost_arr[1][4] = tc[1];
    travel_cost_arr[1][6] = tc[1];
    travel_cost_arr[1][5] = tc[1];
    travel_cost_arr[2][4] = tc[1];
    travel_cost_arr[2][5] = tc[1];
    travel_cost_arr[2][6] = tc[1];
    travel_cost_arr[3][5] = tc[1];
    travel_cost_arr[3][6] = tc[1];
    travel_cost_arr[4][5] = tc[1];
    travel_cost_arr[4][6] = tc[1];
    travel_cost_arr[5][6] = tc[1];
    travel_cost_arr[1][2] = tc[2];
    travel_cost_arr[1][3] = tc[2];
    travel_cost_arr[2][3] = tc[1];
    // fill the other half of the matrix
    for (unsigned short h_rg = 0; h_rg < RG_SIZE; ++h_rg)
    {
        for (unsigned short w_rg = 0; w_rg < h_rg; ++w_rg)
        {
            travel_cost_arr[h_rg][w_rg] = travel_cost_arr[w_rg][h_rg];
        }
        travel_cost_arr[h_rg][h_rg] = 0.0f;
    }
#endif

#if (TRACE_LOAD || SIMULATION)

#ifdef SIMULATION
    if (sim_type == TC_SIM)
#endif
    {
        printf("--------------------------------------------------------------------------------------------------------------------------\n");
        printf("| region |       1       |       2       |       3       |       4       |       5       |       6       |       7       |\n");
        printf("--------------------------------------------------------------------------------------------------------------------------\n");
        for (unsigned short h_rg = 0; h_rg < RG_SIZE; ++h_rg)
        {
            printf("%hu   \t", h_rg + 1);
            for (unsigned short w_rg = 0; w_rg < RG_SIZE; ++w_rg)
            {
                printf("%f\t", travel_cost_arr[h_rg][w_rg]);
            }
            printf("\n");
        }
    }
#endif

    const unsigned short EDU_LEVELS = 3;
    const unsigned short edu_lower[EDU_LEVELS] = {0, 15, 17};
    const unsigned short edu_upper[EDU_LEVELS] = {14, 16, 24};

#ifdef TRACE
    // house distribution table
    unsigned long house_distribution[TYPE_SIZE][RG_SIZE][T]={{{0}}};
    unsigned long house_distribution_count[TYPE_SIZE][T]={{0}};
    unsigned long house_notype_distribution[RG_SIZE][T]={{0}};
    unsigned long house_notype_distribution_count[T]={0};
    unsigned long house_notype_edu_distribution[EDU_LEVELS][RG_SIZE][T]={{{0}}};
    unsigned long house_notype_edu_distribution_count[EDU_LEVELS][T]={{0}};
    // white work region distribution
    unsigned long work_rg_distribution[TYPE_SIZE][RG_SIZE][T]={{{0}}};
    unsigned long work_rg_distribution_count[TYPE_SIZE][T]={{0}};
    unsigned long work_rg_notype_distribution[RG_SIZE][T]={{0}};
    unsigned long work_rg_notype_distribution_count[T]={0};
    // occupation distribution
    unsigned long occ_distribution[TYPE_SIZE][STATE_SIZE][T]={{{0}}};
    unsigned long occ_distribution_count[TYPE_SIZE][T]={{0}};
    unsigned long occ_notype_distribution[STATE_SIZE][T]={{0}};
    unsigned long occ_notype_distribution_count[T]={0};
    unsigned long occ_notype_edu_distribution[EDU_LEVELS][STATE_SIZE][T]={{{0}}};
    unsigned long occ_notype_edu_distribution_count[EDU_LEVELS][T]={{0}};
    // house-work region distribution
    unsigned long house_work_rg_distribution[TYPE_SIZE][RG_SIZE][RG_SIZE]={{{0}}};
    unsigned long house_work_rg_distribution_count[TYPE_SIZE][RG_SIZE]={{0}};
    unsigned long house_work_rg_notype_distribution[RG_SIZE][RG_SIZE]={{0}};
    unsigned long house_work_rg_notype_distribution_count[RG_SIZE]={0};
    // average rent per region
    float rent_rg_sum[TYPE_SIZE][RG_SIZE];
    memset(rent_rg_sum, '\0', sizeof(rent_rg_sum));
    unsigned long rent_rg_count[TYPE_SIZE][RG_SIZE]={{0}};
    float rent_rg_notype_sum[RG_SIZE];
    memset(rent_rg_notype_sum, '\0', sizeof(rent_rg_notype_sum));
    unsigned long rent_rg_notype_count[RG_SIZE]={0};
    // average white wage per region
    float wage_white_rg_sum[TYPE_SIZE][RG_SIZE];
    memset(wage_white_rg_sum, '\0', sizeof(wage_white_rg_sum));
    unsigned long wage_white_rg_count[TYPE_SIZE][RG_SIZE]={{0}};
    float wage_white_notype_rg_sum[RG_SIZE];
    memset(wage_white_notype_rg_sum, '\0', sizeof(wage_white_notype_rg_sum));
    unsigned long wage_white_notype_rg_count[RG_SIZE] = {0};
    // average blue wage per region
    float wage_blue_rg_sum[TYPE_SIZE][RG_SIZE];
    memset(wage_blue_rg_sum, '\0', sizeof(wage_blue_rg_sum));
    unsigned long wage_blue_rg_count[TYPE_SIZE][RG_SIZE]={{0}};
    float wage_blue_rg_notype_sum[RG_SIZE];
    memset(wage_blue_rg_notype_sum, '\0', sizeof(wage_blue_rg_notype_sum));
    unsigned long wage_blue_rg_notype_count[RG_SIZE]={0};
    // average rent
    float rent_sum[TYPE_SIZE];
    memset(rent_sum, '\0', sizeof(rent_sum));
    unsigned long rent_count[TYPE_SIZE]={0};
    float rent_notype_sum = 0.0f;
    unsigned long rent_notype_count = 0;
    // average white wage
    float wage_white_sum[TYPE_SIZE];
    memset(wage_white_sum, '\0', sizeof(wage_white_sum));
    unsigned long wage_white_count[TYPE_SIZE]={0};
    float wage_white_notype_sum = 0.0f;
    unsigned long wage_white_notype_count = 0;
    // average blue wage
    float wage_blue_sum[TYPE_SIZE];
    memset(wage_blue_sum, '\0', sizeof(wage_blue_sum));
    unsigned long wage_blue_count[TYPE_SIZE]={0};
    float wage_blue_notype_sum = 0.0f;
    unsigned long wage_blue_notype_count = 0;
#endif
    unsigned long counter_true = 0;
    unsigned long counter_false = 0;
////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
// The Program Start Here
//////////////////////////////////////////////////////////////////////////////////  
//////////////////////////////////////////////////////////////////////////////////  

    // string likelihood per individual    
    double like_arr[OBSR];
#ifdef TRACE
    float PROB_T0[OBSR];
#endif
    float PROB_T1[OBSR];
    float
        PROB_T2[OBSR];
#ifdef SIMULATION
    double total_max_utility[T];
    memset(total_max_utility, '\0', sizeof(total_max_utility));
    unsigned long total_max_utility_count[T] = {0};
    double total_benefit = 0.0;
#endif

    for (unsigned short I = 0; I < OBSR; ++I)
    {
        // set variables per individual from loaded arrays
        const unsigned short    M = M_arr[I];
        const unsigned short    KIDS = KIDS_arr[I];
        const unsigned short    EXP_U = EXP_U_arr[I];
        const unsigned long     EXP_U_SQ = EXP_U*EXP_U;
        const unsigned char     SCHOOL1 = (SCHOOL_arr[I] < 16) ? 1 : 0;
        const unsigned char     SCHOOL2 = (SCHOOL_arr[I] == 17) ? 1 : 0;
        const unsigned char     SCHOOL3 = (SCHOOL_arr[i] > 17 || SCHOOL_arr[i] == 16) ? 1 : 0;
#ifndef SIMULATION
        const float             WAGE = WAGE_arr[I];
#endif
        const unsigned short    AGE = AGE_arr[I];
        unsigned short  PERIODS = PERIODS_arr[I];
#ifndef SIMULATION
        const unsigned long     RENT_MORT = RENT_MORT_arr[I];
        const unsigned short    D_MORT = D_MORT_arr[I];
#endif
        const unsigned short    REP1 = REP1_arr[I];
        const unsigned short    REP2 = REP2_arr[I];
        const unsigned short    REP3 = REP3_arr[I];
        short int WIFE_EDU_LEVEL = -1;
        float rent[RG_SIZE][TYPE_SIZE];
        float original_rent[RG_SIZE][TYPE_SIZE];
        float prob_nonfired_w[TYPE_SIZE];
        float prob_nonfired_b[TYPE_SIZE];
        float wife[RG_SIZE];
        float moving_cost[TYPE_SIZE];
        float const_lamda_work_2w[TYPE_SIZE];
        float const_lamda_work_2b[TYPE_SIZE];
        float t_const_tmp_w[TYPE_SIZE];
        float t_const_tmp_b[TYPE_SIZE];
        float const_taste[RG_SIZE];
        float rent_for_all_regions[TYPE_SIZE];
#if defined(SIMULATION)
        unsigned short          IND_FILTER = (sim_type == MARRIED_SIM) ? IND_FILTER_arr[I] : 1;
#elif defined(ONLY_MARRIED)
        unsigned short          IND_FILTER = (M_arr[I] != 0);
#endif

    for (unsigned short type = 0; type < TYPE_SIZE; ++type)
    {
        const unsigned short    TYPE1 = (type == 1);
        const unsigned short    TYPE2 = (type == 2);
        const unsigned short    TYPE0 = (type == 0);
        assert(TYPE1 + TYPE2 + TYPE0 == 1);
        if (M != 0 && WIFE_EDU_arr[I] == 99)
        {
            // if education is missing set it to 13 years (average)
             WIFE_EDU_arr[I] = 12;
        }
        if (M == 0 &&  WIFE_EDU_arr[I] != 99)
        {
            // if not married, set to 99
            WIFE_EDU_arr[I] = 99;
        }
        const unsigned short WIFE_EDU = WIFE_EDU_arr[I];
        if (WIFE_EDU != 99)
        {
            for (unsigned int edu_level = 0; edu_level < EDU_LEVELS; ++edu_level)
            {
                if (WIFE_EDU >= edu_lower[edu_level] && WIFE_EDU <= edu_upper[edu_level])
                {
                    WIFE_EDU_LEVEL = edu_level;
                    break;
                }
            }
        }
        const unsigned int wife_edu_level_15 = WIFE_EDU_LEVEL == 1 ? 1 : 0;
        const unsigned int wife_edu_level_17 = WIFE_EDU_LEVEL == 2 ? 1 : 0;
        const float type1_edu = type1 + type1_edu_15*wife_edu_level_15 + type1_edu_17*wife_edu_level_17;
        const float type2_edu = type2 + type2_edu_15*wife_edu_level_15 + type2_edu_17*wife_edu_level_17;

        PROB_T1[I] = expf(type1_edu)/(1.0+(expf(type1_edu)+expf(type2_edu)));
        PROB_T2[I] = expf(type2_edu)/(1.0+(expf(type1_edu)+expf(type2_edu)));
#ifdef TRACE
        PROB_T0[I] = 1.0-PROB_T1[I]-PROB_T2[I];
#endif

        rent_for_all_regions[type] = gama1*M+gama2*KIDS+gama3*TYPE1+gama4*TYPE2; //global rent without gama0 by region

#ifndef WAGE_SELECTION
        moving_cost[type] = TYPE1*alfa2[1]+TYPE2*alfa2[2]+TYPE0*alfa2[0];
#else
        moving_cost[type] = INFINITY; 
#endif
        for (unsigned short rg = 0; rg < RG_SIZE; ++rg)
        {
            const_taste[rg] = teta0[rg]+teta1[rg]*REP1+teta2[rg]*REP2+teta3[rg]*REP3; //taste for housing in a specific region - equation 5 page 13
            // define const taste, move the final calculation into draws
            rent[rg][type] = 6.0f*expf(gama0[rg]+rent_for_all_regions[type])*RENT_REF_PARAM; // full cost of housing rent/mor - equation 6
            wife[rg] = psi1[rg]*M + psi2[rg]*M*WIFE_EDU + psi3[rg]*M*KIDS + psi4[rg]*M*AGE;
#ifdef SIMULATION
            if (sim_type == RENT_SIM)
            {
                if (rg == 4 || rg == 5)
                {
                    rent[rg][type] *= (1.0f - sim_percent);
                }
            }
#endif
            original_rent[rg][type] = rent[rg][type];
        }
#ifndef WAGE_SELECTION
        //probability of not losing your job in white collar - equation 5 page 12
        prob_nonfired_w[type] = 1.0f/(1.0f + expf(TYPE1*aw[1]+TYPE2*aw[2]+TYPE0*aw[0]));
        //probability of not losing your job in blue collar
        prob_nonfired_b[type] = 1.0f/(1.0f + expf(TYPE1*ab[1]+TYPE2*ab[2]+TYPE0*ab[0]));
        const_lamda_work_2w[type] = (lamda21_1*SCHOOL1 + lamda21_2*SCHOOL2 + lamda21_3*SCHOOL3) + 
                                            lamda23*AGE+lamda27*TYPE1+lamda28*TYPE2; //part of the probability of getting job offer in white - page 13
        const_lamda_work_2b[type] = lamda33*AGE+lamda37*TYPE1+lamda38*TYPE2; //part of the probability of getting job offer in blue - page 13
#else
        prob_nonfired_w[type] = 1.0f;
        prob_nonfired_b[type] = 1.0f;
#endif
        t_const_tmp_w[type] = (beta21_1*SCHOOL1 + beta21_2*SCHOOL2 + beta21_3*SCHOOL3) + 
                                            beta22*EXP_U+beta23*EXP_U_SQ+beta27*TYPE1+beta28*TYPE2;  //part of the wage equation  white collar- equation 7 page 14
        t_const_tmp_b[type] = (beta31_1*SCHOOL1 + beta31_2*SCHOOL2 + beta31_3*SCHOOL3) + 
                                            beta32*EXP_U+beta33*EXP_U_SQ+beta37*TYPE1+beta38*TYPE2;  //part of the wage equation  blue collar- equation 7 page 14
        const float ret = 65.0f - (float)AGE + 10.0f;
        const float one_by_beta_pown = powf(one_by_beta, ret);
        const float terminal = (one_by_beta_pown - 1.0f)/(one_by_beta_pown*(one_by_beta - 1.0f));

// in wage selection there is no emax calculation
#ifndef WAGE_SELECTION

        for (unsigned short t = T; t > 0; --t)
        {
            // loop over periods (decdresing)
            const unsigned long t_sq = t*t;
            const unsigned short age40 = (((float)AGE + (float)t/2.0f) > 39.5f);
            //part of the probability of getting job offer in white - page 13 (miss:constant by region +come from unemp)
            const float lamda_work_2w = const_lamda_work_2w[type]+lamda25*t+lamda26*(float)t_sq;
            //part of the probability of getting job offer in blue - page 13(miss:constant by region +come from unemp) 
            const float lamda_work_2b = const_lamda_work_2b[type]+lamda35*t+lamda36*(float)t_sq;
            const float k_const_tmp_w = t_const_tmp_w[type]+beta26*age40;   //part of the wage equation  white collar- equation 7 page 14 (miss:const by region+exp+exp^2)
            const float k_const_tmp_b = t_const_tmp_b[type]+beta36*age40;   //part of the wage equation  blue collar- equation 7 page 14 (miss:const by region+exp+exp^2)

            float prob_work_2w[RG_SIZE];
            float prob_ue_2w[RG_SIZE];
            float prob_work_2b[RG_SIZE];
            float prob_ue_2b[RG_SIZE];

            for (unsigned short rg = 0; rg < RG_SIZE; ++rg)
            {
                float tmp_lamda = lamda_work_2w + lamda20[rg]; // lamda21*SCHOOL+lamda23*AGE+lamda27*TYPE2+lamda28*TYPE3+lamda25*t+lamda26*t_sq+lamda20[rg]
                float tmp_exp = expf(tmp_lamda);
                prob_work_2w[rg] = tmp_exp/(1.0f+tmp_exp); //probability to get job offer in white if come from work
                // lamda21*SCHOOL+lamda23*AGE+lamda27*TYPE2+lamda28*TYPE3+lamda25*t+lamda26*t_sq+lamda20[rg]+lamda22
                tmp_exp = expf(tmp_lamda + lamda22);
                prob_ue_2w[rg] = tmp_exp/(1.0f+tmp_exp); //probability to get job offer in white if come from unemployment

                tmp_lamda = lamda_work_2b+lamda30[rg]; // lamda33*AGE+lamda37*TYPE2+lamda38*TYPE3+lamda35*t+lamda36*t_sq+lamda30[rg]
                tmp_exp = expf(tmp_lamda);
                prob_work_2b[rg] = tmp_exp/(1.0f+tmp_exp); //probability to get job offer in blue if come from work

                // lamda33*AGE+lamda37*TYPE2+lamda38*TYPE3+lamda35*t+lamda36*t_sq+lamda30[rg]+lamda32
                tmp_exp = expf(tmp_lamda + lamda32);
                prob_ue_2b[rg] = tmp_exp/(1.0f+tmp_exp); //probability to get job offer in blue if come from unemployment

                // adjust rent with R
#ifdef SIMULATION
                if (t <  PERIODS_arr[I] - 1)
#else
                if (t <  PERIODS - 1)
#endif
                {
                    rent[rg][type] = rent[rg][type]/(1.0f + R[rg]);
                }
            }

            for (unsigned short k = 0 ; k <= t; ++k)
            {
                // loop over experience
#ifdef PERF_TRACE
                timeval tv = tic();
#endif
                const unsigned long k_sq = k*k;

                //part of the wage equation  white collar- equation 7 page 14 (adding exp and exp^2 still miss:const by region)
                const float rg_const_tmp_w = k_const_tmp_w+beta24*k+beta25*(float)k_sq;
                //part of the wage equation  blue collar- equation 7 page 14 (adding exp and exp^2 still miss:const by region)
                const float rg_const_tmp_b = k_const_tmp_b+beta34*k+beta35*(float)k_sq;

                // initialize to zeros
                float sum_from_ue_max_utility[RG_SIZE];
                float sum_from_b_max_utility[RG_SIZE][D_WAGE];
                float sum_from_w_max_utility[RG_SIZE][RG_SIZE][D_WAGE];
                memset(sum_from_ue_max_utility, '\0', sizeof(sum_from_ue_max_utility));
                memset(sum_from_b_max_utility, '\0', sizeof(sum_from_b_max_utility));
                memset(sum_from_w_max_utility, '\0', sizeof(sum_from_w_max_utility));

                for (unsigned short draw = 0; draw < DRAWS; ++draw)
                {
                    float from_ue_max_utility = -INFINITY;
                    float from_b_max_utility = -INFINITY;
                    float from_w_max_utility = -INFINITY;

                    float choose_ue[RG_SIZE];
                    float choose_w_emax[RG_SIZE][RG_SIZE][D_WAGE];
                    float ue_2b[RG_SIZE];
                    float work_2b[RG_SIZE];
                    float nonfired_2b[RG_SIZE][D_WAGE];
                    float wage_nonfired_2w[RG_SIZE][D_WAGE];
                    float wage_work_2w[RG_SIZE];
                    float wage_ue_2w[RG_SIZE];
                    float taste[RG_SIZE];
                    float wage_b[RG_SIZE];
                    float wage_w[RG_SIZE];

                    for (unsigned short rg = 0; rg < RG_SIZE; ++rg)
                    {
                        float choose_ue_emax;
                        float choose_b_emax[D_WAGE];
                        //float tmp0 = epsilon_r(rg,draw);
                        float tmp1 = epsilon_b(draw,I,t,rg,WHITE);
                        float tmp2 = epsilon_b(draw,I,t,rg,BLUE);
                        float tmp3 = epsilon_b(draw,I,t,rg,UE);
                        taste[rg] = const_taste[rg];// + exp(sgma[3]*tmp0);
                        
                        for (unsigned short dwage = 0; dwage < D_WAGE; ++dwage)
                        {
                            // tmp1_w = beta21*SCHOOL+beta22*EXP_U+beta23*EXP_U_SQ+beta27*TYPE2+beta28*TYPE3+beta26*age40+beta24*k+beta25*k_sq+beta20[rg]  
                            wage_w[dwage] = 6.0f*expf(rg_const_tmp_w + beta20[rg] + sgma[0]*(tmp1 + row_w*P_W_ERROR[dwage]))*WAGE_REF_PARAM;
                            // tmp1_b = beta31*SCHOOL+beta32*EXP_U+beta33*EXP_U_SQ+beta37*TYPE2+beta38*TYPE3+beta36*age40+beta34*k+beta35*k_sq+beta30[rg]
                            wage_b[dwage] = 6.0f*expf(rg_const_tmp_b + beta30[rg] + sgma[1]*(tmp2 + row_b*P_W_ERROR[dwage]))*WAGE_REF_PARAM;           
                            
#ifdef SIMULATION
                            if (sim_type == WAGE_SIM)
                            {
                                if (rg == 4 || rg == 5)
                                {
                                    wage_w[dwage] *= (1.0f + sim_percent);
                                    wage_b[dwage] *= (1.0f + sim_percent);
                                }
                            }
#endif
                            float tmpdw = tmp1 + row_w*P_W_ERROR[dwage];
                            float tmpdb = tmp2 + row_b*P_W_ERROR[dwage];
                            
                            unsigned short D_W_W = get_discrete_index(tmpdw);
                            unsigned short D_W_B = get_discrete_index(tmpdb);
                            
                            wage_nonfired_2w[rg][dwage] = draw_wage(wage_w[dwage], prob_nonfired_w[type]); //equal wage if ind wasn't fired  and -inf if was fired  
                            float wage_nonfired_2b = draw_wage(wage_b[dwage], prob_nonfired_b[type]); //equal wage if ind wasn't fired  and -inf if was fired
                            if (t == T)
                            {
                                choose_ue_emax = 0.0f;//beta*(delta0+delta1*k+delta2*(AGE+t/2.0));
                                choose_b_emax[dwage]= 0.0f;//beta*(delta0+delta1*k+delta2*(AGE+t/2.0)+delta3);
                                for (unsigned short w_rg = 0; w_rg < RG_SIZE; ++w_rg)
                                {
                                    choose_w_emax[rg][w_rg][dwage]= 0.0f;//beta*(delta0+delta1*k+delta2*(AGE+t/2.0)+delta4);
                                }
                            } 
                            else
                            {
                                choose_ue_emax = beta*EMAX(t+1,k,rg,0,UE,0,type);
                                choose_b_emax[dwage] = beta*EMAX(t+1,k+1,rg,0,BLUE,D_W_B,type);

                                for (unsigned short w_rg = 0; w_rg < RG_SIZE; ++w_rg)
                                {
                                    choose_w_emax[rg][w_rg][dwage] = beta*EMAX(t+1,k+1,rg,w_rg,WHITE,D_W_W,type);
                                }
                            }

                            nonfired_2b[rg][dwage] = wage_nonfired_2b + taste[rg] - rent[rg][type] + wife[rg] + choose_b_emax[dwage];
                        } // close dwag                     
                        wage_ue_2w[rg] = draw_wage(wage_w[0], prob_ue_2w[rg]);              //equal wage if ind come fron ue and got an offer and -inf if didn't
                        float wage_ue_2b = draw_wage(wage_b[0], prob_ue_2b[rg]);              //equal wage if ind come fron ue and got an offer and -inf if didn't
                        wage_work_2w[rg] = draw_wage(wage_w[0], prob_work_2w[rg]);          //equal wage if ind come from and got an offer and -inf if didn't
                        float wage_work_2b = draw_wage(wage_b[0], prob_work_2b[rg]);          //equal wage if ind come from and got an offer and -inf if didn't

                        // the equivalent to "wage" when UE is chosen
                        choose_ue[rg] =  taste[rg] - rent[rg][type] + wife[rg] + expf(sgma[2]*tmp3) + choose_ue_emax;
                        float choose_ue_move = choose_ue[rg] - moving_cost[type];

                        float tmp = taste[rg] - rent[rg][type] + wife[rg] + choose_b_emax[0];
                        ue_2b[rg] = wage_ue_2b + tmp;
                        work_2b[rg] = wage_work_2b + tmp;

                        // stay in ue and move housing
                        get_max(from_ue_max_utility, choose_ue_move);
                        // move from ue to blue and move housing
                        get_max(from_ue_max_utility, ue_2b[rg] - moving_cost[type]);
                        // move from blue to ue and move housing
                        get_max(from_b_max_utility, choose_ue_move);
                        // stay in blue and move housing
                        get_max(from_b_max_utility, work_2b[rg] - moving_cost[type]);
                        // move from white to ue and move housing
                        get_max(from_w_max_utility, choose_ue_move);
                        // move from white to blue and move housing
                        get_max(from_w_max_utility, work_2b[rg] - moving_cost[type]);
                    }//close rg

                    float ue_2w[RG_SIZE][RG_SIZE];
                    float work_2w[RG_SIZE][RG_SIZE];
                    float nonfired_2w[RG_SIZE][RG_SIZE][D_WAGE];

                    for (unsigned short h_rg = 0; h_rg < RG_SIZE; ++h_rg)
                    {
                        for (unsigned short w_rg = 0; w_rg < RG_SIZE; ++w_rg)
                        {
                            float tmp = taste[h_rg] - rent[h_rg][type] + wife[h_rg] - travel_cost(h_rg,w_rg) + choose_w_emax[h_rg][w_rg][0];
                            ue_2w[h_rg][w_rg] = wage_ue_2w[w_rg] + tmp;
                            work_2w[h_rg][w_rg] = wage_work_2w[w_rg] + tmp;
                            for (unsigned short dwage = 0; dwage < D_WAGE; ++dwage)
                            {
                                nonfired_2w[h_rg][w_rg][dwage] = 
                                    wage_nonfired_2w[w_rg][dwage] + taste[h_rg] - rent[h_rg][type] + wife[h_rg] - travel_cost(h_rg,w_rg) + choose_w_emax[h_rg][w_rg][dwage];
                            }   
                            // move from ue to white and move housing
                            get_max(from_ue_max_utility, ue_2w[h_rg][w_rg] - moving_cost[type]);
                            // move from blue to white and move housing
                            get_max(from_b_max_utility, work_2w[h_rg][w_rg] - moving_cost[type]);
                            // stay in white in different work region and move housing
                            get_max(from_w_max_utility, work_2w[h_rg][w_rg] - moving_cost[type]);
                        }//close to_w_rg
                    }//close rg

                    /* utility vector
                        -----------------------------------------------------------------------------------------------------------
                        || 0-6       || 7-13    || 14-20    || 21-27    ||28 - 34  || 35 - 41  || 42 - 48 || 49 - 55  || 56 - 62 ||
                        -----------------------------------------------------------------------------------------------------------
                        || UE        || Blue    || White0   ||  White1  || White2  || White3   || White4  || White5   || White6  ||
                        -----------------------------------------------------------------------------------------------------------
                    */

                    for (unsigned short h_rg = 0; h_rg < RG_SIZE; ++h_rg)
                    {
                        float from_b_max_utility_arr[D_WAGE];
                        float from_w_max_utility_arr[RG_SIZE][D_WAGE];  
                        float from_ue_max_utility_h;
                        if(t < T)
                        {
                            float from_b_max_utility_h = from_b_max_utility;
                            float from_w_max_utility_h = from_w_max_utility;
                            from_ue_max_utility_h = from_ue_max_utility;
                            // stay in ue and live in the same region
                            get_max(from_ue_max_utility_h, choose_ue[h_rg]);
                            // move from ue to blue and live in the same region
                            get_max(from_ue_max_utility_h, ue_2b[h_rg]);
                            // move from blue to ue and live in the same region
                            get_max(from_b_max_utility_h, choose_ue[h_rg]);
    
                            for (unsigned short dwage = 0; dwage < D_WAGE; ++dwage)
                            {
                                from_b_max_utility_arr[dwage] = from_b_max_utility_h;
                                // stay in blue and live in the same region
                                get_max(from_b_max_utility_arr[dwage], nonfired_2b[h_rg][dwage]);
                            }
                            // move from white to ue and live in the same region
                            get_max(from_w_max_utility_h, choose_ue[h_rg]);
                            // move from white to blue and live in the same region
                            get_max(from_w_max_utility_h, work_2b[h_rg]);
    
                            for (unsigned short w_rg = 0; w_rg < RG_SIZE; ++w_rg)
                            {
                                // move from ue to white and live in the same region
                                get_max(from_ue_max_utility_h, ue_2w[h_rg][w_rg]);
                                for (unsigned short dwage = 0; dwage < D_WAGE; ++dwage)
                                {
                                    // move from blue to white and live in the same region
                                    get_max(from_b_max_utility_arr[dwage], work_2w[h_rg][w_rg]);
                                }   //close dwage
                            }//close w_rg
    
                            for (unsigned short from_w_rg = 0; from_w_rg < RG_SIZE; ++from_w_rg)
                            {
                                float from_w_max_utility_tmp = -INFINITY;
                                for (unsigned short to_w_rg = 0; to_w_rg < RG_SIZE; ++to_w_rg)
                                {
                                    // stay in white in different work region and live in the same region
                                    get_max(from_w_max_utility_tmp, work_2w[h_rg][to_w_rg]);
                                    for (unsigned short dwage = 0; dwage < D_WAGE; ++dwage)
                                    {
                                        from_w_max_utility_arr[from_w_rg][dwage] = from_w_max_utility_tmp;
                                        // stayed in white in the same work region and move housing
                                        get_max(from_w_max_utility_arr[from_w_rg][dwage], nonfired_2w[to_w_rg][from_w_rg][dwage] - moving_cost[type]);
                                    } //close dwage
                                } // end to_w_rg

                                for (unsigned short dwage = 0; dwage < D_WAGE; ++dwage)
                                {
                                    // stayed in white in the same work region and live in the same region
                                    get_max(from_w_max_utility_arr[from_w_rg][dwage], nonfired_2w[h_rg][from_w_rg][dwage]);

                                    // adding max of unemployment and blue
                                    get_max(from_w_max_utility_arr[from_w_rg][dwage], from_w_max_utility_h);
                                } //close dwage
                            } // end from_w_rg 

                        } //close T<t
                        else
                        {
                            float from_b_max_utility_h = from_b_max_utility + (from_b_max_utility + moving_cost[type])*terminal;
                            float from_w_max_utility_h = from_w_max_utility + (from_w_max_utility + moving_cost[type])*terminal;
                            from_ue_max_utility_h      = from_ue_max_utility + (from_ue_max_utility + moving_cost[type])*terminal;
                            
                            choose_ue[h_rg] *= (1.0f + terminal);
                            ue_2b[h_rg]     *= (1.0f + terminal);
                            work_2b[h_rg]   *= (1.0f + terminal);

                            // stay in ue and live in the same region
                            get_max(from_ue_max_utility_h, choose_ue[h_rg]);
                            // move from ue to blue and live in the same region
                            get_max(from_ue_max_utility_h, ue_2b[h_rg]);
                            // move from blue to ue and live in the same region
                            get_max(from_b_max_utility_h, choose_ue[h_rg]);

                            for (unsigned short dwage = 0; dwage < D_WAGE; ++dwage)
                            {
                                nonfired_2b[h_rg][dwage] *= (1.0f+terminal);
                                from_b_max_utility_arr[dwage] = from_b_max_utility_h;
                                // stay in blue and live in the same region
                                get_max(from_b_max_utility_arr[dwage], nonfired_2b[h_rg][dwage]);
                            }
                            // move from white to ue and live in the same region
                            get_max(from_w_max_utility_h, choose_ue[h_rg]);
                            // move from white to blue and live in the same region
                            get_max(from_w_max_utility_h, work_2b[h_rg]);
    
                            for (unsigned short w_rg = 0; w_rg < RG_SIZE; ++w_rg)
                            {
                                ue_2w[h_rg][w_rg] *= (1.0f + terminal);
                                work_2w[h_rg][w_rg] *= (1.0f + terminal);
                                // move from ue to white and live in the same region
                                get_max(from_ue_max_utility_h, ue_2w[h_rg][w_rg]);
                                for (unsigned short dwage = 0; dwage < D_WAGE; ++dwage)
                                {
                                    // move from blue to white and live in the same region
                                    get_max(from_b_max_utility_arr[dwage], work_2w[h_rg][w_rg]);
                                }   //close dwage
                            }//close w_rg

                            for (unsigned short from_w_rg = 0; from_w_rg < RG_SIZE; ++from_w_rg)
                            {
                                float from_w_max_utility_tmp = -INFINITY;
                                for (unsigned short to_w_rg = 0; to_w_rg < RG_SIZE; ++to_w_rg)
                                {
                                    // stay in white in different work region and live in the same region
                                    get_max(from_w_max_utility_tmp, work_2w[h_rg][to_w_rg]);
                                    for (unsigned short dwage = 0; dwage < D_WAGE; ++dwage)
                                    {
                                        from_w_max_utility_arr[from_w_rg][dwage] = from_w_max_utility_tmp;
                                        // stayed in white in the same work region and move housing
                                        if (h_rg == 0) //since there is a h_rg loop, we should multiply by the terminal value only in the first loop
                                        {
                                            nonfired_2w[to_w_rg][from_w_rg][dwage] *= (1.0f+terminal);
                                        }
                                        get_max(from_w_max_utility_arr[from_w_rg][dwage], nonfired_2w[to_w_rg][from_w_rg][dwage] - moving_cost[type]);
                                    } //close dwage
                                } // end to_w_rg

                                for (unsigned short dwage = 0; dwage < D_WAGE; ++dwage)
                                {
                                    // stayed in white in the same work region and live in the same region
                                    get_max(from_w_max_utility_arr[from_w_rg][dwage], nonfired_2w[h_rg][from_w_rg][dwage]);

                                    // adding max of unemployment and blue
                                    get_max(from_w_max_utility_arr[from_w_rg][dwage], from_w_max_utility_h);
                                } //close dwage
                            } // end from_w_rg

                        } //close T == t
                                        
                        sum_from_ue_max_utility[h_rg] += from_ue_max_utility_h;
                        for (unsigned short dwage = 0; dwage < D_WAGE; ++dwage)
                        {
                            sum_from_b_max_utility[h_rg][dwage] += from_b_max_utility_arr[dwage];
                            for (unsigned short w_rg = 0; w_rg < RG_SIZE; ++ w_rg)
                            {
                                sum_from_w_max_utility[h_rg][w_rg][dwage] += from_w_max_utility_arr[w_rg][dwage];
                            }
                        } //close dwage
                    } // close h_rg
                } // close draws

                for (unsigned short h_rg = 0; h_rg < RG_SIZE; ++h_rg)
                {
                    EMAX(t,k,h_rg,0,UE,0,type) = sum_from_ue_max_utility[h_rg]/(float)DRAWS;
                    for (unsigned short dwage = 0; dwage < D_WAGE; ++dwage)
                    {
                        EMAX(t,k,h_rg,0,BLUE,dwage,type) = sum_from_b_max_utility[h_rg][dwage]/(float)DRAWS;
                        // loop over housing region at t-1
                        for (unsigned short w_rg = 0; w_rg < RG_SIZE; ++w_rg)
                        {
                            // loop over work region at t-1
                            // same value will be inserted to EMAX, regardless of work region at t-1
                            EMAX(t,k,h_rg,w_rg,WHITE,dwage,type) = sum_from_w_max_utility[h_rg][w_rg][dwage]/(float)DRAWS;
                        }   
                    } // close dwage
                }//close h_rg
#ifdef PERF_TRACE
                printf("calculating EMAX for I = %d t = %d and k = %d took: %f seconds (estimated total = %f minutes)\n",
                       I, t, k, toc(tv), toc(tv)*ITER_COUNT/60.0);
#endif
            } // end loop over experience
        } // end loop over periods
#endif // WAGE_SELECTION
    } // enf of loop over types


//////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////                   SOLVING FORWARD                                          /////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef PERF_TRACE
        timeval tv = tic();
#endif

        float p_bar_arr[STATE_VECTOR_SIZE][T];
        short work_rg_arr[T][DRAWS_F];
        short house_rg_arr[T][DRAWS_F];
#ifndef SIMULATION
        short job_arr[T][DRAWS_F];
        unsigned short max_index_arr[T][DRAWS_F];
#endif
        float last_wage[DRAWS_F];
        float last_rent[DRAWS_F];
        memset(p_bar_arr, '\0', sizeof(p_bar_arr));
        unsigned short type = draw_type(PROB_T1[I], PROB_T2[I]);
#ifdef SIMULATION
        const unsigned short draws_f = (unsigned short)((type == 0) ? TYPE_0_OF_1000 : ((type == 1) ? TYPE_1_OF_1000 : TYPE_2_OF_1000));
#else 
        const unsigned short draws_f = DRAWS_F;
#endif // SIMULATION

        for (unsigned short draw = 0; draw < draws_f; ++draw)
        {
#ifdef FULL_TRACE
            printf("%hu %hu %hu %hu ", I, (REP1 ? 1 : (REP2 ? 2 : (REP3 ? 3 : 4))),  type, draw); 
#endif
            unsigned short dwage_b = 0;
            unsigned short dwage_w = 0;
            unsigned short from_state = UE;  
            unsigned short from_h_rg = 0;
            unsigned short from_w_rg = 0;
            unsigned short k = 0;


#ifdef WAGE_SELECTION
            from_state = (unsigned short)((float)STATE_SIZE*(rand()/(RAND_MAX + 1.0f)));
            from_h_rg = (unsigned short)((float)RG_SIZE*(rand()/(RAND_MAX + 1.0f)));
            if (from_state == WHITE)
            {
                from_w_rg = (unsigned short)((float)RG_SIZE*(rand()/(RAND_MAX + 1.0f)));
            }

#endif // WAGE_SELECTION

#ifdef SIMULATION
            if (sim_type != 0)
            {
                PERIODS = T;
            }
#endif

            // winding the rent back, so it would reach the original value in PERIODS
            for (unsigned short rg = 0; rg < RG_SIZE; ++rg)
            {
                rent[rg][type] = original_rent[rg][type];
            }
#ifdef SIMULATION
            for (unsigned short t = 0; t < PERIODS_arr[I] - 1; ++t)
#else
            for (unsigned short t = 0; t < PERIODS - 1; ++t)
#endif
            {
                    for (unsigned short rg = 0; rg < RG_SIZE; ++rg)
                    {
                        rent[rg][type] = rent[rg][type]/(1.0f + R[rg]);
                    }
            }

#ifdef RANDOM_SELECTION
            const short rg = (unsigned short)((float)RG_SIZE*(rand()/(RAND_MAX + 1.0f)));
#endif
            for (unsigned short t = 0; t < PERIODS; ++t)// loop over periods
            {
                short max_index = -1;
                float max_utility = -INFINITY;
                bool w_wage_flag = false;
                bool b_wage_flag = false;
                const unsigned short age40 = (((float)AGE + (float)t/2.0f) > 39.5f);
                const unsigned long k_sq = k*k;
#ifndef WAGE_SELECTION
                const unsigned long t_sq = t*t;
                const float lamda_work_2w = const_lamda_work_2w[type]+lamda25*t+lamda26*(float)t_sq; //part of the probability of getting job offer in white - page 13
                const float lamda_work_2b = const_lamda_work_2b[type]+lamda35*t+lamda36*(float)t_sq; //part of the probability of getting job offer in blue - page 13
#endif
                //part of the wage equation  white collar- equation 7 page 14 (adding exp and exp^2 still miss:const by region)
                const float rg_const_tmp_w = t_const_tmp_w[type]+beta26*age40+beta24*k+beta25*(float)k_sq;
                //part of the wage equation  blue collar- equation 7 page 14 (adding exp and exp^2 still miss:const by region)
                const float rg_const_tmp_b = t_const_tmp_b[type]+beta36*age40+beta34*k+beta35*(float)k_sq;
                
                float work_2b[RG_SIZE];
                float ue_2b[RG_SIZE];
                float nonfired_2b[RG_SIZE];
                float choose_ue[RG_SIZE];
                float choose_w_emax_non_f[RG_SIZE][RG_SIZE];
                float choose_w_emax[RG_SIZE][RG_SIZE];
                float wage_nonfired_2w[RG_SIZE];
                float taste[RG_SIZE];
                float wage_b[RG_SIZE];
                float wage_w[RG_SIZE];
                float wage_b_non_f[RG_SIZE];
                float wage_w_non_f[RG_SIZE];
                float wage_ue_2w[RG_SIZE];
                float wage_work_2w[RG_SIZE];
                unsigned short D_W_W[RG_SIZE];
                unsigned short D_W_B[RG_SIZE];

#ifdef RANDOM_SELECTION
                from_h_rg = rg;
#elif REDUCED_SELECTION
                const short rg = live(I,t);
                if (rg < 0 || rg >= RG_SIZE)
                {
#ifdef FULL_TRACE_INDEX
                    printf("%d ", 999);
#elif FULL_TRACE_WAGE
                    printf("%.3f ",  0.0);
#endif
                    continue;
                }

                from_h_rg = rg;
#else
                for (unsigned short rg = 0; rg < RG_SIZE; ++rg)
#endif
                {
                    // adjust rent with R
                    rent[rg][type] = rent[rg][type]*(1.0f + R[rg]);

                    float prob_work_2w;
                    float prob_ue_2w;
                    float prob_work_2b;
                    float prob_ue_2b;
#ifndef WAGE_SELECTION
                    {
                        const float tmp_lamda = lamda_work_2b + lamda30[rg];        // lamda33*AGE+lamda37*TYPE2+lamda38*TYPE3+lamda35*t+lamda36*t_sq+lamda30[rg]
                        float tmp_exp = expf(tmp_lamda);
                        prob_work_2b = tmp_exp/(1.0f+tmp_exp);          // probability to get job offer in blue if come from work
                        
                        tmp_exp = expf(tmp_lamda + lamda32);            // lamda33*AGE+lamda37*TYPE2+lamda38*TYPE3+lamda35*t+lamda36*t_sq+lamda30[rg]+lamda32
                        prob_ue_2b = tmp_exp/(1.0f+tmp_exp);            // probability to get job offer in blue if come from unemployment

                        if (t == 0)
                        {
                            prob_ue_2b = prob_ue_2b*psai_b;
                        }
                    }
#endif // WAGE_SELECTION
                    taste[rg] = const_taste[rg];// + expf(sgma[3]*tmp0);
                    for (unsigned int w_rg = 0; w_rg < RG_SIZE; ++w_rg)
                    {
                        // tmp1_w = beta21*SCHOOL+beta22*EXP_U+beta23*EXP_U_SQ+beta27*TYPE2+beta28*TYPE3+beta26*age40+beta24*k+beta25*k_sq+beta20[rg]
                        const float tmp_lamda = lamda_work_2w + lamda20[w_rg];  // lamda21*SCHOOL+lamda23*AGE+lamda27*TYPE2+lamda28*TYPE3+lamda25*t+lamda26*t_sq+lamda20[rg]
                        float tmp_exp = expf(tmp_lamda);
                        prob_work_2w = tmp_exp/(1.0f+tmp_exp);          // probability to get job offer in white if come from work
                        
                        tmp_exp = expf(tmp_lamda + lamda22);            // lamda21*SCHOOL+lamda23*AGE+lamda27*TYPE2+lamda28*TYPE3+lamda25*t+lamda26*t_sq+lamda20[rg]+lamda22
                        prob_ue_2w = tmp_exp/(1.0f+tmp_exp);            // probability to get job offer in white if come from unemployment
                        if (t == 0)
                        {
                            prob_ue_2w = prob_ue_2w*psai_w;
                        }
                        const float tmp1 = epsilon_f(draw,I,t,w_rg,WHITE);
                        float tmpdw = tmp1 + row_w*P_W_ERROR[dwage_w];
                        D_W_W[w_rg] = get_discrete_index(tmpdw);
                        wage_w_non_f[w_rg] = 6.0f*expf(rg_const_tmp_w + beta20[w_rg] + sgma[0]*(tmp1+row_w*P_W_ERROR[dwage_w]))*WAGE_REF_PARAM;
                        wage_w[w_rg] = 6.0f*expf(rg_const_tmp_w + beta20[w_rg] + sgma[0]*tmp1)*WAGE_REF_PARAM;
                        wage_ue_2w[w_rg] = draw_wage_f(wage_w[w_rg], prob_ue_2w);       //equal wage if i come fron ue and got an offer and -inf if didn't
                        wage_work_2w[w_rg] = draw_wage_f(wage_w[w_rg], prob_work_2w);   //equal wage if ind come from and got an offer and -inf if didn't
                    }
                    float tmp2 = epsilon_f(draw,I,t,rg,BLUE);
                    wage_b_non_f[rg] = 6.0f*expf(rg_const_tmp_b + beta30[rg] + sgma[1]*(tmp2+row_b*P_W_ERROR[dwage_b]))*WAGE_REF_PARAM;
                    wage_b[rg] = 6.0f*expf(rg_const_tmp_b + beta30[rg] + sgma[1]*tmp2)*WAGE_REF_PARAM;

#ifdef SIMULATION
                    if (sim_type == WAGE_SIM)
                    {
                        if (rg == 4 || rg == 5)
                        {
                            wage_w_non_f[rg] *= (1.0f + sim_percent);
                            wage_w[rg] *= (1.0f + sim_percent);
                            wage_b_non_f[rg] *= (1.0f + sim_percent);
                            wage_b[rg] *= (1.0f + sim_percent);
                            
                        }
                    }
#endif

                    float tmpdb = tmp2 + row_b*P_W_ERROR[dwage_b];
                    D_W_B[rg] = get_discrete_index(tmpdb);

                    // sampling the wage for each of the transitions
                    // note: need to check if it is more efficient to calculate the wage on the fly if needed
                    wage_nonfired_2w[rg] = draw_wage_f(wage_w_non_f[rg], prob_nonfired_w[type]);//equal wage if ind wasn't fired  and -inf if was fired
                    const float wage_nonfired_2b = draw_wage_f(wage_b_non_f[rg], prob_nonfired_b[type]);//equal wage if ind wasn't fired  and -inf if was fired
                    const float wage_ue_2b = draw_wage_f(wage_b[rg], prob_ue_2b);       //equal wage if i come fron ue and got an offer and -inf if didn't
                    const float wage_work_2b = draw_wage_f(wage_b[rg], prob_work_2b);   //equal wage if ind come from and got an offer and -inf if didn't

                    float choose_ue_emax = beta*EMAX(t+1,k,rg,0,UE,0,type);
                    float choose_b_emax_non_f = beta*EMAX(t+1,k+1,rg,0,BLUE,D_W_B[rg],type);
                    float choose_b_emax = beta*EMAX(t+1,k+1,rg,0,BLUE,0,type);

                    for (unsigned short w_rg = 0; w_rg < RG_SIZE; ++w_rg)
                    {
                        choose_w_emax_non_f[rg][w_rg] = beta*EMAX(t+1,k+1,rg,w_rg,WHITE,D_W_W[rg],type);
                        choose_w_emax[rg][w_rg] = beta*EMAX(t+1,k+1,rg,w_rg,WHITE,0,type);
                    }

                    // the equivalent to "wage" when UE is chosen
                    choose_ue[rg] =  taste[rg] - rent[rg][type] + wife[rg] + expf(sgma[2]*epsilon_f(draw,I,t,from_h_rg,UE)) + choose_ue_emax;
                    if (t == 0)
                    {
                        choose_ue[rg] += alfa1[rg];
                    }
                    float tmp = taste[rg] - rent[rg][type] + wife[rg];
                    ue_2b[rg] = wage_ue_2b + tmp + choose_b_emax;
                    work_2b[rg] = wage_work_2b + tmp + choose_b_emax;
                    nonfired_2b[rg] = wage_nonfired_2b + tmp + choose_b_emax_non_f;                  
                } //end rg

//////////////////////////////////////////////////// start maximization    ////////////////////////////////////

                float choices[STATE_VECTOR_SIZE];
                for (unsigned int i = 0; i < STATE_VECTOR_SIZE; ++i)
                {
                    // initialize choices with -inf
                    choices[i] = -INFINITY;
                }

                if (from_state == UE)
                {
#ifndef WAGE_SELECTION
                    if (t > 0)
#endif
                    {
#if !defined(RANDOM_SELECTION) && !defined(REDUCED_SELECTION)
                        for (unsigned short rg = 0; rg < RG_SIZE; ++rg)
#endif
                        {
                            choices[rg] = choose_ue[rg] - moving_cost[type];
                            choices[rg+7] = ue_2b[rg] - moving_cost[type];
                            // stay in ue and move housing
                            get_max_idx(max_utility, max_index, choose_ue[rg] - moving_cost[type], rg);
                            // move from ue to blue and move housing
                            get_max_idx(max_utility, max_index, ue_2b[rg] - moving_cost[type], rg+7);
                        }
                        float ue_2w[RG_SIZE][RG_SIZE];
#if !defined(RANDOM_SELECTION) && !defined(REDUCED_SELECTION)
                        for (unsigned short rg = 0; rg < RG_SIZE; ++rg)
#endif
                        {
                            for (unsigned short w_rg = 0; w_rg < RG_SIZE; ++w_rg)
                            {
                                ue_2w[rg][w_rg] = wage_ue_2w[w_rg] + taste[rg] - rent[rg][type] + wife[rg] - travel_cost(rg,w_rg) + choose_w_emax[rg][w_rg];
                                // move from ue to white and move housing
                                choices[rg+14+7*w_rg] = ue_2w[rg][w_rg] - moving_cost[type];
                                get_max_idx(max_utility, max_index, ue_2w[rg][w_rg] - moving_cost[type], rg+14+7*w_rg);
                            }  //close w_rg
                        }  //close rg

                        // stay in ue and live in the same region
                        choices[from_h_rg] = choose_ue[from_h_rg];
                        get_max_idx(max_utility, max_index, choose_ue[from_h_rg], from_h_rg);
                        // move from ue to blue and live in the same region
                        choices[from_h_rg+7]=ue_2b[from_h_rg];
                        get_max_idx(max_utility, max_index, ue_2b[from_h_rg], from_h_rg+7);
                        for (unsigned short w_rg = 0; w_rg < RG_SIZE; ++w_rg)
                        {
                            // move from ue to white and live in the same region
                            choices[from_h_rg+14+7*w_rg] = ue_2w[from_h_rg][w_rg];
                            get_max_idx(max_utility, max_index, ue_2w[from_h_rg][w_rg], from_h_rg+14+7*w_rg);
                        }

                    }
#ifndef WAGE_SELECTION                  
                    else //t==0
                    {
#if !defined(RANDOM_SELECTION) && !defined(REDUCED_SELECTION)
                        for (unsigned short rg = 0; rg < RG_SIZE; ++rg)
#endif
                        {
                            choices[rg] = choose_ue[rg] - moving_cost[type];
                            choices[rg+7] = ue_2b[rg] - moving_cost[type];
                            // stay in ue and move housing
                            get_max_idx(max_utility, max_index, choose_ue[rg] - moving_cost[type], rg);
                            // move from ue to blue and move housing
                            get_max_idx(max_utility, max_index, ue_2b[rg] - moving_cost[type], rg+7);
                            for (unsigned short w_rg = 0; w_rg < RG_SIZE; ++w_rg)
                            {
                                float ue_2w = wage_ue_2w[w_rg] + taste[rg] - rent[rg][type] + wife[rg] - travel_cost(rg,w_rg) + choose_w_emax[rg][w_rg];
                                // move from ue to white and move housing
                                choices[rg+14+7*w_rg] = ue_2w - moving_cost[type];
                                get_max_idx(max_utility, max_index, ue_2w - moving_cost[type], rg+14+7*w_rg);
                            } // close w_rg
                        }  // close rg
                    }  // close if t>0
#endif // WAGE_SELECTION
                } // close state==UE
                else
                {
                    if (from_state == BLUE)
                    {
#if !defined(RANDOM_SELECTION) && !defined(REDUCED_SELECTION)
                        for (unsigned short rg = 0; rg < RG_SIZE; ++rg)
#endif
                        {
                            choices[rg] = choose_ue[rg] - moving_cost[type];
                            choices[rg+7] = work_2b[rg] - moving_cost[type];
                            // stay in ue and move housing
                            get_max_idx(max_utility, max_index, choose_ue[rg] - moving_cost[type], rg);
                            // move from ue to blue and move housing
                            get_max_idx(max_utility, max_index, work_2b[rg] - moving_cost[type], rg+7);
                        }
                        float work_2w[RG_SIZE][RG_SIZE];
#if !defined(RANDOM_SELECTION) && !defined(REDUCED_SELECTION)
                        for (unsigned short rg = 0; rg < RG_SIZE; ++rg)
#endif
                        {
                            for (unsigned short to_w_rg = 0; to_w_rg < RG_SIZE; ++to_w_rg)
                            {
                                work_2w[rg][to_w_rg] = wage_work_2w[to_w_rg] + taste[rg] - rent[rg][type] + wife[rg] - travel_cost(rg, to_w_rg) + choose_w_emax[rg][to_w_rg];
                                // move from blue to white and move housing
                                choices[rg+14+7*to_w_rg] = work_2w[rg][to_w_rg] - moving_cost[type];
                                get_max_idx(max_utility, max_index, work_2w[rg][to_w_rg] - moving_cost[type], rg+14+7*to_w_rg);
                            }   //to_w_rg

                        }  //rg
                        // move from blue to ue and live in the same region
#ifndef WAGE_SELECTION
                        choices[from_h_rg] = choose_ue[from_h_rg];
                        get_max_idx(max_utility, max_index, choose_ue[from_h_rg], from_h_rg);
#else
                        choices[from_h_rg] = -INFINITY;
#endif
                        // stay in blue and live in the same region
                        choices[from_h_rg+7] = nonfired_2b[from_h_rg];
                        get_max_idx(max_utility, max_index, nonfired_2b[from_h_rg], from_h_rg+7);
                        if(max_index == from_h_rg+7 && t == PERIODS-1)
                        {
                            b_wage_flag = true;
                        }

                        for (unsigned short to_w_rg = 0; to_w_rg < RG_SIZE; ++to_w_rg)
                        {
                            // move from blue to white and live in the same region
                            choices[from_h_rg+14+7*to_w_rg] = work_2w[from_h_rg][to_w_rg]; 
                            get_max_idx(max_utility, max_index, work_2w[from_h_rg][to_w_rg], from_h_rg+14+7*to_w_rg);
                        }
                    }// close if
                    else   // from_state==WHITE
                    {
                        float nonfired_2w[RG_SIZE][RG_SIZE];
#if !defined(RANDOM_SELECTION) && !defined(REDUCED_SELECTION)
                        for (unsigned short rg = 0; rg < RG_SIZE; ++rg)
#endif
                        {
                            choices[rg] = choose_ue[rg]- moving_cost[type];
                            choices[rg+7] = work_2b[rg]- moving_cost[type];
                            // move from white to ue and move housing
                            get_max_idx(max_utility, max_index, choose_ue[rg] - moving_cost[type], rg);
                            // move from white to blue and move housing
                            get_max_idx(max_utility, max_index,  work_2b[rg] - moving_cost[type], rg+7);
                        }//end rg
                        float work_2w[RG_SIZE][RG_SIZE];
#if !defined(RANDOM_SELECTION) && !defined(REDUCED_SELECTION)
                        for (unsigned short rg = 0; rg < RG_SIZE; ++rg)
#endif
                        {
                            for (unsigned short to_w_rg = 0; to_w_rg < RG_SIZE; ++to_w_rg)
                            {
                                float tmp = taste[rg] - rent[rg][type] + wife[rg] - travel_cost(rg, to_w_rg);
                                work_2w[rg][to_w_rg] = wage_work_2w[to_w_rg] + tmp + choose_w_emax[rg][to_w_rg];
                                nonfired_2w[rg][to_w_rg] = wage_nonfired_2w[to_w_rg] + tmp + choose_w_emax_non_f[rg][to_w_rg];
                                // stay in white in different work region and move housing
                                choices[rg+14+7*to_w_rg] = work_2w[rg][to_w_rg] - moving_cost[type];
                                get_max_idx(max_utility, max_index, work_2w[rg][to_w_rg] - moving_cost[type], rg+14+7*to_w_rg);
                            }   //to_w_rg

                        }  //rg
                        // move from white to ue and live in the same region
#ifndef WAGE_SELECTION
                        choices[from_h_rg] = choose_ue[from_h_rg];
                        get_max_idx(max_utility, max_index, choose_ue[from_h_rg], from_h_rg);
#else                   
                        choices[from_h_rg] = -INFINITY;
#endif
                        // move from white to blue and live in the same region
                        choices[from_h_rg+7] = work_2b[from_h_rg];
                        get_max_idx(max_utility, max_index, work_2b[from_h_rg], from_h_rg+7);

                        // FROM_W_RG
                        for (unsigned short to_w_rg = 0; to_w_rg < RG_SIZE; ++to_w_rg)
                        {
                            // stay in white in different work region and live in the same region
                            choices[from_h_rg+14+7*to_w_rg] = work_2w[from_h_rg][to_w_rg];
                            get_max_idx(max_utility, max_index, work_2w[from_h_rg][to_w_rg], from_h_rg+14+7*to_w_rg);
                        } // end to_w_rg
#if !defined(RANDOM_SELECTION) && !defined(REDUCED_SELECTION)
                        for (unsigned short rg = 0; rg < RG_SIZE; ++rg)
#endif
                        {
                            // stayed in white in the same work region and move housing
                            // note: rg stands for the destination housing region, becasue there is no change in work region
                            choices[rg+14+7*from_w_rg] = nonfired_2w[rg][from_w_rg] - moving_cost[type];
                            get_max_idx(max_utility, max_index, nonfired_2w[rg][from_w_rg] - moving_cost[type],
                                    rg+14+7*from_w_rg);
                        } // end rg

                        // stayed in white in the same work region and live in the same region
                        choices[from_h_rg+14+7*from_w_rg] = nonfired_2w[from_h_rg][from_w_rg];               
                        get_max_idx(max_utility, max_index, nonfired_2w[from_h_rg][from_w_rg], from_h_rg+14+7*from_w_rg);
                        if(max_index == from_h_rg+14+7*from_w_rg && t == PERIODS-1)
                        {
                            w_wage_flag = true;
                        }
                    } //end if BLUE
                } //end if UE

                unsigned short tmp_work_rg;
                unsigned short tmp_house_rg;
#ifdef WAGE_SELECTION
                if (max_index == -1)
                {
                    printf("Invalid selection -1\n");
                }
#endif
                {
                    div_t est_house_info = div(max_index,7);
                    tmp_house_rg = (unsigned short)est_house_info.rem;
                    tmp_work_rg = (unsigned short)est_house_info.quot;
                    // work region 0 = unemployment
                    // work region 1 = blue
                    // (work region - 2) = white work region
                }

#ifdef WAGE_SELECTION

                if (tmp_house_rg != from_h_rg)
                {
                    printf("Invalid selection moved housing from %hu to %hu\n", from_h_rg, tmp_house_rg);
                }
                if (tmp_work_rg == 0 && from_state != UE)
                {
                    printf("Invalid selection moved state from %hu to UE\n", from_state);
                }
                if (tmp_work_rg == 1 && from_state != BLUE)
                {
                    printf("Invalid selection moved state from %hu to BLUE\n", from_state);
                }
                if (tmp_work_rg > 1 && from_state != WHITE)
                {
                    printf("Invalid selection moved state from %hu to WHITE\n", from_state);
                }
#endif
 
#if defined(RANDOM_SELECTION) || defined(REDUCED_SELECTION)
                if (rg != tmp_house_rg)
                {
                    printf( "rg = %d tmp_house_rg = %d\n", rg, tmp_house_rg);
                    assert(0);
                }
#endif
#ifdef FULL_TRACE_INDEX
                printf("%hu ", max_index);
#elif FULL_TRACE_WAGE
                {
                    float current_wage = 0.0f;
                    if (tmp_work_rg == 1)
                    {
                        // blue
                        current_wage = ((b_wage_flag == false) ? wage_b[tmp_house_rg] : wage_b_non_f[tmp_house_rg])/6.0f;
                    }
                    else if (tmp_work_rg >= 2)
                    {
                        // white
                        current_wage = ((w_wage_flag == false) ? wage_w[tmp_work_rg-2] : wage_w_non_f[tmp_work_rg-2])/6.0f;
                    }
                    printf("%.3f ",  current_wage);
                }
#elif FULL_TRACE_RENT
                printf("%.3f ", rent[tmp_house_rg][type]/6.0f);   
#endif

#ifdef SIMULATION
                if (sim_type == RENT_SIM)
                {
                    if (tmp_house_rg == 4 || tmp_house_rg == 5)
                    {
                        total_benefit += 6.0f*expf(gama0[tmp_house_rg]+rent_for_all_regions[type])*sim_percent;
                    }
                }
                else if (sim_type == WAGE_SIM)
                {
                    if (tmp_work_rg == 1 && (tmp_house_rg == 4 || tmp_house_rg == 5)) // BLUE
                    {
                        // note: we ignore the difference between non-fired blue and blue
                        total_benefit += 6.0f*expf(rg_const_tmp_b + beta30[tmp_house_rg] + sgma[1]*epsilon_f(draw,I,t,tmp_house_rg,BLUE))*sim_percent; 
                    }
                    else if (tmp_work_rg > 1 && (tmp_work_rg-2 == 4 || tmp_work_rg-2 == 5)) // WHITE
                    {
                        // note: we ignore the difference between non-fired white and white
                        total_benefit += 6.0f*expf(rg_const_tmp_w + beta20[tmp_work_rg-2] + sgma[0]*epsilon_f(draw,I,t,(tmp_work_rg-2),WHITE))*sim_percent;
                    }
                    // else unemployement
                }
                else if (sim_type == TC_SIM)
                {
                    if (tmp_work_rg > 1 && !(tmp_work_rg-2 == 0 || tmp_work_rg-2 == 3 || tmp_work_rg-2 == 6))
                    {
                        // WHITE, not working in 0, 3 or 6
                        total_benefit += sim_percent*travel_cost_arr[tmp_house_rg][tmp_work_rg-2]/(1.0f - sim_percent);
                    }
                }
#endif

                house_rg_arr[t][draw] = tmp_house_rg;
                from_h_rg = tmp_house_rg;
                
                if (tmp_work_rg == 0)
                {
                    // unemployment
#ifndef SIMULATION
                    job_arr[t][draw] = UE;
#endif
                    from_state = UE;
                    dwage_b = 0;
                    dwage_w = 0;
                } 
                else if (tmp_work_rg == 1)
                {
                    // work in blue
#ifndef SIMULATION
                    job_arr[t][draw] = BLUE;
#endif
                    from_state = BLUE;
                    // increase experience
                    ++k;
                    // in blue house_rg equals work_rg
                    work_rg_arr[t][draw] = tmp_house_rg;
                    from_w_rg = tmp_house_rg;
                    dwage_b = D_W_B[tmp_house_rg];
                    dwage_w = 0;
                }
                else
                {
                    // work in white
                    tmp_work_rg -= 2; // now it is a number: 0-6
#ifndef SIMULATION
                    job_arr[t][draw] = WHITE;
#endif
                    from_state = WHITE;
                    // increase experience
                    ++k;
                    work_rg_arr[t][draw] = tmp_work_rg;
                    from_w_rg = tmp_work_rg;
                    dwage_b = 0;
                    dwage_w = D_W_W[tmp_work_rg];
#ifdef TRACE
#if defined(SIMULATION) || defined(ONLY_MARRIED)
                    if (IND_FILTER==1)
#endif
                    {
                        ++work_rg_distribution[type][from_w_rg][t];
                        ++work_rg_distribution_count[type][t];
                        ++work_rg_notype_distribution[from_w_rg][t];
                        ++work_rg_notype_distribution_count[t];
                    
                        ++house_work_rg_distribution_count[type][from_h_rg];
                        ++house_work_rg_distribution[type][from_h_rg][from_w_rg];
                        ++house_work_rg_notype_distribution_count[from_h_rg];
                        ++house_work_rg_notype_distribution[from_h_rg][from_w_rg];
                    }
#endif
                }
#ifndef SIMULATION
                max_index_arr[t][draw] = max_index;
#endif
#ifdef TRACE
#if defined(SIMULATION) || defined(ONLY_MARRIED)
                if (IND_FILTER==1)
#endif
                {
                    ++house_distribution[type][from_h_rg][t];
                    ++house_distribution_count[type][t];
                    ++house_notype_distribution[from_h_rg][t];
                    ++house_notype_distribution_count[t];
                
                    ++occ_distribution[type][from_state][t];
                    ++occ_distribution_count[type][t];
                    ++occ_notype_distribution[from_state][t];
                    ++occ_notype_distribution_count[t];
                    if (WIFE_EDU_LEVEL != -1)
                    {
                        ++house_notype_edu_distribution[WIFE_EDU_LEVEL][from_h_rg][t];
                        ++house_notype_edu_distribution_count[WIFE_EDU_LEVEL][t];
                        ++occ_notype_edu_distribution[WIFE_EDU_LEVEL][from_state][t];
                        ++occ_notype_edu_distribution_count[WIFE_EDU_LEVEL][t];
                    }
                }
#endif
                // calculate last rent and wage
                if (t == PERIODS-1)
                {
#ifndef SIMULATION
                    if (live(I,t) > -1)
#endif
                    {
                        last_rent[draw] = rent[from_h_rg][type]/6.0f;
#ifdef TRACE
#if defined(SIMULATION) || defined(ONLY_MARRIED)
                        if (IND_FILTER==1)
#endif
                        {
                            rent_rg_sum[type][from_h_rg] += last_rent[draw];
                            ++rent_rg_count[type][from_h_rg];
                            rent_rg_notype_sum[from_h_rg] += last_rent[draw];
                            ++rent_rg_notype_count[from_h_rg];
                            rent_sum[type] += last_rent[draw];
                            ++rent_count[type];
                            rent_notype_sum += last_rent[draw];
                            ++rent_notype_count;
                        }
#endif
                    }
#ifndef SIMULATION
                    else
                    {
                        last_rent[draw] = 0.0;
                    }
#endif  
                    if (I!=84 && I!=214 && I!=399 && I!=620 && I!=640)
                    {
                        // all individuals that are not special cases have last_wage calculated here
                        if (from_state == WHITE)
                        {
                            last_wage[draw] = ((w_wage_flag == false) ? wage_w[work_rg_arr[PERIODS-1][draw]] : wage_w_non_f[work_rg_arr[PERIODS-1][draw]])/6.0f;
#ifdef TRACE
#if defined(SIMULATION) || defined(ONLY_MARRIED)
                            if (IND_FILTER==1)
#endif
                            {
                                wage_white_rg_sum[type][from_w_rg] += last_wage[draw];
                                ++wage_white_rg_count[type][from_w_rg];
                                wage_white_sum[type] += last_wage[draw];
                                ++wage_white_count[type];
                                wage_white_notype_rg_sum[from_w_rg] += last_wage[draw];
                                ++wage_white_notype_rg_count[from_w_rg];
                                wage_white_notype_sum += last_wage[draw];
                                ++wage_white_notype_count;
                            }
#endif
                        }
                        else if (from_state == BLUE)
                        {
                            last_wage[draw] = ((b_wage_flag == false) ? wage_b[house_rg_arr[PERIODS-1][draw]] : wage_b_non_f[house_rg_arr[PERIODS-1][draw]])/6.0f;
#ifdef TRACE
#if defined(SIMULATION) || defined(ONLY_MARRIED)
                            if (IND_FILTER==1)
#endif
                            {
                                wage_blue_rg_sum[type][from_h_rg] += last_wage[draw];
                                ++wage_blue_rg_count[type][from_h_rg];
                                wage_blue_rg_notype_sum[from_h_rg] += last_wage[draw];
                                ++wage_blue_rg_notype_count[from_h_rg];
                                wage_blue_sum[type] += last_wage[draw];
                                ++wage_blue_count[type];
                                wage_blue_notype_sum += last_wage[draw];
                                ++wage_blue_notype_count;
                            }
#endif
                        }
                        else
                        {
                            // unemployment
                            last_wage[draw] = 0.0f;
                        }
                    }
                    else
                    {
                        // one of the individuals that were treated in the special cases
                        // doing nothing in the last period
                    }
                }
                else if ((t == 7 && I == 84) || (t == 6 && I == 214) || (t == 1 && I == 399) || (t == 3 && (I == 640 || I == 620)))
                {
                    // not the last period, handle special cases
                    if (from_state == WHITE)
                    {
                        last_wage[draw] = wage_w[from_w_rg]/6.0f;
                    }
                    else if (from_state == BLUE)
                    {
                       last_wage[draw] = wage_b[from_h_rg]/6.0f;
                    }
                    else
                    {
                        last_wage[draw] = 0.0f;
                    }
                }

                {
                    double dvsum = 0.0;
                    double dvtau[STATE_VECTOR_SIZE];

                    for (unsigned short st = 0; st < STATE_VECTOR_SIZE; ++st)
                    {
                        if (choices[st] > -INFINITY)
                        {
                            double dv = choices[st] - max_utility;
                            dvtau[st] = exp(dv/(float)TAU);
                            dvsum += dvtau[st];
                        }
                        else
                        {
                            dvtau[st] = 0.0f;
                        }
                    }
                
                    for (unsigned short st = 0; st < STATE_VECTOR_SIZE; ++st)
                    {
                        if (choices[st] > -INFINITY)
                        {
                            p_bar_arr[st][t] += (float)(dvtau[st]/dvsum);
                        }
                        if (draw == draws_f-1)
                        {
                            p_bar_arr[st][t] = p_bar_arr[st][t]/(float)draws_f;
                        }
                    }
                }

#ifdef SIMULATION
                total_max_utility[t] += max_utility;
                ++total_max_utility_count[t];
#endif

            } // end loop on t
#ifdef FULL_TRACE
            printf("\n");
#endif
        } // end loop of draw  

#ifndef SIMULATION
        like_arr[I] = 0.0;
    
        for (unsigned short draw=0; draw < draws_f; ++draw)
        {
            float p_error[T];
            for (unsigned short t = 0; t < PERIODS ; ++t)// loop over periods
            {
                if (max_index_arr[t][draw] == sample(I,t))
                {
                    // model estimation was correct
                    p_error[t] = error_c + (1.0f - error_c)*p_bar_arr[sample(I,t)][t];
                    ++counter_true;
                }
                else if (sample(I,t) > -1)
                {
                    // model estimation was incorrect
                    p_error[t] = (1.0f - error_c)*p_bar_arr[sample(I,t)][t];    
                    ++counter_false;
                }
                else
                {
                    // missing real information
                    p_error[t] = 0.0f;

                    if (occupation(I,t) == WHITE)
                    {
                        //work in WHITE
                        if (live(I,t) == -1)
                        {
                            // unknown where he live
                            if (work(I,t) == 0)
                            {
                                // unknown where he live & work in 0
                                for (unsigned short rg = 0; rg < RG_SIZE; ++rg)
                                {
                                    if (max_index_arr[t][draw] == 14+rg)
                                    {
                                        p_error[t] += error_c+(1.0f - error_c)*p_bar_arr[14+rg][t];
                                    }
                                    else
                                    {
                                        p_error[t] += (1.0f - error_c)*p_bar_arr[14+rg][t];
                                    }
                                }
                            }
                            else if (work(I,t) == 2)
                            {
                                // unknown where he live & work in 2
                                for (unsigned short rg = 0; rg < RG_SIZE; ++rg)
                                {
                                    if (max_index_arr[t][draw] == 28+rg)
                                    {
                                        p_error[t] += error_c+(1.0f - error_c)*p_bar_arr[28+rg][t];
                                    }
                                    else
                                    {
                                        p_error[t] += (1.0f - error_c)*p_bar_arr[28+rg][t];
                                    }
                                }
                            }
                            else if (work(I,t) == 4)
                            {
                                // unknown where he live & work in 0
                                for (unsigned short rg = 0; rg < RG_SIZE; ++rg)
                                {
                                    if (max_index_arr[t][draw] == 42+rg)
                                    {
                                        p_error[t] += error_c+(1.0f - error_c)*p_bar_arr[42+rg][t];
                                    }
                                    else
                                    {
                                        p_error[t] += (1.0f - error_c)*p_bar_arr[42+rg][t];
                                    }
                                }
                            }
                            else if (work(I,t) == 6)
                            {
                                // unknown where he live & work in 0
                                for (unsigned short rg = 0; rg < RG_SIZE; ++rg)
                                {
                                    if (max_index_arr[t][draw] == 56+rg)
                                    {
                                        p_error[t] += error_c+(1.0f - error_c)*p_bar_arr[56+rg][t];
                                    }
                                    else
                                    {
                                        p_error[t] += (1.0f - error_c)*p_bar_arr[56+rg][t];
                                    }
                                }
                            }
                            else if (work(I,t) == -1)
                            {
                                // work region unknown and unknown where he live
                                for (unsigned short h_rg = 0; h_rg < RG_SIZE; ++h_rg)
                                {
                                    for (unsigned short w_rg = 0; w_rg < RG_SIZE; ++w_rg)
                                    {
                                        if (max_index_arr[t][draw] == 14+h_rg+7*w_rg)
                                        {
                                            p_error[t] += error_c + (1.0f - error_c)*p_bar_arr[14+h_rg+7*w_rg][t];
                                        }
                                        else
                                        {
                                            p_error[t] += (1.0f - error_c)*p_bar_arr[14+h_rg+7*w_rg][t];
                                        }
                                    }
                                }
                            }
                            else
                            {
                                printf("handle_missing_state error - existing invalid work rg (I=%hu t=%hu work=%hu)\n", I, t, work(I,t));
                            }
                        }    
                        else if (work(I,t) == -1)
                        {
                            // work region unknown and known where he live
                            for (unsigned short rg = 0; rg < RG_SIZE; ++rg)
                            {
                                if (job_arr[t][draw] == WHITE && work_rg_arr[t][draw] == rg && live(I,t) == house_rg_arr[t][draw])
                                {
                                    p_error[t] += error_c+(1.0f - error_c)*p_bar_arr[live(I,t) + 7*rg][t];
                                }
                                else
                                {
                                    p_error[t] += (1.0f - error_c)*p_bar_arr[live(I,t)+7*rg][t];
                                }
                            }
                        }
                        else
                        {
                            printf("handle_missing_state error(W): work/live contradict sample (I=%hu t=%hu work=%hu live=%hu sample=%hu)\n",  
                                I, t, work(I,t), live(I,t), sample(I,t));
                        }
                    }
                    else if (occupation(I,t) == UE)
                    {
                        // unemployed
                        if (live(I,t) == -1)
                        {
                            // region of housing is unknown
                            for (unsigned short rg = 0; rg < RG_SIZE; ++rg)
                            {
                                if (max_index_arr[t][draw] == rg)
                                {
                                    p_error[t] += error_c + (1.0f - error_c)*p_bar_arr[rg][t];
                                }
                                else
                                {
                                    p_error[t] += (1.0f - error_c)*p_bar_arr[rg][t];
                                }
                            }
                        }
                        else
                        {
                            printf("handle_missing_state error(UE): live contradict sample (I=%hu t=%hu live=%hu sample=%hu)\n",
                                I, t, live(I,t), sample(I,t));  
                        }
                    }
                    else if (occupation(I,t) == BLUE)
                    {
                        // work in blue
                        if (live(I,t) == -1)
                        {
                            //region of housing is unknown
                            for (unsigned short rg = 0; rg < RG_SIZE; ++rg)
                            {
                                if (max_index_arr[t][draw] == 7+rg)
                                {
                                    p_error[t] += error_c + (1.0f - error_c)*p_bar_arr[7+rg][t];
                                }
                                else
                                {
                                    p_error[t] += (1.0f - error_c)*p_bar_arr[7+rg][t];
                                }
                            }
                        }
                        else
                        {
                            printf("handle_missing_state error(B): live contradict sample (I=%hu t=%hu live=%hu sample=%hu)\n",
                                I, t, live(I,t), sample(I,t));
                        }

                    }
                    else
                    {
                        printf("handle_missing_state error: unknown occupation (I=%hu t=%hu occupation=%hu)\n",
                                I, t, occupation(I,t));
                    }
                }
            } //  close t

            float rent_density = 1.0f;    
            if (D_MORT == 0 && RENT_MORT > 0 && last_rent[draw] > 0.0f)
            {
                float stdnorma = logf((float)RENT_MORT/last_rent[draw]);
                stdnorma = stdnorma/error_h;
                float indexwa = -0.5f*stdnorma*stdnorma;
                rent_density = 0.3989423f*(expf(indexwa)/error_h);
            }
            
            float wage_density = 1.0f;
            if (WAGE > 0.0f && last_wage[draw] > 0.0f)
            {
                float stdnorma;
                if (occupation(I,PERIODS-1) == job_arr[PERIODS-1][draw] && work(I,PERIODS-1) == work_rg_arr[PERIODS-1][draw])
                {
                    stdnorma = logf(WAGE/last_wage[draw]);
                }
                else
                {
                    stdnorma = logf(WAGE/last_wage[draw]*wage_error);
                } 
                stdnorma = stdnorma/error_w;
                float indexwa = -0.5f*stdnorma*stdnorma;
                wage_density = 0.3989423f*(expf(indexwa)/error_w);
            }

            double like = 1.0;
            for (unsigned short t = 0; t < PERIODS; ++t)
            {
                like *= p_error[t];
            }
            
            like *= wage_density*rent_density;
            like_arr[I] += like;

        } // close draws
        like_arr[I] = like_arr[I]/(double)draws_f;

#ifdef PERF_TRACE
        printf("calculating forward for I = %d took: %f seconds (estimated total = %f minutes)\n", I, toc(tv), toc(tv)*OBS/60.0);
#endif

#endif // not SIMULATION

    } // end I

#ifdef INFO
    printf("counter true: %lu counter false: %lu correct: %f%%\n", counter_true, counter_false, 100.0f*(float)counter_true/(float)(counter_true+counter_false));
#endif
    double likelihood = 0.0;
    for (unsigned short I = 0; I < OBSR; ++I)
    {
        double prob = like_arr[I];
        double log_prob;
        if (prob < 1e-300)
        {
#ifdef INFO
            printf("-inf value was calculated for I=%hu\n******************************************************\n", I);
#endif
            log_prob = -110.0;
        }
        else
        {
            log_prob= log(prob);
        }
        likelihood += log_prob;
#ifdef CALC_STDEV
        fprintf(fp, "%f ", log_prob);
#endif
    }
    
#ifdef CALC_STDEV
        fprintf(fp, "%f\n", likelihood);
        fflush(fp);
#endif

#ifdef INFO 
    printf("likelihood = %f\n",likelihood);//should maximize objective function
    fflush(stdout);
#endif

#ifdef SIMULATION
    printf("==============================\n");
    if (sim_type == RENT_SIM)
    {
        printf("rent subsidy simulation\n");
    }
    else if (sim_type == WAGE_SIM)
    {
        printf("wage subsidy simulation\n");
    }
    else if (sim_type == TC_SIM)
    {
        printf("travel cost subsidy simulation\n");
    }
    
    printf("==============================\ntotal cost = %f\n", total_benefit/1000.0);

    printf("------------------------------------------\n");
    printf(" T |   Count   |   Average Max Utility   |\n");
    printf("------------------------------------------\n");
    for (unsigned short t = 0; t < T; ++t)
    {
        printf("%hu\t%lu\t%.0f\n", t, total_max_utility_count[t],  
            (total_max_utility_count[t] > 0) ? total_max_utility[t]/(double)total_max_utility_count[t] : 0.0);
    }
#endif

#ifdef TRACE

    {
        float sum_of_prob_t0 = 0.0;
        float sum_of_prob_t1 = 0.0;
        float sum_of_prob_t2 = 0.0;
        for (unsigned short I = 0; I < OBSR; ++I)
        {
            sum_of_prob_t0 += PROB_T0[I];
            sum_of_prob_t1 += PROB_T1[I];
            sum_of_prob_t2 += PROB_T2[I];
        }
        printf("\ntype 0 probability = %f\t", sum_of_prob_t0/(float)OBSR);
        printf("type 1 probability = %f\t", sum_of_prob_t1/(float)OBSR);
        printf("type 2 probability = %f\t", sum_of_prob_t2/(float)OBSR);
    }

    ////////////////////// Occupation Distribution /////////////////////
    printf("\n\noccupation distribution:\n\n");
    printf("----------------------------------------------------------------\n");
    printf(" T |   count   |     UE      |      WHITE     |       BLUE     |\n");
    printf("----------------------------------------------------------------\n");

#ifdef SIMULATION
    const unsigned short max_T = ((sim_type !=0) ? T : MOMENTS_PERIODS);
#else
    const unsigned short max_T = MOMENTS_PERIODS;
#endif

    for (unsigned short t = 0; t < max_T ; ++t)
    {
        printf("%hu\t%lu\t", t, occ_notype_distribution_count[t]);
        for (unsigned short st = 0; st < STATE_SIZE; ++st)
        {
            printf("%f\t", (float)occ_notype_distribution[st][t]/(float)occ_notype_distribution_count[t]);
        }

        printf("\n");
    }
    printf("\n----------------------------------------------------------\n");
    memset(occ_distribution_count, '\0', sizeof(occ_distribution_count));
    memset(occ_distribution, '\0', sizeof(occ_distribution));
    for (unsigned short I = 0; I < OBSR; ++I)
    {
        unsigned short last_t = PERIODS_arr[I];
        for (unsigned short t = 0; t < last_t; ++t)
        {
            if (occupation(I,t) > -1)
            {
                ++occ_distribution_count[0][t];
                ++occ_distribution[0][occupation(I,t)][t];
            }
        }
    }

    for (unsigned short t = 0; t < MOMENTS_PERIODS ; ++t)
    {
        printf("%hu\t%lu\t", t, occ_distribution_count[0][t]);
        for (unsigned st = 0; st < STATE_SIZE; ++st)
        {
            if (occ_distribution_count[0][t] > 0)
            {
                 printf("%f\t", (float)occ_distribution[0][st][t]/(float)occ_distribution_count[0][t]);
            }
            else
            {
                printf("--------\t");
            }
        }
        printf("\n");
    }

    ////////////////////// Last Whitw Wage  /////////////////////
    printf("\n\naverage white wage in last period:\n\n");
    printf("------------------------------------------------------------------------------------------------------------------------------------\n");
    printf(" ty |      1       |       2       |       3       |       4       |       5       |       6       |       7       |    average    |\n");
    printf("------------------------------------------------------------------------------------------------------------------------------------\n");
    for (unsigned short ty = 0; ty < TYPE_SIZE; ++ty)
    {
        printf("%hu\t", ty);
        unsigned long sum_count = 0;
        float sum_wage = 0.0f;
        for (unsigned short rg = 0; rg  < RG_SIZE; ++rg )
        {
            if (wage_white_rg_count[ty][rg] > 0)
            {
                printf("%.3f\t", wage_white_rg_sum[ty][rg]/(float)wage_white_rg_count[ty][rg]);
                sum_count += wage_white_rg_count[ty][rg];
                sum_wage += wage_white_rg_sum[ty][rg];
            }
            else
            {
                printf("--------\t");
            }
        }
        if (sum_count > 0)
        {
            printf("%.3f\n", sum_wage/(float)sum_count);
        }
        else
        {
            printf("--------\n");
        }
    }
    // average across types per region
    printf("------------------------------------------------------------------------------------------------------------------------------------\n");
    printf("avg\t");
    for (unsigned short rg = 0; rg  < RG_SIZE; ++rg )
    {
        printf("%.3f\t", wage_white_notype_rg_sum[rg]/(float)wage_white_notype_rg_count[rg]);
    }

    // average across all regions
    printf("%.3f\t", wage_white_notype_sum/(float)wage_white_notype_count);

    // real values
    printf("\n------------------------------------------------------------------------------------------------------------------------------------\n");
    memset(wage_white_rg_sum, '\0', sizeof(wage_white_rg_sum));
    memset(wage_white_rg_count, '\0', sizeof(wage_white_rg_count));
    for (unsigned short I = 0; I < OBSR; ++I)
    {
        short st = sample(I, PERIODS_arr[I]-1);
        if (st > 13 && WAGE_arr[I] > 0.0f) // white
        {
            wage_white_rg_sum[0][st/7 - 2] += WAGE_arr[I];
            ++wage_white_rg_count[0][st/7 - 2];
        }
        else if (occupation(I, PERIODS_arr[I]-1) == WHITE && work(I, PERIODS_arr[I]-1) != -1 && WAGE_arr[I] > 0.0f)
        {
            wage_white_rg_sum[0][work(I, PERIODS_arr[I]-1)] += WAGE_arr[I];
            ++wage_white_rg_count[0][work(I, PERIODS_arr[I]-1)];
        }
    }
    {
        printf("real\t");
        unsigned long sum_count = 0;
        float sum_wage = 0.0f;
        for (unsigned short rg = 0; rg  < RG_SIZE; ++rg )
        {
            if (wage_white_rg_count[0][rg] > 0)
            {
                printf("%.3f\t", wage_white_rg_sum[0][rg]/(float)wage_white_rg_count[0][rg]);
                sum_count += wage_white_rg_count[0][rg];
                sum_wage += wage_white_rg_sum[0][rg];
            }
            else
            {
                printf("--------\t");
            }
        }
        if (sum_count > 0)
        {
            printf("%.3f\n", sum_wage/(float)sum_count);
        }
        else
        {
            printf("--------\n");
        }
    }
    
    ////////////////////// Last Blue Wage  /////////////////////
    printf("\n\naverage blue wage in last period:\n\n");
    printf("------------------------------------------------------------------------------------------------------------------------------------\n");
    printf(" ty |      1       |       2       |       3       |       4       |       5       |       6       |       7       |    average    |\n");
    printf("------------------------------------------------------------------------------------------------------------------------------------\n");
    for (unsigned short ty = 0; ty < TYPE_SIZE; ++ty)
    {
        printf("%hu\t", ty);
        unsigned long sum_count = 0;
        float sum_wage = 0.0f;
        for (unsigned short rg = 0; rg  < RG_SIZE; ++rg )
        {
            if (wage_blue_rg_count[ty][rg] > 0)
            {
                printf("%.3f\t", wage_blue_rg_sum[ty][rg]/(float)wage_blue_rg_count[ty][rg]);
                sum_count += wage_blue_rg_count[ty][rg];
                sum_wage += wage_blue_rg_sum[ty][rg];
            }
            else
            {
                printf("--------\t");
            }
        }
        if (sum_count > 0)
        {
            printf("%.3f\n", sum_wage/(float)sum_count);
        }
        else
        {
            printf("--------\n");
        }
    }
    // average across types per region
    printf("------------------------------------------------------------------------------------------------------------------------------------\n");
    printf("avg\t");
    for (unsigned short rg = 0; rg  < RG_SIZE; ++rg )
    {
        printf("%.3f\t",  wage_blue_rg_notype_sum[rg]/(float)wage_blue_rg_notype_count[rg]);
    }
    
    // average across all regions
    printf("%.3f\t", wage_blue_notype_sum/(float)wage_blue_notype_count);

    // real values
    printf("\n------------------------------------------------------------------------------------------------------------------------------------\n");
    memset(wage_blue_rg_sum, '\0', sizeof(wage_blue_rg_sum));
    memset(wage_blue_rg_count, '\0', sizeof(wage_blue_rg_count));
    for (unsigned short I = 0; I < OBSR; ++I)
    {
        short st = sample(I, PERIODS_arr[I]-1);
        if (st > 6 && st < 14 && WAGE_arr[I] > 0.0f) // blue
        {
            wage_blue_rg_sum[0][st-7] += WAGE_arr[I];
            ++wage_blue_rg_count[0][st-7];
        }
        else if (occupation(I, PERIODS_arr[I]-1) == BLUE && live(I, PERIODS_arr[I]-1) != -1 && WAGE_arr[I] > 0.0f)
        {
            wage_blue_rg_sum[0][live(I, PERIODS_arr[I]-1)] += WAGE_arr[I];
            ++wage_blue_rg_count[0][live(I, PERIODS_arr[I]-1)];
        }
    }
    {
        printf("real\t");
        unsigned long sum_count = 0;
        float sum_wage = 0.0;
        for (unsigned short rg = 0; rg  < RG_SIZE; ++rg )
        {
            if (wage_blue_rg_count[0][rg] > 0)
            {
                printf("%.3f\t", wage_blue_rg_sum[0][rg]/(float)wage_blue_rg_count[0][rg]);
                sum_count += wage_blue_rg_count[0][rg];
                sum_wage += wage_blue_rg_sum[0][rg];
            }
            else
            {
                printf("--------\t");
            }
        }
        if (sum_count > 0)
        {
            printf("%.3f\n", sum_wage/(float)sum_count);
        }
        else
        {
            printf("--------\n");
        }
    }

#ifdef RANDOM_SELECTION
    ////////////////////// House Region Distribution //////////////////////
    printf("\n\nhousing region distribution (all types):\n\n");
    printf("-----------------------------------------------------------------------------------------------------------------------------\n");
    printf(" T |   count    |      1       |      2        |      3        |       4        |      5       |       6      |      7      |\n");
    printf("-----------------------------------------------------------------------------------------------------------------------------\n");

    for (unsigned short t = 0; t < max_T; ++t)
    {
        printf("%hu\t%lu\t", t, house_notype_distribution_count[t]);
        for (unsigned rg = 0; rg < RG_SIZE; ++rg)
        {
            printf("%f\t", (float)house_notype_distribution[rg][t]/(float)house_notype_distribution_count[t]);
        }
        printf("\n");
    }

    printf("\n------------------------------------------------------------------------------------------------------------------------------------\n");
    memset(house_distribution_count, '\0', sizeof(house_distribution_count));
    memset(house_distribution, '\0', sizeof(house_distribution));
    for (unsigned short I = 0; I < OBSR; ++I)
    {
        unsigned short last_t = PERIODS_arr[I];
        for (unsigned short t = 0; t < last_t ; ++t)
        {
            if (live(I,t) > -1)
            {
                ++house_distribution_count[0][t];
                ++house_distribution[0][live(I,t)][t];
            }
        }
    }

    for (unsigned short t = 0; t < MOMENTS_PERIODS ; ++t)
    {
        printf("%hu\t%lu\t", t, house_distribution_count[0][t]);
        for (unsigned rg = 0; rg < RG_SIZE; ++rg)
        {
            if (house_distribution_count[0][t] > 0)
            {
                 printf("%f\t", (float)house_distribution[0][rg][t]/(float)house_distribution_count[0][t]);
            }
            else
            {
                printf("--------\t");
            }
        }
        printf("\n");
    }
#endif   

    ////////////////////// Work Region Distribution //////////////////////
    printf("\n\nwork region distribution (all types):\n\n");
    printf("-----------------------------------------------------------------------------------------------------------------------------\n");
    printf(" T |   count    |      1       |      2        |      3        |       4        |      5       |       6      |      7      |\n");
    printf("-----------------------------------------------------------------------------------------------------------------------------\n");
    
    for (unsigned short t = 0; t < max_T; ++t)
    {
        printf("%hu\t%lu\t", t, work_rg_notype_distribution_count[t]);
        for (unsigned short rg = 0; rg < RG_SIZE; ++rg)
        {
            printf("%f\t", (float)work_rg_notype_distribution[rg][t]/(float)work_rg_notype_distribution_count[t]);
        }
        printf("\n");
    }

    printf("\n------------------------------------------------------------------------------------------------------------------------------------\n");
    memset(work_rg_distribution_count, '\0', sizeof(work_rg_distribution_count));
    memset(work_rg_distribution, '\0', sizeof(work_rg_distribution));
    for (unsigned short I = 0; I < OBSR; ++I)
    {
        unsigned short last_t = PERIODS_arr[I];
        for (unsigned short t = 0; t < last_t; ++t)
        {
            if (work(I,t) > -1 && occupation(I,t) > -1 && occupation(I,t) == WHITE)
            {
                ++work_rg_distribution_count[0][t];
                ++work_rg_distribution[0][work(I,t)][t];
            }
        }
    }

    for (unsigned short t = 0; t < MOMENTS_PERIODS ; ++t)
    {
        printf("%hu\t%lu\t", t, work_rg_distribution_count[0][t]);
        for (unsigned rg = 0; rg < RG_SIZE; ++rg)
        {
            if (work_rg_distribution_count[0][t] > 0)
            {
                 printf("%f\t", (float)work_rg_distribution[0][rg][t]/(float)work_rg_distribution_count[0][t]);
            }
            else
            {
                printf("--------\t");
            }
        }
        printf("\n");
    }
    
    printf("\n\nhouse-work region distribution (all types):\n\n");
    printf("----------------------------------------------------------------------------------------------------------------------------\n");
    printf(" rg |   count   |      1       |      2        |      3        |       4        |      5       |       6      |      7      |\n");
    printf("----------------------------------------------------------------------------------------------------------------------------\n");
    
    for (unsigned short h_rg = 0; h_rg < RG_SIZE; ++h_rg)
    {
        printf("%hu\t%lu\t", h_rg+1, house_work_rg_notype_distribution_count[h_rg]);
        for (unsigned short w_rg = 0; w_rg < RG_SIZE; ++w_rg)
        {
            printf("%f\t", (float)house_work_rg_notype_distribution[h_rg][w_rg]/(float)house_work_rg_notype_distribution_count[h_rg]);
        }
        printf("\n");
    }
    printf("\n------------------------------------------------------------------------------------------------------------------------------------\n");
    printf("%hu\t-----\t57.63000\t3.390000\t38.98000\t0.000000\t0.000000\t0.000000\t0.000000\n", 1);
    printf("%hu\t-----\t47.47000\t30.30000\t14.14000\t8.080000\t0.000000\t0.000000\t0.000000\n", 2);
    printf("%hu\t-----\t25.96000\t5.140000\t58.10000\t4.110000\t3.340000\t0.510000\t2.830000\n", 3);
    printf("%hu\t-----\t0.000000\t0.000000\t0.000000\t90.58000\t9.420000\t0.000000\t0.000000\n", 4);
    printf("%hu\t-----\t0.000000\t1.310000\t3.270000\t32.03000\t63.40000\t0.000000\t0.000000\n", 5);
    printf("%hu\t-----\t1.880000\t0.000000\t9.380000\t0.000000\t0.000000\t88.75000\t0.000000\n", 6);
    printf("%hu\t-----\t0.000000\t0.000000\t0.000000\t0.000000\t0.000000\t0.000000\t100.0000\n", 7);

    ////////////////////// House Region Distribution per Education Level//////////////////////
    for (unsigned int edu_level = 0; edu_level < EDU_LEVELS; ++edu_level)
    {
        printf("\n\nhousing region distribution (all types) for education level %hu-%hu:\n\n", edu_lower[edu_level], edu_upper[edu_level]);
        printf("-----------------------------------------------------------------------------------------------------------------------------\n");
        printf(" T |   count    |      1       |      2        |      3        |       4        |      5       |       6      |      7      |\n");
        printf("-----------------------------------------------------------------------------------------------------------------------------\n");

        for (unsigned short t = 0; t < max_T; ++t)
        {
            if (t != 1 && t !=3 && t !=5 && t != 9)
            {
                continue;
            }
            printf("%hu\t%lu\t", t, house_notype_edu_distribution_count[edu_level][t]);
            for (unsigned rg = 0; rg < RG_SIZE; ++rg)
            {
                if (house_notype_edu_distribution_count[edu_level][t] > 0)
                {
                    printf("%f\t", (float)house_notype_edu_distribution[edu_level][rg][t]/(float)house_notype_edu_distribution_count[edu_level][t]);
                }
                else
                {
                    printf("--------\t");
                }

            }
            printf("\n");
        }
        printf("\n------------------------------------------------------------------------------------------------------------------------------------\n");
        memset(house_distribution_count, '\0', sizeof(house_distribution_count));
        memset(house_distribution, '\0', sizeof(house_distribution));
        for (unsigned short I = 0; I < OBSR; ++I)
        {
            unsigned short last_t = PERIODS_arr[I];
            for (unsigned short t = 0; t < last_t ; ++t)
            {
                if (live(I,t) > -1 &&  WIFE_EDU_arr[I] != 99 && WIFE_EDU_arr[I] >= edu_lower[edu_level] && WIFE_EDU_arr[I] <= edu_upper[edu_level])
                {
                    ++house_distribution_count[0][t];
                    ++house_distribution[0][live(I,t)][t];
                }
            }
        }
        for (unsigned short t = 0; t < MOMENTS_PERIODS ; ++t)
        {
            if (t != 1 && t !=3 && t !=5 && t != 9)
            {
                continue;
            }
            printf("%hu\t%lu\t", t, house_distribution_count[0][t]);
            for (unsigned rg = 0; rg < RG_SIZE; ++rg)
            {
                if (house_distribution_count[0][t] > 0)
                {
                    printf("%f\t", (float)house_distribution[0][rg][t]/(float)house_distribution_count[0][t]);
                }
                else
                {
                    printf("--------\t");
                }
            }
            printf("\n");
        }
    }
    
     ////////////////////// Occupation Distribution per Education Level//////////////////////
    for (unsigned int edu_level = 0; edu_level < EDU_LEVELS; ++edu_level)
    {
        printf("\n\noccupation distribution: for education level %hu-%hu:\n\n", edu_lower[edu_level], edu_upper[edu_level]);
        printf("---------------------------------------------------------\n");
        printf(" T |   count   |     UE      |      WHITE     |   BLUE  |\n");
        printf("---------------------------------------------------------\n");
        for (unsigned short t = 0; t < max_T; ++t)
        {
            if (t != 1 && t !=3 && t !=5 && t != 9)
            {
                continue;
            }
            printf("%hu\t%lu\t", t, occ_notype_edu_distribution_count[edu_level][t]);
            for (unsigned st = 0; st < STATE_SIZE; ++st)
            {
                if (occ_notype_edu_distribution_count[edu_level][t] > 0)
                {
                    printf("%f\t", (float)occ_notype_edu_distribution[edu_level][st][t]/(float)occ_notype_edu_distribution_count[edu_level][t]);
                }
                else
                {
                    printf("--------\t");
                }

            }
            printf("\n");
        }

        printf("\n------------------------------------------------------------------------------------------------------------------------------------\n");
        memset(occ_distribution_count, '\0', sizeof(occ_distribution_count));
        memset(occ_distribution, '\0', sizeof(occ_distribution));
        for (unsigned short I = 0; I < OBSR; ++I)
        {
            unsigned short last_t = PERIODS_arr[I];
            for (unsigned short t = 0; t < last_t ; ++t)
            {
                if (occupation(I,t) > -1 && WIFE_EDU_arr[I] != 99 && WIFE_EDU_arr[I] >= edu_lower[edu_level] && WIFE_EDU_arr[I] <= edu_upper[edu_level])
                {
                     const short work_rg = sample(I,t)/7;
                    ++occ_distribution_count[0][t];
                    ++occ_distribution[0][(work_rg == 0) ? UE : ((work_rg == 1) ? BLUE : WHITE)][t];
                }
            }
        }
        for (unsigned short t = 0; t < MOMENTS_PERIODS ; ++t)
        {
            if (t != 1 && t !=3 && t !=5 && t != 9)
            {
                continue;
            }
            printf("%hu\t%lu\t", t, occ_distribution_count[0][t]);
            for (unsigned st = 0; st < STATE_SIZE; ++st)
            {
                if (occ_distribution_count[0][t] > 0)
                {
                    printf("%f\t", (float)occ_distribution[0][st][t]/(float)occ_distribution_count[0][t]);
                }
                else
                {
                    printf("--------\t");
                }
            }
            printf("\n");
        }
    }

#endif

    //percent of correct estimations: -100.0f*(float)counter_true/(float)(counter_true+counter_false);
    return -likelihood;
}

static unsigned short load_dynamic_index(const char* filename, unsigned short* index_arr)
{
    const int COLUMN_NUMBER = 1;
    FILE* fp = fopen(filename,"r");
    if (fp)
    {
        char line[LINE_MAX];
        unsigned short line_idx = 0;
        int col_num = 0;
        unsigned short tmp_param_idx;
        unsigned short dynamic_param_size = MAX_PARAM_LEN;
        while (line_idx <= dynamic_param_size)
        {
            if (fgets(line, LINE_MAX, fp) != 0)
            {
                col_num = sscanf(line, "%hu", &tmp_param_idx);
                if (col_num != COLUMN_NUMBER)
                {
                    printf("wrong format in file %s number of columns: %d\n", filename, col_num);
                    printf("line[%hu]: %s\n", line_idx, line);
                    fclose(fp);
                    return 0;
                }
                else
                {
                    if(line_idx == 0)
                    {
                        if (tmp_param_idx >  dynamic_param_size)
                        {
                            printf("too many (%hu) params in file %s, should be %hu\n", tmp_param_idx, filename, dynamic_param_size);
                            fclose(fp);
                            return 0;
                        }
                        dynamic_param_size = tmp_param_idx;
                    }
                    else
                    {
                        index_arr[line_idx-1] = tmp_param_idx;
#ifdef TRACE_LOAD
                        printf("line[%hu]: %s", line_idx, line);
#endif
                    }
                    ++line_idx;
                }
            } 
            else
            {
                printf("erorr [%s] reading file [%s] at line [%hu]\n", strerror(errno), filename, line_idx);
                fclose(fp);
                return 0;
            }

        }
        fclose(fp);
        return dynamic_param_size;
    } 
    else
    {
        printf("failed to open dynamic index file %s\n", filename);
        return 0;
    }
}
#ifdef WIFE_MODE
static const char* IND_DATA_FILENAME = "wife_ind_data_3.txt";
static const char* MOMENTS_FILENAME = "wife_olim_wide_3.txt";
static const char* WIFE_EDU_FILENAME = "husband_edu.txt"; // this is not a mistake! in "wife mode" the husband is the "wife"
#else // regular mode
static const char* IND_DATA_FILENAME = "ind_data_3.txt";
static const char* MOMENTS_FILENAME = "olim_wide_3.txt";
static const char* WIFE_EDU_FILENAME = "wife_edu.txt";
#endif
static const char* INITIAL_PARAM_FILE = "params.txt";
static const char* PARAM_INDEX_FILE = "params_index.txt";
#ifdef SIMULATION
static const char* FR_FILE_NAME = "fr_params.txt";
static const char* IND_FILTER_FILENAME = "ind_filter.txt";
#endif

// main functions
int main(int argc, char** argv)
{

#ifdef INFO
    timeval tv = tic();
#endif

#ifdef CALC_STDEV
    if (argc < 4)
    {
        fprintf(stderr, "usage: %s <input filename> <output filename> <dx(%%)>\n", argv[0]);
        return -1;
    }
    float dx_percent = (float)atof(argv[3]);
    if (dx_percent == 0.0f || dx_percent > 100.0f || dx_percent < -100.0f)
    {
        fprintf(stderr, "dx(%%) must be between -100 to 100, and cannot be zero\n");
        return -1;
    }
#elif SIMULATION 
    if (argc < 4)
    {
        fprintf(stderr, "usage: %s <input filename> <output filename> <simulation type 0-3,5,6> [percent 0-100]\n", argv[0]);
        return -1;
    }

    unsigned int sim_type = atoi(argv[3]);
    if (sim_type && sim_type != WAGE_SIM && sim_type != RENT_SIM && sim_type != TC_SIM && sim_type != FR_SIM && sim_type != MARRIED_SIM)
    {
        fprintf(stderr, "invalid simulation type value. valid values are:\n0 - none\n1 - rent\n2 - wage\n3 - travel cost\n5 - future interest\n6 - married only\n");
        return -1;
    }

    float sim_percent = 0.0f;
    if (sim_type == WAGE_SIM || sim_type == RENT_SIM || sim_type == TC_SIM)
    {
        if (argc < 5)
        {
            fprintf(stderr, "missing simulation percent value\n");
            return -1;
        }

        sim_percent = (float)atof(argv[4]);
    }

    if (sim_percent < 0.0f || sim_percent > 100.0f)
    {
        fprintf(stderr, "invalid simulation percent value. valid values are: 0 - 100 (inclusive)\n");
        return -1;
    }
    sim_percent /= 100.0f;
#elif REF_PARAM
    if (argc < 5)
    {
        fprintf(stderr, "usage: %s <input filename> <output filename> <rent param> <wage param>\n", argv[0]);
        return -1;
    }
    
    float rent_param = (float)atof(argv[3]);
    float wage_param = (float)atof(argv[4]);
#else // not CALC_STDEV nor SIMULATION
    if (argc < 3)
    {
        fprintf(stderr, "usage: %s <input filename> <output filename> [tag]\n", argv[0]);
        return -1;
    }
#endif

#ifdef SIMULATION
    if (sim_type == MARRIED_SIM)
    {
        if (!load_individuals_filter(IND_FILTER_FILENAME))
        {
            fprintf(stderr, "failed to load individual's filter from file %s - using all individuals\n", IND_FILTER_FILENAME);
            for (unsigned int I = 0; I < OBS; ++I)
            {
                IND_FILTER_arr[I] = 1;
            }
        }
    }
#endif

    if (!load_individuals(IND_DATA_FILENAME))
    {
        fprintf(stderr, "failed to load individual data from file %s\n", IND_DATA_FILENAME);
        return -1;
    }

    if (!load_wife_edu(WIFE_EDU_FILENAME))
    {
        fprintf(stderr, "failed to load wife education data from file %s\n", WIFE_EDU_FILENAME);
        return -1;
    }

    if (!load_moments(MOMENTS_FILENAME))
    {
        fprintf(stderr, "failed to load moments data from file %s\n", MOMENTS_FILENAME);
        return -1;
    }
    init_rand();
    // Starting point
    float x[MAX_PARAM_LEN];
    if (!load_const_params(INITIAL_PARAM_FILE, x))
    {
        fprintf(stderr, "failed to load starting point from file %s\n", INITIAL_PARAM_FILE);
        return -1;
    }
    
    unsigned short dynamic_param_idx[MAX_PARAM_LEN];

    unsigned short dynamic_param_size = load_dynamic_index(PARAM_INDEX_FILE, dynamic_param_idx);
    if (!dynamic_param_size)
    {
        fprintf(stderr, "failed to load dynamic indexes from file %s - using all params\n", PARAM_INDEX_FILE);
    }

    if (!load_dynamic_params(argv[1], x, dynamic_param_idx, dynamic_param_size))
    {
        fprintf(stderr, "failed to load dynamic starting point from file %s\n", argv[1]);
        return -1;
    }

#ifdef INFO
    printf("estimation initialization took: %f seconds\n", toc(tv));
    tv = tic();
#endif

#ifdef CALC_STDEV
    FILE *fp = fopen(argv[2], "w");
    if (fp == 0)
    {
        fprintf(stderr, "failed open output file %s\n", argv[2]);
        return -1;
    }

    for (unsigned short idx = 0; idx < MAX_PARAM_LEN; ++idx)
    {
        float dx = fabs(x[idx])*(dx_percent/100.0f);
        float old_x = x[idx];
        x[idx] += dx;
        double value = estimation(x, fp);
        x[idx] = old_x; 
        printf("done stdev calculation for x[%hu] = %f + %f likelihood = %f\n", idx, x[idx], dx, value);
        fflush(stdout);
    }
#elif SIMULATION

    if (sim_type == FR_SIM && !load_future_intrest(FR_FILE_NAME))
    {
        fprintf(stderr, "failed to load future interest paramsfrom file %s\n", FR_FILE_NAME);
        return -1;
    }
    estimation(x, sim_type, sim_percent);
#elif REF_PARAM
    double value = estimation(x, rent_param, wage_param);
#ifdef INFO
    printf("estimation = %e, took: %f seconds\n", value, toc(tv));
#endif

#else  // not CALC_STDEV nor SIMULATION

    double value = estimation(x);

#ifdef INFO
    printf("estimation = %e, took: %f seconds\n", value, toc(tv));
#endif

    FILE *fp = fopen(argv[2], "w"); 
    if (fp == 0) 
    {
        fprintf(stderr, "failed open output file %s\n", argv[2]);
        return -1;
    }
    
    // write function value to output file
    fprintf(fp, "%e\n", value);
#endif

#if (SIMULATION || REF_PARAM)
#else
    // close output file - either estimation value or stdev output
    fclose(fp);
#endif

    return 0;
}
