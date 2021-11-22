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
#include <algorithm>
#include <iostream>

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

// convert integer state value to string
#define state_to_string(s) ((s) == 0 ? "UE" : ((s) == 1 ? "WHITE" : ((s) == 2 ? "BLUE (FULL)" : ((s) == 3 ? "BLUE (PART)" : "STATE ERROR"))))

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
const unsigned int OBSR         = 412;              // individual 
//const unsigned int OBS          = OBSR*TYPE_SIZE;   // individual multiplay by number of types
const unsigned int DRAWS        = 30;               // draws for emax
const unsigned int RG_SIZE      = 7;                // # of regions
const unsigned int TC_SIZE      = 6;                // travel cost
const unsigned int STATE_SIZE   = 3;                // # of states: white, blue, unemployed
#ifdef TRACE
const unsigned int ALL_STATE_SIZE   = 4;            // # of states: white, blue full, blue part, unemployed
#endif
#ifdef SIMULATION
const unsigned int BASE_DRAWS_F = 1000;
const unsigned int TYPE_0_OF_1000 = BASE_DRAWS_F*0.085;
const unsigned int TYPE_1_OF_1000 = BASE_DRAWS_F*0.605;
const unsigned int TYPE_2_OF_1000 = BASE_DRAWS_F*0.310;
const unsigned int DRAWS_F      = TYPE_1_OF_1000;   // max of: T0 = 85, T1 = 605, T2 = 310
#else
const unsigned int DRAWS_F      = 333;              // draws for forward solving
#endif // SIMULATION
const unsigned int D_WAGE       = 6;                //
const unsigned int TAU          = 50000;
const unsigned int STATE_VECTOR_SIZE = RG_SIZE*RG_SIZE + 3*RG_SIZE; 
// work regions * house regions for white + house regions for blue full + house regions for blue partial + house regions for unemployed
// 7*7 + 7 + 7 + 7 = 70

// random draws
float randn_b_arr[DRAWS][OBSR][T][RG_SIZE][STATE_SIZE][TYPE_SIZE];
#define epsilon_b(draw, I, t, rg, state, type) randn_b_arr[(draw)][(I)][(t)][(rg)][(state)][(type)]

#ifdef SIMULATION
float* randn_f_arr = new float[DRAWS_F*OBSR*T*RG_SIZE*STATE_SIZE*TYPE_SIZE];
#define epsilon_f(draw, I, t, rg, state, type) randn_f_arr[(((((draw)*OBSR + (I))*T + (t))*RG_SIZE + (rg))*STATE_SIZE + (state))*TYPE_SIZE + (type)]
#else
float randn_f_arr[DRAWS_F][OBSR][T][RG_SIZE][STATE_SIZE][TYPE_SIZE];
#define epsilon_f(draw, I, t, rg, state, type) randn_f_arr[(draw)][(I)][(t)][(rg)][(state)][(type)]
#endif // SIMULATION

#define draw_wage(wage,prob) (rand01() < (prob)) ? (wage) : -INFINITY
//TODO what value to put in "full" in case of no work?
inline float draw_blue_wage(float wage, float prob_full, float prob_part, bool& full, float part_wage_factor)
{
    float p = rand01();
    full = true;
    if (p < prob_full)
    {
        full = true;
        return wage;
    }
    else if (p < prob_full + prob_part)
    {
        full = false;
        return wage*part_wage_factor;
    }
    return -INFINITY;
}

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
    
    for (unsigned short i = 0; i < OBSR; ++i)
    {
        for (unsigned short t = 0; t < T; ++t)
        {
            for (unsigned short rg = 0; rg < RG_SIZE; ++rg)
            {
                for (unsigned short s = 0; s < STATE_SIZE; ++s)
                {
                    for (unsigned short type = 0; type < TYPE_SIZE; ++type)
                    {
                        for (unsigned short d = 0; d < DRAWS; ++d)
                        {
                            epsilon_b(d, i, t, rg, s, type) = randn01();
                        }
                        for (unsigned short d = 0; d < DRAWS_F; ++d)
                        {
                            epsilon_f(d, i, t, rg, s, type) = randn01();
                        }
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
const float ITER_COUNT = (float)(OBSR*T*(T+1))/2.0f;
#endif

//////////////////////// profiling - end

// global variables loaded from individual's data file
unsigned short M_arr[OBSR];         // married / single
unsigned short KIDS_arr[OBSR];      //# of children
unsigned short EXP_U_arr[OBSR];     // experience in former USSR
unsigned short SCHOOL_arr[OBSR];    // years of education
float          WAGE_arr[OBSR];      // wage of ind
unsigned short AGE_arr[OBSR];       // age at arrival
unsigned short PERIODS_arr[OBSR];   // number of periods (<T)
unsigned long  RENT_MORT_arr[OBSR]; // cost of housing
unsigned short D_MORT_arr[OBSR];    // dummy
unsigned short REP1_arr[OBSR];      // republic 1 at USSR
unsigned short REP2_arr[OBSR];      // republic 2 at USSR
unsigned short REP3_arr[OBSR];      // republic 3 at USSR
unsigned short TYPE1_arr[OBSR];     // dummy that get 1 of individual is type 1
unsigned short TYPE2_arr[OBSR];     // dummy that get 1 of individual is type 2 
// type 0 is if not 1 or 2
unsigned short HUSBAND_EDU_arr[OBSR];  // husband education array
float USSR_ENG_EXP_arr[OBSR];  // USSR engineer experience array

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
        while (I < OBSR)
        {
            if (fgets(line, LINE_MAX, fp) != 0)
            {
                // formatting the line according to:
                // M KIDS EXP_U SCHOOL WAGE AGE PERIODS RENT_MORT D_MORT REP1 REP2 REP3 TYPE1 TYPE2
                col_num = sscanf(line, "%hu%hu%hu%hu%f%hu%hu%lu%hu%hu%hu%hu%hu%hu",
                                 &(M_arr[I]), &(KIDS_arr[I]), &(EXP_U_arr[I]), &(SCHOOL_arr[I]),
                                 &(WAGE_arr[I]), &(AGE_arr[I]), &(PERIODS_arr[I]), &(RENT_MORT_arr[I]), &(D_MORT_arr[I]), &(REP1_arr[I]), &(REP2_arr[I]),
                                 &(REP3_arr[I]), &(TYPE1_arr[I]), &TYPE2_arr[I]);
                if (col_num != COLUMN_NUMBER)
                {
                    printf("wrong format in file %s number of columns: %d \n", filename, col_num);
                    printf("line[%hu]: %s\n", I, line);
                    printf("format: M\tKIDS\tEXP_U\tSCHOOL\tWAGE\tAGE\tPERIODS\tRENT_MO\tD_MO\tREP1\tREP2\tREP3\tTYPE1\tTYPE2\n");
                    fclose(fp);
                    return false;
                } 
                else
                {
#ifdef TRACE_LOAD
                    printf("line[%hu]: %hu\t%hu\t%hu\t%hu\t%f\t%hu\t%hu\t%lu\t%hu\t%hu\t%hu\t%hu\t%hu\t%hu\n", I, M_arr[I], KIDS_arr[I],
                         EXP_U_arr[I], SCHOOL_arr[I],
                         WAGE_arr[I], AGE_arr[I], PERIODS_arr[I], RENT_MORT_arr[I], D_MORT_arr[I], REP1_arr[I], REP2_arr[I],
                         REP3_arr[I], TYPE1_arr[I], TYPE2_arr[I]);
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
unsigned short IND_FILTER_arr[OBSR];  // filter out individuals for summary

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
        while (I < OBSR)
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
static bool load_ussr_eng_exp(const char* filename)
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
        while (I < OBSR)
        {
            if (fgets(line, LINE_MAX, fp) != 0)
            {
                col_num = sscanf(line, "%f", &(USSR_ENG_EXP_arr[I]));
                if (col_num != COLUMN_NUMBER)
                {
                    printf("wrong format in file %s number of columns: %d \n", filename, col_num);
                    printf("line[%hu]: %s\n", I, line);
                    printf("format: USSR_ENG_EXP\n");
                    fclose(fp);
                    return false;
                } 
                else
                {
#ifdef TRACE_LOAD
                    printf("line[%hu]: %f\n", I, USSR_ENG_EXP_arr[I]);
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
        printf("failed to open USSR engineer experience file %s\n", filename);
        return false;
    }
}

// loading from file to global variables
static bool load_husband_edu(const char* filename)
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
        while (I < OBSR)
        {
            if (fgets(line, LINE_MAX, fp) != 0)
            {
                col_num = sscanf(line, "%hu", &(HUSBAND_EDU_arr[I]));
                if (col_num != COLUMN_NUMBER)
                {
                    printf("wrong format in file %s number of columns: %d \n", filename, col_num);
                    printf("line[%hu]: %s\n", I, line);
                    printf("format: HUSBAND_EDU\n");
                    fclose(fp);
                    return false;
                } 
                else
                {
#ifdef TRACE_LOAD
                    printf("line[%hu]: %hu\n", I, HUSBAND_EDU_arr[I]);
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
        printf("failed to open husband education file %s\n", filename);
        return false;
    }
}

const unsigned short UE    = 0;
const unsigned short BLUE  = 2;
const unsigned short WHITE = 1;

const unsigned short FULL  = 0;
const unsigned short PART  = 1;

const unsigned short MOMENTS_PERIODS = 12;

short occupation_arr[OBSR][MOMENTS_PERIODS]; // real occupation
short live_arr[OBSR][MOMENTS_PERIODS];       // real housing region
short work_arr[OBSR][MOMENTS_PERIODS];       // real work region
short sample_arr[OBSR][MOMENTS_PERIODS];     // real state 0-69
#define occupation(i,j) occupation_arr[(i)][(j)]
#define live(i,j) live_arr[(i)][(j)]
#define work(i,j) work_arr[(i)][(j)]
#define sample(i,j) sample_arr[(i)][(j)]

static bool load_moments(const char* filename)
{
    const int COLUMN_NUMBER = 49; 
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
        while (I < OBSR)
        {
            if (fgets(line, LINE_MAX, fp) != 0)
            {
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
                    for (unsigned short t = 0; t < MOMENTS_PERIODS; ++t)
                    {
                        // fixing data in input file
                        // /////////////////////////
                        // moving regions from 1-7 to 0-6
                        if (live(I,t) != -1)
                        {
                            --live(I,t);
                        }
                        if (work(I,t) != -1)
                        {
                            --work(I,t);
                        }

                        if (sample(I,t) == -1)
                        {
                            if (occupation(I,t) == UE && live(I,t) != -1)
                            {
                                sample(I,t) = live(I,t);
#ifdef TRACE_LOAD
                                printf("line[%hu]: missing state in period = %hu, added = %d\n", I, t, sample(I,t));
#endif
                            }
                            else if (occupation(I,t) == BLUE && live(I,t) != -1)
                            {
                                // TODO: assuming blue in full time
                                sample(I,t) = live(I,t) + 7;
#ifdef TRACE_LOAD
                                printf("line[%hu]: missing state in period = %hu, added = %d\n", I, t, sample(I,t));
#endif
                            }
                            else if (occupation(I,t) == WHITE && live(I,t) != -1 && work(I,t) != -1)
                            {
                                sample(I,t) = live(I,t) + (work(I,t) + 3)*7;
#ifdef TRACE_LOAD
                                printf("line[%hu]: missing state in period = %hu, added = %d\n", I, t, sample(I,t));
#endif
                            }
                        }

                        // verify integrity of input files
                        bool inconsistent = false;
                        if (occupation(I,t) < -1 || occupation(I,t) > 2)
                        {
                            printf("line[%hu]: invalid occupation = %d\n", I, occupation(I,t));
                            inconsistent = true;
                        }
                        if (live(I,t) < -1 || live(I,t) > 6)
                        {
                            printf("line[%hu]: invalid housing region = %d\n", I, live(I,t));
                            inconsistent = true;
                        }
                        if (work(I,t) < -1 || work(I,t) > 6)
                        {
                            printf("line[%hu]: invalid work region = %d\n", I, work(I,t));
                            inconsistent = true;
                        }
                        if (sample(I,t) < -1 || sample(I,t) > 69)
                        {
                            printf("line[%hu]: invalid status = %d\n", I, sample(I,t));
                            inconsistent = true;
                        }
                        short st = sample(I,t);
                        if (st != -1)
                        {
                            // if status is known occupation and house region must be known
                            if (occupation(I,t) == -1 || live(I,t) == -1)
                            {
                                printf("line[%hu]: occupation or housing region data missing for period = %hu\n", I, t);
                                inconsistent = true;
                            }
                            
                            div_t house_info = div(st,7);
                            short house_rg = (short)house_info.rem;
                            short work_rg = (short)house_info.quot;
                            // work region 0 = unemployment
                            // work region 1 = blue, full time
                            // work region 2 = blue, part time
                            // (work region - 3) = white work region
                            if (live(I,t) != -1 && live(I,t) != house_rg)
                            {
                                printf("line[%hu]: housing region inconsistent for period = %hu, %d != %d\n", I, t, house_rg, live(I,t));
                                inconsistent = true;
                            }
                            if (occupation(I,t) != -1 && ((work_rg == 0 && occupation(I,t) != UE) ||
                                ((work_rg == 1 || work_rg == 2) && occupation(I,t) != BLUE) ||
                                (work_rg > 2 && occupation(I,t) != WHITE)))
                            {
                                printf("line[%hu]: occupation inconsistent for period = %hu, %hu != %d\n", I, t,
                                        (work_rg == 0 ? UE : (work_rg == 1 ? BLUE : WHITE)), occupation(I,t));
                                inconsistent = true;
                            }
                            if (work(I,t) != -1 && work_rg > 2 && (work_rg - 3) != work(I,t))
                            {
                                printf("line[%hu]: work region inconsistent for period = %hu, %d != %d\n", I, t, work_rg - 3, work(I,t));
                                inconsistent = true;
                            }
                            if (inconsistent)
                            {
                                printf("Occupation: %d Computed Occupation: %d\n", occupation(I,t), (work_rg == 0 ? UE : (work_rg == 1 || work_rg == 2 ? BLUE : WHITE)));
                                printf("Housing Region: %d Computed Housing Region: %d\n", live(I,t), house_rg);
                                printf("Work Region: %d Computed Work Region: %d\n", work(I,t), (work_rg > 2) ? work_rg - 3 : -1);
                            }
                        }

                        if (inconsistent)
                        {
                            printf("Occupation: %d\n", occupation(I,t));
                            printf("Housing Region: %d\n", live(I,t));
                            printf("Work Region: %d\n", work(I,t));
                        }
                    }
                    ++I;
                }
            } 
            else
            {
                if (errno == 0)
                {
                    printf("invalid number of lines [%hu] in file [%s]\n", I, filename);
                }
                else
                {
                    printf("erorr [%s] reading file [%s] at line [%hu]\n", strerror(errno), filename, I);
                }
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

const unsigned short int MAX_PARAM_LEN = 156;   //# of parameters
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

void print_choices(float* choices)
{
    printf("\n");
    printf("\n");
    printf("---------------------------------------------------------------------------------------------------------------------------------------------\n");
    printf("Utility Vector                                                                                                                             ||\n");
    printf("---------------------------------------------------------------------------------------------------------------------------------------------\n");
    printf("|| 0-6       || 7-13       || 14-20      || 21-27      || 28 - 34    || 35 - 41    || 42 - 48    || 49 - 55    || 56 - 62    || 63 - 69    ||\n");
    printf("---------------------------------------------------------------------------------------------------------------------------------------------\n");
    printf("|| UE        || Blue, full || Blue, part || White0     || White1     || White2     || White3     || White4     || White5     || White6     ||\n");
    printf("---------------------------------------------------------------------------------------------------------------------------------------------\n");

    printf("Unemployment Utility         :");
    for (unsigned short i = 0; i < RG_SIZE; ++i)
    {
        printf(" %6.0f (%hu), ", choices[i], i);
    }
    printf("\n");
    printf("Blue Utility (full time)     :  ");
    for (unsigned short i = 0; i < RG_SIZE; ++i)
    {
        printf("%10.0f (%hu), ", choices[i+7], (unsigned short)(i+7));
    }
    printf("\n");
    printf("Blue Utility (part time)     :");
    for (unsigned short i = 0; i < RG_SIZE; ++i)
    {
        printf("%10.0f (%hu), ", choices[i+14], (unsigned short)(i+7));
    }
    for (unsigned short i = 0; i < RG_SIZE; ++i)
    {
        printf("\n");
        printf("White Utility (work region %hu):", i);
        for (unsigned short j = 0; j < RG_SIZE; ++j)
        {
            printf("%10.0f (%hu), ", choices[(i+3)*7+j], (unsigned short)((i+3)*7+j));
        }
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
// input: array of MAX_PARAM_LEN (156) parameters
// output: likelihood of these params in respect to the individuals' params and the moments

#define RENT_REF_PARAM 1.0
#define WAGE_REF_PARAM 1.0

/* convert k index to experience
index   experience (k)
----------------------
0       0.0
1       0.5
2       1.0
3       1.5
4       2.0
5       2.5
6       3.0
6		3.5
7       4.0
7		4.5
8       5.0
8		5.5
9       6.0
9		6.5
10      7.0
10		7.5
11      8.0
11		8.5
12      9.0
12		9.5
12      10.0
12		10.5
12      11.0
12		11.5
12      12.0
12		12.5
*/
const unsigned short MAXIMUM_K_INDEX=12;
#define index_to_k(I) (((I) >= 6) ? (float)(I) - 3.0 : (float)(I)/2.0)
#define k_to_index(K) (((K) >= 9.0) ? MAXIMUM_K_INDEX : (((K) >= 3.0) ? (int)(K) + 3 : (int)((K)*2.0)))

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
    const float one_by_beta = 1.0/beta;
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
    // regions 1,5,6 must have same value
    beta20[5] = beta20[1];
    beta20[6] = beta20[1];
    float beta21_1 = params[i]/100.0f; ++i;// schooling [48]
    float beta21_3 = params[i]/100.0f; ++i;// schooling [49]
    float beta21_2 = beta21_3;
    float beta22 = params[i]/1000.0f; ++i;// experience in USSR [50]
    float beta23 = params[i]/1000.0f; ++i;// exp^2  in USSR [51]
    float beta24 = params[i]/100.0f; ++i;// experience in Israel [52]
    float beta25 = params[i]/1000.0f; ++i;// exp^2 in Israel [53]
    float beta26 = params[i]/10.0f; ++i;// age over 40 [54]
    float beta27 = params[i]/10.0f; ++i;// type1 [55]
    float beta28 = params[i]/10.0f; ++i;// type2 [56]
    // Wage parameters in blue collar
    float beta30[RG_SIZE]; // constant for every region[57...63] 
    for (j = 0; j < RG_SIZE; ++j, ++i)
    {
         beta30[j] = params[i]/10.0f;
    }
    float beta31_1 = params[i] / 100.0f; ++i; // schooling[64]
    float beta31_3 = params[i] / 100.0f; ++i; //schooling[65]
    float beta31_2 = beta31_3;
    float beta32 = beta22; // experience in USSR - as in white collar
    float beta33 = beta23; // exp^2  in USSR - as in white collar
    set_param(beta34)// experience in Israel[66]
    beta34 = beta34 / 100.0f;
    set_param(beta35)// exp^2 in Israel[67]
    beta35 = beta35 / 1000.0f;
    float beta36 = beta26; // age over 40 - as in white collar
    set_param(beta37)// type1[68]
    beta37 = beta37 / 10.0f;
    set_param(beta38)// type2[69]
    beta38 = beta38 / 10.0f;
    // Travelling Costs
    float tc[TC_SIZE]; //[70..72,153,154,155]
    tc[0] = expf(params[i]);++i;
    tc[1] = expf(params[i]);++i;
    tc[2] = expf(params[i]);++i;
    // tc3, tc4, tc5 (last three) parameter is added at the end
    // probability of Losing Job, by Type
    set_param_array(aw, STATE_SIZE) // Probability of Losing Job in white color - type 0,1,2[73..75]
    float ab[STATE_SIZE]; // Probability of Losing Job in blue color - type 0 - types 1 and 2 are the same as in white[76]
    ab[0] = params[i]; ++i;
    ab[1] = aw[1];
    ab[2] = aw[2];
    // Job offer parameters - white collar
    set_param_array(lamda20, RG_SIZE) // constant for every region[77...83]
    float lamda21_1 = params[i]/10.0f; ++i; // schooling[84]
    float lamda21_3 = params[i]/10.0f; ++i; //schooling[85]
    float lamda21_2 = lamda21_3;
    set_param(lamda22) // unemployment[86]
    set_param(lamda23) // age at arrival[87]
    lamda23 = lamda23/100.0f;
    set_param(lamda25) // time[88]
    set_param(lamda26) // time^2[89]
    set_param(lamda27) // type1[90]
    set_param(lamda28) // type2[91]
    set_param(psai_w) // t==0 - a new one[92]
    set_param(psai_b) // t==0 - a new one[93]

    // Job (full time)  offer parameters - blue collar
    set_param_array(lamda30, RG_SIZE) // constant for every region[94...100]
    //float lamda31 = lamda21; // schooling - the same as in white collar
    float lamda32 = lamda22; // unemployment - the same as in white collar
    set_param(lamda35) // time[101]
    set_param(lamda36) // time^2[102]
    set_param(lamda37) // type1[103]
    set_param(lamda38) // type2[104]
    // type probability
    set_param(type1) // type1 probability[105]
    set_param(type2) // type2 probability[106]
    // Standard deviation of measurement errors
    set_param(error_w) // Standard deviation of measurement errors for wages[107]
    error_w = error_w / 10.0f;
    set_param(error_h) // Standard deviation of measurement errors for houseing[108]
    error_h = expf(error_h / 10.0f);
    set_param(error_c) // Base classification error rate[109]
    error_c = expf(error_c);
    error_c = error_c/(1.0f+error_c);
    // sigma
    set_param_array(sgma, 3) // sigma for white and blue collar and unemployment [110...112]
    sgma[0] = expf(sgma[0]);
    sgma[1] = expf(sgma[1]);
    set_param(wage_error) //punishment for miss the region or occupation in wage [113]
    wage_error = expf(wage_error);      
    set_param(row_w) // serial corellation for white collar [114]
    set_param(row_b) // serial corellation for blue collar  [115]

    set_param_array(R, RG_SIZE) //rent by area [116...122]
    
    set_param_array(psi1, RG_SIZE) // Married by region[123...129]
    set_param_array(psi3, RG_SIZE) // married*kids by region[130...136]

    set_param(lamda29) // kids w [137]
    set_param(lamda39) // kids b full [138]
    set_param(lamda49) // kids b part [139]
    set_param(lamda45) // time for part time blue  [140]
    set_param(lamda46) // time square for part time blue [141]

    set_param(alfa3) // return for leisure at unemployment, 0.5*alfa3 is return for leisure at part time [142]
    set_param(alfa20) // hustband's education influence on moving cost [143]
    set_param(lamda33) // age at arrival - blue full [144]
    lamda33 = lamda33/100.0f;
    set_param(lamda43) // age at arrival - blue part [145]
    lamda43 = lamda43/100.0f;
    set_param(part_wage_factor) // part time wage factor [146]
    set_param(type1_edu_15) // part time wage factor [147]
    set_param(type1_edu_17) // part time wage factor [148]
    set_param(type2_edu_15) // part time wage factor [149]
    set_param(type2_edu_17) // part time wage factor [150]
    set_param(moving_cost_0) // moving cost from region 0 to any region [151]
    moving_cost_0 = -expf(moving_cost_0);
    set_param(lamda233) // return on experience in USSR as engineer [152]
    tc[3] = expf(params[i]); ++i; // travel cost 3 [153]
    tc[4] = expf(params[i]); ++i; // travel cost 4 [154]
    tc[5] = expf(params[i]); ++i; // travel cost 5 [155]
    // P_W_ERROR[0] = 0, P(x<P_W_ERROR[1]) = 10%, P(x<P_W_ERROR[2]) = 30%, P(x<P_W_ERROR[3]) = 50%, P(x<P_W_ERROR[4]) = 70%, P(x<P_W_ERROR[5]) = 90%
    static const float P_W_ERROR[D_WAGE] = {0.0, -1.281551, -0.524401, 0.0, 0.524401, 1.281551};
    // P(x<P_W_ERROR_RNG[1]) = 20%, P(x<P_W_ERROR_RNG[2]) = 40%, P(x<P_W_ERROR_RNG[3]) = 60%, P(x<P_W_ERROR_RNG[4]) = 80%
    static const float P_W_ERROR_RNG[D_WAGE] = {0.0, -0.841621, -0.253347, 0.253347, 0.841621, 0.0};
    float travel_cost_arr[RG_SIZE][RG_SIZE];
#ifdef SIMULATION
    
    travel_cost_arr[0][1] = tc[0]*((sim_type == TC_SIM) ? (1.0f - sim_percent) : 1.0f);
    travel_cost_arr[0][2] = tc[4]*((sim_type == TC_SIM) ? (1.0f - sim_percent) : 1.0f);
    travel_cost_arr[3][4] = tc[3]*((sim_type == TC_SIM) ? (1.0f - sim_percent) : 1.0f);
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
    travel_cost_arr[2][0] = tc[5];
    travel_cost_arr[4][3] = tc[3];
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
    travel_cost_arr[0][2] = tc[4];
    travel_cost_arr[3][4] = tc[3];
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
    // no symmetry between tc[0][2] and tc[2][0]
    travel_cost_arr[2][0] = tc[5];
   // travel_cost_arr[0][1] = 
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
            printf("%d   \t", h_rg + 1);
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
    unsigned long work_rg_notype_edu_distribution[EDU_LEVELS][RG_SIZE][T]={{{0}}};
    unsigned long work_rg_notype_edu_distribution_count[EDU_LEVELS][T]={{0}};
    // occupation distribution
    unsigned long occ_distribution[TYPE_SIZE][ALL_STATE_SIZE][T]={{{0}}};
    unsigned long occ_distribution_count[TYPE_SIZE][T]={{0}};
    unsigned long occ_notype_distribution[ALL_STATE_SIZE][T]={{0}};
    unsigned long occ_notype_distribution_count[T]={0};
    unsigned long occ_notype_edu_distribution[EDU_LEVELS][ALL_STATE_SIZE][T]={{{0}}};
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
    float wage_blue_rg_sum[TYPE_SIZE][RG_SIZE][2];
    memset(wage_blue_rg_sum, '\0', sizeof(wage_blue_rg_sum));
    unsigned long wage_blue_rg_count[TYPE_SIZE][RG_SIZE][2]={{{0}}};
    float wage_blue_rg_notype_sum[RG_SIZE][2];
    memset(wage_blue_rg_notype_sum, '\0', sizeof(wage_blue_rg_notype_sum));
    unsigned long wage_blue_rg_notype_count[RG_SIZE][2]={{0}};
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
    float wage_blue_sum[TYPE_SIZE][2];
    memset(wage_blue_sum, '\0', sizeof(wage_blue_sum));
    unsigned long wage_blue_count[TYPE_SIZE][2]={{0}};
    float wage_blue_notype_sum[2] = {0.0f, 0.0f};
    unsigned long wage_blue_notype_count[2]={0};
#endif
    unsigned long counter_true = 0;
    unsigned long counter_false = 0;

    ///////////////////////////////////////////////////////////////////////////////////
    // The Program Start Here
    //////////////////////////////////////////////////////////////////////////////////  

    // string likelihood per individual    
    double like_arr[OBSR];
#ifdef TRACE
    float PROB_T0[OBSR];
#endif
    float PROB_T1[OBSR];
    float PROB_T2[OBSR];
#ifdef SIMULATION
    double total_max_utility[T];
    memset(total_max_utility, '\0', sizeof(total_max_utility));
    unsigned long total_max_utility_count[T] = {0};
    double total_benefit = 0.0;
#endif

    for (unsigned short I = 0; I < OBSR; ++I)
    {
    const unsigned short    AGE = AGE_arr[I];
    unsigned short  PERIODS = PERIODS_arr[I];
    const unsigned short    REP1 = REP1_arr[I];
    const unsigned short    REP2 = REP2_arr[I];
    const unsigned short    REP3 = REP3_arr[I];
#ifndef SIMULATION
    const float             WAGE = WAGE_arr[I];
    const unsigned long     RENT_MORT = RENT_MORT_arr[I];
    const unsigned short    D_MORT = D_MORT_arr[I];
#else
    unsigned short          IND_FILTER =  (sim_type == MARRIED_SIM) ? IND_FILTER_arr[I] : 1;
#endif
    const float             USSR_ENG_EXP = USSR_ENG_EXP_arr[I];
    float rent[RG_SIZE][TYPE_SIZE];
    float original_rent[RG_SIZE][TYPE_SIZE];
    float prob_nonfired_w[TYPE_SIZE];
    float prob_nonfired_b[TYPE_SIZE];
    float const_lamda_work_2w[TYPE_SIZE];
    float const_lamda_work_2b_full[TYPE_SIZE];
    float const_lamda_work_2b_part[TYPE_SIZE];
    float t_const_tmp_w[TYPE_SIZE];
    float t_const_tmp_b[TYPE_SIZE];
    float const_taste[RG_SIZE];
    float husband[RG_SIZE];
    float moving_cost[TYPE_SIZE];
    short int HUSBAND_EDU_LEVEL = -1;
    float rent_for_all_regions[TYPE_SIZE];

    for (unsigned short type = 0; type < TYPE_SIZE; ++type)
    {
        // set variables per individual from loaded arrays
        const unsigned short    M = M_arr[I];
        const unsigned short    KIDS = KIDS_arr[I];
        const unsigned short    EXP_U = EXP_U_arr[I];
        const unsigned long     EXP_U_SQ = EXP_U*EXP_U;
        const unsigned char     SCHOOL1 = (SCHOOL_arr[I] < 16) ? 1 : 0;
        const unsigned char     SCHOOL2 = (SCHOOL_arr[I] == 17) ? 1 : 0;
        const unsigned char     SCHOOL3 = (SCHOOL_arr[i] > 17 || SCHOOL_arr[i] == 16) ? 1 : 0;
        //const unsigned short    TYPE1 = TYPE1_arr[I];
        //const unsigned short    TYPE2 = TYPE2_arr[I];
        const unsigned short    TYPE1 = (type == 1);
        const unsigned short    TYPE2 = (type == 2);
        const unsigned short    TYPE0 = (type == 0);
        assert(TYPE1 + TYPE2 + TYPE0 == 1);
        if (M != 0 && HUSBAND_EDU_arr[I] == 99)
        {
            // if education is missing set it to 13 years (average)
             HUSBAND_EDU_arr[I] = 12;
        }
        if (M == 0 &&  HUSBAND_EDU_arr[I] != 99)
        {
            // if not married, set to 99
            HUSBAND_EDU_arr[I] = 99;
        }
        const unsigned short    HUSBAND_EDU = HUSBAND_EDU_arr[I];
        if (HUSBAND_EDU != 99)
        {
            for (unsigned int edu_level = 0; edu_level < EDU_LEVELS; ++edu_level)
            {
                if (HUSBAND_EDU >= edu_lower[edu_level] && HUSBAND_EDU <= edu_upper[edu_level])
                {
                    HUSBAND_EDU_LEVEL = edu_level;
                    break;
                }
            }
        }
        const unsigned int husband_edu_level_15 = HUSBAND_EDU_LEVEL == 1 ? 1 : 0;
        const unsigned int husband_edu_level_17 = HUSBAND_EDU_LEVEL == 2 ? 1 : 0;
        const float type1_edu = type1 + type1_edu_15*husband_edu_level_15 + type1_edu_17*husband_edu_level_17;
        const float type2_edu = type2 + type2_edu_15*husband_edu_level_15 + type2_edu_17*husband_edu_level_17;

        PROB_T1[I] = expf(type1_edu)/(1.0+(expf(type1_edu)+expf(type2_edu)));  
        PROB_T2[I] = expf(type2_edu)/(1.0+(expf(type1_edu)+expf(type2_edu)));
#ifdef TRACE
        PROB_T0[I] = 1.0-PROB_T1[I]-PROB_T2[I];
#endif
        rent_for_all_regions[type] = gama1*M+gama2*KIDS+gama3*TYPE1+gama4*TYPE2; //global rent without gama0 by region

        moving_cost[type] = TYPE1*alfa2[1]+TYPE2*alfa2[2]+TYPE0*alfa2[0] + alfa20*HUSBAND_EDU;
        for (unsigned short rg = 0; rg < RG_SIZE; ++rg)
        {
            const_taste[rg] = teta0[rg]+teta1[rg]*REP1+teta2[rg]*REP2+teta3[rg]*REP3; //taste for housing in a specific region - equation 5 page 13
            // define const taste, move the final calculation into draws
            rent[rg][type] = 6.0f*expf(gama0[rg]+rent_for_all_regions[type])*RENT_REF_PARAM; // full cost of housing rent/mor - equation 6
            husband[rg] = psi1[rg]*M  + psi3[rg]*M*KIDS;
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

        //probability of not losing your job in white collar - equation 5 page 12
        prob_nonfired_w[type] = 1.0f/(1.0f + expf(TYPE1*aw[1]+TYPE2*aw[2]+TYPE0*aw[0]));
        //probability of not losing your job in blue collar
        prob_nonfired_b[type] = 1.0f/(1.0f + expf(TYPE1*ab[1]+TYPE2*ab[2]+TYPE0*ab[0]));
        const_lamda_work_2w[type] = (lamda21_1*SCHOOL1 + lamda21_2*SCHOOL2 + lamda21_3*SCHOOL3) + 
                                            lamda23*AGE+lamda27*TYPE1+lamda28*TYPE2+lamda29*KIDS+lamda233*USSR_ENG_EXP; //part of the probability of getting job offer in white - page 13
        const_lamda_work_2b_full[type] = lamda33*AGE+lamda37*TYPE1+lamda38*TYPE2+lamda39*KIDS; //part of the probability of getting job offer in blue - page 13
        const_lamda_work_2b_part[type] = lamda43*AGE+lamda37*TYPE1+lamda38*TYPE2+lamda49*KIDS; //part of the probability of getting job offer in blue - page 13
        t_const_tmp_w[type] = (beta21_1*SCHOOL1 + beta21_2*SCHOOL2 + beta21_3*SCHOOL3) + 
                                            beta22*EXP_U+beta23*EXP_U_SQ+beta27*TYPE1+beta28*TYPE2;  //part of the wage equation  white collar- equation 7 page 14
        t_const_tmp_b[type] = (beta31_1*SCHOOL1 + beta31_2*SCHOOL2 + beta31_3*SCHOOL3) + 
                                            beta32*EXP_U+beta33*EXP_U_SQ+beta37*TYPE1+beta38*TYPE2;  //part of the wage equation  blue collar- equation 7 page 14

        const float ret = 65.0f - (float)AGE + 10.0f;
        const float one_by_beta_pown = powf(one_by_beta, ret);
        const float terminal = (one_by_beta_pown - 1.0f)/(one_by_beta_pown*(one_by_beta - 1.0f));

        for (unsigned short t = T; t > 0; --t)
        {
            // loop over periods (decdresing)
            const unsigned long t_sq = t*t;
            const unsigned short age40 = (((float)AGE + (float)t/2.0f) > 39.5f);
            //part of the probability of getting job offer in white - page 13 (miss:constant by region +come from unemp)
            const float lamda_work_2w = const_lamda_work_2w[type]+lamda25*t+lamda26*(float)t_sq;
            //part of the probability of getting job offer in blue - page 13(miss:constant by region +come from unemp) 
            const float lamda_work_2b_full = const_lamda_work_2b_full[type]+lamda35*t+lamda36*(float)t_sq;
            const float lamda_work_2b_part = const_lamda_work_2b_part[type]+lamda45*t+lamda46*(float)t_sq;
            const float k_const_tmp_w = t_const_tmp_w[type]+beta26*age40;   //part of the wage equation  white collar- equation 7 page 14 (miss:const by region+exp+exp^2)
            const float k_const_tmp_b = t_const_tmp_b[type]+beta36*age40;   //part of the wage equation  blue collar- equation 7 page 14 (miss:const by region+exp+exp^2)

            float prob_work_2w[RG_SIZE];
            float prob_ue_2w[RG_SIZE];
            float prob_work_2b_full[RG_SIZE];
            float prob_work_2b_part[RG_SIZE];
            float prob_ue_2b_full[RG_SIZE];
            float prob_ue_2b_part[RG_SIZE];

            for (unsigned short rg = 0; rg < RG_SIZE; ++rg)
            {
                const float tmp_lamda_w = lamda_work_2w + lamda20[rg]; // lamda21*SCHOOL+lamda23*AGE+lamda27*TYPE1+lamda28*TYPE2+lamda25*t+lamda26*t_sq+lamda20[rg]
                float tmp_exp_w = expf(tmp_lamda_w);
                prob_work_2w[rg] = tmp_exp_w/(1.0+tmp_exp_w); //probability to get job offer in white if come from work
                // lamda21*SCHOOL+lamda23*AGE+lamda27*TYPE1+lamda28*TYPE2+lamda25*t+lamda26*t_sq+lamda20[rg]+lamda22
                tmp_exp_w = expf(tmp_lamda_w + lamda22);
                prob_ue_2w[rg] = tmp_exp_w/(1.0+tmp_exp_w); //probability to get job offer in white if come from unemployment

                const float tmp_lamda_b_full = lamda_work_2b_full+lamda30[rg]; // lamda33*AGE+lamda37*TYPE1+lamda38*TYPE2+lamda35*t+lamda36*t_sq+lamda30[rg]
                float tmp_exp_b_full = expf(tmp_lamda_b_full);
                const float tmp_lamda_b_part = lamda_work_2b_part+lamda30[rg]; // lamda33*AGE+lamda37*TYPE1+lamda38*TYPE2+lamda35*t+lamda36*t_sq+lamda30[rg]
                float tmp_exp_b_part = expf(tmp_lamda_b_part);

                prob_work_2b_full[rg] = tmp_exp_b_full/(1.0+tmp_exp_b_full+tmp_exp_b_part); //probability to get job offer in blue full if come from work
                // lamda33*AGE+lamda37*TYPE1+lamda38*TYPE2+lamda35*t+lamda36*t_sq+lamda30[rg]+lamda32
                prob_work_2b_part[rg] = tmp_exp_b_part/(1.0+tmp_exp_b_part+tmp_exp_b_full); //probability to get job offer in blue full if come from work
                // lamda33*AGE+lamda37*TYPE1+lamda38*TYPE2+lamda35*t+lamda36*t_sq+lamda30[rg]+lamda32

                tmp_exp_b_full = expf(tmp_lamda_b_full + lamda32);
                tmp_exp_b_part = expf(tmp_lamda_b_part + lamda32);
                prob_ue_2b_full[rg] = tmp_exp_b_full/(1.0+tmp_exp_b_full+tmp_exp_b_part); //probability to get job offer in blue full if come from unemployment
                prob_ue_2b_part[rg] = tmp_exp_b_part/(1.0+tmp_exp_b_part+tmp_exp_b_full); //probability to get job offer in blue full if come from unemployment

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
            // TODO: should we run to t? min(t,12) ? or PERIODS?
            //const unsigned short MAX_K = std::min(t,MAXIMUM_K_INDEX);
            for (unsigned short k = 0 ; k <= MAXIMUM_K_INDEX; ++k)
            {
                // loop over experience
#ifdef PERF_TRACE
                timeval tv = tic();
#endif
                const float real_k = index_to_k(k);
                const float k_sq = real_k*real_k;

                //part of the wage equation  white collar- equation 7 page 14 (adding exp and exp^2 still miss:const by region)
                const float rg_const_tmp_w = k_const_tmp_w+beta24*real_k+beta25*k_sq;
                //part of the wage equation  blue collar- equation 7 page 14 (adding exp and exp^2 still miss:const by region)
                const float rg_const_tmp_b = k_const_tmp_b+beta34*real_k+beta35*k_sq;

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
                    float nonfired_2b_full[RG_SIZE][D_WAGE];
                    float nonfired_2b_part[RG_SIZE][D_WAGE];
                    float wage_nonfired_2w[RG_SIZE][D_WAGE];
                    float wage_work_2w[RG_SIZE];
                    float wage_ue_2w[RG_SIZE];
                    float taste[RG_SIZE];
                    float wage_b[RG_SIZE];
                    float wage_w[RG_SIZE];
                    bool ue_2b_full[RG_SIZE];
                    bool work_2b_full[RG_SIZE];

                    for (unsigned short rg = 0; rg < RG_SIZE; ++rg)
                    {
                        float choose_ue_emax;
                        float choose_b_emax[D_WAGE];
                        float tmp1 = epsilon_b(draw,I,t,rg,WHITE,type);
                        float tmp2 = epsilon_b(draw,I,t,rg,BLUE,type);
                        float tmp3 = epsilon_b(draw,I,t,rg,UE,type);
                        taste[rg] = const_taste[rg];// + exp(sgma[3]*tmp0);
                        
                        for (unsigned short dwage = 0; dwage < D_WAGE; ++dwage)
                        {
                            // tmp1_w = beta21*SCHOOL+beta22*EXP_U+beta23*EXP_U_SQ+beta27*TYPE1+beta28*TYPE2+beta26*age40+beta24*k+beta25*k_sq+beta20[rg]  
                            wage_w[dwage] = 6.0f*expf(rg_const_tmp_w + beta20[rg] + sgma[0]*(tmp1 + row_w*P_W_ERROR[dwage]))*WAGE_REF_PARAM;
                            // tmp1_b = beta31*SCHOOL+beta32*EXP_U+beta33*EXP_U_SQ+beta37*TYPE1+beta38*TYPE2+beta36*age40+beta34*k+beta35*k_sq+beta30[rg]
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
                            
                            wage_nonfired_2w[rg][dwage] = draw_wage(wage_w[dwage], prob_nonfired_w[type]);            //equal wage if ind wasn't fired  and -inf if was fired  
                            const float wage_nonfired_2b_full = draw_wage(wage_b[dwage], prob_nonfired_b[type]);      //equal wage if ind wasn't fired  and -inf if was fired
                            const float wage_nonfired_2b_part = draw_wage(wage_b[dwage]/2.0, prob_nonfired_b[type]) + alfa3/2.0; //equal wage if ind wasn't fired  and -inf if was fired
                            if (t == T)
                            {
                                choose_ue_emax = 0.0f;
                                choose_b_emax[dwage]= 0.0f;
                                for (unsigned short w_rg = 0; w_rg < RG_SIZE; ++w_rg)
                                {
                                    choose_w_emax[rg][w_rg][dwage]= 0.0f;
                                }
                            } 
                            else
                            {
                                choose_ue_emax = beta*EMAX(t+1,k,rg,0,UE,0,type);
                                choose_b_emax[dwage] = beta*EMAX(t+1,k_to_index(real_k+1.0),rg,0,BLUE,D_W_B,type);

                                for (unsigned short w_rg = 0; w_rg < RG_SIZE; ++w_rg)
                                {
                                    choose_w_emax[rg][w_rg][dwage] = beta*EMAX(t+1,k_to_index(real_k+1.0),rg,w_rg,WHITE,D_W_W,type);
                                }
                            }

                            nonfired_2b_full[rg][dwage] = wage_nonfired_2b_full + taste[rg] - rent[rg][type] + husband[rg] + choose_b_emax[dwage];
                            nonfired_2b_part[rg][dwage] = wage_nonfired_2b_part + taste[rg] - rent[rg][type] + husband[rg] + choose_b_emax[dwage];
                        } // close dwage             

                        wage_ue_2w[rg] = draw_wage(wage_w[0], prob_ue_2w[rg]);          // equal wage if ind come fron ue and got an offer and -inf if didn't
                        float wage_ue_2b = draw_blue_wage(wage_b[0], prob_ue_2b_full[rg], prob_ue_2b_part[rg], ue_2b_full[rg], part_wage_factor); // equal wage if ind come fron ue and got an offer and -inf if not
                        wage_ue_2b += (ue_2b_full[rg] == false ? alfa3/2.0 : 0.0);      // add part time alfa3 if needed
                        wage_work_2w[rg] = draw_wage(wage_w[0], prob_work_2w[rg]);      // equal wage if ind come from and got an offer and -inf if didn't
                        float wage_work_2b = draw_blue_wage(wage_b[0], prob_work_2b_full[rg], prob_work_2b_part[rg], work_2b_full[rg], part_wage_factor); // equal wage if ind come from and got an offer and -inf
                        wage_work_2b += (work_2b_full[rg] == false ? alfa3/2.0 : 0.0);  // add part time alfa3 if needed

                        // the equivalent to "wage" when UE is chosen
                        choose_ue[rg] =  taste[rg] - rent[rg][type] + husband[rg] + expf(sgma[2]*tmp3) + alfa3 + choose_ue_emax;

                        const float tmp = taste[rg] - rent[rg][type] + husband[rg] + choose_b_emax[0];
                        ue_2b[rg] = wage_ue_2b + tmp;
                        work_2b[rg] = wage_work_2b + tmp;

                    }//close rg

                    float ue_2w[RG_SIZE][RG_SIZE];
                    float work_2w[RG_SIZE][RG_SIZE];
                    float nonfired_2w[RG_SIZE][RG_SIZE][D_WAGE];

                    for (unsigned short h_rg = 0; h_rg < RG_SIZE; ++h_rg)
                    {
                        for (unsigned short w_rg = 0; w_rg < RG_SIZE; ++w_rg)
                        {
                            float tmp = taste[h_rg] - rent[h_rg][type] + husband[h_rg] - travel_cost(h_rg,w_rg) + choose_w_emax[h_rg][w_rg][0];
                            ue_2w[h_rg][w_rg] = wage_ue_2w[w_rg] + tmp;
                            work_2w[h_rg][w_rg] = wage_work_2w[w_rg] + tmp;
                            for (unsigned short dwage = 0; dwage < D_WAGE; ++dwage)
                            {
                                nonfired_2w[h_rg][w_rg][dwage] = 
                                wage_nonfired_2w[w_rg][dwage] + taste[h_rg] - rent[h_rg][type] + husband[h_rg] - travel_cost(h_rg,w_rg) + choose_w_emax[h_rg][w_rg][dwage];
                            }   
                        }//close to_w_rg
                    }//close rg

                    /* utility vector
                        ---------------------------------------------------------------------------------------------------------------------------
                        || 0-6       || 7-13       || 14-20      || 21-27    ||28 - 34  || 35 - 41  || 42 - 48 || 49 - 55  || 56 - 62 || 63 - 69 ||
                        ---------------------------------------------------------------------------------------------------------------------------
                        || UE        || Blue, full || Blue, part || White0   || White1  || White2   || White3  || White4   || White5  || White6  ||
                        ---------------------------------------------------------------------------------------------------------------------------
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
                                get_max(from_b_max_utility_arr[dwage], nonfired_2b_full[h_rg][dwage]);
                                get_max(from_b_max_utility_arr[dwage], nonfired_2b_part[h_rg][dwage]);
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
                                nonfired_2b_full[h_rg][dwage] *= (1.0f+terminal);
                                nonfired_2b_part[h_rg][dwage] *= (1.0f+terminal);
                                get_max(from_b_max_utility_arr[dwage], nonfired_2b_full[h_rg][dwage]);
                                get_max(from_b_max_utility_arr[dwage], nonfired_2b_part[h_rg][dwage]);
                                from_b_max_utility_arr[dwage] = from_b_max_utility_h;
                                // stay in blue and live in the same region
                                get_max(from_b_max_utility_arr[dwage], nonfired_2b_full[h_rg][dwage]);
                                get_max(from_b_max_utility_arr[dwage], nonfired_2b_part[h_rg][dwage]);
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
                printf("calculating EMAX for I = %d type = %d t = %d and k = %d took: %f seconds (estimated total = %f minutes)\n",
                       I, type, t, k, toc(tv), toc(tv)*ITER_COUNT/60.0);
#endif
            } // end loop over experience
        } // end loop over periods
    } // end loop over types
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
        //unsigned short I_id;
        unsigned short type = draw_type(PROB_T1[I], PROB_T2[I]);
        
        //{
        //    div_t I_info = div(I,OBSR);
        //    I_id = (unsigned short)I_info.rem;
        //    I_type = (unsigned short)I_info.quot;
        //}

#ifdef SIMULATION
        const unsigned short draws_f = (unsigned short)((type == 0) ? TYPE_0_OF_1000 : ((type == 1) ? TYPE_1_OF_1000 : TYPE_2_OF_1000));
#else 
        const unsigned short draws_f = DRAWS_F;
#endif // SIMULATION

        for (unsigned short draw = 0; draw < draws_f; ++draw)
        {
#ifdef FULL_TRACE
            printf("%d %d %d %d ", I, (REP1 ? 1 : (REP2 ? 2 : (REP3 ? 3 : 4))), type, draw); 
#endif
            const unsigned short BLUE_STATE_FULL = 0;
            const unsigned short BLUE_STATE_PART = 1;
            unsigned short dwage_b = 0;
            unsigned short dwage_w = 0;
            unsigned short blue_state = BLUE_STATE_FULL; // 0 blue full, 1 blue part
            unsigned short from_state = UE;  
            unsigned short from_h_rg = 0;
            unsigned short from_w_rg = 0;
            float real_k = 0.0;

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
            short rg = (unsigned short)((float)RG_SIZE*(rand()/(RAND_MAX + 1.0f)));
#endif
            for (unsigned short t = 0; t < PERIODS; ++t)// loop over periods
            {
#ifndef RANDOM_SELECTION
                short rg = live(I,t);
#endif
                if (rg < 0)
                {
#ifdef FULL_TRACE_INDEX
                    printf("%d ", 999);
#elif FULL_TRACE_WAGE
                    printf("%.3f ",  0.0);
#endif
                    continue;
                }
                bool w_wage_flag = false;
                bool b_wage_flag = false;
                const unsigned short age40 = (((float)AGE + (float)t/2.0f) > 39.5f);
                const float k_sq = real_k*real_k;
                const unsigned long t_sq = t*t;
                const float lamda_work_2w = const_lamda_work_2w[type]+lamda25*t+lamda26*(float)t_sq; //part of the probability of getting job offer in white - page 13
                const float lamda_work_2b_full = const_lamda_work_2b_full[type]+lamda35*t+lamda36*(float)t_sq; //part of the probability of getting job offer in blue full - page 13
                const float lamda_work_2b_part = const_lamda_work_2b_part[type]+lamda45*t+lamda46*(float)t_sq; //part of the probability of getting job offer in blue part - page 13
                //part of the wage equation  white collar- equation 7 page 14 (adding exp and exp^2 still miss:const by region)
                const float rg_const_tmp_w = t_const_tmp_w[type]+beta26*age40+beta24*real_k+beta25*k_sq;
                //part of the wage equation  blue collar- equation 7 page 14 (adding exp and exp^2 still miss:const by region)
                const float rg_const_tmp_b = t_const_tmp_b[type]+beta36*age40+beta34*real_k+beta35*k_sq;
                
                float work_2b[RG_SIZE];
                float ue_2b[RG_SIZE];
                float nonfired_2b_full[RG_SIZE];
                float nonfired_2b_part[RG_SIZE];
                float ue_2w[RG_SIZE][RG_SIZE];
                float work_2w[RG_SIZE][RG_SIZE];
                float choose_ue[RG_SIZE];
                float nonfired_2w[RG_SIZE][RG_SIZE];
                float taste[RG_SIZE];
                float wage_nonfired_2w[RG_SIZE];
                float wage_b[RG_SIZE];
                float wage_w[RG_SIZE];
                float wage_b_non_f[RG_SIZE];
                float wage_w_non_f[RG_SIZE];
                float wage_ue_2w[RG_SIZE];
                float wage_work_2w[RG_SIZE];
                unsigned short D_W_W[RG_SIZE];
                unsigned short D_W_B[RG_SIZE];
                bool ue_2b_full[RG_SIZE];
                bool work_2b_full[RG_SIZE];

                {
                    assert(rg > -1);

                    // adjust rent with R
                    rent[rg][type] = rent[rg][type]*(1.0 + R[rg]);

                    float prob_work_2w;
                    float prob_ue_2w;
                    float prob_ue_2b_full;
                    float prob_ue_2b_part;
                    float prob_work_2b_full;
                    float prob_work_2b_part;

                    {
                        
                        
                        const float tmp_lamda_b_full = lamda_work_2b_full + lamda30[rg];        // lamda33*AGE+lamda37*TYPE1+lamda38*TYPE2+lamda35*t+lamda36*t_sq+lamda30[rg]
                        float tmp_exp_b_full = expf(tmp_lamda_b_full);
                        const float tmp_lamda_b_part = lamda_work_2b_part + lamda30[rg];        // lamda33*AGE+lamda37*TYPE1+lamda38*TYPE2+lamda35*t+lamda36*t_sq+lamda30[rg]
                        float tmp_exp_b_part = expf(tmp_lamda_b_part);

                        prob_work_2b_full = tmp_exp_b_full/(1.0+tmp_exp_b_full+tmp_exp_b_part); // probability to get job offer in blue (full) if come from work
                        prob_work_2b_part = tmp_exp_b_part/(1.0+tmp_exp_b_part+tmp_exp_b_full); // probability to get job offer in blue (part) if come from work
                        
                        tmp_exp_b_full = expf(tmp_lamda_b_full + lamda32);            // lamda33*AGE+lamda37*TYPE1+lamda38*TYPE2+lamda35*t+lamda36*t_sq+lamda30[rg]+lamda32
                        tmp_exp_b_part = expf(tmp_lamda_b_part + lamda32);            // lamda33*AGE+lamda37*TYPE1+lamda38*TYPE2+lamda35*t+lamda36*t_sq+lamda30[rg]+lamda32

                        prob_ue_2b_full = tmp_exp_b_full/(1.0+tmp_exp_b_full+tmp_exp_b_part);            // probability to get job offer in blue (full) if come from unemployment
                        prob_ue_2b_part = tmp_exp_b_part/(1.0+tmp_exp_b_part+tmp_exp_b_full);            // probability to get job offer in blue (part) if come from unemployment

                        if (t == 0)
                        {
                            prob_ue_2w = prob_ue_2w*psai_w;
                            prob_ue_2b_full = prob_ue_2b_full*psai_b;
                            prob_ue_2b_part = prob_ue_2b_part*psai_b;
                        }
                    }

                    taste[rg] = const_taste[rg];
                    for (unsigned int w_rg = 0; w_rg < RG_SIZE; ++w_rg)
                    {
                        const float tmp1 = epsilon_f(draw,I,t,w_rg,WHITE,type);
                        const float tmp_lamda_w = lamda_work_2w + lamda20[w_rg];  // lamda21*SCHOOL+lamda23*AGE+lamda27*TYPE1+lamda28*TYPE2+lamda25*t+lamda26*t_sq+lamda20[rg]
                        float tmp_exp_w = expf(tmp_lamda_w);
                        prob_work_2w = tmp_exp_w/(1.0+tmp_exp_w);               // probability to get job offer in white if come from work
                        tmp_exp_w = expf(tmp_lamda_w + lamda22);                // lamda21*SCHOOL+lamda23*AGE+lamda27*TYPE1+lamda28*TYPE2+lamda25*t+lamda26*t_sq+lamda20[rg]+lamda22
                        prob_ue_2w = tmp_exp_w/(1.0+tmp_exp_w);                 // probability to get job offer in white if come from unemployment
                        if (t == 0)
                        {
                            prob_ue_2w = prob_ue_2w*psai_w;
                        }
                        const float tmpdw = tmp1 + row_w*P_W_ERROR[dwage_w];
                        D_W_W[w_rg] = get_discrete_index(tmpdw);
                        wage_w_non_f[w_rg] = 6.0f*expf(rg_const_tmp_w + beta20[w_rg] + sgma[0]*(tmp1+row_w*P_W_ERROR[dwage_w]))*WAGE_REF_PARAM;
                        wage_w[w_rg] = 6.0f*expf(rg_const_tmp_w + beta20[w_rg] + sgma[0]*tmp1)*WAGE_REF_PARAM;
                        wage_ue_2w[w_rg] = draw_wage(wage_w[w_rg], prob_ue_2w);                           //equal wage if i come fron ue and got an offer and -inf if didn't
                        wage_work_2w[w_rg] = draw_wage(wage_w[w_rg], prob_work_2w);                       //equal wage if ind come from and got an offer and -inf if didn't
                    }
                    const float tmp2 = epsilon_f(draw,I,t,rg,BLUE,type);
                    wage_b_non_f[rg] = 6.0f*expf(rg_const_tmp_b + beta30[rg] + sgma[1]*(tmp2+row_b*P_W_ERROR[dwage_b]))*WAGE_REF_PARAM;
                    wage_b[rg] = 6.0f*expf(rg_const_tmp_b + beta30[rg] + sgma[1]*tmp2)*WAGE_REF_PARAM;

#ifdef SIMULATION
                    if (sim_type == WAGE_SIM)
                    {
                        assert(0);
                    }
#endif

                    const float tmpdb = tmp2 + row_b*P_W_ERROR[dwage_b];
                    D_W_B[rg] = get_discrete_index(tmpdb);

                    // sampling the wage for each of the transitions
                    
                    wage_nonfired_2w[rg] = draw_wage(wage_w_non_f[rg], prob_nonfired_w[type]);          //equal wage if ind wasn't fired  and -inf if was fired
                    const float wage_nonfired_2b_full = draw_wage(wage_b_non_f[rg], prob_nonfired_b[type]);  //equal wage if ind wasn't fired  and -inf if was fired
                    const float wage_nonfired_2b_part = draw_wage(wage_b_non_f[rg]/2.0, prob_nonfired_b[type]) + alfa3/2.0;  //equal wage if ind wasn't fired  and -inf if was fired
                    float wage_ue_2b = draw_blue_wage(wage_b[rg], prob_ue_2b_full, prob_ue_2b_part, ue_2b_full[rg], part_wage_factor);      //equal wage if i come fron ue and got an offer and -inf if didn't
                    wage_ue_2b += (ue_2b_full[rg] == false ? alfa3/2.0 : 0.0);                      // add part time alfa3 if needed
                    float wage_work_2b = draw_blue_wage(wage_b[rg], prob_work_2b_full, prob_work_2b_part, work_2b_full[rg], part_wage_factor);  //equal wage if ind come from and got an offer and -inf if didn't
                    wage_work_2b += (work_2b_full[rg] == false ? alfa3/2.0 : 0.0);                  // add part time alfa3 if needed

                    const float choose_ue_emax = beta*EMAX(t+1,k_to_index(real_k),rg,0,UE,0,type);
                    const float choose_b_emax_non_f = beta*EMAX(t+1,k_to_index(real_k+1.0),rg,0,BLUE,D_W_B[rg],type);
                    const float choose_b_emax = beta*EMAX(t+1,k_to_index(real_k+1.0),rg,0,BLUE,0,type);

                    const float taste_rent_husband = taste[rg] - rent[rg][type] + husband[rg];
                    
                    // the equivalent to "wage" when UE is chosen
                    choose_ue[rg] =  taste_rent_husband + expf(sgma[2]*epsilon_f(draw,I,t,from_h_rg,UE,type)) + alfa3 + choose_ue_emax;
                    if (t == 0)
                    {
                        // add alfa1 utility for t=0
                        choose_ue[rg] += alfa1[rg];
                    }
                    ue_2b[rg] = wage_ue_2b + taste_rent_husband + choose_b_emax;
                    work_2b[rg] = wage_work_2b + taste_rent_husband + choose_b_emax;
                    nonfired_2b_full[rg] = wage_nonfired_2b_full + taste_rent_husband + choose_b_emax_non_f;                  
                    nonfired_2b_part[rg] = wage_nonfired_2b_part + taste_rent_husband + choose_b_emax_non_f;                  
                    
                } //end rg

                // fill white information
                for (unsigned short rg = 0; rg < RG_SIZE; ++rg)
                {
                    const float taste_rent_husband = taste[rg] - rent[rg][type] + husband[rg];
                    for (unsigned short w_rg = 0; w_rg < RG_SIZE; ++w_rg)
                    {
                        const float choose_w_emax = beta*EMAX(t+1,k_to_index(real_k+1.0),rg,w_rg,WHITE,0,type);
                        const float choose_w_emax_non_f = beta*EMAX(t+1,k_to_index(real_k+1),rg,w_rg,WHITE,D_W_W[rg],type);
                        ue_2w[rg][w_rg] = wage_ue_2w[w_rg] + taste_rent_husband - travel_cost(rg,w_rg) + choose_w_emax;
                        work_2w[rg][w_rg] = wage_work_2w[w_rg] + taste_rent_husband - travel_cost(rg,w_rg) + choose_w_emax; 
                        nonfired_2w[rg][w_rg] = wage_nonfired_2w[rg] + taste_rent_husband + choose_w_emax_non_f;
                    }
                }


                //////////////////////////////////// start maximization ////////////////////////////////////
                float choices[STATE_VECTOR_SIZE];
                for (unsigned int i = 0; i < STATE_VECTOR_SIZE; ++i)
                {
                    // initialize choices with -inf
                    choices[i] = -INFINITY;
                }

                if (from_state == UE)
                {
                    {
                        float tmp_moving_cost = 0.0;
                        if (t != 0 && from_h_rg != rg)
                        {
                            tmp_moving_cost = moving_cost[type];
                            if (from_h_rg == 0)
                            {
                                tmp_moving_cost += moving_cost_0;
                            }
                        }
                        // stay in ue
                        choices[rg] = choose_ue[rg] - tmp_moving_cost;
                        // move to blue accoreding to full/part
                        choices[ue_2b_full[rg]?rg+7:rg+14] = ue_2b[rg] - tmp_moving_cost;
                        // move to white
                        for (unsigned short w_rg = 0; w_rg < RG_SIZE; ++w_rg)
                        {
                            choices[rg+21+7*w_rg] = ue_2w[rg][w_rg] - tmp_moving_cost;
                        }
                    }
                } // close state UE
                else if (from_state == BLUE)
                {
                    {
                        float tmp_moving_cost = 0.0;
                        if (t != 0 && from_h_rg != rg)
                        {
                            tmp_moving_cost = moving_cost[type];
                            if (from_h_rg == 0)
                            {
                                tmp_moving_cost += moving_cost_0;
                            }
                        }
                        // move to ue
                        choices[rg] = choose_ue[rg] - tmp_moving_cost;
						            // stay in blue accoreding to full/part
                        if (from_h_rg != rg)
                        {
                            // move regions
                            choices[work_2b_full[rg]?rg+7:rg+14] = work_2b[rg] - tmp_moving_cost;
                        }
                        else
                        {
                            // stay in same region
                            if (blue_state == BLUE_STATE_FULL)
                            {
                                // full
                                choices[rg+7] = nonfired_2b_full[rg];
                            }
                            else if (blue_state == BLUE_STATE_PART)
                            {
                                // part
                                choices[rg+14] = nonfired_2b_part[rg];
                            }
                            else
                            {
                                assert(0);
                            }
                        }
                        // move to white
                        for (unsigned short w_rg = 0; w_rg < RG_SIZE; ++w_rg)
                        {
                            choices[rg+21+7*w_rg] = work_2w[rg][w_rg] - tmp_moving_cost;
                        }
                    }
                }
                else if (from_state == WHITE)
                {
                    {
                        float tmp_moving_cost = 0.0;
                        if (t != 0 && from_h_rg != rg)
                        {
                            tmp_moving_cost = moving_cost[type];
                            if (from_h_rg == 0)
                            {
                                tmp_moving_cost += moving_cost_0;
                            }
                        }
                        // move to ue
                        choices[rg] = choose_ue[rg]- tmp_moving_cost;
                        // move to blue accoreding to full/part
                        choices[work_2b_full[rg]?rg+7:rg+7] = work_2b[rg]- tmp_moving_cost;
                        // stay in white
                        for (unsigned short w_rg = 0; w_rg < RG_SIZE; ++w_rg)
                        {
                            if (from_w_rg != w_rg)
                            {
                                // move work regions
                                choices[rg+21+7*w_rg] = work_2w[rg][w_rg] - tmp_moving_cost;
                            }
                            else
                            {
                                // stay in same work region
                                choices[rg+21+7*w_rg] = nonfired_2w[rg][w_rg] - tmp_moving_cost;
                            }
                        }
                    }
                }
                else
                {
                    // invalid state
                    assert(0);
                }
   
                // find maximum from choices
                float max_utility = -INFINITY;
                int max_index = -1;
                for (unsigned int i = 0; i < STATE_VECTOR_SIZE; ++i)
                {
                    if (choices[i] > max_utility)
                    {
                        max_index = i;
                        max_utility = choices[i];
                    }
                }

                // check that a maximum was found
                assert(max_index > -1);
                
                if (t == PERIODS-1)
                {
                    if (max_index == from_h_rg+21+7*from_w_rg)
                    {
                        w_wage_flag = true;
                    }
                    else if (max_index == from_h_rg+7 || max_index == from_h_rg+14)
                    {
                        b_wage_flag = true;
                    }
                }

                unsigned short tmp_work_rg;
                unsigned short tmp_house_rg;
                {
                    div_t est_house_info = div(max_index,7);
                    tmp_house_rg = (unsigned short)est_house_info.rem;
                    tmp_work_rg = (unsigned short)est_house_info.quot;
                    // work region 0 = unemployment
                    // work region 1 = blue, full time
                    // work region 2 = blue, part time
                    // (work region - 3) = white work region
                }
                if (rg != tmp_house_rg)
                {
                    std::cout << "rg = " << rg << " tmp_house_rg = " << tmp_house_rg << std::endl;
                    assert(0);
                }
#ifdef FULL_TRACE_INDEX
                printf("%d ", max_index);
#elif FULL_TRACE_WAGE
                {
                    float current_wage;
                    if (tmp_work_rg == 1)
                    {
                        // blue, full
                        current_wage = ((b_wage_flag == false) ? wage_b[tmp_house_rg] : wage_b_non_f[tmp_house_rg])/6.0;
                    }
                    else if (tmp_work_rg == 2)
                    {
                        // blue, part
                        current_wage = ((b_wage_flag == false) ? wage_b[tmp_house_rg]/2.0 : wage_b_non_f[tmp_house_rg])/12.0;
                    }
                    else if (tmp_work_rg > 2)
                    {
                        // white
                        current_wage = ((w_wage_flag == false) ? wage_w[tmp_work_rg-3] : wage_w_non_f[tmp_work_rg-3])/6.0;
                    }
                    else if (tmp_work_rg == 0)
                    {
                        // unemployed
                        current_wage = 0.0;
                    }
                    else
                    {
                        // invalid state
                        assert(0);
                    }
                    if ((max_index < 7 && current_wage > 0.0) || (max_index > 6 && current_wage == 0.0))
                    {
                        printf("Invalid Combination: I=%d draw=%d t=%d selection=%d\n", I, draw, t, max_index);
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
                    if ((tmp_work_rg == 1 || tmp_work_rg == 2) && (tmp_house_rg == 4 || tmp_house_rg == 5)) // BLUE
                    {
                        // note: we ignore the difference between non-fired blue and blue
                        total_benefit += 6.0f*expf(rg_const_tmp_b + beta30[tmp_house_rg] + sgma[1]*epsilon_f(draw,I,t,tmp_house_rg,BLUE,type))*sim_percent; 
                    }
                    else if (tmp_work_rg > 2 && (tmp_work_rg-3 == 4 || tmp_work_rg-3 == 5)) // WHITE
                    {
                        // note: we ignore the difference between non-fired white and white
                        total_benefit += 6.0f*expf(rg_const_tmp_w + beta20[tmp_work_rg-3] + sgma[0]*epsilon_f(draw,I,t,(tmp_work_rg-3),WHITE,type))*sim_percent;
                    }
                    // else unemployment
                }
                else if (sim_type == TC_SIM)
                {
                    if (tmp_work_rg > 2 && !(tmp_work_rg-3 == 0 || tmp_work_rg-3 == 3 || tmp_work_rg-3 == 6))
                    {
                        // WHITE, not working in 0, 3 or 6
                        total_benefit += sim_percent*travel_cost_arr[tmp_house_rg][tmp_work_rg-3]/(1.0f - sim_percent);
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
                else if (tmp_work_rg == 1 || tmp_work_rg == 2)
                {
                    // work in blue
#ifndef SIMULATION
                    job_arr[t][draw] = BLUE;
#endif
                    // 0 blue full, 1 blue part
                    blue_state = ((tmp_work_rg == 1) ? BLUE_STATE_FULL : BLUE_STATE_PART);
                    from_state = BLUE;
                    // increase experience - by 1 for full, by 0.5 for part
                    real_k += ((blue_state == BLUE_STATE_FULL) ? 1.0 : 0.5);
                    // in blue house_rg equals work_rg
                    work_rg_arr[t][draw] = tmp_work_rg;
                    from_w_rg = tmp_work_rg;
                    dwage_b = D_W_B[tmp_house_rg];
                    dwage_w = 0;
                }
                else
                {
                    // work in white
                    tmp_work_rg -= 3; // now it is a number: 0-6
#ifndef SIMULATION
                    job_arr[t][draw] = WHITE;
#endif
                    from_state = WHITE;
                    // increase experience
                    ++real_k;
                    work_rg_arr[t][draw] = tmp_work_rg;
                    from_w_rg = tmp_work_rg;
                    dwage_b = 0;
                    dwage_w = D_W_W[tmp_work_rg];
#ifdef TRACE
#ifdef SIMULATION
                    if (IND_FILTER==1)
#endif
                    {
                        ++work_rg_distribution[type][from_w_rg][t];
                        ++work_rg_distribution_count[type][t];
                        ++work_rg_notype_distribution[from_w_rg][t];
                        ++work_rg_notype_distribution_count[t];
						if (HUSBAND_EDU_LEVEL != -1)
						{
							++work_rg_notype_edu_distribution[HUSBAND_EDU_LEVEL][from_w_rg][t];
							++work_rg_notype_edu_distribution_count[HUSBAND_EDU_LEVEL][t];
						}
                    
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
#ifdef SIMULATION
                if (IND_FILTER==1)
#endif
                {
                    ++house_distribution[type][from_h_rg][t];
                    ++house_distribution_count[type][t];
                    ++house_notype_distribution[from_h_rg][t];
                    ++house_notype_distribution_count[t];
                    ++occ_distribution[type][from_state+blue_state][t];
                    ++occ_distribution_count[type][t];
                    ++occ_notype_distribution[from_state+blue_state][t];
                    ++occ_notype_distribution_count[t];
                    if (HUSBAND_EDU_LEVEL != -1)
                    {
                        ++house_notype_edu_distribution[HUSBAND_EDU_LEVEL][from_h_rg][t];
                        ++house_notype_edu_distribution_count[HUSBAND_EDU_LEVEL][t]; 
                        ++occ_notype_edu_distribution[HUSBAND_EDU_LEVEL][from_state+blue_state][t];
                        ++occ_notype_edu_distribution_count[HUSBAND_EDU_LEVEL][t];
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
#ifdef SIMULATION
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
#ifdef SIMULATION
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
                            if (blue_state == FULL)
                            {
                                last_wage[draw] = (b_wage_flag == false) ? wage_b[house_rg_arr[PERIODS-1][draw]] : wage_b_non_f[house_rg_arr[PERIODS-1][draw]]/6.0;
                            }
                            else if (blue_state == PART)
                            {
                                last_wage[draw] = part_wage_factor*(((b_wage_flag == false) ? wage_b[house_rg_arr[PERIODS-1][draw]] : wage_b_non_f[house_rg_arr[PERIODS-1][draw]])/6.0);
                            }
                            else
                            {
                                // invalid state
                                assert(0);
                            }
#ifdef TRACE
#ifdef SIMULATION
                            if (IND_FILTER==1)
#endif
                            {
                                wage_blue_rg_sum[type][from_h_rg][blue_state] += last_wage[draw];
                                ++wage_blue_rg_count[type][from_h_rg][blue_state];
                                wage_blue_rg_notype_sum[from_h_rg][blue_state] += last_wage[draw];
                                ++wage_blue_rg_notype_count[from_h_rg][blue_state];
                                wage_blue_sum[type][blue_state] += last_wage[draw];
                                ++wage_blue_count[type][blue_state];
                                wage_blue_notype_sum[blue_state] += last_wage[draw];
                                ++wage_blue_notype_count[blue_state];
                            }
#endif
                        }
                        else if (from_state == UE)
                        {
                            // unemployment
                            last_wage[draw] = 0.0;
                        }
                        else
                        {
                            // invalid state
                            assert(0);
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
                        // TODO part/full time in last wage
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
               
                    bool nothing_chosen = true;
                    for (unsigned short st = 0; st < STATE_VECTOR_SIZE; ++st)
                    {
                        if (choices[st] > -INFINITY)
                        {
                            p_bar_arr[st][t] += (float)(dvtau[st]/dvsum);
                            nothing_chosen = false;
                        }
                        if (draw == draws_f-1)
                        {
                            // in the last draw calculate the average
                            p_bar_arr[st][t] = p_bar_arr[st][t]/(float)draws_f;
                        }
                    }
                    assert(!nothing_chosen);
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
                                    if (max_index_arr[t][draw] == 21+rg)
                                    {
                                        p_error[t] += error_c+(1.0f - error_c)*p_bar_arr[21+rg][t];
                                    }
                                    else
                                    {
                                        p_error[t] += (1.0f - error_c)*p_bar_arr[21+rg][t];
                                    }
                                }
                            }
                            else if (work(I,t) == 2)
                            {
                                // unknown where he live & work in 2
                                for (unsigned short rg = 0; rg < RG_SIZE; ++rg)
                                {
                                    if (max_index_arr[t][draw] == 35+rg)
                                    {
                                        p_error[t] += error_c+(1.0f - error_c)*p_bar_arr[35+rg][t];
                                    }
                                    else
                                    {
                                        p_error[t] += (1.0f - error_c)*p_bar_arr[35+rg][t];
                                    }
                                }
                            }
                            else if (work(I,t) == 3)
                            {
                                // unknown where he live & work in 3
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
                            else if (work(I,t) == 4)
                            {
                                // unknown where he live & work in 4
                                for (unsigned short rg = 0; rg < RG_SIZE; ++rg)
                                {
                                    if (max_index_arr[t][draw] == 49+rg)
                                    {
                                        p_error[t] += error_c+(1.0f - error_c)*p_bar_arr[49+rg][t];
                                    }
                                    else
                                    {
                                        p_error[t] += (1.0f - error_c)*p_bar_arr[49+rg][t];
                                    }
                                }
                            }
                            else if (work(I,t) == 6)
                            {
                                // unknown where he live & work in 6
                                for (unsigned short rg = 0; rg < RG_SIZE; ++rg)
                                {
                                    if (max_index_arr[t][draw] == 63+rg)
                                    {
                                        p_error[t] += error_c+(1.0f - error_c)*p_bar_arr[63+rg][t];
                                    }
                                    else
                                    {
                                        p_error[t] += (1.0f - error_c)*p_bar_arr[63+rg][t];
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
                                        if (max_index_arr[t][draw] == 21+h_rg+7*w_rg)
                                        {
                                            p_error[t] += error_c + (1.0f - error_c)*p_bar_arr[21+h_rg+7*w_rg][t];
                                        }
                                        else
                                        {
                                            p_error[t] += (1.0f - error_c)*p_bar_arr[21+h_rg+7*w_rg][t];
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
                                    p_error[t] += (1.0f - error_c)*p_bar_arr[live(I,t) + 7*rg][t];
                                }
                            }
                        }
                        else
                        {
                            printf("handle_missing_state error(W): work/live contradict sample (I=%hu t=%hu work=%hu live=%hu sample=%d)\n",  
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
                            printf("handle_missing_state error(UE): live contradict sample (I=%hu t=%hu live=%hu sample=%d)\n",
                                I, t, live(I,t), sample(I,t));  
                        }
                    }
                    else if (occupation(I,t) == BLUE && sample(I,t) < 14)
                    {
                        // work in blue - full
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
                            printf("handle_missing_state error(B, full): live contradict sample (I=%hu t=%hu live=%hu sample=%d)\n",
                                I, t, live(I,t), sample(I,t));
                        }

                    }
                    else if (occupation(I,t) == BLUE && sample(I,t) > 13)
                    {
                        // work in blue - partial
                        if (live(I,t) == -1)
                        {
                            //region of housing is unknown
                            for (unsigned short rg = 0; rg < RG_SIZE; ++rg)
                            {
                                if (max_index_arr[t][draw] == 14+rg)
                                {
                                    p_error[t] += error_c + (1.0f - error_c)*p_bar_arr[14+rg][t];
                                }
                                else
                                {
                                    p_error[t] += (1.0f - error_c)*p_bar_arr[14+rg][t];
                                }
                            }
                        }
                        else
                        {
                             printf("handle_missing_state error(B, partial): live contradict sample (I=%hu t=%hu live=%hu sample=%d)\n",
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
        printf("calculating forward for I = %d took: %f seconds (estimated total = %f minutes)\n", I, toc(tv), toc(tv)*OBSR/60.0);
#endif

#endif // not SIMULATION

    } // end I

#ifdef INFO
    printf("counter true: %lu counter false: %lu correct: %f%%\n", counter_true, counter_false, 100.0f*(float)counter_true/(float)(counter_true+counter_false));
#endif
    double likelihood = 0.0;
    printf("ID\tMarried\tLikelihood\n");
    for (unsigned short I = 0; I < OBSR; ++I)
    {
        const double prob = like_arr[I];
        double log_prob;
        if (prob < 1e-300 || prob == -INFINITY  || prob == INFINITY || prob == -NAN || prob == NAN)
        {
#ifdef INFO
            printf("failed to calculated value for I=%hu\n******************************************************\n", I);
#endif
            log_prob = -110.0;
        }
        else
        {
            log_prob = log(prob);
        }
        likelihood += log_prob;
        printf("%u\t%d\t%f\n", I, M_arr[I], log_prob);
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

    ////////////////////// Occupation Distribution /////////////////////
    printf("\n\noccupation distribution:\n\n");
    printf("---------------------------------------------------------------------------------\n");
    printf(" T |   count   |     UE      |      WHITE     |   BLUE (FULL)  |   BLUE (PART)  |\n");
    printf("---------------------------------------------------------------------------------\n");

#ifdef SIMULATION
    const unsigned short max_T = ((sim_type !=0) ? T : MOMENTS_PERIODS);
#else
    const unsigned short max_T = MOMENTS_PERIODS;
#endif

    for (unsigned short t = 0; t < max_T ; ++t)
    {
        printf("%hu\t%lu\t", t, occ_notype_distribution_count[t]);
        for (unsigned short st = 0; st < ALL_STATE_SIZE; ++st)
        {
            printf("%f\t", (float)occ_notype_distribution[st][t]/(float)occ_notype_distribution_count[t]);
        }

        printf("\n");
    }
    printf("---------------------------------------------------------------------------------\n");
    memset(occ_distribution_count, '\0', sizeof(occ_distribution_count));
    memset(occ_distribution, '\0', sizeof(occ_distribution));
    for (unsigned short I = 0; I < OBSR; ++I)
    {
        unsigned short last_t = PERIODS_arr[I];
        for (unsigned short t = 0; t < last_t; ++t)
        {
            if (occupation(I,t) > -1)
            {
                const short work_rg = sample(I,t)/7;
                ++occ_distribution_count[0][t];
                ++occ_distribution[0][(work_rg == 0) ? UE : ((work_rg == 1) ? BLUE : ((work_rg == 2) ? BLUE+1 : WHITE))][t];
            }
        }
    }

    for (unsigned short t = 0; t < MOMENTS_PERIODS ; ++t)
    {
        printf("%hu\t%lu\t", t, occ_distribution_count[0][t]);
        for (unsigned st = 0; st < ALL_STATE_SIZE; ++st)
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
    unsigned long sum_count_267 = 0;
    float sum_wage_267 = 0.0f;
    printf("\n\naverage white wage in last period:\n\n");
    printf("---------------------------------------------------------------------------------------------------\n");
    printf(" ty |      1       |       3       |       4       |       5       |     267      |    average    |\n");
    printf("---------------------------------------------------------------------------------------------------\n");
    for (unsigned short ty = 0; ty < TYPE_SIZE; ++ty)
    {
        printf("%hu\t", ty);
        unsigned long sum_count = 0;
        float sum_wage = 0.0f;
        sum_count_267 = 0;
        sum_wage_267 = 0.0f;
        for (unsigned short rg = 0; rg  < RG_SIZE; ++rg )
        {
            if (rg == 1 || rg == 5 || rg == 6) // the 267 region
            {
                sum_count_267 += wage_white_rg_count[ty][rg];
                sum_wage_267 += wage_white_rg_sum[ty][rg];
                if (rg == 6)
                {
                    if (sum_count_267 > 0) 
                    {
                        printf("%.3f\t", sum_wage_267/(float)sum_count_267);
                    }
                    else 
                    {
                        printf("--------\t");
                    }
                }
            }
            else // other regions
            {
                if (wage_white_rg_count[ty][rg] > 0)
                {
                    printf("%.3f\t", wage_white_rg_sum[ty][rg]/(float)wage_white_rg_count[ty][rg]);
                }
                else
                {
                    printf("--------\t");
                }
            }
            sum_count += wage_white_rg_count[ty][rg];
            sum_wage += wage_white_rg_sum[ty][rg];
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
    printf("---------------------------------------------------------------------------------------------------\n");
    printf("avg\t");
    sum_count_267 = 0;
    sum_wage_267 = 0.0f;
    for (unsigned short rg = 0; rg  < RG_SIZE; ++rg )
    {
        if (rg == 1 || rg == 5 || rg == 6) // the 267 region
        {
            sum_count_267 += wage_white_notype_rg_count[rg];
            sum_wage_267 += wage_white_notype_rg_sum[rg];
            if (rg == 6)
            {
                printf("%.3f\t", sum_wage_267/(float)sum_count_267);
            }
        }
        else // other regions
        {
            printf("%.3f\t", wage_white_notype_rg_sum[rg]/(float)wage_white_notype_rg_count[rg]);
        }
    }

    // average across all regions
    printf("%.3f\t", wage_white_notype_sum/(float)wage_white_notype_count);

    // real values
    printf("\n-------------------------------------------------------------------------------------------------\n");
    memset(wage_white_rg_sum, '\0', sizeof(wage_white_rg_sum));
    memset(wage_white_rg_count, '\0', sizeof(wage_white_rg_count));
    for (unsigned short I = 0; I < OBSR; ++I)
    {
        if (occupation(I, PERIODS_arr[I]-1) == WHITE && work(I, PERIODS_arr[I]-1) != -1 && WAGE_arr[I] > 0.0f)
        {
            wage_white_rg_sum[0][work(I, PERIODS_arr[I]-1)] += WAGE_arr[I];
            ++wage_white_rg_count[0][work(I, PERIODS_arr[I]-1)];
        }
    }
    {
        printf("real\t");
        unsigned long sum_count = 0;
        float sum_wage = 0.0f;
        sum_count_267 = 0;
        sum_wage_267 = 0.0f;
        for (unsigned short rg = 0; rg  < RG_SIZE; ++rg )
        {
            if (wage_white_rg_count[0][rg] > 0)
            {
                if (rg == 1 || rg == 5 || rg == 6) // the 267 region
                {
                    sum_count_267 += wage_white_rg_count[0][rg];
                    sum_wage_267 += wage_white_rg_sum[0][rg];
                    if (rg == 6)
                    {
                        printf("%.3f\t", sum_wage_267/(float)sum_count_267);
                    }
                }
                else // other regions
                {
                    printf("%.3f\t", wage_white_rg_sum[0][rg]/(float)wage_white_rg_count[0][rg]);
                }
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
    
    ////////////////////// Last Blue Wage Full /////////////////////
    printf("\n\naverage blue (full time) wage in last period:\n\n");
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
            if (wage_blue_rg_count[ty][rg][FULL] > 0)
            {
                printf("%.4f\t", wage_blue_rg_sum[ty][rg][FULL]/(float)wage_blue_rg_count[ty][rg][FULL]);
                sum_count += wage_blue_rg_count[ty][rg][FULL];
                sum_wage += wage_blue_rg_sum[ty][rg][FULL];
            }
            else
            {
                printf("--------\t");
            }
        }
        if (sum_count > 0)
        {
            printf("%.4f\n", sum_wage/(float)sum_count);
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
        printf("%.4f\t",  wage_blue_rg_notype_sum[rg][FULL]/(float)wage_blue_rg_notype_count[rg][FULL]);
    }
    
    // average across all regions
    printf("%.4f\t", wage_blue_notype_sum[FULL]/(float)wage_blue_notype_count[FULL]);

    // real values
    printf("\n------------------------------------------------------------------------------------------------------------------------------------\n");
    float real_wage_blue_rg_sum[TYPE_SIZE][RG_SIZE];
    memset(real_wage_blue_rg_sum, '\0', sizeof(real_wage_blue_rg_sum));
    unsigned long real_wage_blue_rg_count[TYPE_SIZE][RG_SIZE]={{0}};
    for (unsigned short I = 0; I < OBSR; ++I)
    {
        if (occupation(I, PERIODS_arr[I]-1) == BLUE && live(I, PERIODS_arr[I]-1) != -1 && WAGE_arr[I] > 0.0f && sample(I,PERIODS_arr[I]-1)/7 == 1)
        {
            real_wage_blue_rg_sum[0][live(I, PERIODS_arr[I]-1)] += WAGE_arr[I];
            ++real_wage_blue_rg_count[0][live(I, PERIODS_arr[I]-1)];
        }
    }
    {
        printf("real\t");
        unsigned long sum_count = 0;
        float sum_wage = 0.0;
        for (unsigned short rg = 0; rg  < RG_SIZE; ++rg )
        {
            if (real_wage_blue_rg_count[0][rg] > 0)
            {
                printf("%.4f\t", real_wage_blue_rg_sum[0][rg]/(float)real_wage_blue_rg_count[0][rg]);
                sum_count += real_wage_blue_rg_count[0][rg];
                sum_wage += real_wage_blue_rg_sum[0][rg];
            }
            else
            {
                printf("--------\t");
            }
        }
        if (sum_count > 0)
        {
            printf("%.4f\n", sum_wage/(float)sum_count);
        }
        else
        {
            printf("--------\n");
        }
    }

    ////////////////////// Last Blue Wage Part /////////////////////
    printf("\n\naverage blue (part time) wage in last period:\n\n");
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
            if (wage_blue_rg_count[ty][rg][PART] > 0)
            {
                printf("--------\t");
                sum_count += wage_blue_rg_count[ty][rg][PART];
                sum_wage += wage_blue_rg_sum[ty][rg][PART];
            }
            else
            {
                printf("--------\t");
            }
        }
        if (sum_count > 0)
        {
            printf("%.4f\n", sum_wage/(float)sum_count);
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
        printf("--------\t");
    }
    
    // average across all regions
    printf("%.4f\t", wage_blue_notype_sum[PART]/(float)wage_blue_notype_count[PART]);

    // real values
    printf("\n------------------------------------------------------------------------------------------------------------------------------------\n");
    memset(real_wage_blue_rg_sum, '\0', sizeof(real_wage_blue_rg_sum));
    memset(real_wage_blue_rg_count, '\0', sizeof(real_wage_blue_rg_count));
    for (unsigned short I = 0; I < OBSR; ++I)
    {
        if (occupation(I, PERIODS_arr[I]-1) == BLUE && live(I, PERIODS_arr[I]-1) != -1 && WAGE_arr[I] > 0.0f && sample(I,PERIODS_arr[I]-1)/7 == 2)
        {
            real_wage_blue_rg_sum[0][live(I, PERIODS_arr[I]-1)] += WAGE_arr[I];
            ++real_wage_blue_rg_count[0][live(I, PERIODS_arr[I]-1)];
        }
    }
    {
        printf("real\t");
        unsigned long sum_count = 0;
        float sum_wage = 0.0;
        for (unsigned short rg = 0; rg  < RG_SIZE; ++rg )
        {
            if (real_wage_blue_rg_count[0][rg] > 0)
            {
                printf("%.4f\t", real_wage_blue_rg_sum[0][rg]/(float)real_wage_blue_rg_count[0][rg]);
                sum_count += real_wage_blue_rg_count[0][rg];
                sum_wage += real_wage_blue_rg_sum[0][rg];
            }
            else
            {
                printf("--------\t");
            }
        }
        if (sum_count > 0)
        {
            printf("%.4f\n", sum_wage/(float)sum_count);
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
        //unsigned short last_t = PERIODS_arr[I];
        for (unsigned short t = 0; t < MOMENTS_PERIODS ; ++t)
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
#endif // RANDOM_SELECTION

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
     
    ////////////////////// House Work Distribution //////////////////////
    printf("\n\nhouse-work region distribution (all types):\n\n");
    printf("----------------------------------------------------------------------------------------------------------------------------\n");
    printf(" rg |   count   |      1       |      2        |      3        |       4        |      5       |       6      |      7      |\n");
    printf("----------------|              |               |               |                |              |              |             |\n");
    printf(" house \\ work   |              |               |               |                |              |              |             |\n");
    printf("----------------------------------------------------------------------------------------------------------------------------\n");
    
    for (unsigned short h_rg = 0; h_rg < RG_SIZE; ++h_rg)
    {
        printf("%hu\t%lu\t", (unsigned short)(h_rg+1), house_work_rg_notype_distribution_count[h_rg]);
        for (unsigned short w_rg = 0; w_rg < RG_SIZE; ++w_rg)
        {
            printf("%f\t", (float)house_work_rg_notype_distribution[h_rg][w_rg]/(float)house_work_rg_notype_distribution_count[h_rg]);
        }
        printf("\n");
    }
    printf("\n------------------------------------------------------------------------------------------------------------------------------------\n");
    printf("%d\t-----\t57.63000\t3.390000\t38.98000\t0.000000\t0.000000\t0.000000\t0.000000\n", 1);
    printf("%d\t-----\t47.47000\t30.30000\t14.14000\t8.080000\t0.000000\t0.000000\t0.000000\n", 2);
    printf("%d\t-----\t25.96000\t5.140000\t58.10000\t4.110000\t3.340000\t0.510000\t2.830000\n", 3);
    printf("%d\t-----\t0.000000\t0.000000\t0.000000\t90.58000\t9.420000\t0.000000\t0.000000\n", 4);
    printf("%d\t-----\t0.000000\t1.310000\t3.270000\t32.03000\t63.40000\t0.000000\t0.000000\n", 5);
    printf("%d\t-----\t1.880000\t0.000000\t9.380000\t0.000000\t0.000000\t88.75000\t0.000000\n", 6);
    printf("%d\t-----\t0.000000\t0.000000\t0.000000\t0.000000\t0.000000\t0.000000\t100.0000\n", 7);

#ifdef PRINT_OCC_DIST
    ////////////////////// Occupation Distribution per Education Level//////////////////////
    for (unsigned int edu_level = 0; edu_level < EDU_LEVELS; ++edu_level)
    {
        printf("\n\noccupation distribution: for education level %hu-%hu:\n\n", edu_lower[edu_level], edu_upper[edu_level]);
        printf("---------------------------------------------------------------------------------\n");
        printf(" T |   count   |     UE      |      WHITE     |   BLUE (FULL)  |   BLUE (PART)  |\n");
        printf("---------------------------------------------------------------------------------\n");
        for (unsigned short t = 0; t < max_T; ++t)
        {
            if (t != 1 && t !=3 && t !=5 && t != 9)
            {
                continue;
            }
            printf("%hu\t%lu\t", t, occ_notype_edu_distribution_count[edu_level][t]);
            for (unsigned st = 0; st < ALL_STATE_SIZE; ++st)
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
                if (occupation(I,t) > -1 && HUSBAND_EDU_arr[I] != 99 && HUSBAND_EDU_arr[I] >= edu_lower[edu_level] && HUSBAND_EDU_arr[I] <= edu_upper[edu_level])
                {
                     const short work_rg = sample(I,t)/7;
                    ++occ_distribution_count[0][t];
                    ++occ_distribution[0][(work_rg == 0) ? UE : ((work_rg == 1) ? BLUE : ((work_rg == 2) ? BLUE+1 : WHITE))][t];
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
            for (unsigned st = 0; st < ALL_STATE_SIZE; ++st)
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
#endif // PRINT_OCC_DIST

#ifdef PRINT_EDU_LEVEL
    ////////////////////// House Region Distribution per Education Level//////////////////////
    for (unsigned int edu_level = 0; edu_level < EDU_LEVELS; ++edu_level)
    {
        printf("\n\nhousing region distribution (all types) for education level %hu-%hu:\n\n", edu_lower[edu_level], edu_upper[edu_level]);
        printf("-----------------------------------------------------------------------------------------------------------------------------\n");
        printf(" T |   count    |      1       |      2        |      3        |       4        |      5       |       6      |      7      |\n");
        printf("-----------------------------------------------------------------------------------------------------------------------------\n");

        for (unsigned short t = 0; t < max_T; ++t)
        {
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
                if (live(I,t) > -1 &&  HUSBAND_EDU_arr[I] != 99 && HUSBAND_EDU_arr[I] >= edu_lower[edu_level] && HUSBAND_EDU_arr[I] <= edu_upper[edu_level])
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
    }
	
    ////////////////////// Work Region Distribution per education level//////////////////////
     for (unsigned int edu_level = 0; edu_level < EDU_LEVELS; ++edu_level)
    {
        printf("\n\nwork region distribution (all types) for education level %hu-%hu:\n\n", edu_lower[edu_level], edu_upper[edu_level]);
        printf("-----------------------------------------------------------------------------------------------------------------------------\n");
        printf(" T |   count    |      1       |      2        |      3        |       4        |      5       |       6      |      7      |\n");
        printf("-----------------------------------------------------------------------------------------------------------------------------\n");
    
        for (unsigned short t = 0; t < max_T; ++t)
        {
            printf("%hu\t%lu\t", t, work_rg_notype_edu_distribution_count[edu_level][t]);
            for (unsigned short rg = 0; rg < RG_SIZE; ++rg)
            {
                printf("%f\t", (float)work_rg_notype_edu_distribution[edu_level][rg][t]/(float)work_rg_notype_edu_distribution_count[edu_level][t]);
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
                if (work(I,t) > -1 && occupation(I,t) > -1 && occupation(I,t) == WHITE && 
                        HUSBAND_EDU_arr[I] != 99 && HUSBAND_EDU_arr[I] >= edu_lower[edu_level] && HUSBAND_EDU_arr[I] <= edu_upper[edu_level])
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
    }
	
	////////////////////// Occupation Distribution per Education Level (no period filter)//////////////////////
    for (unsigned int edu_level = 0; edu_level < EDU_LEVELS; ++edu_level)
    {
        printf("\n\noccupation distribution: for education level %hu-%hu:\n\n", edu_lower[edu_level], edu_upper[edu_level]);
        printf("---------------------------------------------------------------------------------\n");
        printf(" T |   count   |     UE      |      WHITE     |   BLUE (FULL)  |   BLUE (PART)  |\n");
        printf("---------------------------------------------------------------------------------\n");
        for (unsigned short t = 0; t < max_T; ++t)
        {
            printf("%hu\t%lu\t", t, occ_notype_edu_distribution_count[edu_level][t]);
            for (unsigned st = 0; st < ALL_STATE_SIZE; ++st)
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
                if (occupation(I,t) > -1 && HUSBAND_EDU_arr[I] != 99 && HUSBAND_EDU_arr[I] >= edu_lower[edu_level] && HUSBAND_EDU_arr[I] <= edu_upper[edu_level])
                {
                     const short work_rg = sample(I,t)/7;
                    ++occ_distribution_count[0][t];
                    ++occ_distribution[0][(work_rg == 0) ? UE : ((work_rg == 1) ? BLUE : ((work_rg == 2) ? BLUE+1 : WHITE))][t];
                }
            }
        }
        for (unsigned short t = 0; t < MOMENTS_PERIODS ; ++t)
        {
            printf("%hu\t%lu\t", t, occ_distribution_count[0][t]);
            for (unsigned st = 0; st < ALL_STATE_SIZE; ++st)
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
#endif // PRINT_EDU_LEVEL
	
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
static const char* IND_DATA_FILENAME = "ind_data_3.txt";
static const char* MOMENTS_FILENAME = "olim_wide_3.txt";
static const char* HUSBAND_EDU_FILENAME = "husband_edu.txt";
static const char* USSR_ENG_EXP_FILENAME = "ussr_eng_exp.txt";
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
            for (unsigned int I = 0; I < OBSR; ++I)
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

    if (!load_husband_edu(HUSBAND_EDU_FILENAME))
    {
        fprintf(stderr, "failed to load husband education data from file %s\n", HUSBAND_EDU_FILENAME);
        return -1;
    }

    if (!load_ussr_eng_exp(USSR_ENG_EXP_FILENAME))
    {
        fprintf(stderr, "failed to load USSR engineer experience data from file %s\n", USSR_ENG_EXP_FILENAME);
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
    if (dynamic_param_size == 0)
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
