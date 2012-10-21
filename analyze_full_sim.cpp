#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

const unsigned short LINE_MAX = 128;
const unsigned short T = 20;
const unsigned short RG_SIZE = 7;
const unsigned short STATE_SIZE = 3;
const unsigned short UE = 0;
const unsigned short WHITE = 1;
const unsigned short BLUE = 2;

int main(int argc, char** argv)
{
    if (argc != 2)
    {
        fprintf(stderr, "usage: %s <full simulation filename>\n", argv[0]);
        return -1;
    }

    const char* filename = argv[1];

    FILE* fp = fopen(filename,"r");
    if (fp)
    {
        char *line = (char*)malloc(LINE_MAX);
        int col_num = 0;
        unsigned long I = 0;
        unsigned short obs;
        unsigned short type;
        unsigned short draw;
        unsigned short t[T];
		unsigned short work_rg;
		unsigned short house_rg;
        unsigned long t_count[T]={0};
		unsigned long work_rg_count[T][RG_SIZE]={{0}};
		unsigned long house_rg_count[T][RG_SIZE]={{0}};
		unsigned long occupation_count[T][STATE_SIZE]={{0}};

        while (fgets(line, LINE_MAX, fp) != 0)
        {
            // formatting the line according to:
            // OBS TYPE DRAW T0 T1 T2 T3 T4 T5 T6 T7 T8 T9 T10 T11 T12 T13 T14 T15 T16 T17 T18 T19
            col_num = sscanf(line, "%hu%hu%hu%hu%hu%hu%hu%hu%hu%hu%hu%hu%hu%hu%hu%hu%hu%hu%hu%hu%hu%hu%hu",
                             &obs, &type, &draw, &t[0], &t[1], &t[2], &t[3], &t[4], &t[5], &t[6], &t[7], &t[8], &t[9],
                             &t[10], &t[11], &t[12], &t[13], &t[14], &t[15], &t[16], &t[17], &t[18], &t[19]);

            for (unsigned short tt = 0; tt < col_num-3; ++tt) 
            {
                div_t info = div(t[tt],RG_SIZE);
                house_rg = (unsigned short)info.rem;
                work_rg = (unsigned short)info.quot;
                // work region 0 = unemployment
                // work region 1 = blue
                // (work region - 2) = white work region
                house_rg_count[tt][house_rg]++;
                t_count[tt]++;
                if (work_rg == 0) 
                {
                    occupation_count[tt][UE]++;
                }
                else if (work_rg == 1) 
                {
                    occupation_count[tt][BLUE]++;
                }
                else
                {
                    work_rg_count[tt][work_rg-2]++;
                    occupation_count[tt][WHITE]++;
                }
            }

            ++I;
        }

        free(line);
        printf("finished reading %lu lines from file %s with status: %s\n", I, filename, strerror(errno));
        fclose(fp);

        ////////////////////// Occupation Distribution /////////////////////
        printf("\n\noccupation distribution:\n\n");
        printf("----------------------------------------------------------------\n");
        printf(" T |   count   |     UE      |      WHITE     |       BLUE     |\n");
        printf("----------------------------------------------------------------\n");

        for (unsigned short tt = 0; tt < T ; ++tt)
        {
            printf("%hu\t%lu\t", tt, t_count[tt]);
            for (unsigned short st = 0; st < STATE_SIZE; ++st) 
            {
                printf("%.3f\t\t", (float)occupation_count[tt][st]/(float)t_count[tt]);
            }
			printf("\n");
        }

        ////////////////////// House Region Distribution //////////////////////
        printf("\n\nhousing region distribution (all types):\n\n");
        printf("-----------------------------------------------------------------------------------------------------------------------------\n");
        printf(" T |   count    |      1       |      2        |      3        |       4        |      5       |       6      |      7      |\n");
        printf("-----------------------------------------------------------------------------------------------------------------------------\n");

        for (unsigned short tt = 0; tt < T; ++tt)
        {
            printf("%hu\t%lu\t", tt, t_count[tt]);
            for (unsigned rg = 0; rg < RG_SIZE; ++rg)
            {
                printf("%.3f\t\t", (float)house_rg_count[tt][rg]/(float)t_count[tt]);
            }
            printf("\n");
        }

        ////////////////////// Work Region Distribution //////////////////////
        printf("\n\nwork region distribution (all types):\n\n");
        printf("-----------------------------------------------------------------------------------------------------------------------------\n");
        printf(" T |   count    |      1       |      2        |      3        |       4        |      5       |       6      |      7      |\n");
        printf("-----------------------------------------------------------------------------------------------------------------------------\n");

        for (unsigned short tt = 0; tt < T; ++tt)
        {
            printf("%hu\t%lu\t", tt, occupation_count[tt][WHITE]);
            for (unsigned rg = 0; rg < RG_SIZE; ++rg)
            {
                printf("%.3f\t\t", (float)work_rg_count[tt][rg]/(float)occupation_count[tt][WHITE]);
            }
            printf("\n");
        }
	
    }
    else
    {
        printf("failed to open file %s\n", filename);
    }

    return 0;
}

